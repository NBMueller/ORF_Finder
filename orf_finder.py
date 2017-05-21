#!/usr/bin/env python3

""" This program takes a FASTA file, reads the sequence and searches for ORF\'s.
    The found ORF\'s can be BLASTed via webBLAST or standalone BLAST.
"""

import argparse
import logging
import os
import sys
import operator
import re
import time
import math
import logging
from urllib.request import urlopen, Request
from urllib.parse import urlencode
from io import StringIO
from datetime import datetime, timedelta
import warnings
import configparser
import tempfile
import subprocess
import queue
import threading
from subprocess import Popen
# External package libraries
import numpy as np
import matplotlib.pyplot as plt
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning) # Silence Biopython warnigns
# from Bio.Blast import NCBIWWW
# from Bio.Blast import NCBIXML
from Bio import SearchIO

__version__ = "1.0.2"
__author__ = "Nico Borgsmueller"
__contact__ = "nicoborgsmueller@web.de"
__credits__ = ["Andreas Andrusch"]
__copyright__ = "Copyright 2017, Nico Borgsmueller"
__license__ = "GNU GPLv3"
__status__ = "Development"
__date__ = "2017.03.24"

# Standalone BLAST Setup for Unix: https://www.ncbi.nlm.nih.gov/books/NBK52640/
# BLAST output specifications: https://www.ncbi.nlm.nih.gov/books/NBK279675/
# BLAST Command Line Applications User Manual: https://www.ncbi.nlm.nih.gov/books/NBK279690/

################################################################################
################################ ARGPARSER #####################################
################################################################################

parser = argparse.ArgumentParser(description='Find orf\'s in a FASTA file.')

parser.add_argument(
    'input_file',
    type=str,
    help="Path to Input file. can be either a FASTA file or a config file."
)
parser.add_argument(
    '-l', '--length',
    type=int,
    default=200,
    help='Minimal ORF length (default=60).'
)
parser.add_argument(
    '-s', '--start',
    type=str,
    default='AUG',
    help='Start codon (default="AUG").',   
)
parser.add_argument(
    '-p', '--stop',
    type=str,
    default='UAA,UAG,UGA',
    help='Stop codon (default="UAA,UAG,UGA").',   
)
parser.add_argument(
    '-oio', '--ORFinORF',
    action='store_true',
    help='Check for all ORFs inside other ORFs.',   
)
parser.add_argument(
    '-o', '--output',
    type=str,
    default=os.getcwd(),
    help='Path for output directory creation (default=<CURRENT_WORKING_DIR>)',   
)

parser.add_argument(
    '-b', '--blast',   
    action='store_true',
    help='Start new web Blast query for found ORFs.',
)
parser.add_argument(
    '-bDB', '--blastDB',   
    type=str,
    default='nt',
    choices=['refseq_rna', 'nt'],
    help='Database used for BLAST query (default="nt").',
)
parser.add_argument(
    '-bProg', '--blastProgram',   
    type=str,
    default='blastn',
    choices=['blastn', 'megablast'],
    help='Program used for BLAST query (default="blastn").',
)
parser.add_argument(
    '-bqn', '--blastQueryNumber',   
    type=int,
    default=20,
    help='Number of queries used for BLASTing results (default=20).'
)
parser.add_argument(
    '-brn', '--blastResultNumber',   
    type=int,
    default=1,
    help='Number of BLAST results saved for each ORF.'
)

################################################################################
################### FIXED VARIABLES AND HELPER CLASSES #########################
################################################################################

# Stop Codon = ?
AA_MAP = {
    "UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
    "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
    "UAU":"Y", "UAC":"Y", "UAA":"?", "UAG":"?",
    "UGU":"C", "UGC":"C", "UGA":"?", "UGG":"W",
    "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
    "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
    "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
    "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G"
}

# URL for running Web BLAST Queries
BLAST_BASE_URL = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi' 

# Logging class for redirecting stderr to logfile
class StreamToLogger(object):
    """
    Fake file-like stream object that redirects writes to a logger instance.
    """
    def __init__(self, logger, log_level=logging.INFO):
      self.logger = logger
      self.log_level = log_level
      self.linebuf = ''

    def write(self, buf):
      for line in buf.rstrip().splitlines():
         self.logger.log(self.log_level, line.rstrip())

################################################################################
################################# FASTA  #######################################
################################################################################

class FASTA:
    def __init__(self, path_to_fasta):
        self.file = path_to_fasta
        self.info, self.raw_seq = self.read_file()
        self.forward_seq, self.reverse_seq = self.process_sequence()


    def read_file(self):
        fasta_content = open(self.file, 'r')

        # First line should start with '>' and contains additional information
        info_content = fasta_content.readline().rstrip()
        if info_content[0] != '>':
            raise IOError(
                'The input FASTA file does\'t start with ">".'
            )
        else:
            info_content.lstrip('>')

        # Rest of the file is expected to be the actual sequence    
        seq = fasta_content.read().replace('\n', '')
        fasta_content.close()
        if seq.count('>') != 0:
            raise IOError(
                'The input FASTA file contains more than one sequence.'
            )

        return info_content, seq


    def process_sequence(self):
        # Convert all lower characters to upper ones
        upper_seq = self.raw_seq.upper()
        no_n_seq = upper_seq.replace('N', '')
        new_seq = no_n_seq.replace('T', 'U').replace(' ', '')
        return new_seq, self._get_compl_seq(new_seq[::-1])


    def _get_compl_seq(self, seq):
        compl_dict = {'A': 'U','C': 'G','G': 'C','U': 'A'}
        return ''.join([compl_dict[i] for i in seq])


################################################################################
################################### ORF ########################################
################################################################################


class ORF:
    def __init__(self, seq, start_pos, stop_pos, strand, frame):
        self.seq = seq
        self.AA_seq = self._get_AA_seq(seq)
        self.start = start_pos
        self.stop = stop_pos
        self.length = self.stop - self.start
        self.strand = strand
        self.frame = frame

        self.blast = []
        self.outsideORF = None
        self.insideORF = []


    def _get_AA_seq(self, seq):
        codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
        return ''.join([AA_MAP[i] for i in codons][:-1])


    def set_insideORF(self, ORF):
        self.insideORF.append(ORF)


    def set_outsideORF(self, ORF):
        self.outsideORF = ORF
        ORF.set_insideORF(self)


    def set_blast_result(self, result):
        self.blast.append(result)


    def _get_insideORF_str(self):
        if not self.insideORF:
            return ''
        else:
            msg = ','.join(
                ['({}-{})'.format(i.start, i.stop) for i in self.insideORF]
            )
            return msg


    def _get_outsideORF_str(self):
        if not self.outsideORF:
            return ''
        else:
            return '({}-{})'.format(self.outsideORF.start, self.outsideORF.stop)



    def _get_blast_result(self):
        if not self.blast:
            return False
        else:
            results_str = '\t'.join([i.to_output() for i in self.blast])
            return results_str


    def get_out_str(self):
        out_str_basic = '{s}\t{f}\t{st}\t{sp}\t{l}\t{seq}\t{prot}\t{encl}\t{embed}' \
            .format(
                s=self.strand,
                f=self.frame,
                st=self.start,
                sp=self.stop,
                l=self.length,
                seq=self.seq.replace('U', 'T'),
                prot=self.AA_seq,
                encl=self._get_outsideORF_str(),
                embed=self._get_insideORF_str(),
            )
        out_str_blast = self._get_blast_result()
        if out_str_blast:
            out_str = '{}\t{}\n'.format(out_str_basic, out_str_blast)
        else:
            out_str = '{}\n'.format(out_str_basic)
        return out_str


    def summary(self):
        if self.outsideORF:
            outside_ORF_msg = '({} - {})' \
                .format(self.outsideORF.start, self.outsideORF.stop)
        else:
            outside_ORF_msg = '-'

        if self.insideORF:
            inside_ORF_msg = ','.join(
                ['({} - {})'.format(i.start, i.stop) for i in self.insideORF]
            )
        else:
            inside_ORF_msg = '-'

        summary_msg = """
        Strand:\t{}
        Frame:\t{}
        Start:\t{}
        Stop:\t{}
        Length:\t{}
        
        Sequence:\t{}
        Protein:\t{}
    
        Enclosing ORF:\t{}
        Contained ORFS:\t{} 
        """.format(
                self.strand, self.frame, self.start, self.stop, self.length,
                self.seq, self.AA_seq,
                outside_ORF_msg,
                inside_ORF_msg,
            )
        print(summary_msg)


################################################################################
################################ ORF FINDER ####################################
################################################################################


class ORFFinder:
    def __init__(self, fasta, min_length, start_codons, stop_codons, ORFinORF):
        self.fasta = fasta
        self.min_length = min_length
        self.start_codons = start_codons
        self.stop_codons = stop_codons
        self.ORFinORF = ORFinORF

        self.statistics = {
            '+': {'orfs': {0: 0, 1: 0, 2: 0},
                'OiO': 0},
            '-': {'orfs': {0: 0,1: 0,2: 0},
                'OiO': 0}
        }

        pretty_out('\tScanning forward sequence...')
        forward_orfs = self._get_orfs(fasta.forward_seq, '+')
        pretty_out('\tScanning forward sequence... done')

        pretty_out('\tScanning reverse sequence...')
        reverse_orfs = self._get_orfs(fasta.reverse_seq, '-')
        pretty_out('\tScanning reverse sequence... done')

        if self.ORFinORF:
            pretty_out('Scanning for ORF\'s inside other ORF\'s...')
            all_orfs_sorted = {
                '+': forward_orfs,
                '-': reverse_orfs
            }
            self.find_OiOs(all_orfs_sorted)
            pretty_out(
                'Scanning for ORF\'s inside other ORF\'s... done', end='\n'
            )

        self.all_orfs = {
            '+': dict(forward_orfs),
            '-': dict(reverse_orfs)
        }


    def _get_orfs(self, seq, strand):
        # Find all start positions and stop positions
        orfs = {0: {}, 1: {}, 2: {}}
        start_pos = {0: [], 1:[], 2:[]}
        for i in range(len(seq)-2):
            # Start codon found
            if seq[i:i+3] in self.start_codons:
                # Safe in right frame list
                start_pos[i%3].append(i)
            # Stop codon found
            elif seq[i:i+3] in self.stop_codons:
                frame = i%3
                outside_orfs = []
                for start_codon_pos in start_pos[frame]:
                    # check if orf is long enough
                    # too short                   
                    if i+3 - start_codon_pos < self.min_length:
                        start_pos[frame] = []
                        break                 
                    new_orf = ORF(
                        seq[start_codon_pos:i+3],
                        start_codon_pos,
                        i+3,
                        strand,
                        '+{}'.format(frame)
                    )                   
                    # Set outside ORF is Flag is true
                    if self.ORFinORF:
                        for outside_orf in outside_orfs:
                            new_orf.set_outsideORF(outside_orf)
                        # Append current ORF to be an possible outside ORF as well
                        outside_orfs.append(new_orf)
                    # Store new ORF in ORF dict
                    orfs[frame][start_codon_pos] = new_orf

        for i,j in orfs.items():
            pretty_out(
                '\t({}{})-strand:\t{:5} ORF\'s'.format(strand, i, len(j))
            )
            self.statistics[strand]['orfs'][i] = len(j)

        # Sort all ORF's by start position
        all_orfs = {}
        for frame_orf in orfs.values():
            all_orfs.update(frame_orf)
        sorted_orfs = sorted(all_orfs.items(), key=operator.itemgetter(0))

        return sorted_orfs


    def get_orfs_list(self):
        orfs_list = []
        for strand, strand_orfs in self.all_orfs.items():
            orfs_list.extend(list(strand_orfs.values()))
        return(orfs_list)



    def find_OiOs(self, all_orfs_sorted):
        for strand, strand_orfs in all_orfs_sorted.items():
            OiO_set = set()
            for orf_idx, orf_details in enumerate(strand_orfs[:-1]):
                orf = orf_details[1]

                next_orf_idx = 1
                next_orf = strand_orfs[orf_idx+next_orf_idx][1]
                # CHeck if stop position from next ORF is smaller than current
                # Stop position --> ORF inside ORF
                while (next_orf.stop <= orf.stop) \
                        and (orf_idx + next_orf_idx < len(strand_orfs)-1):
                    orf.set_insideORF(next_orf)
                    next_orf.set_outsideORF(orf)
                    
                    OiO_set.add(orf)
                    next_orf_idx += 1

                    next_orf = strand_orfs[orf_idx+next_orf_idx][1]

            self.statistics[strand]['OiO'] = len(OiO_set)
            pretty_out(
                '({})-strand: {} ORF\'s inside other ORF\'s found.' \
                    .format(strand, len(OiO_set))
            )


    def write_output(
            self,
            results_dir_path,
            BLASTprogram,
            BLASTdb):
        new_file = os.path.join(results_dir_path, 'ORFs.tsv')
        pretty_out(
            'Writing results to: {}'.format(new_file),
            linebreaks=False
        )
        out_file = open(new_file, 'w')

        out_file.write('# FASTA: {}\n'.format(self.fasta.file))
        out_file.write('# FASTA Description: {}\n'.format(self.fasta.info[1:]))
        out_file.write('# Min ORF length: {}\n'.format(self.min_length))
        out_file.write('# Start Codons: {}\n'.format(','.join(self.start_codons)))
        out_file.write('# Stop Codons: {}\n'.format(','.join(self.stop_codons)))
        out_file.write('# BLAST Program: {}\n'.format(BLASTprogram))
        out_file.write('# BLAST DB: {}\n'.format(BLASTdb))

        all_orfs = []
        for strand, strand_orfs in self.all_orfs.items():
            all_orfs.extend(list(strand_orfs.values()))

        header_basics = 'Strand\tFrame\tStart\tStop\tLength\tSequence\tProtein\tEnclosing_ORF\tEmbedded_ORFs'
        header_blast = '\t'.join([
            'BLAST_accession_{i}\tBLAST_Evalue_{i}\tBLAST_description_{i}' \
                .format(i=idx) \
            for idx in range(1, len(all_orfs[0].blast)+1)
        ])
        out_file.write(
            '{}\t{}\n'.format(header_basics, header_blast)
        )

        length_list = {'AA': [], 'bp': []}

        all_orfs = []
        for strand, strand_orfs in self.all_orfs.items():
            all_orfs.extend(list(strand_orfs.values()))
        for orf in all_orfs:        
            length_list['bp'].append(len(orf.seq))
            length_list['AA'].append(len(orf.AA_seq))
            out_file.write(orf.get_out_str())

        out_file.close()

        # Statistics out
        self._write_statics_out(results_dir_path, length_list)
        self._plot_histogram(results_dir_path, length_list)

        pretty_out('Writing results ...done')


    def _plot_histogram(self, dir_path, length_list):
        for orf_type, type_data in length_list.items():
            new_pic = os.path.join(
                dir_path,
                '{}_histogram.png'.format(orf_type)
            )
            fig = plt.figure()
            ax = fig.add_subplot(111)
            bin_no = int(np.sqrt(len(type_data)))
            ax.hist(type_data, bin_no, facecolor='blue', rwidth=0.9)
            ax.set_title('Length Distribution of ORFs in {}'.format(orf_type))    
            ax.set_xlabel('Length [{}]'.format(orf_type))
            ax.grid()
            ax.set_ylabel('Count')

            fig.savefig(new_pic, dpi=300)



    def _write_statics_out(self, dir_path, length_list):
        new_file = os.path.join(dir_path, 'statistics.txt')
        pretty_out(
            'Writing ORF statistics to: {}'.format(new_file),
            linebreaks=False
        )
        out_file = open(new_file, 'w')

        out_file.write('FASTA length: {}\n\n'.format(len(self.fasta.forward_seq)))
        out_file.write('Total Number of ORFs: {}\n' \
            .format(len(length_list['AA'])))

        for stand, strand_data in self.statistics.items():
            for frame, frame_count in strand_data['orfs'].items():
                out_file.write('({}{})-strand ORFs: {}\n' \
                    .format(stand, frame, frame_count))

        if self.ORFinORF:
            out_file.write('\n(+)-strand ORFs inside ORFs: {}\n' \
                .format(self.statistics['+']['OiO']))
            out_file.write('(-)-strand ORFs inside ORFs: {}\n' \
                .format(self.statistics['-']['OiO']))

        for len_type, type_data in length_list.items():         
            out_file.write('\nShortest ORF in {}: {}\n' \
                .format(len_type, min(type_data)))
            out_file.write('Longest ORF in {}: {}\n' \
                .format(len_type, max(type_data)))
            out_file.write('Median ORF in {}: {}\n' \
                .format(len_type, np.median(type_data)))
            out_file.write('Mean ORF length in {}: {}\n' \
                .format(len_type, np.mean(type_data)))
            out_file.write('STD ORF length in {}: {}\n' \
                .format(len_type, np.std(type_data)))

        out_file.close()



################################################################################
############################## WEB BLAST QUERY #################################
################################################################################

class Web_BLAST:
    def __init__(self, program, db, queryNo, resultNo):
        self.url = BLAST_BASE_URL
        self.db = db
        self.resultNo = resultNo
        self.program = program
        self.queryNo = queryNo


    def query_orfs(self, results_dir_path, orfs_obj, rIDs=None):
        # Sort ORFS by sequence length
        orfs = sorted(orfs_obj.get_orfs_list(), key=lambda orf: len(orf.seq))

        # Redirect sterr to logging file
        pretty_out('Redirecting stderr to logfile!', start='\n', end='\n')
        stderr_logger = logging.getLogger('STDERR')
        sl = StreamToLogger(stderr_logger, logging.ERROR)
        # sys.stderr = sl

        pretty_out(
            'BLASTing {} sequences. This might take a while...' \
                .format(len(orfs))
            )
        pretty_out('\tBLAST program: {}'.format(self.program))
        pretty_out('\tBLAST DB: {}'.format(self.db))
        if len(orfs) < 1:
            pretty_out('BLASTing {} sequences... done'.format(len(orfs)))
            return

        # Create logging file for query
        # Logging file
        log_file = os.path.join(results_dir_path, 'webBLAST.log')
        logging.basicConfig(
            filename=log_file,
            format='%(asctime)s - %(levelname)s - %(message)s',
            level=logging.DEBUG
        )

        # Log query data
        logging.info('FASTA file: {}'.format(orfs_obj.fasta.file))
        logging.info('Min ORF length: {}'.format(orfs_obj.min_length))
        logging.info('Start Codons: {}'.format(','.join(orfs_obj.start_codons)))
        logging.info('Stop Codons: {}'.format(','.join(orfs_obj.stop_codons)))
        logging.info('BLAST program: {}'.format(self.program))
        logging.info('BLAST DB: {}'.format(self.db))

        if not rIDs:
            pretty_out('\tSending queries...')
            logging.info('Posting query to web BLAST')
            logging.info('\tQuery Amount: {}'.format(len(orfs)))
            query_ids = self._send_query(orfs)
            pretty_out(
                '\tIDs: {}'.format(','.join(query_ids.keys())), linebreaks=False
            )
            logging.info('RIDs: {}'.format(','.join(query_ids.keys())))
            logging.info('Posting successful')
            pretty_out('\tSending queries... done')
        else:           
            rIDs_str_raw = '{},' * len(rIDs)
            rIDs_str = rIDs_str_raw[:-1].format(*rIDs)
            pretty_out('Using old query ids: {}'.format(rIDs_str))
            logging.info('Using old query ids: {}'.format(rIDs_str))
            query_ids = {i: None for i in rIDs}

        pretty_out('\tChecking for query results...')
        if rIDs:
            pretty_out(
                'If status for ALL queries is : UNKNOWN\n \
                your queries are deprecated and deleted from the BLAST server.\n \
                Delete this run with CTRL+C and start a new one!',
                start='\n',
                end='\n',
                linebreaks=False
            )
        logging.info('Getting query results from web BLAST')   
        results = self._get_results(query_ids)
        logging.info('Getting successful') 
        pretty_out('\tChecking for query results... done')

        pretty_out('\tStoring results...')
        blast_results_list = {}
        for result_RID, result_obj in results.items():
            blast_qresults = result_obj.read().strip()
            hit = False
            for line in blast_qresults.split('\n'):
                # skip empty lines
                if not line:
                    continue
                if line.startswith('Query='):
                    orf_info = re.split(',|:', line)
                    strand = orf_info[1]
                    start_pos = int(orf_info[3])
                    blast_results_list[(strand, start_pos)] = []
                elif line[0] == '>':
                    hit_details = line[1:].strip().split(' ')
                    accession = hit_details[0]
                    description = ' '.join(hit_details[1:])
                elif line.startswith(' Score'):
                    evalue = float(line.split('=')[-1])
                    single_result = [accession, description, evalue]
                    blast_results_list[(strand, start_pos)].append(single_result)
        
        for identifier, results in blast_results_list.items():
            for res in results:
                result_hit = BLAST_result(res)
                orfs_obj.all_orfs[identifier[0]][identifier[1]] \
                    .set_blast_result(result_hit)
        pretty_out('\tStoring results... done')

        pretty_out('BLASTing {} sequences... done'.format(len(orfs)), end='\n')


    def _send_query(self, orf_list):
        fasta_queries = [''] * self.queryNo

        # Split ORFs into given Query Number FASTA files
        position = 0
        while len(orf_list) > 0:
            for query_no in range(self.queryNo):
                try:
                    orf = orf_list.pop(position)
                # All ORFs added to FASTA already
                except IndexError:
                    break
                new_fasta = '>Strand:{},Start:{}\n{}\n' \
                    .format(orf.strand, orf.start, orf.seq)
                fasta_queries[query_no] += new_fasta
            if position == 0:
                position = -1
            else:
                position = 0
        # Drop empty fasta strings
        actual_queries = [i for i in fasta_queries if i]  
        act_query_no = len(actual_queries)

        pretty_out('\tSplitted into queries: {}'.format(act_query_no)) 
        logging.info('Number of queries: {}'.format(act_query_no)) 

        query_ids = {}
        query = [
            ('CMD', 'Put'),
            ('EXPECT', 10),
            ('PROGRAM', self.program),
            ('DATABASE', self.db),
            ('HITLIST_SIZE', self.resultNo),
        ]
        # Send query for each fasta file
        for idx, query_fasta in enumerate(fasta_queries):
            orfs_length = [
                len(query_fasta.split('\n')[i]) \
                    for i in range(1, len(query_fasta.split('\n')), 2)
            ]
            sum_orf = int(sum(orfs_length))

            # Loop till post successful
            while True:
                cur_query = query + [('QUERY', query_fasta[:-1])]

                query_params = urlencode(cur_query).encode('utf-8')
                request = Request(
                    self.url,
                    query_params,
                    {"User-Agent": "ORFFinderClient"}
                )
                try:
                    req_handle = urlopen(request)
                except OSError:
                    pretty_out('\tInternet connection failed.')
                    time.sleep(3)
                    continue
                else:
                    resp_text = req_handle.read().decode("utf-8")

                    resp_info = self._get_qBLAST_info(resp_text)
                    RID = re.search('RID = (.*)\n', resp_info).group(1)
                    query_ids[RID] = None
                    pretty_out('\t{:2}/{:2} send\t({} ORFs; total bp\'s: {})' \
                        .format(idx+1, act_query_no, len(orfs_length), sum_orf)
                    )
                    time.sleep(3)
                    break

        return query_ids


    def _get_results(self, query_ids):
        query = [
            ('CMD', 'GET'),
            ('FORMAT_TYPE', 'Text'),
        ]
        # Variables to store time of server contact
        last_server_contact = datetime.now()
        last_id_query = {
            i: datetime.now() - timedelta(minutes=1) for i in query_ids
        }
        # Run till there is a result for every query
        while True:
            waiting_queries = 0
            # Loop over single queries
            for query_RID, query_result in query_ids.items():
                # Check if query already successful
                if query_result:
                    continue

                cur_query = query + [('RID', query_RID)]
                query_params = urlencode(cur_query).encode('utf-8')       
                request = Request(
                    self.url,
                    query_params,
                    {"User-Agent": "ORFFinderClient"}
                )
                query_time = datetime.now()
                query_time_diff = query_time - last_id_query[query_RID]
                server_time_diff = query_time - last_server_contact
                # Sleep if last server contact happened less than 3 secs ago
                if server_time_diff.total_seconds() < 3:
                    time.sleep(3 - server_time_diff.total_seconds())
                # Sleep if last server contact for this special ID
                # happened less than 60 secs ago
                if query_time_diff.total_seconds() < 60:
                    time.sleep(60 - query_time_diff.total_seconds())
                try:
                    handle = urlopen(request)
                    last_server_contact = datetime.now()
                # If network connection fails, try again in 60 secs    
                except OSError:
                    pretty_out('\tInternet connection failed.')
                    last_server_contact = datetime.now()
                    continue

                last_id_query[query_RID] = query_time
                resp_text = handle.read().decode("utf-8")
                if resp_text == '\n\n':
                    continue
                request_info = self._get_qBLAST_info(resp_text)

                # Query still running
                try:
                    status = re.search('Status=(.*)\n', request_info).group(1)
                    pretty_out('\tRID: {}\tStatus: {}'.format(query_RID, status))
                    logging.info('\tRID: {},Status: {}'.format(query_RID, status))             
                    if status == 'READY':
                        query_ids[query_RID] = StringIO(resp_text)
                    else:
                        waiting_queries += 1
                # Query successful
                except AttributeError:
                    pretty_out('\tRID: {}\tStatus: READY'.format(query_RID))
                    logging.info('\tRID: {},Status: READY'.format(query_RID))
                    query_ids[query_RID] = StringIO(resp_text)
                # Sleep 3 secs till next check
            pretty_out('\t{dash} {left:2}/{total:2} left  {dash}' \
                .format(dash='-'*13, left=waiting_queries, total=len(query_ids))
            )
            if waiting_queries == 0:
                break

        return query_ids

    


    def _get_qBLAST_info(self, resp_text):
        resp_info_start = resp_text.find('QBlastInfoBegin')
        resp_info_stop = resp_text.find('QBlastInfoEnd', resp_info_start)
        return resp_text[resp_info_start: resp_info_stop-1]


################################################################################
############################## BLAST RESULT ####################################
################################################################################


class BLAST_result:
    def __init__(self, result):
        self.result = result
        for res_details in result:
            try:
                self.eVal = float(res_details)
            except ValueError:
                if re.match('^[a-zA-Z]{2,2}_\d*', res_details):
                    self.accession = res_details
                else:
                    self.description = res_details


    def to_output(self):
        out_str = '{}\t{}\t{}' \
            .format(self.accession, self.eVal, self.description)
        return out_str


################################################################################
########################### STANDALONE BLAST QUERY #############################
################################################################################
class Standalone_BLAST:
    def __init__(self, program_path, db_path, queryNo, resultNo):
        self.program_path = program_path
        self.program = os.path.basename(program_path)
        self.db_path = db_path           
        self.queryNo = queryNo
        self.resultNo = str(resultNo)     


    def query_orfs(self, results_dir_path, orfs_obj):
        # Sort ORFS by sequence length
        orfs = sorted(orfs_obj.get_orfs_list(), key=lambda orf: len(orf.seq))
        # Return if no ORFs found. Pretty unlikely and probably a bug
        if len(orfs) < 1:
            pretty_out('BLASTing {} sequences... done'.format(len(orfs)))
            return

        # Redirect sterr to logging file
        pretty_out('Redirecting stderr to logfile!', start='\n', end='\n')
        stderr_logger = logging.getLogger('STDERR')
        sl = StreamToLogger(stderr_logger, logging.ERROR)
       
        pretty_out('standalone BLAST program: {}'.format(self.program))
        pretty_out('BLAST DB: {}'.format(os.path.basename(self.db_path)))
        
        # Create logging file for query
        # Logging file
        log_file = os.path.join(results_dir_path, 'standaloneBLAST.log')
        logging.basicConfig(
            filename=log_file,
            format='%(asctime)s - %(levelname)s - %(message)s',
            level=logging.DEBUG
        )

        # Log query data
        logging.info('FASTA file: {}'.format(orfs_obj.fasta.file))
        logging.info('Min ORF length: {}'.format(orfs_obj.min_length))
        logging.info('Start Codons: {}'.format(','.join(orfs_obj.start_codons)))
        logging.info('Stop Codons: {}'.format(','.join(orfs_obj.stop_codons)))     
        logging.info(
            'standalone BLAST path: {}' \
                .format(os.path.dirname(self.program_path))
        )  
        logging.info(
            'standalone BLAST program: {}' \
                .format(self.program)
        )
        logging.info('BLAST DB path: {}'.format(os.path.dirname(self.db_path)))
        logging.info('BLAST DB name: {}'.format(os.path.basename(self.db_path)))

        fasta_queries = [''] * self.queryNo

        # Split ORFs into given Query Number FASTA files
        position = 0
        while len(orfs) > 0:
            for query_no in range(self.queryNo):
                try:
                    orf = orfs.pop(position)
                # All ORFs added to FASTA already
                except IndexError:
                    break
                new_fasta = '>Strand:{},Start:{}\n{}\n' \
                    .format(orf.strand, orf.start, orf.seq)
                fasta_queries[query_no] += new_fasta
            if position == 0:
                position = -1
            else:
                position = 0
        # Drop empty fasta strings
        actual_queries = [i for i in fasta_queries if i]  
        act_query_no = len(actual_queries)

        pretty_out('Query Number: {}'.format(act_query_no)) 
        logging.info('Query Number: {}'.format(act_query_no)) 

        # q = queue.Queue()
        # threads = []
        # for idx, query in enumerate(actual_queries):          
        #     pretty_out('\tAdding query as thread: {:5}/{}' \
        #         .format(idx+1, len(actual_queries))
        #     )           
        #     t = threading.Thread(
        #         target=self._run_standlone_blast(query, orfs_obj)
        #     )
        #     t.start()
        #     threads.append(t)

        ## Add all threads  to queue
        # for tread in threads:
        #     q.put(tread)
        # # Block until all tasks are done
        # q.join()

        pretty_out('Running BLAST. Time for number crunching...') 
        running_queries = []
        temp_files = []
        for query_str in actual_queries:
            temp = tempfile.NamedTemporaryFile() 
            temp.write(query_str.encode('utf-8')  )
            temp.flush()
            temp_files.append(temp)
            # Run standalone BLAST query as subprocess
            single_query = Popen(
                [self.program_path,
                '-query', temp.name,
                '-db', self.db_path,
                '-task', self.program,
                '-dust', 'no',
                '-outfmt', "6 qseqid sacc evalue stitle",
                '-max_target_seqs', self.resultNo],
                stdout=subprocess.PIPE
            )
            running_queries.append(single_query)
        # Sleep 5 secs to allow writing of temp files to filesystem
        # BLAST takes longer anyhow...
        time.sleep(5)
        for tfile in temp_files:
            tfile.close()

        while running_queries:
            for blast_proc in running_queries:
                poll_code = blast_proc.poll()
                if poll_code is not None:
                    running_queries.remove(blast_proc)
                    results_list = str(blast_proc.stdout.read(), 'utf-8') \
                        .strip().split('\n')
                    # Add results to orf object
                    for result in results_list:
                        orf_info = result.split('\t')[0]
                        strand = orf_info.split(',')[0][-1]
                        start_pos = int(orf_info.split(':')[-1])
                        result_list = result.split('\t')[1:]
                        orfs_obj.all_orfs[strand][start_pos] \
                            .set_blast_result(BLAST_result(result_list))
                    # process results
                else:
                    time.sleep(30.0 / act_query_no)
            pretty_out(
                '\t{} / {} queries finished...' \
                    .format(act_query_no-len(running_queries), act_query_no)
            ) 

       





    def _run_standlone_blast(self, fasta_str, orfs_obj):
        # Write orf in temporary fasta file
        temp = tempfile.NamedTemporaryFile() 
        temp.write(fasta_str.encode('utf-8')  )
        temp.flush()
        # Run standalone BLAST query as subprocess
        single_query = subprocess.run(
            [self.program_path,
            '-query', temp.name,
            '-db', self.db_path,
            '-task', self.program,
            '-dust', 'no',
            '-outfmt', "6 qseqid sacc evalue stitle",
            '-max_target_seqs', '1'],
            stdout=subprocess.PIPE
        )
        temp.close()
        # Get results
        results_list = str(single_query.stdout, 'utf-8').strip().split('\n')
        # Add results to orf object
        for result in results_list:
            orf_info = result.split('\t')[0]
            strand = orf_info.split(',')[0][-1]
            start_pos = int(orf_info.split(':')[-1])
            result_str = '\t'.join(result.split('\t')[1:])
            orfs_obj.all_orfs[strand][start_pos] \
                .set_blast_result(BLAST_result([result_str]))





################################################################################
############################ HELPER FUNCTIONS ##################################
################################################################################

def pretty_out(msg, start='', end='', linebreaks=True):
    if linebreaks:
        msg_formatted = '\n\t\t\t'.join(
            msg[i:i+60] for i in range(0, len(msg), 60)
        )
    else:
        msg_formatted = msg
    print('{s}{t}\t{m}{e}'.format(
        t=datetime.strftime(datetime.now(), '%Y.%m.%d %H.%M.%S'),
        m=msg_formatted,
        s=start,
        e=end,
        )
    )

################################################################################
################################## MAIN ########################################
################################################################################

if __name__ == '__main__':
    pretty_out('ORF_finder start')
    args = parser.parse_args()

    # Input file is a config file
    input_file_type = args.input_file.split('.')[-1]

    # Input is config
    if input_file_type == 'cfg':
        pretty_out('Reading input from config file...')
        # Parse config
        config = configparser.ConfigParser()
        config.read(args.input_file) 
        # Parse arguments
        input_file = config['INPUT']['InputPath']
        output_dir = config.get(
            'OUTPUT',
            'OutputPath',
            fallback=os.getcwd(),
        )

        start_codons = config.get(
            'ORF_Finder',
            'StartCodon',
            fallback='AUG'
        ).strip().split(',')
        stop_codons = config.get(
            'ORF_Finder',
            'StopCodons',
            fallback='UAA,UAG,UGA'
        ).strip().split(',')
        min_length = config. getint(
            'ORF_Finder',
            'MinLength',
            fallback=200
        )
        OrfInOrf = config.getboolean(
            'ORF_Finder',
            'ORFinORF',
            fallback=False
        )

        oldWebBlastIDs = False # Not implemented for config yet
        blastQueryNo = config. getint(
            'BLAST',
            'NumberOfQueries',
            fallback=20
        )
        blastResultNo = config. getint(
            'BLAST',
            'ResultsPerORF',
            fallback=1
        )
        blastProgram = config.get(
            'BLAST',
            'BlastProgram',
            fallback='blastn'
        )
        blastDB = config.get(
            'BLAST',
            'DBName',
            fallback='nt'
        )

        stAlBlastResults = config.getboolean(
            'BLAST',
            'RunStandaloneBlast',
            fallback=False
        )
        if stAlBlastResults:
            stAlBlastPath = os.path.join(
                config['BLAST']['BlastPath'],
                'bin',
                config['BLAST']['BlastProgram']
            )
            stAlDBPath = os.path.join(
                config['BLAST']['DBPath'],
                config['BLAST']['DBName']
            )
        webBlastResults = config.getboolean(
            'BLAST',
            'RunWebBlast',
            fallback=False
        )
        if webBlastResults and stAlBlastResults:
            stAlBlastResults = False
            pretty_out(
                'Flags for Standalone and web BLAST set to True. \
                Web version is used...'
            )      
    # Input is fasta
    elif input_file_type in  ['fa', 'mpfa', 'fna', 'fsa', 'fasta']:
        pretty_out('Reading input from command line...')
        # Parse arguments
        input_file = args.input_file
        output_dir = args.output

        start_codons = args.start.strip().split(',')
        stop_codons = args.stop.strip().split(',')
        min_length = args.length
        OrfInOrf = args.ORFinORF

        blastQueryNo = args.blastQueryNumber
        blastResultNo = args.blastResultNumber
        blastProgram = args.blastProgram
        blastDB = args.blastDB

        stAlBlastResults = False # Not implemented for command line yet
        webBlastResults = args.blast      
    # Input is unknown
    elif input_file_type == 'log':
        pretty_out('Reading input from log file...')
        log_file = open(args.input_file, 'r')
        for line in log_file:
            content = line.split(' - ')
            in_data = content[2].split(':')
            option = in_data[0].strip()
            try:
                option_value = in_data[1].strip()
            except IndexError:
                continue
            print(option)
            if option == 'FASTA file':
                input_file = option_value
            elif option == 'Min ORF length':
                min_length = int(option_value)
            elif option == 'Start Codons':
                start_codons = option_value.split(',')
            elif option == 'Stop Codons':
                stop_codons = option_value.split(',')
            elif option == 'BLAST program':
                blastProgram = option_value
            elif option == 'BLAST DB':
                blastDB = option_value
            elif option == 'Number of queries':
                blastQueryNo = int(option_value)
            elif option == 'RIDs' :
                oldWebBlastIDs = option_value.split(',')
            elif option == 'Using old query ids':                
                # Deprecated: Only needed for "old" configs.
                # Can be removed in future
                option_value = option_value.lstrip('[').rstrip(']')
                oldWebBlastIDs = [
                    i.strip()[1:-1] for i in option_value.split(',')
                ]
            elif option == 'RID':
                break
        log_file.close()

        # If webBLAST.log file is already a run started with RIDs,
        # no query number is written to config
        if not 'blastQueryNo' in locals():
            blastQueryNo = None

        output_dir = args.output
        OrfInOrf = args.ORFinORF
        blastResultNo = args.blastResultNumber
        stAlBlastResults = False
        webBlastResults = True
    else:
        raise IOError('Unknown input type. Use: .cfg or .fa/.fna/.fasta/.fsa')

    pretty_out('Reading file: {}'.format(input_file))
    fasta = FASTA(input_file)
    pretty_out('\tbp: {}'.format(len(fasta.forward_seq)))
    pretty_out('Reading file... done')

    pretty_out('Using start codons: {}'.format(', '.join(start_codons)))
    pretty_out('Using stop codons: {}'.format(', '.join(stop_codons)))
    pretty_out('Report ORFs with min length: {}'.format(min_length), end='\n')

    pretty_out('Seaching for ORF\'s...')
    orfs = ORFFinder(
        fasta,
        min_length,
        start_codons,
        stop_codons,
        args.ORFinORF
    )
    pretty_out('Seaching for ORF\'s... done', end='\n')

    # Create directory to store results in
    timestamp = datetime.strftime(datetime.now(), '%Y.%m.%dT%H.%M.%S')
    if webBlastResults or oldWebBlastIDs or stAlBlastResults:
        results_dir_name = '{}_ORF_Finder_results_l{}_BLASTed_{}' \
            .format(timestamp, min_length, blastDB) 
    else:
        results_dir_name = '{}_ORF_Finder_results_l{}' \
            .format(timestamp, min_length) 
    results_dir_path = os.path.join(output_dir, results_dir_name)
    os.mkdir(results_dir_path)
    pretty_out(
        'Results Directory: {}'.format(results_dir_path),
        linebreaks=False
    )

    # Run BLAST if argument given
    # Run standalone BLAST
    if stAlBlastResults:
        standaloneBlast = Standalone_BLAST(
            stAlBlastPath,
            stAlDBPath,
            blastQueryNo,
            blastResultNo
        )
        standaloneBlast.query_orfs(results_dir_path, orfs)
        orfs.write_output(results_dir_path, blastProgram, blastDB)
    # Run web BLAST
    elif webBlastResults or oldWebBlastIDs:
        webBlast = Web_BLAST(
            blastProgram,
            blastDB,
            blastQueryNo,
            blastResultNo
        )        
        if oldWebBlastIDs:
            webBlast.query_orfs(
                results_dir_path,
                orfs,
                rIDs=oldWebBlastIDs
            )                       
        else:
            webBlast.query_orfs(results_dir_path, orfs)
        orfs.write_output(results_dir_path, blastProgram, blastDB)
    # Don not run any BLAST
    else:
        orfs.write_output(results_dir_path, 'not used', 'not used')

    pretty_out('ORF_finder exit')