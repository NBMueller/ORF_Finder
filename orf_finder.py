#!/usr/bin/env python3

""" This program takes a FASTA file, reads the sequence and searches for ORF\'s.
    In Addition some functions are run on the found ORF\'s.
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
from datetime import datetime
# External package libraries
import numpy as np
import matplotlib.pyplot as plt
# from Bio.Blast import NCBIWWW
# from Bio.Blast import NCBIXML
from Bio import SearchIO

__version__ = "1.0.0"
__author__ = "Nico Borgsmueller"
__contact__ = "nicoborgsmueller@web.de"
__credits__ = ["Andreas Andrusch"]
__copyright__ = "Copyright 2017, Nico Borgsmueller"
__license__ = "GNU GPLv3"
__status__ = "Development"
__date__ = "2017.03.24"

################################################################################
################################ ARGPARSER #####################################
################################################################################

parser = argparse.ArgumentParser(description='Find orf\'s in a FASTA file.')

parser.add_argument(
    'fasta_file',
    type=str,
    help="Path to FASTA file."
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
    '-bRID', '--blastRID',   
    type=str,
    help='Use old web Blast query RIDs ("," separated).',
)
parser.add_argument(
    '-bqn', '--blastQueryNumber',   
    type=int,
    default=20,
    help='Number of queries used for BLASTing results (default=20).')

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

        self.blast = None
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
        self.blast = result


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
            return '-\t-\t-'
        else:
            return self.blast.to_output()


    def get_out_str(self):
        out_str= '{s}\t{f}\t{st}\t{sp}\t{l}\t{seq}\t{prot}\t{encl}\t{embed}\t{b}\n' \
            .format(
                s=self.strand,
                f=self.frame,
                st=self.start,
                sp=self.stop,
                l=self.length,
                seq=self.seq,
                prot=self.AA_seq,
                encl=self._get_outsideORF_str(),
                embed=self._get_insideORF_str(),
                b=self._get_blast_result()
            )
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

        pretty_out('Start scanning forward sequence...', start='\n')
        forward_orfs = self._get_orfs(fasta.forward_seq, '+')
        pretty_out('...done', end='\n')

        pretty_out('Start scanning reverse sequence...')
        reverse_orfs = self._get_orfs(fasta.reverse_seq, '-')
        pretty_out('...done', end='\n')

        if self.ORFinORF:
            pretty_out('Scanning for ORF\'s inside other ORF\'s...')
            all_orfs_sorted = {
                '+': forward_orfs,
                '-': reverse_orfs
            }
            self.find_OiOs(all_orfs_sorted)
            pretty_out('...done', end='\n')

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
                '({}{})-strand:\t{:5} ORF\'s'.format(strand, i, len(j))
            )
            self.statistics[strand]['orfs'][i] = len(j)

        # Sort all ORF's by start position
        all_orfs = {}
        for frame_orf in orfs.values():
            all_orfs.update(frame_orf)
        sorted_orfs = sorted(all_orfs.items(), key=operator.itemgetter(0))

        return sorted_orfs



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


    def write_output(self, results_dir_path, BLASTprogram, BLASTdb):
        new_file = os.path.join(results_dir_path, 'ORFs.tsv')
        pretty_out('Writing results to: {}'.format(new_file))
        out_file = open(new_file, 'w')

        out_file.write('# FASTA: {}\n'.format(self.fasta.file))
        out_file.write('# FASTA Description: {}\n'.format(self.fasta.info[1:]))
        out_file.write('# Min ORF length: {}\n'.format(self.min_length))
        out_file.write('# Start Codons: {}\n'.format(','.join(self.start_codons)))
        out_file.write('# Stop Codons: {}\n'.format(','.join(self.stop_codons)))
        out_file.write('# BLAST Program: {}\n'.format(BLASTprogram))
        out_file.write('# BLAST DB: {}\n'.format(BLASTdb))

        out_file.write(
            'Strand\tFrame\tStart\tStop\tLength\tSequence\tProtein\tEnclosing_ORF\tEmbedded_ORFs\tBLAST_accession\tBLAST_Evalue\tBLAST_description\n'
        )

        all_orfs_list = []
        for strand, strand_orfs in self.all_orfs.items():
            all_orfs_list.extend(list(strand_orfs.values()))

        length_list = {'AA': [], 'bp': []}
        for orf in all_orfs_list:        
            length_list['bp'].append(len(orf.seq))
            length_list['AA'].append(len(orf.AA_seq))
            out_file.write(orf.get_out_str())
        out_file.close()

        # Statistics out
        self._write_statics_out(results_dir_path, length_list)
        self._plot_histogram(results_dir_path, length_list)

        pretty_out('...done')


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
        pretty_out('Writing ORF statistics to: {}'.format(new_file))
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


 # # Stdout Functions for Development

 #    def get_all_inside_OiO(self, strand='both'):
 #        return self._get_all_OiO(strand, 'inside')


 #    def get_all_outside_OiO(self, strand='both'):
 #        return self._get_all_OiO(strand, 'outside')


 #    def _get_all_OiO(self, strand, OiO_type):
 #        OiO_list = []
 #        if strand =='both':
 #            OiO_list.extend(self._get_OiO('forward', OiO_type))
 #            OiO_list.extend(self._get_OiO('reverse', OiO_type))
 #        else:
 #            OiO_list.extend(self._get_OiO(strand, OiO_type))
 #        return OiO_list


 #    def _get_OiO(self, strand, OiO_type):
 #        true_orfs = []
 #        for orf in self.all_orfs[strand]:
 #            if OiO_type == 'inside':
 #                if orf[1].outsideORF:
 #                    true_orfs.append(orf[1])
 #            elif OiO_type == 'outside':
 #                if orf[1].insideORF:
 #                    true_orfs.append(orf[1])
 #        return true_orfs



################################################################################
############################## WEB BLAST QUERY #################################
################################################################################

class Web_BLAST:
    def __init__(self, program, db, queryNo, hitlist_size=1):
        self.url = BLAST_BASE_URL
        self.db = db
        self.hitlist_size = hitlist_size
        self.program = program
        self.queryNo = queryNo


    def query_orfs(self, results_dir_path, orfs_obj, rIDs=None):
        orfs = []
        for strand, strand_data in orfs_obj.all_orfs.items():
            orfs.extend(strand_data.values())

        pretty_out(
            'BLASTing {} sequences. This might take a while...' \
                .format(len(orfs))
            )
        pretty_out('\tBLAST program: {}'.format(self.program))
        pretty_out('\tBLAST DB: {}'.format(self.db))
        if len(orfs) < 1:
            pretty_out('...done')
            return

        # Create logging file for query
        # Logging file
        log_file = os.path.join(results_dir_path, 'webBLAST.log')
        logging.basicConfig(
            filename=log_file,
            format='%(asctime)s - %(levelname)s - %(message)s',
            level=logging.DEBUG
        )
        # Redirect sterr to logging file
        pretty_out('Redirecting stderr to logfile!', start='\n', end='\n')
        stderr_logger = logging.getLogger('STDERR')
        sl = StreamToLogger(stderr_logger, logging.ERROR)
        # sys.stderr = sl

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
            pretty_out('\tIDs: {}'.format(','.join(query_ids.keys())))
            logging.info('RIDs: {}'.format(','.join(query_ids.keys())))
            logging.info('Posting successful')
            pretty_out('\t...done')
        else:
            pretty_out('Using old query ids: {}'.format(','.join(rIDs)))
            logging.info('Using old query ids: {}'.format(rIDs))
            query_ids = {i: None for i in rIDs}

        pretty_out('\tChecking for query results...')
        logging.info('Getting query results from web BLAST')   
        results = self._get_results(query_ids)
        logging.info('Getting successful') 
        pretty_out('\t...done')

        pretty_out('\tStoring results...')
        blast_results_list = []
        for result_RID, result_obj in results.items():
            blast_qresults = SearchIO.parse(result_obj, 'blast-xml')
            blast_results_list.extend([i for i in blast_qresults])

        success_flag = False
        for result in blast_results_list:
            orf_info = result.id.split(',')
            strand = orf_info[0].split(':')[1]
            start_pos = int(orf_info[1].split(':')[1])
            if not success_flag and len(result) > 0:
                success_flag = True
            orfs_obj.all_orfs[strand][start_pos].set_blast_result(BLAST_result(result))
        pretty_out('\t...done')

        if success_flag:
            pretty_out('BLASTing successful!', start='\n', end='\n') 
        else:
            pretty_out('BLASTing failed! Probably too many queries...') 

        pretty_out('...done')


    def _send_query(self, orf_list):
        fasta_queries = [''] * self.queryNo

        # Split ORFs into given Query Number FASTA files
        for orf_idx, orf in enumerate(orf_list):
            query_no = orf_idx % self.queryNo
            new_fasta = '>Strand:{},Start:{}\n{}\n' \
                .format(orf.strand, orf.start, orf.seq)
            fasta_queries[query_no] += new_fasta
              
        # Drop empty fasta strings
        actual_queries = [i for i in fasta_queries if i]  
        actual_query_no = len(actual_queries)

        pretty_out('\tSplitted into queries: {}'.format(actual_query_no)) 
        logging.info('Number of queries: {}'.format(actual_query_no)) 

        query_ids = {}
        query = [
            ('CMD', 'Put'),
            ('EXPECT', 10),
            ('PROGRAM', self.program),
            ('DATABASE', self.db),
            ('HITLIST_SIZE', self.hitlist_size),
        ]
        # Send query for each fasta file
        for idx, query_fasta in enumerate(fasta_queries):
            # Loop till post successful
            while True:
                cur_query = query + [('QUERY', query_fasta[:-1])]
                # import pdb; pdb.set_trace()
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
                    pretty_out('\t{}/{} send'.format(idx+1, actual_query_no))
                    time.sleep(3)
                    break

        return query_ids


    def _get_results(self, query_ids):
        query = [
            ('CMD', 'GET'),
            ('FORMAT_TYPE', 'XML'),
        ]
        # Variables to store time of server contact
        last_server_contact = datetime.now()
        last_id_query = {i: datetime.now() for i in query_ids}
        # RUn till there is a result for every query
        while True:
            all_query_successfull = True 
            # Loop over single queries
            for query_RID, query_result in query_ids.items():
                # Check if query already successful
                if query_result:
                    continue
                else:
                    all_query_successfull = False

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
                # Query successful
                except AttributeError:
                    pretty_out('\tRID: {}\tStatus: READY'.format(query_RID))
                    logging.info('\tRID: {},Status: READY'.format(query_RID))
                    query_ids[query_RID] = StringIO(resp_text)
                # Sleep 3 secs till next check
            pretty_out('\t' + '-'*40)
            if all_query_successfull:
                break

        return query_ids

    


    def _get_qBLAST_info(self, resp_text):
        resp_info_start = resp_text.find('QBlastInfoBegin')
        resp_info_stop = resp_text.find('QBlastInfoEnd', resp_info_start)
        return resp_text[resp_info_start: resp_info_stop-1]


class BLAST_result:
    def __init__(self, result_obj):
        self.result = result_obj
        if result_obj.hits:
            self.description = result_obj[0].description
            self.accession = result_obj[0].accession
            self.HSPs = result_obj[0].hsps[0].evalue
        else:
            self.description = 'no_hit'
            self.accession = 'no_hit'
            self.HSPs = 'no_hit'


    def to_output(self):
        out_str = '{}\t{}\t{}' \
            .format(self.accession, self.HSPs, self.description)
        return out_str


################################################################################
############################ HELPER FUNCTIONS ##################################
################################################################################

def pretty_out(msg, start='', end=''):
    msg_formatted = '\n\t\t\t'.join(msg[i:i+60] for i in range(0, len(msg), 60))
    print('{s}{t}\t{m}{e}'.format(
        t=datetime.strftime(datetime.now(), '%Y.%m.%dT%H.%M.%S'),
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

    pretty_out('Reading file: {}'.format(args.fasta_file))
    fasta = FASTA(args.fasta_file)
    pretty_out('\tbp: {}'.format(len(fasta.forward_seq)))
    pretty_out('...done')

    start_codons = args.start.strip().split(',')
    stop_codons = args.stop.strip().split(',')

    pretty_out('Using start codons: {}'.format(', '.join(start_codons)))
    pretty_out('Using stop codons: {}'.format(', '.join(stop_codons)))
    pretty_out('Report ORFs with min length: {}'.format(args.length))

    pretty_out('Seaching for ORF\'s...')
    orfs = ORFFinder(
        fasta,
        args.length,
        start_codons,
        stop_codons,
        args.ORFinORF
    )
    pretty_out('...done')

    # Create directory to store results in
    timestamp = datetime.strftime(datetime.now(), '%Y.%m.%dT%H.%M.%S')
    old_path = os.path.split(fasta.file)
    if args.blast or args.blastRID:
        results_dir_name = '{}_ORF_Finder_results_l{}_BLASTed_{}' \
            .format(timestamp, args.length, args.blastDB) 
    else:
        results_dir_name = '{}_ORF_Finder_results_l{}' \
            .format(timestamp, args.length) 
    results_dir_path = os.path.join(old_path[0], results_dir_name)
    os.mkdir(results_dir_path)
    pretty_out('...done')
    # Run BLAST if argument given
    if args.blast or args.blastRID:
        webBlast = Web_BLAST(
            args.blastProgram,
            args.blastDB,
            args.blastQueryNumber
        )
        if args.blast:
            webBlast.query_orfs(results_dir_path, orfs)
        else:
            webBlast.query_orfs(
                results_dir_path,
                orfs,
                rIDs=args.blastRID.strip().split(',')
            )
        orfs.write_output(results_dir_path, args.blastProgram, args.blastDB)
    else:
        orfs.write_output(results_dir_path, 'not used', 'not used')

    pretty_out('ORF_finder exit')
