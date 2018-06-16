#!/usr/bin/env python3

from ProcessInputFiles import ParseInputFiles
from DemultiplexRun import RunDemultiplex
import time
import argparse


def launch_demultiplex(*args, qseq_directory='path', sample_key='path', file_label='',
                       output_directory=None, barcode_reverse=False, hamming_distance=1, mixed_length=False):
    """Call demultiplex objects"""
    start_time = time.time()
    demultiplex_input = ParseInputFiles(*args,
                                        qseq_directory=qseq_directory,
                                        sample_file=sample_key,
                                        barcode_reverse=barcode_reverse,
                                        file_label=file_label)
    demultiplex_run = RunDemultiplex(input_object=demultiplex_input,
                                     output_directory=output_directory,
                                     ham_dist=hamming_distance,
                                     mixed_length=mixed_length)
    run_metrics = demultiplex_run.run_stats

    end_time = time.time()
    print('Total reads:%s' % str(run_metrics['reads']))
    print('Reads passing filter:%s' % str(run_metrics['reads_passing_filter']))
    print('Indexed reads:%s' % str(run_metrics['indexed_reads']))
    print('Unmatched reads:%s' % str(run_metrics['unmatched_reads']))
    print('Total time:%s minutes' % str(round((end_time - start_time) / 60.0, 2)))


parser = argparse.ArgumentParser(description='Demultiplexing script. Script demultiplexes Illumnina qseq lane files '
                                             'outputing sample fastq files. Works with .gz and uncompressed qseq files.'
                                             ' Options for single and dual indexes'
                                             '\n\n Usage; python3 Demultiplex.py -D directory -S sample_key'
                                             ' -BR False -L file_labels -H hamming_distance_threshold -O '
                                             'output_directory -I input_file_1 input_file_2 ...')

parser.add_argument('-D', type=str, help='/path/ to qseq directory')
parser.add_argument('-S', type=str, help='/path/sample_file.txt file should be formatted as \''
                                         'barcode tab sample_name\' for single index and '
                                         '\'barcode tab barcode tab sample_name\' '
                                         'for dual indexes ')
parser.add_argument('-BR', action="store_true", default=False, help='Consider Barcodes Reverse Complements')
parser.add_argument('-L', type=str, help='string of r and b character to designate input files as '
                                         'barcode or read files, should be the same order as input'
                                         'file')
parser.add_argument('-O', type=str, help='Path to Output Directory')
parser.add_argument('-I', type=str, nargs='*', help='qseq file prefix and suffix separated'
                                                    'by ^, ie. -I s_1_^.qseq.txt '
                                                    's_2_^.qseq.txt ')
parser.add_argument('-H', type=int, default=0,
                    help='Minimum hamming distance threshold for a sequencing barcode to be considered, default=0')
parser.add_argument('-M', type=int, help='If the reference barcodes for an index contain barcodes of different length,'
                    'hamming distance includes the difference in reference sequencing barcode length')
arguments = parser.parse_args()

if __name__ == "__main__":
    try:
        print('Working')
        launch_demultiplex(*arguments.I,
                           qseq_directory=arguments.D,
                           sample_key=arguments.S,
                           output_directory=arguments.O,
                           file_label=arguments.L,
                           barcode_reverse=arguments.BR,
                           hamming_distance=arguments.H,
                           mixed_length=arguments.M)
    except TypeError:
        print('python3 Demultiplex.py --help, for usage')
