#!/usr/bin/env python3

from futures_demulti import *
from demultiplex_input_process import ProcessDemultiplexInput
import time
import argparse


def launch_demultiplex(*args, directory='path', sample_key='path', file_label='', barcode_1=None,
                       barcode_2=None, output_directory=None, b1_reverse=False, b2_reverse=False, gnu_zipped=False,
                       workers=2):
    """Simple function to initialize DemultiplexClass"""
    start_time = time.time()
    demultiplex_input = ProcessDemultiplexInput(*args, directory=directory,
                                                sample_key=sample_key,
                                                barcode_1=barcode_1, barcode_2=barcode_2, file_label=file_label
                                                )
    demultiplex_input.run(b1_reverse=b1_reverse, b2_reverse=b2_reverse)
    output_run = output_objects(sample_list=demultiplex_input.sample_list, output_directory=output_directory,
                                read_count=demultiplex_input.read_count)
    run_metrics = iterate_through_qseq(demultiplex_instance=demultiplex_input,
                                       output_dictionary=output_run, gnu_zipped=gnu_zipped, workers=workers)

    end_time = time.time()
    print('Total reads:' + str(sum(run_metrics.reads)))
    print('Reads passing filter:' + str(sum(run_metrics.reads_pass_filter)))
    print('Indexed reads:' + str(sum(run_metrics.indexed_reads)))
    print('Unmatched reads:' + str(sum(run_metrics.unmatched_reads)))
    print('Total time:' + str(round((end_time - start_time) / 60.0, 2)) + ' minutes')


parser = argparse.ArgumentParser(description='Demultiplexing script. Script demultiplexes Illumnina qseq lane files '
                                             'outputing sample fastq files. Works with .gz and uncompressed qseq files.'
                                             ' Options for single and dual indexes'
                                             '\n\n Usage; demultiplex -D directory -S sample_key'
                                             ' -B1 barcode_1 -B2 barcode_2 -L file_labels -M mismatch_number -O '
                                             'output_directory -I input_file_1 input_file_2 ...')

parser.add_argument('-D', type=str, help='/path/ to qseq directory')
parser.add_argument('-S', type=str, help='/path/sample_file.txt file should be formatted as \''
                                         'barcode tab sample_name\' for single index and '
                                         '\'barcode tab barcode tab sample_name\' '
                                         'for dual indexes ')
parser.add_argument('-B1', type=str, help='/path/barcode_1_file, barcode \t index key')
parser.add_argument('-B2', type=str, default=None, help='/path/barcode_2_file, barcode \t index key')
parser.add_argument('-B1R', action="store_true", default=False, help='Consider Barcode1 Reverse Complement')
parser.add_argument('-B2R', action="store_true", default=False, help='Consider Barcode2 Reverse Complement')
parser.add_argument('-L', type=str, help='string of r and b character to designate input files as '
                                         'barcode or read files, should be the same order as input'
                                         'file')
parser.add_argument('-W', type=int, default=2, help='Number of cores available, default = 2')

parser.add_argument('-O', type=str, help='path to output directory')
parser.add_argument('-Z', action="store_true", default=False, help='if qseq files gzipped, slows processing')
parser.add_argument('-I', type=str, nargs='*', help='qseq file prefix and suffix separated'
                                                    'by ^, ie. -I s_1_^.qseq.txt '
                                                    's_2_^.qseq.txt ')
arguments = parser.parse_args()


try:
    print('Working')
    launch_demultiplex(*arguments.I, directory=arguments.D, barcode_1=arguments.B1, barcode_2=arguments.B2,
                       sample_key=arguments.S, output_directory=arguments.O, workers=arguments.W, gnu_zipped=arguments.Z,
                       file_label=arguments.L, b1_reverse=arguments.B1, b2_reverse=arguments.B2)
except TypeError:
    print('python3 Demultiplex.py --help, for usage')
