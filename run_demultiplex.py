#!/usr/bin/env python3

from futures_demulti import *
from demultiplex_input_process import ProcessDemultiplexInput

test = ProcessDemultiplexInput('1_test.^.qseq.txt', '2_test.^.qseq.txt', directory='tests/test_qseq/',
sample_key='tests/test_sample_files/single_index_test.txt', barcode_1='tests/test_sample_files/N700_nextera_barcodes.txt',
file_label='rb')
test.run(b1_reverse=False, b2_reverse=False)
output_run = output_objects(sample_list=test.sample_list, output_directory='tests/test_output/', read_count=test.read_count)

output_dictionary=output_run
barcode_1=test.barcode_1
barcode_2=test.barcode_2
sample_key=test.sample_key
file_label=test.file_label
directory=test.directory
gnu_zipped=False

iterate_through_qseq(file_list=test.file_list)

import DemultiplexClass
import time
import argparse


def launch_demultiplex(*args, directory='path', sample_key='path', mismatch=3, file_label='', barcode_1=None,
                       barcode_2=None, output_directory=None, gnu_zipped=False):
    """Simple function to initialize DemultiplexClass"""
    start_time = time.time()
    demultiplex = DemultiplexClass.Demuliplex(*args, directory=directory, barcode_1=barcode_1, barcode_2=barcode_2,
                                              sample_key=sample_key, file_label=file_label, mismatch=mismatch,
                                              gnu_zipped=gnu_zipped)

    demultiplex.run(output_directory=output_directory)
    end_time = time.time()
    print('Total reads:' + str(demultiplex.reads))
    print('Reads passing filter:' + str(demultiplex.reads_pass_filter))
    print('Indexed reads:' + str(demultiplex.indexed_reads))
    print('Unmatched reads:' + str(demultiplex.unmatched_read))
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
parser.add_argument('-L', type=str, help='string of r and b character to designate input files as '
                                         'barcode or read files, should be the same order as input'
                                         'file')
parser.add_argument('-M', type=int, default=2, help='number of barcode mismatches to consider, default = 2')
parser.add_argument('-O', type=str, help='path to output directory')
parser.add_argument('-Z', action="store_true", default=False, help='if qseq files gzipped, slows processing')
parser.add_argument('-I', type=str, nargs='*', help='qseq file prefix and suffix separated'
                                                    'by ^, ie. -I s_1_^.qseq.txt '
                                                    's_2_^.qseq.txt ')
arguments = parser.parse_args()

print('Started Job')

launch_demultiplex(*arguments.I, directory=arguments.D, barcode_1=arguments.B1, barcode_2=arguments.B2,
                   sample_key=arguments.S, output_directory=arguments.O, mismatch=arguments.M, gnu_zipped=arguments.Z,
                   file_label=arguments.L)