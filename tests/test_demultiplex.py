#!/usr/bin/env python3

import subprocess
import os
import unittest

from ProcessInputFiles import ParseInputFiles
from DemultiplexRun import RunDemultiplex

test_directory = os.path.dirname(os.path.realpath(__file__))
root_directory = '/'.join(test_directory.split('/')[0:-1]) + '/'

test_single_index_demultiplex = ParseInputFiles('1_test.^.qseq.txt.gz', '2_test.^.qseq.txt',
                                                qseq_directory=test_directory + '/test_qseq/',
                                                sample_file=test_directory + '/test_sample_files/single_index_test.txt',
                                                barcode_reverse=False,
                                                file_label='rb')


test_dual_index_demultiplex = ParseInputFiles('1_test.^.qseq.txt.gz', '2_test.^.qseq.txt',
                                              '3_test.^.qseq.txt', '4_test.^.qseq.txt',
                                              qseq_directory=test_directory + '/test_qseq/',
                                              sample_file=test_directory + '/test_sample_files/dual_index_test.txt',
                                              barcode_reverse=True,
                                              file_label='rbbr')

test_single_run = RunDemultiplex(input_object=test_single_index_demultiplex,
                                 output_directory=test_directory + '/test_output/',
                                 ham_dist=1)
test_single_metrics = dict(test_single_run.run_stats)

test_dual_run = RunDemultiplex(input_object=test_dual_index_demultiplex,
                               output_directory=test_directory + '/test_output/',
                               ham_dist=1)
test_dual_metrics = dict(test_dual_run.run_stats)
print('fuck')

class TestDemultiplex(unittest.TestCase):
    def setUp(self):
        pass

    def test_single_filter_pass(self):
        self.assertEqual(test_single_metrics['reads_passing_filter'], 19614)

    def test_single_unmatched_reads(self):
        self.assertEqual(test_single_metrics['unmatched_reads'], 2774)

    def test_single_index_reads(self):
        self.assertEqual(test_single_metrics['indexed_reads'], 16840)

    def test_single_samples(self):
        x = test_single_index_demultiplex.file_list[0][0].replace('1_test', '').replace('.gz', '')
        y = test_single_index_demultiplex.file_list[0][1].replace('2_test', '')
        self.assertEqual(x, y)

    def test_dual_filter_pass(self):
        self.assertEqual(test_dual_metrics['reads_passing_filter'], 19264)

    def test_total_reads(self):
        self.assertEqual(test_dual_metrics['reads'], test_single_metrics['reads'])

    def test_dual_unmatched_reads(self):
        self.assertEqual(test_dual_metrics['unmatched_reads'], 9420)

    def test_dual_index_reads(self):
        self.assertEqual(test_dual_metrics['indexed_reads'], 9844)

    def test_dual_samples(self):
        x = test_dual_index_demultiplex.file_list[0][2].replace('3_test', '')
        y = test_dual_index_demultiplex.file_list[0][3].replace('4_test', '')
        self.assertEqual(x, y)

    def test_command_line_single_index(self):
        parser_open = subprocess.run(['python3', root_directory + 'Demultiplex.py', '-D',
                                      test_directory + '/test_qseq/', '-S',
                                      test_directory + '/test_sample_files/single_index_test.txt',
                                      '-L', 'rb', '-O', test_directory + '/test_output/',
                                      '-I', '1_test.^.qseq.txt.gz', '2_test.^.qseq.txt'],
                                     stdout=subprocess.PIPE)
        output = parser_open.stdout
        output_categories = (output.decode()).split('\n')
        filtered = int(output_categories[2].split(':')[1])
        indexed = int(output_categories[3].split(':')[1])
        unmatched = int(output_categories[4].split(':')[1])
        self.assertEqual(filtered, 19614)
        self.assertEqual(indexed, 16840)
        self.assertEqual(unmatched, 2774)

    def test_command_line_dual_index(self):
        parser_open = subprocess.run(['python3', root_directory + 'Demultiplex.py', '-D',
                                      test_directory + '/test_qseq/', '-S',
                                      test_directory + '/test_sample_files/dual_index_test.txt', '-L', 'rbbr', '-O',
                                      test_directory + '/test_output/', '-I', '1_test.^.qseq.txt.gz',
                                      '2_test.^.qseq.txt', '3_test.^.qseq.txt', '4_test.^.qseq.txt'],
                                     stdout=subprocess.PIPE)
        output = parser_open.stdout
        print(output)
        output_categories = (output.decode()).split('\n')
        filtered = int(output_categories[2].split(':')[1])
        indexed = int(output_categories[3].split(':')[1])
        unmatched = int(output_categories[4].split(':')[1])
        self.assertEqual(filtered, 19264)
        self.assertEqual(indexed, 9844)
        self.assertEqual(unmatched, 9420)


if __name__ == '__main__':
    unittest.main()
