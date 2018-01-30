#!/usr/bin/env python3

import subprocess
import unittest

from demultiplex_input_process import ProcessDemultiplexInput
from futures_demulti import *

test_single_index_demultiplex = ProcessDemultiplexInput('1_test.^.qseq.txt', '2_test.^.qseq.txt',
                                                        directory='tests/test_qseq/',
                                                        sample_key='tests/test_sample_files/single_index_test.txt',
                                                        barcode_1='tests/test_sample_files/'
                                                                  'N700_nextera_barcodes.txt',
                                                        file_label='rb')
test_single_index_demultiplex.run(b1_reverse=True, b2_reverse=True)
test_single_metrics = iterate_through_qseq(demultiplex_instance=test_single_index_demultiplex,
                                           output_directory='tests/test_output/', gnu_zipped=False)

test_dual_index_demultiplex = ProcessDemultiplexInput('1_test.^.qseq.txt', '2_test.^.qseq.txt',
                                                      '3_test.^.qseq.txt', '4_test.^.qseq.txt',
                                                      directory='tests/test_qseq/',
                                                      sample_key='tests/test_sample_files/dual_index_test.txt',
                                                      barcode_1='tests/test_sample_files/N700_nextera_barcodes.txt',
                                                      barcode_2='tests/test_sample_files/N500_nextera_barcodes.txt',
                                                      file_label='rbbr'
                                                      )
test_dual_index_demultiplex.run(b1_reverse=True, b2_reverse=True)

test_dual_metrics = iterate_through_qseq(demultiplex_instance=test_dual_index_demultiplex,
                                         output_directory='tests/test_output/', gnu_zipped=False)


class TestDemultiplex(unittest.TestCase):
    def setUp(self):
        pass

    def test_single_filter_pass(self):
        self.assertEqual(test_single_metrics[3], 19614)

    def test_single_unmatched_reads(self):
        self.assertEqual(test_single_metrics[1], 2774)

    def test_single_index_reads(self):
        self.assertEqual(test_single_metrics[2], 16840)

    def test_single_samples(self):
        x = test_single_index_demultiplex.file_list[0][0].replace('1_test', '')
        y = test_single_index_demultiplex.file_list[1][0].replace('2_test', '')
        self.assertEqual(x, y)

    def test_dual_filter_pass(self):
        self.assertEqual(test_dual_metrics[3], 19264)

    def test_total_reads(self):
        self.assertEqual(test_dual_metrics[0], test_single_metrics[0])

    def test_dual_unmatched_reads(self):
        self.assertEqual(test_dual_metrics[1], 9420)

    def test_dual_index_reads(self):
        self.assertEqual(test_dual_metrics[2], 9844)

    def test_dual_samples(self):
        x = test_dual_index_demultiplex.file_list[2][0].replace('3_test', '')
        y = test_dual_index_demultiplex.file_list[3][0].replace('4_test', '')
        self.assertEqual(x, y)

    def test_command_line_single_index(self):
        parser_open = subprocess.run(['python3', 'Demultiplex.py', '-D', 'tests/test_qseq/', '-S',
                                      'tests/test_sample_files/single_index_test.txt', '-B1',
                                      'tests/test_sample_files/N700_nextera_barcodes.txt', '-W', '2',
                                      '-L', 'rb', '-O', 'tests/test_output/',
                                      '-I', '1_test.^.qseq.txt', '2_test.^.qseq.txt'],
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
        parser_open = subprocess.run(['python3', 'Demultiplex.py', '-D', 'tests/test_qseq/', '-S',
                                      'tests/test_sample_files/dual_index_test.txt', '-B1',
                                      'tests/test_sample_files/N700_nextera_barcodes.txt', '-B2',
                                      'tests/test_sample_files/N500_nextera_barcodes.txt', '-L', 'rbbr', '-O',
                                      'tests/test_output/', '-I', '1_test.^.qseq.txt', '2_test.^.qseq.txt',
                                      '3_test.^.qseq.txt', '4_test.^.qseq.txt'],
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
