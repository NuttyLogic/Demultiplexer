#!/usr/bin/env python3

import unittest
import DemultiplexClass
import subprocess


test_single_index_demultiplex = DemultiplexClass.Demuliplex('1_test.^.qseq.txt', '2_test.^.qseq.txt',
                                                            directory='tests/test_qseq/',
                                                            sample_key='tests/test_sample_files/single_index_test.txt',
                                                            barcode_1='tests/test_sample_files/'
                                                                      'N700_nextera_barcodes.txt',
                                                            file_label='rb',
                                                            )
test_single_index_demultiplex.run(output_directory='tests/test_output/')


test_dual_index_demultiplex = DemultiplexClass.Demuliplex('1_test.^.qseq.txt', '2_test.^.qseq.txt',
                                                          '3_test.^.qseq.txt', '4_test.^.qseq.txt',
                                                          directory='tests/test_qseq/',
                                                          sample_key='tests/test_sample_files/dual_index_test.txt',
                                                          barcode_1='tests/test_sample_files/N700_nextera_barcodes.txt',
                                                          barcode_2='tests/test_sample_files/N500_nextera_barcodes.txt',
                                                          file_label='rbbr'
                                                          )
test_dual_index_demultiplex.run(output_directory='tests/test_output/')


class TestDemultiplex(unittest.TestCase):

    def setUp(self):
        pass

    def test_single_filter_pass(self):
        self.assertEqual(test_single_index_demultiplex.reads_pass_filter, 19614)

    def test_single_unmatched_reads(self):
        self.assertEqual(test_single_index_demultiplex.unmatched_read, 2774)

    def test_single_index_reads(self):
        self.assertEqual(test_single_index_demultiplex.indexed_reads, 16840)

    def test_single_samples(self):
        x = test_single_index_demultiplex.file_list[0][0].replace('1_test', '')
        y = test_single_index_demultiplex.file_list[0][1].replace('2_test', '')
        self.assertEqual(x, y)

    def test_dual_filter_pass(self):
        self.assertEqual(test_dual_index_demultiplex.reads_pass_filter, 19264)

    def test_total_reads(self):
        self.assertEqual(test_single_index_demultiplex.reads, test_dual_index_demultiplex.reads)

    def test_dual_unmatched_reads(self):
        self.assertEqual(test_dual_index_demultiplex.unmatched_read, 9420)

    def test_dual_index_reads(self):
        self.assertEqual(test_dual_index_demultiplex.indexed_reads, 9844)

    def test_dual_samples(self):
        x = test_dual_index_demultiplex.file_list[1][2].replace('3_test', '')
        y = test_dual_index_demultiplex.file_list[1][3].replace('4_test', '')
        self.assertEqual(x, y)

    def test_command_line_single_index(self):
        parser_open = subprocess.run(['python3', 'DemultiplexRun.py', '-D', 'tests/test_qseq/', '-S',
                            'tests/test_sample_files/single_index_test.txt', '-B1',
                            'tests/test_sample_files/N700_nextera_barcodes.txt', '-L', 'rb', '-O', 'tests/test_output/',
                            '-I', '1_test.^.qseq.txt', '2_test.^.qseq.txt'],
                            stdout=subprocess.PIPE)
        output = parser_open.stdout
        print(output)
        output_categories = (output.decode()).split('\n')
        filter = int(output_categories[2].split(':')[1])
        indexed = int(output_categories[3].split(':')[1])
        unmatched = int(output_categories[4].split(':')[1])
        self.assertEqual(filter, 19614)
        self.assertEqual(indexed, 16840)
        self.assertEqual(unmatched, 2774)

    def test_command_line_dual_index(self):
        parser_open = subprocess.run(['python3', 'DemultiplexRun.py', '-D', 'tests/test_qseq/', '-S',
                                      'tests/test_sample_files/dual_index_test.txt', '-B1',
                                      'tests/test_sample_files/N700_nextera_barcodes.txt', '-B2',
                                      'tests/test_sample_files/N500_nextera_barcodes.txt', '-L', 'rbbr', '-O',
                                      'tests/test_output/', '-I', '1_test.^.qseq.txt', '2_test.^.qseq.txt',
                                      '3_test.^.qseq.txt', '4_test.^.qseq.txt'],
                                     stdout=subprocess.PIPE)
        output = parser_open.stdout
        print(output)
        output_categories = (output.decode()).split('\n')
        filter = int(output_categories[2].split(':')[1])
        indexed = int(output_categories[3].split(':')[1])
        unmatched = int(output_categories[4].split(':')[1])
        self.assertEqual(filter, 19264)
        self.assertEqual(indexed, 9844)
        self.assertEqual(unmatched, 9420)


if __name__ == '__main__':
    unittest.main()
