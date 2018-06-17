#!/usr/bin/env python3

import os
import sys
from BarcodeParser import BarcodeFileParser


class ParseInputFiles:
    """
    Parse external sample key, and store args, input file is formatted as 'barcode\tbarcode\tsample_id' where barcode
    the number of barcodes can be value >= 1
    Arguments:
        *args (str): strs specifying read file format ^ is used for read pair identifier
        qseq_directory (str): path to input files
        sample_file (str): path to file listing barcodes and sample ids
        barcode_reverse (bool): consider reverse complement
        file_label (str): str labeling input reads as sequencing reads or barcodes 'r' for seq read 'b' for barcode
                          labels should match to input order of the *args
    Attributes:
        self.file_label (str): string describing input files
        self.barcode_count (int): count of barcode in file_lablel
        self.file_list (list): list of paired qseq files
        self.file_description (list): input *args split by '^'
        self.qseq_directory (str): path to qseq directory
        self.sample_file (str): path to sample input file
        self.barcode_reverse (bool): consider reverse complement
        self.file_label (str): str describing input files as reads or barcodes
        self.barcodes (list): list of index barcodes
        self.sample_ids (dict): dict of sample name linked to barcode combination
        self.read_count (int):  count of input seq files
        self.barcode_count (int) count of input barcode files

    """

    def __init__(self, *args, qseq_directory=None, sample_file=None,
                 barcode_reverse=False, file_label=None):
        self.file_list = []
        self.file_description = []
        for arg in args:
            self.file_description.append(arg.split('^'))
        self.qseq_directory = qseq_directory
        self.sample_file = sample_file
        self.barcode_reverse = barcode_reverse
        self.file_label = file_label
        self.barcodes = None
        self.sample_ids = {}
        self.read_count = None
        self.barcode_count = None
        self.run()

    def run(self):
        self.get_demultiplex_info()
        self.set_barcode_objects()
        self.get_directory_lists()
        self.set_file_label()
        self.set_iteration_list()

    def get_demultiplex_info(self):
        with open(self.sample_file, 'r') as sample:
            while True:
                line = self.process_line(sample.readline())
                if not line[0]:
                    break
                if not self.barcodes:
                    self.barcodes = [[] for _ in range(len(line) - 1)]
                for count, barcode in enumerate(line[:-1]):
                    self.barcodes[count].append(barcode.upper())
                sample_id = '_'.join(line[:-1])
                self.sample_ids[sample_id] = line[-1]

    def set_barcode_objects(self):
        barcode_objects = []
        assert(isinstance(self.barcodes, list))
        for barcodes in self.barcodes:
            if self.barcode_reverse:
                barcode_objects.append(BarcodeFileParser(barcode_list=barcodes, rev=True))
            else:
                barcode_objects.append(BarcodeFileParser(barcode_list=barcodes, rev=False))
        self.barcodes = barcode_objects

    def get_directory_lists(self):
        """Link to directory, pull  list of files, combine files with same unique id. Function only works with files
        that have an integer as a unique id
        """
        file_list = os.listdir(self.qseq_directory)
        # initialize list to hold sample names (ie. coupled file IDs)
        sample_names = [[] for _ in range(len(self.file_description))]
        # key to sort files bases on proper ID
        sorting_key = [[] for _ in range(len(self.file_description))]
        for count, file_title in enumerate(self.file_description):
            for file_name in file_list:
                # order files based on read type, ie all read 1 go together etc.
                if file_title[0] in file_name:
                    sample_names[count].append(file_name)
                    # store unique file ID
                    try:
                        sort_id = int((file_name.replace(file_title[0], '')).replace(file_title[1], ''))
                    except ValueError:
                        print('File ID must be integer')
                        sys.exit()
                    # store sort id in list for file type
                    sorting_key[count].append(sort_id)
        for count, seq_file_list in enumerate(sample_names):
            # sort list on unique ID and append sorted list of file, one list per file prefix
            self.file_list.append([x for x, y in sorted(zip(sample_names[count], sorting_key[count]))])

    def set_file_label(self):
        """Parses string describing input files, the file_label should be formatted as r for read and b for barcode,
        so a string 'rbbr' would describe a four qseq file input with read barcode barcode read.  Used to properly
        parse the qseq files
       """
        label_list = []
        for character in self.file_label:
            if character.lower() == 'r':
                label_list.append('read')
            elif character.lower() == 'b':
                label_list.append('barcode')
        self.file_label = label_list
        self.barcode_count = label_list.count('barcode')
        self.read_count = label_list.count('read')

    def set_iteration_list(self):
        self.file_list = list(map(list, zip(*self.file_list)))

    @staticmethod
    def process_line(line):
        return line.replace('\n', '').split('\t')
