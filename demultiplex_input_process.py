#!/usr/bin/env python3

from BarcodeParser import BarcodeFileParser
import os
import sys


class ProcessDemultiplexInput:
    """Process Input Files"""

    def __init__(self, *args, directory='path', sample_key='path', file_label='', barcode_1=None,
                 barcode_2=None):
        # store file description
        self.file_description = []
        for arg in args:
            self.file_description.append(arg.split('^'))
        self.directory = directory
        self.sample_key = sample_key
        self.file_label = file_label
        self.file_list = []
        self.sample_list = []
        self.barcode_1 = barcode_1
        self.barcode_2 = barcode_2
        self.barcode_count = None
        self.read_count = None
        self.left_offset = 0
        self.right_offset = 0

    def run(self, b1_reverse=False, b2_reverse=False):
        self.get_directory_lists()
        self.process_barcodes(barcode1_reverse=b1_reverse, barcode2_reverse=b2_reverse)
        self.process_file_label()
        self.get_sample_labels()

    def get_barcodes(self, barcode1_reverse=False, barcode2_reverse=False):
        """If barcode file supplied process files and store values in dictionary
        -----------------------------------------------------
        opens self.barcode*, a path to text file with a new barcode on each line
        returns self.barcode*, a dictionary hashing barcodes to an Illumina ID
        returns self.left_offset, set barcode read offset if barcode length differs from read length
        returns self.right_offset, default barcode length equal to read length"""
        self.barcode_1 = BarcodeFileParser(barcode_file_path=self.barcode_1, rev=barcode1_reverse)
        if self.barcode_2:
            self.barcode_2 = BarcodeFileParser(barcode_file_path=self.barcode_2, rev=barcode2_reverse)

    def get_sample_labels(self):
        """Takes sample label file and processes barcode sample IDs.  Note this function assumes 'barcode1 \t barcode 2
        \t sample_name \n' .  If only one barcode is used then the sample file should be formatted as
        'barcode1 \t sample_name \n' and no barcode_2 file should be supplied. Function will fail if sample key is not
        in this format.
        -----------------------------------------------------
        opens self.sample_key; parses files and hashes to sample name based on unique barcode ID
        returns self.sample_key; a dictionary hashing to sample name"""
        # initialize dict
        sample_dict = {}
        # loop through text file
        for line in open(self.sample_key):
            # replace new line indicator
            line_replace = line.replace('\n', '')
            # split on tabs
            line_split = line_replace.split('\t')
            # if barcode_2 supplied sample id is a combination of 2 barcodes
            if self.barcode_2:
                # int(), to check if file in is proper format
                try:
                    sample_dict['key' + str(int(line_split[0])) + 'key' + str(int(line_split[1]))] = line_split[2]
                except ValueError:
                    print('Sample Key not in proper format\nPlease format file Barcode1 tab Barcode2 tab SampleName')
                    sys.exit()
                self.sample_list.append(line_split[2])
            # else id is only one barcode
            else:
                try:
                    sample_dict['key' + str(int(line_split[0]))] = line_split[1]
                except ValueError:
                    print('Sample Key not in proper format\nPlease format file Barcode tab SampleName')
                    sys.exit()
                self.sample_list.append(line_split[1])
        self.sample_key = sample_dict

    def get_directory_lists(self):
        """Link to directory, pull  list of files, combine files with same unique id. Function only works with files
        that have an integer as a unique id
        -----------------------------------------------------
        self.directory: path to illumina sequencing lane directory
        self.file_description: list of lists containing sequencing file prefix and suffix, sequencing file designated
        by prefix*suffix
        returns; sorted list of relevant files names in directory"""
        file_list = os.listdir(self.directory)
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

    def process_file_label(self):
        """Parses string describing input files, the file_label should be formatted as r for read and b for barcode,
        so a string 'rbbr' would describe a four qseq file input with read barcode barcode read.  Used to properly
        parse the qseq files
        -----------------------------------------------------
        self.file_label: string describing input files
        self.barcode_count: returns count of barcode in file_lablel as a downstream control"""
        label_list = []
        for character in self.file_label:
            if character.lower() == 'r':
                label_list.append('read')
            elif character.lower() == 'b':
                label_list.append('barcode')
        self.file_label = label_list
        self.barcode_count = label_list.count('barcode')
        self.read_count = label_list.count('read')
