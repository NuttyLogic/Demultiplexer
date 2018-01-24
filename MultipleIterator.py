#!/usr/bin/env python3

import gzip
import io


class MultipleSequencingFileIterator:
    """Open multiple fastq or gseq files together and iterate over them as a group"""

    def __init__(self, input_file, directory='path', gnu_zipped=False):
        """Initiate iteration object, yield line in gseq files
         -----------------------------------------------------
         *args='path_to_gesq': returns an iterator object for paired sequencing files
         gnu_zipped=False: if gnu_zipped=True will process files with python gzip library, increases processing time
         unix gunzip is faster"""
        file_list = []
        # store files in list
        for file in input_file:
            file_list.append(directory + file)

        # iterator initializer
        def iteration_call(iter_file):
            if gnu_zipped:
                # wrap encoded line in buffer to stream text file
                with io.TextIOWrapper(io.BufferedReader(gzip.open(iter_file, 'rb'))) as seq:
                    for line in seq:
                        # yield a line split by tabs and stripped of line identifier, '\n'
                        yield ((line.replace('\n', '')).split('\t'))
            else:
                with open(iter_file) as seq:
                    for line in seq:
                        # yield a line split by tabs and stripped of line identifier, '\n'
                        yield ((line.replace('\n', '')).split('\t'))
        self.iter_list = []
        # append iterator object to list, object is iterator for individual files
        for file in file_list:
            self.iter_list.append(iteration_call(file))

    # zip files together to iterate over all of the files at once and yield one line at a time for looping
    def iterator_zip(self):
        for line in zip(*self.iter_list):
            yield line
