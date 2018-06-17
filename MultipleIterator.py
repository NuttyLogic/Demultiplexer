#!/usr/bin/env python3

from QseqIterator import QseqIterator


class MultipleQseqIterator:
    """Iterate through group of qseq files
    Arguments:
        input_files (list): list of files to iterate over
        directory (str): path to qseq directory
    Attributes:
        self.iter_list (list): list of QseqIterator objects
        """

    def __init__(self, input_files=None, directory=None):
        """Initiate iteration object, yield line in gseq files
         -----------------------------------------------------
         *args='path_to_gesq': returns an iterator object for paired sequencing files
         """
        assert(isinstance(input_files, list))
        assert(isinstance(directory, str))
        self.iter_list = []
        # store files in list
        for file in input_files:
            self.iter_list.append(QseqIterator(qseq='%s%s' % (directory, file)))

    def __iter__(self):
        for line in zip(*self.iter_list):
            yield line
