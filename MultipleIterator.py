#!/usr/bin/env python3

from QseqIterator import QseqIterator


class MultipleQseqIterator:
    """Open qseq files together and iterate over them as a group"""

    def __init__(self, input_files=None, directory=None):
        """Initiate iteration object, yield line in gseq files
         -----------------------------------------------------
         *args='path_to_gesq': returns an iterator object for paired sequencing files
         """
        assert(input_files, list)
        assert(directory, str)
        self.iter_list = []
        # store files in list
        for file in input_files:
            self.iter_list.append(QseqIterator(qseq='%s%s' % (directory, file)))

    def __iter__(self):
        for line in zip(*self.iter_list):
            yield line
