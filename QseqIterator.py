#!/usr/bin/env python3

import gzip
import io


class QseqIterator:
    """Initializes a qseq file iterator, opens gziped or plain text
    Arguments:
        qseq (str): path to qseq file
    Attributes:
        self.qseq (object): iteration object
        self.process_line (func): return split qseq line
        """

    def __init__(self, qseq=None):
        if qseq.endswith(".gz"):
            self.qseq = gzip.open(qseq, 'rb')
        else:
            self.qseq = open(qseq, 'r')

    def __iter__(self):
        with self.qseq as qseq:
            while True:
                line = qseq.readline()
                if not line:
                    break
                yield self.process_line(line)

    @staticmethod
    def process_line(line):
        if isinstance(line, bytes):
            return line.decode('utf-8').replace('\n', '').split('\t')
        else:
            return line.replace('\n', '').split('\t')
