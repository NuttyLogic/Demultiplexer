#!/usr/bin/env python3


class HammingDistance:
    """Calculates the similarity between a reference barcode and a sequencing read
     - because barcode sequences should be continuous """

    def __init__(self, reference_barcode=None, seq_barcode=None):
        self.reference_barcode = reference_barcode
        self.seq_barcode = seq_barcode
        self.hamming = None
        if len(self.reference_barcode) == len(self.seq_barcode):
            self.get_same_hamming()
        else:
            self.get_diff_hamming()

    def get_same_hamming(self):
        self.hamming = self.hamming_dist(self.reference_barcode, self.seq_barcode)

    def get_diff_hamming(self):
        if self.reference_barcode < self.seq_barcode:
            difference = len(self.seq_barcode) - len(self.reference_barcode)
            hamming1 = self.hamming_dist(self.reference_barcode, self.seq_barcode[difference:])
            hamming2 = self.hamming_dist(self.reference_barcode, self.seq_barcode[:-difference])
            self.hamming = sorted([hamming1, hamming2])[0] + abs(difference)
        else:
            difference = len(self.reference_barcode) - len(self.seq_barcode)
            hamming1 = self.hamming_dist(self.reference_barcode[difference:], self.seq_barcode)
            hamming2 = self.hamming_dist(self.reference_barcode[:-difference], self.seq_barcode)
            self.hamming = sorted([hamming1, hamming2])[0] + abs(difference)

    @staticmethod
    def hamming_dist(string1, string2):
        return sum(chr1 != chr2 for chr1, chr2 in zip(string1, string2))
