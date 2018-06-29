#!/usr/bin/env python3


class HammingDistance:
    """Calculates the similarity between a reference barcode and a sequencing read, if the reference barcode and the
    sequencing barcode are different lengths hamming distance is calculated for equal length strings relative to the
    beginning and the end of the longer string
    Arguments:
        reference_barcode (str): reference barcode string
        seq_barcode (str): sequence barcode string
        mixed_length (bool): if a list of barcodes is composed of different length string the different in length
                            between the reference and sequencing barcode is added to the hamming distance
    Attributes:
        self.reference_barcode (str): reference barcode str
        self.seq_barcode (str): sequencing barcode str
        self.hamming (int): store calculated hamming distance
        self.mixed_length (bool): set hamming distance plus len difference
        self.get_same_hamming (func): set hamming distance for strs of same length
        self.get_diff_hamming (func): set hamming distance for different length str
        self.hamming_distance (func): returns count of str differences
      """

    def __init__(self, reference_barcode=None, seq_barcode=None, mixed_length=False):
        self.reference_barcode = str(reference_barcode)
        self.seq_barcode = str(seq_barcode)
        self.hamming = None
        self.mixed_length = mixed_length
        if len(self.reference_barcode) == len(self.seq_barcode):
            self.get_same_hamming()
        else:
            self.get_diff_hamming()

    def get_same_hamming(self):
        self.hamming = self.hamming_dist(self.reference_barcode, self.seq_barcode)

    def get_diff_hamming(self):
        hamming_list = []
        if len(self.reference_barcode) < len(self.seq_barcode):
            difference = len(self.seq_barcode) - len(self.reference_barcode)
            # set str length starting at difference len
            for count in range(difference + 1):
                start = count
                end = len(self.seq_barcode) - difference + count
                hamming_list.append(self.hamming_dist(self.reference_barcode,
                                                      self.seq_barcode[start:end]))
            if self.mixed_length:
                self.hamming = sorted(hamming_list)[0] + abs(difference)
            else:
                self.hamming = sorted(hamming_list)[0]
        else:
            difference = len(self.reference_barcode) - len(self.seq_barcode)
            for count in range(difference + 1):
                start = count
                end = len(self.reference_barcode) - difference + count
                hamming_list.append(self.hamming_dist(self.reference_barcode[start:end],
                                                      self.seq_barcode))
            if self.mixed_length:
                self.hamming = sorted(hamming_list)[0] + abs(difference)
            else:
                self.hamming = sorted(hamming_list)[0]

    @staticmethod
    def hamming_dist(string1, string2):
        # count str differences
        return sum(chr1 != chr2 for chr1, chr2 in zip(string1, string2))
