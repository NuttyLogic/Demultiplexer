#!/usr/bin/env python3


from Barcode import Barcode


class BarcodeFileParser:
    """Class to organize input barcodes into hamming objects and reference objects
    Arguments:
        barcode_list (list): list of str containing all barcodes for an index
        rev (bool): consider reverse complement of barcodes
    Attributes:
        self.barcode_dict (dict): barcodes hashing to id, expanded with sequencing barcodes
        self.hamming_dict (dict): barcodes hashing to id, not expanded used as reference
        self.rev (bool): consider reverse complement of barcodes
        self.barcode_list (list): list of input barcodes
        self.get_barcodes (func): set barcode_dict and hamming dict
        """

    def __init__(self, barcode_list=None, rev=False):
        self.barcode_dict = {}
        self.hamming_dict = {}
        self.rev = rev
        self.barcode_list = barcode_list
        self.get_barcodes()

    def get_barcodes(self):
        for bc in self.barcode_list:
            barcode = Barcode(bc, bc)
            if self.rev:
                self.barcode_dict[str(barcode)] = barcode.get_id()
                self.barcode_dict[reversed(barcode)] = barcode.get_id()
                self.hamming_dict[str(barcode)] = barcode.get_id()
                self.hamming_dict[reversed(barcode)] = barcode.get_id()
            else:
                self.barcode_dict[str(barcode)] = barcode.get_id()
                self.hamming_dict[str(barcode)] = barcode.get_id()
