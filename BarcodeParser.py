#!/usr/bin/env python3


from Barcode import Barcode


class BarcodeFileParser:

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
