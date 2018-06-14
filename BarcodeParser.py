#!/usr/bin/env python3


from Barcode import Barcode


class BarcodeFileParser:

    def __init__(self, barcode_file_path=None, rev=False):
        self.barcode_dict = {}
        self.hamming_dict = {}
        self.rev = rev
        self.get_barcodes(barcode_file_path=barcode_file_path)

    def get_barcodes(self, barcode_file_path=None):
        with open(barcode_file_path, 'r') as bar_file:
            for line in bar_file:
                line_info = self.process_line(line)
                barcode = Barcode(line_info[0], line_info[1])
                if self.rev:
                    self.barcode_dict[str(barcode)] = barcode.get_number()
                    self.barcode_dict[reversed(barcode)] = barcode.get_number()
                    self.hamming_dict[str(barcode)] = barcode.get_number()
                    self.hamming_dict[reversed(barcode)] = barcode.get_number()
                else:
                    self.barcode_dict[str(barcode)] = [barcode]
                    self.hamming_dict[str(barcode)] = barcode.get_number()

    @staticmethod
    def process_line(line):
        return line.replace('\n', '').split('\t')
