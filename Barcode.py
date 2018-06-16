#!/usr/bin/env python3


class Barcode:
    """Wrapper for barcodes inputs
    Arguments:
        barcode (str): str listing barcode sequence
        id (str/int): id for reference barcode
        """

    def __init__(self, barcode, id):
        self.barcode = barcode
        self.id = id
        self.complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

    def __reversed__(self):
        """Returns reverse complement of barcode -> str"""
        # reverse string
        reversed_string = []

        rev_barcode = self.barcode[::-1]
        # complementary bp lookup dictionary
        bases = list(rev_barcode)
        # iterate over string list
        for i in bases:
            reversed_string.append(self.complement[i])

        # return joined string
        return ''.join(reversed_string)

    def __str__(self):
        """Returns barcode sequence"""
        return self.barcode

    def get_id(self):
        """Returns barcode id """
        return self.id
