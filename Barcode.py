#!/usr/bin/env python3


class Barcode:

    def __init__(self, barcode, number):
        self.barcode = barcode
        self.number = number
        self.complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

    def __reversed__(self):
        """Simple reverse complement function used to initialize a barcode dictionary, (original sequence and reverse
                complement both link to same hash in dictionary).
                -----------------------------------------------------
                string='string': string must be composed of BP ATGC
                returns; reverse complement string"""
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

    def reverse(self):
        """Simple reverse complement function used to initialize a barcode dictionary, (original sequence and reverse
        complement both link to same hash in dictionary).
        -----------------------------------------------------
        string='string': string must be composed of BP ATGC
        returns; reverse complement string"""
        # reverse string
        reversed_string = []

        self.barcode = self.barcode[::-1]
        # complementary bp lookup dictionary
        bases = list(self.barcode)
        # iterate over string list
        for i in bases:
            reversed_string.append(self.complement[i])

        # return joined string
        self.barcode = ''.join(reversed_string)
        
        return self.barcode

    def get(self):
        return self.barcode

    def get_number(self):
        return self.number
