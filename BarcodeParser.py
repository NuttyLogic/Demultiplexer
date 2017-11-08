#!/usr/bin/env python3
from Barcode import Barcode
from collections import Counter


class BarcodeFileParser:

    def __init__(self,
                 path=None,
                 mismatch=1):

        self.path = path
        self.mismatch = 0 if mismatch <= 0 else mismatch
        self.mismatch_list = ['A', 'T', 'G', 'C', '.']
        self.barcode_list = []

    def read_barcodes(self):
        for index, line in enumerate(open(self.path)):
            barcode = self.read_barcode_from_line(line=line)
            self.barcode_list.append(barcode)

    def get_barcodes(self):
        if not self.path:
            return None
        else:
            self.read_barcodes()
            mismatched_barcodes = self.get_possible_mismatches()
            barcode_dict = self.get_barcodes_from_mismatch_lists(mismatched_barcodes)
            return barcode_dict

    @staticmethod
    def mismatch(barcode_string, mismatches):
        mismatched_barcodes = []
        for possible_mismatch in mismatches:
            # loop over every position in barcode
            for index in range(len(barcode_string)):
                barcode_list = list(barcode_string)
                # set character in to mismatch
                barcode_list[index] = possible_mismatch
                # if mismatched barcode not already in list then add the barcode to the list
                mismatched_barcodes.append(''.join(barcode_list))
        return mismatched_barcodes

    def read_barcode_from_line(self, line=None):
        parsed = (line.replace('\n', '')).split('\t')
        return Barcode(parsed[0], int(parsed[1]))

    def get_possible_mismatches(self):
        barcode_list = []
        barcode_higher = []
        barcode_lower = []

        for count, barcode in enumerate(self.barcode_list):
            barcode_list.append([[barcode.get(), barcode.reverse()]])
            barcode_higher.append(str(barcode.get()))
            barcode_higher.append(str(barcode.reverse()))
            forward = BarcodeFileParser.mismatch(barcode.get(), self.mismatch_list)
            reverse = BarcodeFileParser.mismatch(barcode.reverse(), self.mismatch_list)
            all_barcodes = list(set(forward + reverse))
            barcode_list[count].append(all_barcodes)
            barcode_lower = barcode_lower + all_barcodes

        counted_lower = list(Counter(barcode_lower).items())
        validation_list = []

        for value in counted_lower:
            if value[1] == 1 and value[0] not in barcode_higher:
                validation_list.append(value[0])

        for barcode_set in barcode_list:
            for barcode in barcode_set[1]:
                if barcode not in validation_list:
                    barcode_set[1].remove(barcode)

        barcode_higher = barcode_higher + validation_list

        for mismatch_count in range(2):
            barcode_lower = []
            for barcode_set in barcode_list:
                barcode_set.append([])
                for barcode in barcode_set[1 + mismatch_count]:
                    mismatches = BarcodeFileParser.mismatch(barcode, self.mismatch_list)
                    barcode_set[2 + mismatch_count] = barcode_set[2 + mismatch_count] + mismatches
                    barcode_lower = barcode_lower + mismatches

            counted_lower = list(Counter(barcode_lower).items())
            validation_list = []

            for value in counted_lower:
                if value[1] == 1 and value[0] not in barcode_higher:
                    validation_list.append(value[0])

            for barcode_set in barcode_list:
                for barcode in barcode_set[2 + mismatch_count]:
                    if barcode not in validation_list:
                        barcode_set[2 + mismatch_count].remove(barcode)
            barcode_higher = barcode_higher + validation_list

        return barcode_list

    def get_barcodes_from_mismatch_lists(self, barcode_list):
        barcode_dict = {}
        cleaned_barcodes = []
        len_of_barcode_list = len(barcode_list[0]) - 1

        for barcode_set in barcode_list:
            cleaned_barcode_set = []
            for count in range(len_of_barcode_list):
                cleaned_barcode_set = cleaned_barcode_set + barcode_set[count]
            cleaned_barcodes.append(cleaned_barcode_set)

        for barcode_info in zip(cleaned_barcodes, self.barcode_list):
            barcode_index = barcode_info[1].get_number()
            for barcode in barcode_info[0]:
                barcode_dict[barcode] = barcode_index

        return barcode_dict