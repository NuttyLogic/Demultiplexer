from MultipleIterator import *
from sequence_helpers import *
import sys


class BarcodeOffset:

    def __init__(self, input_list, demultiplex_instance, gnu_zipped):
        self.input_list = input_list
        self.class_instance = demultiplex_instance
        self.gnu_zipped = gnu_zipped
        self.run()

    def run(self):
        self.set_barcode_offset()

    def set_barcode_offset(self):
        """Iterate through the first 1,000,000 lines of the first qseq set to set the barcode length, if the barcode
        and read length are equal the function exits
        ---------------------------------------------
        sets class_instance (demultiplex_input instance) self.left_offset/right_offset"""
        iterator = MultipleSequencingFileIterator(self.input_list, directory=self.class_instance.directory,
                                                  gnu_zipped=self.gnu_zipped)
        # get position of barcode files
        barcode_indexes = duplicates(self.class_instance.file_label, 'barcode')
        # get position of read files
        # set barcode list, for looping
        barcode_list = [self.class_instance.barcode_1, self.class_instance.barcode_2]
        # loop through grouped files
        left_offset = 0
        right_offset = int(self.class_instance.right_offset)
        barcode_offset_set = False
        offset_list = None
        for count, line in enumerate(iterator.iterator_zip()):
            if count < 100000:
                # set string with Illumina quality control information
                combined_filter = ''.join([qual[-1] for qual in line])
                # initialize empty sample_id value
                # if all reads don't pass filter don't consider
                if '0' not in combined_filter:
                    # loop through barcode_indexes, get sample key
                    for index_count, index in enumerate(barcode_indexes):
                        # get sequence location in qseq file, and set barcode offset if needed
                        barcode_read = line[barcode_indexes[index_count]][8]
                        if self.class_instance.right_offset != len(barcode_read):
                            if self.class_instance.right_offset > len(barcode_read):
                                print('Barcode length greater than read length')
                                sys.exit()
                            elif self.class_instance.right_offset < len(barcode_read):
                                offset = len(barcode_read) - self.class_instance.right_offset
                                if not offset_list:
                                    offset_list = [0 for _ in range(offset + 1)]
                                for offset_count, offset_range in enumerate(range(offset + 1)):
                                    try:
                                        _ = (barcode_list[index_count][line[barcode_indexes[index_count]][8]
                                        [offset_range:(self.class_instance.right_offset + offset_range)]])
                                        offset_list[offset_count] = offset_list[offset_count] + 1
                                    except KeyError:
                                        continue
                        else:
                            barcode_offset_set = True
                            break
                if barcode_offset_set:
                    break
            else:
                break
        if not barcode_offset_set:
            left_offset = offset_list.index(max(offset_list))
            right_offset = self.class_instance.right_offset + left_offset
        self.class_instance.left_offset = left_offset
        self.class_instance.right_offset = right_offset
