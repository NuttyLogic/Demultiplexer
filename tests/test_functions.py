import subprocess
import unittest

from ProcessInputFiles import ParseInputFiles
from DemultiplexRun import RunDemultiplex
from Hamming import HammingDistance


barcode_list = ['TAAGGCGA',	'CGTACTAG',	'AGGCAGAA',	'TCCTGAGC',	'GGACTCCT',	'TAGGCATG']
barcode_test_hams = ['TAAGGCGA', 'TAAGGC', 'AGGCGA', 'TxxGGCGA', 'TAAGGCxx', 'xxxxxxx']


def test_ham(barcode_list, barcode_test_hams):
    hamming_list = []
    for ref_barcode, test_barcode in zip(barcode_list, barcode_test_hams):
        ham_score = HammingDistance(test_barcode, ref_barcode).hamming
        hamming_list.append((ham_score, ref_barcode))
    hamming_list.sort()
    print(hamming_list)

test_ham(barcode_list, barcode_test_hams)


class TestDemultiplex(unittest.TestCase):
    def setUp(self):
        pass




if __name__ == '__main__':
    unittest.main()