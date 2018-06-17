import unittest

from DemultiplexHelpers import duplicates
from Hamming import HammingDistance


barcode_list = ['TAAGGCGA',	'CGTACTAG',	'AGGCAGAA',	'TCCTGAGC',	'GGACTCCT',	'TAGGCATG']
barcode_test_hams = ['TAAGGCGA', 'TAAGGC', 'AGGCGA', 'TxxGGCGA', 'TAAGGCxx', 'xxxxxxxx']


def test_ham(barcode_list, barcode_test_hams, mixed_length=False):
    hammed_list = []
    for test_barcode in barcode_test_hams:
        hamming_list = []
        for ref_barcode in barcode_list:
            ham_score = HammingDistance(reference_barcode=ref_barcode,
                                        seq_barcode=test_barcode,
                                        mixed_length=mixed_length).hamming
            print(ref_barcode, test_barcode, ham_score)
            hamming_list.append((ham_score, ref_barcode))
        hamming_list.sort()
        hammed_list.append(hamming_list)
    hammed_list.sort()
    print(hammed_list)
    return hammed_list


hamming_test_list = test_ham(barcode_list, barcode_test_hams)
print('Hamming with mixed length')
hamming_test_list_mixed = test_ham(barcode_list, barcode_test_hams, mixed_length=True)

test_indices = duplicates(['read', 'barcode', 'read', 'blah', 'read', 'blah', 'barcodes'], 'read')


class TestDemultiplex(unittest.TestCase):

    def setUp(self):
        pass

    def test_hamming_same(self):
        self.assertEqual(hamming_test_list[0][0][0], 0)

    def test_hamming_left_shorter(self):
        self.assertEqual(hamming_test_list[1][0][0], 0)

    def test_hamming_right_shorter(self):
        self.assertEqual(hamming_test_list[2][0][0], 0)

    def test_hamming_middle_two(self):
        self.assertEqual(hamming_test_list[3][0][0], 2)

    def test_hamming_right_two(self):
        self.assertEqual(hamming_test_list[4][0][0], 2)

    def test_hamming_all(self):
        self.assertEqual(hamming_test_list[5][0][0], 8)

    def test_mixed_length_left(self):
        self.assertEqual(hamming_test_list_mixed[1][0][0], 2)

    def test_mixed_length_right(self):
        self.assertEqual(hamming_test_list_mixed[2][0][0], 2)

    def test_duplicate(self):
        self.assertEqual(test_indices, [0, 2, 4])


if __name__ == '__main__':
    unittest.main()