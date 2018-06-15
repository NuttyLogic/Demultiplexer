#!/usr/bin/env python3

from ProcessInputFiles import ParseInputFiles
from MultipleIterator import MultipleQseqIterator
from Hamming import HammingDistance
from FastqOutput import FastqOut
from DemultiplexHelpers import duplicates
from DemultiplexHelpers import qseq_fastq_conversion


class RunDemultiplex:

    def __init__(self, input_object=None, output_directory=None, ham_dist=1):
        """

        :type input_object: object
        """
        assert(isinstance(input_object, ParseInputFiles))
        self.input_object = input_object
        self.output_objects = None
        self.output_directory = output_directory
        self.ham_dist = ham_dist
        self.run_stats = {'reads': 0, 'indexed_reads': 0, 'unmatched_reads': 0, 'reads_passing_filter': 0}
        self.barcode_indices = None
        self.read_indices = None
        self.demultiplex_barcode_dicts = None
        self.hamming_dicts = None
        self.run()

    def run(self):
        self.set_indices()
        self.set_barcode_objects()
        self.get_output_objects()
        self.iterate_through_qseq()

    def set_indices(self):
        self.barcode_indices = duplicates(self.input_object.file_label, 'barcode')
        self.read_indices = duplicates(self.input_object.file_label, 'read')

    def get_output_objects(self):
        self.output_objects = FastqOut(sample_dict=self.input_object.sample_ids,
                                       read_count=self.input_object.read_count,
                                       output_directory=self.output_directory).output_dict

    def set_barcode_objects(self):
        self.demultiplex_barcode_dicts = [dict(x.barcode_dict) for x in self.input_object.barcodes]
        self.hamming_dicts = [dict(x.hamming_dict) for x in self.input_object.barcodes]

    def close_output_object(self):
        for out in self.output_objects:
            out.close()

    def get_seq_barcodes(self, line):
        return [line[index][8] for index in self.barcode_indices]

    def get_seq_reads(self, line):
        return [line[index] for index in self.read_indices]

    def set_new_barcode(self, barcode, index):
        hamming_list = []
        for ref_barcode, sample_id in self.hamming_dicts[index].items():
            ham_score = HammingDistance(barcode, ref_barcode).hamming
            hamming_list.append((ham_score, sample_id))
        hamming_list.sort()
        if hamming_list[0][0] != hamming_list[1][0] and hamming_list[0][0] <= self.ham_dist:
            hamming_id = hamming_list[0][1]
            self.demultiplex_barcode_dicts[index][barcode] = hamming_id
        else:
            hamming_id = 'unmatched'
            self.demultiplex_barcode_dicts[index][barcode] = hamming_id
        return hamming_id

    def output_read(self, read_id, seq_reads):
        try:
            outs = self.output_objects[self.input_object.sample_ids[read_id]]
        except KeyError:
            outs = self.output_objects['unmatched']
        for read, out in zip(seq_reads, outs):
            out.write(qseq_fastq_conversion(read))

    def iterate_through_qseq(self):
        for grouping in self.input_object.file_list:
            for qseq_lines in MultipleQseqIterator(grouping, self.input_object.qseq_directory):
                self.run_stats['reads'] += 1
                if self.get_read_quality(qseq_lines):
                    self.run_stats['reads_passing_filter'] += 1
                    barcodes = self.get_seq_barcodes(qseq_lines)
                    sample_id = []
                    for index, barcode in enumerate(barcodes):
                        try:
                            sample_id.append(self.demultiplex_barcode_dicts[index][barcode])
                        except KeyError:
                            sample_id.append(self.set_new_barcode(barcode, index))
                    if 'unmatched' in sample_id:
                        self.run_stats['unmatched_reads'] += 1
                    else:
                        self.run_stats['indexed_reads'] += 1
                    seq_reads = self.get_seq_reads(qseq_lines)
                    self.output_read('_'.join(sample_id), seq_reads)

    @staticmethod
    def get_read_quality(line):
        combined_filter = ''.join([qual[-1] for qual in line])
        if '0' not in combined_filter:
            return True
        else:
            return False
