#!/usr/bin/env python3

from ProcessInputFiles import ParseInputFiles
from MultipleIterator import MultipleQseqIterator
from Hamming import HammingDistance
from FastqOutput import FastqOut
from DemultiplexHelpers import duplicates
from DemultiplexHelpers import qseq_fastq_conversion


class RunDemultiplex:
    """Parsers qseq files and outputs formatted sample fastq files
    Arguments:
        input_object (object): instance of ParseInputFiles class with all class attributes
        output_directory (str): path to directory to output sample fastq files
        ham_dist (int): hamming threshold to consider sequencing barcode
        mixed_length (bool): if barcodes for an index are mixed length include distance difference in hamming value
    Attributes:
        self.input_object (ParseInputFiles): instance of ParseInputFiles
        self.output_object (dict): dict of objects used to write formatted fastq lines
        self.output_directory (str): path to output directory
        self.ham_dict (int): hamming threshold
        self.mixed_length (bool): barcodes of mixed length
        self.run_stat (dict): counts of reads, indexed_reads, unmatched_reads, reads_passing_filter
        self.read_indices (list): indices of sequence read files
        self.barcode_indices (list): indices of barcode read files
        self.demultiplex_barcode_dicts (list): list of barcode hashing to barcode ids
        self.hamming_dicts (list): list of barcode dicts used to calculate hamming distance
    """

    def __init__(self, input_object=None, output_directory=None, ham_dist=1, mixed_length=False):
        assert(isinstance(input_object, ParseInputFiles))
        self.input_object = input_object
        self.output_objects = None
        self.output_directory = output_directory
        self.ham_dist = ham_dist
        self.mixed_length = mixed_length
        self.run_stats = {'reads': 0, 'indexed_reads': 0, 'unmatched_reads': 0, 'reads_passing_filter': 0}
        self.barcode_indices = None
        self.read_indices = None
        self.demultiplex_barcode_dicts = None
        self.hamming_dicts = None
        self.run()

    def run(self):
        """Run demultiplex
        Attributes:
            self.set_indices (func): retrieve read and barcode indices
            self.set_barcode_objects (func): copy demultiplex barcode dicts
            self.get_output_objects (func): initialize output objects
            self.iterate_through_qseq (func): process qseq files
            self.close_output_objects (func): close output objects
            """
        self.set_indices()
        self.set_barcode_objects()
        self.get_output_objects()
        self.iterate_through_qseq()
        self.close_output_object()

    def set_indices(self):
        """run DemultiplexHelpers duplicats to se list of indices"""
        self.barcode_indices = duplicates(self.input_object.file_label, 'barcode')
        self.read_indices = duplicates(self.input_object.file_label, 'read')

    def get_output_objects(self):
        """Initialize lists of output objects linked to sample_ids"""
        self.output_objects = FastqOut(sample_dict=self.input_object.sample_ids,
                                       read_count=self.input_object.read_count,
                                       output_directory=self.output_directory).output_dict

    def set_barcode_objects(self):
        """Copy barcode dicts from self.input_object to simplify calling dicts"""
        self.demultiplex_barcode_dicts = [dict(x.barcode_dict) for x in self.input_object.barcodes]
        self.hamming_dicts = [dict(x.hamming_dict) for x in self.input_object.barcodes]

    def close_output_object(self):
        """Close output objects"""
        for _, outs in self.output_objects.items():
            # iterate through list of writers
            for out in outs:
                out.close()

    def get_seq_barcodes(self, line):
        """Retrieve barcode sequences from barcode files based on index"""
        return [line[index][8] for index in self.barcode_indices]

    def get_seq_reads(self, line):
        """Retrieve read sequences from read files based on index"""
        return [line[index] for index in self.read_indices]

    def set_new_barcode(self, barcode, index):
        """Add barcode to seq barcode dict, set hashing ID based on hamming distance threshold and if hamming distance
        value is unique
        Arugments:
            barcode (str): sequencing barcode observed
            index (int): index of barcode list, which index is the barcode for
        Attributes:
            HammingDistance (class): calculates hamming distance between two strings,
            bool to adjust for length difference
        Returns:
            hamming_id (str): used to form sample_id to pull output object
            """
        # store hamming distance results
        hamming_list = []
        # iterate over all reference sequences
        for ref_barcode, sample_id in self.hamming_dicts[index].items():
            ham_score = HammingDistance(seq_barcode=barcode,
                                        reference_barcode=ref_barcode,
                                        mixed_length=self.mixed_length).hamming
            hamming_list.append((ham_score, sample_id))
        hamming_list.sort()
        # select sample_id with lowest hamming distance, check to make sure value is unique
        if hamming_list[0][0] != hamming_list[1][0] and hamming_list[0][0] <= self.ham_dist:
            hamming_id = hamming_list[0][1]
            # add barcode to dict matched to sample_id hamming_id
            self.demultiplex_barcode_dicts[index][barcode] = hamming_id
        else:
            # if hamming distance fails checks set sample match to unmatched
            hamming_id = 'unmatched'
            self.demultiplex_barcode_dicts[index][barcode] = hamming_id
        return hamming_id

    def output_read(self, read_id, seq_reads):
        """ Retrieve outpupt object based on read_id and ouput formatted line"""
        try:
            outs = self.output_objects[self.input_object.sample_ids[read_id]]
        except KeyError:
            outs = self.output_objects['unmatched']
        for read, out in zip(seq_reads, outs):
            # format and write line
            out.write(qseq_fastq_conversion(read))

    def iterate_through_qseq(self):
        """Iterate through qseq files and write out lines to fastq sample files
        Attributes:
            self.input_object.file_list
            """
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
