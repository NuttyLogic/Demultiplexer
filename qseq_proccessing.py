from MultipleIterator import *
from sequence_helpers import *


class ProcessQseq:

    def __init__(self, input_list):
        self.input_files = input_list[:-1]
        class_arguments = input_list[-1]
        self.directory = class_arguments[0]
        self.file_label = class_arguments[1]
        self.gnu_zipped = class_arguments[2]
        self.barcode_1 = class_arguments[3]
        self.barcode_2 = class_arguments[4]
        self.left_offset = class_arguments[5]
        self.right_offset = class_arguments[6]
        self.sample_key = class_arguments[7]
        self.read_count = class_arguments[8]
        self.queue = class_arguments[9]
        self.metrics = class_arguments[10]
        self.run()

    def run(self):
        self.process_qseq()

    def process_qseq(self):
        """Iterate through groups of ID'd qseq files and output demultiplexed fastq files
        -----------------------------------------------------
        self.file_list = sorted list of input files;
        self.barcode*: dictionary of barcodes
        self.file_label: list of barcode/ read file positions
        self.output_dict: dict of output objects, 1 per read qseq files
        returns;
        self.read = number of reads processes across all of the groups of files
        self.reads_pass_filter = number of reads that pass the Illumina filter
        self.indexed_reads = reads matched to sample index
        self.unmatched_reads = number of unmatched reads
        """
        # initialize iterator object for sorted group of files
        # directory, gnu_zipped, file_label, barcode_1, bacode_2, left_offset, right_offset, sample_key
        reads = 0
        indexed_reads = 0
        unmatched_read = 0
        reads_pass_filter = 0
        iterator = MultipleSequencingFileIterator(self.input_files, directory=self.directory,
                                                  gnu_zipped=self.gnu_zipped)
        # get position of barcode files
        barcode_indexes = duplicates(self.file_label, 'barcode')
        # get position of read files
        read_indexes = duplicates(self.file_label, 'read')
        # set barcode list, for looping
        barcode_list = [self.barcode_1, self.barcode_2]
        # loop through grouped files
        for count, line in enumerate(iterator.iterator_zip()):
            reads += 1
            # set string with Illumina quality control information
            combined_filter = ''.join([qual[-1] for qual in line])
            # initialize empty sample_id value
            sample_id = ''
            # if all reads don't pass filter don't consider
            if '0' not in combined_filter:
                reads_pass_filter += 1
                # loop through barcode_indexes, get sample key
                for index_count, index in enumerate(barcode_indexes):
                    try:
                        # get sequence location in qseq file, and set barcode offset if needed
                        key = barcode_list[index_count][line[barcode_indexes[index_count]][8][self.left_offset:self.right_offset]]
                    except KeyError:
                        # if barcode sequence not in barcode dictionary set key to 'x'
                        key = 'x'
                    sample_id = '{0}key{1}'.format(sample_id, str(key))
                # if barcode matches with key proceed
                if 'x' not in sample_id:
                    try:
                        # look up sample, if matched set sample name
                        sample = self.sample_key[sample_id]
                        indexed_reads += 1
                    except KeyError:
                        # if sample unmatched write to unmatched reads
                        unmatched_read += 1
                        sample = 'unmatched'
                    # retrieve list of output objects
                    out = [sample, []]
                    # write line to file
                    for out_count in range(self.read_count):
                        # convert qseq line to fatq format
                        out[1].append(qseq_fastq_conversion(line[read_indexes[out_count]]))

                else:
                    # if barcode sequence not in dictionary write to unmatched
                    unmatched_read += 1
                    sample = 'unmatched'
                    out = [sample, []]
                    # write line to file
                    for out_count in range(self.read_count):
                        # convert qseq line to fatq format
                        out[1].append(qseq_fastq_conversion(line[read_indexes[out_count]]))
                self.queue.put(out)
        self.metrics[0] += reads
        self.metrics[1] += unmatched_read
        self.metrics[2] += indexed_reads
        self.metrics[3] += reads_pass_filter
