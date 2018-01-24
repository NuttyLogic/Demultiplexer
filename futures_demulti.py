#!/usr/bin/env python3

from MultipleIterator import MultipleSequencingFileIterator
import os
from multiprocessing import Pool
import time
from collections import namedtuple


def duplicates(lst, item):
    """Python index lookup in a list returns the first index by default, this function returns all indexes of an item
    in a list. Loops through list and returns count if item = item in list.
       -----------------------------------------------------
       lst=List
       item=object in list
       returns; list with indices of object"""
    return [i for i, x in enumerate(lst) if x == item]


def qseq_fastq_conversion(qseq_list):
    """Convert a .qseq file to .fastq before output
       -----------------------------------------------------
       qseq_list: '\t' split qseq string
       returns; an output string with four new line indicators in fastq format"""
    # pull fastq header information
    fastq_id = '@%s:%s:%s:%s:%s:%s#%s/%s' % (qseq_list[0], qseq_list[1], qseq_list[2], qseq_list[3], qseq_list[4],
                                          qseq_list[5], qseq_list[6], qseq_list[7])
    # mask any missing base calls
    seq = qseq_list[8].replace('.', 'N')
    # line 3 is used by some tools to store additional information, so the output will be static
    line_3 = '+'
    # The quality output is a single value in qseq list
    quality = qseq_list[9]
    fastq_out = fastq_id + '\n' + seq + '\n' + line_3 + '\n' + quality + '\n'
    return fastq_out


def output_objects(sample_list='path', output_directory='path', read_count=None):
    """Initialized objects to output reads in fastq format, will generate a file for every 'read' labeled file in
    the file_label plus a file for unmatched reads
    -----------------------------------------------------
    output_directory = path to write files; folder must already exist
    self.sample_list: list of input samples
    self.read_count: number of read files labeled in file label
    returns self.output_dict; hashes to output object based on sample name"""
    # initialize output objects for all samples
    output_dict = {}
    for sample in sample_list:
        object_list = []
        for count in range(read_count):
            object_list.append(open(output_directory + sample + '_' +
                                    str(count + 1) + '.fastq', 'w'))
        output_dict[sample] = object_list
    object_list = []
    # initialize output objects for unmatched reads
    for count in range(read_count):
        if os.path.exists(output_directory + 'unmatched' + '_' + str(count + 1) + '.fastq'):
            object_list.append(open(output_directory + 'unmatched' + '_' + str(count + 1) + '.fastq', 'a'))
        else:
            object_list.append(open(output_directory + 'unmatched' + '_' +
                                    str(count + 1) + '.fastq', 'w'))
    output_dict['unmatched'] = object_list
    return output_dict


def set_barcode_offset(input_list, class_instance='test', gnu_zipped=False):
    """Iterate through the first 1,000,000 lines of the first qseq set to set the barcode length, if the barcode
    and read length are equal the fucntion exits
    ---------------------------------------------
    sets class_instance (demultiplex_input instance) self.left_offset/right_offset"""
    iterator = MultipleSequencingFileIterator(input_list, directory=class_instance.directory, gnu_zipped=gnu_zipped)
    # get position of barcode files
    barcode_indexes = duplicates(class_instance.file_label, 'barcode')
    # get position of read files
    read_indexes = duplicates(class_instance.file_label, 'read')
    # set barcode list, for looping
    barcode_list = [class_instance.barcode_1, class_instance.barcode_2]
    # loop through grouped files
    left_offset = 0
    right_offset = int(class_instance.right_offset)
    barcode_offset_set = False
    offset_list = None
    for count, line in enumerate(iterator.iterator_zip()):
        if count < 100000:
            # set string with Illumina quality control information
            combined_filter = ''.join([qual[-1] for qual in line])
            # initialize empty sample_id value
            sample_id = ''
            # if all reads don't pass filter don't consider
            if '0' not in combined_filter:
                # loop through barcode_indexes, get sample key
                for index_count, index in enumerate(barcode_indexes):
                    # get sequence location in qseq file, and set barcode offset if needed
                    barcode_read = line[barcode_indexes[index_count]][8]
                    if class_instance.right_offset != len(barcode_read):
                        if class_instance.right_offset  > len(barcode_read):
                            print('Barcode length greater than read length')
                            os.exit()
                        elif class_instance.right_offset < len(barcode_read):
                            offset = len(barcode_read) - class_instance.right_offset
                            if not offset_list:
                                offset_list = [0 for _ in range(offset + 1)]
                            for offset_count, offset_range in enumerate(range(offset + 1)):
                                try:
                                    (barcode_list[index_count]
                                     [line[barcode_indexes[index_count]][8]
                                        [offset_range:(class_instance.right_offset + offset_range)]])
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
        right_offset = class_instance.right_offset + left_offset
    class_instance.left_offset = left_offset
    class_instance.right_offset = right_offset


def iterate_through_qseq(workers=8, demultiplex_instance='input_class', output_dictionary='out_dict',
                         gnu_zipped=False):
    """function to wrap process_qseq and set_barcode_offset
    -------------------------------------------------------
    inputs: ProcessDemultiplexInput.instance, Output_dictionary, gnu_zipped status
    returns: NameTuple with read value"""
    # transpose iterator list
    iterator_file_list = list(map(list, zip(*demultiplex_instance.file_list)))
    file_metrics = namedtuple('file_set_metrics', ('reads', 'indexed_reads', 'unmatched_reads', 'reads_pass_filter'))
    metrics = file_metrics([], [], [], [])
    # loop through lists of files
    pool = Pool(processes=workers)
    jobs = []
    set_barcode_offset(iterator_file_list[0], class_instance=demultiplex_instance, gnu_zipped=gnu_zipped)
    for iterator_file in iterator_file_list:
        proc = pool.apply_async(func=process_qseq(iterator_file, output_dictionary=output_dictionary,
                                                  class_instance=demultiplex_instance, shared_metrics=metrics,
                                                  gnu_zipped=gnu_zipped))
        jobs.append(proc)
    # Wait for jobs to complete before exiting
    while not all([p.ready() for p in jobs]):
        time.sleep(5)
    # Safely terminate the pool
    pool.close()
    pool.join()
    for sample in output_dictionary.values():
        for out_object in sample:
            out_object.close()
    return metrics


def process_qseq(input_list, class_instance='test', shared_metrics='named_tuple', output_dictionary='test',
                 gnu_zipped=False):
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
    reads = 0
    indexed_reads = 0
    unmatched_read = 0
    reads_pass_filter = 0
    iterator = MultipleSequencingFileIterator(input_list, directory=class_instance.directory, gnu_zipped=gnu_zipped)
    # get position of barcode files
    barcode_indexes = duplicates(class_instance.file_label, 'barcode')
    # get position of read files
    read_indexes = duplicates(class_instance.file_label, 'read')
    # set barcode list, for looping
    barcode_list = [class_instance.barcode_1, class_instance.barcode_2]
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
                    key = barcode_list[index_count][line[barcode_indexes[index_count]][8][class_instance.left_offset:class_instance.right_offset]]
                except KeyError:
                    # if barcode sequence not in barcode dictionary set key to 'x'
                    key = 'x'
                sample_id = '{0}key{1}'.format(sample_id, str(key))
            # if barcode matches with key proceed
            if 'x' not in sample_id:
                try:
                    # look up sample, if matched set sample name
                    sample = class_instance.sample_key[sample_id]
                    indexed_reads += 1
                except KeyError:
                    # if sample unmatched write to unmatched reads
                    unmatched_read += 1
                    sample = 'unmatched'
                # retrieve list of output objects
                out = output_dictionary[sample]
                # write line to file
                for out_count, output_object in enumerate(out):
                    # convert qseq line to fatq format
                    output_object.write(qseq_fastq_conversion(line[read_indexes[out_count]]))
            else:
                # if barcode sequence not in dictionary write to unmatched
                unmatched_read += 1
                sample = 'unmatched'
                out = output_dictionary[sample]
                for out_count, output_object in enumerate(out):
                    output_object.write(qseq_fastq_conversion(line[read_indexes[out_count]]))
    shared_metrics.reads.append(reads)
    shared_metrics.unmatched_reads.append(unmatched_read)
    shared_metrics.indexed_reads.append(indexed_reads)
    shared_metrics.reads_pass_filter.append(reads_pass_filter)
