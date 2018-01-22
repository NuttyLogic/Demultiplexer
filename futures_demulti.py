#!/usr/bin/env python3

from MultipleIterator import MultipleSequencingFileIterator
from demultiplex_input_process import ProcessDemultiplexInput
import os
import sys
from concurrent import futures
import threading as th
from multiprocessing import Process
from multiprocessing import Pool
import time


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
    fastq_id = '@%s:%s:%s:%s:%s#%s/%s' % (qseq_list[0], qseq_list[2], qseq_list[3], qseq_list[4],
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


def iterate_through_qseq(file_list=None, workers=8):
    # transpose iterator list
    iterator_file_list = list(map(list, zip(file_list)))
    # loop through lists of files
    start = time.time()
    pool = Pool(processes=workers)
    jobs = []
    print(iterator_file_list)
    for iterator_file in iterator_file_list:
        qseq_input = '?'.join(iterator_file[0])
        print(qseq_input)
        proc = pool.apply_async(func=process_qseq(qseq_input), args=(qseq_input,))
        jobs.append(proc)
    # Wait for jobs to complete before exiting
    while (not all([p.ready() for p in jobs])):
        time.sleep(5)
    # Safely terminate the pool
    pool.close()
    pool.join()
    end = time.time()
    print(end - start)
    for sample in output_dictionary.values():
        for out_object in sample:
            out_object.close()


def process_qseq(input_list):
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
    iterator = MultipleSequencingFileIterator(input_list, directory=directory, gnu_zipped=False)
    # get position of barcode files
    barcode_indexes = duplicates(file_label, 'barcode')
    # get position of read files
    read_indexes = duplicates(file_label, 'read')
    # set barcode list, for looping
    barcode_list = [barcode_1, barcode_2]
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
                    # get sequence location in qseq file
                    key = barcode_list[index_count][line[barcode_indexes[index_count]][8]]
                except KeyError:
                    # if barcode sequence not in barcode dictionary set key to 'x'
                    key = 'x'
                sample_id = '{0}key{1}'.format(sample_id, str(key))
            # if barcode matches with key proceed
            if 'x' not in sample_id:
                try:
                    # look up sample, if matched set sample name
                    sample = sample_key[sample_id]
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
    print([reads, unmatched_read, indexed_reads, reads_pass_filter])
    return [reads, unmatched_read, indexed_reads, reads_pass_filter]

def test_func(input_file, output_dict):
    print(output_dict)
    print('ahhhhh')
    print(input_file)

test = ProcessDemultiplexInput('1_test.^.qseq.txt', '2_test.^.qseq.txt', directory='tests/test_qseq/',
sample_key='tests/test_sample_files/single_index_test.txt', barcode_1='tests/test_sample_files/N700_nextera_barcodes.txt',
file_label='rb')
test.run()
output_run = output_objects(sample_list=test.sample_list, output_directory='tests/test_output/', read_count=test.read_count)

output_dictionary=output_run
barcode_1=test.barcode_1
barcode_2=test.barcode_2
sample_key=test.sample_key
file_label=test.file_label
directory=test.directory
gnu_zipped=False

iterate_through_qseq(file_list=test.file_list)
