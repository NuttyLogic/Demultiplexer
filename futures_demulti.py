#!/usr/bin/env python3

from set_read_offset import BarcodeOffset
from qseq_proccessing import *
from fastq_output import *
import multiprocessing as mp


def iterate_through_qseq(workers=8, demultiplex_instance='input_class', output_directory='path', gnu_zipped=False):
    """function to wrap process_qseq and set_barcode_offset
    -------------------------------------------------------
    inputs: ProcessDemultiplexInput.instance, Output_dictionary, gnu_zipped status
    returns: NameTuple with read value"""
    # transpose iterator list
    iterator_file_list = list(map(list, zip(*demultiplex_instance.file_list)))
    fastq_manager = mp.Manager()
    fastq_queue = fastq_manager.Queue()
    BarcodeOffset(iterator_file_list[0], demultiplex_instance, gnu_zipped)
    read_stats = fastq_manager.list([0, 0, 0, 0])
    input_arguments = (demultiplex_instance.directory,
                       demultiplex_instance.file_label,
                       gnu_zipped,
                       demultiplex_instance.barcode_1,
                       demultiplex_instance.barcode_2,
                       demultiplex_instance.left_offset,
                       demultiplex_instance.right_offset,
                       demultiplex_instance.sample_key,
                       demultiplex_instance.read_count,
                       fastq_queue,
                       read_stats)
    for sample in iterator_file_list:
        sample.append(input_arguments)
    pool = mp.Pool(processes=workers - 1)
    pool.map(ProcessQseq, iterator_file_list)
    # Wait for jobs to complete before exiting
    fastq_output_tupple = (fastq_queue,
                           demultiplex_instance.sample_list,
                           output_directory,
                           demultiplex_instance.read_count)
    p = mp.Process(target=FastqOut, args=(fastq_output_tupple,))
    p.start()
    # Safely terminate the pool
    pool.close()
    pool.join()
    p.join()
    return read_stats
