import os


class FastqOut:
    """Initialized objects to output reads in fastq format, will generate a file for every 'read' labeled file in
        the file_label plus a file for unmatched reads
    Arguments:
        sample_dict (dict): dict of sample_id's
        read_count (int): number of read files to output
        output_directory (str): path to directory to write out file
    Attributes:
        self.output_dict (dict): dict hashing sample_ids to output objects
        self.output_objects (func): set output_dict
    """

    def __init__(self, sample_dict=None, read_count=None, output_directory=None):
        self.output_dict = {}
        self.output_objects(sample_dict=sample_dict, read_count=read_count, output_directory=output_directory)

    def output_objects(self, sample_dict=None, read_count=None, output_directory=None):
        # initialize output objects for all samples
        for _, sample in sample_dict.items():
            object_list = []
            for count in range(read_count):
                output_path = '%s%s_%s.fastq' % (output_directory, sample, str(count + 1))
                object_list.append(open(output_path, 'w'))
            self.output_dict[sample] = object_list
        # initialize output objects for unmatched reads
        for count in range(read_count):
            output_path = '%sunmatched_%s.fastq' % (output_directory, str(count + 1))
            if os.path.exists(output_path):
                self.output_dict['unmatched'] = [open(output_path, 'a')]
            else:
                self.output_dict['unmatched'] = [open(output_path, 'w')]
