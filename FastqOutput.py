import os


class FastqOut:

    def __init__(self, sample_dict=None, read_count=None, output_directory=None):
        self.output_dict = {}
        self.output_objects(sample_dict=sample_dict, read_count=read_count, output_directory=output_directory)

    def output_objects(self, sample_dict=None, read_count=None, output_directory=None):
        """Initialized objects to output reads in fastq format, will generate a file for every 'read' labeled file in
        the file_label plus a file for unmatched reads
        -----------------------------------------------------
        output_directory = path to write files; folder must already exist
        self.sample_list: list of input samples
        self.read_count: number of read files labeled in file label
        returns self.output_dict; hashes to output object based on sample name"""
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
