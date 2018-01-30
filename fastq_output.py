import os
import time


class FastqOut:

    def __init__(self, input_list):
        self.input_queue = input_list[0]
        self.sample_list = input_list[1]
        self.output_directory = input_list[2]
        self.read_count = input_list[3]
        self.output_dict = {}
        self.run()

    def run(self):
        self.output_objects()
        self.output()

    def output_objects(self):
        """Initialized objects to output reads in fastq format, will generate a file for every 'read' labeled file in
        the file_label plus a file for unmatched reads
        -----------------------------------------------------
        output_directory = path to write files; folder must already exist
        self.sample_list: list of input samples
        self.read_count: number of read files labeled in file label
        returns self.output_dict; hashes to output object based on sample name"""
        # initialize output objects for all samples
        for sample in self.sample_list:
            object_list = []
            for count in range(self.read_count):
                object_list.append(open(self.output_directory + sample + '_' +
                                        str(count + 1) + '.fastq', 'w'))
            self.output_dict[sample] = object_list
        object_list = []
        # initialize output objects for unmatched reads
        for count in range(self.read_count):
            if os.path.exists(self.output_directory + 'unmatched' + '_' + str(count + 1) + '.fastq'):
                object_list.append(open(self.output_directory + 'unmatched' + '_' + str(count + 1) + '.fastq', 'a'))
            else:
                object_list.append(open(self.output_directory + 'unmatched' + '_' +
                                        str(count + 1) + '.fastq', 'w'))
        self.output_dict['unmatched'] = object_list

    def output(self):
        time.sleep(1)
        while not self.input_queue.empty():
            sample_fastq = self.input_queue.get()
            output = self.output_dict[sample_fastq[0]]
            for sample in zip(sample_fastq[1], output):
                sample[1].write(sample[0])
        for output_object in self.output_dict.values():
            for out in output_object:
                out.close()
