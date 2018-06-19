#!/usr/bin/env python3


def duplicates(lst, item):
        """Python index lookup in a list returns the first index by default, this function returns all indexes of an item
        in a list. Loops through list and returns count if item = item in list.
         Arguments:
           lst (list): list of strings
           item (str): object to return indices
        Returns:
            list with indices of object
            """
        return [i for i, x in enumerate(lst) if x == item]


def qseq_fastq_conversion(qseq_list):
    """Convert a .qseq file to .fastq before output
    Arguments:
        qseq_list (list): '\t' split qseq string
    Returns:
       formatted string in fastq format
       """
    # pull fastq header information
    fastq_id = '@%s:%s:%s:%s:%s:%s#%s/%s' % (qseq_list[0], qseq_list[1], qseq_list[2], qseq_list[3], qseq_list[4],
                                             qseq_list[5], qseq_list[6], qseq_list[7])
    # mask any missing base calls
    seq = qseq_list[8].replace('.', 'N')
    # line 3 is used by some tools to store additional information, so the output will be static
    line_3 = '+'
    # The quality output is a single value in qseq list
    quality = qseq_list[9]
    fastq_out = '%s\n%s\n%s\n%s\n' % (fastq_id, seq, line_3, quality)
    return fastq_out
