# Demultiplexer

Script to demultiplex .qseq files to .fastq files. The current version supports single and dual index barcodes. The script is hash based,
and calculates 3 mismatches by defualt.  

## Usage

```python
python3 Demultiplex -W cores -D directory -S sample_key -B1 barcode_1 -B2 barcode_2 -L file_labels -O output_directory -I input_file_1 input_file_2 ...
```

## Inputs

- -D, /path/ to qseq directory
- -W, Number of cores available, default = 2
- -S, /path/sample_file.txt; file should be formatted as 'barcode tab sample_name' for single index and 'barcode tab barcode tab sample_name for dual indexes, see ~/tests/single_index_test or ~/tests/dual_index_test for an example
- -B1, /path/barcode_1_file, line separated list of barcodes
- -B2, /path/barcode_2_file, line separated list of barcodes
- -B1R, Consider Barcode1 Reverse Complement
- -B2R, Consider Barcode2 Reverse Complement
- -L, string of r and b character to designate input files as barcode or read files, should be the same order as input files
- -O, path to output directory
- -Z, designate is inpute qseq files are gzipped, slows processing
- -I, qseq file prefix and suffix separated by \^, ie. -I s_1_.\^.qseq.txt s_2_.\^.qseq.txt

## Examples

### Single Index Demultiplex

```python
python3 Demultiplex -D ~/Demultiplexer/tests/test_qseq -W 2 -S ~/Demultiplexer/tests/test_sample_files/single_index_test.txt -B1 ~/Demultiplexer/tests/test_sample_files/N700_nextera_bacrodes.txt -L 'rb' -M 1 -O ~/Demultiplexer/tests/test_output/ -I 1_test.^.qseq.txt 2_test.^.qseq.txt
```
### Dual Index Demultiplex

```python
python3 Demultiplex -D ~/Demultiplexer/tests/test_qseq -W 2 -S ~/Demultiplexer/tests/test_sample_files/single_index_test.txt -B1 ~/Demultiplexer/tests/test_sample_files/N700_nextera_bacrodes.txt -B1R -B2 ~/Demultiplexer/tests/test_sample_files/N500_nextera_bacrodes.txt -B2R -L 'rbbr'  -O ~/Demultiplexer/tests/test_output/ -I 1_test.^.qseq.txt 2_test.^.qseq.txt 3_test.^.qseq.txt 4_test.^.qseq.txt
```

### Multiple Read Files with Single Index

```python
python3	Demultiplex.py	-D	~/tests/test_qseq/	-S ~/tests/test_sample_files/single_index_test.txt	-B1	~/tests/test_sample_files/N700_nextera_barcodes.txt	-W	2	-L	rrb	-O	~/tests/test_output/	-I	1_test.^.qseq.txt	4_test.^.qseq.txt 2_test.^.qseq.txt
```

## Setup/Requirements
- Download a [release](https://github.com/NuttyLogic/Demultiplexer/release), extract, and run. Demultiplex will work with python > 3.4.
