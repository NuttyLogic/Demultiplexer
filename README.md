# Demultiplexer

Package to demultplex sequencing lane .qseq files to sample .fastq files. 

## Usage

```bash
python3 Demultiplex.py -D directory -S sample_key -H 1 -BR -L file_labels -O output_directory -I input_file_1 input_file_2 ...
```

## Inputs

  - -D            /path/ to qseq directory
  - -S           /path/sample_file.txt file should be formatted as 'barcode
                  tab sample_name' for single index and 'barcode tab barcode
                  tab sample_name' for dual indexes
  - -BR             Consider Barcodes Reverse Complements
  - -L            string of r and b character to designate input files as
                  barcode or read files, should be the same order as inputfile
  - -O            Path to Output Directory
  - -I [I [I ...]]  qseq file prefix and suffix separatedby ^, ie. -I
                  s_1_^.qseq.txt s_2_^.qseq.txt
  - -H             Minimum hamming distance threshold for a sequencing barcode
                  to be considered, default=0
  - -M             If the reference barcodes for an index contain barcodes of
                  different length, return hamming distance plus length difference

## Examples

### Single Index Demultiplex

```bash
python3 Demultiplex.py -D tests/test_qseq/ -S tests/test_sample_files/single_index_test.txt -L rb -O tests/test_output/ -I 1_test.^.qseq.txt.gz 2_test.^.qseq.txt
```
### Dual Index Demultiplex

```bash
python3 Demultiplex.py -D tests/test_qseq/ -S tests/test_sample_files/dual_index_test.txt -L rbbr -O tests/test_output/ -I 1_test.^.qseq.txt.gz 2_test.^.qseq.txt 3_test.^.qseq.txt 4_test.^.qseq.txt
```

### Multiple Read Files with Single Index

```bash
python3 Demultiplex.py -D tests/test_qseq/ -S tests/test_sample_files/dual_index_test.txt -L rbb -O tests/test_output/ -I 1_test.^.qseq.txt.gz 2_test.^.qseq.txt 3_test.^.qseq.txt 
```

## Sample Key Formatting Example
- Single Index Sample Key
```
TAAGGCGA	sample_1
CGTACTAG	sample_2
AGGCAGAA	sample_3
TCCTGAGC	sample_4
GGACTCCT	sample_5
TAGGCATG	sample_6
```
- Dual Index Sample Key
```
TAAGGCGA	TAGATCGC	sample_1
CGTACTAG	TAGATCGC	sample_2
AGGCAGAA	TAGATCGC	sample_3
TCCTGAGC	TAGATCGC	sample_4
GGACTCCT	TAGATCGC	sample_5
TAGGCATG	TAGATCGC	sample_6
```
## Setup/Requirements
- Download a [release](https://github.com/NuttyLogic/Demultiplexer/release), extract, and run. Demultiplex will work with python > 3.4.
