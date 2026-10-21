# demuxolotl
Demultiplex BLink-seq barcodes from raw FASTQ files.

This is a small program, originally part of Harpy, that demultiplexes the linked-read barcodes out of raw BLink-seq FASTQ files.
It's intended to work on data that's already been demultiplexed at the sample level, which should have been done by the sequencing
facility (NovaSeq).

## Install
- requires Go and Make
```bash
git clone https://github.com/BLinkseq/demuxolotl.git
cd demoxolotl
make
```
The binary gets written to `bin/demoxolotl`. You can use it from there or copy/move it to a location in your PATH.


## Run
```bash
demuxolotl <options> PREFIX R1.fq[.gz] R2.fq[.gz] > stats.txt
```

## Arguments
| Name | What it is |
|:------------|:----|
| `PREFIX` | The output filename prefix. Writes PREFIX.R1.fq.gz and PREFIX.R2.fq.gz |
| `R1.fq.gz` | Input R1 (forward) FASTQ file. Should be fresh from the sequencer, no trimming. |
| `R2.fq.gz` | Input R2 (reverse) FASTQ file. Should be fresh from the sequencer, no trimming. |
| `stats.txt` | Demultiplexing stats that get written to `stdout` |

### Output stats
A few output diagnostics are provided:
```
Total read pairs processed:          5000
Total unique barcodes:               4815
Total reads with barcode failure:    64
Total reads with valid barcodes:     9830
Total reads with corrected barcodes: 334
Discarded (ME sequence not found):   53
```

## Options
| Option | Short Option | Default | Description |
|:------------|:---:|-----:|:------|
| `--threads` | `-t` | 4 | Threads to use |
| `--me` | `-m` | `AGATGTGTATAAGAGACAG` | ME sequence to search for |
| `--max-mismatch` | `-x` | 2 | Maximum allowable ME mismatch score |
| `--qx-rx` | `-q` |  | Output tags QX:Z (barcode quality) and RX:Z (corrected barcode) |
| `--schema` | `-s` | `use embedded` | [Pheniqs-style JSON schema](demux/schema.json) file with the barcodes and decoder settings |

**TL;DR**: In most cases, you will only be touching the `--threads` options. Don't modify the other options unless you have a really good reason to do so.

### ME Sequence
The BLink-seq chemistry requires ME sequences and the demultiplexing identifies the ME sequence as an anchor point to pad the variable-end stagger at the end to a fixed length. Unless you were given explicit instructions to use a different ME sequence, do not modify this value.

### Max Mismatch
When looking for the ME sequence, this is the maximum mismatch score you're willing to tolerate before deciding it's a bad match. Each ATCG nucleotide mismatch counts as 1 point, whereas N nucleotides in the read sequence count as 0.3 points. Since the ME sequence is 19 bases long, it's really not recommended to increase this value above the default of `2`.

### QX and RX tags
The `RX:Z` tag records the _corrected_ nucleotide barcode. The original observed barcode is only given if an unrecognized barcode couldn't be corrected. The `QX:Z` tag is the PHRED score of the barcode. Both tags are optional and aren't used by any downstream software we know if-- they are for diagnostic purposes.

### Schema
This program ports over the PAMLD algorithm and decoding logic from [Pheniqs](https://github.com/biosails/pheniqs), stripping all the extra features BLink-seq doesn't need. Given the Pheniqs origins, we continue to use a pheniqs configuration file for decoding barcodes. If using the default BLink-seq chemistry, you will not need to include your own schema file.
