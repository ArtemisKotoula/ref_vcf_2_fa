# ref_vcf_2_fa - CDS Reconstructor from SNP Variants

A Perl script for reconstructing sample-specific coding sequences (CDS) by introducing SNPs from a VCF-like variant file into reference sequences.

## Overview

This script reconstructs haplotype-specific coding sequences for multiple samples by applying SNPs from a tab-delimited variant file to reference CDS sequences.

For each sequence, the script:

- Retrieves the reference CDS from a one-line FASTA file.
- Applies SNPs according to each sample's genotype.
- Reconstructs two haplotypes for every sample.
- Writes the reconstructed sequences to a FASTA file.

Each sequence generates one output FASTA file containing all reconstructed haplotypes and the reference sequence.

## Input Files

Examples for the input files can be found in the repository, within the `example` directory.

### 1. Reference FASTA

The reference FASTA must contain one sequence per line (no wrapped sequences).

The transcript IDs in the FASTA headers must match those in the variant file.

### 2. Variant File

The script expects a tab-delimited VCF file.

## Usage

1. Edit the input file paths near the beginning of the script:

>$file  = "vcf.txt";
>
>$fasta = "one_line.fa";

2. Run the script:

> perl ref_vcf_2_fa.pl

## Output

The script creates a `ref_out/` directory and one FASTA file per transcript within it.

