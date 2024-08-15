# Nanopore Data Assembly Tutorial

This tutorial will guide you through the process of filtering, assembling, circularizing, and annotating nanopore sequencing data using various bioinformatics tools. The tools used in this tutorial are:

- **NanoFilt** for filtering low-quality and short reads.
- **Flye** for de novo assembly of long reads.
- **DNAapler** for reorienting the circularized assembly.
- **Prokka** for annotating the assembled genome.

## Prerequisites

Before starting, ensure that you have the following tools installed:

- NanoFilt
- Flye
- DNAapler
- Prokka

Additionally, ensure that you have your sequencing data in `.fastq.gz` format.

## Step 1: Filter Low-Quality and Short Reads

Reads under 1kb in length can introduce noise in the assembly process. We will filter out reads shorter than 1000 bp and with a quality score lower than 10.

```bash
mkdir -p results/filtered
gunzip -c data/P2-Rapid010_Proferment_barcode09.fastq.gz | NanoFilt -l 1000 -q 10 > results/filtered/P2-Rapid010_Proferment_barcode09_filtered.fastq
```

