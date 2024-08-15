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

The first step, even before processing any data is to prepare the working environment. In bioinformatics, an organized workspace is vital, so when you come after some time to your project, you can find and understand whant you were doing, rather thatn spend hours searching through weirldy named directories. It is inportant to always create three directories:

- scripts: all the scripts will be stored here, with meaningful names
- data: Raw data will go in here and, if you want and need, databases
- results: Create a sub directory for every different process you do. If you run a process multiple times with different parameters, include them in the directory name, so you will differenciate them in the future.

```bash
mkdir scripts
mkdir data
mkdir results
```

The dta that we will use in this tutorial is publicly available in NCBI


## Step 1: Filter Low-Quality and Short Reads

Filtering out short and low-quality reads is a crucial step in the assembly process. Reads that are under 1,000 base pairs in length can introduce noise, leading to a less accurate assembly. Additionally, reads with a quality score lower than 10 may contain sequencing errors, which can further compromise the accuracy of the assembly.

To achieve this, we use NanoFilt to filter the reads. However, NanoFilt requires the input file to be in an uncompressed format. Since our sequencing data is in a compressed .fastq.gz file, we first need to unzip the file. We use the gunzip -c command to decompress the file while streaming its content directly to NanoFilt. This approach ensures that the original compressed file remains intact, as the -c option tells gunzip to write the output to standard output rather than replacing the original file. NanoFilt then processes the uncompressed data, filtering out any reads shorter than 1,000 base pairs or with a quality score below 10, resulting in a cleaner dataset for downstream analysis.

```bash
# Create target directory
mkdir -p results/filtered
# Run the filtering
gunzip -c data/P2-Rapid010_Proferment_barcode09.fastq.gz | NanoFilt -l 1000 -q 10 > results/filtered/P2-Rapid010_Proferment_barcode09_filtered.fastq
```

 ## Step 2: Assemble the reads

We will use Flue, a *de novo* assembler for long reads. It is designed for a wide range of datasets, and it has several parameters that will have to specify to obtain the most optimal assembly. First of all, we have to indicate the type of input we are using, `nano-raw` in our case, since we have regular uncorrected nanopore reads. Also, it is important to specify genome size, so flye know's what to expect and can estimate the correct depth of sequencing and act accordingly.

```bash 
flye --nano-raw results/filtered/P2-Rapid010_Proferment_barcode09_filtered.fastq \
    -o results/assembly_brk09_basic \
    --genome-size 430000
```

There is another useful option in Flye, which is `asm-coverage`, where you can indicate the desired covereage you want in your assembly and it will automatically subset your input to meet the requirement. This might be useful in cases where there is too much information (over 100x), which will slow down the process and might lead Flye to produce errors.

Flye produces as ouptut multiple files, such as the final assembly (assembly.fasta) and the assembly graphs. The assembly graphs are indicators of how good the assembly was done, and if it went accordingly. The graphs can be visualized using Bandage, a software tool. If we are assembling an isolate, ideally we would like to see one big circular fragment, the chromosome, with some smaller circular or linear fragments (plasmids)

![Successful assembly with one circular chromosome and one plasmid](Images/Assembly_graph_good.png)

![Failed assembly, since no circular chromosome has been found, only multiple contigs](Images/Assembly_graph_bad.png)
 
 ## Step 3: Circularize the assembly

Once the assembly is done, the first contig will contain our chromosome. By consensus, the bacterial chromosomes have to be oriented so the first gene present is the origin of replication (dnaA in most cases). That way, different assemblies of the same genome can be compared between them (otherwise the genes would have different positons!). For that purpose, we use a circularizer tool, such as dnaapler. It will search for any origin of replication genes of both genomes and plasmids and then alter the contigs so these are the first positions.

```bash
dnaapler all -i results/assembly_brk09_basic/assembly.fasta \
    -o results/assembly_brk09_basic/circularization \
    -t 15

mv results/assembly_brk09_basic/circularization/dnaapler_reoriented.fasta results/assembly_brk09_basic/assembly_circularized.fasta
```

## Step 4: Annotate the assembly

The final step in this pipeline is the annotation. This process consists of two disticnt phases: calling the Open Reading Frames (ORFs) and then assigning a function to them. An ORF is found by determine regions in the genome that begin with an start codon and are closed with a stop codon. Each of this sequences will be putative proteins, which are compared against a gene database to assign known fucntions to each one of them. Depending on the database used, we will recover some genes and miss others, so it is a good practice to try and run our annotation on a more specific database if we are looking for a specific function (such as a CAZy database or ncRNA).

```bash
prokka --outdir results/annotation \
    --prefix brk09 \
    --genus Bacillus results/assembly_brk09_basic/assembly_circularized.fasta
```