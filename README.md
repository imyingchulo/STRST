# STRST: STAR and Trinity RNA-seq Tool

Pipelines for RNA-seq analysis: **STAR** and **Trinity**.
This package is made easy to use two approaches for RNA-seq:
a) STAR (alignment) and b) Trinity (_de novel_ assembly)

The package includes the data cleaning, alignment/assembly, DE/Annotation analysis and data visualization for GO.
please find more detail for STAR and Trinity-wiki in the links below:
STAR: https://github.com/alexdobin/STAR
Trinity: https://github.com/trinityrnaseq/trinityrnaseq/wiki

p.s. For **GO** visualization script is not included in the package. But the template is provide in the readme.

## Environment setup
It is better to run this package in virtual environment.
Set up the virtual environment with following steps:

0. Install virtualenv

```shell
$ pip3 install virtualenv 
```

1. Setup virtualenv

```shell
$ python3 -m venv venv 
```

2. Activate virtualenv

```shell
$ source ./venv/bin/activate
```

3. Install package

```shell
$ pip install -e .
```

## Usage guidance

:octocat: Functions include in the package:

:memo: Explanation of the function - `function name`

1. Pre-alignment quality control - `fastqc`
2. Trimmed adapter - `trimmed_adapter`
3. Prepare genome for STAR alignment - `star_genome_prep`
4. STAR alignment - `star_alignment` 
5. Post STAR alignment - `bamfiles_postqc`
6. Prepare files for post STAR alignment - `qualimap_qc`
7. Count the Transcripts - `salmon_quantification` 
8. Count the Transcripts - `htseq_count` 
9. Count the Transcripts - `rsem_count` 
10. Trinity assembly - `trinity_assemble`
11. Estimate the expression after Trinity assembly - `trinity_estimate_expression`
12. Calculate the transcripts and transfer into a matrix - `trinity_abundance_estimates_to_matrix`
13. Trinity DE analysis - `trinity_run_de_analysis `
14. RNA-seq data visualization - `trinity_run_glimma`

## STEP1: DATA QC and PREPARE INPUTS
###  Pre-alignment quality control (fastqc)
For both the STAR and Trinity pipelines, the first step is always to check the raw data quality.

:mag_right: Generate fasta QC report.
```
strst fastqc --fastqc [put the information of "which fastqc"] --data [an absolute directory to rna-seq data'] --outputfolder [output dir]
```

###  Trimmed adapter
The adapters information can find from `fastqc` report.

:mag_right: Trimmed adapter by cutadapt.
```
strst trimmed-adapter --cutadapt [indicate the path where cutadapt installed] --adapter_fwd [adapter forward sequence] --adapter_rev [adapter reverse sequence] --data [fasta files] --extention [fastq.gz or fastq] 
```

## STEP2: SEQUENCE ALIGNMENT or ASSEMBLY
### Using SATR for RNA-seq alignment:

1. Prepare genome reference
```
strst star_genome_prep --star [indicate the path where star installed] --gtf_dir [absolute directory to gft file; example: "/home/ref/gft_dir"] --gtf_file [full path to gft file; example: "/home/ref/gft_dir/sample.gtf"] --ref_fasta [full path to reference fasta file; example: "/home/fasta_dir/sample.fasta"]
```
2. Run STAR alignment
```
strst star_alignment --star [indicate the path where star installed] --data [absolute directory to gft file; example: "/home/ref/gft_dir"] --gtf_dir [absolute directory to gft file; example: "/home/ref/gft_dir"] --gtf_file [full path to gft file; example: "/home/ref/gft_dir/sample.gtf"] --sep [charaters after R1/R2, i.e., file named: 1-G1_HNFVVDRXX_L1_R1.fastq.gz, put "."; file named: 1-G1_HNFVVDRXX_L1_R1.trimmed.fastq.gz, put ".trimmed] --gz [if RNA-seq data is .gz, put "true", if not put "false"] --threads [number of threads, default=15]
```
3. QC alignment output
```
strst bamfiles_postqc --samtools [indicate the path where samtools installed] --bam_dir [absolute path to bam files (the files after STAR), default="result/alignment/"]
```
```
strst qualimap_qc --qualimap [indicate the path where qualimap installed] --bam_dir [absolute path to bam files (the files after STAR), default="result/alignment/"] --gff [An absolute directory to gff or gtf] --output_dir [output folder; i.e. /home/output]
```
4. Quantification ( Three methods are provide in this package: Salmon, htseq and RSEM)
```
strst salmon_quantification --salmon [indicate the path where salmon installed] --fasta [absolute path to fasta] --output_dir [the dir for output file] --bam_dir [absolute path to bam files]
```
```
strst htseq_count --htseq [indicate the path where htseq installed] --gtf_file [absolute path to gtf_file] --bam_dir [absolute path to bam files] --output_dir [output dir, default="result/gene_counts"]
```
```
strst rsem_count --rsem [indicate the path where rsem installed] --fasta_file [absolute path to fasta] --bam_dir [absolute path to bam files] --output_dir [output dir, default="result/gene_counts"]
```

### Using Trinity for _de novel_ assembly RNA-seq:
To obtain more information about Trynity, please visit the developer's website: https://github.com/trinityrnaseq/trinityrnaseq/wiki

In strst package followed the steps from the website provide above:
1. Assembly
```
strst trinity_assemble --trinity [indicate the path where trinity installed] --data [absolute path to RNA-seq data] --sample_list [sample_list; default="/home/ubuntu/sample_input.txt"] --sep [charaters after R1/R2, i.e., file named: 1-G1_HNFVVDRXX_L1_R1.fastq.gz, put "."; file named: 1-G1_HNFVVDRXX_L1_R1.trimmed.fastq.gz, put ".trimmed."] --cpu [request cpu] --memory [request memory] --output_prefix [output dir, default="result/trinity"]
```
2. Estimate expression
```
strst trinity_estimate_expression --trinity [indicate the path where trinity installed] --method [salmon or RSEM] --seqtype [fq or fa] --fasta [absolute path to meta-fasta] --samplelist [An absolute directory to samplelist (generated from trinity_assemble)] --outputdir [Full path to output]
```
3. Abundance estimates to matrix
```
strst trinity_abundance_estimates_to_matrix --trinity_dir [An absolute directory to trinity_dir] --data_dir [Path to access folders] --gene_map [Path to gene_map, includes file name and extension] --output_name [output file name]
```
4. DE analysis
```
strst trinity_run_de_analysis --trinity_dir [An absolute directory to trinity_dir] --method [edgeR or DESeq2] --matrix [Path to matrix, includes file name and extension] --output_name [output file name] --sample_file [An absolute directory to sample_file (generated from trinity_assemble)]
```
5. Glimma visuliation
```
strst trinity_run_glimma --trinity_dir [An absolute directory to trinity_dir] --edger_dir [edgeR DE analysis output dir] --sample_file [An absolute directory to sample_file (generated from trinity_assemble)]
```
