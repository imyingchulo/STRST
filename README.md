# STRST: STAR and Trinity RNA-seq Tool

Pipelines for RNA-seq analysis: STAR and Trinity.
This package using two approaches for RNA-seq a) STAR and b) Trinity. 
The package includes the data cleaning, alignment/assembly, DE/Annotation analysis and data visualization for GO.

please find more detail for STAR and Trinity-wiki in the links below:
STAR: https://github.com/alexdobin/STAR
Trinity: https://github.com/trinityrnaseq/trinityrnaseq/wiki

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
After setting up virtualenv, install 

3. Install package

```shell
$ pip install -e .
```

## Usage guidance

:octocat: Functions include in the package:

:memo: Function explain - `function name`

1. Pre-alignment quality control - `fastqc`
2. Trimmed adapter - `trimmed_adapter`
3. Prepare genome for STAR alignment - `star_genome_prep`
4. STAR alignment - `star_alignment` :construction:
5. Post STAR alignment - `quality control` :construction: 
6. Prepare files for post STAR alignment - `prepare_bamfiles_postqc` :construction: 
7. Count the Transcripts - `htseq_count` :construction:  
8. Trinity assembly - `trinity_assemble`
9. Estimate the expression after Trinity assembly - `trinity_estimate_expression`
10. Calculate the transcripts and transfer into a matrix - `trinity_abundance_estimates_to_matrix`
11. Trinity DE analysis - `trinity_run_de_analysis `
12. RNA-seq data visualization - `trinity_run_glimma` :construction:

###  Pre-alignment quality control (fastqc)
:construction: `under construction`

### RNA-seq alignment/assembly
:construction: `under construction`

###  post-alignment quality control
:construction: `under construction`

###  RNA-seq data visualization
:construction: `under construction`

###  Functional analysis
:construction: `under construction`

