import click
from strst.strst_main import STRST


@click.group()
def main():
    pass

# @click.command()
# def test():
#     print("hello world")

@click.command()
@click.option( '--fastqc', metavar='<str>', help='An absolute directory to fastqc. to check:type "which fastqc" get the full path to fastqc, then paste the path after "--fastqc"' )
@click.option( '--data', metavar='<str>', required=True, help='An absolute directory to rna-seq data' )
@click.option('--outputfolder', metavar='<str>', help='fast file QC output folder; default:pre_ran_seq_qc')
def fastqc(fastqc, data, outputfolder):
    strst = STRST()
    strst.fastqc(fastqc=fastqc, data=data, outputfolder=outputfolder)

@click.command()
@click.option( '--cutadapt', metavar='<str>', help='An absolute directory to cutadapt. to check:type "which cutadapt get the full path to cutadapt, then paste the path after "--cutadapt"' )
@click.option( '--adapter_fwd', metavar='<str>', required=True, help='forward adapter sequence' )
@click.option( '--adapter_rev', metavar='<str>', help='reverse adapter sequence' )
@click.option( '--data', metavar='<str>', required=True, help='An absolute directory to rna-seq data' )
@click.option('--extention', metavar='<str>',  required=True, help='fastq.gz or fastq')
def trimmed_adapter(cutadapt, adapter_fwd, adapter_rev, data, extention):
    strst = STRST()
    strst.trimmed_adapter(cutadapt=cutadapt, adapter_fwd=adapter_fwd, adapter_rev=adapter_rev, data=data, extention=extention)

@click.command()
@click.option( '--star', metavar='<str>', required=True, help='An absolute directory to star. to check:type "which star get the full path to star, then paste the path after "--star"' )
@click.option( '--gtf_dir', metavar='<str>', required=True, help='absolute directory to gft file; example: "/home/ref/gft_dir"' )
@click.option( '--gtf_file', metavar='<str>', required=True, help='full path to gft file; example: "/home/ref/gft_dir/sample.gtf"' )
@click.option( '--ref_fasta', metavar='<str>', required=True, help='full path to reference fasta file; example: "/home/fasta_dir/sample.fasta"' )
def star_genome_prep(star, gtf_dir, gtf_file, ref_fasta):
    strst = STRST()
    strst.star_genome_prep(star=star, gtf_dir=gtf_dir, gtf_file=gtf_file, ref_fasta=ref_fasta)

@click.command()
@click.option( '--star', metavar='<str>', required=True, help='An absolute directory to STAR. to check:type "which STAR get the full path to star, then paste the path after "--star"' )
@click.option( '--data', metavar='<str>', required=True, help='RNA-seq data' )
@click.option( '--gtf_dir', metavar='<str>', help='absolute directory to gft file; example: "/home/ref/gft_dir"' )
@click.option( '--gtf_file', metavar='<str>', required=True, help='gft file name; example: "sample.gtf"' )
@click.option( '--sep', metavar='<str>', required=True, help='charaters after R1/R2, i.e., file named: 1-G1_HNFVVDRXX_L1_R1.fastq.gz, put "."; file named: 1-G1_HNFVVDRXX_L1_R1.trimmed.fastq.gz, put ".trimmed."' )
@click.option( '--gz', metavar='<str>', required=True, help='if RNA-seq data is .gz, put "true", if not put "false"' )
@click.option( '--threads', metavar='<int>', required=True, help='number of threads, default=15' )
def star_alignment(star, data, gtf_dir, gtf_file, sep, threads, gz=" "):
    strst = STRST()
    strst.star_alignment(star=star, data=data, gtf_dir=gtf_dir, gtf_file=gtf_file, sep=sep, threads=threads, gz=gz)

@click.command()
@click.option( '--samtools', metavar='<str>', required=True, help='An absolute directory to samtools. to check:type "which samtools get the full path to samtools, then paste the path after "--samtools"' )
@click.option( '--bam_dir', metavar='<str>', required=True, help='absolute path to bam files (the files after STAR), default="result/alignment/"' )
def bamfiles_postqc(samtools, bam_dir):
    strst = STRST()
    strst.bamfiles_postqc(samtools=samtools, bam_dir=bam_dir)

@click.command()
@click.option('--qualimap', metavar='<str>', help='An absolute directory to qualimap. to check:type "which qualimap" get the full path to qualimap, then paste the path after" --qualimap"' )
@click.option('--bam_dir', metavar='<str>', required=True, help='absolute path to bam files (the files after STAR), default="result/alignment/"' )
@click.option('--gff', metavar='<str>', required=True, help='An absolute directory to gff or gtf' )
@click.option('--output_dir', metavar='<str>', help='output folder; i.e. /home/output')
def qualimap_qc(qualimap, bam_dir, gff, output_dir):
    strst = STRST()
    strst.qualimap_qc(qualimap=qualimap, bam_dir=bam_dir, gff=gff, output_dir=output_dir)

@click.command()
@click.option( '--salmon', metavar='<str>', required=True, help='An absolute directory to salmon. to check:type "which salmon get the full path to salmon, then paste the path after "--salmon"' )
@click.option('--fasta', metavar='<str>', required=True, help='absolute path to fasta')
@click.option( '--output_dir', metavar='<str>', required=True, help='the dir for output file')
@click.option( '--bam_dir', metavar='<str>', required=True, help='absolute path to bam files' )
def salmon_quantification(salmon, fasta, bam_dir,  output_dir):
    strst = STRST()
    strst.salmon_quantification(salmon=salmon, fasta=fasta, bam_dir=bam_dir, output_dir=output_dir)

@click.command()
@click.option('--htseq', metavar='<str>', required=True, help='An absolute directory to htseq. to check:type "which htseq get the full path to htseq, then paste the path after "--htseq"')
@click.option('--gtf_file', metavar='<str>', required=True, help='absolute path to gtf_file')
@click.option('--bam_dir', metavar='<str>', required=True, help='absolute path to bam files')
@click.option('--output_dir', metavar='<str>', required=True, help='output dir, default="result/gene_counts"')
def htseq_count(htseq, gtf_file, bam_dir, output_dir="result/gene_counts"):
    strst = STRST()
    strst.htseq_count(htseq=htseq, gtf_file=gtf_file, bam_dir=bam_dir, output_dir=output_dir)

@click.command()
@click.option('--rsem', metavar='<str>', required=True, help='An absolute directory to rsem-calculate-expression')
@click.option('--fasta_file', metavar='<str>', required=True, help='absolute path to fasta')
@click.option('--bam_dir', metavar='<str>', required=True, help='absolute path to bam files')
@click.option('--output_dir', metavar='<str>', required=True, help='output dir, default="result/gene_counts"')
def rsem_count(rsem, fasta_file, bam_dir, output_dir="result/gene_counts"):
    strst = STRST()
    strst.rsem_count(rsem=rsem, fasta_file=fasta_file, bam_dir=bam_dir, output_dir=output_dir)



@click.command()
@click.option('--trinity', metavar='<str>', required=True, help='An absolute directory to trinity')
@click.option('--data', metavar='<str>', required=True, help='absolute path to RNA-seq data')
@click.option('--sample_list', metavar='<str>', required=True, help='sample_list; default="/home/ubuntu/sample_input.txt"')
@click.option('--sep', metavar='<str>', required=True, help='charaters after R1/R2, i.e., file named: 1-G1_HNFVVDRXX_L1_R1.fastq.gz, put "."; file named: 1-G1_HNFVVDRXX_L1_R1.trimmed.fastq.gz, put ".trimmed."' )
@click.option('--cpu', metavar='<int>', required=True, help='request cpu')
@click.option('--memory', metavar='<int>', required=True, help='request memory')
@click.option('--output_prefix', metavar='<str>', help='output dir, default="result/trinity"')
def trinity_assemble(trinity, data, sample_list, sep, cpu, memory, output_prefix):
    strst = STRST()
    strst.trinity_assemble(trinity=trinity, data=data, sample_list=sample_list, sep=sep, cpu=cpu, memory=memory, output_prefix=output_prefix)

@click.command()
@click.option('--trinity', metavar='<str>', required=True, help='An absolute directory to trinity (./Trinity)')
@click.option('--method', metavar='<str>', required=True, help='salmon or RSEM')
@click.option('--seqtype', metavar='<str>', required=True, help='fq or fa')
@click.option('--fasta', metavar='<str>', required=True, help='absolute path to meta-fasta')
@click.option('--samplelist', metavar='<str>', required=True, help='An absolute directory to samplelist (generated from trinity_assemble)')
@click.option('--outputdir', metavar='<str>', required=True, help='Full path to output')
def trinity_estimate_expression(trinity, method, seqtype, fasta, samplelist, outputdir):
    strst = STRST()
    strst.trinity_estimate_expression(trinity=trinity, method=method, seqtype=seqtype, fasta=fasta, samplelist=samplelist, outputdir=outputdir)

@click.command()
@click.option('--trinity_dir', metavar='<str>', required=True, help='An absolute directory to trinity_dir')
@click.option('--data_dir', metavar='<str>', required=True, help='Path to access folders')
@click.option('--gene_map', metavar='<str>', required=True, help='Path to gene_map, includes file name and extension')
@click.option('--output_name', metavar='<str>', required=True, help='output file name')
def trinity_abundance_estimates_to_matrix(trinity_dir, data_dir, gene_map, output_name):
    strst = STRST()
    strst.trinity_abundance_estimates_to_matrix(trinity_dir=trinity_dir, data_dir=data_dir, gene_map=gene_map, output_name=output_name)

@click.command()
@click.option('--trinity_dir', metavar='<str>', required=True, help='An absolute directory to trinity_dir')
@click.option('--method', metavar='<str>', required=True, help='edgeR or DESeq2')
@click.option('--matrix', metavar='<str>', required=True, help='Path to matrix, includes file name and extension')
@click.option('--output_name', metavar='<str>', required=True, help='output file name')
@click.option('--sample_file', metavar='<str>', required=True, help='An absolute directory to sample_file (generated from trinity_assemble)')
def trinity_run_de_analysis(trinity_dir, method, matrix, output_name, sample_file):
    strst = STRST()
    strst.trinity_run_de_analysis(trinity_dir=trinity_dir, method=method, matrix=matrix, output_name=output_name, sample_file=sample_file)

@click.command()
@click.option('--trinity_dir', metavar='<str>', required=True, help='An absolute directory to trinity_dir')
@click.option('--edger_dir', metavar='<str>', required=True, help='edgeR DE analysis output dir')
@click.option('--sample_file', metavar='<str>', required=True, help='An absolute directory to sample_file (generated from trinity_assemble)')
def trinity_run_glimma(trinity_dir, edger_dir, sample_file):
    strst = STRST()
    strst.trinity_run_glimma(trinity_dir=trinity_dir, edger_dir=edger_dir, sample_file=sample_file)


main.add_command(fastqc)
main.add_command(trimmed_adapter)
main.add_command(star_genome_prep)
main.add_command(star_alignment)
main.add_command(bamfiles_postqc)
main.add_command(qualimap_qc)
main.add_command(salmon_quantification)
main.add_command(htseq_count)
main.add_command(rsem_count)
main.add_command(trinity_assemble)
main.add_command(trinity_estimate_expression)
main.add_command(trinity_abundance_estimates_to_matrix)
main.add_command(trinity_run_de_analysis)
main.add_command(trinity_run_glimma)
