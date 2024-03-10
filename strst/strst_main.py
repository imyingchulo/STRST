import os
from pathlib import Path
from subprocess import call
import csv

# This pipeline is focused on "pair-end" RNA-seq data

# 1. Pre-alingment quality control (fastqc)
# 2. STAR (RNA-seq alingment)
# 3. post-alingment quality control (stat)
# 4. RNA-seq data visulization
# 5. Functional analysis

class STRST(object):
    def __init__(self, root_of_result = None):
        """
        # self means => class STRST
        At the beginning of the package,
        result and sub-folders will generate automatically when the user first run the package.
        """
        self.root_of_result = root_of_result if root_of_result else Path('{}/result'.format(os.getcwd())).resolve()
        self.result_dirs = ['pre_ran_seq_qc', 'pred_rna_seq', 'post_align_qc', 'alignment', 'visulization', 'functional_analysis_result', 'gene_counts', 'trimmed','assemble','trinity']
        self.init()

    def init(self):  # The setup_dir function is automatically create folders
        result = Path('{}'.format(self.root_of_result)).resolve()
        if not os.path.exists(result):
            os.mkdir(result)
        for dir in self.result_dirs:
            dest = "{}/{}".format(result, dir)
            if not os.path.exists(dest):
                os.mkdir(dest)

    def print_msg_box(msg, indent=1, width=None, title=None):
        """Print message-box with optional title."""
        lines = msg.split('\n')
        space = " " * indent
        if not width:
            width = max(map(len, lines))
        box = f'╔{"═" * (width + indent * 2)}╗\n'  # upper_border
        if title:
            box += f'║{space}{title:<{width}}{space}║\n'  # title
            box += f'║{space}{"-" * len(title):<{width}}{space}║\n'  # underscore
        box += ''.join([f'║{space}{line:<{width}}{space}║\n' for line in lines])
        box += f'╚{"═" * (width + indent * 2)}╝'  # lower_border
        print(box)

    def fastqc(self, fastqc, data):
        for i in os.listdir(data):
            if i.endswith('.fastq.gz'):
                os.system("{} {}/{} -o result/{}".format(fastqc, data, i, "pre_ran_seq_qc"))

    def trimmed_adapter(self, cutadapt, adapter_fwd, adapter_rev, data):
        # read files and create file name dictionary
        dict = {}
        for i in os.listdir(data):
            if 'fastq' and "R1" in i:
                file_name_r1 = i.strip().split("R1")[0]
                if file_name_r1 not in dict:
                    dict[i] = file_name_r1

        for i in os.listdir(data):
            if i.endswith('.fastq.gz') and "R1" in i:
                visited=set()
                identifier = dict[i]
                if identifier not in visited:
                    r1 = "{}/{}R1{}".format(data, identifier, ".fastq.gz")
                    r2 = "{}/{}R2{}".format(data, identifier, ".fastq.gz")
                    r1_output = "result/trimmed/{}R1.trimmed.fastq.gz".format(identifier)
                    r2_output = "result/trimmed/{}R2.trimmed.fastq.gz".format(identifier)
                    os.system("{} --minimum-length 1 -a {} -A {} -o {} -p {} {} {}".format(cutadapt, adapter_fwd, adapter_rev, r1_output, r2_output, r1, r2))
                visited.add(identifier)

    def star_genome_prep(self, star, gtf_dir, gtf_file, ref_fasta, extract_commands=''):
        os.system("{} --runMode genomeGenerate --genomeDir {} "
                  "--genomeFastaFiles {} --sjdbGTFfile {} --sjdbOverhang 49 {}".format(star, gtf_dir, ref_fasta, gtf_file, extract_commands))

    def star_alignment(self, star, data, gtf_dir, gtf_file, sep, threads="15", gz=" "):
        # read files name and creat the file dictionary
        dict = {}
        for i in os.listdir(data):
            if 'fastq' and "R1" in i:
                file_name_r1 = i.strip().split("R1")[0]
                if file_name_r1 not in dict:
                    dict[i] = file_name_r1

        if gz == "false":
            print("start to run STAR-alignment, RNA-seq data is not gzip")
            for i in os.listdir(data):
                if i.endswith('fastq') and "R1" in i:
                    identifier = dict[i]
                    visited = set()
                    if identifier not in visited:
                        r1 = "{}/{}R1{}{}".format(data, identifier, sep, "fastq")
                        r2 = "{}/{}R2{}{}".format(data, identifier, sep, "fastq")
                        output = identifier[:-1]
                        print("read1:{}\nread2:{}".format(r1, r2))
                        os.system("{} --runThreadN {} --genomeDir {} --sjdbGTFfile {} --genomeLoad NoSharedMemory --readMapNumber 10000"
                                  " --outSAMtype BAM Unsorted --quantMode GeneCounts "
                                  " --readFilesIn {} {} --outFileNamePrefix result/alignment/{}".format(
                            star, threads, gtf_dir, gtf_file, r2, r1, output))
                        print("Finished {} and {} STAR alignment! searching for next pair...".format(r1, r2))
                        visited.add(identifier)
            print("ALL works done!")

        elif gz == "true":
            print("start to run STAR-alignment, RNA-seq data is gzip")
            for i in os.listdir(data):
                if i.endswith('fastq.gz') and "R1" in i:
                    identifier = dict[i]
                    visited = set()
                    if identifier not in visited:
                        r1 = "{}/{}R1{}{}".format(data, identifier, sep, "fastq.gz")
                        r2 = "{}/{}R2{}{}".format(data, identifier, sep, "fastq.gz")
                        output = identifier[:-1]
                        print("read1:{}\nread2:{}".format(r1, r2))
                        print("{} --runThreadN {} --genomeDir {} --sjdbGTFfile {} --genomeLoad NoSharedMemory"
                                  "--outSAMtype BAM Unsorted --quantMode TranscriptomeSAM GeneCounts --outFilterScoreMinOverLread 0.3 --readMapNumber 10000"
                                  "--readFilesIn {} {} --outFileNamePrefix result/alignment/{} --readFilesCommand zcat".format(star, threads,
                                                                                                                               gtf_dir, gtf_file, r1, r2, output))
                        print("Finished {}R1{}{} and {}R2{}{} STAR alignment! searching for next pair...".format(identifier, sep, "fastq.gz", identifier, sep, "fastq.gz"))
                        visited.add(identifier)
            print("ALL works done!")

    def prepare_bamfiles_postqc(self, samtools, bam_dir="result/alignment/", ):
        for i in os.listdir(bam_dir):
            if i.endswith("bam"):
                output_name = i.strip().split(".bam")
                os.system("{} sort {}.bam {}_sorted.bam".format(samtools, output_name, output_name))
                os.system("{} index {}_sorted.bam".format(samtools, output_name))
                os.system("{} flagstat {}_sorted.bam {}_sorted.flagstat".format(samtools, output_name, output_name))

    def htseq_count(self, htseq, gtf_file, bam_file, output, commands, output_dir="result/gene_counts"):
        # -a 10 -r pos -t gene -i gene_id -m intersection-nonempty
        os.system("{} {} -f bam {} {} > {}/{}".format(htseq, commands, bam_file, gtf_file, output_dir, output))

    def trimmomatic_pairend(self,trimmomatic, data, sep, ends="fastq.gz", filter_condition="ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36"):
        # read files and create file name dictionary
        dict = {}
        for i in os.listdir(data):
            if i.endswith(ends) and "R1" in i:
                file_name_r1 = i.strip().split("R1")[0]
                if file_name_r1 not in dict:
                    dict[i] = file_name_r1

        for i in os.listdir(data):
            if i.endswith(ends) and "R1" in i:
                identifier = dict[i]
                visited = set()
                if identifier not in visited:
                    right = "{}/{}R1{}".format(data, identifier, sep)
                    left = "{}/{}R2{}".format(data, identifier, sep)
                    output = identifier[:-1]
                    print("start to assemble:\nright:{}\nleft:{}".format(right, left))
                    os.system("{} PE -trimlog {}.log {}{} {}{} {}pair.trimmed.{} {}unpair.trimmed.{} {}pair.trimmed.{} {}unpair.trimmed.{} {}".format(trimmomatic, output, right, ends, left, ends, right, ends, right, ends, left, ends, left, ends, filter_condition))

    # Trinity
    def trinity_assemble(self, trinity, data, sep, cpu, memory, sample_list='/home/ubuntu/sample_input.txt', output_prefix='result/trinity'):
        # read files and create file name dictionary
        dict = {}
        for i in os.listdir(data):
            if 'fastq' and "R1" in i:
                file_name_r1 = i.strip().split("R1")[0]
                if file_name_r1 not in dict:
                    dict[i] = file_name_r1

        # creat sample file
        tmp =[]
        sample_input=[]
        for i in os.listdir(data):
            if 'fastq' and "R1" in i:
                identifier = dict[i]
                visited = set()
                if identifier not in visited:
                    right = "{}/{}R1{}{}".format(data, identifier, sep, "fastq.gz")
                    left = "{}/{}R2{}{}".format(data, identifier, sep, "fastq.gz")
                    condition = identifier[:-1]
                    tmp.append(("{}\t{}\t{}\t{}".format(condition, condition, right, left)))
                    sample_input = '\n'.join(tmp)
                visited.add(identifier)
        with open('{}'.format(sample_list), 'w') as fout:
            fout.write(sample_input)

        print("Start to run Trinity")
        # print("{} --seqType fq --max_memory {} --samples_file {} --output {} --CPU {}".format(trinity, memory, sample_list, output_prefix, cpu))
        call("{} --seqType fq --max_memory {} --samples_file {} --min_kmer_cov 2 --no_parallel_norm_stats --output {} --CPU {}".format(trinity, memory, sample_list, output_prefix, cpu), shell=True)
        print("ALL works done!")

    # Estimate expression levels for each transcript in each sample
    def trinity_estimate_expression(self, trinity, method, seqtype, fasta, samplelist, outputdir):
        if method == "salmon":
            print('start to estimate expression with salmon')
            print("salmon indexing")
            os.system("salmon index --index quasi --type quasi --transcripts {}".format(fasta))
            print("finished index")
            print("Run trinity-estimate-expression")
            with open(samplelist, 'r') as fin:
                for i in fin.readlines():
                    input1 = i.split("\t")[2]
                    input2 = i.split("\t")[3]
                    output = "{}/{}".format(outputdir, i.split("\t")[1])
                    print("working on {} and {}".format(input1, input2))
                    os.system("{}/util/align_and_estimate_abundance.pl --seqType {} --transcripts {} --est_method {} --trinity_mode --output_dir {} --left {} --right {} ".format(trinity, seqtype, fasta, method, output, input1, input2))
                    print("finished {} and {}, finding for next pa``ir".format(input1, input2))
            print("All works done!")

        elif method == "RSEM":
            print('start to estimate expression with RSEM')
            with open(samplelist, 'r') as fin:
                for i in fin.readlines():
                    input1 = i.split("\t")[2]
                    input2 = i.split("\t")[3]
                    output = "{}/{}".format(outputdir, i.split("\t")[1])
                    print("working on {} and {}".format(input1, input2))
                    os.system("{}/util/align_and_estimate_abundance.pl --prep_reference --seqType {} --transcripts {} --est_method {} --trinity_mode --aln_method bowtie2 --output_dir {} --left {} --right {} ".format(trinity, seqtype, fasta, method, output, input1, input2))
                    print("finished {} and {}, finding for next pair".format(input1, input2))
            print("All works done!")


    def trinity_abundance_estimates_to_matrix(self, trinity_dir, data_dir, gene_map, output_name):
        dict= {}
        for i in os.listdir(data_dir):
            if i.endswith("_L1") or i.endswith("_L2"):
                for j in os.listdir("{}/{}".format(data_dir, i)):
                    if j.endswith(".isoforms.results"):
                        value = '{}/{}/{}'.format(data_dir, i, j)
                        # print(value)
                        key = "key"
                        if key not in dict:
                            dict[key] = [] #initialize as an empty list
                        dict[key].append(value)
        call("{}/util/abundance_estimates_to_matrix.pl --est_method RSEM {} --gene_trans_map {} --out_prefix {}".format(trinity_dir,str(dict[key]).replace("'", "").replace("]", "").replace("[", "").replace(",", ""),gene_map, output_name), shell=True)


    def trinity_run_de_analysis(self, trinity_dir, method, matrix, output_name, sample_file):
        if method == "edgeR":
            call("{}/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix {} --method edgeR --dispersion 0.1 --output {}".format(trinity_dir, matrix, output_name), shell=True)
        elif method == "DESeq2":
            call("{}/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix {} --method DESeq2 --samples_file {} --output {}".format(trinity_dir, matrix, sample_file,output_name), shell=True)

    def trinity_run_glimma(self, trinity_dir, edger_dir, sample_file):
        for i in os.listdir(edger_dir):
            if i.endswith(".edgeR.DE_results"):
                de = i.strip().split(".edgeR.DE_results")[0]
                call("{}/Analysis/DifferentialExpression/Glimma.Trinity.Rscript --samples_file {} --DE_results {}{}.edgeR.DE_results --counts_matrix {}{}.edgeR.count_matrix".format(trinity_dir, sample_file, edger_dir, de, edger_dir, de),shell=True)


    # def trinotate_gene_annotate(self, transDecoder_dir, fasta, trinotate_dir, working_dir):
    #     print("preparing the database...")
    #     call( " -t {}".format(transDecoder_dir, "TransDecoder.LongOrfs", fasta), shell=True)
    #
    #
    #     call("cp ../data/trinotate_data/Trinotate.boilerplate.sqlite  {}/Trinotate.sqlite".format(trinotate_dir, working_dir),shell=True)
    #     call("chmod 644 Trinotate.sqlite", shell=True
    #
    #     "TRINOTATE_HOME/Trinotate Trinotate.sqlite init \
    #      --gene_trans_map ../trinity_out_dir/Trinity.fasta.gene_trans_map \
    #      --transcript_fasta ../trinity_out_dir/Trinity.fasta \
    #      --transdecoder_pep Trinity.fasta.transdecoder.pep"





    # after assemble QC
    def trinity_access_read_counts(self, bowtie2, fasta, sample_list, outputdir):
        print("build a bowtie2 index for the transcriptome")
        os.system("bowtie2-build {} {}".format(fasta, fasta))
        print("bowtie2 index done")

        print("perform the alignment to just capture the read alignment statistics")
        with open(sample_list,'r') as fin:
            for i in fin.readlines():
                input1=i.split("\t")[2]
                input2=i.split("\t")[3]
                stat_output="{}/{}".format(outputdir, i.split("\t")[1])
                print("working on {} and {}".format(input1, input2))
                os.system("{} --local --no-unal -x {} -q -1 {} -2 {} | samtools sort -o {}.coordSorted.bam".format(bowtie2, fasta, input1, input2, stat_output))
                print("finished {} and {}, finding for next pair".format(input1, input2))
        print("All works done!")









    # def reference_prep(self, makeblastdb, fasta):

        #makeblastdb -in gtf/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa -dbtype nucl






    #brew install brewsci/bio/trimmomatic
    #brew install brewsci/bio/trimmomatic

    # def trinity_prep_n_alignment(self, gmap, genome_fasta, genome_fasta_name, data = 'result/assemble'):
    #     print("start to prepare ref genome:{}".format(genome_fasta))
    #     os.system("gmap_build -d {} -D . -k 13 {}".format(gmap, genome_fasta_name, genome_fasta)) #gmap_build -d Arabidopsis_thaliana.TAIR10.dna.toplevel -D gtf
    #     print("Genome prepared!")
    #     print("start to run alignment with GMAP")
    #     for i in data:
            # gmap -D . -d Arabidopsis_thaliana.TAIR10.dna.toplevel result/assemble/CL2_HT7M2DRXX_L2_both.fa -f samse > CL2_HT7M2DRXX_L2_both.sam


if __name__ == '__main__':
    # r = STRST()
    STRST.print_msg_box('test')













