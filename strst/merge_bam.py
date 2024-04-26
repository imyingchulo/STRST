import os


def merge_bam_n_stat(self, samtools, bam_dir, output_name):
    def bamfiles_postqc(bam_dir):
        for i in os.listdir(bam_dir):
            if i.endswith("_merged.bam"):
                bam_name = i.strip().split("_merged.bam")[0]
                print("{} sort -o {}/{}_merged_sorted.bam {}/{}_merged.bam".format(samtools, bam_dir, bam_name, bam_dir, bam_name))
                print("{} index {}/{}_merged_sorted.bam".format(samtools, bam_dir, bam_name))
                print("{} flagstat {}/{}_merged_sorted.bam".format(samtools, bam_dir, bam_name))

    tmp = []
    for i in os.listdir(bam_dir):
        if i.endswith(".bam"):
            tmp.append("{}/{}".format(bam_dir, i))
    bam_files = ' '.join(tmp)
    print("{} merge -o {}/{}_merged.bam {}".format(samtools, bam_dir, output_name, bam_files))
    bamfiles_postqc(bam_dir)

if __name__ == '__main__':
    merge_bam_n_stat()
