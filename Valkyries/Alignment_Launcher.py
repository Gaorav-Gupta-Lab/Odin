"""

@author: Dennis A. Simpson
         RTP Genomics
         Chapel Hill, NC  27517
@copyright: 2020

"""

import subprocess
from Valkyries import Tool_Box

__author__ = 'Dennis A. Simpson'
__version__ = "0.5.0"


class AlignmentLauncher:
    """
    Builds the commands and launches the aligner.
    """

    def __init__(self, args, log, paired_end):

        self.args = args
        self.paired_end = paired_end
        self.log = log

    def run_bwa_aligner(self, fq1_name, fq2_name, sam_file, bam_file):
        threads = int(self.args.Spawn)
        aligner_options = getattr(self.args, "Aligner_Options", "")

        if self.args.BWA_Method == "mem":
            self.log.info("Begin alignment of {0} and {1} with BWA mem.".format(fq1_name, fq2_name))
            cmd = "bwa mem -t {0} {1} {2} {3} {4} > {5}" \
                .format(threads, aligner_options, self.args.Aligner_RefSeq, fq1_name, fq2_name, sam_file)
            self.log.debug(cmd)
            subprocess.run([cmd], shell=True)

            self.log.info("Alignment complete, begin SAM to BAM conversion.")
            cmd = "samtools view -bh {0} -o {1}".format(sam_file, bam_file)
            subprocess.run([cmd], shell=True)
            self.log.debug(cmd)

        elif self.args.BWA_Method == "aln":
            self.log.info("Begin alignment of {0} and {1} with BWA aln.".format(fq1_name, fq2_name))
            sai_file1 = fq1_name.replace(".fastq", ".sai")
            sai_file2 = fq2_name.replace(".fastq", ".sai")

            cmd1 = "bwa aln {0} {1} {2} {3} > {4}" \
                .format(threads, aligner_options, self.args.Aligner_RefSeq, fq1_name, sai_file1)
            cmd2 = "bwa aln {0} {1} {2} {3} > {4}" \
                .format(threads, aligner_options, self.args.Aligner_RefSeq, fq2_name, sai_file2)

            subprocess.run([cmd1], shell=True)
            subprocess.run([cmd2], shell=True)
            self.log.info("Alignment complete, begin SAI to SAM conversion.")

            if self.paired_end:
                cmd3 = "bwa sampe {0} {1} {2} {3} {4} > {5}" \
                    .format(self.args.Aligner_RefSeq, sai_file1, sai_file2, fq1_name, fq2_name, sam_file)

            elif not self.paired_end:
                cmd3 = "bwa {0} samse {1} {2} > {3}".format(self.args.Aligner_RefSeq, sai_file1, fq1_name,
                                                            sam_file)

            subprocess.run([cmd3], shell=True)

            self.log.info("SAI to SAM conversion complete. Begin SAM to BAM conversion.")
            cmd = "samtools view -bh {0} -o {1}".format(sam_file, bam_file)
            subprocess.run([cmd], shell=True)

        else:
            self.log.error("No valid BWA alignment method provided")
            raise SystemExit(1)

    def run_bowtie(self, fq1_name, fq2_name, sam_file, bam_file):
        if self.paired_end:
            sub_cmd = " -1 {0} -2 {1}".format(fq1_name, fq2_name)
            message = "file {0} and {1}".format(fq1_name, fq2_name)
        else:
            sub_cmd = " -U {0}".format(fq1_name)
            message = "file {0}".format(fq1_name)

        if self.args.local:
            bowtie2_local = " --local --ma {0}".format(self.args.ma)
        else:
            bowtie2_local = ''

            self.log.info("Begin alignment of {0} with Bowtie2".format(message))

        # Construct command for aligner and call subprocess to run.
        cmd = "bowtie2 {0} --trim5 {1} --trim3 {2} -x {3} {4} -S {5}" \
            .format(bowtie2_local, self.args.trim5, self.args.trim3, self.args.Aligner_RefSeq, sub_cmd,
                    sam_file)
        subprocess.run([cmd], shell=True)

        self.log.info("Alignment complete. Converting SAM format to BAM format")
        # Construct command for samtools and call subprocess to run.
        cmd = "samtools view -bh {0} -o {1}".format(sam_file, bam_file)
        subprocess.run([cmd], shell=True)

    def run_aligner(self, outfile_list_dict):
        """
        This will run our aligner of choice.  Final output is a BAM file.
        :param outfile_list_dict:
        :return:
        """
        for sample_index in outfile_list_dict:
            bam_file = "{}{}_{}.bam".format(self.args.WorkingFolder, self.args.JobName, sample_index)
            fq1_name = outfile_list_dict[sample_index][0]
            fq2_name = None
            if self.paired_end:
                fq2_name = outfile_list_dict[sample_index][1]

            sam_file = "{}{}_{}.sam".format(self.args.WorkingFolder, self.args.JobName, sample_index)
            if bam_file is None:
                bam_file = fq1_name.replace(".fastq", ".bam")

            if self.args.Aligner == "BWA":
                self.run_bwa_aligner(fq1_name, fq2_name, sam_file, bam_file)

            elif self.args.Aligner == "Bowtie2":
                self.run_bowtie(fq1_name, fq2_name, sam_file, bam_file)

            else:
                self.log.error("\033[1;31m***Warning:\033[m {0} not allowed.  Currently only BWA or Bowtie2 supported."
                               .format(self.args.Aligner))
                raise SystemExit(1)

            self.log.info("File conversion process complete.  Removing temporary files.")
            if not self.args.Demultiplex:
                Tool_Box.delete([fq1_name, fq2_name, sam_file, "{0}{1}_{2}.log"
                                .format(self.args.WorkingFolder, self.args.JobName, fq1_name)])
            else:
                Tool_Box.delete([sam_file, "{0}{1}_{2}.log"
                                .format(self.args.WorkingFolder, self.args.JobName, fq1_name)])

