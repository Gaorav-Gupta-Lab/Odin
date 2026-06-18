"""
June 15, 2019
    Dennis A. Simpson
    Version 0.20.1
    Refactor code according to pylint recommendations.  Specifically remove all eval functions.
September 3, 2018
    Dennis A. Simpson
    Version 0.18.6
    Minor bug fixes.  Changes in output formatting.  Removed the <> for insertions in the filtering section.  Done to
    make the output compatible with snpEff and other annotators.
July 24, 2018
    Dennis A. Simpson
    Version 0.18.1
    Changed read depth filter in variant calling such that both the forward and reverse reads must exceed the minimum
    read depth provided by user in the parameter file.
July 21, 2018
    Dennis A. Simpson
    Version 0.17.0
    Added filtering based on error matrix.  The pickle file containing the error matrix will automatically be placed
    in the Mimir directory no matter where the directory resides.  Still needs some work on the output so that
    annotation programs will read it.
June 24, 2018
    Dennis A. Simpson
    Version 0.17.0
    Insertions were not scoring correctly.  They are now.  Reports the consensus insertion.  Altered output format of
    raw data.  Insertions marked as <insertion> and deleted sequence is not in ALT column.  Changed name of output to
    <Job_Name>_raw_variants.txt
June 9, 2018
    Dennis A. Simpson
    Version 0.16.1
    Introduced a bug in Column_Cruncher that was not scoring SNV's.
June 8, 2018
    Dennis A. Simpson
    Version 0.16.0
    Major changes.  Changed how read overlaps are identified, changed main variant dictionary to Pandas DataFrame,
    changed column cruncher function to Cython.  Added a filter that removes the first aligned position from every
    read regardless of position in read.  Now runs at about 250,000 positions per second.  This is a 70x
    improvement in performance.
May 9, 2018
    Dennis A. Simpson
    Version 0.12.0
    Major bug in the region parsing function that was causing Skadi to skip some exons on some genes.  Also added a
    section to filter VCF files using control data, strand count, and minimum allele count.
May 2, 2018
    Dennis A. Simpson
    Version 0.11.0
    Made changes that improved efficiency by finally getting the pool to run the largest datasets first.  Scraped all
    references to the control BAM file from code.  For my purposes the control file will be the same for a single
    sequencing run so it makes sense to process the control file separate and then use the VCF file from that to filter
    variants.
April 20, 2018
    Dennis A. Simpson
    Version 0.10.2
    Made a few bug fixes that was causing some regions to loop.  Also trying out globals for several of the file objects
April 11, 2018
    Dennis A. Simpson
    Version 0.9.1
    Deadlocking is solved.  Code is cleaner but v0.9.0 was even slower on Longleaf.  I am attempting to have the largest
    computational blocks processed first so the multiprocessor pool stays full longer.  I was using the size of the
    region as a proxy for how much information it contains.  This does not work.  This version actually counts the
    number of reads during the preanalysis steps and build a dictionary that has the key being the number of reads.  If
    two regions share the same number of reads then the region count is incremented by one.

April 3, 2018
    Dennis A. Simpson
    Still deadlocking on Longleaf.  This refactor is yet another attempt to solve this.  It is also about 40% faster on
    my development machine.
March 29, 2018
    Dennis A. Simpson
    This is a variant caller that can use UMI's and read strand counts.  Allows simple definition of target regions as
    genes, exons, or user defined.  While it does work it is also slow.  On the UNC Longleaf cluster it has a bad habit
    deadlocking.  So far I have not been able to reproduce this bug on other platforms.


@author: Dennis A. Simpson
         RTP Genomics, LLC
         Durham, NC
@copyright: 2019
"""

import csv
import itertools
import subprocess
from contextlib import suppress
import datetime
import os
import collections
from collections import Counter
import pysam
from natsort import natsort

import vcf
import pathos
import pandas
import scipy
from scipy.stats import gmean
import dill
from Odin import Column_Cruncher
from Odin import Utilities
import Valkyries.Ensembl as Ensembl
import Valkyries.Tool_Box as Tool_Box
from Valkyries import BamTools

__author__ = 'Dennis A. Simpson'
__version__ = "0.20.1"
__package__ = 'Odin'


class ReportGenerator:
    def __init__(self, args):
        self.args = args
        self.output_dict = collections.defaultdict(lambda: collections.defaultdict(str))
        self.chromosome_key = {
            "chr1": "NC_000001.11", "chr2": "NC_000002.12", "chr3": "NC_000003.12", "chr4": "NC_000004.12",
            "chr5": "NC_000005.10", "chr6": "NC_000006.12", "chr7": "NC_000007.14", "chr8": "NC_000008.11",
            "chr9": "NC_000009.12", "chr10": "NC_000010.11", "chr11": "NC_000011.10", "chr12": "NC_000012.12",
            "chr13": "NC_000013.11", "chr14": "NC_000014.09", "chr15": "NC_000015.10", "chr16": "NC_000016.10",
            "chr17": "NC_000017.11", "chr18": "NC_000018.10", "chr19": "NC_000019.10", "chr20": "NC_000020.11",
            "chr21": "NC_000021.09", "chr22": "NC_000022.11", "chrX": "NC_000023.11", "chrY": "NC_000024.10"
        }

    def generate_report(self, annotated_vcf):
        """

        :param annotated_vcf:
        :return:
        """
        column_header = True
        line_count = 0
        location = ""
        with open(annotated_vcf) as f:
            for l in f:
                variant_dict = collections.defaultdict(list)
                # Break each line into a list.
                l_list = [x for x in l.strip("\n").split("\t")]

                if not column_header:
                    # This is the VCF data.  Key is chr:position.

                    for d in l_list[7].split(";"):
                        k = d.split("=")[0]
                        new_d = d.split("=")[1].split(",")
                        variant_dict[k] = new_d
                    lx = self.parse_data(variant_dict, l_list)

                    if len(lx) > 0:
                        location += "{}\n".format(lx)

                elif column_header:
                    if l_list[0] == "#CHROM":
                        # These are the column names and define the end of the header section.
                        column_header = False

                    elif line_count == 3:
                        # Get the sample name from the VCF file.
                        sample_name = l_list[0].split("=")[1]

                    elif line_count != 1 or line_count != 4:
                        pass
                        # Line 1 is the date, line 4 is the source.  This builds a header without those or the sample

            out_string = ""
        for position in natsort.natsorted(self.output_dict['HIGH']):
            out_string += "{}\n".format(self.output_dict['HIGH'][position])
        out_string += "\n"

        for position in natsort.natsorted(self.output_dict['MODERATE']):
            out_string += "{}\n".format(self.output_dict['MODERATE'][position])
        out_string += "\n"

        for position in natsort.natsorted(self.output_dict['LOW']):
            out_string += "{}\n".format(self.output_dict['LOW'][position])

        outfile = open("{}temp_tables_{}.txt".format(self.args.WorkingFolder, self.args.Job_Name), 'w')
        outfile.write(out_string)

        return

    def parse_data(self, variant_dict, l_list):
        """

        :param variant_dict:
        :param l_list:
        :return:
        """
        min = 0.003
        position_string = ""
        position = "{}:{}".format(l_list[0], int(l_list[1]))

        # Insertions
        ins_position = "g.{}_{}ins" .format(l_list[1], int(l_list[1]) + 1)

        # Deletions
        del_position = ""
        snv_list = ["", "", "", ""]
        block_list = variant_dict["ANN"]

        for block in block_list:
            data = block.split("|")
            # I don't want to deal with this feature.
            if data[1] == 'sequence_feature':
                continue

            d10 = ""
            l1 = []
            variant = []
            effect = []
            freq_list = []
            insertion = False
            deletion = False
            if len(data[10]) > 4:
                d10 = data[10]
            if float(l_list[32]) > min and 'ins' in data[9]:
                for alt in l_list[4].split(","):
                    if len(alt) > len(l_list[3]):
                        insertion = True
                        ins_position += "{}".format(alt[1:])
                variant.append(data[9])
                effect.append(d10)
                freq = l_list[32]

            if float(l_list[31]) > min and 'del' in data[9]:
                deletion = True
                variant.append(data[9])
                effect.append(d10)
                freq = l_list[31]

                if len(l_list[3]) > 2:
                    del_position = "g.{}{}del".format(l_list[1], "_{}".format(int(l_list[1]) + len(l_list[3]) - 1))
                else:
                    del_position = "g.{}del".format(l_list[1])

            if data[0][0] == "G" and float(l_list[27]) > min:
                snv_list[0] = data[0]
                variant.append(data[9])
                effect.append(d10)
                freq = l_list[27]
            elif data[0][0] == "A" and float(l_list[28]) > min:
                snv_list[1] = data[0]
                variant.append(data[9])
                effect.append(d10)
                freq = l_list[28]
            elif data[0][0] == "T" and float(l_list[29]) > min:
                snv_list[2] = data[0]
                variant.append(data[9])
                effect.append(d10)
                freq = l_list[29]
            elif data[0][0] == "C" and float(l_list[30]) > min:
                snv_list[3] = data[0]
                variant.append(data[9])
                effect.append(d10)
                freq = l_list[30]

            snv = False
            s1 = ""
            snv_position = ""
            for s in snv_list:
                if s != "":
                    snv = True
                    l1.append("{}>{}".format(l_list[3], s))
            if len(l1) > 1:
                s1 = ",".join(l1)

            elif len(l1) == 1:
                s1 = l1[0]

            if snv:
                snv_position = "g.{}{}".format(l_list[1], s1)
            if snv and deletion and insertion:
                position_string = "{}:[{},{}, {}]".format(self.chromosome_key[l_list[0]], snv_position, del_position, ins_position)
            elif snv and insertion:
                position_string = "{}:[{},{}]".format(self.chromosome_key[l_list[0]], snv_position, ins_position)
            elif snv and deletion:
                position_string = "{}:[{},{}]".format(self.chromosome_key[l_list[0]], snv_position, del_position)
            elif deletion and insertion:
                position_string = "{}:[{},{}]".format(self.chromosome_key[l_list[0]], del_position, ins_position)
            elif snv:
                position_string = "{}:{}".format(self.chromosome_key[l_list[0]], snv_position)
            elif insertion:
                position_string = "{}:{}".format(self.chromosome_key[l_list[0]], ins_position)
            elif deletion:
                position_string = "{}:{}".format(self.chromosome_key[l_list[0]], del_position)

            if snv or deletion or insertion:
                gene = data[3]
                transcript = data[6]
                effect = list(set(effect))
                variant_string = ",".join(variant)
                effect_string = ",".join(effect)

                for i in range(27, 33):
                    if float(l_list[i]) > min:
                        freq_list.append(l_list[i])

                if data[2] == 'MODIFIER' or data[2] == 'LOW':
                    self.output_dict['LOW'][position] = "{}\t{}\t{}\t{}\t{}\t{}".format(position_string, gene, transcript, variant_string, effect_string, freq)
                else:
                    self.output_dict[data[2]][position] = "{}\t{}\t{}\t{}\t{}\t{}".format(position_string, gene, transcript, variant_string, effect_string, freq)

        return position_string


class RegionGenerator:
    """
    Returns a dictionary of contiguous regions and a sorted list of the number of reads in each region.
    """
    def __init__(self, args, log, ensembl_data):
        self._args = args
        self._log = log
        self._ensembl_data = ensembl_data
        self.sub_region_dict = collections.defaultdict(list)
        self.region_dictionary = collections.defaultdict(list)
        self. region_read_counts = []
        self.qbam = pysam.AlignmentFile(self._args.BAM_File)
        self.previous_ctg = ""
        self.previous_start = 0
        self.previous_stop = 1000000000

    def __sub_regions(self, locus, gene, exon=None):
        """
        Gathers all the sub regions in the target space.
        :param locus:
        :param gene:
        :param exon:
        :return:
        """
        current_locus_end = locus.end + int(self._args.Boundary_Padding)
        current_locus_start = locus.start - int(self._args.Boundary_Padding)
        primary_key = "chr{}:{}-{}".format(locus.contig, current_locus_start, current_locus_end)
        self.sub_region_dict[primary_key] = [gene, exon]

    def __exon_generator(self, gene):
        try:
            exon_list = self._ensembl_data.exon_ids_of_gene_name(gene)

        except ValueError:
            Tool_Box.debug_messenger("Gene Name Incorrect!!!!")
            self._log.warning("{} not found in Ensembl database.".format(gene))
            return

        track_list = []
        for exon in exon_list:
            exon_locus = self._ensembl_data.locus_of_exon_id(exon)

            # My method to eliminate duplicate exons.
            key = (exon_locus.start, exon_locus.end)
            if key not in track_list:
                track_list.append(key)
                yield exon

    def region_parsing(self):
        """
        Coordinates main work of class.  Returns dictionary and list when complete.
        :return:
        """
        self._log.info("Parsing regions")
        # Because the parser was written for a multicolumn bed file it returns a list of lists.
        # This makes that a list.
        target_list = [x[0] for x in Tool_Box.FileParser.indices(self._log, self._args.Target_File)]

        for gene in target_list:
            if self._args.Target == "Genes":
                locus = self._ensembl_data.locus_of_gene_id("".join(self._ensembl_data.gene_ids_of_gene_name(gene)))
                self.__sub_regions(locus, gene)
            elif self._args.Target == "Exons":
                for exon in self.__exon_generator(gene):
                    locus = self._ensembl_data.locus_of_exon_id(exon)
                    self.__sub_regions(locus, gene, exon)

        # Combine overlapping regions.
        new_start = True
        regions_processed = 0
        gene = ""
        for key in natsort.natsorted(self.sub_region_dict):
            region_count = len(self.sub_region_dict)
            ctg = key.split(":")[0]
            start = int(key.split(":")[1].split("-")[0])
            stop = int(key.split("-")[1])
            regions_processed += 1

            if new_start:
                new_start = False
                self.previous_start = start
                self.previous_stop = stop
                self.previous_ctg = ctg
                gene = self.sub_region_dict[key][0]

            elif start <= self.previous_stop and self.previous_ctg == ctg:
                if start < self.previous_start:
                    self.previous_start = start
                if self.previous_stop < stop:
                    self.previous_stop = stop

            elif not new_start or not self.previous_ctg == ctg:
                self.store_region(ctg, gene)
                self.previous_start = start
                self.previous_stop = stop
                self.previous_ctg = ctg
                gene = self.sub_region_dict[key][0]

            if region_count == regions_processed:
                regions_processed = 0
                self.store_region(ctg, gene)

        self.qbam.close()
        self.region_read_counts.sort(reverse=True)

        return self.region_dictionary, self.region_read_counts

    def store_region(self, ctg, gene):
        """
        Filters contiguous regions for number of reads and fills region dictionary and region count list.
        :param ctg:
        :param gene:
        :return:
        """
        new_key = "{}:{}-{}".format(self.previous_ctg, self.previous_start, self.previous_stop)
        read_count = self.qbam.count(contig=self.previous_ctg, start=self.previous_start, end=self.previous_stop)

        if read_count > int(self._args.Minimum_Read_Depth):
            # Prevent read counts equal read counts for two regions from clashing.
            while read_count in self.region_dictionary:
                read_count += 1
            # this list is gene, region, and read counts.  The counts are used to id the temp files.
            self.region_dictionary[read_count] = [gene, new_key, read_count]
            self.region_read_counts.append(read_count)
        self.previous_ctg = ctg
        new_start = True
        return new_start


def variant_search(args, log):
    """
    Entry point for variant caller.  Get here from Odin.py.
    :param args:
    :param log:
    :return:
    """

    if not os.path.isfile(args.BAM_File):
        log.error("BAM File {} Not Found.  Check File Name and Path.".format(args.BAM_File))
        raise SystemExit(1)

    ensembl_data = Ensembl.initialize(args, log)
    target_regions = RegionGenerator(args, log, ensembl_data)
    region_dictionary, region_read_counts = target_regions.region_parsing()

    indices = None
    # Currently the UNC folks are giving me demultiplexed FASTQ files.
    if getattr(args, "Index_File", None) is not None:
        indices = Tool_Box.FileParser.indices(log, args.Index_File)

    log.info("Spawning {} jobs to process {} regions".format(args.Spawn, len(region_read_counts)))
    p = pathos.multiprocessing.Pool(int(args.Spawn))

    # chunksize=1 is to ensure that the largest regions are loaded before the smaller ones.  Keeps the pool full.
    argvs = args, ensembl_data, log, indices, region_dictionary
    multiprocessor_data = \
        p.starmap(gene_processing, zip(itertools.repeat(argvs), region_read_counts), chunksize=1)

    log.info("Parallel Jobs Complete.")

    tumor_sample = args.BAM_File.split("/")[-1]
    ref_file = args.RefSeq.split("/")[-1]

    variant_data_dict = collections.defaultdict(lambda: collections.defaultdict(list))
    files = []

    log.debug("Read temporary VCF files into dictionary")
    # Read all the temporary vcf files into a dictionary then delete them.
    for k, v in region_dictionary.items():
        first_line = True
        with suppress(FileNotFoundError):
            with open("{}{}_{}_variants.txt".format(args.WorkingFolder, args.Job_Name, v[2]), 'r') as vcf_file:
                gene_vcf = csv.reader(vcf_file, delimiter='\t')
                for line in gene_vcf:
                    if first_line:
                        first_line = False
                        continue

                    chrom = line[0]
                    position = line[1]
                    variant_data_dict[chrom][position] = line
            files.append("{}{}_{}_variants.txt".format(args.WorkingFolder, args.Job_Name, v[2]))
    Tool_Box.delete(files)

    outstring = ""
    log.debug("Sort variant data dict and create outstring")
    # Sort the data and format it for writing vcf file.
    for chrom in natsort.natsorted(variant_data_dict):
        for position in natsort.natsorted(variant_data_dict[chrom]):
            outdata = "\t".join(variant_data_dict[chrom][position])
            outstring += "\n{}".format(outdata)

    log.debug("Format metatdata for VCF file")
    # Format metadata for vcf file.
    file_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    dp = '<ID=DP,Number=1,Type=Integer,Description="Read Depth">'
    svtype = '<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">'
    svlen = '<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">'
    allel_freq = '<ID=AF,Number=A,Type=Float,Description="Query Allele Frequency">'
    db = '<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 151, Annotation Release 108">'

    metadata = "##fileformat=VCFv4.3\n##filedate={0}\n##reference=file:{1}\n##sample={2}\n" \
               "##source=Skadi_v{3}\n##Ensembl_Release={10}\n##species={5}\n##INFO={5}\n##INFO={6}\n##INFO={7}\n" \
               "##FORMAT={8}\n##FORMAT={9}" \
        .format(file_date, ref_file, tumor_sample, __version__, args.Species, svtype, svlen, db,
                allel_freq, dp, Ensembl.__release__)

    query_orientation_meta = \
        "TotalForward\tTotalReverse\tfG\tfA\tfT\tfC\tfN\trG\trA\trT\trC\trN\tfDel\trDel\tfIns\trIns"

    metadata += "\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tGENE\tDP_Query\t{}"\
        .format(query_orientation_meta)

    metadata += outstring
    vcf_outfile = "{}{}_raw_variants.txt".format(args.WorkingFolder, args.Job_Name)
    outfile = open(vcf_outfile, "w")
    outfile.write(metadata)
    outfile.close()

    """
        Each multiprocessor data block is a list of two dictionaries.
        Query Family Size Data = block[0]
        Control Family Size Data = block[1]
    """
    qfamilies = Utilities.FamilyDistribution(log, args, args.BAM_File.split("/")[-1].split(".")[0])
    target_counts = collections.defaultdict(lambda: collections.defaultdict(int))
    for block in multiprocessor_data:
        qfamilies.data_collection("query", block["query"])
        target_counts["query"]["total"] += block["target_counts"]["total"]
        target_counts["query"]["reverse"] += block["target_counts"]["reverse"]

    qfamilies.data_processing()

    query_targeting = target_counts["query"]["total"] / BamTools.total_align_count(args.BAM_File)
    log.info("{} On-Target={}".format(tumor_sample, query_targeting))

    return vcf_outfile


class TargetProcessing:
    def __init__(self, argvs):
        self.args, self.ensembl_data, self.log, self.indices, self.region_dict = argvs
        self.qbam = pysam.AlignmentFile(self.args.BAM_File)
        self.refseq = pysam.FastaFile(self.args.RefSeq)
        # Dictionary for Family data
        self.qfamilies = collections.defaultdict(list)

        # Dataframe for variant data
        self.position_data_df = pandas.DataFrame(index=["snv", "del", "ins", "depth"])

        # Dictionary for targeting data
        self.qtargeted_count = collections.defaultdict(int)

        self.previous_locus_end = 0
        self.locus_count = 0

    def stacker(self, key):
        """

        :param key:
        :return:
        """
        target = self.region_dict[key][1]
        region_start = int(target.split(":")[1].split("-")[0])
        region_stop = int(target.split("-")[1])
        ctg = target.split(":")[0]
        gene = self.region_dict[key][0]
        region_length = region_stop-region_start

        # This is here in case something is very wrong with our coordinates.  This should never occur but does.
        if region_length < 1:
            self.log.warning("{} for {} has no length".format(target, gene))
            return

        # Get family data and gene targeting data for tumor and control.
        self.log.debug("Tumor Family Size and Target Counts for {}".format(target))

        self.qfamilies, self.qtargeted_count = \
            Utilities.family_counts(self.qbam.fetch(region=target), ctg, self.qfamilies, self.previous_locus_end,
                                    self.qtargeted_count, index="query")

        sample = "{}|{}".format(gene, target)
        self.pileup_processor(target, sample, region_length)

        return

    def read_generator(self, pileup_column):
        duplicate_tracking = []
        for read, quality, name in zip(pileup_column.pileups, pileup_column.get_query_qualities(),
                                       pileup_column.get_query_names()):
            umt = name.split(":")[0]
            if umt in duplicate_tracking or quality < int(self.args.Minimum_Base_Quality) or read.is_head:
                continue

            duplicate_tracking.append(umt)
            yield read, name

    def pileup_processor(self, target, sample, locus_length):
        """
        Process a Samtools pileup one position at a time.
        :param target:
        :param sample:
        :param locus_length:
        :return:
        """

        ctg = target.split(":")[0]
        progress_check = int(locus_length / 10)
        pileup = self.qbam.pileup(region=target, truncate=self.args.Strict_Boundaries, max_depth=100000000)

        # If a region is <6 nt the int is 0.  Tha is a big problem.  Actually should never occur but does.
        if progress_check < 1:
            progress_check = 1

        del_tracking = collections.defaultdict(list)
        column_count = 0
        for pileup_column in pileup:
            if column_count % progress_check == 0:
                self.log.debug("{}% of {} positions for {} complete"
                               .format(round((column_count / locus_length) * 100, 3), locus_length, sample))
            column_count += 1

            # Skip columns that don't contain enough reads.
            if len(pileup_column.get_query_names()) < int(self.args.Minimum_Read_Depth):
                continue

            # self.log.debug("Processing {} reads at {}".format(len(pileup_column.get_query_names()), position_key))
            position = "{}:{}".format(ctg, pileup_column.reference_pos)
            # start_time = time.time()

            column_data, del_tracking = \
                Column_Cruncher.column_cruncher(int(self.args.Minimum_Base_Quality), pileup_column, ctg[3:],
                                                pileup_column.reference_pos, self.refseq, del_tracking)
            self.position_data_df.loc[:, position] = pandas.Series(column_data, index=self.position_data_df.index)
            # self.log.debug("{} reads per second"
            #                .format(int(len(pileup_column.get_query_names()) / (time.time() - start_time))))

    def close_files(self):
        self.refseq.close()
        self.qbam.close()


def gene_processing(argvs, key):
    """
    :param argvs:
    :param key:
    :return:
    """

    args, ensembl_data, log, indices, region_dictionary = argvs
    region = region_dictionary[key][1]
    gene = region_dictionary[key][0]
    log.info("Begin: {}|{} containing {} reads".format(gene, region, key))

    # Branch off based on targeting choice.
    tp = TargetProcessing(argvs)
    tp.stacker(key)

    qfamilies = tp.qfamilies
    qfamilies["target_counts"] = tp.qtargeted_count

    variants = VariantCaller(log, args, region_dictionary[key], tp.refseq)
    variant_block_writer(args, region_dictionary[key], variants.process_locus(tp.position_data_df))

    log.info("Completed: {}|{}".format(gene, region))

    return qfamilies


def snv_processor(snv_list, strand_list, min_allele_count):
    snv_data_dict = {"G": [0, []],
                     "A": [0, []],
                     "T": [0, []],
                     "C": [0, []],
                     "N": [0, []],
                     "strand": []}

    snv_counts = Counter(snv_list)
    for snv, strand in zip(snv_list, strand_list):
        if snv_counts[snv] > min_allele_count:
            snv_data_dict[snv][0] += 1
            snv_data_dict[snv][1].append(strand)
            if snv != "N":
                snv_data_dict["strand"].append(strand)

    return snv_data_dict


def strand_counter(data_list):
    return data_list.count("forward"), data_list.count("reverse")


class VariantCaller:
    """
    This class gathers up and processes the variant data.
    """
    def __init__(self, log, args, region, refseq):
        self._log = log
        self._args = args
        self._region = region
        self._refseq = refseq
        self._variant_gene_data = {}
        self.vcf_snp = vcf.Reader(filename=args.dbSNP)

        # Data Frame for VCF output file.
        self.vcf_data_df = pandas.DataFrame(index=["snv", "del", "ins", "depth"])

    @property
    def variant_gene_data(self):
        return self._variant_gene_data

    def insertions(self, position_data_df, column_name, ref_base, position):
        """
        :param position:
        :param ref_base:
        :param position_data_df:
        :param column_name:
        :return:
        """

        gapped_alignment_dict = collections.defaultdict(list)
        ins_list = []
        fasta_data = ""
        count = 0
        consensus_ins_seq = ""
        query_orientation_ins_count = Counter(position_data_df.loc["ins", column_name][0])

        if len(position_data_df.loc["ins", column_name][0]) >= int(self._args.Minimum_Allele_Count):
            # Format insertion data as FASTA
            for strand, insertion in zip(position_data_df.loc["ins", column_name][0],
                                         position_data_df.loc["ins", column_name][1]):

                fasta_data += ">{}\n{}\n".format(count, insertion)
                count += 1

            # Create gapped alignment file in FASTA format using MUSCLE
            cmd = ['muscle', "-quiet", "-maxiters", "1", "-diags"]
            muscle = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
            muscle.stdin.write(fasta_data)
            muscle.stdin.close()
            gapped_alignment_dict.clear()

            # Put gapped sequences into dictionary for consensus calling
            for line in muscle.stdout:
                if line[0] != ">":
                    seq = line.strip("\n")
                    for i in range(len(seq)):
                        gapped_alignment_dict[i].append(seq[i])

            # Determine Consensus Sequence and Format it for Output
            for k in gapped_alignment_dict:
                algn_count = len(gapped_alignment_dict[k])
                c = collections.Counter(gapped_alignment_dict[k]).most_common(1)[0]
                if c[1] / algn_count > 0.51 and c[0] != "-":
                    consensus_ins_seq += c[0]
                elif c[0] != "-":
                    consensus_ins_seq += "N"
            ins_list.append("{}{}".format(ref_base, consensus_ins_seq))

        return ins_list, (query_orientation_ins_count["forward"], query_orientation_ins_count["reverse"])

    def process_locus(self, position_data_df):
        gene = self._region[0]
        region = self._region[1]
        self._log.debug("Processing Variant Data for {}|{}".format(gene, region))
        min_allele_count = int(self._args.Minimum_Allele_Count)
        region_length = position_data_df.shape[1]
        column_names = position_data_df.columns.tolist()
        position_count = 0

        vcf_columns = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "GENE", "DP_Query",
                       "TotalForward", "TotalReverse", "fG", "fA", "fT", "fC", "fN", "rG", "rA", "rT", "rC", "rN",
                       "fDel", "rDel", "fIns", "rIns"]

        vcf_data_df = pandas.DataFrame(columns=vcf_columns)

        for i in range(region_length):
            position_count += 1
            snpid = "."
            qread_depth = len(position_data_df.loc["depth", column_names[i]])

            # Because of filtering this can happen.
            if qread_depth == 0:
                continue

            query_orientation = Counter(position_data_df.loc["depth", column_names[i]])

            # Remove any position where there are not enough filtered reads.
            if query_orientation["forward"] <= int(self._args.Minimum_Read_Depth) or \
                    query_orientation["reverse"] <= int(self._args.Minimum_Read_Depth):
                continue

            contig = column_names[i].split(":")[0]
            position = int(column_names[i].split(":")[1])
            ref_base = self._refseq.fetch(reference=contig[3:], start=position, end=position + 1)

            # Identify SNP's in query.
            for record in self.vcf_snp.fetch(contig[3:], position - 2, position):
                snp_alt = []
                for snp in record.ALT:
                    snp_alt.append(str(snp))

                if len(record.ID) > 1 and set(map(tuple, snp_alt)) & \
                        set(map(tuple, position_data_df.loc["snv", column_names[i]][1])):
                    snpid = record.ID

            # Identify and Count SNV's
            qsnv_data_dict = None
            if position_data_df.loc["snv", column_names[i]][1]:
                qsnv_data_dict = snv_processor(position_data_df.loc["snv", column_names[i]][1],
                                               position_data_df.loc["snv", column_names[i]][0],
                                               min_allele_count)

            # Filter Deletions
            if len(position_data_df.loc["del", column_names[i]][0]) > min_allele_count:
                # if the deletion is the first position I ignore it.
                exception = False
                try:
                    refbase = vcf_data_df.loc[position-1, "REF"]
                    alt_list = vcf_data_df.loc[position-1, "ALT"].split(",")
                except KeyError:
                    exception = True

                if not exception:
                    del_columns = ["REF", "ALT", "fDel", "rDel"]
                    query_orientation_del_count = Counter(position_data_df.loc["del", column_names[i]][0])
                    del_qstrands = (query_orientation_del_count["forward"], query_orientation_del_count["reverse"])

                    for x in range(len(alt_list)):
                        alt_list[x] = "{}{}".format(alt_list[x], ref_base)

                    alt_list.append(refbase)
                    alt_string = ",".join(alt_list)
                    refbase += ref_base

                    vcf_data_df.loc[position-1, del_columns] = [refbase, alt_string, del_qstrands[0], del_qstrands[1]]

            # Filter and Format insertions
            variant_list, ins_qstrands = self.insertions(position_data_df, column_names[i], ref_base, position)

            # Filter SNV's
            snv_count = 0
            snv_forward_counts = [0, 0, 0, 0, 0]
            snv_reverse_counts = [0, 0, 0, 0, 0]

            if qsnv_data_dict and qsnv_data_dict["G"][0] > min_allele_count:
                variant_list += "G"
                snv_forward_counts[0], snv_reverse_counts[0] = strand_counter(qsnv_data_dict["G"][1])
                snv_count += qsnv_data_dict["G"][0]

            if qsnv_data_dict and qsnv_data_dict["A"][0] > min_allele_count:
                variant_list += "A"
                snv_count += qsnv_data_dict["A"][0]
                snv_forward_counts[1], snv_reverse_counts[1] = strand_counter(qsnv_data_dict["A"][1])

            if qsnv_data_dict and qsnv_data_dict["T"][0] > min_allele_count:
                variant_list += "T"
                snv_count += qsnv_data_dict["T"][0]
                snv_forward_counts[2], snv_reverse_counts[2] = strand_counter(qsnv_data_dict["T"][1])

            if qsnv_data_dict and qsnv_data_dict["C"][0] > min_allele_count:
                variant_list += "C"
                snv_count += qsnv_data_dict["C"][0]
                snv_forward_counts[3], snv_reverse_counts[3] = strand_counter(qsnv_data_dict["C"][1])

            if qsnv_data_dict and qsnv_data_dict["N"][0] > min_allele_count:
                snv_forward_counts[4], snv_reverse_counts[4] = strand_counter(qsnv_data_dict["N"][1])

            if variant_list:
                snv_string = ",".join(variant_list)
                variant_count_list = [0]*10
                if snv_count > 0:
                    # Forward G, A, T, C, N; Reverse G, A, T, C, N
                    variant_count_list = [snv_forward_counts[0], snv_forward_counts[1], snv_forward_counts[2],
                                          snv_forward_counts[3], snv_forward_counts[4], snv_reverse_counts[0],
                                          snv_reverse_counts[1], snv_reverse_counts[2], snv_reverse_counts[3],
                                          snv_reverse_counts[4]]

                # Forward/Reverse Deletions; Forward/Reverse Insertions
                indel_list = [0, 0, ins_qstrands[0], ins_qstrands[1]]
                # Python uses 0 based numbering while most bioinformatic packages use 1 based numbering.
                d_list = [contig, position+1, snpid, ref_base, snv_string, ".", ".", ".", ".",
                          gene, qread_depth, query_orientation["forward"], query_orientation["reverse"]]
                d_list.extend(variant_count_list)
                d_list.extend(indel_list)
                vcf_data_df.loc[position, vcf_columns] = d_list

        return vcf_data_df


def variant_block_writer(args, labels, vcf_data_df):
    region_id = labels[2]
    index_list = vcf_data_df.index.values
    column_names = vcf_data_df.columns.tolist()
    outstring = ""

    for idx in index_list:
        line_list = []
        for colmn in column_names:
            line_list.append(str(vcf_data_df.loc[idx, colmn]))
        line = "\t".join(line_list)
        outstring += "{}\n".format(line)

    outfile = open("{}{}_{}_variants.txt".format(args.WorkingFolder, args.Job_Name, region_id), "w")
    outfile.write(outstring)
    outfile.close()


class VCFfilter:
    """
    Applies error correction matrix filtering to raw VCF file.
    """

    def __init__(self, args, log):
        self.args = args
        self.log = log
        self.refseq = pysam.FastaFile(self.args.RefSeq)
        self.tumor_sample = None
        self.column_names = None
        self.vcf_header = ""
        self.ks_dataout = ""
        self.forward_columns = ["fG", "fA", "fT", "fC", "fDel", "fIns"]
        self.reverse_columns = ["rG", "rA", "rT", "rC", "rDel", "rIns"]
        self.error_columns = ["G_dcMAF", "A_dcMAF", "T_dcMAF", "C_dcMAF", "Del_dcMAF", "Ins_dcMAF"]

        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        os.chdir("..")
        pickle_file = "{0}{1}pickles{1}error_correction_df.pkl".format(os.getcwd(), os.sep)
        if not os.path.isfile(pickle_file):
            log.error("Error Correction Matrix Not Found.  Have You Created the Error Matrix?")
            return
        with open(pickle_file, 'rb') as file:
            self.error_data_df = dill.load(file)

        pickle_file2 = "{0}{1}pickles{1}KS_error_correction_df.pkl".format(os.getcwd(), os.sep)
        if not os.path.isfile(pickle_file2):
            log.error("Error Correction Matrix Not Found.  Have You Created the Error Matrix?")
            return
        with open(pickle_file2, 'rb') as file:
            self.ks_error_data_df = dill.load(file)

    def dataframe_generator(self, vcf_input_file):
        data_vcf_dict = collections.OrderedDict()
        column_header = True
        line_count = 0
        with open(vcf_input_file) as f:
            for l in f:
                # Break each line into a list.
                l_list = [x for x in l.strip("\n").split("\t")]

                if not column_header:
                    # This is the VCF data.  Key is chr:position.
                    key = "{}:{}".format(l_list[0], l_list[1])
                    data_vcf_dict[key] = l_list
                elif column_header:
                    if l_list[0] == "#CHROM":
                        # These are the column names and define the end of the header section.
                        column_header = False
                        self.column_names = l_list
                    elif line_count == 3:
                        # Get the sample name from the VCF file.
                        sample_name = l_list[0].split("=")[1]

                    elif line_count != 1 or line_count != 4:
                        # Line 1 is the date, line 4 is the source.  This builds a header without those or the sample
                        self.vcf_header += "{}\n".format(l_list[0])
                line_count += 1

        # Convert dictionary to dataframe and transpose
        vcf_df = pandas.DataFrame(data_vcf_dict, index=self.column_names)
        return vcf_df.T, sample_name

    @staticmethod
    def triplet_key_generator(df_row, variant_df):
        # ToDo: confirm that this function is used here.
        if df_row.ALT != ".":
            back_index = "{}:{}".format(df_row._1, int(df_row.POS) - 1)
            forward_index = "{}:{}".format(df_row._1, int(df_row.POS) + 1)
            try:
                back_ref = variant_df.loc[back_index, "REF"]
                forward_ref = variant_df.loc[forward_index, "REF"]
                back_chrom = variant_df.loc[back_index, "#CHROM"]
                forward_chrom = variant_df.loc[forward_index, "#CHROM"]
            except KeyError:
                return False
            if back_chrom != forward_chrom:
                return False

            ref = df_row.REF
            if len(df_row.REF) == 3:
                ref = df_row.REF[1]
            if len(back_ref) == 3:
                back_ref = back_ref[1]
            if len(forward_ref) == 3:
                forward_ref = forward_ref[1]

            return "{}{}{}".format(back_ref, ref, forward_ref)

    def filter_vcf_file(self):
        """
        Main entry point for filtering class.
        :return:
        """
        def mutant_allele_frequency(forward_count, reverse_count, strand_counts):
            """
            Called once for each G A T C Del and Ins count.  I am leaving this here for future additional filters.
            :param forward_count:
            :param reverse_count:
            :param strand_counts:
            :return:
            """
            forward_strand_counts = strand_counts[0]
            reverse_strand_counts = strand_counts[1]
            raw_forward_freq = 0
            raw_reverse_freq = 0

            # Check for minimum alt allele counts.
            if forward_count <= int(self.args.Minimum_Allele_Count) or \
                    reverse_count <= int(self.args.Minimum_Allele_Count):
                return 0, 0

            # Check for minimum read depth
            if forward_strand_counts <= int(self.args.Minimum_Read_Depth) or \
                    reverse_strand_counts <= int(self.args.Minimum_Read_Depth):
                return 0, 0

            # Determine our raw mutant allele fraction.
            with suppress(ZeroDivisionError):
                raw_forward_freq = forward_count/forward_strand_counts
            with suppress(ZeroDivisionError):
                raw_reverse_freq = reverse_count/reverse_strand_counts

            adjusted_forward_freq = 0
            adjusted_reverse_freq = 0

            # Adjust raw MAF based on triplet error frequency.  Currently this is not used.  Triplet freq is 0.
            if raw_forward_freq > forward_triplet_frequency:
                adjusted_forward_freq = raw_forward_freq-forward_triplet_frequency
            if raw_reverse_freq > reverse_triplet_frequency:
                adjusted_reverse_freq = raw_reverse_freq-reverse_triplet_frequency

            return adjusted_forward_freq, adjusted_reverse_freq

        self.log.info("Reading Tumor VCF File")
        tumor_vcf_df, self.tumor_sample = self.dataframe_generator(self.args.Tumor_VCF)
        tumor_freq_list = []
        ks_freq_list = []
        dataout = ""
        row_count = 0
        for row in tumor_vcf_df.itertuples():
            row_count += 1
            if row_count % 10000 == 0:
                self.log.info("{} rows completed".format(row_count))

            # Only analyze positions that are in the tumor and error correction dataframes.
            try:
                self.error_data_df[row.Index]
            except KeyError:
                continue

            tumor_strand_counts = (int(tumor_vcf_df.loc[row.Index, "TotalForward"]),
                                   int(tumor_vcf_df.loc[row.Index, "TotalReverse"]))

            # There are 12 frequencies for each row. Make sure we are starting with an empty list.
            tumor_freq_list.clear()
            ks_freq_list.clear()
            for forward_column, reverse_column, error_column in zip(self.forward_columns, self.reverse_columns,
                                                                    self.error_columns):
                # Currently triplet frequencies are not used.
                # forward_triplet_frequency = statistics.mean(self.triplet_data_df.loc[triplet_key, forward_column])
                forward_triplet_frequency = 0
                reverse_triplet_frequency = 0
                # reverse_triplet_frequency = statistics.mean(self.triplet_data_df.loc[triplet_key, reverse_column])

                freq_forward_adjusted, freq_reverse_adjusted = \
                    mutant_allele_frequency(int(tumor_vcf_df.loc[row.Index, forward_column]),
                                            int(tumor_vcf_df.loc[row.Index, reverse_column]),
                                            tumor_strand_counts)

                # Get the error from the appropriate source
                error_freq, error_stdev = self.error_data_df[row.Index][error_column]

                expected_error = self.ks_error_data_df[row.Index][error_column]

                ks_freq_list.append(self.ks_allele_frequency(freq_forward_adjusted, freq_reverse_adjusted,
                                                             expected_error))

                # Get our filtered, allele frequency.
                tumor_freq_list.append(self.allele_frequency(freq_forward_adjusted, freq_reverse_adjusted,
                                                             error_freq, error_stdev))

            self.ks_output(tumor_vcf_df, row, ks_freq_list, tumor_strand_counts)

            alt_allele = []
            del_found = False
            alt_found = False

            position_frequencies = ["0", "0", "0", "0", "0", "0"]
            alt_list = tumor_vcf_df.loc[row.Index, "ALT"].split(",")
            rbase = tumor_vcf_df.loc[row.Index, "REF"]
            v_list = ["G", "A", "T", "C", "Del", "Ins"]

            if tumor_freq_list[4] > 0:
                del_found = True

            for i in reversed(range(len(v_list))):
                allel_freq = tumor_freq_list[i]
                v = v_list[i]
                if allel_freq > 0:
                    position_frequencies[i] = str(allel_freq)

                    if v == "Ins":
                        for item in alt_list:
                            if len(item) > len(rbase):
                                alt_allele.append(item)

                    elif v == "Del":
                        for item in alt_list:
                            if len(item) < len(rbase):
                                del_found = True
                                alt_allele.append(item)

                    elif v != "Del" or v != "Ins":
                        if del_found:
                            alt_allele.append("{}{}".format(v, rbase[1:]))
                        else:
                            alt_allele.append(v)

                    alt_found = True

            filtered_freq_string = "\t".join(position_frequencies)
            alt = ",".join(alt_allele)

            filter_pass = self.strand_filter(tumor_strand_counts[0], tumor_strand_counts[1])
            row_data = ""

            for column in self.column_names:
                data = tumor_vcf_df.loc[row.Index, column]
                if column == "FILTER" and not filter_pass:
                    row_data += "FAIL\t"
                elif column == "ALT" and alt_found:
                    row_data += "{}\t".format(alt)
                elif column == "ALT" and not alt_found:
                    row_data += ".\t"

                elif column == "ID" and not alt_found:
                    row_data += ".\t"
                elif column == "REF":
                    # ctg = row.Index.split(":")[0]
                    # position = int(row.Index.split(":")[1])
                    # ref_base = self.refseq.fetch(reference=ctg[3:], start=position-1, end=position)
                    # row_data += "{}\t".format(ref_base)
                    if del_found:
                        row_data += "{}\t".format(rbase)
                    else:
                        ctg = row.Index.split(":")[0]
                        position = int(row.Index.split(":")[1])
                        row_data += "{}\t".format(self.refseq.fetch(reference=ctg[3:], start=position-1, end=position))
                else:
                    row_data += "{}\t".format(data)

            row_data += "{}\n".format(filtered_freq_string)

            # Allows user to choose to output all date, including normal positions or just alt positions.
            if alt_found or self.args.Include_All:
                dataout += row_data

        self.write_vcf(dataout)
        self.write_ks_vcf(self.ks_dataout)

    @staticmethod
    def strand_filter(forward, reverse):
        """

        :param forward:
        :param reverse:
        :return:
        """
        if reverse == 0:
            return False
        strand_ratio = forward / reverse
        if 0.25 <= strand_ratio <= 4:
            return True

        return False

    def ks_output(self, tumor_vcf_df, row, ks_freq_list, tumor_strand_counts):
        """

        :param tumor_vcf_df:
        :param row:
        :param ks_freq_list:
        :param tumor_strand_counts:
        """
        alt_allele = []
        del_found = False
        alt_found = False

        position_frequencies = ["0", "0", "0", "0", "0", "0"]
        alt_list = tumor_vcf_df.loc[row.Index, "ALT"].split(",")
        rbase = tumor_vcf_df.loc[row.Index, "REF"]
        v_list = ["G", "A", "T", "C", "Del", "Ins"]

        if ks_freq_list[4] > 0:
            del_found = True

        for i in reversed(range(len(v_list))):
            allel_freq = ks_freq_list[i]
            v = v_list[i]
            if allel_freq > 0:
                position_frequencies[i] = str(allel_freq)

                if v == "Ins":
                    for item in alt_list:
                        if len(item) > len(rbase):
                            alt_allele.append(item)

                elif v == "Del":
                    for item in alt_list:
                        if len(item) < len(rbase):
                            del_found = True
                            alt_allele.append(item)

                elif v != "Del" or v != "Ins":
                    if del_found:
                        alt_allele.append("{}{}".format(v, rbase[1:]))
                    else:
                        alt_allele.append(v)

                alt_found = True

        filtered_freq_string = "\t".join(position_frequencies)
        alt = ",".join(alt_allele)

        filter_pass = self.strand_filter(tumor_strand_counts[0], tumor_strand_counts[1])
        row_data = ""

        for column in self.column_names:
            data = tumor_vcf_df.loc[row.Index, column]
            if column == "FILTER" and not filter_pass:
                row_data += "FAIL\t"
            elif column == "ALT" and alt_found:
                row_data += "{}\t".format(alt)
            elif column == "ALT" and not alt_found:
                row_data += ".\t"

            elif column == "ID" and not alt_found:
                row_data += ".\t"
            elif column == "REF":
                if del_found:
                    row_data += "{}\t".format(rbase)
                else:
                    ctg = row.Index.split(":")[0]
                    position = int(row.Index.split(":")[1])
                    row_data += "{}\t".format(self.refseq.fetch(reference=ctg[3:], start=position - 1, end=position))
            else:
                row_data += "{}\t".format(data)

        row_data += "{}\n".format(filtered_freq_string)

        # Allows user to choose to output all date, including normal positions or just alt positions.
        if alt_found or self.args.Include_All:
            self.ks_dataout += row_data

    def ks_allele_frequency(self, forward_freq, reverse_freq, expected_error):
        """

        :param forward_freq:
        :param reverse_freq:
        :param expected_error:
        :return:
        """
        try:
            allele_freq_ratio = forward_freq / reverse_freq
        except ZeroDivisionError:
            allele_freq_ratio = 0

        freq = 0
        if allele_freq_ratio < 0.25 or allele_freq_ratio > 4:
            freq = min(forward_freq, reverse_freq)
        elif 0.25 <= allele_freq_ratio <= 4:
            freq = gmean([forward_freq, reverse_freq])
        else:
            Tool_Box.debug_messenger("No duplex frequency generated")

        # If the observed frequency is below the lower limit for Odin then return a 0.
        if freq < float(self.args.Minimum_Allele_Freq):
            return 0

        # Do a Kolmogorov-Smirnov test
        ks_result = scipy.stats.kstest(expected_error, scipy.stats.norm.cdf, args=(freq,))
        if ks_result[1] < 0.00015:
            return freq
        return 0

    def allele_frequency(self, forward_freq, reverse_freq, error_freq, error_stdev):
        """

        :param forward_freq:
        :param reverse_freq:
        :param error_freq:
        :param error_stdev:
        :return:
        """
        try:
            allele_freq_ratio = forward_freq / reverse_freq
        except ZeroDivisionError:
            allele_freq_ratio = 0

        freq = 0
        if allele_freq_ratio < 0.25 or allele_freq_ratio > 4:
            freq = min(forward_freq, reverse_freq)
        elif 0.25 <= allele_freq_ratio <= 4:
            freq = gmean([forward_freq, reverse_freq])
        else:
            self.log.error("No duplex frequency generated")

        # If the observed frequency is below the lower limit for Odin then return a 0.
        if freq <= float(self.args.Minimum_Allele_Freq):
            freq = 0
        elif freq <= error_freq+(error_stdev*float(self.args.Min_Fold_Increase)):
            freq = 0
        # Checks if the error reported in the error correction matrix is likely from SNP's
        # if error_freq <= float(self.args.SNP_Error_Freq):
        #     # Check if the observed ALT freq is greater than the minimum fold increase re the error matrix.
        #     if freq <= error_freq*float(self.args.Min_Fold_Increase):
        #         freq = 0
        #
        # elif freq <= error_freq*0.9:
        #     # ToDo: The 0.9 factor here is arbitrary.  Should use a value based on some logic.
        #     # if the error correction value is from a SNP it is possible our observed value is real if it is greater
        #     freq = 0

        return freq

    def write_vcf(self, dataout):
        """

        :param dataout:
        """
        file_header = "##fileformat=VCFv4.3\n##filedate={}\n"\
            .format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        file_header += "##tumorsample={}\n##tumorVCFfile={}\n"\
            .format(self.tumor_sample, self.args.Tumor_VCF.split("/")[-1])
        file_header += "{}\n".format(self.vcf_header)
        file_header += "\t" .join(self.column_names)
        file_header += "\tG_dcMAF\tA_dcMAF\tT_dcMAF\tC_dcMAF\tDel_dcMAF\tIns_dcMAF\n"

        data_string = "{}{}".format(file_header, dataout)
        outfile = open("{}{}_filtered_variants.txt".format(self.args.WorkingFolder, self.args.Job_Name), 'w')
        outfile.write(data_string)
        outfile.close()

    def write_ks_vcf(self, dataout):
        """

        :param dataout:
        """
        file_header = "##fileformat=VCFv4.3\n##filedate={}\n"\
            .format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        file_header += "##tumorsample={}\n##tumorVCFfile={}\n"\
            .format(self.tumor_sample, self.args.Tumor_VCF.split("/")[-1])
        file_header += "{}\n".format(self.vcf_header)
        file_header += "\t" .join(self.column_names)
        file_header += "\tG_dcMAF\tA_dcMAF\tT_dcMAF\tC_dcMAF\tDel_dcMAF\tIns_dcMAF\n"

        data_string = "{}{}".format(file_header, dataout)
        outfile = open("{}{}_KS_filtered_variants.txt".format(self.args.WorkingFolder, self.args.Job_Name), 'w')
        outfile.write(data_string)
        outfile.close()


class TripletFiltering:
    """
    This uses the triplet context of errors to filter
    """
    def __init__(self, args, log):
        self.args = args
        self.log = log
        self.control_sample = None
        self.tumor_sample = None
        self.column_names = None
        self.vcf_header = ""
        self.forward_columns = ["fG", "fA", "fT", "fC", "fDel", "fIns"]
        self.reverse_columns = ["rG", "rA", "rT", "rC", "rDel", "rIns"]
        # dir_path = os.path.dirname(os.path.abspath(__file__))
        # pickle_file = "{}{}error_correction_df.pkl".format(dir_path, os.sep)
        with open('triplet_data_df.pkl', 'rb') as file:
            self.triplet_data_df = dill.load(file)

    def dataframe_generator(self, vcf_input_file):
        """

        :param vcf_input_file:
        :return:
        """
        data_vcf_dict = collections.OrderedDict()
        column_header = False
        line_count = 0
        with open(vcf_input_file) as f:
            for l in f:
                l_list = [x for x in l.strip("\n").split("\t")]

                if column_header:
                    key = "{}:{}".format(l_list[0], l_list[1])
                    data_vcf_dict[key] = l_list
                elif not column_header:
                    if l_list[0] == "#CHROM":
                        column_header = True
                        self.column_names = l_list
                    elif line_count == 3:
                        sample_name = l_list[0].split("=")[1]
                    elif line_count != 1 or line_count != 4:
                        self.vcf_header += "{}\n".format(l_list[0])
                line_count += 1

        # Convert dictionary to dataframe and transpose
        vcf_df = pandas.DataFrame(data_vcf_dict, index=self.column_names)
        return vcf_df.T, sample_name

    @staticmethod
    def triplet_key_generator(df_row, variant_df):
        """

        :param df_row:
        :param variant_df:
        :return:
        """
        if df_row.ALT != ".":
            back_index = "{}:{}".format(df_row._1, int(df_row.POS) - 1)
            forward_index = "{}:{}".format(df_row._1, int(df_row.POS) + 1)
            try:
                back_ref = variant_df.loc[back_index, "REF"]
                forward_ref = variant_df.loc[forward_index, "REF"]
                back_chrom = variant_df.loc[back_index, "#CHROM"]
                forward_chrom = variant_df.loc[forward_index, "#CHROM"]
            except KeyError:
                return False
            if back_chrom != forward_chrom:
                return False

            ref = df_row.REF
            if len(df_row.REF) == 3:
                ref = df_row.REF[1]
            if len(back_ref) == 3:
                back_ref = back_ref[1]
            if len(forward_ref) == 3:
                forward_ref = forward_ref[1]

            return "{}{}{}".format(back_ref, ref, forward_ref)

    def filter_vcf_file(self):
        """

        :return:
        """
        def triplet_frequency_adjustment(forward_count, reverse_count, strand_counts):
            forward_strand_counts = strand_counts[0]
            reverse_strand_counts = strand_counts[1]
            raw_forward_freq = 0
            raw_reverse_freq = 0

            with suppress(ZeroDivisionError):
                raw_forward_freq = forward_count/forward_strand_counts
            with suppress(ZeroDivisionError):
                raw_reverse_freq = reverse_count/reverse_strand_counts

            adjusted_forward_freq = 0
            adjusted_reverse_freq = 0
            if raw_forward_freq > forward_triplet_frequency:
                adjusted_forward_freq = raw_forward_freq-forward_triplet_frequency
            if raw_reverse_freq > reverse_triplet_frequency:
                adjusted_reverse_freq = raw_reverse_freq-reverse_triplet_frequency

            return adjusted_forward_freq, adjusted_reverse_freq

        self.log.info("Reading Control VCF File")
        control_vcf_df, self.control_sample = self.dataframe_generator(self.args.Control_VCF)

        self.log.info("Reading Tumor VCF File")
        tumor_vcf_df, self.tumor_sample = self.dataframe_generator(self.args.Tumor_VCF)
        tumor_freq_list = []
        control_freq_list = []
        dataout = ""
        row_count = 0
        for row in tumor_vcf_df.itertuples():
            row_count += 1
            if row_count % 10000 == 0:
                self.log.info("{} rows completed".format(row_count))

            # Only analyze positions that are in the tumor and control and we can get a triplet key.
            try:
                control_strand_counts = (
                    int(control_vcf_df.loc[row.Index, "TotalForward"]),
                    int(control_vcf_df.loc[row.Index, "TotalReverse"]))
            except KeyError:
                continue
            triplet_key = self.triplet_key_generator(row, tumor_vcf_df)
            if not triplet_key:
                continue

            tumor_strand_counts = (int(tumor_vcf_df.loc[row.Index, "TotalForward"]),
                                   int(tumor_vcf_df.loc[row.Index, "TotalReverse"]))

            # There are 12 frequencies for each row.
            control_freq_list.clear()
            tumor_freq_list.clear()
            for forward_column, reverse_column in zip(self.forward_columns, self.reverse_columns):
                # forward_triplet_frequency = statistics.mean(self.triplet_data_df.loc[triplet_key, forward_column])
                forward_triplet_frequency = 0
                reverse_triplet_frequency = 0
                # reverse_triplet_frequency = statistics.mean(self.triplet_data_df.loc[triplet_key, reverse_column])

                control_forward_adjusted, control_reverse_adjusted = \
                    triplet_frequency_adjustment(int(control_vcf_df.loc[row.Index, forward_column]),
                                                 int(control_vcf_df.loc[row.Index, reverse_column]),
                                                 control_strand_counts)
                tumor_forward_adjusted, tumor_reverse_adjusted = \
                    triplet_frequency_adjustment(int(tumor_vcf_df.loc[row.Index, forward_column]),
                                                 int(tumor_vcf_df.loc[row.Index, reverse_column]),
                                                 tumor_strand_counts)

                control_freq_list.append(self.allele_frequency(control_forward_adjusted, control_reverse_adjusted))

                tumor_freq_list.append(self.allele_frequency(tumor_forward_adjusted, tumor_reverse_adjusted))

            alt_allele = []
            alt_found = False
            filter_pass = True
            position_frequencies = ["0", "0", "0", "0", "0", "0"]

            for i in range(len(["G", "A", "T", "C", "Del", "Ins"])):
                # c = frequency_data_dict["Control"][i]
                c = control_freq_list[i]
                t = tumor_freq_list[i]
                # t = frequency_data_dict["Tumor"][i]
                v = ["G", "A", "T", "C", "Del", "Ins"][i]
                if t > 0 and t >= c * int(self.args.Min_Fold_Increase):
                    position_frequencies[i] = str(t)

                    if v == "Ins":
                        ins_seq = "<{}>".format(tumor_vcf_df.loc[row.Index, "ALT"].split("<")[1].split(">")[0])
                        alt_allele.append(ins_seq)
                    elif v != "Del":
                        alt_allele.append(v)

                    alt_found = True

            filtered_freq_string = "\t".join(position_frequencies)
            alt = ",".join(alt_allele)

            if not self.strand_filter(tumor_strand_counts[0], tumor_strand_counts[1]):
                filter_pass = False

            for column in self.column_names:
                data = tumor_vcf_df.loc[row.Index, column]
                if column == "FILTER" and not filter_pass:
                    dataout += "FAIL\t"
                elif v == "ALT" and alt_found:
                    dataout += "{}\t".format(alt)
                elif v == "ALT" and not alt_found:
                    dataout += ".\t"
                else:
                    dataout += "{}\t".format(data)

            dataout += "{}\n".format(filtered_freq_string)
        self.write_vcf(dataout)

    @staticmethod
    def strand_filter(forward, reverse):
        """

        :param forward:
        :param reverse:
        :return:
        """
        if reverse == 0:
            return False
        strand_ratio = forward / reverse
        if 0.25 <= strand_ratio <= 4:
            return True

        return False

    @staticmethod
    def allele_frequency(forward_variant_count, reverse_variant_count):
        """

        :param forward_variant_count:
        :param reverse_variant_count:

        :return:
        """
        # forward_counts = strand_counts[0]
        # reverse_counts = strand_counts[1]

        # if forward_counts == 0 or forward_variant_count <= min_allele_count:
        #     forward_freq = 0
        # else:
        #     forward_freq = forward_variant_count / forward_counts
        #
        # if reverse_counts == 0 or reverse_variant_count <= min_allele_count:
        #     reverse_freq = 0
        #     # strand_ratio = 0
        #
        # else:
        #     # strand_ratio = forward_counts/reverse_counts
        #     reverse_freq = reverse_variant_count / reverse_counts

        forward_freq = forward_variant_count
        reverse_freq = reverse_variant_count
        try:
            maf_strand_ratio = forward_freq / reverse_freq
        except ZeroDivisionError:
            maf_strand_ratio = 0
        freq = 0
        if 0 < maf_strand_ratio < 0.25:
            # freq = forward_variant_count / sum(strand_counts)
            freq = min(forward_freq, reverse_freq)
        elif 0.25 <= maf_strand_ratio <= 4:
            # freq = (forward_variant_count+reverse_variant_count)/sum(strand_counts)
            # freq = (forward_freq + reverse_freq) / 2
            freq = gmean([forward_freq, reverse_freq])
        elif maf_strand_ratio > 4:
            # freq = reverse_variant_count / sum(strand_counts)
            freq = min(forward_freq, reverse_freq)

        return freq

    def write_vcf(self, dataout):
        """
        Writes our VCF file
        :param dataout:
        """
        file_header = "##filedate={}\n".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        file_header += "##controlsample={}\n##controlVCFfile={}\n".format(self.control_sample, self.args.Control_VCF.split("/")[-1])
        file_header += "##tumorsample={}\n##tumorVCFfile={}\n".format(self.tumor_sample, self.args.Tumor_VCF.split("/")[-1])
        file_header += "{}\n".format(self.vcf_header)
        file_header += "\t" .join(self.column_names)
        file_header += "\tG_dcMAF\tA_dcMAF\tT_dcMAF\tC_dcMAF\tDel_dcMAF\tIns_dcMAF\n"

        data_string = "{}{}".format(file_header, dataout)
        outfile = open("{}{}_filtered_variants.txt".format(self.args.WorkingFolder, self.args.Job_Name), 'w')
        outfile.write(data_string)
        outfile.close()
