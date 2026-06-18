"""
March 22, 2019
    Version 0.4.0
    Dennis A. Simpson
    Per position error matrix now includes global error for each type and standard deviations.  Both in the output file
    and in the pickle file.

This library is intended to produce a Matrix of mutation types

@author: Dennis A. Simpson
         Asystbio Laboratories, LLC
         Durham, NC
@copyright: 2019
"""
import collections
import datetime
import os
from scipy.stats import gmean
import statistics
import pandas
import dill
from natsort import natsort
from Valkyries import Tool_Box

__author__ = 'Dennis A. Simpson'
__version__ = "0.4.0"
__package__ = 'Odin'


class ErrorMatrix:
    def __init__(self, args, log):
        self.args = args
        self.log = log
        self.file_names = []
        self.data_vcf_dict, self.file_count = self.vcf_file_parsing()
        self.outstring = ""
        self.error_correction_df = None

        self.error_correction_snp_dict = collections.defaultdict(lambda: collections.defaultdict(float))

        if os.path.isfile(getattr(args, "SNP_Table", None)):
            self.log.info("Reading Error Correction SNP File")
            with open(self.args.SNP_Table) as f:
                first_line = True
                for l in f:
                    mylist = l.split("\t")
                    if first_line:
                        first_line = False
                        continue
                    if mylist[5] == "G":
                        column_key = "G_dcMAF"
                    elif mylist[5] == "A":
                        column_key = "A_dcMAF"
                    elif mylist[5] == "T":
                        column_key = "T_dcMAF"
                    elif mylist[5] == "C":
                        column_key = "C_dcMAF"

                    self.error_correction_snp_dict["{}:{}".format(mylist[0], mylist[1])][column_key] = 0.0000

    def vcf_file_parsing(self):
        # Gather names of remaining text files and define some objects.
        vcf_files_path = [x for x in os.listdir(self.args.Working_Folder) if x.endswith(".txt")]
        data_vcf_dict = collections.defaultdict(lambda: collections.defaultdict(list))
        column_name_dict = collections.defaultdict()

        for vcf_file_name in vcf_files_path:
            column_header = False
            self.file_names.append(vcf_file_name)
            self.log.info("Reading {}".format(vcf_file_name))
            vcf_file = "{}{}".format(self.args.Working_Folder, vcf_file_name)
            column_names = []
            dc_names_list = ["G_dcMAF", "A_dcMAF", "T_dcMAF", "C_dcMAF", "Del_dcMAF", "Ins_dcMAF"]
            raw_counts_labels = [("fG", "rG"), ("fA", "rA"), ("fT", "rT"), ("fC", "rC"), ("fN", "rN"),
                                 ("fDel", "rDel"), ("fIns", "rIns")]
            # Read VCF file into dictionary
            with open(vcf_file) as f:
                for l in f:
                    l_list = [x for x in l.strip("\n").split("\t")]

                    if column_header:
                        chrom = l_list[0]
                        position = l_list[1]

                        for i, v in enumerate(l_list):
                            column_name_dict[column_names[i]] = v

                        key = "{}:{}".format(chrom, position)
                        data_vcf_dict[key]["GENE"] = column_name_dict["GENE"]
                        for c, dc_name in enumerate(dc_names_list):
                            alt_forward_count = raw_counts_labels[c][0]
                            alt_reverse_count = raw_counts_labels[c][1]
                            dc_freq = self.mutant_allele_frequency(int(column_name_dict["TotalForward"]),
                                                                   int(column_name_dict["TotalReverse"]),
                                                                   int(column_name_dict[alt_forward_count]),
                                                                   int(column_name_dict[alt_reverse_count]))

                            data_vcf_dict[key][dc_name].append(dc_freq)

                    elif l_list[0] == "#CHROM":
                        column_header = True
                        for name in l_list:
                            column_name_dict[name] = ""
                        column_names = l_list

        return data_vcf_dict, len(vcf_files_path)

    def process_data(self):
        self.log.info("Begin Processing Data")
        column_names = ["#CHROM", "POS", "GENE", "G_dcMAF", "A_dcMAF", "T_dcMAF", "C_dcMAF", "Del_dcMAF", "Ins_dcMAF"]
        count = 0
        avg_g_list = []
        avg_a_list = []
        avg_t_list = []
        avg_c_list = []
        avg_del_list = []
        avg_ins_list = []

        # Geometric means cannot be calculated for 0.
        gmean_zero_correction = 0.0000000001

        total_rows = len(self.data_vcf_dict)
        error_correction_dict = collections.OrderedDict()
        for key in natsort.natsorted(self.data_vcf_dict):
            if count % 10000 == 0:
                self.log.debug("{} rows processed of {}".format(count, total_rows))
            count += 1

            # Only score positions that are in at least 25% of the files.  Doesn't matter which we use to check.
            if len(self.data_vcf_dict[key]["G_dcMAF"]) <= 0.25*self.file_count:
                continue
            avg_g = statistics.mean(self.data_vcf_dict[key]["G_dcMAF"])
            avg_a = statistics.mean(self.data_vcf_dict[key]["A_dcMAF"])
            avg_t = statistics.mean(self.data_vcf_dict[key]["T_dcMAF"])
            avg_c = statistics.mean(self.data_vcf_dict[key]["C_dcMAF"])
            avg_del = statistics.mean(self.data_vcf_dict[key]["Del_dcMAF"])
            avg_ins = statistics.mean(self.data_vcf_dict[key]["Ins_dcMAF"])

            sd_g = statistics.pstdev(self.data_vcf_dict[key]["G_dcMAF"])
            sd_a = statistics.pstdev(self.data_vcf_dict[key]["A_dcMAF"])
            sd_t = statistics.pstdev(self.data_vcf_dict[key]["T_dcMAF"])
            sd_c = statistics.pstdev(self.data_vcf_dict[key]["C_dcMAF"])
            sd_del = statistics.pstdev(self.data_vcf_dict[key]["Del_dcMAF"])
            sd_ins = statistics.pstdev(self.data_vcf_dict[key]["Ins_dcMAF"])

            avg_g_list.append(avg_g+gmean_zero_correction)
            avg_a_list.append(avg_a+gmean_zero_correction)
            avg_t_list.append(avg_t+gmean_zero_correction)
            avg_c_list.append(avg_c+gmean_zero_correction)
            avg_del_list.append(avg_del+gmean_zero_correction)
            avg_ins_list.append(avg_ins+gmean_zero_correction)

            chrom = key.split(":")[0]
            position = int(key.split(":")[1])
            gene = self.data_vcf_dict[key]["GENE"]
            self.outstring += "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n"\
                .format(chrom, position, gene, avg_g, sd_g, avg_a, sd_a, avg_t, sd_t, avg_c, sd_c, avg_del, sd_del,
                        avg_ins, sd_ins)

            # error_correction_dict[key] = (chrom, position, gene, avg_g, avg_a, avg_t, avg_c, avg_del, avg_ins)
            error_correction_dict[key] = \
                [chrom, position, gene, (avg_g, sd_g), (avg_a, sd_a), (avg_t, sd_t), (avg_c, sd_c), (avg_del, sd_del),
                 (avg_ins, sd_ins)]

        total_avg_g = gmean(avg_g_list)-gmean_zero_correction
        total_avg_a = gmean(avg_a_list)-gmean_zero_correction
        total_avg_t = gmean(avg_t_list)-gmean_zero_correction
        total_avg_c = gmean(avg_c_list)-gmean_zero_correction
        total_avg_del = gmean(avg_del_list)-gmean_zero_correction
        total_avg_ins = gmean(avg_ins_list)-gmean_zero_correction
        stdev_total_avg_g = statistics.pstdev(avg_g_list)
        stdev_total_avg_a = statistics.pstdev(avg_a_list)
        stdev_total_avg_t = statistics.pstdev(avg_t_list)
        stdev_total_avg_c = statistics.pstdev(avg_c_list)
        stdev_total_avg_del = statistics.pstdev(avg_del_list)
        stdev_total_avg_ins = statistics.pstdev(avg_ins_list)

        error_correction_dict["Totals"] = \
            ["N.A.", "N.A.", "N.A.", (total_avg_g, stdev_total_avg_g), (total_avg_a, stdev_total_avg_a), (total_avg_t, stdev_total_avg_t),
             (total_avg_c, stdev_total_avg_c), (total_avg_del, stdev_total_avg_del),
             (total_avg_ins, stdev_total_avg_ins)]

        # Convert dictionary to dataframe and transpose
        self.error_correction_df = pandas.DataFrame(error_correction_dict, index=column_names)

        # Set the error for any position that contains SNP's to the average error for that type.
        for k1 in self.error_correction_snp_dict:
            for k2 in self.error_correction_snp_dict[k1]:
                if k2 == "G_dcMAF":
                    self.error_correction_df.at[k2, k1] = (total_avg_g, stdev_total_avg_g)
                elif k2 == "A_dcMAF":
                    self.error_correction_df.at[k2, k1] = (total_avg_a, stdev_total_avg_a)
                elif k2 == "T_dcMAF":
                    self.error_correction_df.at[k2, k1] = (total_avg_t, stdev_total_avg_t)
                elif k2 == "C_dcMAF":
                    self.error_correction_df.at[k2, k1] = (total_avg_t, stdev_total_avg_t)

    def data_output(self):
        # Pickle and save the variant dataframe in the pickle directory
        os.chdir(os.path.dirname(os.path.abspath(__file__)))
        os.chdir("..")

        pickle_file = "{0}{1}pickles{1}error_correction_df.pkl".format(os.getcwd(), os.sep)
        with open(pickle_file, 'wb') as file:
            dill.dump(self.error_correction_df, file, protocol=-1)

        # Pickle the dcMAF dictionary for use in error correction
        column_names = ["#CHROM", "POS", "GENE", "G_dcMAF", "A_dcMAF", "T_dcMAF", "C_dcMAF", "Del_dcMAF", "Ins_dcMAF"]
        ks_data_df = pandas.DataFrame(self.data_vcf_dict, index=column_names)
        pickle_file2 = "{0}{1}pickles{1}KS_error_correction_df.pkl".format(os.getcwd(), os.sep)
        with open(pickle_file2, 'wb') as file:
            dill.dump(ks_data_df, file, protocol=-1)

        # Now write this into a nice file.
        file_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        files = ",".join(self.file_names)
        group_labels = "#CHROM\tPOS\tGENE\tG_dcMAF\t\tA_dcMAF\t\tT_dcMAF\t\tC_dcMAF\t\tDel_dcMAF\t\tIns_dcMAF"
        column_labels = "\t\t\tAverage\tStandard Deviation\tAverage\tStandard Deviation\tAverage\tStandard Deviation\t" \
                        "Average\tStandard Deviation\tAverage\tStandard Deviation\tAverage\tStandard Deviation"
        header = "##Odin Error Correction File.\n" \
                 "##Generated with Mutation_Matrix v{} module ErrorMatrix.\n" \
                 "##File Date: {}\n" \
                 "##VCF Source Data File(s): {}\n\n{}\n{}\n"\
            .format(__version__, file_date, files, group_labels, column_labels)
        outfile = open("{}Odin_{}_Error_Correction_Data.bed".format(self.args.Working_Folder, self.args.JobName), 'w')
        outfile.write(header+self.outstring)
        outfile.close()

    def mutant_allele_frequency(self, forward_count, reverse_count, alt_forward_count, alt_reverse_count):
        """
        Called once for each G A T C Del and Ins count.  I am leaving this here for future additional filters.
        :param alt_reverse_count:
        :param alt_forward_count:
        :param forward_count:
        :param reverse_count:
        :return:
        """
        # Check for minimum alt allele counts.
        # if alt_forward_count <= int(self.args.Minimum_Allele_Count):
        #     if alt_reverse_count <= int(self.args.Minimum_Allele_Count):
        #         return 0

        try:
            raw_forward_freq = alt_forward_count / forward_count
        except ZeroDivisionError:
            raw_forward_freq = 0
        try:
            raw_reverse_freq = alt_reverse_count / reverse_count
        except ZeroDivisionError:
            raw_reverse_freq = 0
        try:
            allele_freq_ratio = raw_forward_freq / raw_reverse_freq
        except ZeroDivisionError:
            allele_freq_ratio = 0

        freq = 0
        if 0 < allele_freq_ratio < 0.25:
            freq = min(raw_forward_freq, raw_reverse_freq)
        elif 0.25 <= allele_freq_ratio <= 4:
            freq = gmean([raw_forward_freq, raw_reverse_freq])
        elif allele_freq_ratio > 4:
            freq = min(raw_forward_freq, raw_reverse_freq)

        return freq


class MutationMatrix:
    def __init__(self, args, log):
        self.args = args
        self.log = log
        self.triplet_data_df = self.triplet_matrix_dataframe()
        self.file_names = []
        self.forward_columns = ["fG", "fA", "fT", "fC", "fDel", "fIns"]
        self.reverse_columns = ["rG", "rA", "rT", "rC", "rDel", "rIns"]

    def variant_parser(self, vcf_file):

        column_header = False
        data_vcf_dict = collections.OrderedDict()
        self.file_names.append(vcf_file.split(os.sep)[-1])
        self.log.info("Reading {}".format(vcf_file.split(os.sep)[-1]))
        column_names = False

        # Read VCF file into dictionary
        with open(vcf_file) as f:
            for l in f:
                l_list = [x for x in l.strip("\n").split("\t")]

                if column_header:
                    key = "{}:{}".format(l_list[0], l_list[1])
                    data_vcf_dict[key] = l_list

                elif l_list[0] == "#Chr":
                    column_header = True
                    l_list[0] = "Chr"
                    column_names = l_list

        # Convert dictionary to dataframe and transpose
        if not column_names:
            return pandas.DataFrame()
        variant_df = pandas.DataFrame(data_vcf_dict, index=column_names)
        return variant_df.T

    @staticmethod
    def count_normalization(read_depth, variant_count):
        """
        Normalize variant counts
        :param read_depth:
        :param variant_count:
        :return:
        """
        try:
            # adjusted_count = (variant_count / read_depth) * 10000
            adjusted_count = variant_count / read_depth
        except ZeroDivisionError:
            adjusted_count = 0
        return adjusted_count

    def triplet_matrix_dataframe(self):
        self.log.debug("Triplet Matrix Dataframe")
        triplet_data_dict = collections.defaultdict(list)
        codons = ["G", "A", "T", "C"]
        column_names = ["fG", "fA", "fT", "fC", "fDel", "fIns", "rG", "rA", "rT", "rC", "rDel", "rIns"]
        for v1 in codons:
            for v2 in codons:
                for v3 in codons:
                    key = "{}{}{}".format(v1, v2, v3)
                    triplet_data_dict[key] = [[], [], [], [], [], [], [], [], [], [], [], []]

        # Convert dictionary to dataframe and transpose
        triplet_data_df = pandas.DataFrame(triplet_data_dict, index=column_names)
        triplet_data_dict.clear()
        return triplet_data_df.T

    def triplet_matrix_processing(self, vcf_file):
        variant_df = self.variant_parser("{}{}".format(self.args.Working_Folder, vcf_file))
        if variant_df.empty:
            return
        self.log.debug("Triplet Matrix Processing past variant_df")

        # Go through the variant dataframe by row
        for row in variant_df.itertuples():

            if row.ALT != ".":
                back_index = "{}:{}".format(row.Chr, int(row.POS)-1)
                forward_index = "{}:{}".format(row.Chr, int(row.POS)+1)
                try:
                    back_ref = variant_df.loc[back_index, "REF"]
                    forward_ref = variant_df.loc[forward_index, "REF"]
                    back_chrom = variant_df.loc[back_index, "Chr"]
                    forward_chrom = variant_df.loc[forward_index, "Chr"]
                except KeyError:
                    continue
                if back_chrom != forward_chrom:
                    continue
                ref = row.REF
                if len(row.REF) == 3:
                    ref = row.REF[1]
                if len(back_ref) == 3:
                    back_ref = back_ref[1]
                if len(forward_ref) == 3:
                    forward_ref = forward_ref[1]

                triplet_key = "{}{}{}".format(back_ref, ref, forward_ref)
                for forward_column, reverse_column in zip(self.forward_columns, self.reverse_columns):
                    forward_variant_count = int(variant_df.loc[row.Index, forward_column])
                    forward_depth = int(variant_df.loc[row.Index, "TotalForward"])
                    reverse_variant_count = int(variant_df.loc[row.Index, reverse_column])
                    reverse_depth = int(variant_df.loc[row.Index, "TotalReverse"])

                    forward_adjusted_count = self.count_normalization(forward_depth, forward_variant_count)
                    reverse_adjusted_count = self.count_normalization(reverse_depth, reverse_variant_count)

                    # This prevents SNP's from being included in results
                    if forward_adjusted_count > 0.25 or 0.25 < reverse_adjusted_count:
                        continue

                    self.triplet_data_df.loc[triplet_key, forward_column].append(forward_adjusted_count)
                    self.triplet_data_df.loc[triplet_key, reverse_column].append(reverse_adjusted_count)

    def write_output(self):
        outstring = ""
        for file in self.file_names:
            outstring += "#DataSource:{}\n".format(file)
        outstring += "\nTriplet\tfG\tfA\tfT\tfC\tfDel\tfIns\trG\trA\trT\trC\trDel\trIns\n"

        for row in self.triplet_data_df.itertuples():
            f = ""
            r = ""
            for forward_column, reverse_column in zip(self.forward_columns, self.reverse_columns):
                fmean = statistics.mean(self.triplet_data_df.loc[row.Index, forward_column])
                rmean = statistics.mean(self.triplet_data_df.loc[row.Index, reverse_column])

                f += "\t{}".format(fmean)
                r += "\t{}".format(rmean)

            outstring += "{}{}{}\n".format(row.Index, f, r)
        outfile = open("{}{}_triplet_variants.txt".format(self.args.Working_Folder, self.args.JobName), 'w')
        outfile.write(outstring)
        outfile.close()

        # Pickle and save the variant dataframe
        # self.triplet_data_df.to_pickle("triplet_data_df.pkl")
        with open('triplet_data_df.pkl', 'wb') as file:
            dill.dump(self.triplet_data_df, file, protocol=-1)
