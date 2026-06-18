"""


@author: Dennis A. Simpson
         RTP Genomics
         Chapel Hill, NC  27517
@copyright: 2020
"""

import os
import gzip
import collections
from Levenshtein import distance
import time
import natsort
from Valkyries import Tool_Box as Tool_Box

__author__ = 'Dennis A. Simpson'
__version__ = "0.20.0"


def batch_reader(repeater, fq1_batch):
    """
    This is debugging code used to test the multiprocessing pool setup.
    :param repeater:
    :param fq1_batch:
    :return:
    """
    pool_id = os.getpid()
    Tool_Box.debug_messenger([repeater, "Batch Reader", pool_id])
    for fq_read in fq1_batch:
        print(fq_read[0])
    time.sleep(1.5)
    print("Repeater: {}; PID {} Finished".format(repeater, pool_id))
    return


class FastqProcessing:
    def __init__(self, args, log, paired_end=False):
        self.log = log
        self.args = args
        self.fq1_batch = None
        self.fq2_batch = None
        self.fastq3_file = None
        self.paired_end = paired_end
        self.master_index_dict = {}
        self.sample_manifest_list = None
        self.sample_manifest_dictionary = collections.defaultdict(list)
        self.demultiplex_file_names = collections.defaultdict(list)
        self.summary_data = collections.defaultdict(list)

    def dataframe_build(self):
        """
        Build the dataframes containing the worker ID, indices, and file names.
        """
        self.log.info("Begin building dataframes")
        with open(self.args.MasterIndexFile) as f:
            for l in f:
                if "#" in l or not l:
                    continue
                l_list = [x for x in l.strip("\n").split("\t")]
                self.master_index_dict[l_list[0]] = [l_list[1], l_list[2]]

        self.sample_manifest_list = Tool_Box.FileParser.indices(self.log, self.args.SampleManifest)

        outfile_list_dict = collections.defaultdict(list)

        '''
        self.demultiplex_file_names["Unknown"] = \
            [Writer(self.log, "{}{}_Unknown_R1.fastq.gz".format(self.args.WorkingFolder, self.args.JobName), "Unknown"),
             Writer(self.log, "{}{}_Unknown_R2.fastq.gz".format(self.args.WorkingFolder, self.args.JobName), "Unknown")]

        self.demultiplex_file_names["GhostIndex"] = \
            [Writer(self.log, "{}{}_GhostIndex_R1.fastq.gz".format(self.args.WorkingFolder, self.args.JobName), "GhostIndex"),
             Writer(self.log, "{}{}_GhostIndex_R2.fastq.gz".format(self.args.WorkingFolder, self.args.JobName), "GhostIndex")]
        '''
        # Initialize our summary dataframe
        self.summary_data["Unknown"] = ["", "", 0]
        self.summary_data["GhostIndex"] = ["", "", 0]
        self.summary_data['total_reads'] = ["", "", 0]

        for i in range(int(self.args.Spawn)):
            outfile_list_dict["Unknown"].append(
                ["{}{}_Unknown_R1.fastq.gz".format(self.args.WorkingFolder, i),
                 "{}{}_Unknown_R2.fastq.gz".format(self.args.WorkingFolder, i)])

            '''
            I noticed with Nextera Libraries we occasionally found an index in a read that was valid but never used
            in the library prep.  I don't know, yet, how prevalent these Ghost Indices are so I am capturing them for
            later analysis.
            '''
            outfile_list_dict["GhostIndex"].append(
                ["{}{}_GhostIndex_R1.fastq.gz".format(self.args.WorkingFolder, i),
                 "{}{}_GhostIndex_R2.fastq.gz".format(self.args.WorkingFolder, i)])

            tmp1 = gzip.open("{}{}_Unknown_R1.fastq.gz".format(self.args.WorkingFolder, i), "wb")
            tmp2 = gzip.open("{}{}_Unknown_R2.fastq.gz".format(self.args.WorkingFolder, i), "wb")
            tmp1.close()
            tmp2.close()
            tmp1 = gzip.open("{}{}_GhostIndex_R1.fastq.gz".format(self.args.WorkingFolder, i), "wb")
            tmp2 = gzip.open("{}{}_GhostIndex_R2.fastq.gz".format(self.args.WorkingFolder, i), "wb")
            tmp1.close()
            tmp2.close()

        for sample in self.sample_manifest_list:
            sample_index = sample[0]
            self.sample_manifest_dictionary[sample_index] = [sample[1], sample[2], sample[3]]

            # Check for errors in the sample manifest
            if sample_index not in self.master_index_dict:
                self.log.error("Sample index {} for sample {}, replicate {} is not in Master Index File.".
                               format(sample_index, sample[1], sample[2]))
                raise SystemExit(1)

            for i in range(int(self.args.Spawn)):
                sample_key = "{}|{}".format(i, sample_index)

                # Check for more errors in the sample manifest
                if sample_key in self.demultiplex_file_names:
                    self.log.error("The index {} is duplicated.  Correct the error in {} and try again."
                                   .format(sample_index, self.args.SampleManifest))
                    raise SystemExit(1)

                # Get a list of the temporary FASTQ file names generated.
                outfile_list_dict[sample_index].append(
                    ["{}{}_{}_R1.fastq.gz".format(self.args.WorkingFolder, i, sample_index),
                     "{}{}_{}_R2.fastq.gz".format(self.args.WorkingFolder, i, sample_index)])

                # Initialize FASTQ file objects
                tmp1 = gzip.open("{}{}_{}_R1.fastq.gz".format(self.args.WorkingFolder, i, sample_index), "wb")
                tmp2 = gzip.open("{}{}_{}_R2.fastq.gz".format(self.args.WorkingFolder, i, sample_index), "wb")
                tmp1.close()
                tmp2.close()
                '''
                self.demultiplex_file_names[sample_key] = \
                    [Writer(self.log,
                            "{}{}_{}_R1.fq.gz".format(self.args.WorkingFolder, i, sample_index),
                            sample_index),
                     Writer(self.log,
                            "{}{}_{}_R2.fq.gz".format(self.args.WorkingFolder, i, sample_index),
                            sample_index)]
            '''
            # Add sample tracking framework
            self.summary_data[sample_index] = [sample[1], sample[2], 0]

        self.log.info("Dataframes built")

        return outfile_list_dict

    def file_writer(self, worker_id, batch_count, fq1_batch=None, fq2_batch=None, fq3_batch=None):
        """
        Called by multiprocessor in Odin.py.  Process FASTQ file(s) and write temp version(s).
        :param worker_id:
        :param batch_count:
        :param fq1_batch:
        :param fq2_batch:
        :param fq3_batch:
        :return:
        """
        # pool_id = os.getpid()
        # read_count = 0
        # run_start = time.time()

        temp_data_dict = collections.defaultdict(list)
        # self.log.info("Processing Batch {} using Pool {}.".format(batch_count, pool_id))
        if not fq3_batch:
            fq3_batch = fq1_batch

        for fq1_read, fq2_read, fq3_read in zip(fq1_batch, fq2_batch, fq3_batch):
            fq1_name = fq1_read[0]
            fq1_seq = fq1_read[1]
            fq1_qual = fq1_read[3]
            f1_seq_length = len(fq1_seq)

            fq2_name = fq2_read[0]
            fq2_seq = fq2_read[1]
            fq2_qual = fq2_read[3]
            f2_seq_length = len(fq2_seq)

            fq3_name = fq3_read[0]
            fq3_seq = fq3_read[1]
            fq3_qual = fq3_read[3]
            f3_seq_length = len(fq3_seq)

            # Apply Filters
            trim_5 = int(self.args.Trim5)
            trim_3 = int(self.args.Trim3)
            min_length = int(self.args.Minimum_Length) + trim_5 + trim_3

            # Filter reads based on length and number of N's.
            if (f1_seq_length < min_length or f2_seq_length < min_length
                    or fq1_seq.count("N") / f1_seq_length >= float(self.args.N_Limit)
                    or fq2_seq.count("N")/f2_seq_length >= float(self.args.N_Limit)):
                continue
            # ToDo: HaloPLEX not setup correctly yet.
            # Add the UMT's to the header.
            if self.args.HaloPLEX:
                header1 = "{0}|{0}:{1}".format(fq3_seq, fq1_name)
                header2 = "{0}|{0}:{1}".format(fq3_seq, fq2_name)

                # Fixme: This needs to be exposed to the user.
                # Short HaloPLEX reads have issues.  Found that reads <= 100 all show a 3' -1 or -2 error
                if f1_seq_length <= 100:
                    fq1_seq, fq1_qual = read_trim(fq1_seq, fq1_qual, trim5=0, trim3=3)
                if f2_seq_length <= 100:
                    fq2_seq, fq2_qual = read_trim(fq2_seq, fq2_qual, trim5=0, trim3=3)

            elif self.args.ThruPLEX:
                # header1 = "{0}|{1}".format(fastq1_read.name.split(":")[-1], fastq1_read.name)
                umt1 = fq1_seq[:6]
                umt2 = fq2_seq[:6]
                header1 = "{0}|{1}:{2}".format(umt1, umt2, fq1_name)
                header2 = "{0}|{1}:{2}".format(umt1, umt2, fq2_name)
                fq1_seq, fq1_qual = read_trim(fq1_seq, fq1_qual, trim5=len(umt1), trim3=0)
                fq2_seq, fq2_qual = read_trim(fq2_seq, fq2_qual, trim5=len(umt2), trim3=0)

            elif self.args.AlignOnly:
                header1 = fq1_name
                header2 = fq2_name

            else:
                self.log.error("Only HaloPLEX or ThruPLEX currently enabled.")
                raise SystemExit(1)

            # Trim sequences from ends if needed.
            if trim_5 > 0 or trim_3 > 0:
                fq1_seq, fq1_qual = read_trim(fq1_seq, fq1_qual, trim_5, trim_3)
                fq2_seq, fq2_qual = read_trim(fq2_seq, fq2_qual, trim_5, trim_3)

            fq1_name = header1
            fq2_name = header2

            # Sort the reads into their respective files and count.
            sample_index = self.index_search(fq1_name)

            try:
                self.sample_manifest_dictionary[sample_index]
            except IndexError:
                sample_index = "GhostIndex"

            try:
                temp_data_dict[sample_index][0] += "{}\n{}\n+\n{}".format(fq1_name, fq1_seq, fq1_qual)
                temp_data_dict[sample_index][1] += "{}\n{}\n+\n{}".format(fq2_name, fq2_seq, fq2_qual)
            except IndexError:
                temp_data_dict[sample_index] = ["{}\n{}\n+\n{}".format(fq1_name, fq1_seq, fq1_qual),
                                                "{}\n{}\n+\n{}".format(fq2_name, fq2_seq, fq2_qual)]

            # read_count += 1

        # self.summary_data['total_reads'][2] += read_count
        self.output_process(worker_id, temp_data_dict)

        return

    def output_process(self, worker_id, temp_data_dict):
        """
        Write data to compressed file.
        :param worker_id:
        :param temp_data_dict:
        """

        for sample_index in temp_data_dict:
            fq1 = gzip.open("{}{}_{}_R1.fastq.gz".format(self.args.WorkingFolder, worker_id, sample_index), "ab")
            fq2 = gzip.open("{}{}_{}_R2.fastq.gz".format(self.args.WorkingFolder, worker_id, sample_index), "ab")

            fq1.write(temp_data_dict[sample_index][0].encode())
            fq2.write(temp_data_dict[sample_index][1].encode())

            fq1.close()
            fq2.close()

    def index_search(self, fastq_name):
        def match_maker(query, unknown):
            """
            This little ditty gives us some wiggle room in identifying our indices and any other small targets.
            :param query
            :param unknown
            :return:
            """
            query_mismatch = distance(query, unknown)

            # Unknown length can be longer than target length.  Need to adjust mismatch index to reflect this.
            adjusted_query_mismatch = query_mismatch - (len(unknown) - len(query))

            return adjusted_query_mismatch

        temp_dict = collections.defaultdict(str)

        # The indices are after the last ":" in the header.
        # left_query = fastq_read.split(":")[-1].split("+")[0]
        # right_query = fastq_read.split(":")[-1].split("+")[1]

        left_query = fastq_name.split(":")[-1].split("+")[0]
        right_query = fastq_name.split(":")[-1].split("+")[1]

        for index_key in self.master_index_dict:
            left_index = self.master_index_dict[index_key][0]
            right_index = self.master_index_dict[index_key][1]

            left_match = match_maker(left_index, left_query)
            right_match = match_maker(right_index, right_query)

            if left_match == 0 == right_match:
                return index_key
            elif left_match <= 1 and right_match <= 1:
                temp_dict[left_match, right_match] = index_key

        if temp_dict:
            natsort.natsorted(temp_dict)
            return next(iter(temp_dict.values()))

        else:
            return "Unknown"

class OLDFastqProcessing:
    def __init__(self, args, log, fastq1, fastq2, fastq3=None, paired_end=False):
        self.log = log
        self.args = args
        self.fastq1_file = fastq1
        self.fastq2_file = fastq2
        self.fastq3_file = fastq3
        self.paired_end = paired_end
        self.master_index_dict = {}
        self.sample_manifest_list = None
        self.demultiplex_file_names = collections.defaultdict(list)
        self.summary_data = collections.defaultdict(list)

    def dataframe_build(self):
        """
        If we are demultiplexing the FASTQ files then build the dataframes containing the indices and file names.
        """
        self.log.info("Begin building dataframes")
        with open(self.args.MasterIndexFile) as f:
            for l in f:
                if "#" in l or not l:
                    continue
                l_list = [x for x in l.strip("\n").split("\t")]
                self.master_index_dict[l_list[0]] = [l_list[1], l_list[2]]

        self.sample_manifest_list = Tool_Box.FileParser.indices(self.log, self.args.SampleManifest)

        outfile_list_dict = collections.defaultdict(list)
        if self.args.AnalyzeUnknowns:
            outfile_list_dict["Unknown"] = \
                ["{}{}_unknown_R1.fastq".format(self.args.WorkingFolder, self.args.JobName),
                 "{}{}_unknown_R2.fastq".format(self.args.WorkingFolder, self.args.JobName)]

            '''
            I noticed with Nextera Libraries we occasionally found an index in a read that was valid but never used
            in the library prep.  I don't know, yet, how prevalent these Ghost Indices are so I am capturing them for
            later analysis.
            '''
            outfile_list_dict["GhostIndex"] = \
                ["{}{}_GhostIndex_R1.fastq".format(self.args.WorkingFolder, self.args.JobName),
                 "{}{}_GhostIndex_R2.fastq".format(self.args.WorkingFolder, self.args.JobName)]

        self.demultiplex_file_names["Unknown"] = \
            [Writer(self.log, "{}{}_unknown_R1.fastq".format(self.args.WorkingFolder, self.args.JobName), "Unknown"),
             Writer(self.log, "{}{}_unknown_R2.fastq".format(self.args.WorkingFolder, self.args.JobName), "Unknown")]

        self.demultiplex_file_names["GhostIndex"] = \
            [Writer(self.log, "{}{}_GhostIndex_R1.fastq".format(self.args.WorkingFolder, self.args.JobName), "GhostIndex"),
             Writer(self.log, "{}{}_GhostIndex_R2.fastq".format(self.args.WorkingFolder, self.args.JobName), "GhostIndex")]

        # Initialize our summary dataframe
        self.summary_data["Unknown"] = ["", "", 0]
        self.summary_data["GhostIndex"] = ["", "", 0]

        for sample in self.sample_manifest_list:
            sample_index = sample[0]

            # Check for errors in the sample manifest
            if sample_index not in self.master_index_dict:
                self.log.error("Sample index {} for sample {}, replicate {} is not in Master Index File.".
                               format(sample_index, sample[1], sample[2]))
                raise SystemExit(1)
            if sample_index in self.demultiplex_file_names:
                self.log.error("The index {} is duplicated.  Correct the error in {} and try again."
                               .format(sample[0], self.args.SampleManifest))
                raise SystemExit(1)

            # Add sample tracking framework
            self.summary_data[sample_index] = [sample[1], sample[2], 0]

            # Get a list of FASTQ file names generated.
            outfile_list_dict[sample_index] = \
                ["{}{}_{}_R1.fastq".format(self.args.WorkingFolder, self.args.JobName, sample_index),
                 "{}{}_{}_R2.fastq".format(self.args.WorkingFolder, self.args.JobName, sample_index)]

            # Initialize FASTQ file objects
            self.demultiplex_file_names[sample_index] = \
                [Writer(self.log, "{}{}_{}_R1.fastq".format(self.args.WorkingFolder, self.args.JobName, sample_index), sample_index),
                 Writer(self.log, "{}{}_{}_R2.fastq".format(self.args.WorkingFolder, self.args.JobName, sample_index), sample_index)]

        self.log.info("Dataframes built")

        return outfile_list_dict

    def file_writer(self):
        """
        Process FASTQ file(s) and write new version(s) suitable for aligners.  Return a dictionary of file names.
        :return:
        """

        outfile_list_dict = self.dataframe_build()

        self.log.info("Begin writing temporary FASTQ files.")
        current_read_count = 0
        eof = False

        # This generator returns objects not lines.
        while not eof:
            try:
                fastq1_read = next(self.fastq1_file.seq_read())
                fastq2_read = next(self.fastq2_file.seq_read())
                if self.fastq3_file is not None:
                    fastq3_read = next(self.fastq3_file.seq_read())
            except StopIteration:
                eof = True
                continue

            current_read_count += 1

            # Apply Filters
            trim_5 = int(self.args.Trim5)
            trim_3 = int(self.args.Trim3)
            min_length = int(self.args.Minimum_Length) + trim_5 + trim_3

            # Filter reads based on length and number of N's.
            if (len(fastq1_read.seq) < min_length or len(fastq2_read.seq) < min_length
                    or fastq1_read.seq.count("N") / len(fastq1_read.seq) >= float(self.args.N_Limit)
                    or fastq2_read.seq.count("N")/len(fastq2_read.seq) >= float(self.args.N_Limit)):
                continue

            # Add the UMT's to the header.
            if self.args.HaloPLEX:
                header1 = "{0}|{0}:{1}".format(fastq3_read.seq, fastq1_read.name)
                header2 = "{0}|{0}:{1}".format(fastq3_read.seq, fastq2_read.name)

                # Fixme: This needs to be exposed to the user.
                # Short HaloPLEX reads have issues.  Found that reads <= 100 all show a 3' -1 or -2 error
                if len(fastq1_read.seq) <= 100:
                    read_trim(fastq1_read, trim5=0, trim3=3)
                if len(fastq2_read.seq) <= 100:
                    read_trim(fastq2_read, trim5=0, trim3=3)

            elif self.args.ThruPLEX:
                # header1 = "{0}|{1}".format(fastq1_read.name.split(":")[-1], fastq1_read.name)
                umt1 = fastq1_read.seq[:6]
                umt2 = fastq2_read.seq[:6]
                header1 = "{0}|{1}:{2}".format(umt1, umt2, fastq1_read.name)
                header2 = "{0}|{1}:{2}".format(umt1, umt2, fastq2_read.name)
                read_trim(fastq1_read, trim5=len(umt1), trim3=0)
                read_trim(fastq2_read, trim5=len(umt2), trim3=0)

            elif self.args.AlignOnly:
                header1 = fastq1_read.name
                header2 = fastq2_read.name
            else:
                self.log.error("Only HaloPLEX or ThruPLEX currently enabled.")
                raise SystemExit(1)

            # Trim sequences from ends if needed.
            if trim_5 > 0 or trim_3 > 0:
                read_trim(fastq1_read, trim_5, trim_3)
                read_trim(fastq2_read, trim_5, trim_3)

            fastq1_read.name = header1
            fastq2_read.name = header2

            # Sort the reads into their respective files and count.
            sample_index = self.index_search(fastq1_read.name)

            self.output_process(sample_index, fastq1_read, fastq2_read)

            if current_read_count % 500000 == 0:
                self.log.info("Processed {} total reads".format(current_read_count))

        # Send signal to writer to write remaining data and close open files.
        for sample_key in self.demultiplex_file_names:
            self.log.info("Finalizing sample index {}".format(sample_key))

            try:
                fq1 = self.demultiplex_file_names[sample_key][0]
                fq2 = self.demultiplex_file_names[sample_key][1]
                fq1.fastq_write("", True)
                fq2.fastq_write("", True)
                fq1.close()
                fq2.close()
            except IndexError:
                break

        self.summary_data['total_reads'] = ["", "", current_read_count]
        self.summary_data["manifest"] = [self.sample_manifest_list, "", 0]

        self.log.info("Modified FASTQ file(s) written")

        return outfile_list_dict, self.summary_data

    def output_process(self, sample_index, r1, r2):
        """
        Send file object and sequence read objects to writer.
        :param sample_index:
        :param r1:
        :param r2:
        """
        try:
            fq1 = self.demultiplex_file_names[sample_index][0]
            fq2 = self.demultiplex_file_names[sample_index][1]
            self.summary_data[sample_index][2] += 1
        except IndexError:
            # Process and count Ghost Indices
            fq1 = self.demultiplex_file_names["GhostIndex"][0]
            fq2 = self.demultiplex_file_names["GhostIndex"][1]
            self.summary_data["GhostIndex"][2] += 1

        fq1.fastq_write(r1)
        fq2.fastq_write(r2)

    def index_search(self, fastq_read):
        def match_maker(query, unknown):
            """
            This little ditty gives us some wiggle room in identifying our indices and any other small targets.
            :param query
            :param unknown
            :return:
            """
            query_mismatch = distance(query, unknown)

            # Unknown length can be longer than target length.  Need to adjust mismatch index to reflect this.
            adjusted_query_mismatch = query_mismatch - (len(unknown) - len(query))

            return adjusted_query_mismatch

        temp_dict = collections.defaultdict(str)

        # The indices are after the last ":" in the header.
        left_query = fastq_read.split(":")[-1].split("+")[0]
        right_query = fastq_read.split(":")[-1].split("+")[1]
        for index_key in self.master_index_dict:
            left_index = self.master_index_dict[index_key][0]
            right_index = self.master_index_dict[index_key][1]

            left_match = match_maker(left_index, left_query)
            right_match = match_maker(right_index, right_query)

            if left_match == 0 == right_match:
                return index_key
            elif left_match <= 1 and right_match <= 1:
                temp_dict[left_match, right_match] = index_key

        if temp_dict:
            natsort.natsorted(temp_dict)
            return next(iter(temp_dict.values()))

        else:
            return "Unknown"


class Writer:
    """
    Write new FASTQ file.
    """
    __slots__ = ['log', 'file', 'outstring', 'counter', 'sample_index']

    def __init__(self, log, out_file_string, sample_index):
        """
        :param log:
        :param out_file_string:
        """
        #self.file = open(out_file_string, "w")
        self.file = gzip.open(out_file_string, "wb")
        self.log = log
        self.outstring = ""
        self.counter = 0
        self.sample_index = sample_index

    def fastq_write(self, read, cleanup=False):
        """
        Writes new FASTQ file by being given one read at a time.
        :param read:
        :param cleanup:
        :return:
        """

        self.counter += 1

        if not cleanup:
            try:
                assert len(read.seq) == len(read.qual)
            except AssertionError:
                self.log.error("Sequence and quality scores of different lengths! Read Name {}; Seq Length {}; Qual "
                               "Length {}".format(read.name, len(read.seq), len(read.qual)))
                raise SystemExit(1)

            self.outstring += "@{0}\n{1}\n{2}\n{3}\n".format(read.name, read.seq, read.index, read.qual)
            # self.file.write("@{0}\n{1}\n{2}\n{3}\n".format(read.name, read.seq, read.index, read.qual))
        if self.counter % 50 == 0 or cleanup:
            self.file.write(self.outstring)
            self.outstring = ""

    def fastq_list_write(self, read_list):
        """
        Writes our new FASTQ file.
        :param read_list:
        :return:
        """
        outstring = ""
        for read in read_list:
            try:
                assert len(read.seq) == len(read.qual)
            except AssertionError:
                self.log.error("Sequence and quality scores of different lengths! Read Name {0}; Seq Length {1}; Qual "
                               "Length {2}".format(read.name, len(read.seq), len(read.qual)))
                raise SystemExit(1)
            outstring += "@{0}\n{1}\n{2}\n{3}\n".format(read.name, read.seq, read.index, read.qual)

        self.file.write(outstring)
        read_list.clear()

        return True

    def close(self):
        """
        Closes FASTQ file
        :return:
        """
        self.file.close()
        return True


def read_trim(fq_seq, fq_qual, trim5, trim3):
    """
    Provide additional trimming to reads beyond the adaptor trim.
    :param fastq_read:
    :param trim5:
    :param trim3:
    """
    return fq_seq[trim5:-trim3], fq_qual[trim5:-trim3]
    # fastq_read.seq = fastq_read.seq[trim5:-trim3]
    # fastq_read.qual = fastq_read.qual[trim5:-trim3]
    #
    # if trim5 and trim3:
    #     fastq_read.seq = fastq_read.seq[trim5:-trim3]
    #     fastq_read.qual = fastq_read.qual[trim5:-trim3]
    # elif trim5:
    #     fastq_read.seq = fastq_read.seq[trim5:]
    #     fastq_read.qual = fastq_read.qual[trim5:]
    # elif trim3:
    #     fastq_read.seq = fastq_read.seq[:-trim3]
    #     fastq_read.qual = fastq_read.qual[:-trim3]
