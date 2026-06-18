#!/usr/bin/env python3

"""

@author: Dennis A. Simpson
         RTP Genomics, LLC
         Durham, NC
@copyright: 2020
"""

import argparse
import collections
import fileinput
import gzip
import itertools
import pathlib
import sys
import os
import time
import subprocess
import datetime
from distutils.util import strtobool
import pathos
from Valkyries import Tool_Box, FASTQ_Tools, Alignment_Launcher, BamTools, Version_Dependencies as VersionDependencies
from Valkyries import FASTQReader
from Odin import Skadi, Mutation_Matrix, Process_Duplicates
from Odin.Mutation_Matrix import ErrorMatrix

__author__ = "Dennis A. Simpson"
__version__ = "0.10.0"
__package__ = "Odin"


def main():
    VersionDependencies.python_check()
    VersionDependencies.pysam_check()
    parser = argparse.ArgumentParser(description="A program package to find knowledge and wisdom in DNA sequences.\n{0}"
                                     .format(__version__), formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--options_file', action='store', dest='options_file', required=True,
                        help='File containing program parameters.')

    options_parser = Tool_Box.options_file(parser)
    args = options_parser.parse_args()
    start_time = time.time()
    run_start = datetime.datetime.today().strftime("%a %b %d %H:%M:%S %Y")
    log = Tool_Box.Logger(args)
    Tool_Box.log_environment_info(log, args, sys.argv)
    args, options_parser = string_to_boolean(args, options_parser)
    fastq1, fastq2 = config_error_check(args, log)

    # Define our library method.
    method = "ThruPLEX"
    if args.HaloPLEX:
        method = "HaloPLEX"

    if args.HaloPLEX or args.ThruPLEX:
        log.info("{0} v{1} {2} Analysis Beginning".format(__package__, __version__, method))
        if args.Atropos_Trim:
            fastq1, fastq2 = atropos_trim(args, log, method)

        sorted_deduplicated_bamfile_name = pangea(args, log, method, fq1=fastq1, fq2=fastq2)

        # Add some parameters to our options parser object.
        options_parser.add_argument("--BAM_File", dest="BAM_File", default=sorted_deduplicated_bamfile_name)
        options_parser.add_argument("--Variant_Search", dest="Variant_Search", default=True)
        args = options_parser.parse_args()

        skadi(args, log, options_parser)

    elif args.AlignOnly:
        method = "Other"
        if args.Atropos_Trim:
            fastq1, fastq2 = atropos_trim(args, log, method)
            
        summary_data_dict = pangea(args, log, method, fq1=fastq1, fq2=fastq2)
        align_summary_output(args, run_start, fastq1, fastq2, summary_data_dict)

    elif args.FASTQ_Quality:
        log.error("Function Not Implemented")

    elif args.Skadi:
        skadi(args, log, options_parser)

    elif args.Create_Error_Matrix:
        log.info("{0} v{1} Mutation Matrix (v{2}) Building Error Correction DataFrame"
                 .format(__package__, __version__, Mutation_Matrix.__version__))

        error_matrix = ErrorMatrix(args, log)
        error_matrix.process_data()
        error_matrix.data_output()

    elif args.Annotate_Results:
        annotate_vcf(args, log)

    else:
        log.error('No module selected to run.')

    warning = "\033[1;31m **See warnings above**\033[m" if log.warning_occurred else ''
    elapsed_time = int(time.time() - start_time)
    log.info("****Odin complete ({:,}  seconds, {:,}  Mb peak memory).****  \n{}"
             .format(elapsed_time, Tool_Box.peak_memory(), warning))


def annotate_vcf(args, log):
    """
    Right now the only annotator in the pipeline is snpEff.
    :param args:
    :param log:

    """
    # "java -Xmx12g -jar snpEff.jar  GRCh38.86 test.txt > test_ann.vcf"
    stats_file = "{}{}.html".format(args.WorkingFolder, args.JobName)
    annotated_vcf = "{}{}_ann.vcf.txt".format(args.WorkingFolder, args.JobName)

    cmd = "java -Xmx{0}g -jar {1}/snpEff.jar -c {1}/snpEff.config -stats {2} -canon {3} {4} > {5}"\
        .format(args.Java_Mem, args.snpEff_Folder, stats_file, args.Dataset, args.Input_VCF, annotated_vcf)

    log.info("Running command | {}".format(cmd))
    subprocess.run([cmd], shell=True)
    report = Skadi.ReportGenerator(args)
    report.generate_report(annotated_vcf)


def skadi(args, log, options_parser, vcf_raw_data_file=None):
    """
    Allows user to run the Skadi variant search and/or VCF file filtering submodules
    :param vcf_raw_data_file:
    :param options_parser:
    :param args:
    :param log:
    :return:
    """
    if args.Variant_Search:
        log.info("{0} v{1} Skadi Variant Search (v{2}) Beginning".format(__package__, __version__, Skadi.__version__))
        vcf_raw_data_file = Skadi.variant_search(args, log)

    if args.Filter_VCF_File:
        if vcf_raw_data_file is None:
            vcf_raw_data_file = args.Tumor_VCF
        else:
            options_parser.add_argument("--Tumor_VCF", dest="Tumor_VCF", default=vcf_raw_data_file)
            args = options_parser.parse_args()

        if not os.path.isfile(vcf_raw_data_file):
            log.error("{} not found".format(vcf_raw_data_file))
            raise SystemExit(1)

        log.info("{0} v{1} Skadi (v{2}) VCF File Filtering Beginning"
                 .format(__package__, __version__, Skadi.__version__))

        skadi_filter = Skadi.VCFfilter(args, log)
        skadi_filter.filter_vcf_file()


def align_summary_output(args, run_start, fastq1, fastq2, summary_data_dict):
    orig_fq = ""
    trim_fq = "## FASTQ1: {}\n## FASTQ2: {}\n".format(args.FASTQ1, args.FASTQ2)
    run_end = datetime.datetime.today().strftime("%a %b %d %H:%M:%S %Y")
    if args.Atropos_Trim:
        orig_fq = "## Original FASTQ1: {}\n## Original FASTQ2: {}\n".format(args.FASTQ1, args.FASTQ2)
        trim_fq = "## Trimmed FASTQ1: {}\n## Trimmed FASTQ2: {}\n".format(fastq1, fastq2)

    total_reads = summary_data_dict["total_reads"][2]
    outdata = "## Run Start: {}\n## Run End: {}\n{}{}## Total Reads: {}\n## Unknown Reads: {}\n" \
              "## Ghost Indexed Reads: {}\n\n# Index\tSample Name\tReplicate\tCount\n".\
        format(run_start, run_end, orig_fq, trim_fq, total_reads, summary_data_dict["Unknown"][2],
               summary_data_dict["GhostIndex"][2])

    sample_manifest_list = summary_data_dict["manifest"][0]

    for sample in sample_manifest_list:
        sample_index = sample[0]
        sample_name = sample[1]
        sample_replicate = sample[2]
        sample_read_count = summary_data_dict[sample_index][2]
        fraction_total = sample_read_count/total_reads
        outdata += \
            "{}\t{}\t{}\t{}\t{}\n".format(sample_index, sample_name, sample_replicate, sample_read_count, fraction_total)

    outfile = open("{}{}_Summary.txt".format(args.WorkingFolder, args.JobName), 'w')
    outfile.write(outdata)
    outfile.close()


def pangea(args, log, method, fq1=None, fq2=None):
    """
    This will run the primary Odin FASTQ processing and alignment tools.
    :param args:
    :param log:
    :param method:
    :param fq1:
    :param fq2:
    :return:
    """

    bamfile_name = getattr(args, "BAM_File", None)

    if bamfile_name is None:
        paired_end = False
        fastq1 = FASTQReader.Reader(fq1, log, int(args.BatchSize))

        if os.path.isfile(args.FASTQ2):
            log.info("Paired End Mode")
            paired_end = True
            fastq2 = FASTQReader.Reader(fq2, log, int(args.BatchSize))
        else:
            log.info("Single End Mode")
            fastq2 = FASTQReader.Reader(fq1, log, int(args.BatchSize))

        if method == "HaloPLEX":
            fastq3 = FASTQReader.Reader(args.FASTQ3, log, int(args.BatchSize))
        else:
            fastq3 = None
        '''
        # get number of lines
        with gzip.open(args.FASTQ1, 'rb') as f:
            for i, l in enumerate(f):
                pass
        batch_limit = int((0.25*i)/int(args.BatchSize))
        '''
        # log.info("Processing {} Reads in {} Batches".format(round(int(0.25*i), 0), batch_limit))

        p = pathos.multiprocessing.Pool(int(args.Spawn))
        worker = []
        counter = []
        for i in range(int(args.Spawn)):
            worker.append(i)
            counter.append(i)

        fastq_data = FASTQ_Tools.FastqProcessing(args, log, paired_end)
        outfile_list_dict = fastq_data.dataframe_build()
        batch_count = 0
        count_increment = 0
        previous_count = 0
        run_start = time.time()

        log.info("Begin Demultiplexing FASTQ File")

        while fq1:
            fq_1 = []
            fq_2 = []

            for i in worker:
                fq_1.append(next(fastq1.grouper()))
                fq_2.append(next(fastq2.grouper()))
                count_increment += 1

            p.starmap(fastq_data.file_writer, zip(worker, itertools.count(batch_count), fq_1, fq_2))

            if count_increment % (int(args.Spawn)*2) == 0:
                elapsed_time = int(time.time() - run_start)
                increment_completed = count_increment-previous_count
                log.info("{} Batches completed in {} seconds.  Total Batches Completed: {}".
                         format(increment_completed, elapsed_time, count_increment))
                previous_count = count_increment
                run_start = time.time()

            batch_count = count_increment

        p = pathos.multiprocessing.Pool(int(args.Spawn))
        log.info("Begin Combining FASTQ Files")
        p.starmap(combine_files, zip(outfile_list_dict, itertools.repeat(outfile_list_dict), itertools.repeat(args)))

        return SystemExit("Dragons")

        # log.info("Sending FASTQ file(s) to splitter.")
        # fastq_data = FASTQ_Tools.FastqProcessing(args, log, fastq1, fastq2, fastq3, paired_end)

        # Preprocess FASTQ files.
        # One file each read for use with multiple threads on BWA aligner.
        # outfile_list_dict, summary_data_dict = fastq_data.file_writer()
        '''
        Tool_Box.debug_messenger("Code hacked to run single FASTQ")
        log.warning("Code hacked to run single FASTQ")
        outfile_list_dict = collections.defaultdict(list)
        outfile_list_dict["sgMre11_T1"] = \
            ["{}".format(args.FASTQ1), "{}".format(args.FASTQ2)]
        '''
        log.info("Send FASTQ data to aligner.")
        aligner = Alignment_Launcher.AlignmentLauncher(args, log, paired_end)

        # Running sequential alignment.  Not the best but hey.
        # Aligner using multiple threads option
        aligner.run_aligner(outfile_list_dict)
        
        log.info("Sort and Index BAM file.")
        for sample_index in outfile_list_dict:
            bam_file = "{}{}_{}.bam".format(args.WorkingFolder, args.JobName, sample_index)

            # Sort and Index BAM file
            sorted_bamfile_name = BamTools.bamfile_sort(bam_file, int(args.Spawn), args.CompressionLevel)

        if args.Demultiplex:
            log.info("Sending FASTQ Files to pigz")
            for f1, f2 in outfile_list_dict.values():
                Tool_Box.compress_files(args, log, f1)
                Tool_Box.compress_files(args, log, f2)

        if args.AlignOnly:
            return summary_data_dict

    else:
        # Old method in here to take an existing sorted BAM file and deduplicate.
        if not os.path.isfile("{}".format(args.BAM_File)):
            log.error("BAM File {} not found".format(args.BAM_File))
            raise SystemExit(1)

        sorted_bamfile_name = "{}".format(args.BAM_File)

    log.info("Processing Duplexes in BAM file with {} protocol".format(method))
    sorted_deduped_bamfile_name = Process_Duplicates.bam_processing(args, log, sorted_bamfile_name, method)

    return sorted_deduped_bamfile_name


def combine_files(sample_index, outfile_list_dict, args):
    r1_list = []
    r2_list = []

    for file_list in outfile_list_dict[sample_index]:
        r1_list.append(str(file_list[0]))
        r2_list.append(str(file_list[1]))

    with gzip.open("{}{}_{}_R1.fq.gz".format(args.WorkingFolder, args.JobName, sample_index), "wb") as outfile:
        hook = fileinput.hook_compressed
        for line in fileinput.input(files=r1_list, openhook=hook):
            outfile.write(line)
    with gzip.open("{}{}_{}_R2.fq.gz".format(args.WorkingFolder, args.JobName, sample_index), 'wb') as outfile:
        hook = fileinput.hook_compressed
        for line in fileinput.input(files=r2_list, openhook=hook):
            outfile.write(line)

    r1_list.extend(r2_list)
    Tool_Box.delete(r1_list)


def atropos_trim(args, log, method):
    """
    Trim adapters from reads with Atropos.  This creates the Atropos config file and runs Atropos.
    :return:
    :param args:
    :param log:
    :param method:
    :return:
    """

    fastq1_trimmed = "{}{}_Trim.R1.fastq.gz".format(args.WorkingFolder, args.JobName)
    trim_report = "--report-file {}Atropos_{}_Trim_Report.txt".format(args.WorkingFolder, args.JobName)

    if method == "HaloPLEX":
        fastq2_trimmed = "{}{}_Trim.R3.fastq.gz".format(args.WorkingFolder, args.JobName)
    else:
        fastq2_trimmed = "{}{}_Trim.R2.fastq.gz".format(args.WorkingFolder, args.JobName)

    if args.NextSeq_Trim:
        nextseq_trim = "--nextseq-trim 1\n"
        op_order = "--op-order GAWCQ\n"
    else:
        nextseq_trim = "\n"
        op_order = "--op-order AWCQ\n"

    module = "trim\n"
    aligner = "--aligner {}\n".format(args.Atropos_Aligner)
    threads = "--threads {}\n".format(args.Threads)
    additional_adapters = ""
    r1_adapters_5p = "-G file:{}\n".format(args.Anchored_Adapters_5p)
    r2_adapters_5p = "-g file:{}\n".format(args.Anchored_Adapters_5p)
    r1_adapters_3p = "-A file:{}\n".format(args.Anchored_Adapters_3p)
    r2_adapters_3p = "-a file:{}\n".format(args.Anchored_Adapters_3p)
    fq1_trimmed = "-o {}\n".format(fastq1_trimmed)
    fq2_trimmed = "-p {}\n".format(fastq2_trimmed)
    fastq1 = "-pe1 {}\n".format(args.FASTQ1)
    fastq2 = "-pe2 {}\n".format(args.FASTQ2)
    error_rate = "--error-rate {}\n".format(args.Adapter_Mismatch_Fraction)
    stats_method = "--stats both \n"
    batch_size = "--batch-size {}\n".format(args.BatchSize)
    read_queue = ""
    results_queue = ""
    if args.Read_Queue_Size:
        read_queue = "--read-queue-size {}\n".format(args.Read_Queue_Size)
    if args.Result_Queue_Size:
        results_queue = "--result-queue-size {}\n".format(args.Result_Queue_Size)

    config_block = "{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}"\
            .format(module, aligner, threads, additional_adapters, r1_adapters_5p, r2_adapters_5p, r1_adapters_3p,
                    r2_adapters_3p, fq1_trimmed, fq2_trimmed, fastq1, fastq2, nextseq_trim, op_order, error_rate,
                    stats_method, batch_size, read_queue, results_queue, trim_report)

    config_file = open("{}{}_Atropos_Config.txt".format(args.WorkingFolder, args.JobName), "w")
    config_file.write(config_block)
    config_file.close()

    log.info("Beginning Atropos Trim of {} library".format(method))
    subprocess.run("atropos --config {}{}_Atropos_Config.txt".format(args.WorkingFolder, args.JobName), shell=True)

    return fastq1_trimmed, fastq2_trimmed


def string_to_boolean(args, options_parser):
    """

    :param args:
    :param options_parser:
    :return:
    """
    options_parser.set_defaults(HaloPLEX=bool(strtobool(args.HaloPLEX)))
    options_parser.set_defaults(ThruPLEX=bool(strtobool(args.ThruPLEX)))
    options_parser.set_defaults(FASTQ_Quality=bool(strtobool(args.FASTQ_Quality)))
    options_parser.set_defaults(Skadi=bool(strtobool(args.Skadi)))
    options_parser.set_defaults(AlignOnly=bool(strtobool(args.AlignOnly)))
    options_parser.set_defaults(Create_Error_Matrix=bool(strtobool(args.Create_Error_Matrix)))
    options_parser.set_defaults(Annotate_Results=bool(strtobool(args.Annotate_Results)))
    options_parser.set_defaults(Atropos_Trim=bool(strtobool(args.Atropos_Trim)))
    if bool(strtobool(args.HaloPLEX)) or bool(strtobool(args.ThruPLEX)) or bool(strtobool(args.Skadi)):
        if not bool(strtobool(args.Skadi)):

            options_parser.set_defaults(NextSeq_Trim=bool(strtobool(args.NextSeq_Trim)))
        elif bool(strtobool(args.Skadi)):
            options_parser.set_defaults(Inclde_All=bool(strtobool(args.Include_All)))
            options_parser.set_defaults(Variant_Search=bool(strtobool(args.Variant_Search)))

            if bool(strtobool(getattr(args, "Strict_Boundaries", "False"))):
                options_parser.set_defaults(Strict_Boundaries=bool(strtobool(args.Strict_Boundaries)))
            else:
                options_parser.set_defaults(Strict_Boundaries=False)

        options_parser.set_defaults(Demultiplex=bool(strtobool(getattr(args, "Demultiplex", "False"))))
        options_parser.set_defaults(Filter_VCF_File=bool(strtobool(args.Filter_VCF_File)))
        options_parser.set_defaults(AnalyzeUnknowns=bool(strtobool(args.AnalyzeUnknowns)))

    args = options_parser.parse_args()

    return args, options_parser


def config_error_check(args, log):
    """
    Get FASTQ files and check options file for common errors.
    :param args:
    :param log:
    :return:
    """
    fastq1 = getattr(args, "FASTQ1", None)
    fastq2 = getattr(args, "FASTQ2", None)
    fastq3 = getattr(args, "FASTQ3", None)

    # Handle potential configuration file errors.
    if not fastq1:
        log.error("FASTQ file parameter missing from options file. Correct error and try again.")
        raise SystemExit(1)

    if not pathlib.Path(fastq1).is_file():
        log.error("FASTQ file {} not found.  Correct error and run again.".format(fastq1))
        raise SystemExit(1)

    if fastq2 and not pathlib.Path(fastq1).is_file():
        log.error("FASTQ file {} not found.  Correct error and run again.".format(fastq2))
        raise SystemExit(1)

    if fastq3 and not pathlib.Path(fastq3).is_file():
        log.error("FASTQ file {} not found.  Correct error and run again.".format(fastq3))
        raise SystemExit(1)

    if fastq1 and fastq1 == fastq2:
        log.error("FASTQ1 file name is the same as FASTQ2")
        raise SystemExit(1)
    if fastq1 and fastq1 == fastq3:
        log.error("FASTQ1 file name is the same as FASTQ3")
        raise SystemExit(1)
    if fastq3 and fastq2 == fastq3:
        log.error("FASTQ2 file name is the same as FASTQ3")
        raise SystemExit(1)

    if args.HaloPLEX and args.ThruPLEX:
        log.error("More than one library type selected")
        raise SystemExit(1)

    return fastq1, fastq2


if __name__ == '__main__':
    main()
