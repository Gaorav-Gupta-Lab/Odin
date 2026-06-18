"""
v0.1.8
    August, 23, 2019
    Dennis Simpson
    Made a few tweaks so that the debug mode will only process 300,000 reads

Deduplicates BAM file based on UMI's or UMT's and map position.

Based on the program "Connor".  Connor is copyright 2014 Bioinformatics Core, University of Michigan
Licensed under the Apache License, Version 2.0 (the "License");

"""
import collections
from collections import defaultdict, Counter
from copy import deepcopy
from functools import partial
import sys
import traceback
import itertools
import pathos
from sortedcontainers import SortedSet
import Odin.Family_Handler as FamilyHandler
import Valkyries.Tool_Box as Tool_Box
from Valkyries import BamTools
import Valkyries.Sequence_Magic as Sequence_Magic

__version__ = "0.1.8"


class _TagFamily:
    umi_sequence = 0

    def __init__(self, umt, alignments, inexact_match_count, consensus_threshold, family_filter=lambda x: None):

        # Count how many tag families are found.
        self.umi_sequence = _TagFamily.umi_sequence
        _TagFamily.umi_sequence += 1

        self._umt = umt
        (self.distinct_cigar_count, majority_cigar) = _TagFamily._get_dominant_cigar_stats(alignments)

        self.align_pairs = alignments
        self._mark_minority_cigar(majority_cigar)
        self.inexact_match_count = inexact_match_count
        self.consensus_threshold = consensus_threshold
        self.consensus = self._build_consensus(umt, self.align_pairs)
        self.included_pair_count = sum([1 for p in self.align_pairs if not p.filter_value])
        self.filter_value = family_filter(self)

    # Formats the UMT that ends up in the BAM file tag.
    def umt(self, format_string=None):
        """

        :param format_string:
        :return:
        """
        if format_string:
            return format_string.format(left=self._umt[0], right=self._umt[1])

        return self._umt

    @staticmethod
    def _get_cigarstring_tuple(paired_alignment):
        return paired_alignment.left.cigarstring, paired_alignment.right.cigarstring

    def _mark_minority_cigar(self, majority_cigar):
        for pair in self.align_pairs:

            if _TagFamily._get_cigarstring_tuple(pair) != majority_cigar:
                pair.left.filter_value = "minority CIGAR"
                pair.right.filter_value = "minority CIGAR"

    def _generate_consensus_sequence(self, alignment_pairs):

        left_alignments = []
        right_alignments = []
        for align_pair in alignment_pairs:
            left_alignments.append(align_pair.left)
            right_alignments.append(align_pair.right)

        left_consensus_seq = self._consensus_sequence(left_alignments)
        right_consensus_seq = self._consensus_sequence(right_alignments)

        return left_consensus_seq, right_consensus_seq

    def _consensus_sequence(self, alignments):
        """
        The original method was placing families of more than 1 alignment into the one alignment consensus.
        :param alignments:
        :return:
        """
        # seq_counter = Counter([a.query_sequence for a in alignments])
        query_list = []
        for a in alignments:
            query_list.append(a.query_sequence)
        seq_counter = Counter(query_list)

        consensus_seq = None
        if len(query_list) == 1:
            consensus_seq = query_list[0]

        elif len(seq_counter) == 2:
            majority_seq, majority_count = seq_counter.most_common(1)[0]
            total = len(alignments)

            if majority_count / total > self.consensus_threshold:
                consensus_seq = majority_seq

        if consensus_seq is None:
            consensus = []
            for i in range(len(alignments[0].query_sequence)):
                counter = Counter([s.query_sequence[i:i + 1] for s in alignments])
                base = counter.most_common(1)[0][0]
                freq = counter[base] / sum(counter.values())

                if freq >= self.consensus_threshold:
                    consensus.append(base)
                else:
                    consensus.append("N")

            consensus_seq = "".join(consensus)

        return consensus_seq

    @staticmethod
    def _select_template_alignment_pair(alignment_pairs):
        top_alignment_pair = None
        best_template = (0, None)
        for alignment_pair in alignment_pairs:
            query_name = alignment_pair.left.query_name
            qual_sum = alignment_pair.left.mapping_quality + alignment_pair.right.mapping_quality

            if (-qual_sum, query_name) < best_template:
                best_template = (-qual_sum, query_name)
                top_alignment_pair = alignment_pair

        return top_alignment_pair

    @staticmethod
    def _get_dominant_cigar_stats(alignments):
        counter = Counter([_TagFamily._get_cigarstring_tuple(s) for s in alignments])
        number_distict_cigars = len(counter)
        top_two_cigar_count = counter.most_common(2)
        dominant_cigar = top_two_cigar_count[0][0]
        dominant_cigar_count = top_two_cigar_count[0][1]

        if number_distict_cigars > 1 and dominant_cigar_count == top_two_cigar_count[1][1]:
            dominant_cigar = sorted(counter.most_common(), key=lambda x: (-x[1], x[0]))[0][0]

        return number_distict_cigars, dominant_cigar

    def is_consensus_template(self, alignment_object):
        """

        :param alignment_object:
        :return:
        """
        return self.consensus.left.query_name == alignment_object.query_name

    def _build_consensus(self, umt, align_pairs):
        included_pairs = [p for p in align_pairs if not p.filter_value]
        template_pair = _TagFamily._select_template_alignment_pair(included_pairs)

        left_align = deepcopy(template_pair.left)
        right_align = deepcopy(template_pair.right)

        (left_sequence, right_sequence) = self._generate_consensus_sequence(included_pairs)

        left_align.query_sequence = left_sequence
        right_align.query_sequence = right_sequence
        left_align.query_qualities = template_pair.left.query_qualities
        right_align.query_qualities = template_pair.right.query_qualities
        consensus_pair = BamTools.PairedAlignment(left_align, right_align)

        if not (umt[0] or umt[1]) or (len(umt[1]) != len(umt[0])):
            msg = "Each UMT must match tag_length ({})"
            raise ValueError(msg.format(len(umt[0])))

        return consensus_pair


def _build_coordinate_pairs(alignment_object):
    """
    This generator yields mates in a family.
    :param alignment_object:
    """
    missing_mate_filter = 'read mate was missing or excluded'
    coords = defaultdict(dict)

    # Block added by Dennis Simpson to deal with preprocessed FASTQ headers.
    for alignment in alignment_object:
        key = (alignment.reference_id, alignment.next_reference_start)
        query_name = alignment.query_name

        if alignment.query_name.count(":") > 7:
            # this is here to handle headers that have not been preprocessed to move the index to the front.
            query_name = ":" .join(alignment.query_name.split(":")[:-4])

        if alignment.orientation == 'left':
            coords[key][query_name] = alignment

        elif alignment.orientation == 'neither':

            if key in coords and query_name in coords[key]:
                align1 = coords[key].pop(query_name)

                yield BamTools.PairedAlignment(align1, alignment)

            else:
                coords[key][query_name] = alignment

        else:
            key = (alignment.reference_id, alignment.reference_start)
            coord = coords[key]
            l_align = coord.pop(query_name, None)

            # Clear empty coordinate dict[key]
            if not coord:
                del coords[key]

            if l_align:
                yield BamTools.PairedAlignment(l_align, alignment)

            else:
                alignment.filter_value = missing_mate_filter

    for aligns in coords.values():
        for align in aligns.values():
            align.filter_value = missing_mate_filter


class CoordinateFamilyHolder:
    """
    Encapsulates how stream of paired aligns are iteratively released as
    sets of pairs which share the same coordinate (coordinate families)
    """
    def __init__(self):
        self._coordinate_family = defaultdict(partial(defaultdict, list))
        self._right_coords_in_progress = defaultdict(SortedSet)
        self.pending_pair_count = 0
        self.pending_pair_peak_count = 0

    def _add(self, pair):
        def _start(align):
            return align.reference_name, align.reference_start

        self._right_coords_in_progress[pair.right.reference_name].add(pair.right.reference_start)
        right_coord = self._coordinate_family[_start(pair.right)]
        right_coord[_start(pair.left)].append(pair)
        self.pending_pair_count += 1
        self.pending_pair_peak_count = max(self.pending_pair_count, self.pending_pair_peak_count)

    def _completed_families(self, reference_name, rightmost_boundary):
        """returns one or more families whose end < rightmost boundary"""
        in_progress = self._right_coords_in_progress[reference_name]
        while in_progress:
            right_coord = in_progress[0]
            if right_coord < rightmost_boundary:
                in_progress.pop(0)
                left_families = self._coordinate_family.pop((reference_name, right_coord), {})

                for family in sorted(left_families.values(), key=lambda x: x[0].left.reference_start):
                    family.sort(key=lambda x: x.query_name)
                    self.pending_pair_count -= len(family)

                    yield family
            else:
                break

    def _remaining_families(self):
        for left_families in self._coordinate_family.values():
            for family in left_families.values():
                self.pending_pair_count -= len(family)
                yield family

            left_families.clear()
        self._coordinate_family.clear()

    def build_coordinate_families(self, paired_aligns):
        """
        Given a stream of paired aligns, return a list of pairs that share
        same coordinates (coordinate family).  Flushes families in progress
        when any of:
        a) incoming right start > family end
        b) incoming chrom != current chrom
        c) incoming align stream is exhausted
        :param paired_aligns:
        """
        rightmost_start = None
        current_chrom = None

        for pair in paired_aligns:

            if rightmost_start is None:
                rightmost_start = pair.right.reference_start
                current_chrom = pair.right.reference_name

            # New Chromosome
            if current_chrom != pair.right.reference_name:
                # Since I now do this on a single chromosome at a time in parallel, this should never be called.
                self._right_coords_in_progress[current_chrom].clear()
                rightmost_start = None
                current_chrom = None

                for family in self._remaining_families():
                    yield family

            # New Alignment Start.
            elif pair.right.reference_start != rightmost_start:
                # A new alignment start means we are done with the family.
                right = pair.right
                for family in self._completed_families(right.reference_name, right.reference_start):
                    yield family

            self._add(pair)

        for family in self._remaining_families():
            yield family


def _build_tag_families(tagged_paired_aligns, ranked_tags, mismatch_threshold, consensus_threshold,
                        family_filter=lambda x: None):

    """
    Partition paired aligns into families.
    Each read is considered against each ranked tag until all reads are
    partitioned into families.
    :param tagged_paired_aligns:
    :param ranked_tags:
    :param mismatch_threshold:
    :param consensus_threshold:
    :param family_filter:
    :return:
    """
    tag_aligns = defaultdict(set)
    tag_inexact_match_count = defaultdict(int)

    for paired_align in tagged_paired_aligns:
        (left_umi, right_umi) = paired_align.umt

        for best_tag in ranked_tags:
            if paired_align.umt == best_tag:
                tag_aligns[best_tag].add(paired_align)
                break

            # It is possible that the UMI's only differ because of sequencing errors.
            elif Sequence_Magic.match_maker(left_umi, best_tag[0]) <= mismatch_threshold \
                    or Sequence_Magic.match_maker(right_umi, best_tag[1]) <= mismatch_threshold:

                tag_aligns[best_tag].add(paired_align)
                tag_inexact_match_count[best_tag] += 1
                break

    tag_families = []
    for tag in sorted(tag_aligns):
        tag_family = _TagFamily(tag, tag_aligns[tag], tag_inexact_match_count[tag], consensus_threshold, family_filter)
        tag_families.append(tag_family)

    return tag_families


def _rank_tags(tagged_paired_aligns):
    """
    Return the list of tags ranked from most to least popular.
    :param tagged_paired_aligns:
    :return:
    """
    tag_count_dict = defaultdict(int)
    for paired_align in tagged_paired_aligns:
        tag_count_dict[paired_align.umt] += 1

    tags_by_count = Tool_Box.sort_dict(tag_count_dict)
    ranked_tags = [tag_count[0] for tag_count in tags_by_count]

    return ranked_tags


def _progress_logger(args, base_generator, total_rows, log, chromosome, coordinate_holder):
    """

    :param args:
    :param base_generator:
    :param total_rows:
    :param log:
    :param chromosome:
    :param coordinate_holder:
    """
    row_count = 0
    breakpoint_counter = 0
    completed = 0

    for item in base_generator:
        if breakpoint_counter == row_count:
            log.info("{}; {}% ({}/{}) alignments processed".format(chromosome, completed, row_count, total_rows))
            log.debug("{}; {} Mb peak memory".format(chromosome, Tool_Box.peak_memory()))
            log.debug("{}; {} pending alignment pairs; {} peak pairs"
                      .format(chromosome, coordinate_holder.pending_pair_count,
                              coordinate_holder.pending_pair_peak_count))

            breakpoint_counter += int(total_rows/10)
            completed = round(100*(row_count/total_rows), 2)
        row_count += 1

        if args.Verbose == "DEBUG":
            # This is some debug code that allows faster testing of the deduplication.
            limit = 300000
            if row_count == limit:
                Tool_Box.debug_messenger("Limiting output here for debugging")
                log.warning("Odin read limit , {}; {} alignments; {} Mb peak memory"
                            .format(chromosome, row_count, Tool_Box.peak_memory()))
                break
        yield item


def _build_family_filter(args):
    """

    :param args:
    :return:
    """
    min_family_size = int(args.Minimum_Family_Size)
    too_small_msg = 'family too small (<{})'.format(min_family_size)

    def family_size_filter(family):
        """

        :param family:
        :return:
        """
        if family.included_pair_count < min_family_size:
            return too_small_msg
        return None

    return family_size_filter


def deduplicate_alignments(args, input_bamfile, consensus_writer, log, chromosome):
    """

    :param args:
    :param input_bamfile:
    :param consensus_writer:
    :param log:
    :param chromosome:
    """
    log.info('reading input bam [{}] for {}', input_bamfile, chromosome)
    total_aligns = BamTools.total_align_count(input_bamfile, chromosome)
    family_filter = _build_family_filter(args)
    handlers = FamilyHandler.build_family_handlers(args, consensus_writer, log, chromosome)

    bamfile = BamTools.alignment_file(input_bamfile, 'rb')

    coord_family_holder = CoordinateFamilyHolder()
    pysam_alignments = _progress_logger(args, bamfile.fetch(chromosome), total_aligns, log, chromosome, coord_family_holder)
    filtered_aligns_gen = BamTools.filter_alignments(pysam_alignments)
    paired_align_gen = _build_coordinate_pairs(filtered_aligns_gen)
    coord_family_gen = coord_family_holder.build_coordinate_families(paired_align_gen)

    for coord_family in coord_family_gen:

        ranked_tags = _rank_tags(coord_family)
        tag_families = _build_tag_families(coord_family, ranked_tags, int(args.UMT_Distance_Threshold),
                                           float(args.Consensus_Freq_Threshold), family_filter)

        for handler in handlers:
            for tag_family in tag_families:
                handler.handle(tag_family)

    for handler in handlers:
        handler.end()

    bamfile.close()


def entry_point(argvs, chromosome):
    """
    This is where the parallel jobs enter.  Jobs run for each chromosome.
    :param argvs:
    :param chromosome:
    """

    args, log, input_bamfile, analysis = argvs
    log.info("Working on Alignments on {}".format(chromosome))

    bam_tags = BamTools.build_bam_tags(analysis)
    dedup_bam_outfile = "{0}{1}_{2}_dedup.bam".format(args.Folder, args.JobName, chromosome)
    consensus_writer = BamTools.build_writer(log, input_bamfile, dedup_bam_outfile, bam_tags, analysis, __version__)

    deduplicate_alignments(args, input_bamfile, consensus_writer, log, chromosome)

    consensus_writer.close()
    log.info("Completed {}".format(chromosome))


def bam_processing(args, log, input_bamfile, analysis):
    """
    Deduplicate BAM file.  Done in parallel, one chromosome per CPU/thread.
    :param args:
    :param log:
    :param input_bamfile:
    :param analysis:
    :return:
    """
    try:
        log.info('{} Deduplication Beginning (v{})'.format(analysis, __version__))
        chromosome_list = Sequence_Magic.chromosomes(args.Species, log, True)

        argvs = args, log, input_bamfile, analysis
        p = pathos.multiprocessing.Pool(int(args.Spawn))
        p.starmap(entry_point, zip(itertools.repeat(argvs), chromosome_list))
        dedup_list = []

        log.debug("Combining Frequency Data")
        frequency_dict = collections.defaultdict(int)
        total_count = 0

        for chromosome in chromosome_list:
            dedup_list.append("{0}{1}_{2}_dedup.bam".format(args.WorkingFolder, args.JobName, chromosome))
            infile = open("{}{}_{}_family_size_freq.txt".format(args.WorkingFolder, args.JobName, chromosome))
            first_line = True
            for line in infile:
                if first_line:
                    first_line = False
                    continue

                family_size = int(line.split("\t")[0])
                family_count = int(line.split("\t")[1])
                total_count += family_count
                frequency_dict[family_size] += family_count

            infile.close()

        outstring = "Family Size\tFamily Count\tFrequency"
        for family_size in sorted(frequency_dict):
            family_count = frequency_dict[family_size]
            outstring += "\n{}\t{}\t{}".format(family_size, family_count, family_count/total_count)

        outfile = open("{}{}_family_size_freq.txt".format(args.WorkingFolder, args.JobName), "w")
        outfile.write(outstring)
        outfile.close()
        deduped_bamfile_name = "{0}{1}_dedup.bam".format(args.WorkingFolder, args.JobName)
        deduped_sorted_bamfile = BamTools.merge_bam(args, log, dedup_list, merged_bamfile=deduped_bamfile_name)

        return deduped_sorted_bamfile

    except Tool_Box.UsageError as usage_error:
        message = "{} usage problem: {}".format(analysis, usage_error)
        log.error(message)
        raise SystemExit(1)

    except Exception:
        if log:
            show = log.error
        else:
            show = partial(sys.stderr)
        show("Well, that did not work so well.  Hopefully the traceback will tell you why.")
        show(traceback.format_exc())
        raise SystemExit(1)
