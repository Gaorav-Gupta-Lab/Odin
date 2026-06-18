"""

BAM_Tools v0.7.0
    May 22, 2019
    Dennis A. Simpson
    Refactor based on pylint.
BAM_Tools v0.6.0
    March 16, 2018
    Dennis A. Simpson
    Gave the user the ability to control the compression level of the sorted BAM files.

@author: Dennis A. Simpson
         RTP_Genomics LLC
         Chapel Hill, NC
@copyright: 2019
"""

import collections
from copy import deepcopy
import subprocess
import pysam
from Odin import Utilities
import Valkyries.Tool_Box as Tool_Box

__author__ = "Dennis A. Simpson"
__version__ = "0.8.0"


class AlignWriter:
    class _NullWriter:
        def write(self, family, paired_align, alignment):
            """
            This doesn't do anything so needs to be mapped and then removed.
            :param family:
            :param paired_align:
            :param alignment:
            """

        def close(self):
            """
            Again, doesn't do anything
            """
            pass

    NULL = _NullWriter()

    def __init__(self, header, bam_path, tags=None):
        if tags is None:
            self._tags = []
        else:
            self._tags = sorted(tags)

        new_header = deepcopy(header)
        if 'CO' not in new_header:
            new_header['CO'] = []

        new_header['CO'].extend([tag.header_comment for tag in self._tags])
        self._bam_path = bam_path
        self._bam_file = pysam.AlignmentFile(bam_path, "wb", header=new_header)

    @property
    def bam_file_path(self):
        """

        :return:
        """
        return self._bam_path

    def _add_bam_tags(self, family, paired_align, alignment_object):
        """

        :param family:
        :param paired_align:
        :param alignment_object:
        """
        for tag in self._tags:
            tag.set_tag(family, paired_align, alignment_object)

    def write(self, family, paired_align, alignment_object):
        """

        :param family:
        :param paired_align:
        :param alignment_object:
        """
        self._add_bam_tags(family, paired_align, alignment_object)
        self._bam_file.write(alignment_object.pysam_align_segment)

    def close(self):
        """

        """
        self._bam_file.close()


class LoggingWriter:

    class UnplacedFamily:
        def __init__(self):
            self.filter_value = 'unplaced'
            self.umi_sequence = -1

    UNPLACED_FAMILY = UnplacedFamily()

    def __init__(self, base_writer, log, chromosome):
        self._base_writer = base_writer
        self.log = log
        self._align_filter_stats = collections.defaultdict(int)
        self._family_filter_stats = collections.defaultdict(set)
        self.chromosome = chromosome

    def write(self, family, paired_align, alignment):
        """

        :param family:
        :param paired_align:
        :param alignment:
        """
        if not family:
            family = LoggingWriter.UNPLACED_FAMILY

        self._align_filter_stats[(family.filter_value, alignment.filter_value)] += 1
        self._family_filter_stats[family.filter_value].add(family.umi_sequence)

        if family == LoggingWriter.UNPLACED_FAMILY:
            self._base_writer.write(None, paired_align, alignment)
        else:
            self._base_writer.write(family, paired_align, alignment)

    @staticmethod
    def _log_line(text, count, total, filter_name):
        """

        :param text:
        :param count:
        :param total:
        :param filter_name:
        :return:
        """
        line = '{:.2f}% ({}/{}) {}: {}'
        return line.format(100 * count / total, count, total, text, filter_name)

    @property
    def _unplaced_aligns(self):
        """

        :return:
        """
        unplaced_aligns = {}
        for (fam_filter, align_filter), cnt in self._align_filter_stats.items():
            if fam_filter == LoggingWriter.UNPLACED_FAMILY.filter_value and align_filter:
                unplaced_aligns[align_filter] = cnt

        return collections.OrderedDict(Tool_Box.sort_dict(unplaced_aligns))

    @staticmethod
    def _discarded_filter_value(fam_filter, align_filter):
        """

        :param fam_filter:
        :param align_filter:
        :return:
        """
        if fam_filter == LoggingWriter.UNPLACED_FAMILY.filter_value:
            return None

        filter_values = []
        if fam_filter:
            filter_values.append(fam_filter)
        if align_filter:
            filter_values.append(align_filter)

        return "; ".join(filter_values)

    @property
    def _discarded_aligns(self):
        discarded_aligns = {}
        for (fam_filter, align_filter), cnt in self._align_filter_stats.items():
            filter_value = LoggingWriter._discarded_filter_value(fam_filter,
                                                                 align_filter)
            if filter_value:
                discarded_aligns[filter_value] = cnt

        return collections.OrderedDict(Tool_Box.sort_dict(discarded_aligns))

    @property
    def _family_stats(self):
        family_filter_stats = dict(self._family_filter_stats)
        family_filter_stats.pop(LoggingWriter.UNPLACED_FAMILY.filter_value, None)
        included_count = len(family_filter_stats.pop(None, []))
        discarded_count = 0
        filter_counts = collections.OrderedDict()

        for name, fam_ids in family_filter_stats.items():
            align_count = len(fam_ids)
            discarded_count += align_count
            filter_counts[name] = align_count

        total_count = included_count + discarded_count

        return included_count, total_count, collections.OrderedDict(Tool_Box.sort_dict(filter_counts))

    @property
    def _align_stats(self):
        """

        :return:
        """
        included_filter = (None, None)
        included_count = self._align_filter_stats[included_filter]
        excluded_count = sum([count for fam_align_filter,
                              count in self._align_filter_stats.items() if fam_align_filter != included_filter])
        total_count = included_count + excluded_count
        return included_count, excluded_count, total_count

    @staticmethod
    def _percent_stat_str(count, total):
        """

        :param count:
        :param total:
        :return:
        """
        return '{:.2f}% ({}/{})'.format(100 * count / total, count, total)

    @staticmethod
    def _log_filter_counts(filter_counts, log_method, msg_format, total, chromosome):
        for name, count in filter_counts.items():
            percent = LoggingWriter._percent_stat_str(count, total)
            log_method(msg_format.format(chromosome=chromosome, filter_name=name, percent_stat=percent))

    def _log_results(self):
        (included_align_count, excluded_align_count, tot_align_count) = self._align_stats
        (included_fam_count, total_fam_count, discarded_fam_filter_counts) = self._family_stats

        self.log.info('{}: {} alignments unplaced or discarded'
                      .format(self.chromosome, LoggingWriter._percent_stat_str(excluded_align_count, tot_align_count)))

        LoggingWriter._log_filter_counts(self._unplaced_aligns, self.log.debug,
                                         '{chromosome}: alignments unplaced: {percent_stat} {filter_name}',
                                         tot_align_count, self.chromosome)

        LoggingWriter._log_filter_counts(self._discarded_aligns, self.log.debug,
                                         '{chromosome} alignments discarded: {percent_stat} {filter_name}',
                                         tot_align_count, self.chromosome)

        LoggingWriter._log_filter_counts(discarded_fam_filter_counts, self.log.info,
                                         '{chromosome} families discarded: {percent_stat} {filter_name}',
                                         total_fam_count, chromosome=self.chromosome)

        percent_stat = LoggingWriter._percent_stat_str(included_align_count, tot_align_count)
        self.log.info('{}: {} alignments included in {} families'
                      .format(self.chromosome, percent_stat, included_fam_count))

        if included_align_count == 0:
            self.log.warning("{}: No alignments passed filters. (Was input BAM downsampled?)".format(self.chromosome))
        else:
            percent_dedup = 100 * (1 - (included_fam_count / included_align_count))
            msg = '{} {:.2f}% deduplication rate (1 - {} families/{} included alignments)'

            self.log.info(msg.format(self.chromosome, percent_dedup, included_fam_count, included_align_count))

    def close(self):
        """

        :rtype:
        """
        if self._align_filter_stats:
            self._log_results()
        self._base_writer.close()


class BamTag:

    class _NullObject:
        """Returns None for all method calls"""

        def __init__(self):
            self.included_pair_count = None
            self.filter_value = None
            self.umi_sequence = None
            self.umt = lambda *args: None
            self.is_consensus_template = lambda *args: None
            self.positions = lambda *args: None
            self.cigars = lambda *args: None

    _NULL_OBJECT = _NullObject()

    def __init__(self, tag_name, tag_type, description, get_value, analysis):
        self._tag_name = tag_name
        self._tag_type = tag_type
        self._get_value = get_value
        self._description = description

        self.header_comment = "{}\tBAM tag\t{}: {}".format(analysis, tag_name, description)

    def __lt__(self, other):
        return (self._tag_name, self._description) < (other._tag_name, other._description)

    def set_tag(self, family, paired_align, align_obj):
        """

        :param family:
        :param paired_align:
        :param align_obj:
        :return:
        """
        family = family if family else BamTag._NULL_OBJECT
        paired_align = paired_align if paired_align else BamTag._NULL_OBJECT
        value = self._get_value(family, paired_align, align_obj)
        align_obj.set_tag(self._tag_name, value, self._tag_type)


class BamFlag:
    PAIRED = 1
    PROPER_PAIR = 2
    UNMAP = 4
    MUNMAP = 8
    REVERSE = 16
    MREVERSE = 32
    READ1 = 64
    READ2 = 128
    SECONDARY = 256
    QCFAIL = 512
    DUP = 1024
    SUPPLEMENTARY = 2048


class PairedAlignment:
    """Represents the left and right align pairs from an single sequence."""
    def __init__(self, left_alignment, right_alignment):

        if left_alignment.query_name.count(":") > 7:
            left_query_name = ":".join(left_alignment.query_name.split(":")[:-4])
            right_query_name = ":".join(right_alignment.query_name.split(":")[:-4])
        else:
            left_query_name = left_alignment.query_name
            right_query_name = right_alignment.query_name

        if left_query_name != right_query_name:
            msg = 'Inconsistent query names ({} != {})'
            raise ValueError(msg.format(left_query_name, right_query_name))

        self.query_name = left_alignment.query_name
        self.left = left_alignment
        self.right = right_alignment

        left_umt = left_alignment.query_name.split(":")[0].split("|")[1]
        right_umt = right_alignment.query_name.split(":")[0].split("|")[1]

        if len(right_alignment.query_name.split(":")[0].split("|")) == 3:
            right_umt = right_alignment.query_name.split(":")[0].split("|")[2]
            left_umt = left_alignment.query_name.split(":")[0].split("|")[1]

        self.umt = (left_umt, right_umt)
        self._tag_length = len(left_umt)

    @property
    def filter_value(self):
        """

        :return:
        """
        if self.left.filter_value or self.right.filter_value:
            return self.left.filter_value, self.right.filter_value

        return None

    def cigars(self, format_string=None):
        """
        This is the CIGAR used in the tags.
        :param format_string:
        :return:
        """
        if format_string:
            return format_string.format(left=self.left.cigarstring, right=self.right.cigarstring)

        return self.left.cigarstring, self.right.cigarstring

    def positions(self, format_string=None):
        """
        Make Samtools positions.
        :param format_string:
        :return:
        """
        left_value = self.left.reference_start + 1
        right_value = self.right.reference_end + 1

        if format_string:
            return format_string.format(left=left_value, right=right_value)

        return left_value, right_value

    def __eq__(self, other):
        return self.__dict__ == other.__dict__

    def __hash__(self):
        return hash(self.left) * hash(self.right)

    def __repr__(self):
        return ("Pair({}|{}|{}, "
                "{}|{}|{})").format(self.left.query_name,
                                    self.left.reference_start,
                                    self.left.query_sequence,
                                    self.right.query_name,
                                    self.right.reference_start,
                                    self.right.query_sequence)


class OdinAlign:
    """
    This doesn't really align anything.  It is a class to manipulate certain Pysam values.  From the original Connor
    notes:
    cgates: FYI, you can use dynamic delegation via __setattr__ and __getattr__ but it's awkward and about twice as slow
    """

    def __init__(self, pysam_align_segment, filter_value=None):
        self.pysam_align_segment = pysam_align_segment
        self.filter_value = filter_value

    def __eq__(self, other):
        return other.__dict__ == self.__dict__

    # cgates: the native pysam hashing is not performant for ultradeep pileups
    def __hash__(self):
        return hash(self.filter_value) ^ \
               hash(self.pysam_align_segment.query_name) ^ \
               self.pysam_align_segment.reference_start

    @property
    def cigarstring(self):
        """

        :return:
        """
        return self.pysam_align_segment.cigarstring

    @cigarstring.setter
    def cigarstring(self, value):
        """

        :param value:
        """
        self.pysam_align_segment.cigarstring = value

    @property
    def flag(self):
        """

        :return:
        """
        return self.pysam_align_segment.flag

    @flag.setter
    def flag(self, value):
        """

        :param value:
        """
        self.pysam_align_segment.flag = value

    def get_tag(self, name, with_value_type=False):
        """

        :param name:
        :param with_value_type:
        :return:
        """
        return self.pysam_align_segment.get_tag(name, with_value_type)

    def get_tags(self, with_value_type=False):
        """

        :param with_value_type:
        :return:
        """
        return self.pysam_align_segment.get_tags(with_value_type)

    @property
    def mapping_quality(self):
        """

        :return:
        """
        return self.pysam_align_segment.mapping_quality

    @mapping_quality.setter
    def mapping_quality(self, value):
        """

        :rtype: object
        """
        self.pysam_align_segment.mapping_quality = value

    @property
    def next_reference_start(self):
        """

        :return:
        """
        return self.pysam_align_segment.next_reference_start

    @next_reference_start.setter
    def next_reference_start(self, value):
        """

        :param value:
        :return:
        """
        self.pysam_align_segment.next_reference_start = value

    @property
    def orientation(self):
        """
        If the next reference is to the left then we are on the opposite strand.  If it is to the right we are on the
        same strand and if it is neither then the next reference is possibly part of the same family
        :return:
        """
        if self.reference_start < self.next_reference_start:
            return 'left'
        if self.reference_start > self.next_reference_start:
            return 'right'

        return 'neither'

    @property
    def query_name(self):
        """

        :return:
        """
        return self.pysam_align_segment.query_name

    @query_name.setter
    def query_name(self, value):
        """

        :param value:
        """
        self.pysam_align_segment.query_name = value

    @property
    def query_sequence(self):
        """

        :return:
        """
        return self.pysam_align_segment.query_sequence

    @query_sequence.setter
    def query_sequence(self, value):
        """

        :param value:
        :return:
        """
        self.pysam_align_segment.query_sequence = value

    @property
    def query_qualities(self):
        """

        :return:
        """
        return self.pysam_align_segment.query_qualities

    @query_qualities.setter
    def query_qualities(self, value):
        """

        :param value:
        """
        self.pysam_align_segment.query_qualities = value

    @property
    def reference_end(self):
        """

        :return:
        """
        return self.pysam_align_segment.reference_end

    @property
    def reference_id(self):
        """

        :return:
        """
        return self.pysam_align_segment.reference_id

    @reference_id.setter
    def reference_id(self, value):
        """

        :param value:
        """
        self.pysam_align_segment.reference_id = value

    @property
    def reference_name(self):
        """

        :return:
        """
        return self.pysam_align_segment.reference_name

    @property
    def reference_start(self):
        """

        :return:
        """
        return self.pysam_align_segment.reference_start

    @reference_start.setter
    def reference_start(self, value):
        """

        :param value:
        """
        self.pysam_align_segment.reference_start = value

    def set_tag(self, tag_name, tag_value, value_type):
        """

        :param tag_name:
        :param tag_value:
        :param value_type:
        """
        self.pysam_align_segment.set_tag(tag_name, tag_value, value_type)

    @property
    def template_length(self):
        """

        :return:
        """
        return self.pysam_align_segment.template_length

    @template_length.setter
    def template_length(self, value):
        """

        :param value:
        """
        self.pysam_align_segment.template_length = value


def filter_alignments(pysam_alignments):
    """
    :param pysam_alignments:

    """
    filters = {'cigar unavailable': lambda a: a.cigarstring is None,
               'mapping quality < 1': lambda a: a.mapping_quality < 1,
               'not in proper pair': lambda a: a.flag & BamFlag.PROPER_PAIR == 0,
               'qc failed': lambda a: a.flag & BamFlag.QCFAIL != 0,
               'secondary alignment': lambda a: a.flag & BamFlag.SECONDARY != 0,
               'supplementary alignment': lambda a: a.flag & BamFlag.SUPPLEMENTARY != 0,
               }

    generator = Utilities.FilteredGenerator(filters)
    for pysam_align, filter_value in generator.filter(pysam_alignments):
        # If the filters pass then yield the alignment.
        if not filter_value:
            yield OdinAlign(pysam_align, filter_value)


def alignment_file(filename, mode, template=None):
    """

    :param filename:
    :param mode:
    :param template:
    :return:
    """
    return pysam.AlignmentFile(filename, mode, template)


def bamfile_sort(input_bam_filepath, threads, compression_level):
    """
    Sort and index the BAM file using pysam samtools.
    :param input_bam_filepath:
    :param threads:
    :param compression_level:
    :return:
    """
    sorted_bamfile = input_bam_filepath.replace(".bam", "_sorted.bam")
    pysam.samtools.sort("-@", str(threads), "-l", compression_level, '-o', sorted_bamfile, input_bam_filepath,
                        catch_stdout=False)

    Tool_Box.delete([input_bam_filepath])

    pysam.samtools.index(sorted_bamfile, catch_stdout=False)

    return sorted_bamfile


def merge_bam(args, log, file_list, merged_bamfile=None):
    """
    Merge BAM files that resulted from alignment of the temporary FASTQ files.  Then sort and index merged BAM.
    :param args:
    :param log:
    :param file_list:
    :param merged_bamfile:
    :return:
    """

    log.info("Merging tmp BAM files.")
    if not merged_bamfile:
        merged_bamfile = "{0}{1}_merged.bam".format(args.WorkingFolder, args.JobName)
    bam_files = " ".join(file_list)

    log.debug("Cannot figure out how to get pysam.merge() to work so I use Samtools directly.")
    cmd = "samtools merge -f -O BAM {0} {1}".format(merged_bamfile, bam_files)
    subprocess.run([cmd], shell=True)
    # pysam.samtools.merge("-f", merged_bamfile, bam_files, catch_stdout=False)

    log.info("Temporary BAM files merged.  Sort and index merged BAM file.")
    sorted_bamfile = bamfile_sort(merged_bamfile, int(args.Spawn), args.CompressionLevel)

    log.info("BAM file sorted and indexed.  Deleting Temporary BAM files..")
    Tool_Box.delete(file_list)

    return sorted_bamfile


def build_writer(log, bamfile, output_bam, tags, analysis, version):
    """

    :param log:
    :param bamfile:
    :param output_bam:
    :param tags:
    :param analysis:
    :param version:
    :return:
    """
    if not output_bam:
        log.error("No output bam file provided to BamTools.build_writer.")
        raise SystemExit(1)

    input_bam = alignment_file(bamfile, 'rb')
    hd = input_bam.header
    input_bam.close()
    header = hd.to_dict()

    header_pg_key = 'PG'
    pg_headers = header.get(header_pg_key, [])
    pg_headers.append({"ID": analysis, "PN": analysis, "VN": version})
    header[header_pg_key] = pg_headers

    return AlignWriter(header, output_bam, tags)


def total_align_count(input_bam, chromosome=None):
    """Returns count of all mapped alignments in input BAM (based on index)
    :param input_bam:
    :param chromosome:
    :return:
    """
    count = 0

    for line in pysam.samtools.idxstats(input_bam).split('\n'):
        if line:
            chrom, _, mapped, unmapped = line.strip().split('\t')

            if chrom != '*' and not chromosome:
                count += int(mapped) + int(unmapped)
            elif chrom == chromosome:
                count += int(mapped) + int(unmapped)
    return count


def build_bam_tags(analysis):
    """

    :param analysis:
    :return:
    """
    def combine_filters(fam, pair, align):
        """

        :param fam:
        :param pair:
        :param align:
        :return:
        """
        filters = [x.filter_value for x in [fam, align] if x and x.filter_value]

        if filters:
            return ";".join(filters).replace('; ', ';')

        return None

    boolean_tag_value = {True: 1}
    tags = [
        BamTag("X0", "Z", "filter (why the alignment was excluded)", combine_filters, analysis),
        BamTag("X1", "Z", "leftmost~rightmost matched pair positions",
               lambda fam, pair, align: pair.positions('{left}~{right}'), analysis),
        BamTag("X2", "Z", "L~R CIGARs",
               lambda fam, pair, align: pair.cigars('{left}~{right}'), analysis),
        BamTag("X3", "i", "unique identifier for this alignment family",
               lambda fam, pair, align: fam.umi_sequence, analysis),
        BamTag("X4", "Z", "L~R UMT barcodes for this alignment family; because of fuzzy matching the family UMT may "
                          "be distinct from the UMT of the original alignment",
               lambda fam, pair, align: fam.umt('{left}~{right}'), analysis),
        BamTag("X5", "i", "family size (number of align pairs in this family)",
               lambda fam, pair, align: fam.included_pair_count, analysis),
        BamTag("X6", "i", "presence of this tag signals that this alignment would be the template for the consensus "
                          "alignment",
               lambda fam, pair, align: boolean_tag_value.get(fam.is_consensus_template(align), None), analysis)]
    return tags
