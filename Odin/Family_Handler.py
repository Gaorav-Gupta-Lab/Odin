"""
This was taken from Connor.  Needs work.
"""

import math
import collections
import Valkyries.Tool_Box as Tool_Box


def build_family_handlers(args, consensus_writer, logger, chromosome):
    # _WriteAnnotatedAlignsHandler(annotated_writer),
    handlers = [_FamilySizeStatHandler(logger, args, chromosome),
                _MatchStatHandler(args, logger, chromosome),
                _WriteConsensusHandler(consensus_writer)]

    return handlers


class _WriteConsensusHandler:
    def __init__(self, consensus_writer):
        self._writer = consensus_writer

    def handle(self, family):
        if not family.filter_value:
            self._writer.write(family, family.consensus, family.consensus.left)
            self._writer.write(family, family.consensus, family.consensus.right)

    def end(self):
        pass


class _FamilySizeStatHandler:
    def __init__(self, logger, args, chromosome):
        self.collection = []
        self.min = None
        self.quartile_1 = None
        self.median = None
        self.mean = None
        self.quartile_3 = None
        self.max = None
        self.log = logger
        self.args = args
        self.chromosome = chromosome

    def handle(self, tag_family):
        """

        :param tag_family:
        """
        self.collection.append(len(tag_family.align_pairs))

    @staticmethod
    def _percentile(collection, percent):
        """

        :param collection:
        :param percent:
        :return:
        """
        if not collection:
            return None
        fractional_index = (len(collection)-1) * percent
        floor_index = int(math.floor(fractional_index))
        ceiling_index = int(math.ceil(fractional_index))
        index = int(fractional_index)
        if floor_index == ceiling_index:
            value = collection[index]
        else:
            lower = collection[floor_index] * (ceiling_index - fractional_index)
            upper = collection[ceiling_index] * (fractional_index - floor_index)
            value = lower + upper
        return value

    @property
    def summary(self):
        return (self.min,
                self.quartile_1,
                self.median,
                self.mean,
                self.quartile_3,
                self.max)

    def end(self):
        """
        This closes out the duplicate search.
        :return:
        """
        percentile = _FamilySizeStatHandler._percentile
        self.collection.sort()

        if not self.collection:
            Tool_Box.debug_messenger("WARNING:  {}; No Duplicates Found.".format(self.chromosome))
            self.log.warning("{}: No Duplicates Found.".format(self.chromosome))
            return

        self.min = self.collection[0]
        self.quartile_1 = percentile(self.collection, 0.25)
        self.median = percentile(self.collection, 0.50)
        self.quartile_3 = percentile(self.collection, 0.75)
        self.max = self.collection[-1]
        total_counts = sum(self.collection)
        self.mean = total_counts / len(self.collection)

        self.log.debug(('{} family_stat|family size distribution (original pair counts: min, 1Q, median, mean, 3Q, '
                        'max): {}').format(self.chromosome, ', '.join(map(lambda xz: str(round(xz, 2)), self.summary))))

        outstring = "Family Size\tCount\tFrequency"
        for i in range(self.collection[-1]+1):
            x = dict(collections.Counter(self.collection))
            try:
                outstring += "\n{}\t{}\t{}".format(i, x[i], x[i]/total_counts)
            except KeyError:
                outstring += "\n{}\t0".format(i)

        outfile = open("{}{}_{}_family_size_freq.txt"
                       .format(self.args.Working_Folder, self.args.JobName, self.chromosome), "w")
        outfile.write(outstring)
        outfile.close()


class _MatchStatHandler:
    def __init__(self, args, logger, chromosome):
        self.log = logger
        self.mismatch_allowance = int(args.UMT_Distance_Threshold)
        self.total_inexact_match_count = 0
        self.total_pair_count = 0
        self.chromosome = chromosome

    def handle(self, tag_family):
        """

        :param tag_family:
        """
        self.total_inexact_match_count += tag_family.inexact_match_count
        self.total_pair_count += len(tag_family.align_pairs)

    def end(self):
        exact_count = self.total_pair_count - self.total_inexact_match_count
        self.log.debug('{}; {:.2f}% ({}/{}) original pairs matched UMT exactly'
                       .format(self.chromosome, 100 * (1 - self.percent_inexact_match), exact_count,
                               self.total_pair_count))

        self.log.debug(('{}; {:.2f}% ({}/{}) original pairs matched by Levenshtein mismatch threshold (<={}) on left '
                        'or right UMI').format(self.chromosome, 100 * self.percent_inexact_match,
                                               self.total_inexact_match_count, self.total_pair_count,
                                               self.mismatch_allowance))

    @property
    def percent_inexact_match(self):
        """

        :return:
        """
        try:
            v = self.total_inexact_match_count / self.total_pair_count
        except ZeroDivisionError:
            v = 0
        return v
