# -*- coding: utf-8 -*-
"""
Sequence_Magic.py v0.2.0
    Jan. 12, 2017
    Dennis A. Simpson
    Added function to generate list of chromosome names.
Sequence_Magic.py v0.1.0
    Dec. 1, 2016
    Dennis A. Simpson
    Changed versioning to conform to semantic versioning (http://semver.org/).  Added author, version, package data.
Created on Fri Dec 18 12:59:12 2015

@author: Dennis A. Simpson
         University of North Carolina at Chapel Hill
         Chapel Hill, NC  27599
@copyright: 2016
"""

from Levenshtein import distance

__author__ = 'Dennis A. Simpson'
__version__ = '0.2.1'


def chromosomes(species, log, include_chrY):
    """
    Little ditty to generate a list of chromosome names.
    :param include_chrY:
    :param species:
    :param log:
    :return:
    """
    chrom_list = ["chrX", "chrM"]
    count = 0
    if species == "Mouse":
        chrom_list = ["chrX", "chrMT"]
        count = 20
    elif species == "Human":
        chrom_list = ["chrX", "chrM"]
        count = 23
    else:
        log.error("Species Must be Mouse or Human")

    for i in range(1, count):
        chrom_list.append("chr{0}".format(i))
    if include_chrY:
        chrom_list.append("chrY")

    return chrom_list


def rcomp(seq):
    """
    reverse complement our sequence
    :param seq:
    :return:
    """

    def _complement(rseq):
        """
        This is code copied from BioPython.  It is here because Python 3.3 does not give the same result as Python 3.4 when
        called from the Biopython Seq module.
        :param rseq:
        :return:
        """
        ambiguous_dna_complement = {
            "A": "T",
            "C": "G",
            "G": "C",
            "T": "A",
            "M": "K",
            "R": "Y",
            "W": "W",
            "S": "S",
            "Y": "R",
            "K": "M",
            "V": "B",
            "H": "D",
            "D": "H",
            "B": "V",
            "X": "X",
            "N": "N",
        }
        complement_mapping = ambiguous_dna_complement
        before = ''.join(complement_mapping.keys())
        after = ''.join(complement_mapping.values())
        before += before.lower()
        after += after.lower()
        ttable = str.maketrans(before, after)
        comp_seq = rseq.translate(ttable)

        return comp_seq
    
    _complement(''.join(reversed(seq)))


def match_maker(query, unknown):
    """
    This little ditty gives us some wiggle room in identifying our indices and any other small targets.
    :param query
    :param unknown
    :return:
    """

    query_mismatch = distance(query, unknown)

    # Unknown length can be longer than target length.  Need to adjust mismatch index to reflect this.
    adjusted_query_mismatch = query_mismatch-(len(unknown) - len(query))

    return adjusted_query_mismatch
