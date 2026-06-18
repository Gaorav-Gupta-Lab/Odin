"""
Cython version of the column cruncher

@author: Dennis A. Simpson
         Asystbio Laboratories, LLC
         Durham, NC
@copyright: 2019
"""

# cython: language_level=3

__author__ = 'Dennis A. Simpson'
__version__ = "0.3.1"
__package__ = 'Odin'

cpdef object column_cruncher(int min_base_qual, object pileup_column, str ctg, int position, object refseq, object del_tracking):
    """
    Process pileup column.
    :param del_tracking: 
    :param min_base_qual:
    :param pileup_column:
    :param ctg:
    :param position
    :param refseq:
    
    """

    cdef int quality
    cdef str name, strand, qbase, position_filter, query_alignment_start, alignment_end, ins_string, ref_base, del_string
    cdef object read
    cdef set duplicate_tracking = set()
    duplicate_tracking_add = duplicate_tracking.add

    cdef list depth = []
    cdef list snv_strand = []
    cdef list snv_qbase = []
    cdef list del_strand = []
    cdef list del_qbase = []
    cdef list ins_strand = []
    cdef list ins_qbase = []

    ref_base = refseq.fetch(reference=ctg, start=position, end=position + 1)
    for read, quality, name in zip(pileup_column.pileups,
                                   pileup_column.get_query_qualities(),
                                   pileup_column.get_query_names()):
            umt = name.split(":")[0]
            position_filter = str(read.query_position)
            query_alignment_start = str(read.alignment.query_alignment_start)
            alignment_end = str(read.alignment.query_alignment_end)

            if position_filter == query_alignment_start or position_filter == alignment_end:
                continue

            if umt in duplicate_tracking or (read.indel == 0 and quality < min_base_qual):
                continue

            if read.indel != 0 and quality < 13:
                continue

            duplicate_tracking_add(umt)
            # Get strand orientation.
            strand = "forward"
            if read.alignment.is_reverse:
                strand = "reverse"

            depth.append(strand)

            if read.is_del or read.is_refskip:
                del_string = ""
                del_tracking[position].append(umt)
                if umt in del_tracking[position-1]:
                    del_string = "*"


                elif read.indel < 0:
                    for i in range(-1*read.indel):
                        del_string += refseq.fetch(reference=ctg, start=position+i, end=position+i+1)
                else:
                    del_string = ref_base

                del_strand.append(strand)
                del_qbase.append(del_string)

            elif read.indel > 0:
                ins_string = ""
                for i in range(read.indel):
                    ins_string += read.alignment.query_sequence[read.query_position+i]
                ins_qbase.append(ins_string)
                ins_strand.append(strand)

            elif not read.alignment.query_sequence[read.query_position] == ref_base:
                snv_strand.append(strand)
                snv_qbase.append(read.alignment.query_sequence[read.query_position])

    return [[snv_strand, snv_qbase], [del_strand, del_qbase], [ins_strand, ins_qbase], depth], del_tracking