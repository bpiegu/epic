import logging
from sys import platform
from re import search, IGNORECASE
from io import BytesIO
from subprocess import check_output
from argparse import Namespace

import pysam

import pandas as pd

from epic.config import logging_settings

__author__ = "Endre Bakken Stovner https://github.com/endrebak/"
__license__ = "MIT"


def find_readlength(args):
    # type: (Namespace) -> int
    """Estimate length of reads based on 10000 first."""

    bed_file = args.treatment[0]
    nseq = 10000
    readlengths = None

    if args.bam:
        bam = pysam.AlignmentFile(bed_file, "rb")
        iter = bam.fetch()
        count_seq = 0
        bam_reads_length = []
        for aln in iter:
           bam_reads_length.append(aln.query_length)
           readlengths = pd.Series(bam_reads_length)
           count_seq +=1
           if count_seq == nseq: break
    else:
        filereader = "cat "
        if bed_file.endswith(".gz") and search("linux", platform, IGNORECASE):
            filereader = "zcat "
        elif bed_file.endswith(".gz") and search("darwin", platform, IGNORECASE):
            filereader = "gzcat "
        elif bed_file.endswith(".bz2"):
            filereader = "bzgrep "

        command = filereader + ("{} | head -{}").format(bed_file, nseq)
        output = check_output(command, shell=True)

        df = pd.read_table(
            BytesIO(output),
            header=None,
            usecols=[1, 2],
            sep="\t",
            names=["Start", "End"])

        readlengths = df.End - df.Start

    mean_readlength = readlengths.mean()
    median_readlength = readlengths.median()
    max_readlength = readlengths.max()
    min_readlength = readlengths.min()

    logging.info((
        "Used first {} reads of {} to estimate a median read length of {}\n"
        "Mean readlength: {}, max readlength: {}, min readlength: {}.").format(
            nseq, bed_file, median_readlength, mean_readlength, max_readlength,
            min_readlength))

    return median_readlength


def get_closest_readlength(estimated_readlength):
    # type: (int) -> int
    """Find the predefined readlength closest to the estimated readlength.

    In the case of a tie, choose the shortest readlength."""

    readlengths = [36, 50, 75, 100]
    differences = [abs(r - estimated_readlength) for r in readlengths]
    min_difference = min(differences)
    index_of_min_difference = [i
                               for i, d in enumerate(differences)
                               if d == min_difference][0]

    return readlengths[index_of_min_difference]
