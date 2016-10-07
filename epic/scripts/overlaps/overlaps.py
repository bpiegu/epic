from io import StringIO
from subprocess import check_output
from os.path import basename
import logging

import pandas as pd
from joblib import Parallel, delayed


def overlap_matrix_nucleotides(all_files, window_size, nb_cpu):

    nucleotide_matrixes = Parallel(n_jobs=nb_cpu)(
        delayed(_create_overlap_matrix_nucleotides)(bed_file, all_files, window_size)
        for bed_file in all_files)

    nucleotides_df = pd.concat(nucleotide_matrixes).reset_index(drop=True)

    return nucleotides_df


def _compute_nucleotide_overlap(df, window_size):

    main_file = df.Main.ix[0,0]

    nucleotide_overlap = df.drop_duplicates("Chromosome Start End File".split()).groupby("Chromosome Start End".split()).size().value_counts().to_frame().reset_index()

    nucleotide_overlap.columns = "Other Overlaps".split()
    nucleotide_overlap.Overlaps = nucleotide_overlap.Overlaps * window_size

    nucleotide_overlap.insert(0, "Main", main_file)

    return nucleotide_overlap


def _create_overlap_matrix_nucleotides(bed_file, all_files, window_size):

    all_files_str = " ".join(all_files)

    base_bed = basename(bed_file).split(".")[0]
    logging.info("Processing {} at nucleotide level".format(base_bed))
    command = "bedtools makewindows -w {} -b {} | bedtools intersect -wo -filenames -a - -b {}".format(window_size, bed_file, all_files_str)
    output = check_output(command, shell=True).decode()

    df = pd.read_table(StringIO(output), usecols=[0, 1, 2, 3], header=None, names="Chromosome Start End File".split())

    df.File = df.File.apply(basename).str.split(".", expand=True).ix[:, 0]

    df.insert(0, "Main", base_bed)

    return _compute_nucleotide_overlap(df, window_size)


def overlap_matrix_regions(all_files, nb_cpu):

    regions_matrixes = Parallel(n_jobs=nb_cpu)(
        delayed(_create_overlap_matrix_regions)(bed_file, all_files)
        for bed_file in all_files)

    regions_df = pd.concat(regions_matrixes).reset_index(drop=True)

    return regions_df



def _create_overlap_matrix_regions(bed_file, all_files):

    all_files_str = " ".join(all_files)

    base_bed = basename(bed_file).split(".")[0]
    logging.info("Processing {}".format(base_bed))
    command = "bedtools intersect -wo -a {} -b {} -filenames".format(bed_file, all_files_str)
    output = check_output(command, shell=True).decode()

    df = pd.read_table(StringIO(output), usecols=[0, 1, 2, 3], header=None, names="Chromosome Start End File".split())
    df.File = df.File.apply(basename).str.split(".", expand=True).ix[:, 0]

    df.insert(0, "Main", base_bed)

    return _compute_region_overlap(df)



def _compute_region_overlap(df):

    main_file = df.Main.ix[0,0]

    nb_overlap = df.drop_duplicates("Chromosome Start End File".split()).groupby("Chromosome Start End".split()).size().value_counts().to_frame().reset_index()

    # nb_overlap = df.drop_duplicates("Chromosome Start End File".split()).File.value_counts().to_frame().reset_index()
    nb_overlap.columns = "Other Overlaps".split()

    nb_overlap.insert(0, "Main", main_file)

    return nb_overlap
