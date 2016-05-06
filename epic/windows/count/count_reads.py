from os.path import expanduser, abspath


def count_reads(bed_file,
                genome,
                fragment_size,
                window_size,
                keep_duplicates,
                paired_end,
                nb_cpus=1):

    # # making paths absolute, hence equal, to trigger cache on every
    # # opportunity
    # bed_file = abspath(expanduser(bed_file))

    # from epic.windows.count.count_reads_in_windows import count_reads_in_windows
    # bed_dfs = count_reads_in_windows(bed_file, genome, fragment_size,
    #                                  window_size, keep_duplicates, paired_end,
    #                                  nb_cpus)
    # # else:
    # #     from epic.windows.count.pandas_count_reads_in_windows import pandas_count_reads_in_windows
    # #     bed_dfs = pandas_count_reads_in_windows(bed_file, genome,
    # #                                             fragment_size, window_size,
    # #                                             keep_duplicates, nb_cpus)

    # return bed_dfs
