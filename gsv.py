import pysam, argparse

# TODO: check chromosome length to make sure query bounds don't go too far
def get_query_bounds (svtype, position, conf_int, fetch_flank):
    """
    Returns the boudns to query the BAM file around the given position/confidence interval for the specified svtype and fetch_flank.
    >>> get_query_bounds('DEL', 10, (-4, 10), 1)
    (5, 22)
    >>> get_query_bounds('DEL', 100, (-10, 0), 20)
    (70, 121)
    """
    if svtype == 'DEL':
        return (max(position + conf_int[0] - fetch_flank, 0), position + conf_int[1] + fetch_flank + 1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('test', action='store_true', help='run test cases')
    args = parser.parse_args()
    if args.test:
        import doctest
        doctest.testmod()
