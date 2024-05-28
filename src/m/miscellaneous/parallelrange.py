#!/usr/bin/env python
def parallelrange(rank, numprocs, globalsize):
    """
    PARALLELRANGE - from a rank, and a number of processors, figure out a range, for parallel tasks.

       Usage:
          i1, i2 = parallelrange(rank, numprocs, globalsize)
    """

    #We use floor. we under distribute rows. The rows left are then redistributed, therefore resulting in a more even distribution.
    num_local_rows = [int(globalsize / numprocs) for i in range(numprocs)]

    #There may be some rows left. Distribute evenly.
    row_rest = globalsize - numprocs * int(globalsize / numprocs)

    for i in range(row_rest):
        num_local_rows[i] = num_local_rows[i] + 1

    i1 = 0
    for i in range(rank - 1):
        i1 += num_local_rows[i]
    i2 = i1 + num_local_rows[rank - 1] - 1

    return i1, i2
