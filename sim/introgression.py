import pandas as pd
from collections import defaultdict


def get_introgressed(ts, from_pop, to_pop):
    """Detect fragments introgressed from a population 'from_pop'
    into a set of individuals 'to_inds'."""
    source, dest = to_pop, from_pop

    source_inds = ts.samples(source)

    # dict of introgresed segments per individual
    segments = defaultdict(list)
    for migr in ts.migrations():
        if migr.source == source and migr.dest == dest:
            # find the trees corresponding to this migration
            for tree in ts.trees(leaf_lists=True):
                interval = tree.get_interval()
                if migr.left > interval[0]:  # migration within tree segment
                    continue
                if migr.right <= interval[0]:  # migration outside tree segment
                    break
                for leaf in tree.leaves(migr.node):
                    if leaf in source_inds:
                        segments[leaf].append(interval)
    segments = {ind: merge_segments(segments[ind]) for ind in source_inds}

    return segments


def merge_segments(segments):
    """Merge several adjacent tree segments into one continuous block.

    Recombination events can create several adjacent gene trees for
    a single migration event detected for an individual. This function
    merges segments belonging to those trees into a single haplotype.
    """
    segments = sorted(segments)

    merged = []
    merged_start, prev_end = segments[0]  # 'open' the first segment
    for curr_start, curr_end in segments[1:]:
        # the following segment is not directly adjacent
        if curr_start != prev_end:
            # 'close' the current segment and start a new one
            merged.append((merged_start, prev_end))
            merged_start = curr_start
        # remember the end of the current segment for the next iteration
        prev_end = curr_end
    merged.append((merged_start, prev_end))  # 'close' the last segment

    return pd.DataFrame(merged, columns=["start", "end"], dtype=float)


# segments = get_introgressed(ts, pop_params["neand"]["id"], pop_params["eur"]["id"])
# segments = pd.Series(sum(x.end - x.start)/args.seq_len for x in segments.values())
# segments.mean()
