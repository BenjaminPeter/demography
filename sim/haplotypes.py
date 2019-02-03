import numpy as np
from collections import defaultdict


def get_tracts(ts, source, dest):
    tracts = []
    for m in ts.migrations():
        if m.source == source and m.dest == dest:
            tracts.append((m.left, m.right))
    return np.array(tracts)

tracts = get_tracts(ts, pop_params["eur"]["id"], pop_params["neand"]["id"])
tracts_total = np.sum(tracts[:, 1] - tracts[:, 0])
print(tracts_total, tracts_total / args.seq_len)


i = list(enumerate(all_inds(pop_params)))

def draw(tree):
    print(tree.draw(format="unicode"))


source = pop_params["eur"]["id"]
dest = pop_params["neand"]["id"]

test_inds = ts.samples(source)

# dict of introgresed segments per individual
segments = defaultdict(list)
for migr in ts.migrations():
    if migr.source == source and migr.dest == dest:
        # find the trees corresponding to this migration
        for tree in ts.trees(leaf_lists=True):
            interval = tree.get_interval()
            if migr.left > interval[0]: # migration within tree segment
                continue
            if migr.right <= interval[0]: # migration outside tree segment
                break
            for leaf in tree.leaves(migr.node):
                if leaf in test_inds:
                    segments[leaf].append(interval)

names = dict(zip(test_inds, pop_inds("eur", pop_params)))
segments = {names[ind] : merge_segments(segments[ind]) for ind in test_inds}

def merge_segments(segments):
    segments = sorted(segments)
    merged = []
    merged_start, prev_end = segments[0]
    for curr_start, curr_end in segments[1:]:
        if curr_start != prev_end:
            merged.append((merged_start, prev_end))
            merged_start = curr_start
        prev_end = curr_end
    merged.append((merged_start, prev_end))
    return pd.DataFrame(merged, columns = ["start", "end"], dtype = float)

pd.Series(sum(x.end - x.start)/args.seq_len for x in segments.values())
_.mean()



















def combine_segs(segs, get_segs = False):
    merged = np.empty([0, 2])
    if len(segs) == 0:
        if get_segs:
            return([])
        else:
            return(0)
    sorted_segs = segs[np.argsort(segs[:, 0]), :]
    for higher in sorted_segs:
        if len(merged) == 0:
            merged = np.vstack([merged, higher])
        else:
            lower = merged[-1, :]
            if higher[0] <= lower[1]:
                upper_bound = max(lower[1], higher[1])
                merged[-1, :] = (lower[0], upper_bound)
            else:
                merged = np.vstack([merged, higher])
    if get_segs:
        return(merged)
    else:
        return(np.sum(merged[:, 1] - merged[:, 0])/seq_len)


de_seg = defaultdict(list)
for mr in ts.migrations():
    if mr.source == source and mr.dest == dest:
        for tree in ts.trees(leaf_lists=True):
            if mr.left > tree.get_interval()[0]:
                continue
            if mr.right <= tree.get_interval()[0]:
                break
            for l in tree.leaves(mr.node):
                if l in test_inds:
                    # print(l, mr)
                    de_seg[l].append(tree.get_interval())

true_de_segs = [combine_segs(np.array(de_seg[i]), True) for i in sorted(de_seg.keys())]

pd.Series(sum(x[:, 1] - x[:, 0])/args.seq_len for x in true_de_segs)
_.mean()
