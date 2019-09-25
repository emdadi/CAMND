__author__ = 'Abolfazl'


def discard_small(partitioning, threshold):
    res = {}
    discard = set()
    for part_name, part in partitioning.iteritems():
        if len(part) >= threshold:
            res[part_name] = part
        else:
            discard = discard.union(set(part))
    return res, discard