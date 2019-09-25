from operator import itemgetter
# matplotlib.use("PDF")
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as cluster


figure_number = 1


def get_cut_mods(cut):
    cut_mods = {}
    for c, mod in enumerate(cut):
        if mod not in cut_mods:
            cut_mods[mod] = []
        cut_mods[mod].append(c)
    return cut_mods


def get_mod_labels_in_cut(first_cut_mod, second_cut):
    return map(lambda fcm: second_cut[fcm], first_cut_mod)


def get_super_in_compressed_tree(sub, sub_super_map, super_sub_map):
    super = sub_super_map[sub]
    if len(super_sub_map[super]) > 1 or super not in sub_super_map:
        return super
    else:
        return get_super_in_compressed_tree(super, sub_super_map, super_sub_map)


def collapse_more_than_two_way_branches(super_sub_map):
    keys_snapshot = list(super_sub_map.keys())
    for super in keys_snapshot:
        subs = super_sub_map[super]

        if len(subs) >= 3:
            cur_sub = subs[0]
            for i in range(1, len(subs) - 1):
                nex_sub = subs[i]

                new_clust = (super[0], super[1] - 0.5 + 0.0001 * i)  # NOTE: this causes the final sort of branch points for building merge occurrences work properly
                super_sub_map[new_clust] = [cur_sub, nex_sub]

                cur_sub = new_clust

            super_sub_map[super] = [cur_sub, subs[-1]]

    print 'collapse to two-way branches done!'
    return super_sub_map


def build_merge_occurrences(super_sub_map, cluster_id_map, len_map, observation_count):
    super_sub_keys_sorted_by_depth = sorted(super_sub_map.keys(), cmp=lambda x, y: cmp(x[1], y[1]) if cmp(x[0], y[0]) == 0 else cmp(y[0], x[0]), reverse=False)  # first sort by depth
    for s in super_sub_map.keys():
        super_sub_map[s] = sorted(super_sub_map[s], key=itemgetter(1))  # sort children of each module by number (depth of all are the same)

    merge_occurrences = {}
    tree_height = max([s[0] for s in super_sub_keys_sorted_by_depth])

    for superi, super in enumerate(super_sub_keys_sorted_by_depth):
        subs = super_sub_map[super]

        new_cluster_id = (observation_count + 1) + superi
        cluster_id_map[super] = new_cluster_id
        len_map[super] = len_map[subs[0]] + len_map[subs[1]]

        merge_occurrences[new_cluster_id] = (cluster_id_map[subs[0]], cluster_id_map[subs[1]], tree_height - super[0] + 1, len_map[super])

    print 'build merge occurrences done!'
    return merge_occurrences


def save_dendogram(tree, tree_height, thresholds, file_path):
    global figure_number

    # draw dendogram
    plt.figure(figure_number, figsize=(50, 40), dpi=400)
    plt.title("dendogram")

    cluster.dendrogram(tree, show_leaf_counts=False, no_labels=True)
    # plt.show()
    # fig = plt.gcf()
    # fig.set_size_inches(30, 12)

    # draw threshold lines
    axes = plt.gca()
    xmin, xmax = axes.get_xlim()

    thresh_handles = []
    for thr in sorted(thresholds, reverse=True):
        line_height = thr * tree_height
        thresh_handle, = plt.plot([xmin, xmax], [line_height, line_height], label='%.2f%%' % (100 * thr), linewidth=2.0, linestyle='--')
        thresh_handles.append(thresh_handle)

    plt.legend(handles=thresh_handles, loc=2, prop={'size': 30})

    # save figure
    file_format = 'png'
    plt.savefig('%s.%s' % (file_path, file_format), format=file_format)

    figure_number += 1

