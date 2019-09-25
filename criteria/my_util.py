import random
import os
import errno
import time

from scipy import io as spio


def allPairtShortestPath(adjacency):
    v = len(adjacency)
    inf = 9999999
    dist = [[inf if x == 0 else x for x in line] for line in adjacency]    #-1 is infinity
    for i in xrange(v): dist[i][i] = 0
    for k in xrange(v):
        for i in xrange(v):
            for j in xrange(v):
                if dist[i][k] + dist[k][j] < dist[i][j]:
                    dist[i][j] = dist[i][k] + dist[k][j]
    return dist


def getAdjacencyForStoichPresent(S, revs):
    rc = len(S[0])
    mc = len(S)
    exts = mc
    reacts = [{} for i in xrange(rc)]
    for r in xrange(rc):
        prod = False
        cons = False
        for m in xrange(mc):
            f = S[m][r]
            if f:
                reacts[r][m] = f
                if f > 0:
                    prod = True
                if f < 0:
                    cons = True
        if not prod:
            reacts[r][exts] = 1
            exts += 1
        if not cons:
            reacts[r][exts] = -1
            exts += 1

            #        if len(reacts[r]) == 2: #only for non-hyper arcs (for hyper arc another method is taken, follow on to see ...)
            #            if revs[r] == 1:    #reversible! add reverse to reacts list!
            #                reacts.append({})
            #                for x in reacts[r]:
            #                    reacts[-1][x] = -reacts[r][x]

    adj = [[0 for j in xrange(exts)] for i in xrange(exts)]

    dummsLabels = []
    dumms = exts
    for r, re in enumerate(reacts):
        cs = dict((m, re[m]) for m in re if re[m] < 0) #{m:re[m] for m in re if re[m]<0}
        ps = dict((m, re[m]) for m in re if re[m] > 0) #{m:re[m] for m in re if re[m]>0}
        lb = str(r) #if r < rc else ''
        if len(re) == 2:
            c = cs.keys()[0]
            p = ps.keys()[0]
            adj[c][p] = (lb, cs[c], ps[p], revs[r] == 1)   #label, source factor, dest factor!
        elif len(re) > 2:
            for x in adj:   # add a new column and a new row to adj
                x.append(0)
            adj.append([0 for i in xrange(dumms + 1)])

            for x in cs:    # connect from consumers to a new dummy node for hyperarc
                adj[x][dumms] = ('', cs[x], 1, revs[r] == 1)
                #                if revs[r] == 1:
            #                    adj[dumms][x] = ('',-1,1)

            for x in ps:    # connect from the new dummy node to products of hyperarc
                adj[dumms][x] = ('', -1, ps[x], False)

            dummsLabels.append(str(r))
            dumms += 1
        else:
            print 'BAD CONVERT FROM s TO adjacency! A REACTIONS WITH LESS THAN TWO SIDES!'

    return adj, mc, exts, dumms, dummsLabels


def getAdjacencyFromStoichiometric(S, undirected=True, completeNotBipartite=True, addExternals=False):
    mc = len(S)
    rc = len(S[0])
    adj = [[0] * mc for i in xrange(mc)]
    #    for i in xrange(mc):
    #        for j in xrange(mc):
    #            for c in xrange(rc):
    #                if S[i][c] != 0 and S[j][c] != 0 and i != j:
    #                    adj[i][j] = 1

    exts = mc
    for c in xrange(rc):
        ins = []
        outs = []
        for r in xrange(mc):
            if S[r][c] < 0:
                ins.append(r)
            elif S[r][c] > 0:
                outs.append(r)

        if addExternals:
            if ins or outs:
                add = False
                if not ins:
                    ins.append(exts)
                    add = True
                elif not outs:
                    outs.append(exts)
                    add = True
                if add:
                    exts += 1
                    for r in adj:
                        r.append(0)
                    adj.append([0] * exts)

        for i in ins:
            for j in outs:
                adj[i][j] = 1
                if undirected:
                    adj[j][i] = 1

        if completeNotBipartite:
            for gees in [ins, outs]:
                for i in gees:
                    for j in gees:
                        adj[i][j] = 1

    return adj


def removeReaction(reactName, S, rxns, revs):
    idx = rxns.index(reactName)
    rxns.pop(idx)
    revs.pop(idx)
    for m in S:
        m.pop(idx)


def removeMetabolite(metabName, S, mets, metNames, compartNames):
    bmIdx = mets.index(metabName)
    mets.pop(bmIdx)
    metNames.pop(bmIdx)
    S.pop(bmIdx)
    if compartNames is not None:
        compartNames.pop(bmIdx)


def adjacencyMatrixToList(adjm):
    return dict([(i, [[j for j, x in enumerate(xs) if x], False, None]) for i, xs in enumerate(adjm)])


def getTotalConnectedComponentsCount(adjacency, ans, groupOfMets):
    """
    :param graph: {  nodeId: [[x, y, z], False, partLabelHolder], ...  }  -> means node 'nodeId' is connected to x, y, and z.
        partLabelHolder is important, it will be filled from 'ans' so that each node partition is determined and used when counting
        components count!
    :return:
    """
    parts = groupMetabolitesByLabel(ans, groupOfMets)
    for lb, part in parts.iteritems():
        for m in part:
            adjacency[m][2] = lb
    return getConnectedComponentCount(adjacency)


def getConnectedComponentCount(graph):
    """
    Only computes count!
    :param graph: {  nodeId: [[x, y, z], False, partLabel], ...  }  -> means node 'nodeId' is connected to x, y, and z
    :return:
    """
    for v in graph.values():
        v[1] = False

    componentsCount = 0
    for nodeID in graph:
        if not graph[nodeID][1]:
            componentsCount += 1
            recursivelyMark(nodeID, graph)
    return componentsCount


def recursivelyMark(nodeID, nodes):
    (connections, visited, nodeLb) = nodes[nodeID]
    if visited:
        return

    connections = filter(lambda n2id: nodeLb == nodes[n2id][2], connections)

    nodes[nodeID][1] = True
    for connectedNodeID in connections:
        recursivelyMark(connectedNodeID, nodes)


def getConnectedComponents(graph):
    """
    :param graph: { nodeId: [[x, y, z], False, ignored], ...  }  -> means node 'nodeId' is connected to x, y, and z
    :return: a GENERATOR which generates connected components (CC). Each CC, in turn, is a generator of nodeIds of the CC
    """
    seen = set()

    def component(node):
        nodes = set([node])
        while nodes:
            node = nodes.pop()
            seen.add(node)
            nodes |= set(graph[node][0]) - seen
            yield node

    for node in graph:
        if node not in seen:
            yield component(node)


def sampleUST(graph):
    """
    :param graph: { nodeId: [[x, y, z], False, ignored], ...  }  -> means node 'nodeId' is connected to x, y, and z
    :return: a spanning tree over the graph (with the same format)
    """
    for v in graph.values():
        v[1] = False

    tree = {}

    ccs = getConnectedComponents(graph)
    for cc in ccs:
        c = dict([(n, graph[n]) for n in cc])
        tree.update(sampleUST4CC(c))

    return tree


def sampleUST4CC(graph):
    tree = {}
    sz = len(graph)

    todo = set(graph.keys())
    firstPath = True
    while todo:
        x = todo.pop()
        if graph[x][1]:
            continue

        walk = [x]
        walknodes = set(walk)

        while True:
            yi = int(len(graph[x][0]) * random.random())
            y = graph[x][0][yi]

            if y in walknodes:
                p = walk.pop()
                while y != p:
                    walknodes.remove(p)
                    p = walk.pop()

            walk += [y]
            walknodes.add(y)
            if graph[y][1]:
                break
            else:
                x = y

            if len(walk) == sz:
                break

            if firstPath:
                # this is checked on every elongation, so when len(walk)=sz/2 probability of n successfull enlongation
                # without breaking is (3/4)^n which shrinks rapidly as n grows! 0.5 is here not to let the probability
                # be (2/4)^n which is much smaller :)
                if random.random() < 0.5 * len(walk) / sz:
                    break

        firstPath = False

        for w in walk:
            graph[w][1] = True

        if len(walk) == 1:
            tree[walk[0]] = [[], False, None]
        else:
            for i in xrange(len(walk) - 1):
                if walk[i] not in tree:
                    tree[walk[i]] = [[], False, None]
                tree[walk[i]][0] += [walk[i + 1]]
                if walk[i + 1] not in tree:
                    tree[walk[i + 1]] = [[], False, None]   #TODO:TC: any change in graph structure makes problem here!
                tree[walk[i + 1]][0] += [walk[i]]

    return tree


def good_split(l, c):
    parts = []

    partst = 0
    parted = 0
    indc = False
    for z in l:
        if z == '"':
            if indc:
                indc = False
            else:
                indc = True
        if z == c and not indc:
            p = l[partst: parted]
            if p.startswith('"') and p.endswith('"'):
                p = p[1:-1]
            parts.append(p)
            partst = parted + 1

        parted += 1

    p = l[partst: parted]
    if p.startswith('"') and p.endswith('"'):
        p = p[1:-1]
    parts.append(p)

    return parts


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def graph_by_explode_reactions_to_complete_bipartite(stoich, mets):
    grph = []

    for i in range(len(stoich)):
        for j in range(i + 1, len(stoich)):
            for z in range(len(stoich[0])):
                if stoich[i][z] != 0 and stoich[j][z] != 0:
                    grph.append((mets[i], mets[j]))

    return list(set(grph))


def graph_by_explode_metabolites_to_complete_bipartite(stoich, rxns):
    grph = []

    for i in range(len(stoich[0])):
        for j in range(i + 1, len(stoich[0])):
            for z in range(len(stoich)):
                if stoich[z][i] != 0 and stoich[z][j] != 0:
                    grph.append((rxns[i], rxns[j]))

    return list(set(grph))


def write_line(f, l):
    f.write(l)
    f.write('\n')


def prepare_matlab_file_and_exec_and_wait_finish(script_file, result_file, no_prepare=False):
    if not no_prepare:
        f = open(script_file + '.m', 'r')
        lines = []
        for l in f:
            lines.append(l.strip())
        f.close()

        f = open(script_file + '.m', 'w')

        write_line(f, 'try')

        for l in lines:
            write_line(f, '    %s' % l)

        write_line(f, 'catch')
        write_line(f, '    zzz = lasterror;')
        write_line(f, '    save %s zzz' % result_file)
        write_line(f, 'end')
        write_line(f, 'quit')

        f.close()

    try:
        os.remove(result_file)
    except:
        pass

    # res = os.system('LD_PRELOAD="/usr/lib64/libstdc++.so.6" matlab -nodisplay -nosplash -nodesktop -r "%s"' % script_file)
    res = os.system('matlab -nodisplay -nosplash -nodesktop -r "%s"' % script_file)

    while not os.path.isfile(result_file):
        time.sleep(5)
    while file_is_locked_for_writing(result_file):
        time.sleep(5)


def file_is_locked_for_writing(f):
    try:
        mf = open(f, 'a+b')
    except IOError:
        return True

    mf.close()
    return False


def try_load_matlab_result(file):
    f = spio.loadmat(file)
    if 'zzz' in f:
        raise Exception(f['zzz'][0].tolist()[0])
    return f


def get_remove_reaction_matlab_script(model_name, react):
    ret = []

    ret.append("")
    ret.append("idx = find(not(cellfun('isempty', strfind(%s.rxns, '%s'))));" % (model_name, react))
    ret.append("%s.rxns = %s.rxns([1:idx-1, idx+1:end]);" % (model_name, model_name))
    ret.append("%s.rev = %s.rev([1:idx-1, idx+1:end]);" % (model_name, model_name))
    ret.append("%s.lb = %s.lb([1:idx-1, idx+1:end]);" % (model_name, model_name))
    ret.append("%s.ub = %s.ub([1:idx-1, idx+1:end]);" % (model_name, model_name))
    ret.append("%s.c = %s.c([1:idx-1, idx+1:end]);" % (model_name, model_name))
    ret.append("%s.rxnNames = %s.rxnNames([1:idx-1, idx+1:end]);" % (model_name, model_name))
    ret.append("%s.S = [%s.S(:, 1:idx-1), %s.S(:, idx+1:end)];" % (model_name, model_name, model_name))
    ret.append("")

    return ret


def read_standard_written_2d_list(str, st_idx, ed_idx, convert):
    res = []

    idx = str.find('[', st_idx + 1)
    while idx != -1:
        parts = str[idx + 1: str.index(']', idx)].split(',')
        res.append(map(convert, parts))
        idx = str.find('[', idx + 1)

    return res


def get_organism(species):
    return species[: species.index('_')]


def small_to_large_sort(alist, anorder):
    def sort_order_compare(x, y):
        try:
            xi = anorder.index(x)
        except:
            xi = len(alist) + len(anorder) + 1

        try:
            yi = anorder.index(y)
        except:
            yi = len(alist) + len(anorder) + 1

        return cmp(xi, yi)

    return sorted(alist, cmp=sort_order_compare)
