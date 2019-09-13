import networkx as nx
import math
import csv
import random as rand
import sys

outf = None

def print_or_write(s):
    global outf
    if outf is not None:
        outf.write(str(s))
        outf.write('\n')
        outf.flush()
    else:
        print s

#this method just reads the graph structure from the file
def buildG(G, file_, delimiter_):
    #construct the weighted version of the contact graph from cgraph.dat file
    #reader = csv.reader(open("/home/kazem/Data/UCI/karate.txt"), delimiter=" ")
    reader = csv.reader(open(file_), delimiter=delimiter_)
    for line in reader:
        if float(line[2]) != 0.0:
            G.add_edge(int(line[0]),int(line[1]),weight=float(line[2]))

#keep removing edges from Graph until one of the connected components of Graph splits into two
#compute the edge betweenness
def CmtyGirvanNewmanStep(G):
    #print "call CmtyGirvanNewmanStep"
    init_ncomp = nx.number_connected_components(G)    #no of components
    ncomp = init_ncomp
    while ncomp <= init_ncomp:
        bw = nx.edge_betweenness_centrality(G, weight='weight')    #edge betweenness for G
        #find the edge with max centrality
        max_ = max(bw.values())
        #find the edge with the highest centrality and remove all of them if there is more than one!
        for k, v in bw.iteritems():
            if float(v) == max_:
                G.remove_edge(k[0],k[1])    #remove the central edge
        ncomp = nx.number_connected_components(G)    #recalculate the no of components

#compute the modularity of current split
def _GirvanNewmanGetModularity(G, deg_, m_):
    New_A = nx.adj_matrix(G)
    New_deg = {}
    New_deg = UpdateDeg(New_A, G.nodes())
    #Let's compute the Q
    comps = nx.connected_components(G)    #list of components
    #print 'no of comp: %d' % len(comps)
    Mod = 0    #Modularity of a given partitionning
    len_comps = 0
    for c in comps:
        len_comps += 1
        EWC = 0    #no of edges within a community
        RE = 0    #no of random edges
        for u in c:
            EWC += New_deg[u]
            RE += deg_[u]        #count the probability of a random edge
        Mod += ( float(EWC) - float(RE*RE)/float(2*m_) )
    print_or_write('no of comp: %d' % len_comps)
    Mod = Mod/float(2*m_)
    #print "Modularity: %f" % Mod
    return Mod

def UpdateDeg(A, nodes):
    deg_dict = {}
    n = len(nodes)
    B = A.sum(axis = 1)
    for i in range(n):
        deg_dict[nodes[i]] = B[i, 0]
    return deg_dict

#run GirvanNewman algorithm and find the best community split by maximizing modularity measure
def runGirvanNewman(G, Orig_deg, m_):
    #let's find the best split of the graph
    BestQ = 0.0
    Q = 0.0
    while True:
        CmtyGirvanNewmanStep(G)
        Q = _GirvanNewmanGetModularity(G, Orig_deg, m_)
        print_or_write("current modularity: %f" % Q)
        if Q > BestQ:
            BestQ = Q
            Bestcomps = nx.connected_components(G)    #Best Split
            print_or_write("START_COMP")
            for best_comp in Bestcomps:
                print_or_write(best_comp)
            print_or_write("END_COMP")
        if G.number_of_edges() == 0:
            break
    if BestQ > 0.0:
        print_or_write("Best Q: %f" % BestQ)
        for best_comp in Bestcomps:
            print_or_write(best_comp)
    else:
        print_or_write("Best Q: %f" % BestQ)

def main(argv):
    global outf

    if len(argv) < 2:
        sys.stderr.write("Usage: %s <input graph> [<output file>]\n" % (argv[0],))
        return 1

    if len(argv) == 3:
        outf = open(argv[2], 'w')

    graph_fn = argv[1]
    G = nx.Graph()  #let's create the graph first
    buildG(G, graph_fn, ',')

    print_or_write(G.nodes())
    print_or_write(G.number_of_nodes())

    n = G.number_of_nodes()    #|V|
    A = nx.adj_matrix(G)    #adjacenct matrix

    m_ = 0.0    #the weighted version for number of edges
    for i in range(0,n):
        for j in range(0,n):
            m_ += A[i,j]
    m_ = m_/2.0
    print_or_write("m: %f" % m_)

    #calculate the weighted degree for each node
    Orig_deg = {}
    Orig_deg = UpdateDeg(A, G.nodes())

    #run Newman alg
    runGirvanNewman(G, Orig_deg, m_)

    if len(argv) == 3:
        outf.close()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
