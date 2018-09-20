import random
import numpy as np
from week_1 import debruijn_graph_from_kmers


def preprocess_adj_list(path):
    adj_list = {}
    with open(path, 'r') as f:
        lines = f.readlines()
        for line in lines:
            node = int(line.split("->")[0])
            neighbours = line.split("->")[1].split(",")
            adj_list[node] = list(map(int, neighbours))
    return adj_list


def write_eulerian_cycle(cycle, filename="res.txt"):
    with open(filename, 'w') as f:
        f.write('->'.join(map(str, cycle)))


def read_kmers(path):
    with open(path, 'r') as f:
        k = f.readline()
        lines = f.readlines()
    return [l.strip() for l in lines]


def get_edge(node, edges):
    """ Returns edge with node as the starting node """
    for e in edges:
        if e[0] == node:
            return e
    return None


def get_cycle(cur_node, edges):
    """ Returns cycle starting from cur_node from given edges """
    local_cycle = [cur_node]
    while get_edge(cur_node, edges):
        cur_edge = get_edge(cur_node, edges)
        cur_node = cur_edge[1]
        edges.remove(cur_edge)
        if cur_node == local_cycle[0]:
            break
        local_cycle.append(cur_node)
    return local_cycle + [local_cycle[0]]


def find_eulerian_cycle(adj_list):
    """ Find eulerian cycle from given adjacency list """

    # Get list of edges from adjecency list
    virgin_edges = []
    for node in adj_list:
        virgin_edges += [(node, nb) for nb in adj_list[node]]
    # Get random node and build initial cycle
    cur_node = random.sample(virgin_edges, 1)[0][0]
    eulerian_cycle = get_cycle(cur_node, virgin_edges)[:-1]

    # Extend current cycle with other cycles until all edges go out
    while len(virgin_edges):
        i_end = 0
        while i_end < len(eulerian_cycle):
            if get_edge(eulerian_cycle[i_end], virgin_edges):
                new_cycle = get_cycle(eulerian_cycle[i_end], virgin_edges)
                eulerian_cycle_updated = eulerian_cycle[:i_end] + new_cycle + eulerian_cycle[i_end+1:]
                eulerian_cycle = eulerian_cycle_updated[:]
                break
            i_end += 1
    return eulerian_cycle + [eulerian_cycle[0]]


def find_eulerian_path(adj_list):
    """ Find eulerian path from given adjacency list """

    def get_unbalanced_nodes(adj_list):
        edges = [(node, nb) for node in adj_list for nb in adj_list[node]]
        out_nodes = []    # nodes with more outgoing edges
        in_nodes = []   # nodes with more ingoind edges
        for node in list(set(adj_list.keys()).union(*[set(v) for v in adj_list.values()])):
            n_in = sum([e[1] == node for e in edges])
            n_out = sum([e[0] == node for e in edges])
            if n_out > n_in:
                out_nodes.append(node)
            elif n_in > n_out:
                in_nodes.append(node)
        return in_nodes, out_nodes

    # Find unbalanced nodes and connect them to form a cycle
    in_nodes, out_nodes = get_unbalanced_nodes(adj_list)
    if in_nodes[0] in adj_list.keys():
        adj_list[in_nodes[0]].append(out_nodes[0])
    else:
        adj_list[in_nodes[0]] = [out_nodes[0]]

    euler_cycle = find_eulerian_cycle(adj_list)[:-1]
    # Shift eulerian cycle for out_node to be placed in the beggining
    out_node_ind = euler_cycle.index(out_nodes[0])
    in_node_ind = euler_cycle.index(in_nodes[0])

    #euler_path = euler_cycle[out_node_ind:] + euler_cycle[:out_node_ind]
    euler_path = euler_cycle[in_node_ind+1:] + euler_cycle[:in_node_ind+1]
    out_node_ind = euler_path.index(out_nodes[0])
    in_node_ind = euler_path.index(in_nodes[0])
    return euler_path


def genome_reconstruction(kmers):
    """ Sequence genome from givel list of kmers """
    adj_list = debruijn_graph_from_kmers(kmers)
    euler_path = find_eulerian_path(adj_list)
    genome = ''.join([euler_path[0]] + [euler_path[i][-1] for i in range(1, len(euler_path))])
    return genome


def k_univ_circular_string(k):
    """ Build univercal circular string from all binary values of size k """

    def get_binary_strings(k):
        a = np.zeros((2**k, k), np.int16)
        for j in range(k):
            a[:, j] = np.tile( [0]* 2**(k-j-1) + [1]* 2**(k-j-1) , reps=2**j)
        return [''.join(map(str, list(l))) for l in a]

    kmers = get_binary_strings(k)
    adj_list = debruijn_graph_from_kmers(kmers)
    euler_path = find_eulerian_cycle(adj_list)[:-1]
    genome = ''.join([euler_path[i][-1] for i in range(len(euler_path))])
    return genome


def process_kdmers(path):
    with open(path, 'r') as f:
        first_line = f.readline()
        k = int( first_line.split(' ')[0] )
        d = int( first_line.split(' ')[1] )
        # Read first line and get kdmer length (due to the obscure newline symbol logic in the dataset)
        kdmers = []
        lines = f.readlines()
        for line in lines:
            kdmers.append([kdmer.strip() for kdmer in line.split('|')])
        return k, d, kdmers


def debruijn_graph_from_kdmers(kdmers):
    """ Construc DeBruijn graph from list of kmers """
    nodes = set()
    for i, kdmer in enumerate(kdmers):
        nodes.add((kdmer[0][1:], kdmer[1][1:]))
        nodes.add((kdmer[0][:-1], kdmer[1][:-1]))
    adj_list = dict()
    for nd in nodes:
        for kdmer in kdmers:
            if kdmer[0][:-1] == nd[0] and kdmer[1][:-1] == nd[1]:
                if nd not in adj_list.keys():
                    adj_list[nd] = [(kdmer[0][1:], kdmer[1][1:])]
                else:
                    adj_list[nd].append((kdmer[0][1:], kdmer[1][1:]))
    return adj_list


def string_spelled_by_gapped_patterns(gapped_patterns, k, d):
    """ Return genome string from kdmers. """
    adj_list = debruijn_graph_from_kdmers(gapped_patterns)
    path = find_eulerian_path(adj_list)
    kdmers_path = []
    for i in range(len(path) - 1):
        kdmers_path.append((path[i][0] + path[i+1][0][-1], path[i][1] + path[i+1][1][-1]))

    prefix_kmers = [kdmer[0] for kdmer in kdmers_path]
    suffix_kmers = [kdmer[1] for kdmer in kdmers_path]

    prefix_string = ''.join([km[0] for km in prefix_kmers[:-1]] + [prefix_kmers[-1]])
    suffix_string = ''.join([km[0] for km in suffix_kmers[:-1]] + [suffix_kmers[-1]])
    ind_overlap = k+d
    if prefix_string[ind_overlap:] == suffix_string[:-ind_overlap]:
        return prefix_string[:ind_overlap] + suffix_string


def maximal_nonbarnching_paths(adj_list):
    """ Return nonbranching paths from adjacent list """

    def get_out_deg(adj_list, nd):
        if nd not in adj_list.keys():
            return 0
        return len(adj_list[nd])

    def get_in_deg(adj_list, nd):
        return sum([nd in adj_list[src] for src in adj_list.keys()])

    def get_out_node(adj_list, nd):
        if nd not in adj_list or not adj_list[nd]:
            return None
        return adj_list[nd][0]

    def one_in_out(adj_list, nd):
        return get_out_deg(adj_list, nd) == 1 and get_in_deg(adj_list, nd) == 1

    paths = []
    isolated_nodes = []
    # (1) Collect nonbarnching pahts
    for nd in adj_list.keys():
        if not one_in_out(adj_list, nd):
            if get_out_deg(adj_list, nd) > 0:
                # For each outoging edge
                for nd_out in adj_list[nd]:
                    nonbranch_path = [nd, nd_out]
                    while one_in_out(adj_list, nd_out):
                        nd_out = get_out_node(adj_list, nd_out)
                        nonbranch_path.append(nd_out)
                    paths.append(nonbranch_path)
        elif nd not in isolated_nodes:
            # Check for isolated path
            nonbranch_path = [nd]
            nd_out = nd
            isolated = False
            while not isolated:
                isolated_nodes.append(nd_out)
                nd_out = get_out_node(adj_list, nd_out)
                if nd_out == None: break
                if nd_out in nonbranch_path:
                    isolated = True
                nonbranch_path.append(nd_out)
            if isolated:
                paths.append(nonbranch_path)
    return paths


def contig_generation(kmers):
    adj_list = debruijn_graph_from_kmers(kmers)
    print(adj_list)
    paths = maximal_nonbarnching_paths(adj_list)
    print("paths: ", len(paths))
    contigs = []
    for pth in paths:
        contigs.append(''.join([pth[0]] + [p[-1] for p in pth[1:]]))
    return contigs

if __name__ == "__main__":
    # (1) Test find_eulerian_cycle
    #adj_list = preprocess_adj_list("test_ds.txt")
    #eulerian_cycle = find_eulerian_cycle(adj_list)
    #write_eulerian_cycle(eulerian_cycle)

    # (2) Test find_eulerian_path
    #write_eulerian_cycle(find_eulerian_path(adj_list))

    # (3) Test genome_reconstruction
    #kmers = read_kmers('test_ds.txt')
    #print(genome_reconstruction(kmers))

    # (4) Test k_univ_circular_string
    #print(k_univ_circular_string(8))

    # (5) Test string_spelled_by_gapped_patterns
    #k, d, kdmers = process_kdmers('test_ds.txt')
    #print(string_spelled_by_gapped_patterns(kdmers, k, d))

    # (6) Test maximum nonbarnching paths
    #adj_list = preprocess_adj_list('test_ds.txt')
    #paths = maximal_nonbarnching_paths(adj_list)
    #    for path in paths:
    #        f.write('->'.join(map(str, path)) + '\n')

    # (6) Test contig generation problem
    # TODO: Current logic do not working with duplicate kmers
    #kmers = read_kmers('test_ds.txt')
    #paths = contig_generation(kmers)
    #with open('res.txt', 'w') as f:
    #    for p in paths:
    #        f.write(p + '\n')

    # Quiz questions

    # Problem 1
    kmers = ['AAAT', 'AATG', 'ACCC','ACGC', 'ATAC', 'ATCA', 'ATGC',
    'CAAA', 'CACC', 'CATA', 'CATC', 'CCAG', 'CCCA', 'CGCT', 'CTCA',
    'GCAT', 'GCTC', 'TACG', 'TCAC', 'TCAT', 'TGCA']
    #print(genome_reconstruction(kmers))

    # Problem 3
    k, d, kdmers = process_kdmers('test_ds.txt')
    print(string_spelled_by_gapped_patterns(kdmers, k, d))
