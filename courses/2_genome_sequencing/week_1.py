def read_from_file(file_path):
    with open(file_path, 'r') as f:
        return f.read().split()

def write_to_file(res_list, file_name='res.txt'):
    with open(file_name, 'w') as f:
        f.write('\n'.join(res_list))


def decompose_string(str, k):
    """ Decompose string str to the k-mers """
    kmers = set()
    for i in range(len(str) - k + 1):
        kmers.add(str[i:i+k])
    return sorted(kmers)


def string_spelled(genome_kmer_path):
    """ Construct genome string from consequtive kmers """
    str = genome_kmer_path[0]
    for kmer in genome_kmer_path[1:]:
        str += kmer[-1]
    return str


def graph_overlap(kmers):
    """ Build adjecency list from list of kmers """

    def get_postfixes(kmer, kmers_list):
        postfixes = []
        for km in kmers_list:
            if kmer[1:] == km[:-1]:
                postfixes.append(km)
        return postfixes

    adj_list = dict()
    for i, kmer in enumerate(kmers):
        postfixes = get_postfixes(kmer, kmers[:i] + kmers[i+1:])
        if len(postfixes):
            adj_list[kmer] = postfixes
    return adj_list


def de_bruijn_graph(str, k):
    """ Return ajdacency list of k-1 mers as nodes """

    nodes = list(set([str[i:i+(k-1)] for i in range(len(str)-k+2)]))
    #adj_list = {k: [] for k in nodes}
    adj_list = {}
    for node in nodes:
        edges = [str[i:i+k] for i in range(len(str)+k-1) if str[i:i+k][:-1] == node]
        if len(edges):
            if node not in adj_list:
                adj_list[node] = []
            for e in edges:
                adj_list[node] += [e[1:]]
    return adj_list


def debruijn_graph_from_kmers(kmers):
    """ Construc DeBruijn graph from list of kmers """
    nodes = set()
    for kmer in kmers:
        nodes.add(kmer[1:])
        nodes.add(kmer[:-1])
    adj_list = dict()
    for nd in nodes:
        for kmer in kmers:
            if kmer[:-1] == nd:
                if nd not in adj_list.keys():
                    adj_list[nd] = [kmer[1:]]
                else:
                    adj_list[nd].append(kmer[1:])
    return adj_list

if __name__ == "__main__":
    # (1) Test decompose string
    #str = "CAATCCAAC"
    #print(decompose_string(str, 5))

    # (2) Test string spelled
    #l = [
    #'ACCGA',
    #'CCGAA',
    #'CGAAG',
    #'GAAGC',
    #'AAGCT']
    #print(string_spelled(l))

    # (3) Test adjacent list
    def format_adj_list(adj_list_dict):
        res = []
        for k, v in adj_list_dict.items():
            res.append(k + " -> " + ",".join(sorted(v)))
        return res

    #l = read_from_file("dataset_198_10.txt")
    #write_to_file(format_adj_list(graph_overlap(l)))

    # (4) Test de_bruijn_graph
    #str = 'AAGATTCTCTAAGA'
    #print("\n".join(format_adj_list(de_bruijn_graph(str, 4))))

    # (5) Test debruijn_graph_from_kmers
    l = ['GAGG',
    'CAGG',
    'GGGG',
    'GGGA',
    'CAGG',
    'AGGG',
    'GGAG']
    print("\n".join(format_adj_list(debruijn_graph_from_kmers(l))))
