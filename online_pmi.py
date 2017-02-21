from collections import defaultdict
import itertools as it
import sys, distances, igraph, utils
import numpy as np
import random

MAX_ITER = 5
tolerance = 0.001
infomap_threshold = 0.5
min_batch = 250000

fname = sys.argv[1]
dataname = sys.argv[2]
coding_option = sys.argv[3]

def igraph_clustering(matrix, threshold, method='infomap'):
    """
    Method computes Infomap clusters from pairwise distance data.
    """

    G = igraph.Graph()
    vertex_weights = []
    for i in range(len(matrix)):
        G.add_vertex(i)
        vertex_weights += [0]
    
    # variable stores edge weights, if they are not there, the network is
    # already separated by the threshold
    weights = None
    for i,row in enumerate(matrix):
        for j,cell in enumerate(row):
            if i < j:
                if cell <= threshold:
                    G.add_edge(i, j, weight=1-cell, distance=cell)
                    weights = 'weight'

    if method == 'infomap':
        comps = G.community_infomap(edge_weights=weights,
                vertex_weights=None)
        #comps = G.community_label_propagation(weights=weights,
        #        initial=None, fixed=None)

    elif method == 'ebet':
        dg = G.community_edge_betweenness(weights=weights)
        oc = dg.optimal_count
        comps = False
        while oc <= len(G.vs):
            try:
                comps = dg.as_clustering(dg.optimal_count)
                break
            except:
                oc += 1
        if not comps:
            print('Failed...')
            comps = list(range(len(G.sv)))
            input()
    elif method == 'multilevel':
        comps = G.community_multilevel(return_levels=False)
    elif method == 'spinglass':
        comps = G.community_spinglass()

    D = {}
    for i,comp in enumerate(comps.subgraphs()):
        vertices = [v['name'] for v in comp.vs]
        for vertex in vertices:
            D[vertex] = i+1

    return D
    
def infomap_concept_evaluate_scores(d, lodict, gop, gep):
    #fout = open("output.txt","w")
    average_fscore = []
    f_scores = defaultdict(list)
    bin_mat, n_clusters = None, 0
    for concept in d:
        ldn_dist_dict = defaultdict(lambda: defaultdict(float))
        langs = list(d[concept].keys())
        scores, cognates = [], []
        ex_langs = list(set(lang_list) - set(langs))
        for l1, l2 in it.combinations(langs, r=2):
            if d[concept][l1].startswith("-") or d[concept][l2].startswith("-"): continue
            w1, w2 = d[concept][l1], d[concept][l2]
            score = distances.nw(w1, w2, lodict=lodict, gp1=gop, gp2=gep)[0]
            score = 1.0 - (1.0/(1.0+np.exp(-score)))
            ldn_dist_dict[l1][l2] = score
            ldn_dist_dict[l2][l1] = ldn_dist_dict[l1][l2]
        distMat = np.array([[ldn_dist_dict[ka][kb] for kb in langs] for ka in langs])
        clust = igraph_clustering(distMat, infomap_threshold, method='infomap')
        predicted_labels = defaultdict()
        predicted_labels_words = defaultdict()
        for k, v in clust.items():
            predicted_labels[langs[k]] = v
            predicted_labels_words[langs[k]+":"+d[concept][langs[k]]] = v
        print(concept, "\n")
        print(predicted_labels_words)
        n_clusters += len(set(clust.values()))
        t = utils.dict2binarynexus(predicted_labels, ex_langs, lang_list)
        #print(concept, "\n",t)
        if bin_mat is None:
            bin_mat = t
        else:
            bin_mat += t
        
        print("No. of clusters ", n_clusters)
    return bin_mat

def calc_pmi(alignment_dict, char_list, pmidict, initialize=False):
    sound_dict = defaultdict(float)
    relative_align_freq = 0.0
    relative_sound_freq = 0.0
    count_dict = defaultdict(float)
    
    if initialize == True:
        for c1, c2 in product(char_list, repeat=2):
            if c1 == "-" or c2 == "-":
                continue
            count_dict[c1,c2] += 0.001
            count_dict[c2,c1] += 0.001
            sound_dict[c1] += 0.001
            sound_dict[c2] += 0.001
            relative_align_freq += 0.001
            relative_sound_freq += 0.002

    for alignment in alignment_dict:
        for a1, a2 in alignment:
            if a1 == "-" or a2 == "-":
                continue
            count_dict[a1,a2] += 1.0
            #count_dict[a2,a1] += 1.0
            sound_dict[a1] += 1.0
            sound_dict[a2] += 1.0
            relative_align_freq += 1.0
            relative_sound_freq += 2.0

    for a in count_dict.keys():
        m = count_dict[a]
        assert m>0

        num = np.log(m)-np.log(relative_align_freq)
        denom = np.log(sound_dict[a[0]])+np.log(sound_dict[a[1]])-(2.0*np.log(relative_sound_freq))
        val = num - denom
        count_dict[a] = val
    
    return count_dict

word_list = [line.strip().split(" ") for line in open(fname).readlines()]
char_list = [line.strip() for line in open("IPA.characters.list").readlines()]

data_dict = defaultdict(lambda: defaultdict())
concept_list, lang_list = [], []

lines = open(dataname).readlines()
for line in lines[1:]:
    line = line.replace("\n","")
    arr = line.split("\t")
    lang, concept = arr[0], arr[1]
    word = None
    if arr[4].startswith("-"): continue
    if coding_option.upper() == "SCA":
        word = arr[4]
    elif coding_option.upper() == "ASJP":
        word = utils.cleanASJP(arr[3])
    elif coding_option.upper() == "IPA":
        word = utils.cleanIPA(arr[2])
    data_dict[concept][lang] = word.split(",")[0]
    if lang not in lang_list: lang_list.append(lang)
    concept_list.append(concept)

pmidict = None
n_examples, n_updates, alpha = len(word_list), 0.0, 0.0

for n_iter in range(MAX_ITER):
    print("Iteration ", n_iter)
    algn_list = []
    for w1, w2 in word_list:
        s, alg = distances.nw(w1, w2, lodict=pmidict)
        algn_list.append(alg)
    pmidict = calc_pmi(algn_list, char_list, pmidict, initialize=False)

for k, v in pmidict.items():
    print(k, v)

#for n_iter in range(MAX_ITER):
#    random.shuffle(word_list)
#    print("Iteration ", n_iter)
#    for idx in range(0, n_examples, min_batch):
#        wl = word_list[idx:idx+min_batch]
#        eta = np.power(1.0/(n_updates+2), alpha)
#        algn_list = []
#        for w1, w2 in wl:
#            s, alg = distances.nw(w1, w2, lodict=pmidict)
#            algn_list.append(alg)
#        mb_pmi_dict = calc_pmi(algn_list, char_list, pmidict, initialize=False)
        #for k, v in mb_pmi_dict.items():
        #    pmidict_val = pmidict[k]
        #    pmidict[k] = (eta*v) + ((1.0-eta)*pmidict_val)
        #n_updates += 1
#    bin_mat = infomap_concept_evaluate_scores(data_dict, pmidict, -2.5, -1.75)

bin_mat = infomap_concept_evaluate_scores(data_dict, pmidict, -2.5, -1.75)
nchar, nlangs = np.array(bin_mat).shape

print("begin data;")
print("   dimensions ntax=", str(nlangs), "nchar=", str(nchar), ";\nformat datatype=restriction interleave=no missing= ? gap=-;\nmatrix\n")

for row, lang in zip(np.array(bin_mat).T, lang_list):
    #print(row,len(row), "\n")
    rowx = "".join([str(x) for x in row])
    print(lang, "\t", rowx.replace("2","?"))
print(";\nend;")

    
    
    
    
