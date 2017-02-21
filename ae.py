# -*- coding: utf-8 -*-
from keras.models import Model
from collections import defaultdict
import glob, codecs, utils
import numpy as np
from keras import regularizers
from keras.layers import Input, Dense, Convolution2D, MaxPooling2D, UpSampling2D, Dropout, Reshape, Flatten, LSTM, RepeatVector, GRU, Bidirectional, SimpleRNN
from keras.layers.wrappers import TimeDistributed
#from keras.callbacks import TensorBoard
import multiprocessing
from sklearn import metrics
np.random.seed(1337)  # for reproducibility
from keras.preprocessing import sequence
import itertools as it
import sys, igraph, DPmeans
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn import mixture, preprocessing

unique_chars = []
max_word_len = 10
concept_list, lang_list = [], []


def pad_word(x):
    if len(x) > max_word_len:
        return x[:max_word_len]
    else:
        return x.center(max_word_len,"0")

def wrd_to_onehot(w):
    w2d = []
    n_chars = len(unique_chars)+1
    for x in w:
        temp = n_chars*[0]
        if x == "0":
            w2d.append(temp)
        else:
            i = unique_chars.index(x)+1
            temp[i] = 1
            w2d.append(temp)
    return np.array(w2d)

def writeNexus(dm,f):
    l=len(dm)
    f.write("#nexus\n"+
                  "\n")
    f.write("BEGIN Taxa;\n")
    f.write("DIMENSIONS ntax="+str(l)+";\n"
                "TAXLABELS\n"+
                "\n")
    i=0
    for ka in dm:
        f.write("["+str(i+1)+"] '"+ka+"'\n")
        i=i+1
    f.write(";\n"+
               "\n"+
               "END; [Taxa]\n"+
               "\n")
    f.write("BEGIN Distances;\n"
            "DIMENSIONS ntax="+str(l)+";\n"+
            "FORMAT labels=left;\n")
    f.write("MATRIX\n")
    i=0
    for ka in dm:
        row="["+str(i+1)+"]\t'"+ka+"'\t"
        for kb in dm:
            row=row+str(dm[ka][kb])+"\t"
        f.write(row+"\n")
    f.write(";\n"+
    "END; [Distances]\n"+
    "\n"+
    "BEGIN st_Assumptions;\n"+
    "    disttransform=NJ;\n"+
    "    splitstransform=EqualAngle;\n"+
    "    SplitsPostProcess filter=dimension value=4;\n"+
    "    autolayoutnodelabels;\n"+
        "END; [st_Assumptions]\n")
    f.flush()

def read_data(fname):
    lines = open(fname).readlines()
    data_dict = defaultdict(lambda: defaultdict())
    ld = defaultdict(lambda: defaultdict())
    
    lwords = []
    for line in lines[1:]:
        line = line.replace("\n","")
        arr = line.split("\t")
        lang, concept = arr[0], arr[1]
        word = None
        if arr[4].startswith("-"):
            #data_dict[concept][lang] = "NA"
            continue
        if coding_option.upper() == "IPA":
            ipa_words = utils.cleanIPA(arr[2]).split(",")
            lwords += ipa_words
            data_dict[concept][lang] = ipa_words[0]
            ld[lang][concept] = data_dict[concept][lang]
            for x in ipa_words:
                for ch in x:
                    if ch not in unique_chars:
                        unique_chars.append(ch)
        elif coding_option.upper() == "ASJP":
            asjp_words = arr[3].split(", ")
            clean_asjp_words = [utils.cleanASJP(z) for z in asjp_words]
            lwords += clean_asjp_words
            data_dict[concept][lang] = clean_asjp_words[0]
            ld[lang][concept] = data_dict[concept][lang]
            for x in clean_asjp_words:
                for ch in x:
                    if ch not in unique_chars:
                        unique_chars.append(ch)
        if lang not in lang_list: lang_list.append(lang)
        concept_list.append(concept)
    print(unique_chars)
    return data_dict, ld, lwords

train_words = []
test_words = []
test_concepts = []

dataset = sys.argv[1]
coding_option = sys.argv[2]
print(dataset)

d, ld, train_words = read_data(dataset)
#train_words = list(set(train_words))

for concept in d:
    test_concepts.append(concept)
    test_words.append(d[concept])


train_phonetic, test_phonetic = [], []

temp = [x.replace("\n","").split(" ") for x in open("gondi_ldn_ipa1.txt").readlines()]
source_phonetic, target_phonetic = [], []

for src, tgt in temp:
    #print(src, tgt)
    source_phonetic.append(wrd_to_onehot(src))
    target_phonetic.append(wrd_to_onehot(tgt))

for x in train_words:
    onehotp1 = wrd_to_onehot(x)
    train_phonetic.append(onehotp1)
    
n_dim = len(unique_chars)+1
latent_dim = 32

train_phonetic = np.array(train_phonetic)

train_phonetic = sequence.pad_sequences(train_phonetic, maxlen=max_word_len)

source_phonetic = sequence.pad_sequences(source_phonetic, maxlen=max_word_len)
target_phonetic = sequence.pad_sequences(target_phonetic, maxlen=max_word_len)

print(train_phonetic.shape)

inputs = Input(shape=(max_word_len, n_dim))
#encoded = GRU(latent_dim)(inputs)
encoded = Bidirectional(GRU(latent_dim), merge_mode="concat")(inputs)

decoded = RepeatVector(max_word_len)(encoded)
decoded = GRU(n_dim, return_sequences=True)(decoded)
decoded = TimeDistributed(Dense(n_dim, activation="softmax"))(decoded)

sequence_autoencoder = Model(inputs, decoded)
encoder = Model(input=inputs, output=encoded)

sequence_autoencoder.compile(optimizer='rmsprop', loss='categorical_crossentropy')
encoder.compile(optimizer='rmsprop', loss='categorical_crossentropy')

sequence_autoencoder.summary()
sequence_autoencoder.fit(train_phonetic, train_phonetic,
                nb_epoch=50,
                batch_size=32)

dm = defaultdict(lambda: defaultdict(float))
repr_phonetic = encoder.predict(train_phonetic)
repr_dict = defaultdict()

for w, wv in zip(train_words, repr_phonetic):
    repr_dict[w] = wv
    #print(w, wv)

langs_list = list(set(lang_list))
concepts = list(set(concept_list))

print("\n", langs_list, "\n", len(langs_list))

fout = open(dataset+"_"+coding_option+"_"+str(latent_dim)+".word_vectors.txt","w")
for c in concepts:
    lang_vector, vocabulary = [], []
    for lang in langs_list:
        if concept not in d:
            continue
        if lang not in d[c]:
            continue
        w = d[c][lang]
        if w == "NA":
            continue
        #print(c, lang, w, repr_dict[w])
        word_vec = []
        for x in repr_dict[w]:
            word_vec.append(str(x))
        fout.write(c+"\t"+lang+"\t"+w+"\t"+"\t".join(word_vec)+"\n")
        lang_vector.append(repr_dict[w])
        vocabulary.append(w)
    wv = np.array(lang_vector)
    #wv = preprocessing.scale(wv)
    #pca = PCA(n_components=2)
    #mds = MDS(n_components=2)
    #np.set_printoptions(suppress=True)
    #Y = wv
    #Y = pca.fit_transform(wv)

    #plt.scatter(Y[:, 0], Y[:, 1])
    #for label, x, y in zip(vocabulary, Y[:, 0], Y[:, 1]):
    #    plt.annotate('$%s$' %label, xy=(x, y), xytext=(0, 0), textcoords='offset points')
    #plt.show()



for l1, l2 in it.combinations(langs_list, r=2):
    #print(l1, l2)
    n_concepts = 0.0
    num_sim = 0.0
    for concept in concepts:
        if concept not in ld[l1] or concept not in ld[l2]:
            continue
        
        x1 = ld[l1][concept]
        x2 = ld[l2][concept]
        if d[concept][l1] == "NA" or d[concept][l2] == "NA":
            continue
        else:
            r1, r2 = repr_dict[x1], repr_dict[x2]
            r1, r2 = r1.flatten(), r2.flatten()
            cos_sim = np.dot(r1, r2)/(np.linalg.norm(r1)*np.linalg.norm(r2))
            print(concept, l1, l2, x1, x2, cos_sim)
            #cos_sim = np.exp(-np.sum(np.square(r1-r2)))
            #cos_sim = cos_sim/(np.linalg.norm(r1) **2 + np.linalg.norm(r2) **2 - cos_sim)
            dm[l1][l2] += 1.0-cos_sim
            dm[l2][l1] += 1.0-cos_sim
            n_concepts += 1.0

    denom_sim = 0.0
    n_denom = 0.0
    dm[l1][l2] = dm[l1][l2]/n_concepts
    dm[l2][l1] = dm[l1][l2]
    continue
    for c1, c2 in it.product(concepts, concepts):
        if c1 not in ld[l1] or c2 not in ld[l2]:
            continue
        x1 = ld[l1][c1]
        x2 = ld[l2][c2]
        if d[c1][l1] == "NA" or d[c2][l2] == "NA":
            continue
        else:
            r1, r2 = repr_dict[x1], repr_dict[x2]
            r1, r2 = r1.flatten(), r2.flatten()
            cos_sim = np.dot(r1, r2)/(np.linalg.norm(r1)*np.linalg.norm(r2))
            denom_sim += 1.0-(cos_sim+1.0)/2.0
            n_denom += 1.0
    
    denom_sim = denom_sim/n_denom
    dm[l1][l2] = num_sim/denom_sim
    dm[l2][l1] = dm[l1][l2]
    
fout = codecs.open(dataset+"_"+coding_option+".srctgt.LSTM_concat_bidir.nex","w","utf-8")
writeNexus(dm, fout)
fout.close()

