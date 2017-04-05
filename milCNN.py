import sys
import os
import numpy as np
import pdb
from keras.models import Sequential
from keras.layers.core import Dense, Dropout, Activation, Flatten, Merge
from keras.layers.normalization import BatchNormalization
from keras.layers.advanced_activations import PReLU, LeakyReLU
from keras.optimizers import SGD, RMSprop, Adadelta, Adagrad, Adam
from keras.layers import normalization
from keras.layers.convolutional import Convolution2D, MaxPooling2D
from keras.layers import LSTM, Bidirectional 
from keras.layers.embeddings import Embedding
from keras.layers.convolutional import Convolution2D, MaxPooling2D,Convolution1D, MaxPooling1D
from keras.callbacks import ModelCheckpoint, EarlyStopping
from keras.constraints import maxnorm
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from keras import objectives
from keras import backend as K
from keras.utils import np_utils
from sklearn.metrics import roc_curve, auc, roc_auc_score
_EPSILON = K.epsilon()
import random
import gzip
import pickle

def padding_sequence_new(seq, max_len = 101, repkey = 'N'):
    seq_len = len(seq)
    new_seq = seq
    if seq_len < max_len:
        gap_len = max_len -seq_len
        new_seq = seq + repkey * gap_len
    return new_seq

def read_rna_dict(rna_dict = 'rna_dict'):
    odr_dict = {}
    with open(rna_dict, 'r') as fp:
        for line in fp:
            values = line.rstrip().split(',')
            for ind, val in enumerate(values):
                val = val.strip()
                odr_dict[val] = ind
    
    return odr_dict

def get_6_trids():
    nucle_com = []
    chars = ['A', 'C', 'G', 'U']
    base=len(chars)
    end=len(chars)**6
    for i in range(0,end):
        n=i
        ch0=chars[n%base]
        n=n/base
        ch1=chars[n%base]
        n=n/base
        ch2=chars[n%base]
        n=n/base
        ch3=chars[n%base]
        n=n/base
        ch4=chars[n%base]
        n=n/base
        ch5=chars[n%base]
        nucle_com.append(ch0 + ch1 + ch2 + ch3 + ch4 + ch5)
    return  nucle_com 
def get_RNA_seq_concolutional_array(seq, motif_len = 4):
    seq = seq.replace('U', 'T')
    alpha = 'ACGT'
    #for seq in seqs:
    #for key, seq in seqs.iteritems():
    row = (len(seq) + 2*motif_len - 2)
    new_array = np.zeros((row, 4))
    for i in range(motif_len-1):
        new_array[i] = np.array([0.25]*4)
    
    for i in range(row-3, row):
        new_array[i] = np.array([0.25]*4)
        
    #pdb.set_trace()
    for i, val in enumerate(seq):
        i = i + motif_len-1
        if val not in 'ACGT':
            new_array[i] = np.array([0.25]*4)
            continue
        #if val == 'N' or i < motif_len or i > len(seq) - motif_len:
        #    new_array[i] = np.array([0.25]*4)
        #else:
        try:
            index = alpha.index(val)
            new_array[i][index] = 1
        except:
            pdb.set_trace()
        #data[key] = new_array
    return new_array

def get_embed_dim_new(embed_file):
    with open(embed_file) as f:
        pepEmbedding = pickle.load(f)
        
    embedded_dim = pepEmbedding[0].shape
    print embedded_dim
    n_aa_symbols, embedded_dim = embedded_dim
    print n_aa_symbols, embedded_dim
    # = embedded_dim[0]
    embedding_weights = np.zeros((n_aa_symbols + 1,embedded_dim))
    embedding_weights[1:,:] = pepEmbedding[0]
    
    return embedded_dim, pepEmbedding[0], n_aa_symbols

def split_overlap_seq(seq):
    window_size = 101
    overlap_size = 20
    #pdb.set_trace()
    bag_seqs = []
    seq_len = len(seq)
    if seq_len >= window_size:
        num_ins = (seq_len - 101)/(window_size - overlap_size) + 1
        remain_ins = (seq_len - 101)%(window_size - overlap_size)
    else:
        num_ins = 0
    bag = []
    end = 0
    for ind in range(num_ins):
        start = end - overlap_size
        if start < 0:
            start = 0
        end = start + window_size
        subseq = seq[start:end]
        bag_seqs.append(subseq)
    if num_ins == 0:
        seq1 = seq
        pad_seq = padding_sequence_new(seq1)
        bag_seqs.append(pad_seq)
    else:
        if remain_ins > 10:
            #pdb.set_trace()
            #start = len(seq) -window_size
            seq1 = seq[-window_size:]
            pad_seq = padding_sequence_new(seq1)
            bag_seqs.append(pad_seq)
    return bag_seqs
            
def get_6_nucleotide_composition(tris, seq, ordict):
    seq_len = len(seq)
    tri_feature = []
    k = len(tris[0])
    #tmp_fea = [0] * len(tris)
    for x in range(len(seq) + 1- k):
        kmer = seq[x:x+k]
        if kmer in tris:
            ind = tris.index(kmer)
            tri_feature.append(ordict[str(ind)])
        else:
            tri_feature.append(-1)
    #tri_feature = [float(val)/seq_len for val in tmp_fea]
        #pdb.set_trace()        
    return np.asarray(tri_feature)

def read_seq_graphprot(seq_file, label = 1):
    seq_list = []
    labels = []
    seq = ''
    with open(seq_file, 'r') as fp:
        for line in fp:
            if line[0] == '>':
                name = line[1:-1]
            else:
                seq = line[:-1].upper()
                seq = seq.replace('T', 'U')
                seq_list.append(seq)
                labels.append(label)
    
    return seq_list, labels

def get_RNA_concolutional_array(seq, motif_len = 4):
    seq = seq.replace('U', 'T')
    alpha = 'ACGT'
    #for seq in seqs:
    #for key, seq in seqs.iteritems():
    row = (len(seq) + 2*motif_len - 2)
    new_array = np.zeros((row, 4))
    for i in range(motif_len-1):
        new_array[i] = np.array([0.25]*4)
    
    for i in range(row-3, row):
        new_array[i] = np.array([0.25]*4)
        
    #pdb.set_trace()
    for i, val in enumerate(seq):
        i = i + motif_len-1
        if val not in 'ACGT':
            new_array[i] = np.array([0.25]*4)
            continue
        try:
            index = alpha.index(val)
            new_array[i][index] = 1
        except:
            pdb.set_trace()
        #data[key] = new_array
    return new_array

def load_graphprot_data(protein, train = True, path = '../data/GraphProt_CLIP_sequences/'):
    data = dict()
    tmp = []
    listfiles = os.listdir(path)
    
    key = '.train.'
    if not train:
        key = '.ls.'
    mix_label = []
    mix_seq = []
    mix_structure = []    
    for tmpfile in listfiles:
        if protein not in tmpfile:
            continue
        if key in tmpfile:
            if 'positive' in tmpfile:
                label = 1
            else:
                label = 0
            seqs, labels = read_seq_graphprot(os.path.join(path, tmpfile), label = label)
            #pdb.set_trace()
            mix_label = mix_label + labels
            mix_seq = mix_seq + seqs
    
    data["seq"] = mix_seq
    data["Y"] = np.array(mix_label)
    
    return data

def loaddata_graphprot(protein, train = True, ushuffle = True):
    #pdb.set_trace()
    data = load_graphprot_data(protein, train = train)
    label = data["Y"]
    rna_array = []
    #trids = get_6_trids()
    #nn_dict = read_rna_dict()
    for rna_seq in data["seq"]:
        #rna_seq = rna_seq_dict[rna]
        rna_seq = rna_seq.replace('T', 'U')
        
        seq_array = get_RNA_seq_concolutional_array(seq)
        #tri_feature = get_6_nucleotide_composition(trids, rna_seq_pad, nn_dict)
        rna_array.append(seq_array)
    
    return np.array(rna_array), label

def get_bag_data(data):
    bags = []
    seqs = data["seq"]
    labels = data["Y"]
    for seq in seqs:
        #pdb.set_trace()
        bag_seqs = split_overlap_seq(seq)
        #flat_array = []
        bag_subt = []
        for bag_seq in bag_seqs:
            tri_fea = get_RNA_seq_concolutional_array(bag_seq)
            bag_subt.append(tri_fea)

        bags.append(np.array(bag_subt))
    return bags, labels
    #for data in pairs.iteritems():
    #    ind1 = trids.index(key)
    #    emd_weight1 = embedding_rna_weights[ord_dict[str(ind1)]]


def custom_objective(y_true, y_pred):
    '''Just another crossentropy'''
    #y_true = K.clip(y_true, _EPSILON, 1.0-_EPSILON)
    #y_true = max(y_true)
    #y_armax_index = numpy.argmax(y_pred)
    y_new = K.clip(y_pred, _EPSILON, 1.0-_EPSILON)
    #y_new = max(y_pred)
    '''
    if y_new >= 0.5:
        y_new_label = 1
    else:
        y_new_label = 0
    cce = abs(y_true - y_new_label)
    '''
    cce = - (y_true * K.log(y_new) + (1 - y_true)* K.log(1-y_new))
    return cce

def set_cnn_model(input_dim = 4, input_length = 107):
    nbfilter = 16
    model = Sequential()
    #model.add(brnn)
    model.add(Convolution1D(input_dim=input_dim,input_length=input_length,
                            nb_filter=nbfilter,
                            filter_length=10,
                            border_mode="valid",
                            #activation="relu",
                            subsample_length=1))
    model.add(Activation('relu'))
    model.add(MaxPooling1D(pool_length=3))
    model.add(Flatten())
    model.add(Dropout(0.5))
    model.add(Dense(nbfilter, activation='relu'))
    model.add(Dropout(0.5))

    return model

def set_cnn_embed(n_aa_symbols, input_length, embedded_dim, embedding_weights, nb_filter = 16):
    #nb_filter = 64
    filter_length = 10
    dropout = 0.5
    model = Sequential()
    #pdb.set_trace()
    model.add(Embedding(input_dim=n_aa_symbols+1, output_dim = embedded_dim, weights=[embedding_weights], input_length=input_length, trainable = True))
    print 'after embed', model.output_shape
    model.add(Convolution1D(nb_filter, filter_length, border_mode='valid', init='glorot_normal'))
    model.add(Activation(LeakyReLU(.3)))
    model.add(MaxPooling1D(pool_length=3))
    model.add(Dropout(dropout))
    
    return model

def get_cnn_network_graphprot(rna_len = 501, nb_filter = 16):
    print 'configure cnn network'
    embedded_rna_dim, embedding_rna_weights, n_nucl_symbols = get_embed_dim('rnaEmbedding25.pickle')
    print 'symbol', n_nucl_symbols
    model = set_cnn_embed(n_nucl_symbols, rna_len, embedded_rna_dim, embedding_rna_weights, nb_filter = nb_filter)
    
    #model.add(Bidirectional(LSTM(2*nbfilter)))
    #model.add(Dropout(0.10))
    model.add(Flatten())
    model.add(Dense(nb_filter*50, activation='relu')) 
    model.add(Dropout(0.50))
    model.add(Dense(nb_filter*10, activation='sigmoid')) 
    model.add(Dropout(0.50))
    print model.output_shape
    
    return model
        
def get_all_embedding(protein):
    
    data = load_graphprot_data(protein)
    #pdb.set_trace()
    train_bags, label = get_bag_data(data)
    #pdb.set_trace()
    test_data = load_graphprot_data(protein, train = False)
    test_bags, true_y = get_bag_data(test_data) 
    
    return train_bags, label, test_bags, true_y

def run_network(model, total_hid, train_bags, test_bags, y_bags):
    model.add(Dense(1))
    model.add(Activation('softmax'))
    #categorical_crossentropy, binary_crossentropy
    #sgd = SGD(lr=0.01, decay=1e-6, momentum=0.9, nesterov=True)
    model.compile(loss='binary_crossentropy', optimizer='rmsprop')
    print 'model training'
    nb_epos= 10
    for iterate in range(nb_epos):
        print 'train epoch', iterate
        for training, y in zip(train_bags, y_bags):
            tmp_size = len(training)
            #pdb.set_trace()
            ys = np.array(tmp_size *[y])
            #model.fit(training, ys, batch_size = tmp_size, nb_epoch=1) np_utils.to_categorical(ys)
            ys = np_utils.to_categorical(ys)
            model.train_on_batch(training, ys)
            
    predictions = []
    for testing in test_bags:
        pred = model.predict_proba(testing)[:,1]
        predictions.append(max(pred))
    return predictions

def run_milcnn():
    data_dir = '../data/GraphProt_CLIP_sequences/'
    #trids =  get_6_trids()
    #ordict = read_rna_dict()
    #embedded_rna_dim, embedding_rna_weights, n_nucl_symbols = get_embed_dim_new('rnaEmbedding25.pickle')
    for protein in os.listdir(data_dir):
        
        protein = protein.split('.')[0]
        print protein
        train_bags, train_labels, test_bags, test_labels = get_all_embedding(protein)
        net =  set_cnn_model()
        
        #seq_auc, seq_predict = calculate_auc(seq_net)
        hid = 16
        predict = run_network(net, hid, train_bags, test_bags, train_labels)
        
        auc = roc_auc_score(test_labels, predict)
        print auc
        #run_mil_classifier(train_bags, train_labels, test_bags, test_labels)
run_milcnn()
#seq= 'TTATCTCCTAGAAGGGGAGGTTACCTCTTCAAATGAGGAGGCCCCCCAGTCCTGTTCCTCCACCAGCCCCACTACGGAATGGGAGCGCATTTTAGGGTGGTTACTCTGAAACAAGGAGGGCCTAGGAATCTAAGAGTGTGAAGAGTAGAGAGGAAGTACCTCTACCCACCAGCCCACCCGTGCGGGGGAAGATGTAGCAGCTTCTTCTCCGAACCAA'
#print len(seq)
#split_overlap_seq(seq)
