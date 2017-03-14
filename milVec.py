import os
import pdb
import gzip
import numpy as np
#from pandas.core import window
import misvm       

def padding_sequence_new(seq, max_len = 101, repkey = 'N'):
    seq_len = len(seq)
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
    num_ins = len(seq)/window_size
    remain_ins = len(seq)%window_size
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
               seq_list.append(seq)
               labels.append(label)

   return seq_list, labels

def load_graphprot_data(protein, train = True, path = '/home/panxy/eclipse/rna-protein/data/GraphProt_CLIP_sequences/'):
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

def get_bag_data(data, tris, ordict):
    bags = []
    seqs = data["seq"]
    
    for seq in seqs:
        bag_seqs = split_overlap_seq(seq)
        flat_array = []
        for bag_seq in bag_seqs:
            tri_fea = get_6_nucleotide_composition(tris, bag_seq, ordict)
        flat_array = np.ndarray.flatten(tri_fea)
        bags.append(flat_array)
    return bags

def get_all_embedding(protein):
    trids =  get_6_trids()
    ord_dict = read_rna_dict()
    embedded_rna_dim, embedding_rna_weights, n_nucl_symbols = get_embed_dim_new('rnaEmbedding25.pickle')
    data, label = loaddata_graphprot(protein)
   
    test_data, true_y = loaddata_graphprot(protein, train = False)
      
    for key, val in pairs.iteritems():
        ind1 = trids.index(key)
        emd_weight1 = embedding_rna_weights[ord_dict[str(ind1)]]

def run_mil_classifier():
    classifiers = {}
    classifiers['MissSVM'] = misvm.MissSVM(kernel='linear', C=1.0, max_iters=10)
    classifiers['sbMIL'] = misvm.sbMIL(kernel='linear', eta=0.1, C=1.0)
    classifiers['SIL'] = misvm.SIL(kernel='linear', C=1.0)

    # Train/Evaluate classifiers
    accuracies = {}
    for algorithm, classifier in classifiers.items():
        classifier.fit(train_bags, train_labels)
        predictions = classifier.predict(test_bags)
        accuracies[algorithm] = np.average(test_labels == np.sign(predictions))

    for algorithm, accuracy in accuracies.items():
        print '\n%s Accuracy: %.1f%%' % (algorithm, 100 * accuracy)

def run_milvec():
    data_dir = '/home/panxy/eclipse/rna-protein/data/GraphProt_CLIP_sequences/'
    for protein in os.listdir(data_dir):
        
        protein = protein.split('.')[0]
        print protein
       get_all_embedding(protein)
#seq= 'TTATCTCCTAGAAGGGGAGGTTACCTCTTCAAATGAGGAGGCCCCCCAGTCCTGTTCCTCCACCAGCCCCACTACGGAATGGGAGCGCATTTTAGGGTGGTTACTCTGAAACAAGGAGGGCCTAGGAATCTAAGAGTGTGAAGAGTAGAGAGGAAGTACCTCTACCCACCAGCCCACCCGTGCGGGGGAAGATGTAGCAGCTTCTTCTCCGAACCAA'
#print len(seq)
#split_overlap_seq(seq)
