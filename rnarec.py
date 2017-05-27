def padding_sequence(seq, max_len = 2695, repkey = 'N'):
    seq_len = len(seq)
    if seq_len < max_len:
        gap_len = max_len -seq_len
        new_seq = seq + repkey * gap_len
    else:
        new_seq = seq[:max_len]
    return new_seq

def local_ushuffle(seq, dishuff = True):
    '''
    shuffle dinucletide
    '''

    if dishuff:
        return_seq = ushuffle.shuffle(seq, len(seq), 2)
    else:
        l = list(seq)
        random.shuffle(l)
        return_seq = ''.join(l)        
        
    return return_seq

def get_rna_seqs_rec(rnas, rna_seq_dict, rna_nax_len, trids, nn_dict):

    label = []
    rna_array = []
    for rna in rnas:
        rna_seq = rna_seq_dict[rna]
        rna_seq = rna_seq.replace('T', 'U')
        rna_seq_pad = padding_sequence(rna_seq, max_len = rna_nax_len, repkey = 'N')
        tri_feature = get_6_nucleotide_composition(trids, rna_seq_pad, nn_dict)
        rna_array.append(tri_feature)
        label.append(1)
        if ushuffle:
            shuffle_rna_seq = local_ushuffle(rna_seq)
            shuffle_rna_seq_pad = padding_sequence(shuffle_rna_seq, max_len = rna_nax_len, repkey = 'N')
            #onehot_rna = get_RNA_concolutional_array(shuffle_rna_seq_pad)
            tri_feature_shu = get_6_nucleotide_composition(trids, shuffle_rna_seq_pad, nn_dict)
            label.append(0)
            rna_array.append(tri_feature_shu)
    
    return np.array(rna_array), np.array(label)

def read_fasta_file(fasta_file):
    seq_dict = {}    
    fp = open(fasta_file, 'r')
    name = ''
    for line in fp:
        #let's discard the newline at the end (if any)
        line = line.rstrip()
        #distinguish header from sequence
        if line[0]=='>': #or line.startswith('>')
            #it is the header
            name = line[1:] #discarding the initial >
            seq_dict[name] = ''
        else:
            #it is sequence
            seq_dict[name] = seq_dict[name] + line.upper()
    fp.close()
    
    return seq_dict

def run_rnacomend(datadir = 'data/', ushuffle = True):
    fw = open('result_file_rbp67', 'w')
    pair_file = datadir + 'interactions_HT.txt'
    #rbp_seq_file = datadir + 'rbps_HT.fa'
    rna_seq_file = datadir + 'utrs.fa'
     
    rna_seq_dict = read_fasta_file(rna_seq_file)
    rna_len = [] # 2695
    for val in rna_seq_dict.values():
        rna_len.append(len(val))
    
    rna_len.sort()
    rna_nax_len = rna_len[int(len(rna_len)*0.8)] 
    seq_hid = 16
    label = []
    rna_array = []
    protein_pair = {}
    trids = get_6_trids()
    nn_dict = read_rna_dict()
    with open(pair_file, 'r') as fp:
        for line in fp:
            values = line.rstrip().split()
            protein = values[0]
            rna = values[1]
            protein_pair.setdefault(protein, []).append(rna)
    
    for protein, rnas in protein_pair.iteritems():
        print protein
        #pdb.set_trace()
        fw.write(protein + '\t')
        if len(rnas) < 2000:
            continue
        data, label = get_rna_seqs_rec(rnas, rna_seq_dict, rna_nax_len, trids, nn_dict)
        seq_net = get_cnn_network_graphprot(rna_len = rna_nax_len - 5, nb_filter = seq_hid)

        training_val_indice, train_val_label, test_indice, test_label = split_training_validation(label)
        train_val = data[training_val_indice]
        train_val_label = label[training_val_indice]  
        test_data = data[test_indice]
        test_label = label[test_indice]
        #pdb.set_trace()
        training_indice, training_label, val_indice, val_label = split_training_validation(train_val_label)
        cnn_train = train_val[training_indice]
        training_label = train_val_label[training_indice]  
        cnn_validation = train_val[val_indice]
        validation_label = train_val_label[val_indice]
                 
        y, encoder = preprocess_labels(training_label)
        val_y, encoder = preprocess_labels(validation_label, encoder = encoder) 
        print 'predicting'    
        model_name = 'model/' + protein +'.pickle'
        seq_auc, seq_predict = calculate_auc(seq_net, seq_hid, cnn_train, test_data, test_label, y, validation = cnn_validation,
                                              val_y = val_y, model_name = model_name)
        print str(seq_auc)
        fw.write( str(seq_auc) +'\n')
        mylabel = "\t".join(map(str, test_label))
        myprob = "\t".join(map(str, seq_predict))  
        fw.write(mylabel + '\n')
        fw.write(myprob + '\n')
        
    fw.close()
