import os
import pdb
import gzip
#from pandas.core import window
       

def padding_sequence_new(seq, max_len = 101, repkey = 'N'):
    seq_len = len(seq)
    if seq_len < max_len:
        gap_len = max_len -seq_len
        new_seq = seq + repkey * gap_len
    return new_seq

        
def split_overlap_seq(seq):
    window_size = 101
    overlap_size = 20
    #pdb.set_trace()
    num_ins = len(seq)/window_size
    remain_ins = len(seq)%window_size
    bag = []
    for ind in range(num_ins):
        start = ind *window_size - overlap_size
        if start < 0:
            start = 0
        end = start + window_size
        subseq = seq[start:end]
        print subseq
    if num_ins == 0:
        seq1 = seq
        pad_seq = padding_sequence_new(seq1)
    else:
        #if remain_ins >= overlap_size:
            #pdb.set_trace()
        start = len(seq) -window_size
        seq1 = seq[start:]
        pad_seq = padding_sequence_new(seq1)
        print pad_seq

        
seq= 'TTATCTCCTAGAAGGGGAGGTTACCTCTTCAAATGAGGAGGCCCCCCAGTCCTGTTCCTCCACCAGCCCCACTACGGAATGGGAGCGCATTTTAGGGTGGTTACTCTGAAACAAGGAGGGCCTAGGAATCTAAGAGTGTGAAGAGTAGAGAGGAAGTACCTCTACCCACCAGCCCACCCGTGCGGGGGAAGATGTAGCAGCTTCTTCTCCGAACCAA'
print len(seq)
split_overlap_seq(seq)
