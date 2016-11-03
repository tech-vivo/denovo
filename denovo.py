# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 10:38:36 2016

@author: will mcfadden

This is my package that does all the denovo alignment.
I have tried to comment the functions pretty explicitly.
See the accompanying readme for details on how to use this.
"""

from math import floor
import sys, os.path


def readFastaReads(filename):
    """ Reads Fasta files with different reads seperated by a > lines
     
        returns: list of DNA sequences
    """
    reads = []
    read_frag = ''
    with open(filename, 'r') as f:
        for line in f:
            # ignore header lines
            if not line[0] == '>':
                read_frag += line.rstrip()
            else:
                if(len(read_frag)>0):
                    reads.append(read_frag)
                read_frag = ''
    if(len(read_frag)>0):
        reads.append(read_frag)
    return reads

def overlap(a, b, overlap_n):
    """ Returns length of longest suffix of 'a' matching
        a prefix of 'b' that is at least 'min_length'
        characters long.  If no such overlap exists,
        return 0. 
        
        Start from left of a and return first occurence the first n letters of b
        See if the rest of a matches the start of b
        If not move to next index and repeat
    """
    start = 0  
    while True:
        start = a.find(b[:overlap_n], start) 
        if start == -1:  
            return 0
        if b.startswith(a[start:]):
            return len(a)-start
        start += 1  
        
def make_kmers(seq,k):
    """ Builds a kmer dictionary of all the sequences of length k
        
        For every k-mer contained in the list of sequences seq
        we save the set of all sequences and the corresponding
        index of that sequence.
    
    """
    kmers=dict()
    for i,s in enumerate(seq):
        for j in range(len(s)-k):
            if(s[j:j+k] in kmers):
                kmers[s[j:j+k]].add((i,s))
            else:
                kmers[s[j:j+k]]=set()
                kmers[s[j:j+k]].add((i,s))
    return kmers

def overlap_all_pairs_tricky(seq):
    """ This is an attempt to do a tricky kmer indexing and lookup before comparing
        each sequence.  See the ipython notebook for a performance comparison 
        with overlap_all_pairs_easy(seq)
        
        Rationale: This will allow us to reduce the work from O(N^2) to O(N+n^2)
        Where N is the total number of sequences and n < N is just the 
        sequences where we already know that there is some overlap
        
        returns:
        start - the index of the first read (ie no read maps to it)
        followers - a dictionary with key equal to the read index
        and value equal to a tuple of (next_read, location_of_overlap)
            next_read is the immediate following read 
            location_of_overlap is how far from the end the next_read begins
    """
    min_n = int(floor(min([len(s) for s in seq])/2))
    k=min_n-1
    kmers = make_kmers(seq,k)
    followers = {}
    no_preceder=set(range(len(seq)))
    for i,a in enumerate(seq):
        suff = a[-1-k:-1]
        max_spt = 0
        for j,b in kmers[suff]:
            n = int(floor(min(len(a),len(b))/2))
            if(a!=b):
                spt = overlap(a,b,n)
                if(spt>max_spt):
                    followers[i]=(j,spt)
                    if j in no_preceder:
                        no_preceder.remove(j)
    start = list(no_preceder)[0]
    return start,followers
    
def overlap_all_pairs_easy(seq):
    """    
    Brute force through all pairs in 'seq' and calls overlap() on each to 
    determine if there is an overlap between them
    
    returns:
        start - the index of the first read (ie no read maps to it)
        followers - a dictionary with key equal to the read index
        and value equal to a tuple of (next_read, location_of_overlap)
            next_read is the immediate following read 
            location_of_overlap is how far from the end the next_read begins
    """
    
    followers = {}
    no_preceder=set(range(len(seq)))
    for i,a in enumerate(seq):
        max_spt = 0        
        for j,b in enumerate(seq):
            n = int(floor(min(len(a),len(b))/2))
            if(a!=b):
                spt = overlap(a,b,n)
                if(spt>max_spt):
                    followers[i]=(j,spt)
                    if j in no_preceder:
                        no_preceder.remove(j)
    start = list(no_preceder)[0]
    return start,followers
    
def conjoiner(start,followers,reads):
    """    
    Starts at the index 'start' and then walks through the 'followers' dictionary
    adding the appropriate fragment of that read to the 'conjoined' string
    
    returns:
        conjoined - the string that consists of all conjoined reads
    """
    
    i = start
    conjoined = reads[i]
    while True:
        if i in followers:
            next_up = followers[i]
            new_read = reads[next_up[0]]
            to_add = new_read[next_up[1]:]
            conjoined+=to_add
            i = next_up[0]
        else:
            break
    return conjoined
            
def assemble(filename,use_tricky=0):
    """ wrapper that calls the subfunctions above to go from filename to assembly
    
    """
    reads = readFastaReads(filename)
    if(use_tricky):
        start,followers = overlap_all_pairs_easy(reads)
    else:
        start,followers = overlap_all_pairs_easy(reads)
    return conjoiner(start,followers,reads)










"""
  This is just the logic for the command line call
"""
if __name__ == "__main__":
    if len(sys.argv)<2 or len(sys.argv)>3:
        print "\nThat's not gonna work. \n   usage: $denovo inputfile [outputfile optional] \n\n"
        sys.exit()
    
    file_in = sys.argv[1]
    
    if os.path.isfile(file_in):
        output = assemble(file_in)
    else:
        print "\nSorry input file not found! \n\n"
        sys.exit()

    if len(sys.argv)==3:
        file_out= sys.argv[2]
        with open(file_out, "w") as text_file:
            text_file.write(output)
    else:
        print output
        