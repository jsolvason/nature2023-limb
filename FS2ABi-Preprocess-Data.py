#!/usr/bin/env python
# coding: utf-8
# 1_20220902-Compare-ETS-Pbm-And-Chip-1.py

print('importing libraries...')

from bx.intervals.io import GenomicIntervalReader
from bx.bbi.bigwig_file import BigWigFile
import numpy as np
import time
import sys
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
import pickle 

import helpher_functions as hf

inDir = '../../data-downloaded/'

fn2genome={
    'GSM2218592_B_Ets1_R1_bow':'mm9', # b cells, Saelee P et al
    'GSM4110116_Ets1_S0_1_129929':'hg19', # t cells, McCarter AC et al
    'GSM3520734_NK_ETS1_chip':'hg38', # nk cells, Taveirne S et al
}

print('loading genomes...')
genome2chr2seq={
    'hg38':hf.faLoadGenome('hg38.fa'),
    'hg19':hf.faLoadGenome('hg19.fa'),
    'mm9' :hf.faLoadGenome('mm9.fa')
}

seq2aff=hf.loadAff('parsed_Ets1_8mers.txt')

# plot fast 
beta=False

for fn,genome in fn2genome.items():
    
    print(f'reading {fn}...')
    
    outDir=fn.split('.')[0]
    hf.mkdir_if_dir_not_exists('data-outputted/'+outDir)
        
    # read chip bigwig
    bw = BigWigFile( open( inDir+fn ,'rb') )

    # generate all possible ets 8mers
    allEtsKmers=hf.IupacToAllPossibleSequences("NNGGAWNN")

    fwdEts=set(['GGAA','GGAT'])
    revEts=set(['TTCC','ATCC'])

    doneList=[]

    for chrom,seq in genome2chr2seq[genome].items():
        if chrom in doneList: continue

        print(f'\t{chrom}',end=',  ')

        Kmer2Data={kmer:{'pbm-aff':seq2aff[kmer], #record affinity for that kmer 
                         'chip-signal-list':[]  , #list of chip signals since one kmer might occur many times throughout the chromosome 
                         'chip-nan-count':0} for kmer in allEtsKmers} #create dictionary for all possible ets 8mers

        lc=0

        # iterate over each kmer in the genome
        for start,kmer in enumerate(hf.get_kmers(seq,8)): #find all 8mers in the chromosome

            # if it is an unknown region of chromosome (ends) then don't use it
            if 'N' in kmer: continue 

            # if rev kmer, make forward
            if kmer[2:6] in revEts: 
                lc+=1
                kmer=hf.revcomp(kmer)

            # now that it is fwd ets, find it in the chip data and record signal as result
            if kmer[2:6] in hf.etsCores: 
                
                # get mean chip signal over 8mer
                result=bw.query(chrom, start, start+8, 1)
                if result==None: 
                    continue

                result=result[0]['mean']

                if np.isnan(result):
                    Kmer2Data[kmer]['chip-nan-count']+=1 #if kmer has no signal, then increment count of nan
                else:
                    Kmer2Data[kmer]['chip-signal-list'].append(round(result,5)) #add signal to the list of signals for that kmer

            if beta:
                if lc>10:
                    break

        bn=fn.split('/')[-1].split('.')[0]
        with open(f'data-outputted/{outDir}/{bn}__Kmer2Data__beta={beta}__chrom={chrom}.pydict.pickle','wb') as f:
            pickle.dump(Kmer2Data,f,protocol=pickle.HIGHEST_PROTOCOL) #dump dictionary into pickle file for that chromosome

        if beta: break

    print(f'\nDone with {fn}!!!')

print('done!')