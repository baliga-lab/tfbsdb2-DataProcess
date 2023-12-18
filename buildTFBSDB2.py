#################################################################
# @Program: buildTFBSDB.py                                      #
# @Version: 1                                                   #
# @Author: Christopher L Plaisier, PhD                          #
# @Sponsored by:                                                #
# Nitin Baliga, ISB                                             #
# Institute for Systems Biology                                 #
# 401 Terry Ave North                                           #
# Seattle, Washington  98109-5234                               #
# (216) 732-2139                                                #
# @Also Sponsored by:                                           #
# Luxembourg Systems Biology Grant                              #
# American Cancer Society Postdoctoral Fellowship               #
#                                                               #
# If this program is used in your analysis please mention who   #
# built it. Thanks. :-)                                         #
#                                                               #
# Copyrighted by Chris Plaisier  2/15/2013                      #
# Modifyed    by Yaqiao Li       12/09/2022                     #
#################################################################

###############
### IMPORTS ###
###############
from pssm import pssm
import gzip, os, sys, re, os, math, shutil
from copy import deepcopy
from subprocess import *
from random import sample
from multiprocessing import Pool, cpu_count, Manager
import time
import _pickle as cPickle
import versionRecord

versionRecord.getProgramsVersion()

# Boot a manager for the shared memory objects
mgr = Manager()

## Multipocessing function for FIMO
def runFimo(motif):
    if not os.path.exists('tmp/pssm_'+motif.replace('$','_')+'.meme'):
        outFile = open('tmp/pssm_'+motif.replace('$','_')+'.meme','w')
        # Header crap
        memeHeader = ''
        memeHeader += 'MEME version 3.0\n\n'
        memeHeader += 'ALPHABET= ACGT\n\n'
        # Here is where we tell it what strand: for miRNAs this would just be '+'
        memeHeader += 'strands: + -\n\n'
        memeHeader += 'Background letter frequencies (from dataset with add-one prior applied):\n'
        memeHeader += 'A 0.250 C 0.250 G 0.250 T 0.250\n\n'
        outFile.write(memeHeader)
        outFile.write(writeMotifs[motif].getMemeFormatted())
        outFile.close()

    for chrNum in ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']:
    #for chrNum in ['chrY']:
        if not os.path.exists('motifHits2/motifHits_'+motif.replace('$','_')+'_'+str(chrNum)+'.pkl'):
            tmp = fimo(queryFile='tmp/pssm_'+motif.replace('$','_')+'.meme', seqFile='footprintSeqs/'+chrNum+'.fa')
            pklFile = open('motifHits2/motifHits_'+motif.replace('$','_')+'_'+str(chrNum)+'.pkl','wb')
            cPickle.dump(tmp, pklFile)
            pklFile.close()

## Run FIMO on specified files
def fimo(queryFile=None, seqFile=None):
    fimoArgs = '--verbosity 4 --bgfile bgFile.meme -text --thresh 1e-5 '
    fimoArgs += str(queryFile)+' '+str(seqFile)
    print(fimoArgs)
    errOut = open('stderr.out','a')
    fimoProc = Popen("/users/yli/software/anaconda/envs/tfbsdb2022/bin/fimo " + fimoArgs, shell=True,stdout=PIPE,stderr=errOut)
    errOut.close()
    stdoutput = fimoProc.communicate()[0]
    output = stdoutput.decode('utf-8').split('\n')
    # Write out to a file
    #outFile = open('fimo_out.txt','w')
    #outFile.write('\n'.join(output))
    #outFile.close()
    # Return results like cMonkey
    res1 = []
    a = output.pop(0)
    # a is the first line, header of fimo output. 
    #print(output[0])#yli: if fimo found no hits in this chr, this line will cause error: index out of range. Because the output is empty.
    #pattern name    sequence name    start    stop    strand    score    p-value    q-value    matched sequence
    for line1 in output:
        splitUp = line1.split('\t')
        if len(splitUp)==10:
            splitUp2 = splitUp[1].split('.')
            """ For running on footprints only
            chrom = splitUp2[0]
            splitUp3 = splitUp2[1].split('_')
            start = int(splitUp3[0])+(int(splitUp[2])-1)
            stop = int(splitUp3[0])+(int(splitUp[3])-1)
            """
            chrom = seqFile.split('/')[1].rstrip('.fa')
            start = int(splitUp[3])
            stop = int(splitUp[4])
            res1.append({'chr': chrom, 'start': start, 'stop': stop, 'strand': splitUp[5], 'score':splitUp[6], 'p.value':float(splitUp[7]), 'match.sequence':splitUp[9]})
    #print(res1[0])
    return res1

#########################################################
# Load up motifs and turn them into MEME formatted file #
#########################################################
# JASPAR
#pklFile = open('PSSMs/jasparCoreVertebrata.pkl','rb')
#pklFile = open('PSSMs/test-JASPAR2022_CORE_vertebrates_non-redundant_pfms.pkl','rb')
pklFile = open('PSSMs/JASPAR2022_CORE_vertebrates_non-redundant_pfms.pkl','rb')
jaspar = cPickle.load(pklFile)
pklFile.close()

# TransFac
pklFile = open('PSSMs/transfac_2012.1_PSSMs_vertabrate.pkl','rb')
transfac = cPickle.load(pklFile)
pklFile.close()

# SELEX
#pklFile = open('PSSMs/selexPSSMsNonRedundant.pkl','rb')
pklFile = open('PSSMs/SELEX_all-selex-motifs.pkl','rb')
selex = cPickle.load(pklFile)
pklFile.close()

# Uniprobe
pklFile = open('PSSMs/uniprobePSSMsNonRedundant.pkl','rb')
uniprobe = cPickle.load(pklFile)
pklFile.close()

compendium = dict(jaspar.items())
compendium.update(dict(transfac.items()))
compendium.update(dict(selex.items()))
compendium.update(dict(uniprobe.items()))

#compendium = dict(transfac.items())
#compendium = dict(jaspar.items() + transfac.items() + selex.items() + uniprobe.items())
writeMotifs = mgr.dict(compendium)
del compendium
print(len(jaspar))
print(len(transfac))
print(len(selex))
print(len(uniprobe))
print(len(writeMotifs))

# Prep for run
if not os.path.exists('tmp'):
    os.mkdir('tmp')

# Run fimo to get the target sites
print('Running FIMO...')

# Multiprocessing version
cpus = 64
print('There are '+ str(cpus)+ ' CPUs avialable.')
pool = Pool(processes=cpus)
pool.map(runFimo, writeMotifs.keys())
#doEm = ['ERR188698']
#pool.map(runFimo, doEm)
