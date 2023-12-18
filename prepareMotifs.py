from pssm import pssm
import cPickle, gzip, os, sys, re, os, math, shutil
from copy import deepcopy
from subprocess import *
from random import sample
from multiprocessing import Pool, cpu_count, Manager
import time
import sys
import readmeme as rmm
from versionRecord import getProgramsVersion

'''
### get version of softwares
getProgramsVersion()
'''


'''
### read in JASPAR: all .meme files in one directory, each .meme file has only one motif, transform them into one single .pkl file
dmpjaspar={}
#memeDir = "./Motifs/JASPAR/"
memeDir = "./Motifs/test/"
filenames = os.listdir(memeDir)
for filename in filenames:
    filePath = memeDir + filename
    rt = rmm.readSingleMotifMeme(filePath)
    dmpjaspar[rt[0]]=rt[1]
print("Motifs Number in JASPAR: " + str(len(dmpjaspar)))
#pklFile = open('./PSSMs/'+'JASPAR2022_CORE_vertebrates_non-redundant_pfms.pkl','wb')
pklFile = open('./PSSMs/'+'test-JASPAR2022_CORE_vertebrates_non-redundant_pfms.pkl','wb')
cPickle.dump(dmpjaspar, pklFile)
pklFile.close()
'''


### read in SELEX: a single .meme file which has more than one motifs, transform them into one single .pkl file
dmpselex={}
memeFile = "./Motifs/SELEX/all-selex-motifs.meme"
names = []
contents = []
(names, contents) = rmm.readMultiMotifMeme(memeFile)
for i in range(0,len(names)):
    dmpselex[names[i]] = contents[i]
print("Motifs Number in SELEX: " + str(len(names)))
pklFile = open('./PSSMs/'+'SELEX_all-selex-motifs.pkl','wb')
cPickle.dump(dmpselex, pklFile)
pklFile.close()

print(dmpselex.items())



'''
###-------Test for output a single meme from pkl-------BEGIN
pkl =  open('./test.pkl', 'rb')
test = cPickle.load(pkl)
pkl.close()

mgr = Manager()
writeMotifs = mgr.dict(test)
motif = 'MA2001.1'
meme = open(motif + '.meme', 'w')
memeHeader = ''
memeHeader += 'MEME version 3.0\n\n'
memeHeader += 'ALPHABET= ACGT\n\n'
memeHeader += 'strands: + -\n\n'
memeHeader += 'Background letter frequencies (from dataset with add-one prior applied):\n'
memeHeader += 'A 0.250 C 0.250 G 0.250 T 0.250\n\n'
meme.write(memeHeader)
meme.write(writeMotifs[motif].getMemeFormatted())
meme.close()
###-------Test for output a single meme from pkl---------END
'''






