import os, sys, re, os, math, shutil
from copy import deepcopy
from subprocess import *
from random import sample
import time
import sys
from pssm import pssm

def readMultiMotifMeme(filePath=None):
    memeFile = open(filePath, 'rb')

    allname = []
    allpp = []

    while 1:
        line = memeFile.readline()
        if not line:
            break

        if line.startswith("MEME"):
            continue

        if line.startswith("ALPHABET"):
            continue

        if line.startswith("MOTIF"):
            name = line.strip().split(" ")[1]
            allname.append(name)
            continue

        if line.startswith("URL"):
            continue

        if line.startswith("letter"):
            width = 1000000
            # width can be any number that couldn't be the length of a motif
            matrix = []
            l = line.strip().split(" ")
            length = l[3]
            if len(l) > 8:
                eValue = l[9]
            else:
                eValue = '0'
            nsites = l[7]
            width = int(l[5])
            genes = []
            print length, eValue, nsites, width
            for w in range(0,width):
                line=memeFile.readline()
                matrix.append([float(i) for i in line.strip().split()])
            pp = pssm(pssmFileName=None, biclusterName=name, nsites=nsites, eValue=eValue, pssm=matrix, genes=genes)
            allpp.append(pp)
            continue

    memeFile.close()
    return(allname, allpp)





def readMultiMotifMeme2(filePath=None):
    memeFile = open(filePath, 'rb')

    m = t = 1
    width = 10000000 
    # width can be any number that couldn't be the length of a motif
    allname = []
    allpp = []
    n=0

    for line in memeFile:
        n+=1
        if line.startswith("MEME"):
            continue
        elif line.startswith("ALPHABET"):
            continue
        elif line.startswith("MOTIF"):
            name = line.strip().split(" ")[1]
            allname.append(name)
            m = t = 1
        elif line.startswith("letter"):
            matrix = []
            l = line.strip().split(" ")
            length = l[3]
            eValue = l[9]
            nsites = l[7]
            width = l[5]
            genes = []
            m = 0
        elif m == 0 and t <= int(width):
            t += 1
            matrix.append([float(i) for i in line.strip().split()])

        if m == 0 and t == int(width) + 1:
            pp = pssm(pssmFileName=None, biclusterName=name, nsites=nsites, eValue=eValue, pssm=matrix, genes=genes)
            allpp.append(pp)
            m = t = 1

    memeFile.close()
    return(allname, allpp)


def readSingleMotifMeme(filePath=None):
    memeFile = open(filePath, 'rb')
    b = 0
    m = 0
    t = 1
    width = 0
    matrix=[]
    for line in memeFile:
        if b == 1:
            background = line.strip()
            b = 0
        elif line.startswith("MEME"):
            version = line.strip().split(" ")[2]
        elif line.startswith("ALPHABET"):
            alphabet = line.strip().split(" ")[1]
        elif line.startswith("strands"):
            strands = line.strip().split(" ")[1:]
        elif line.startswith("Background"):
            b = 1
        elif line.startswith("MOTIF"):
            name = line.strip().split(" ")[1]
            #name = '_'.join(name)
            # motif may have alterative names, just use the original name, no need to join them together
            m=0
            t=1
        elif line.startswith("letter"):
            l = line.strip().split(" ")
            length = l[3]
            eValue = l[9]
            nsites = l[7]
            width = l[5]
            genes = []
            m = 1
        elif m==1 and t <= int(width):
            t += 1
            #matrix += [float(i) for i in line.strip().split()]
            matrix.append([float(i) for i in line.strip().split()])
        else:
            continue

    #print(matrix)
    memeFile.close()
    pp = pssm(pssmFileName=None, biclusterName=name, nsites=nsites, eValue=eValue, pssm=matrix, genes=genes)
    #meme={}
    #meme['name'] = name
    #meme['nsites'] = nsites
    #meme['genes'] = genes
    #meme['matrix'] = matrix
    #meme['eValue'] = eValue
    return(name,pp)

