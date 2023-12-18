#################################################################
# @Program: compileSequences.py                                 #
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
# Copyrighted by Chris Plaisier  1/15/2013                      #
#################################################################


import os, gzip
from copy import deepcopy
import tarfile

################################################
# 1. Read in footprints (8,374,968 footprints) #
################################################
entries = 0
footprints = {} # chr -> [footprint, ... ] - ordered by genomic region
inFile = open('combined.fps','r')
# Read in line by line and build footprints dictionary
while 1:
    inLine = inFile.readline()
    if not inLine:
        break
    splitUp = inLine.strip().split()
    if not splitUp[0] in footprints:
        footprints[splitUp[0]] = []
    footprints[splitUp[0]].append([splitUp[0], int(splitUp[1]), int(splitUp[2])])
    entries += 1

inFile.close()
print 'Footprints =',entries

"""
## Merge overlapping sequences 6bp
# For each chromosome
entMerg = 0
fpMerged ={}
for i in footprints:
    fpMerged[i] = []
    prev = ''
    fpLen = len(footprints[i])
    cur = 0
    while 1:
        cur += 1
        if not cur<fpLen:
            break
        if not prev=='' and prev[2]<=(footprints[i][cur][1]-6):
            prev[2] = footprints[i][cur][2]
        else:
            if not prev=='':
                fpMerged[i].append(prev)
                entMerg += 1
            prev = footprints[i][cur]

print 'Footprints =',entries,'; Merged = ',entMerg,'; Merge Length = 6'

##########################
# 2. Parse out sequences #
##########################
print '  Extracting the sequence data...'
# Unzip sequences for extraction
#tar = tarfile.open('sequences/fasta/chromFa.tar.gz')
#tar.extractall(path='sequences/fasta')
#tar.close()

# 6. Extract the sequences
for chrom in fpMerged:
    footprintSeqFile = open('footprintSequences_'+chrom+'.fasta','w')
    if os.path.exists('sequences/fasta/'+str(chrom)+'.fa'):
        chrSeqFile = open('sequences/fasta/'+str(chrom)+'.fa','r')
    elif os.path.exists('sequences/fasta/'+str(chrom).lstrip('chr').replace('_random','')+'/'+str(chrom)+'.fa'):
            chrSeqFile = open('sequences/fasta/'+str(chrom).lstrip('chr').replace('_random','')+'/'+str(chrom)+'.fa','r')
    else:
        print 'FATAL ERROR!!!! Arghhh',chrom,'(',str(chrom).lstrip('chr').replace('_random',''),')does not have a seqeunce file!'
        break
    chrSeqFile.readline() # Get rid of header
    chrSeq = [x.strip().upper() for x in chrSeqFile.readlines()]
    chrSeq = ''.join(chrSeq)
    print '  ',chrom,' (',len(chrSeq),')...'
    for fp1 in fpMerged[chrom]:
        footprintSeq = chrSeq[(fp1[1]-1):(fp1[2]-1)]
        footprintSeqFile.write('>'+str(fp1[0])+':'+str(fp1[1])+'-'+str(fp1[2])+'\n'+str(footprintSeq)+'\n')
    footprintSeqFile.close()
"""

## Merge overlapping sequences 6bp
# For each chromosome
entMerg = 0
fpMergedMask ={}
outFile = open('merges.tsv','w')
for i in footprints:
    fpMergedMask[i] = []
    prev = ''
    fpLen = len(footprints[i])
    merged = []
    for cur in range(len(footprints[i])):
        if not prev=='' and prev[2]>=footprints[i][cur][1]:
            #print prev, footprints[i][cur]
            if prev[2]<footprints[i][cur][2]:
                prev[2] = footprints[i][cur][2]
            merged.append(footprints[i][cur])
        else:
            if not prev=='':
                fpMergedMask[i].append(prev)
                entMerg += 1
                if len(merged)>1:
                    outFile.write('\n'+'\t'.join([str(j) for j in merged]))
            merged = [footprints[i][cur]]
            prev = footprints[i][cur][:]

if not prev=='':
    fpMergedMask[i].append(prev)
    entMerg += 1
    if len(merged)>1:
        outFile.write('\n'+'\t'.join([str(j) for j in merged]))
outFile.close()

print 'Footprints =',entries,'; Merged = ',entMerg,'; Merge Length = 1'


#######################################
# 3. Mask out non-footprint sequences #
#######################################
print '  Maksing out sequence data...'

# Unzip sequences for extraction
#tar = tarfile.open('sequences/fasta/chromFa.tar.gz')
#tar.extractall(path='sequences/fasta')
#tar.close()

# Make masked chomosomal fasta sequences
lenFile = open('seqLens.csv','w')
lenFile.write('chromosome,len,fpLen,%')
for chrom in footprints:
    # Read in sequences
    if os.path.exists('footprintSeqs/'+str(chrom)+'.fa'):
        chrSeqFile = open('footprintSeqs/'+str(chrom)+'.fa','r')
    elif os.path.exists('footprintSeqs/'+str(chrom).lstrip('chr').replace('_random','')+'/'+str(chrom)+'.fa'):
            chrSeqFile = open('footprintSeqs/'+str(chrom).lstrip('chr').replace('_random','')+'/'+str(chrom)+'.fa','r')
    else:
        print 'FATAL ERROR!!!! Arghhh',chrom,'(',str(chrom).lstrip('chr').replace('_random',''),')does not have a seqeunce file!'
        break
    header = chrSeqFile.readline() # Get rid of header
    chrSeq = [x.strip().upper() for x in chrSeqFile.readlines()]
    chrSeq = ''.join(chrSeq)
    sumLen = sum([i[2]-i[1] for i in fpMergedMask[chrom]])
    print '  ',chrom,' (',len(chrSeq),'); merged:  '+str(sumLen)+'...'
    lenFile.write('\n'+chrom+','+str(len(chrSeq))+','+str(sumLen)+','+str(float(sumLen)/float(len(chrSeq))))
    # Make N'd version of chrSeq
    maskedChrSeq = 'N'*(fpMergedMask[chrom][0][1])+chrSeq[(fpMergedMask[chrom][0][1]-1):(fpMergedMask[chrom][0][2]-1)]
    maskedChrSeq += ''.join(['N'*((fpMergedMask[chrom][i][1])-(fpMergedMask[chrom][i-1][2]))+chrSeq[(fpMergedMask[chrom][i][1]-1):(fpMergedMask[chrom][i][2]-1)] for i in range(1,len(fpMergedMask[chrom]))])
    maskedChrSeq += 'N'*(len(chrSeq)-(fpMergedMask[chrom][len(fpMergedMask[chrom])-1][2]))
    print '  After masking length:',len(maskedChrSeq)
    print '  Writing out sequence...'
    outFile = open('footprintSeqs2/'+str(chrom)+'.fa','w')
    outFile.write(header.strip())
    divs = len(maskedChrSeq)/100
    lo1 = len(maskedChrSeq)-(divs*100)
    [outFile.write('\n'+maskedChrSeq[((i*100)-100):(i*100)]) for i in range(1,divs+1)]
    #[outFile.write('\n'+maskedChrSeq[((i*100)-100):(i*100)]) for i in range(1,divs)]
    #yli: I changed "range(1,divs)" to "range(1,divs+1)". 
    #yli: In the 2016 version data, masked reference was 100bp shorter than the original reference, the second last line was missing. 
    #yli: I think this bug actually would not affect the final result because the last couple of lines for a chr sequences are usually N bases.
    if lo1>0:
        outFile.write('\n'+maskedChrSeq[(divs*100):])
    outFile.close()
lenFile.close()



