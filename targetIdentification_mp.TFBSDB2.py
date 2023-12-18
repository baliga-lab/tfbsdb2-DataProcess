#################################################################
# Program: targetIdentification_mp.v2022.py                     #
# Version: v2022                                                #
# Author: Christopher L Plaisier, PhD                           #
# Modified by: Yaqiao Li                                        #
# Baliga Lab, ISB                                               #
# Institute for Systems Biology                                 #
#################################################################

import gzip, os, sys, re, os, math, shutil, errno
from copy import deepcopy
from subprocess import *
from ftplib import FTP
import tarfile
from random import sample
from multiprocessing import Pool, cpu_count, Manager
import time
from sys import stdout
import _pickle as cPickle
import argparse

DESCRIPTION = """targetIdentification_mp.v2022.py - get statistic information of optimal promoter region"""

# Creat directories
def create_dirs(dir_list):
    for dir_name in dir_list:
        try:
            os.mkdir(dir_name)
        except OSError as exc:
            if exc.errno != errno.EEXIST:
                raise
            pass


# Return start as the key for sorting
def startKey1(hit):
    return hit[0][0]

# Compound sorting
def sortMergedSets(hits):
        return sorted(hits, key=startKey1)

# Return start as the key for sorting
def startKey2(hit):
    return hit['start']

# Compound sorting
def sortHits(hits):
        return sorted(hits, key=startKey2)

# Uniquify
def uniquify(seq):
    # Not order preserving
    return {}.fromkeys(seq).keys()

# Complement
def complement(seq):
    complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    complseq = [complement[base] for base in seq]
    return complseq

# Reverse complement
def reverseComplement(seq):
    seq = list(seq)
    seq.reverse()
    return ''.join(complement(seq))

# Function to retreive boundaries for 3pUTR
def getPromoter(geneCoords,upstream):
    if geneCoords['strand']=='+':
        return [(geneCoords['txStart'] - upstream[0]), (geneCoords['txStart'] - upstream[1])]
    elif geneCoords['strand']=='-':
        return [(geneCoords['txEnd'] + upstream[1]), (geneCoords['txEnd'] + upstream[0])]

# Merge overlaps and give back sequences [[5pStart,5pEnd], [[3pStart,3pEnd], ... ]]
# !!! - Assumes that the mergeDem entries come from the same chromosome
def mergeSeqs(mergeDem,upstream): #,min5pUTR,min3pUTR):
    orig = mergeDem[0]
    strand = orig['strand']
    # Grab the starting promoter, 5' UTR, coding and 3' UTR sequences
    promoter = getPromoter(orig,upstream)
    """p5utr = get5pUTR(orig,min5pUTR)
    cds = getCDS(orig)
    p3utr = get3pUTR(orig,min3pUTR)"""
    # Now iterate through the rest and merge
    for i in range(1,len(mergeDem)):
        mergeMe = mergeDem[i]
        # Merge promoter sequences
        promoterM = getPromoter(mergeMe,upstream)
        if promoterM[0] < promoter[0] <= promoterM[1] <= promoter[1]:
            promoter[0] = promoterM[0]
        elif promoter[0] <= promoterM[0] <= promoter[1] < promoterM[1]:
            promoter[1] = promoterM[1]
        elif not promoter==promoterM:
            # COMPROMISE HERE: Then will take the one for the longest transcript
            if strand=='-':
                if promoterM[0] > promoter[1]:
                    promoter = promoterM
            elif strand=='+':
                if promoterM[1] < promoter[0]:
                    promoter = promoterM
        """# Merge 5' UTR
        p5utrM = get5pUTR(mergeMe,min5pUTR)
        #print p5utrM,'; ',p5utr
        if not ((p5utrM[0][0]==p5utr[0][0]) and (p5utrM[len(p5utrM)-1][1]==p5utr[len(p5utr)-1][1])):
            #print "Need merging!"
            #print p5utrM,p5utr
            p5utr = mergeICS(p5utr,p5utrM)
            #print "Merged: ",p3utr
        # Merge coding sequences
        cdsM = getCDS(mergeMe)
        #print cdsM,'; ',cds
        if not ((cdsM[0][0]==cds[0][0]) and (cdsM[len(cdsM)-1][1]==cds[len(cds)-1][1])):
            #print "Need merging!"
            #print cdsM,cds
            cds = mergeICS(cds,cdsM)
            #print "Merged: ",cds
        # Merge 3' UTR
        p3utrM = get3pUTR(mergeMe,min3pUTR)
        #print p3utrM,'; ',p3utr
        if not ((p3utrM[0][0]==p3utr[0][0]) and (p3utrM[len(p3utrM)-1][1]==p3utr[len(p3utr)-1][1])):
            #print "Need merging!"
            #print p3utrM,p3utr
            p3utr = mergeICS(p3utr,p3utrM)
            #print "Merged: ",p3utr"""
    return [promoter] #, p5utr, cds, p3utr]

def overlap(a,b):
    """
    Possible overlaps between two genomic features
    <-a->      : b.start <= a.stop <= b.stop (4) *
        <-b->  : a.start <= b.start <= a.stop (3)
    or
    <-b->      : a.start <= b.stop <= a.stop (3) *
        <-a->  : b.start <= a.start <= b.stop (4)
    or
    <--a-->    : a.start <= b.start <= b.stop <= a.stop
     <-b->     :
    or
     <-a->     : b.start <= a.start <= a.stop <= b.stop
    <--b-->    :
    """
    if (b['start'] <= a['stop'] <= b['stop']) or (a['start'] <= b['stop'] <= a['stop']):
        return True
    return False


#------------------------------------------------------------------------------------ Main -----------#




parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                    description=DESCRIPTION)
parser.add_argument('--outputDir', help='Output dir of statistic results', type=str)
parser.add_argument('--promoterStart', help='Start position of promoter from TSS', nargs='?', const=0, type=int)
parser.add_argument('--promoterEnd', help='End position of promoter from TSS', nargs='?', const=0, type=int)
parser.add_argument('--proximalStart', help='Start position of proximal site from TSS', nargs='?', const=0, type=int)
parser.add_argument('--proximalEnd', help='End position of proximal site from TSS', nargs='?', const=0, type=int)
parser.add_argument('--distalStart', help='Start position of distal site from TSS', nargs='?', const=0, type=int)
parser.add_argument('--distalEnd', help='End position of distal site from TSS', nargs='?', const=0, type=int)

args = parser.parse_args()

output_dir     = args.outputDir
promoter_start = args.promoterStart
promoter_end   = args.promoterEnd
proximal_start = args.proximalStart
proximal_end   = args.proximalEnd
distal_start   = args.distalStart
distal_end     = args.distalEnd
 
promoterSeq = [promoter_start, promoter_end]
proximalSeq = [proximal_start, proximal_end]
distalSeq   = [distal_start, distal_end]

    #promoterSeq = [5000,-5000]
    #proximalSeq = [2500,-500]
    #distalSeq = [5000,2500]

dir_list = []
dir_list.append(output_dir)
geneHits_dir = output_dir + '/geneHitsDB'
dir_list.append(geneHits_dir)

create_dirs(dir_list)









#------------------------------------------------------------------------------------ Prepare promoter region-----------#

# Download gene identifier conversion table from NCBI if not already done
print('Downloading conversion table for Entrez IDs to RefSeq IDs...')
if not os.path.exists('gene2refseq.gz'):
    ftp1 = FTP('ftp.ncbi.nih.gov')
    ftp1.login()
    ftp1.cwd('/gene/DATA/')
    outFile = open('gene2refseq.gz','wb')
    ftp1.retrbinary('RETR gene2refseq.gz',outFile.write)
    outFile.close()
    ftp1.quit()

# Setup so can run on other species if wanted (specifically mouse)
org = 'Homo_sapiens'
orgDict = { 'Homo_sapiens': { 'orgId':9606 } }

# Download the genomic coordinates
print('Starting on ' + str(org) + '...')
print('  Downloading genomic coordinates...')
# Download genome information for organism if not already done
if not os.path.exists(str(org)+'/refGene.txt.gz'):
    # Download genome information from UCSC
    ftp1 = FTP('hgdownload.cse.ucsc.edu')
    ftp1.login()
    # Get gene coordinate data
    ftp1.cwd('/goldenPath/currentGenomes/'+str(org)+'/database/')
    outFile = open(str(org)+'/refGene.txt.gz','wb')
    ftp1.retrbinary('RETR refGene.txt.gz',outFile.write)
    outFile.close()
    ftp1.quit()

print('  Building dictionaries...')
# 1. Read in refSeqCoords
inFile = gzip.open(str(org)+'/refGene.txt.gz','r')
refseqCoords = {}
chrs = []
while 1:
    line = inFile.readline()
    if not line:
        break
    line = str(line, encoding="utf-8")
    splitUp = line.strip().split('\t')
    if len(splitUp)>=13 and splitUp[13]=='cmpl':
        if not splitUp[1] in refseqCoords:
            if not splitUp[2] in chrs:
                chrs.append(splitUp[2])
            refseqCoords[splitUp[1]] = {'chr':splitUp[2], 'strand':splitUp[3], 'txStart':int(splitUp[4]), 'txEnd':int(splitUp[5]), 'cdsStart':int(splitUp[6]), 'cdsEnd':int(splitUp[7]), 'exonCount':int(splitUp[8]), 'exonStarts':splitUp[9], 'exonEnds':splitUp[10], 'geneName':splitUp[12], 'exonFrames':[int(x) for x in splitUp[15].split(',') if x]}
    elif not len(splitUp)>=13:
        if not splitUp[0] in refseqCoords:
            if not splitUp[1] in chrs:
                chrs.append(splitUp[1])
            # Build the exonFrames determine which 
            refseqCoords[splitUp[0]] = {'chr':splitUp[1], 'strand':splitUp[2], 'txStart':int(splitUp[3]), 'txEnd':int(splitUp[4]), 'cdsStart':int(splitUp[5]), 'cdsEnd':int(splitUp[6]), 'exonCount':int(splitUp[7]), 'exonStarts':splitUp[8], 'exonEnds':splitUp[9]}
inFile.close()

print('  Building Entrez ID to RefSeq ID dictionary...')
# 2. Make a dictionary of EntrezIDs to RefSeqIds
if not os.path.exists('entrezID2refSeq.pkl'):

    inFile = gzip.open('gene2refseq.gz','r')
    inFile.readline() # skip header
    entrezId2refSeq = {}
    while 1:
        line = inFile.readline()
        if not line:
            break
        # Only add those that have the correct NCBI organism ID
        line = str(line, encoding="utf-8")
        splitUp = line.strip().split('\t')
        if int(splitUp[0])==orgDict[org]['orgId']:
            #print splitUp[3],splitUp[3].split('.')[0]
            # Check that the nucleotide ID is not a '-' and that it has genomic coordiantes assocaited with it
            if not splitUp[3]=='-' and splitUp[3].split('.')[0] in refseqCoords:
                if not int(splitUp[1]) in entrezId2refSeq:
                    entrezId2refSeq[int(splitUp[1])] = [splitUp[3].split('.')[0]]
                else:
                    entrezId2refSeq[int(splitUp[1])].append(splitUp[3].split('.')[0])
    inFile.close()
    pklFile = open('entrezID2refSeq.pkl','wb')
    cPickle.dump(entrezId2refSeq,pklFile)
    pklFile.close()
else:
    pklFile = open('entrezID2refSeq.pkl','rb')
    entrezId2refSeq = cPickle.load(pklFile)
    pklFile.close()
print(' ' + "entrezIdnum" + "refseqCoordsnum")
print(' ' + str(len(entrezId2refSeq)) + str(len(refseqCoords)))


print('  Now collapsing and merging RefSeq IDs into Entrez IDs...')
# 3. Merege multiple refseq IDs corresponding to a single entrezID
# Now merge the data
#chrs = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY'] # This needs to come from the chromsome fasta masked files
mergedSet = {}
for chr in chrs:
    mergedSet[chr] = {}
baddies = []
chrNoMatch = 0
gotGenes = []
for entrezId in entrezId2refSeq:
    chr = refseqCoords[entrezId2refSeq[entrezId][0]]['chr']
    strand = refseqCoords[entrezId2refSeq[entrezId][0]]['strand']
    if chr in mergedSet:
        if len(entrezId2refSeq[entrezId])>1:
            # There are duplicates so build a list of their refseqCoords and merge them
            mergeDem = []
            for refseq in entrezId2refSeq[entrezId]:
                tot = 0
                negOne = 0
                goodOrBad = 1
                if 'exonFrames' in refseqCoords[refseq]:
                    for i in refseqCoords[refseq]['exonFrames']:
                        if i == -1:
                            negOne += 1
                        tot += 1
                    if len(refseqCoords[refseq]['exonFrames'])>=5:
                        goodOrBad = 1-float(negOne)/float(tot)
                #print goodOrBad
                if chr=='' and goodOrBad>=0.5:
                    chr = refseqCoords[refseq]['chr']
                    mergeDem.append(refseqCoords[refseq])
                elif chr==refseqCoords[refseq]['chr'] and refseqCoords[refseq]['strand'] and goodOrBad>=0.5:
                    mergeDem.append(refseqCoords[refseq])
                elif goodOrBad<0.5:
                    #print 'Baddie taken out: ',refseq
                    baddies.append(refseq)
                else:
                    #print 'Uh oh! Chr and Strand don\'t match! EntrezID = ',entrezId,'; refseqID = ',refseq
                    chrNoMatch += 1
            if len(mergeDem)>1:
                mergedSet[chr][entrezId] = mergeSeqs(mergeDem,promoterSeq) + [strand] #,min5pUTR,min3pUTR) + [strand]
                gotGenes.append(entrezId)
                #print entrez, len(entrez2refseq[entrez]), lenICS(mergedSet[chr][entrez][1])
            else:
                promoter = getPromoter(refseqCoords[(entrezId2refSeq[entrezId])[0]],promoterSeq)
                """p5utr = get5pUTR(refseqCoords[(entrezId2refSeq[entrezId])[0]],min5pUTR)
                cds = getCDS(refseqCoords[(entrezId2refSeq[entrezId])[0]])
                p3utr = get3pUTR(refseqCoords[(entrezId2refSeq[entrezId])[0]],min3pUTR)"""
                mergedSet[chr][entrezId] = [promoter,strand]
                gotGenes.append(entrezId)
        else:
            promoter = getPromoter(refseqCoords[(entrezId2refSeq[entrezId])[0]],promoterSeq)
            """p5utr = get5pUTR(refseqCoords[(entrezId2refSeq[entrezId])[0]],min5pUTR)
            cds = getCDS(refseqCoords[(entrezId2refSeq[entrezId])[0]])
            p3utr = get3pUTR(refseqCoords[(entrezId2refSeq[entrezId])[0]],min3pUTR)"""
            mergedSet[chr][entrezId] = [promoter,strand]
            gotGenes.append(entrezId)
badFile = open(output_dir + '/promoter_' + str(promoter_start) + '_' + str(promoter_end) + '_baddies.txt','w')
badFile.write('\n'.join(baddies))
badFile.close()
del refseqCoords
del entrezId2refSeq

# 4. Calculate total regions searched for proximal and distal
proximalSeqLen = 0
distalSeqLen = 0
#chrSeqs = {}
for chrom in ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']:
        # Read in sequences
        if os.path.exists('footprintSeqs2/'+str(chrom)+'.fa'):
            chrSeqFile = open('footprintSeqs2/'+str(chrom)+'.fa','r')
        chrSeqFile.readline() # Get rid of header
        chrSeq = [x.strip().upper() for x in chrSeqFile.readlines()]
        chrSeq = ''.join(chrSeq)
        #chrSeqs[chrom] = chrSeq
        stdout.write(chrom+' ')
        stdout.flush()
        # For each gene in mergedSet determine if binding site in promoter
        for gene in mergedSet[chrom]:
            if mergedSet[chrom][gene][1]=='+':
                prox1 = [((mergedSet[chrom][gene][0][0]-1)-proximalSeq[0]), ((mergedSet[chrom][gene][0][0]-1)-proximalSeq[1])]
                dist1 = [((mergedSet[chrom][gene][0][0]-1)-distalSeq[0]), ((mergedSet[chrom][gene][0][0]-1)-distalSeq[1])]
            else:
                prox1 = [((mergedSet[chrom][gene][0][0]-1)+proximalSeq[1]), ((mergedSet[chrom][gene][0][0]-1)+proximalSeq[0])]
                dist1 = [((mergedSet[chrom][gene][0][0]-1)-distalSeq[1]), ((mergedSet[chrom][gene][0][0]-1)-distalSeq[0])]
            tmp1 = chrSeq[prox1[0]:prox1[1]]
            proximalSeqLen += len(tmp1)-tmp1.count('N')
            tmp1 = chrSeq[dist1[0]:dist1[1]]
            distalSeqLen += len(tmp1)-tmp1.count('N')
stdout.write('\n')
stdout.flush()
print("proximalSeqLen" + "\t" + "distalSeqLen")
print(str(proximalSeqLen) + "\t" + str(distalSeqLen))
del chrSeq

# 6. Determine the footprints inside the defined promoter region
"""
print '  Find footprints in promoters...'
if not os.path.exists('geneFootprints.pkl'):
    # 5. Read in footprints (8,374,968 footprints) #
    print '  Loading footprints...'
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
    print '  Footprints =',entries

    # Create overlap
    geneFootprints = {}
    stdout.write('  ')
    for chrom in ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']:
        # Read in sequences
        stdout.write(chrom+' ')
        stdout.flush()
        # Create segments for searching in the footprints
        numSeg2 = 500
        lmh2 = len(footprints)
        seg2 = lmh2/numSeg2
        segments2 = []
        for i in range(1,numSeg2):
            segments2.append([int(footprints[chrom][i*seg2][1]), i*seg2])
        # Prep storage for footprints
        geneFootprints[chrom] = {}
        # For each gene in mergedSet determine if binding site in promoter
        for gene in mergedSet[chrom]:
            stdout.write('.')
            stdout.flush()
            curStart2 = 0
            for i in range(len(segments2)):
                if mergedSet[chrom][gene][0][1] <= segments2[i][0]:
                    break
                curStart2 = segments2[i][1]
            for fp1 in footprints[chrom][curStart2:]:
                if mergedSet[chrom][gene][0][1] < int(fp1[2]):
                    break
                if overlap({'start':mergedSet[chrom][gene][0][0],'stop':mergedSet[chrom][gene][0][1]}, {'start':fp1[1], 'stop':fp1[2]}):
                    if not gene in geneFootprints[chrom]:
                        geneFootprints[chrom][gene] = []
                    geneFootprints[chrom][gene].append(fp1)
    print '  Done.\n'
    pklFile = open('geneFootprints.pkl','wb')
    cPickle.dump(geneFootprints, pklFile)
else:
    pklFile = open('geneFootprints.pkl','rb')
    geneFootprints = cPickle.load(pklFile)
pklFile.close()
"""

# 7. Get the names of the motifs to screen through
tfs = uniquify(['_'.join(('.'.join(i.split('.')[:-1])).split('_')[:-1]) for i in os.listdir('./motifHits2') if i.count('pkl')==1])
print('TFs:' + str(len(tfs)))

# For multicore processing 
def overlapFIMO(chrom):
   # Read in sequences
    pklFile = open('./motifHits2/'+stats['tf']+'_'+chrom+'.pkl','rb')
    motifHits = cPickle.load(pklFile)
    motifHits = sortHits(motifHits)
    stats['hits'] += len(motifHits)
    pklFile.close()
    if len(motifHits)>0:
        # Create segments for searching in the motifHits
        numSeg = 50
        lmh = len(motifHits)
        seg = lmh//numSeg
        segments = []
        for i in range(1,numSeg):
            segments.append([int(motifHits[i*seg]['start']), i*seg])
        # For each gene in mergedSet determine if binding site in promoter
        for gene in geneFootprints[chrom]:
            locations = []
            strands = []
            for fp1 in geneFootprints[chrom][gene]:
                curStart = 0
                for i in range(len(segments)):
                    if mergedSet[chrom][gene][0][1] <= segments[i][0]:
                        break
                    curStart = segments[i][1]
                for hit in motifHits[curStart:]:
                    if fp1[2] < int(hit['start']):
                        break
                    if overlap({'start':fp1[1], 'stop':fp1[2]},{'start':int(hit['start']), 'stop':int(hit['stop'])}):
                        locations.append([hit['start'],hit['stop']])
                        strands.append(hit['strand'])
                        dist = 0
                        if mergedSet[chrom][gene][1]=='+':
                            if hit['strand']=='+':
                                stats['plus'] += 1
                            else:
                                stats['minus'] += 1
                            dist = int(mergedSet[chrom][gene][0][1])-int(hit['stop'])
                        else:
                            if hit['strand']=='-':
                                stats['minus'] += 1
                            else:
                                stats['plus'] += 1
                            dist = int(hit['start'])-int(mergedSet[chrom][gene][0][0])
                        if dist<=stats['proximalSeq'][0] and dist>stats['proximalSeq'][1]:
                            stats['proximal'] += 1
                        elif dist<=stats['distalSeq'][0] and dist>stats['distalSeq'][1]:
                            stats['distal'] += 1
            if len(locations)>0:
                writeMe.append('\n'+str(gene)+','+str(mergedSet[chrom][gene][0])+','+chrom+','+str(len(locations))+','+';'.join(['-'.join([str(j) for j in i]) for i in locations])+','+';'.join([str(i) for i in strands]))
    stdout.write(chrom+' ')
    stdout.flush()

# For multicore processing 
def overlapFIMOFaster(chrom):
    # Temporary variables
    hits = 0
    proximal = 0
    distal = 0
    plus = 0
    minus = 0
    proximalSeq = stats['proximalSeq']
    distalSeq = stats['distalSeq']
    writeMeTmp = []
    # Read in sequences
    chrSeqFile = open('footprintSeqs2/'+str(chrom)+'.fa','r')
    chrSeqFile.readline() # Get rid of header
    chrSeq = [x.strip().upper() for x in chrSeqFile.readlines()]
    chrSeq = ''.join(chrSeq)
    chrSeqFile.close()
    #global chrSeqs
    # Load up the TFBS predicted by FIMO
    pklFile = open('./motifHits2/'+stats['tf']+'_'+chrom+'.pkl','rb')
    motifHits = cPickle.load(pklFile)
    motifHits = sortHits(motifHits)
    hits += len(motifHits)
    pklFile.close()
    # If the motifHits is empty then skip this and return empty
    if len(motifHits)>0:
        # Create segments for searching in the motifHits
        numSeg = 50
        lmh = len(motifHits)
        seg = lmh//numSeg
        segments = []
        for i in range(1,numSeg):
            segments.append([int(motifHits[i*seg]['start']), i*seg])
        # For each gene in mergedSet determine if binding site in promoter
        for gene in mergedSet[chrom]:
            locations = []
            strands = []
            pValues = []
            matchSequences = []
            curStart = 0
            for i in range(len(segments)):
                if mergedSet[chrom][gene][0][1] <= segments[i][0]:
                    break
                curStart = segments[i][1]
            for hit in motifHits[curStart:]:
                if mergedSet[chrom][gene][0][1] < int(hit['start']):
                    break
                if overlap({'start':mergedSet[chrom][gene][0][0], 'stop':mergedSet[chrom][gene][0][1]},{'start':int(hit['start']), 'stop':int(hit['stop'])}):
                    locSeq = chrSeq[(int(hit['start'])-1):(int(hit['stop'])-1)]
                    if (len(locSeq)-locSeq.count('N')) > 0:
                        locations.append([hit['start'],hit['stop']])
                        strands.append(hit['strand'])
                        pValues.append(hit['p.value'])
                        matchSequences.append(hit['match.sequence'])
                        dist = 0
                        if mergedSet[chrom][gene][1]=='+':
                            if hit['strand']=='+':
                                plus += 1
                            else:
                                minus += 1
                            dist = int(mergedSet[chrom][gene][0][1])-int(hit['stop'])
                        else:
                            if hit['strand']=='-':
                                minus += 1
                            else:
                                plus += 1
                            dist = int(hit['start'])-int(mergedSet[chrom][gene][0][0])
                        if dist<=proximalSeq[0] and dist>proximalSeq[1]:
                            proximal += 1
                        elif dist<=distalSeq[0] and dist>distalSeq[1]:
                            distal += 1
            if len(locations)>0:
                writeMeTmp.append(str(gene)+','+str('-'.join([str(i) for i in mergedSet[chrom][gene][0]]))+','+chrom+','+str(len(locations))+','+';'.join(['-'.join([str(j) for j in i]) for i in locations])+','+';'.join([str(i) for i in strands])+','+';'.join([str(i) for i in pValues])+','+';'.join([str(i) for i in matchSequences]))
    stdout.write(chrom+' ')
    stdout.flush()
    stats['hits'] += hits
    stats['proximal'] += proximal
    stats['distal'] += distal
    stats['plus'] += plus
    stats['minus'] += minus
    writeMe[chrom] = writeMeTmp
    
# 8. Load a FIMO output for a TF and determine if binding site is in promoter
# Retain: Entrez ID,Instances,Locations,Strandednesses
#tfs = ['motifHits/DistalBias/motifHits_Tcf7l2.1']
#tfs = ['motifHits/DistalBias/motifHits_Tcf7l2.1','motifHits/DistalBias/motifHits_SOX2_HMG_full_dimeric_17_1','motifHits/DistalBias/motifHits_V_OCT4_01_M01125']
#tfs = ['motifHits/ProximalBias/motifHits_ERG_ETS_full_dimeric_14_1','motifHits/ProximalBias/motifHits_V_NFYC_Q5_M02107','motifHits/ProximalBias/motifHits_E2F1_E2F_DBD_dimeric_12_1']
#tfs = ['motifHits/DistalBias/motifHits_Tcf7l2.1','motifHits/DistalBias/motifHits_Tcf7l2.2','motifHits/DistalBias/motifHits_SOX2_HMG_full_dimeric_17_1','motifHits/DistalBias/motifHits_V_OCT4_01_M01125','motifHits/ProximalBias/motifHits_ERG_ETS_full_dimeric_14_1','motifHits/ProximalBias/motifHits_V_NFYC_Q5_M02107','motifHits/ProximalBias/motifHits_E2F1_E2F_DBD_dimeric_12_1']
print('  Starting overlap analysis...')
cpus = cpu_count()
print('There are' + str(cpus) + 'CPUs avialable.')
for tf in tfs:
    mgr = Manager()
    stats = mgr.dict()
    stats['tf'] = tf
    stats['hits'] = 0
    stats['proximal'] = 0
    stats['distal'] = 0
    stats['plus'] = 0
    stats['minus'] = 0
    stats['proximalSeq'] = proximalSeq
    stats['distalSeq'] = distalSeq
    writeMe = mgr.dict()
    stdout.write('  ')
    #overlapFIMO('chr1')
    chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
    geneHitscsv = geneHits_dir + '/fullGenome_'+tf+'.csv'
    if not os.path.exists(geneHitscsv):
        print(tf + " not exists\n")
        if not os.path.exists(geneHits_dir):# add by yli
            os.mkdir(geneHits_dir)# add by yli
        #if os.path.getsize('./motifHits2/'+tf+'_chr1.pkl')<=83886080:
        pool = Pool(processes=cpus)
        pool.map(overlapFIMOFaster,chroms)
        pool.close()
        pool.join()
        #else:
        #    print 'Files big going single core...'
        #    for chrom in chroms:
        #        print chrom
        #        overlapFIMOFaster(chrom)
        outFile = open(geneHitscsv, 'w')
        outFile.write('Entrez ID,Chr,Instances,Locations,Strands')
        outFile.write('\n'+'\n'.join(['\n'.join(writeMe[chrom]) for chrom in chroms]))
        outFile.close()
        if not stats['distal']==0 and not stats['proximal']==0:
            percDistal = ((float(stats['distal'])/float(distalSeqLen))/((float(stats['proximal'])/float(proximalSeqLen))+(float(stats['distal'])/float(distalSeqLen))))
        else:
            percDistal = 'NA'
        if not stats['plus']==0 and not stats['minus']==0:
            percPlus = float(stats['plus'])/float(stats['plus']+stats['minus'])
        else:
            percPlus = 'NA'
        stdout.write('\n')
        stdout.flush()
        print('TF: ' + tf + '; Hits: ' + str(stats['hits']) + '; Proximal: ' + str(stats['proximal']) + '; Distal: ' + str(stats['distal']) + '; % Distal: ' + str(percDistal) + '; Minus: ' + str(stats['minus']) + '; Plus: ' + str(stats['plus']) + '; % Plus: ' + str(percPlus))
        out_txt = output_dir + '/tf_promoter_' + str(promoter_start) + '-' + str(promoter_end) + '_proximal_' + str(proximal_start) + '-' + str(proximal_end) + '_distal_' + str(distal_start) + '-' + str(distal_end) + '.txt'
        if not os.path.exists(out_txt):
            outFile2 = open(out_txt,'w')
            outFile2.write('TF,Hits,Proximal,ProximalSeqLen,Distal,DistalSeqLen,% Distal,Minus,Plus,% Plus')
        else:
            outFile2 = open(out_txt,'a')
        outFile2.write('\n'+tf+','+str(stats['hits'])+','+str(stats['proximal'])+','+str(proximalSeqLen)+','+str(stats['distal'])+','+str(distalSeqLen)+','+str(percDistal)+','+str(stats['minus'])+','+str(stats['plus'])+','+str(percPlus))
        outFile2.close()
    del mgr
print('  Done.')
