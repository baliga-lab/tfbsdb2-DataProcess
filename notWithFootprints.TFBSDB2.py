#################################################################
# @Program: targetIdentification.py                             #
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
# Copyrighted by Chris Plaisier 2/6/2013                        #
#################################################################

###############
### IMPORTS ###
###############
#from pssm import pssm
import gzip, os, sys, re, os, math, shutil
from copy import deepcopy
from subprocess import *
from ftplib import FTP
import tarfile
from random import sample
from multiprocessing import Pool, cpu_count, Manager
import time
from sys import stdout
import _pickle as cPickle
import statistics

#########################
## Important Variables ##
#########################

promoter_start = int(sys.argv[1])
promoter_end = int(sys.argv[2])
output_dir = sys.argv[3]

if not os.path.exists(output_dir):
    os.mkdir(output_dir)

promoterSeq = [promoter_start, promoter_end]
proximalSeq = [2500,-500]
distalSeq = [5000,2500]

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
print('Starting on '+str(org)+'...')
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
        line = str(line, encoding="utf-8")
        # Only add those that have the correct NCBI organism ID
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
print('---------------------------------------------------')
print(' '+str(len(entrezId2refSeq))+str(len(refseqCoords)))


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
badFile = open(str(org)+'/baddies.txt','w')
badFile.write('\n'.join(baddies))
badFile.close()
del refseqCoords
del entrezId2refSeq

# 7. Get the names of the motifs to screen through
tfs = uniquify(['_'.join(('.'.join(i.split('.')[:-1])).split('_')[:-1]) for i in os.listdir('./motifHits2') if i.count('pkl')==1])
print('TFs: ' + str(len(tfs)))

# For multicore processing 
def overlapFIMOFaster(chrom):
    # Temporary variables
    writeMeTmp = []
    # Load up the TFBS predicted by FIMO
    pklFile = open('./motifHits2/'+stats['tf']+'_'+chrom+'.pkl','rb')
    motifHits = cPickle.load(pklFile)
    motifHits = sortHits(motifHits)
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
            pvalues = []
            curStart = 0
            for i in range(len(segments)):
                if mergedSet[chrom][gene][0][1] <= segments[i][0]:
                    break
                curStart = segments[i][1]
            for hit in motifHits[curStart:]:
                if mergedSet[chrom][gene][0][1] < int(hit['start']):
                    break
                if overlap({'start':mergedSet[chrom][gene][0][0], 'stop':mergedSet[chrom][gene][0][1]},{'start':int(hit['start']), 'stop':int(hit['stop'])}):
                    locations.append([hit['start'],hit['stop']])
                    strands.append(hit['strand'])
                    pvalues.append(hit['p.value'])#--add by yli
            if len(locations)>0:
                #--anno by yli writeMeTmp.append(str(gene)+','+str(mergedSet[chrom][gene][0])+','+chrom+','+str(len(locations))+','+';'.join(['-'.join([str(j) for j in i]) for i in locations])+','+';'.join([str(i) for i in strands]))
                writeMeTmp.append(str(gene)+','+str(mergedSet[chrom][gene][0])+','+chrom+','+str(len(locations))+','+';'.join(['-'.join([str(j) for j in i]) for i in locations])+','+';'.join([str(i) for i in strands])+','+';'.join([str(i) for i in pvalues]))#--modified by yli
    stdout.write(chrom+' ')
    stdout.flush()
    writeMe[chrom] = writeMeTmp
    
# 8. Load a FIMO output for a TF and determine if binding site is in promoter
# Retain: Entrez ID,Instances,Locations,Strandednesses
#tfs = ['motifHits/DistalBias/motifHits_Tcf7l2.1']
#tfs = ['motifHits/DistalBias/motifHits_Tcf7l2.1','motifHits/DistalBias/motifHits_SOX2_HMG_full_dimeric_17_1','motifHits/DistalBias/motifHits_V_OCT4_01_M01125']
#tfs = ['motifHits/ProximalBias/motifHits_ERG_ETS_full_dimeric_14_1','motifHits/ProximalBias/motifHits_V_NFYC_Q5_M02107','motifHits/ProximalBias/motifHits_E2F1_E2F_DBD_dimeric_12_1']
#tfs = ['motifHits/DistalBias/motifHits_Tcf7l2.1','motifHits/DistalBias/motifHits_Tcf7l2.2','motifHits/DistalBias/motifHits_SOX2_HMG_full_dimeric_17_1','motifHits/DistalBias/motifHits_V_OCT4_01_M01125','motifHits/ProximalBias/motifHits_ERG_ETS_full_dimeric_14_1','motifHits/ProximalBias/motifHits_V_NFYC_Q5_M02107','motifHits/ProximalBias/motifHits_E2F1_E2F_DBD_dimeric_12_1']
print('  Starting overlap analysis...')
cpus = cpu_count()
print('There are', cpus,'CPUs avialable.')
for tf in tfs:
    mgr = Manager()
    stats = mgr.dict()
    stats['tf'] = tf
    writeMe = mgr.dict()
    stdout.write('  ')
    #overlapFIMO('chr1')
    chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
    if not os.path.exists(output_dir + '/fullGenome_'+tf+'.csv'):
        print(tf)
        pool = Pool(processes=cpus)
        pool.map(overlapFIMOFaster,chroms)
        pool.close()
        pool.join()
        outFile = open(output_dir + '/fullGenome_'+tf+'.csv', 'w')
        outFile.write('Entrez ID,PromoterRegion,Chr,Instances,Locations,Strands,motifHitsPvalue') #--modified by yli
        outFile.write('\n'+'\n'.join(['\n'.join(writeMe[chrom]) for chrom in chroms]))
        outFile.close()
    del mgr
    stdout.write('\n')
print('  Done.')

