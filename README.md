# tfbsdb2-DataProcess

Data process for constructing the transcription factor (TF) to target gene interactions database.

## Dependencies

Python version: 3.7.15

FIMO version: 5.0.5

## Data Preparation
### Digital footprints

wget https://resources.altius.org/~jvierstra/projects/footprinting.2020/consensus.index/consensus_footprints_and_collapsed_motifs_hg38.bed.gz

### PSSMs
#### SELEX

wget https://downloads.wenglab.org/factorbook-download/all-selex-motifs.meme.gz

#### JASPAR
#### TRANSFAC
#### UniPROBE

### mask reference

python maskSequences.py

## Step 1: build database

python buildTFBSDB.TFBSDB2.py

## Step 2: target identification

python targetIdentification_mp.TFBSDB2.py --outputDir targetIdentification_tests/targetIdentification_1000_1000 --promoterStart 1000 --promoterEnd -1000 --proximalStart 2500 --proximalEnd -500 --distalStart 5000 --distalEnd 2500
