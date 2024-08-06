#!/bin/bash
# Get the FASTA file suggested by the user
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/GRCh38.p13.genome.fa.gz
gzip -d GRCh38.p13.genome.fa.gz