#!/bin/bash

set -e

# Download chromomsome sizes

hg38='http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes'
hg19='http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes'
hg18='http://hgdownload.cse.ucsc.edu/goldenPath/hg18/bigZips/hg18.chrom.sizes'
mm10='http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes'
mm9='http://hgdownload.cse.ucsc.edu/goldenPath/mm9/bigZips/mm9.chrom.sizes'

wget -P references $hg38
wget -P references $hg19
wget -P references $hg18
wget -P references $mm10
wget -P references $mm9