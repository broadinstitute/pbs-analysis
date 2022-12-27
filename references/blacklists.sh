#!/bin/bash

set -e

# Download blacklist files
# Website (URL): https://sites.google.com/site/anshulkundaje/projects/blacklists

hg38='https://www.encodeproject.org/files/ENCFF419RSJ/@@download/ENCFF419RSJ.bed.gz'
hg19='https://www.encodeproject.org/files/ENCFF001TDO/@@download/ENCFF001TDO.bed.gz'
mm10='https://www.encodeproject.org/files/ENCFF547MET/@@download/ENCFF547MET.bed.gz'
mm9='http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/mm9-mouse/mm9-blacklist.bed.gz'

localize () {
    # Expects file name first and file prefix second
    name=${1##*/}
    wget $1

    zcat $name > references/${2}_blacklist.bed

    rm $name
}

localize $hg38 hg38
localize $hg19 hg19
localize $mm10 mm10
localize $mm9 mm9