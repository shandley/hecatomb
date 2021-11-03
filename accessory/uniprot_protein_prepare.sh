# Manually Donwload via web interface: https://www.uniprot.org/uniprot/?query=taxonomy:%22Viruses%20[10239]%22&fil=reviewed%3Ano&sort=score
# Transfer to $DIR and rename to uniprot-viral.fa.gz

# Dependencies:
# CD-HIT: http://weizhongli-lab.org/cd-hit/
# mmseqs2: https://github.com/soedinglab/MMseqs2
# seqkit: https://bioinf.shenwei.me/seqkit/

# Set variables
DIR=/mnt/data1/databases/hecatomb/proteins/uniprot
TAX=/mnt/data1/databases/hecatomb/tax/taxonomy

# Timestamp file
now=`date +"%Y-%m-%d"`

echo The last time this script was run was: ${now} > ${DIR}/uniprot_protein_prepare.timestamp

### Fitler very short sequences
# Unzip as CD-HIT requires unzipped files
gunzip ${DIR}/uniprot-viral.fa.gz

# Remove proetien entries shorter than 30 amino acids
seqkit seq -m 30 ${DIR}/uniprot-viral.fa > ${DIR}/uniprot-viral-size-filtered.fa

# Cluster DB @ 99% ID
# Change -c flag to cluster at a different %ID
cd-hit -i ${DIR}/uniprot-viral-size-filtered.fa -o ${DIR}/virus_uniprot_c99.fa -c 0.99 -M 64000 -T 64;

# Compress clustered file
pigz ${DIR}/virus_uniprot_c99.fa

# Create mmseqs2 database
mmseqs createdb ${DIR}/virus_uniprot_c99.fa.gz ${DIR}/sequencdDB --dbtype 1;
mmseqs createtaxdb ${DIR}/sequenceDB ${DIR}/tmp --ncbi-tax-dump ${TAX};
