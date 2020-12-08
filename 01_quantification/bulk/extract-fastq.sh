#$ -pe local 6
#$ -cwd
#$ -o log/
#$ -e log/
#$ -l mem_free=10G,h_vmem=15G,h_fsize=100G

a=/fastscratch/myscratch/shicks1
b=/fastscratch/myscratch/shicks1/RNAseq_sleep/01_quantification

for file in $(<$b/bulk/data/SRR_files.txt)
do
    echo "$file.sra"
    fasterq-dump -O $a/RNAseq_sleep/fastq_files/ \
        -f --threads 6 \
        $a/sra/$file.sra
done