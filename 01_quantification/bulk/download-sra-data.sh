#$ -cwd
#$ -o log/
#$ -e log/
#$ -l mem_free=10G,h_vmem=10G,h_fsize=10G

d=/fastscratch/myscratch/shicks1/RNAseq_sleep/01_quantification/bulk

prefetch --option-file $d/data/SRR_files.txt