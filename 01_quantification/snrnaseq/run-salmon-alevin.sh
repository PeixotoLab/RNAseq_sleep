#$ -pe local 6
#$ -cwd
#$ -o log/
#$ -e log/
#$ -l mem_free=10G,h_vmem=10G

a=/fastscratch/myscratch/shicks1/RNAseq_sleep

# Use this if you are copying files from Dario's server
mice_samples=(1C 2E 3C 4E 5C 8E)
for c in "${mice_samples[@]}"
do
  ls -1 $a/fastq_files/peixoto_snRNA_SD_fastq/$c/*R1*.gz | sort -u > $a/01_quantification/snrnaseq/ID_fastq1
  ls -1 $a/fastq_files/peixoto_snRNA_SD_fastq/$c/*R2*.gz | sort -u > $a/01_quantification/snrnaseq/ID_fastq2
  ID_fastq1=`cat $a/01_quantification/snrnaseq/ID_fastq1`
  ID_fastq2=`cat $a/01_quantification/snrnaseq/ID_fastq2`
  salmon alevin --libType ISR \
      --index $a/salmon_index_files/gencode.vM25-salmon-index-v1.0.0-mouse-withdecoys-annotation-expanded \
      -1 $ID_fastq1 \
      -2 $ID_fastq2 \
      --tgMap $a/salmon_index_files/gencode.vM25.annotation.expanded.tx2gene.tsv \
      --chromiumV3 \
      --threads 6 \
      --output $a/01_quantification/snrnaseq/salmon_quants/${c}_quant \
      --forceCells 10000 \
      --dumpFeatures --dumpBfh \
      --numCellBootstraps 30
done

rm $a/01_quantification/snrnaseq/ID_fastq1
rm $a/01_quantification/snrnaseq/ID_fastq2