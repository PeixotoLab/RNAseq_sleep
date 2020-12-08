#$ -pe local 6
#$ -cwd
#$ -o log/
#$ -e log/
#$ -l mem_free=10G,h_vmem=10G

a=/fastscratch/myscratch/shicks1/RNAseq_sleep

# Use this if you are copying files from Dario's server
for file in `ls $a/fastq_files/peixoto_bulk_SD_fastq/*.fastq.gz | sed 's/[12].fastq.gz//' | sort -u`;
do 
    fastq_1=`echo $file | sed "s/$/1.fastq.gz/"`
    fastq_2=`echo $file | sed "s/$/2.fastq.gz/"`
    output_file=`basename $file | sed "s/_R//"`
    salmon quant -l A \
        --index $a/salmon_index_files/gencode.vM25-salmon-index-v1.0.0-mouse-withdecoys \
        -1 $fastq_1 \
        -2 $fastq_2 \
        --threads 6 \
        --output $a/01_quantification/bulk/salmon_quants/${output_file}_quant \
        --numBootstraps 30
done

# # Use this if you are downloading SRA files
# for file in $(<$a/01_quantification/bulk/data/SRR_files.txt)
# do
#     fastq_1=`echo $file | sed "s/$/.sra_1.fastq/"`
#     fastq_2=`echo $file | sed "s/$/.sra_2.fastq/"`
#     salmon quant -l A \
#         --index $a/salmon_index_files/gencode.vM25-salmon-index-v1.0.0-mouse-withdecoys \
#         -1 $a/fastq_files/$fastq_1 \
#         -2 $a/fastq_files/$fastq_2 \
#         --threads 6 \
#         --output $a/01_quantification/bulk/salmon_quants/${file}_quant \
#         --numBootstraps 30
# done

