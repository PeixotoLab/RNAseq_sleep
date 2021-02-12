#$ -pe local 4
#$ -cwd
#$ -o log/
#$ -e log/
#$ -l mem_free=10G,h_vmem=10G

###############################################
########## Files here are for bulk ############
# to build salmon index
# Use
# #$ -l mem_free=10G,h_vmem=10G
###############################################

# create salmon index with decoys (this process takes ~1 hour for mouse)
salmon index -t /fastscratch/myscratch/shicks1/RNAseq_sleep/salmon_index_files/gentrome_transcripts_mouse.fa.gz \
             -d /fastscratch/myscratch/shicks1/RNAseq_sleep/salmon_index_files/decoys_mouse.txt \
             -i /fastscratch/myscratch/shicks1/RNAseq_sleep/salmon_index_files/gencode.vM25-salmon-index-v1.0.0-mouse-withdecoys \
             --gencode --threads 4 -k 31

# create salmon index without decoys (this process takes <30 mins for mouse)
# salmon index -t /fastscratch/myscratch/shicks1/RNAseq_sleep/salmon_index_files/gencode.vM25.transcripts.fa.gz \
#              -i /fastscratch/myscratch/shicks1/RNAseq_sleep/salmon_index_files/gencode.vM25-salmon-index-v1.0.0-mouse-nodecoys \
#              --gencode --threads 4 -k 31


###############################################
######## Files here are for snRNAseq ##########
# Following https://combine-lab.github.io/alevin-tutorial/2020/alevin-velocity/
# to build salmon index
# Use
# #$ -l mem_free=10G,h_vmem=10G
###############################################

# cd /fastscratch/myscratch/shicks1/RNAseq_sleep/salmon_index_files

# # create salmon index with decoys (this process takes ~1 hour for mouse)
# salmon index -t <(cat gencode.vM25.annotation.expanded.fa GRCm38.primary_assembly.genome.fa) \
#              -d GRCm38.primary_assembly.genome.chrnames.txt \
#              -i gencode.vM25-salmon-index-v1.0.0-mouse-withdecoys-annotation-expanded \
#              --gencode --threads 4 -k 31
