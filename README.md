# RNA-seq (bulk and snRNA) of sleep deprivation in wildtype and Shank3 mutant mice

This repository contains the code and analyses for the analysis of bulk RNA-seq and snRNA-seq analysis of sleep deprivation in wildtype and Shank3 mutant mice

## Authors

- Stephanie Hicks (shicks19@jhu.edu)
- Davide Risso (drisso@gmail.com)
- Katie Ford (kaitlyn.ford@wsu.edu)
- Elena Zuin (elena.zuin3@gmail.com)

## Data 

### bulk RNA-seq

- 2020-11-28: The FASTQ files for the WT and Shank 3 mutant bulk RNA-seq samples are available on SRA here: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA453921
- 2020-12-04: Downloaded all previous fastq files from Dario's server, but these also include nap mice. I noticed the ordering is slightly different and not sure what this means. 

Used the following code: 
```
# First login to JHPCE server using the ssh command:
$ ssh user@serverA

# Next, sftp to Dario's server:
$ sftp user@serverB

mget peixoto_bulk_SD_fastq/*.fastq.gz
mget peixoto_snRNA_SD_fastq/*
exit
```

If you use this approach, skip over next 3 sections, and start with the salmon sections.

#### Download GEO metadata

First, we manually downloaded the `SraRunInfo.csv` file: 

- Go here: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA453921
- Then click on `SRA Experiments`. Click `Send to`. Choose `File`.
- Change `Format` to `RunInfo`. Click `Create File`.
- This will download a file called `SraRunInfo.csv`.

This file has already been downloaded for you and is available at `01_quantification/bulk/data/SraRunInfo.csv`. 

Using this file, you can download and create GEO metadata with the shell script `01_quantification/bulk/download-geo-data.sh`.

#### Download SRA files

[Install](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit) the [SRA toolkit](https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/). 

**Notes**:

- You will need to add the Toolkit functions to your PATH variable (https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit). This can be done by adding `~/.../sratoolkit.2.10.6-centos_linux64/bin` to the PATH variable in `~/.bash_profile`. Check that it worked with `which fastq-dump`.
- You will also want to specify the default location to download the SRA files to using the toolkit.

Once you have installed the toolkit, run the shell script `01_quantification/bulk/download-sra-data.sh`. 
An example of a description of an SRA file is https://www.ncbi.nlm.nih.gov/sra/SRX4003606. 

#### Extract FASTQ files

References:

- https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump
- https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=fastq-dump

After downloading the SRA files, convert them to FASTQ files using the shell script `01_quantification/bulk/extract-fastq.sh`.

### snRNA-seq 

I downloaded all snRNA-seq fastq files from Dario's server. 

See `Notes.txt`: 

```
10k cells per sample
3 controls + 3 sleep deprived

20191220_Peixoto_1C_fastq.zip 
20200106_Peixoto_3C_fastq.zip
20200108_Peixoto_5C_fastq.zip
20200103_Peixoto_2E_fastq.zip
20200108_Peixoto_4E_fastq.zip
20200110_Peixoto_8E_fastq.zip

400M reads total
```

Used the following code: 
```
# First login to JHPCE server using the ssh command:
$ ssh user@serverA

# Next, sftp to Dario's server:
$ sftp user@serverB

mget peixoto_snRNA_SD_fastq/1C/*
mget peixoto_snRNA_SD_fastq/2E/*
mget peixoto_snRNA_SD_fastq/3C/*
mget peixoto_snRNA_SD_fastq/4E/*
mget peixoto_snRNA_SD_fastq/5C/*
mget peixoto_snRNA_SD_fastq/8E/*
exit
```

We use these fastq files directly in `salmon alevin` (see section below).


## Preprocessing 

This [website](https://combine-lab.github.io/alevin-tutorial/2018/setting-up-resources/) is helpful for setting up the resources needed to run `salmon` (bulk RNA-seq) and `salmon alevin` (snRNA-seq). 

### Download GENCODE files

Run the shell script `01_quantification/download-gencode-files.R`. This downloads GENCODE files and creates the files necessary for `salmon`. Reference files were obtained for mouse release version M25. 

1. `GRCm38.primary_assembly.genome.fa.gz` - nucleotide (DNA) sequences of the **GRCm38 primary genome assembly**.
2. `gencode.vM25.transcripts.fa.gz` - nucleotide (DNA) sequences of **all transcripts** on reference chromosomes.
3. `gencode.vM25.annotation.gtf.gz` - gene annotation on the reference chromosomes (i.e. for humans, these are chromosomes 1 to 22, X, and Y), i.e. locations of genes and other information about the genes, gene structure

> Gene transfer format (GTF) is a file format used to hold information about gene structure. It is a tab-delimited text format based on the general feature format (GFF), but contains some additional conventions specific to gene information.

The specific locations of where the files were pulled from are from here:

- `ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.primary_assembly.genome.fa.gz`
- `ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.transcripts.fa.gz`
- `ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz`


### Create decoys for `salmon` index

The decoy sequence is going to be the whole genome sequence (`GRCm38.primary_assembly.genome.fa.gz`). You can read more about decoy sequences in Salmon below:

* https://salmon.readthedocs.io/en/latest/salmon.html#preparing-transcriptome-indices-mapping-based-mode
* https://github.com/COMBINE-lab/SalmonTools/blob/master/README.md
* https://combine-lab.github.io/alevin-tutorial/2020/alevin-velocity/

Source for code: https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/

In general: 

> After creating the fasta file with transcript and intron sequences as above, 
we can index it using Salmon. Here, we add the full genome sequence as decoy 
sequences (see Srivastava et al., 2019 for more details). This is recommended 
in order to avoid reads truly originating from intergenic regions being assigned 
to a suboptimal transcriptome location. However, the effect of including decoys 
is typically smaller when both transcripts and introns are being quantified than 
for ‘regular’ gene expression quantification, since in the former case a larger 
fraction of the genome is already covered by the features of interest.

We use decoys for both bulk RNA-seq and snRNA-seq analyses.

To use a decoy, we need to create two files:

1. `decoys.txt` is the names of the genome targets (decoys), will be used in the `-d` parameter in `build-index-salmon.sh`
2. `gentrome_transcripts_mouse.fa.gz` is a concatenated FASTA transcriptome, will be used in the `-t` parameter in `build-index-salmon.sh` (see below). Note that you need to recreate this once per time you set up salmon to for quantification.

These two files are created in the `01_quantification/create-decoys-salmon.sh` and will both be used `01_quantification/build-index-salmon.sh` file to build the salmon index (see next section).

### Install and build `salmon index` 

This part will have to be done for each user. 
I installed the salmon 1.3.0 binary in my home directory here `/users/shicks1/src/`. 

To install salmon v1.3.0: 
```{bash}
cd /users/shicks1/src/
wget https://github.com/COMBINE-lab/salmon/releases/download/v1.3.0/salmon-1.3.0_linux_x86_64.tar.gz
tar xzvf salmon-1.3.0_linux_x86_64.tar.gz
rm salmon-1.3.0_linux_x86_64.tar.gz
```

Also, make sure this is in your `.bash_profile` file
```{bash}
PATH=$PATH:/users/shicks1/src/salmon-latest_linux_x86_64/bin
```

You can check to make sure salmon has been upgraded correctly using `salmon -h` inside terminal (or help with specific parts of using salmon using e.g. `salmon index -h` for help with the index step). 

OK, we are ready to use `salmon`. 

- The `-t` argument is the input transcripts file. 
- The `-i` argument is the index file to create. 
- The `-d` argument is the decoy sequence. 
- The `--keepDuplicates` argument forces all duplicate transcripts (for example, multiple unspliced transcript of the same gene that are identical for example) that appear in the input will be retained and quantified separately. If you keep the duplicates they will be assigned identical expression levels since salmon can’t tell them apart. When you aggregate on the gene level, this will not make a difference any more. Therefore, I do not keep the duplicates as we are interested in gene level aggregation. 
- The `--gencode` flag will handle the composite fasta headers in GENCODE transcript fasta files and split the transcript name at the first '|' character. 
- The `--threads` argument says how many threads to use when building the index. 

The salmon index is built using the `01_quantification/build-index-salmon.sh` file (used 4 cores). The bulk RNA-seq index uses the `decoys.txt` and is built from the combined FASTA file (`gentrome_transcripts_mouse.fa.gz`). The snRNA-seq index uses the `GRCm38.primary_assembly.genome.chrnames.txt` and is built from the combined FASTA file.

There are two indexes created: 

- bulk RNA-seq: `salmon_index_files/gencode.vM25-salmon-index-v1.0.0-mouse-withdecoy`
- snRNA-seq: `salmon_index_files/gencode.vM25-salmon-index-v1.0.0-mouse-withdecoys-annotation-expanded`


### bulk RNA-seq

#### Running `salmon quant`

We use the index (`gencode.vM25-salmon-index-v1.0.0-mouse-withdecoy`) created by `build-index-salmon.sh` for the bulk RNA-seq analysis.
See the `01_quantification/bulk/run-salmon-quant.sh` file.
This script quantifies reads at the transcript level, which can be summarized later at the gene level using the `tximeta::summarizeToGene()` function. 

#### Create `SummarizedExperiment` object

Here we import the `quant.sf` files produced by `salmon quant` into R/Bioconductor using the `tximeta` package with the `01_quantification/bulk/run-tximeta.R` script. 
We save the `SummarizedExperiment` object at the transcript level (and gene level) as a `.RDS` file: `data/se_mouse_sleep_complete.rds` (and `data/gse_mouse_sleep_complete.rds`)

- Helpful vignette: https://bioconductor.org/packages/devel/bioc/vignettes/tximeta/inst/doc/tximeta.html

### snRNA-seq 

#### Running `salmon alevin`

We use the index (`gencode.vM25-salmon-index-v1.0.0-mouse-withdecoys-annotation-expanded`) created by `build-index-salmon.sh` for the snRNA-seq analysis.
See the `01_quantification/snrnaseq/run-salmon-alevin.sh` file.
This script quantifies reads at the gene level (this is the only option available).

#### Create `SingleCellExperiment` object

Here we import the `quants_mat.gz` files produced by `salmon alevin` into R/Bioconductor using the `tximeta` package with the `01_quantification/snrnaseq/run-tximeta.R` script. 
We quantified reads for both spliced mRNA and introns using the `eisaR::getFeatureRanges()` function. 

We ran the `salmon alevin` output through `alevinQC` and produced a set of HTML files (one HTML file per sample) to report the QC for each quantification. 
These can be found at `01_quantification/snrnaseq/alevinQC/*.html`. 

We combined each sample into one `SingleCellExperiment` object. Here is the output: 

```
> sce
class: SingleCellExperiment 
dim: 87308 65989 
metadata(6): tximetaInfo quantInfo ... txomeInfo txdbInfo
assays(3): counts variance mean
rownames(87308): ENSMUSG00000102693.1 ENSMUSG00000064842.1 ...
  ENSMUSG00000099399.6-I ENSMUSG00000095366.2-I
rowData names(2): gene_id SYMBOL
colnames(65989): TTCGCTGCAAGTCGTT TATGTTCCACTATGTG ... GGGCTCATCGAAGTGG
  CCGGTGACACTTGAAC
colData names(2): whitelist sample_id
reducedDimNames(0):
altExpNames(0):
```

This `sce` object is not small: 

```
> pryr::object_size(sce)
6.81 GB
```

Gene symbols were added: 

```
> mcols(sce)
DataFrame with 87308 rows and 2 columns
                                      gene_id      SYMBOL
                                  <character> <character>
ENSMUSG00000102693.1     ENSMUSG00000102693.1          NA
ENSMUSG00000064842.1     ENSMUSG00000064842.1     Gm26206
ENSMUSG00000102851.1     ENSMUSG00000102851.1     Gm18956
ENSMUSG00000089699.1     ENSMUSG00000089699.1          NA
ENSMUSG00000103147.1     ENSMUSG00000103147.1      Gm7341
...                                       ...         ...
ENSMUSG00000102011.1-I ENSMUSG00000102011.1-I          NA
ENSMUSG00000100964.1-I ENSMUSG00000100964.1-I          NA
ENSMUSG00000099619.6-I ENSMUSG00000099619.6-I          NA
ENSMUSG00000099399.6-I ENSMUSG00000099399.6-I          NA
ENSMUSG00000095366.2-I ENSMUSG00000095366.2-I          NA
```

The meta data is sparse at this point, so this is all that is in the `colData` slot: 

```
> colData(sce)
DataFrame with 65989 rows and 2 columns
                 whitelist   sample_id
                 <logical> <character>
TTCGCTGCAAGTCGTT      TRUE          1C
TATGTTCCACTATGTG      TRUE          1C
AGGGTCCTCCAGCCTT      TRUE          1C
CACATGACATCGGTTA      TRUE          1C
GCAGTGGTATCAACGC      TRUE          1C
...                    ...         ...
CCTCACAGTCTGTAGT     FALSE          8E
GAGAAATCAAGTGACG     FALSE          8E
GCAGCTGAGGCACTAG     FALSE          8E
GGGCTCATCGAAGTGG     FALSE          8E
CCGGTGACACTTGAAC     FALSE          8E
```

We see there are roughly 9K-6K valid cell barcodes per sample. 

```
> table(sce$whitelist, sce$sample_id)
       
          1C   2E   3C   4E   5C   8E
  FALSE 3650 4033 4393 3780 1835 3654
  TRUE  7350 6964 6606 7214 9165 7345
```

We save the combined `SingleCellExperiment` object at the gene-level as a `.RDS` file: `data/sce_mouse_sleep_snrnaseq_complete.rds`. 

- Helpful vignette: https://combine-lab.github.io/alevin-tutorial/2020/alevin-velocity/

#### Data setting, quality control, normalization and doublets removal
To analyze the data the UMI counts of exons and the introns with the same ensembl IDs were added up. To identify mitochondrial genes, mouse ensembl IDs were converted to chromosome identifier with the EnsDb.Mmusculus.v79 (v2.99.0) package. We split the data into six SingleCellExperiment objects, one for each mouse.

For each sample, the Bioconductor scuttle (v1.8.4) package was used to detect low quality and damaged droplets. Particularly, the perCellQCMetrics function to calculate useful per-cell QC metrics: the sum of counts and the number of detected features. In this case, these statistics were computed also for mitochondrial genes. Then, we used the logNormCounts function of the scuttle package to apply the normalization factors and obtain the log-normalized counts.

Lastly, for each sample the doublets were removed with the scDblFinder (v1.12.0) package, using the computeDoubletDensity function to calculate the scores and the doubletThresholding function to set the doublet scores threshold with griffiths method.

#### Cell-type annotation
To identify cell types we used the Allen Whole Cortex & Hippocampus - 10x genomics (v2021) as reference dataset (cite https://www.sciencedirect.com/science/article/pii/S0092867421005018). This dataset was imported by the loadHDF5SummarizedExperiment function of the HDF5Array (v1.26.0) package and was converted to a SingleCellExperiment object available at https://github.com/drighelli/AllenInstituteBrainData. 
We then selected the “Non-Neuronal”, “Glutamatergic” and “GABAergic” clusters coming from the Visual Cortex (VIS, VISl, VISm, VISp) to annotate our dataset. For computational issues, we selected a random subset of 100,000 cortical cells.

Cell annotation was computed using two methods: Azimuth and SingleR. For the first method, the reference data was converted into a Seurat object and into a Azimuth compatible object, using the AzimuthReference function of the Azimuth (v0.4.6) package. Then query samples were merged and were converted into a Seurat object. Cell annotation was computed using the RunAzimuth function of the Azimuth package. The t-SNE and the UMAP projections were computed using the RunTSNE and RunUMAP functions of the Seurat (v 4.3.0) package with seed.use = 1.

For the second method, the reference dataset was aggregated across groups of cell type and was normalized, using the aggregateAcrossCells and the logNormCounts functions of the scuttle (v1.8.4) package, and the gene symbols were converted to ensembl with the function getBM of the biomaRt (v2.54.0) package. Then, cell annotation was computed using the SingleR function of SingleR (v2.0.0) package.

For each side, to visualize the assigned labels in two dimensions, the t-SNE and the UMAP projections were computed using the DimPlot function of Seurat package, with option reduction = "integrated_dr”, where "integrated_dr” is the supervised principal component analysis obtained by the Azimuth method.

Also, a pseudo-bulk level Multidimensional Scaling (MDS) plot was created with the pbMDS function of muscat (v1.12.1) package. Each point represents one subpopulation-sample instance; points are colored by subpopulation and shaped by group ID.

#### Differential expression analysis
For each neuronal cell-type with more than 500 cells, the differential gene expression analysis was carried out with two methods: a negative binomial generalized linear model (GLM) on pseudo-bulk samples and a negative binomial mixed effect model (NBMM) at the cell level.

For the GLM model, we created the pseudo-bulk samples with the function aggregateAcrossCells of the scuttle package. In other words, we computed sum counts values for each feature across cell-type and mouse groups. We made a Remove Unwanted Variation (RUV) normalization with k=1 on each neuronal cell-type, using the RUVs function of the RUVSeq (v1.32.0) package. We used the negative control genes coming from microarray analysis to estimate the factor of unwanted variation. The 10% negative control genes were randomly selected. The remaining control genes were used to fit RUV normalization.

To visualize the principal component analysis of RUV normalization, we used the plotPCA function of the EDASeq (v2.32.0) package.

We then used the Bioconductor edgeR (v3.40.2) package. Before the differential gene expression analysis, the genes were filtered with the function filterByExpr (with default parameters). The factor of unwanted variation was added in the design matrix. The differential gene expression analysis was computed with the function glmLRT by specifying “SD-HC” (Sleep Deprived vs Home Cage Control) as contrast and offset term equal to zero.

The negative binomial mixed effect model was computed with the function mmDS of the muscat (v1.12.1) package with option method=”nbinom”.

To visualize the differential expressed genes the volcano plot was made for each cell type and each model, using the ggplot2 (v3.3.6) package. Also, the p-value histogram was made for each cell-type and model, to check if the distribution was uniformly distributed between 0 and 1.

We used the negative and positive controls defined in the ‘Bulk genome-wide gene expression (RNA-seq)’ section to evaluate the concordance between the bulk and single-nuclear differential expression results.
