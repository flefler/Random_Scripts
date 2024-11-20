# Cyanobacterial genomic assembly script
### Last updated November 15 2024 by Forrest W. Lefler

Cyanobacterial cultures are, more often than not, a uni-algal culture with a various assemblage of other (potentially beneficial) microbes, primarily bacteria. Because of this, cyanobacterial genomes can be difficult to recover.
Read more [here](https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2017.00809/full), and [here](https://www.sciencedirect.com/science/article/abs/pii/S2211926424000274)


To accomplish this goal we will take the following steps:

First, we will remove low quality reads with [fastP](https://github.com/OpenGene/fastp), then use [HoCoRT](https://github.com/ignasrum/hocort) to remove potential contaminate reads like human, mouse, and phiX, followed by analyses with [fastQC](https://github.com/s-andrews/FastQC), and tie it all together with [multiQC](https://github.com/MultiQC/MultiQC) 

Second, we will look at the taxonomic characteristics of the raw reads using [singleM](https://github.com/wwood/singlem). We can refer back to this if our bins lack cyanobacteria. Also, they just look cool.

Third, we will assemble contigs with [metaSPAdes](https://github.com/ablab/spades) and build a report of the output using [quast](https://github.com/ablab/quast) and [multiQC](https://github.com/MultiQC/MultiQC) and make BAM files with [CoverM](https://github.com/wwood/CoverM)

Fourth, we will identify viral and plasmid contigs using [geNomad](https://github.com/apcamargo/genomad).

Fifth, we will bin using a few binners ([SemiBin2](https://github.com/BigDataBiology/SemiBin), [METABAT2](https://pubmed.ncbi.nlm.nih.gov/31388474/), [Comebin](https://github.com/ziyewang/COMEBin)], and [metacoag](https://github.com/metagentools/MetaCoAG)) and the reconcile those bins using [binette](https://github.com/genotoul-bioinfo/Binette).

Finally, we will use [compareM2](https://github.com/cmkobel/comparem2) to analyze our bins with [GTDBtk](https://github.com/Ecogenomics/GTDBTk), [checkM2](https://github.com/chklovski/CheckM2), [SeqKit](https://github.com/shenwei356/seqkit) and get a nice HTML file. and then annotate with [metacerberus](https://github.com/raw-lab/metacerberus)



## Other info

This script is for illumina data, primarily NovoSeq data. 

```00_Reads``` contains our raw reads, with one folder per sample. ```~/00_Reads/${SAMPLE}/*.fq.gz```

```99_logs``` contains [logs](https://en.wikipedia.org/wiki/Logging#/media/File:Mexico_logs.jpg)

I am submitting jobs via [slurm](https://en.wikipedia.org/wiki/Slurm_Workload_Manager), so these scripts include [loops](https://www.geeksforgeeks.org/loops-programming/) and sbatch commands.

samples_1.txt are the names of my samples, it looks like this:

Sample names\
XYZ123\
ABC098

# Assembly portion

## QA/QC

### fastP
We want to deduplicate the reads to reduce redundancy in assembly, less is more. we also want the -g option as these were sequenced with NovaSeq. Adjust according to your data.
```
mamba activate hocort
SAMPLES=`cut -f 1 samples.txt | sed '1d'`
for SAMPLE in $SAMPLES; do

    N=${SAMPLE}_fastp
    mkdir 00_CLEANREADS/fastP
    mkdir 00_CLEANREADS/fastP/HTML
    mkdir 00_CLEANREADS/fastP/JSON
    mkdir 00_CLEANREADS/fastP/${SAMPLE}

    F=00_Reads/${SAMPLE}/*_1.fq.gz
    R=00_Reads/${SAMPLE}/*_2.fq.gz
    F_cp=00_CLEANREADS/fastP/${SAMPLE}/${SAMPLE}_clean_1.fq.gz
    R_cp=00_CLEANREADS/fastP/${SAMPLE}/${SAMPLE}_clean_2.fq.gz
    HTML=00_CLEANREADS/fastP/HTML/${SAMPLE}.html
    JSON=00_CLEANREADS/fastP/JSON/${SAMPLE}.json

    CMD="fastp -i ${F} -I ${R} -o ${F_cp} -O ${R_cp} -D -g -j ${JSON} -h ${HTML} -w 16"

    sbatch -A ${ACCOUNT} -J ${N} -c 16 --mem=100G -o 99_logs/${N}_o.txt -e 99_logs/${N}_e.txt --export=ALL --mail-type=FAIL --mail-user=${EMAIL} -t 24:00:00 --wrap="${CMD}"

done
```
### hocort
Remove human, mouse, and phix reads with Hocort, remove low quality and duplicate reads with fastP
```
mamba activate hocort
SAMPLES=`cut -f 1 samples.txt | sed '1d'`
for SAMPLE in $SAMPLES; do

    mkdir 00_CLEANREADS/hocort
    mkdir 00_CLEANREADS/hocort/${SAMPLE}

    N=${SAMPLE}_decontam
    F_cp=00_CLEANREADS/fastP/${SAMPLE}/${SAMPLE}_clean_1.fq.gz
    R_cp=00_CLEANREADS/fastP/${SAMPLE}/${SAMPLE}_clean_2.fq.gz
    REF=/blue/${ACCOUNT}/flefler/contam_databases/homo_mus_phix/homo_mus_phix
    F_ch=00_CLEANREADS/hocort/${SAMPLE}/${SAMPLE}_clean_1.fq.gz
    R_ch=00_CLEANREADS/hocort/${SAMPLE}/${SAMPLE}_clean_2.fq.gz

    CMD="hocort map bowtie2 -x ${REF} -i ${F_cp} ${R_cp} -o ${F_ch} ${R_ch} -c=--very-fast"

    sbatch -A ${ACCOUNT} --qos=${ACCOUNT}-b -J ${N} -c 20 --mem=150G -o ${N}_o.txt -e ${N}_e.txt --export=ALL --mail-type=ALL --mail-user=${EMAIL} -t 72:00:00 --wrap="${CMD}"

done
```

### fastQC
```
mamba activate Assemble_Bin

mkdir 01_QC/FASTQC

SAMPLES=`cut -f 1 samples_1.txt | sed '1d'`
for SAMPLE in $SAMPLES; do

    N=${SAMPLE}_FASTQC

    CMD="fastqc --quiet 00_CLEANREADS/hocort/${SAMPLE}/*.fq.gz -o 01_QC/FASTQC"

    sbatch -A ${ACCOUNT} -J ${N} -c 8 --mem=32G -o 99_logs/${N}_o.txt -e 99_logs/${N}_e.txt --export=ALL --mail-type=FAIL --mail-user=${EMAIL} -t 2:00:00 --wrap="${CMD}"

done
```
### multiQC
```
multiqc 01_QC -o 01_QC --interactive
```

## Now lets look at the data using singleM
This helps understand who is there, its useful if we do not recover a cyano genome
```
mkdir 02_singleM

mamba activate Assemble_Bin
SAMPLES=`cut -f 1 samples_1.txt | sed '1d'`
for SAMPLE in $SAMPLES; do

    N=${SAMPLE}_singleM
    
    F_ch=00_CLEANREADS/hocort/${SAMPLE}/${SAMPLE}_clean_1.fq.gz
    R_ch=00_CLEANREADS/hocort/${SAMPLE}/${SAMPLE}_clean_2.fq.gz

	CMD="singlem pipe \
	--forward ${F_ch} --reverse ${R_ch} \
	--taxonomic-profile-krona 02_singleM/${SAMPLE}_krona.html \
	--threads 8 --quiet"

    sbatch -A ${ACCOUNT} -J ${N} -c 8 --mem=32G -o 99_logs/${N}_o.txt -e 99_logs/${N}_e.txt --export=ALL --mail-type=FAIL --mail-user=${EMAIL} -t 24:00:00 --wrap="${CMD}"

done
```

## Assembly
### metaSPAdes
```
mamba activate Assemble_Bin
mkdir 03_ASSEMBLIES

SAMPLES=`cut -f 1 samples_1.txt | sed '1d'`

for SAMPLE in $SAMPLES; do

    N=${SAMPLE}_metaSPAdes

    F_ch=00_CLEANREADS/fastP/${SAMPLE}/${SAMPLE}_clean_1.fq.gz
    R_ch=00_CLEANREADS/fastP/${SAMPLE}/${SAMPLE}_clean_2.fq.gz

    CMD="spades.py --meta -1 ${F_ch} -2 ${R_ch} -o 03_ASSEMBLIES/${SAMPLE} -m 250 -t 32"
    
    sbatch -A ${ACCOUNT} -J ${N} -c 32 --mem=250G -o 99_logs/${N}_o.txt -e 99_logs/${N}_e.txt --export=ALL --mail-type=FAIL --mail-user=${EMAIL} -t 48:00:00 --wrap="${CMD}"

done
```


### gzip all fasta files and subset contigs >1000bp
Spades does not have a minumum contig size, we are going to only use those >1000 bp for further analyses. Furthermore, spades produces quite a few output files, this is really annoying so were going to zip them up.
```
SAMPLES=`cut -f 1 samples_1.txt | sed '1d'`
for SAMPLE in $SAMPLES; do
    N=${SAMPLE}_gzip
    CMD="seqkit seq -m 1000 03_ASSEMBLIES/${SAMPLE}/scaffolds.fasta -o 03_ASSEMBLIES/${SAMPLE}/cleaned_scaffolds.fasta.gz &&\
    gzip 03_ASSEMBLIES/${SAMPLE}/*.fasta &&\
    gzip 03_ASSEMBLIES/${SAMPLE}/*/*.fasta"
    sbatch -A ${ACCOUNT} -J ${N} -c 4 --mem=10G -o 99_logs/${N}_o.txt -e 99_logs/${N}_e.txt --export=ALL --mail-type=FAIL --mail-user=${EMAIL} -t 48:00:00 --wrap="${CMD}"
done
```

### Assess Assemblies with Quast
```
mkdir 01_QC/QUAST

mamba activate Assemble_Bin
SAMPLES=`cut -f 1 samples_1.txt | sed '1d'`

for SAMPLE in $SAMPLES; do

    N=${SAMPLE}_quast

    CMD="metaquast 03_ASSEMBLIES/${SAMPLE}/cleaned_scaffolds.fasta.gz \
        --output-dir 01_QC/QUAST/${SAMPLE} \
        --max-ref-number 0 -L \
        --threads 8 --fast --silent"
    
    sbatch -A ${ACCOUNT} -J ${N} -c 8 --mem=30G -o 99_logs/${N}_o.txt -e 99_logs/${N}_e.txt --export=ALL --mail-type=FAIL --mail-user=${EMAIL} -t 2:00:00 --wrap="${CMD}"

done
```
multiqc 01_QC -o 01_QC --interactive


## Identify viral and plasmid contigs with geNomad
### run geNomad
I dont want to identify proviruses and i want hits to have 3 genes to be considered a plasmid/virus
```
mamba activate genomad
mkdir 04_geNomad
SAMPLES=`cut -f 1 samples_1.txt | sed '1d'`
for SAMPLE in $SAMPLES; do

    N=${SAMPLE}_geNomad
    
    CMD="genomad end-to-end --cleanup --disable-find-proviruses --quiet --enable-score-calibration --composition metagenome --threads 16 \
    --min-score 0.8 --min-virus-hallmarks 3 --min-plasmid-hallmarks 3 --min-plasmid-hallmarks-short-seqs 3 --min-virus-hallmarks-short-seqs 3 \
    03_ASSEMBLIES/${SAMPLE}/cleaned_scaffolds.fasta.gz 04_geNomad/${SAMPLE} /blue/${ACCOUNT}/flefler/genomad_db"

    sbatch -A ${ACCOUNT} -J ${N} -c 16 --mem=100G -o 99_logs/${N}_o.txt -e 99_logs/${N}_e.txt --export=ALL --mail-type=FAIL --mail-user=${EMAIL} -t 196:00:00  --wrap="${CMD}"

done
```

### subset viral and plasmid contigs

#### Extract viral contigs
```
mamba activate Assemble_Bin
SAMPLES=`cut -f 1 samples_1.txt | sed '1d'`
for SAMPLE in $SAMPLES; do

    N=${SAMPLE}_viral

    CMD="seqkit grep -f <(cut -f1 04_geNomad/${SAMPLE}/cleaned_scaffolds_summary/scaffolds_virus_summary.tsv) 03_ASSEMBLIES/${SAMPLE}/cleaned_scaffolds.fasta.gz -o 04_geNomad/${SAMPLE}/virus/viral_contigs.fasta.gz"

    sbatch -A ${ACCOUNT} -J ${N} -c 4 --mem=10G -o 99_logs/${N}_o.txt -e 99_logs/${N}_e.txt --export=ALL --mail-type=ALL --mail-user=${EMAIL} -t 196:00:00  --wrap="${CMD}"

done
```

#### Extract plasmid contigs
```
mamba activate Assemble_Bin
SAMPLES=`cut -f 1 samples_1.txt | sed '1d'`
for SAMPLE in $SAMPLES; do

    N=${SAMPLE}_plasmid

    CMD="seqkit grep -f <(cut -f1 04_geNomad/${SAMPLE}/cleaned_scaffolds_summary/scaffolds_plasmid_summary.tsv) 03_ASSEMBLIES/${SAMPLE}/cleaned_scaffolds.fasta.gz \
    | gzip > 04_geNomad/${SAMPLE}/plasmid/plasmid_contigs.fasta.gz"

    sbatch -A ${ACCOUNT} -J ${N} -c 4 --mem=10G -o 99_logs/${N}_o.txt -e 99_logs/${N}_e.txt --export=ALL --mail-type=FAIL --mail-user=${EMAIL} -t 196:00:00  --wrap="${CMD}"

done
```

#### Combine viral and plasmid contig names
```
SAMPLES=`cut -f 1 samples_1.txt | sed '1d'`
for SAMPLE in $SAMPLES; do
    cat 04_geNomad/${SAMPLE}/cleaned_scaffolds_summary/cleaned_scaffolds_virus_summary.tsv 04_geNomad/${SAMPLE}/cleaned_scaffolds_summary/cleaned_scaffolds_plasmid_summary.tsv \
    | cut -f 1 | grep -v "seq_name" |  sort | uniq > 04_geNomad/${SAMPLE}/viral_and_plasmid_contigs.txt
done
```

#### Extract contigs that are not in the viral_and_plasmid_contigs.txt file
```
SAMPLES=`cut -f 1 samples_1.txt | sed '1d'`
for SAMPLE in $SAMPLES; do

    N=${SAMPLE}_removevirusplasmid

    CMD="seqkit grep -v -f 04_geNomad/${SAMPLE}/viral_and_plasmid_contigs.txt 03_ASSEMBLIES/${SAMPLE}/cleaned_scaffolds.fasta.gz \
    -o 03_ASSEMBLIES/${SAMPLE}/non_viral_non_plasmid_contigs.fasta.gz"

    sbatch -A ${ACCOUNT} -J ${N} -c 4 --mem=10G -o 99_logs/${N}_o.txt -e 99_logs/${N}_e.txt --export=ALL --mail-type=FAIL --mail-user=${EMAIL} -t 196:00:00  --wrap="${CMD}"
    
done
```

## Make BAM files
```
mkdir 05_BIN
conda activate reassemble
SAMPLES=`cut -f 1 samples_1.txt | sed '1d'`
for SAMPLE in $SAMPLES; do

    N=${SAMPLE}_BAM

    CONTIGS=03_ASSEMBLIES/${SAMPLE}/non_viral_non_plasmid_contigs.fasta.gz
    F_cp=00_CLEANREADS/fastP/${SAMPLE}/${SAMPLE}_clean_1.fq.gz
    R_cp=00_CLEANREADS/fastP/${SAMPLE}/${SAMPLE}_clean_2.fq.gz

    CMD="coverm make -r ${CONTIGS} -1 ${F_cp} -2 ${R_cp} -o 05_BIN/${SAMPLE} -t 16"

    sbatch -A ${ACCOUNT} -J ${N} -c 16 --mem=100G -o 99_logs/${N}_o.txt -e 99_logs/${N}_e.txt --export=ALL --mail-type=FAIL --mail-user=${EMAIL} -t 24:00:00 --wrap="${CMD}"

done
```

### bam file stats
```
SAMPLES=`cut -f 1 samples_1.txt | sed '1d'`
for SAMPLE in $SAMPLES; do

    mkdir 01_QC/Bamstats/${SAMPLE}
    N=${SAMPLE}_Bamstats

    CMD="bamtools stats -in 05_BIN/${SAMPLE}/*.bam > 01_QC/Bamstats/${SAMPLE}stats.txt"
    sbatch -A ${ACCOUNT} -J ${N} -c 2 --mem=10G -o 99_logs/${N}_o.txt -e 99_logs/${N}_e.txt --export=ALL --mail-type=FAIL --mail-user=${EMAIL} -t 24:00:00 --wrap="${CMD}"
done
```

## MultiQC
```
multiqc 01_QC -o 01_QC --interactive 
```

## Taxonomy of contigs
I did this for no reason. Its cool to see but not all that useful.
```
mamba activate metabuli
SAMPLES=`cut -f 1 samples_3.txt | sed '1d'`
for SAMPLE in $SAMPLES; do

    mkdir 06_METABULI/${SAMPLE}

    N=${SAMPLE}_METABULI
    CONTIGS=03_ASSEMBLIES/${SAMPLE}/non_viral_non_plasmid_contigs.fasta.gz
    RESULTS=06_METABULI/${SAMPLE}
    DB=/blue/${ACCOUNT}/flefler/metabuliGTDB/gtdb
    TMP=/blue/${ACCOUNT}/flefler/BLCC_Genomes/06_METABULI/${SAMPLE}/tmp

    CMD="metabuli classify --seq-mode 3 ${CONTIGS} ${DB} ${RESULTS} ${SAMPLE}"

    sbatch -A ${ACCOUNT} -J ${N} -c 2 --mem=100G \
    -o 99_logs/${N}_o.txt -e 99_logs/${N}_e.txt --export=ALL --mail-type=FAIL --mail-user=${EMAIL} -t 24:00:00 --wrap="${CMD}"

done
```

# BINNING STEPS
## BIN

### binners
#### METABAT2
```
conda activate metabat2
SAMPLES=`cut -f 1 samples_1.txt | sed '1d'`
for SAMPLE in $SAMPLES; do

    mkdir 05_BIN/${SAMPLE}/METABAT

    N=METABAT_${SAMPLE}

    CONTIGS=03_ASSEMBLIES/${SAMPLE}/non_viral_non_plasmid_contigs.fasta.gz
    BAM=05_BIN/${SAMPLE}/*.bam

    CMD="jgi_summarize_bam_contig_depths ${BAM} --referenceFasta ${CONTIGS} --outputDepth 05_BIN/${SAMPLE}/METABAT/depth.txt &&\
    metabat2 -i ${CONTIGS} -a 05_BIN/${SAMPLE}/METABAT/depth.txt -o 05_BIN/${SAMPLE}/METABAT/BINS"

    sbatch -A ${ACCOUNT} -J ${N} -c 16 --mem=100G -o 99_logs/${N}_o.txt -e 99_logs/${N}_e.txt --export=ALL --mail-type=FAIL --mail-user=${EMAIL} -t 196:00:00  --wrap="${CMD}"

done

for SAMPLE in $SAMPLES; do
    gzip 05_BIN/${SAMPLE}/METABAT/*.fa
done
```

#### SEMIBIN2
```
conda activate Assemble_Bin
SAMPLES=`cut -f 1 samples_1.txt | sed '1d'`
for SAMPLE in $SAMPLES; do

    N=${SAMPLE}_SemiBin
    CONTIGS=03_ASSEMBLIES/${SAMPLE}/non_viral_non_plasmid_contigs.fasta.gz
    BAM=05_BIN/${SAMPLE}/*.bam

    CMD="SemiBin2 single_easy_bin \
    --input-fasta ${CONTIGS} --input-bam ${BAM} \
    --environment wastewater --quiet --output 05_BIN/${SAMPLE}/SEMIBIN"

    sbatch -A ${ACCOUNT} -J ${N} -c 16 --mem=100G -o 99_logs/${N}_o.txt -e 99_logs/${N}_e.txt --export=ALL --mail-type=FAIL --mail-user=${EMAIL} -t 196:00:00 --wrap="${CMD}"
done
```

#### comebin dont run
I typically dont run this as it is very time consuming
```
#conda activate comebin_env
#SAMPLES=`cut -f 1 samples_1.txt | sed '1d'`
#for SAMPLE in $SAMPLES; do

    #gunzip 03_ASSEMBLIES/${SAMPLE}/non_viral_non_plasmid_contigs.fasta.gz

    #N=comebin_${SAMPLE}
    #CONTIGS=03_ASSEMBLIES/${SAMPLE}/non_viral_non_plasmid_contigs.fasta

    #CMD="run_comebin.sh -a ${CONTIGS} \
    #-o 05_BIN/${SAMPLE}/comebin \
    3-p 05_BIN/${SAMPLE} -t 16"
    #sbatch -A ${ACCOUNT} -J ${N} -c 16 --mem=100G -o 99_logs/${N}_o.txt -e 99_logs/${N}_e.txt \
    #--export=ALL --mail-type=FAIL --mail-user=${EMAIL} -t 200:00:00 --wrap="${CMD}"

#done
```
#### metacoag
This one uses the unfiltered contigs
#### Prepare the bam files
```
mamba activate reassemble
SAMPLES=`cut -f 1 samples_1.txt | sed '1d'`
SAMPLES=F202
for SAMPLE in $SAMPLES; do

    #gunzip 03_ASSEMBLIES/${SAMPLE}/scaffolds.fasta.gz

    N=${SAMPLE}_covermcontig
    F_cp=00_CLEANREADS/fastP/${SAMPLE}/${SAMPLE}_clean_1.fq.gz
    R_cp=00_CLEANREADS/fastP/${SAMPLE}/${SAMPLE}_clean_2.fq.gz
    CONTIGS=03_ASSEMBLIES/${SAMPLE}/scaffolds.fasta
    ABUNDANCE=05_BIN/${SAMPLE}/abundance.tsv

    CMD="coverm contig -1 ${F_cp} -2 ${R_cp} -r ${CONTIGS} -o ${ABUNDANCE} -t 8 
    sed -i '1d' ${ABUNDANCE}"

    sbatch -A ${ACCOUNT} -J ${N} -c 8 --mem=50G -o 99_logs/${N}_o.txt -e 99_logs/${N}_e.txt \
    --export=ALL --mail-type=FAIL --mail-user=${EMAIL} -t 24:00:00 --wrap="${CMD}"

done
```
#### Run the binner
```
mamba activate gbintk
SAMPLES=`cut -f 1 samples_1.txt | sed '1d'`
for SAMPLE in $SAMPLES; do

    mkdir 05_BIN/${SAMPLE}/metacoag

    N=${SAMPLE}_metacoag
    CONTIGS=03_ASSEMBLIES/${SAMPLE}/scaffolds.fasta
    GRAPH=03_ASSEMBLIES/${SAMPLE}/assembly_graph_with_scaffolds.gfa
    PATHS=03_ASSEMBLIES/${SAMPLE}/scaffolds.paths
    ABUNDANCE=05_BIN/${SAMPLE}/abundance.tsv

    CMD="gbintk metacoag --assembler spades --graph ${GRAPH} --contigs ${CONTIGS} --paths ${PATHS} --abundance ${ABUNDANCE} --output 05_BIN/${SAMPLE}/metacoag"

    sbatch -A ${ACCOUNT} -J ${N} -c 8 --mem=50G -o 99_logs/${N}_o.txt -e 99_logs/${N}_e.txt \
    --export=ALL --mail-type=FAIL --mail-user=${EMAIL} -t 24:00:00 --wrap="${CMD}"

done
```

### binette
```
conda activate binette

SAMPLES=`cut -f 1 samples_1.txt | sed '1d'`
for SAMPLE in $SAMPLES; do

    N=${SAMPLE}_binette

    CONTIGS=03_ASSEMBLIES/${SAMPLE}/scaffolds.fasta.gz

    CMD="binette --bin_dirs \
    05_BIN/${SAMPLE}/SEMIBIN/output_bins \
    05_BIN/${SAMPLE}/METABAT \
    05_BIN/${SAMPLE}/metacoag/bins \
    --contigs ${CONTIGS} \
    -m 80 -t 8 \
    -o 05_BIN/${SAMPLE}/binette"

    sbatch -A ${ACCOUNT} -J ${N} -c 8 --mem=100G -o 99_logs/${N}_o.txt -e 99_logs/${N}_e.txt \
    --export=ALL --mail-type=FAIL --mail-user=${EMAIL} -t 24:00:00 --wrap="${CMD}"

done
```

SAMPLES=`cut -f 1 samples_1.txt | sed '1d'`
for SAMPLE in $SAMPLES; do
    rm -rf 05_BIN/${SAMPLE}/binette
done

### Zip it all up
```
SAMPLES=`cut -f 1 samples_1.txt | sed '1d'`
for SAMPLE in $SAMPLES; do

    N=${SAMPLE}_gzip

    CMD="gzip 05_BIN/${SAMPLE}/SEMIBIN/output_bins/*.fa
    gzip 05_BIN/${SAMPLE}/METABAT/*.fa
    gzip 05_BIN/${SAMPLE}/metacoag/bins/*.fasta
    gzip 05_BIN/${SAMPLE}/binette/final_bins/*.fa
    gzip 05_BIN/${SAMPLE}/metacoag/low_quality_bins/*.fasta
    gzip 03_ASSEMBLIES/${SAMPLE}/*.fasta"

    sbatch -A ${ACCOUNT} -J ${N} -c 2 --mem=20G -o 99_logs/${N}_o.txt -e 99_logs/${N}_e.txt \
    --export=ALL --mail-type=FAIL --mail-user=${EMAIL} -t 24:00:00 --wrap="${CMD}"

done
```

### move and rename bins based on sample name
```
SAMPLES=`cut -f 1 samples_1.txt | sed '1d'`
for SAMPLE in $SAMPLES; do 
    for file in 05_BIN/${SAMPLE}/binette/final_bins/*.fa.gz; do 
        cp "$file" "06_GENOMES/${SAMPLE}_$(basename "$file")"; 
    done
done
```

# Annotate bins

## Simple stats
```
seqkit stats -a 06_GENOMES/\*\.fa.gz -T | csvtk tab2csv | csvtk cut -f -2,-3,-9,-10,-11,-12,-15,-16,-17 -o Genome_Info.csv
```
## compareM2
```
mamba activate comparem2
N=comparem2
CMD="comparem2 --config input_genomes=*.fa.gz --until sequence_lengths checkm2 gtdbtk"
sbatch -A ${ACCOUNT} --qos=${ACCOUNT}-b -p bigmem -J ${N} -c 128 --mem=1000G -o ${N}_o.txt -e ${N}_e.txt --export=ALL --mail-type=ALL --mail-user=${EMAIL} -t 24:00:00 --wrap="${CMD}"
```
## annotate genes
```
why is this missing
```

# ultra secret bonus code
We can map the raw reads back to the assembly and reassemble those with spades (or whatever you like) to try and imporve the assembly. sometimes its not worth the effort though.
```
SAMPLES=`cut -f 1 samples_1.txt | sed '1d'`
for SAMPLE in $SAMPLES; do 

    mkdir 07_Reassemble
    mkdir 08_REASSEMBLIES
    mkdir 09_FINALGENOMES
    mkdir 07_Reassemble/${SAMPLE}

    for file in 05_BIN/${SAMPLE}/binette/final_bins/*.fa.gz; do 

        base_name=${SAMPLE}_$(basename "$file" .fa.gz)
        N=${base_name}_fastq

        F_cp=00_CLEANREADS/fastP/${SAMPLE}/${SAMPLE}_clean_1.fq.gz
        R_cp=00_CLEANREADS/fastP/${SAMPLE}/${SAMPLE}_clean_2.fq.gz
        BAM1=07_Reassemble/${SAMPLE}/${base_name}/*.bam
        BAM2=07_Reassemble/${SAMPLE}/${base_name}/${base_name}_pairs.bam

        CMD="coverm make -r ${file} -1 ${F_cp} -2 ${R_cp} -o 07_Reassemble/${SAMPLE}/${base_name} -t 16 --discard-unmapped &&\
        coverm filter -b ${BAM1} -o ${BAM2} --proper-pairs-only --threads 16
        samtools fastq ${BAM2} \
        -1 07_Reassemble/${SAMPLE}/${base_name}/${base_name}_1.fastq.gz \
        -2 07_Reassemble/${SAMPLE}/${base_name}/${base_name}_2.fastq.gz \
        -n -@ 16 &&\
        spades.py -1 07_Reassemble/${SAMPLE}/${base_name}/${base_name}_1.fastq.gz \
        -2 07_Reassemble/${SAMPLE}/${base_name}/${base_name}_2.fastq.gz \
        -o 08_REASSEMBLIES/${base_name} -m 125 -t 16 ;\
        seqkit seq -m 1000 08_REASSEMBLIES/${base_name}/scaffolds.fasta -o 09_FINALGENOMES/${base_name}_scaffolds.fasta.gz ;\
        gzip 08_REASSEMBLIES/*/*.fasta ; gzip 08_REASSEMBLIES/*/*/*.fasta"

        sbatch -A ${ACCOUNT} -J ${N} -c 16 --mem=125G -o 99_logs/${N}_o.txt -e 99_logs/${N}_e.txt \
        --export=ALL --mail-type=ALL --mail-user=${EMAIL} -t 24:00:00 --wrap="${CMD}"

    done
done
```
