### This script has been written my Miriam Llorian
### To analyse RNA Seq data for PE in order to run rMATS
### Need to first check the length distribution after trimming
### and keep as many reads as possible of unique length

### need to use and modify the parameters CROP= and MINLEN= in Trimmomatic
#### IMPORTANT NEED TO CHANGE CROP and MINLEN in Trimmomatic command

# module load R/3.5.1-foss-2016b-BABS
# R

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2),
        sep = "", collapse = " ")
}

#### IMPORTANT THINGS TO MODIFY IN THIS SCRIPT: 

#fastqDir,
#PI, 
#project, 
#scientist, 
#library type,
#star index, 
#REFERENCE_HOME, 
#REF,
#GTF,
#refflat for CollectRnaSeqMetrics
#riboint
#strand specificity for picard "STRAND_SPECIFICITY=" 
# RNA-SeQC: needs to be told if it's single end
### change rRNA interval list, refflat and  accordingly (see below)

species <- "homo_sapiens"
assembly <- "GRCh38"
## to keep consistency with previous analysiss
release.no <- "89"
read.length <- "100"

### Reference genome ###
REFERENCE_HOME <- paste0("/camp/svc/reference/Genomics/babs/", species, "/ensembl/",assembly,"/release-",release.no,"/")
REF <- paste0(REFERENCE_HOME,"genome_idx/bowtie2/", simpleCap(species),".",assembly,".dna_sm.primary_assembly")  ### problem with this link
ref <- paste0(REFERENCE_HOME,"genome/", simpleCap(species),".",assembly,".dna_sm.primary_assembly")  ### problem with this link
## check to change in CollectAlignmentMetrics and rnaseqc commads when corrected
GTF <- paste0(REFERENCE_HOME,"gtf/", simpleCap(species),".",assembly, ".", release.no, ".gtf")


###  VERY IMPORTANT set library type:
rsem.library.type  <- "--forward-prob 0"

#### star indexes change depending on the length of the reads, change accordingly
star.idx      <- paste0(REFERENCE_HOME,"genome_idx/rsem/star/", read.length,"bp/genome")

### strand specificity For strand-specific library prep. For unpaired reads, use FIRST_READ_TRANSCRIPTION_STRAND if the reads are expected to be on the transcription strand.
s <- "SECOND_READ_TRANSCRIPTION_STRAND" ## change in CollectRnaSeqMetrics

### Ribosomal interval file
riboInt <- paste0(REFERENCE_HOME, "gtf/", simpleCap(species),".",assembly, ".", release.no, ".rRNA.interval_list")

### reflat
refflat <- paste0(REFERENCE_HOME, "gtf/", simpleCap(species),".",assembly, ".", release.no, ".refflat")

###gtf_rnaseqc
gtf_rnaseqc <-  paste0(REFERENCE_HOME, "gtf/", simpleCap(species),".",assembly, ".", release.no,".rnaseqc.gtf")

####
rRNA <-  paste0(REFERENCE_HOME, "gtf/", simpleCap(species),".",assembly, ".", release.no, ".rRNA.list")

bed <-  paste0(REFERENCE_HOME, "gtf/", simpleCap(species),".",assembly, ".", release.no, ".bed")
#######################################
### Create a directory structure  #####
#######################################
projectDir <- "/camp/stp/babs/working/sopenam/projects/"
### Change the following entries for each project accordingly
PI <- "tybulewiczv/"
scientist <- "eva.lana-elola/"
project <- "12_RNASeq_HumanEmbryonicHearts_DownSyndrome/"

### Read sample names from fastqfiles ####
###  data is in three sequencing directories 

SequencingDir <- "/camp/stp/babs/inputs/sequencing/fastq/"

run <- c("201001_K00371_0407_BHJ5YHBBXY")

ProjNum <- "RN20085/"

SeqDir <- paste0( SequencingDir, run[1], "/fastq/", ProjNum)

samples <- dir(SeqDir, pattern="_L00[1-9]_R[1]_001.fastq.gz")
samples <- gsub("_R[1]_001.fastq.gz","", samples) # replace suffix with empty space
samples


shortsamplename <- unique(gsub("_S[0-9]*_L00[1-9]", "", samples))
shortsamplename <- unique(shortsamplename)
shortsamplename

#### Define directories in the project tree
dataDir <- paste0(projectDir, PI,scientist,project,"data/")
fastqDir <- paste0(dataDir, "fastq/")
trimDir <- paste0(dataDir,"trimmed/")
trim.fastQCDir <- paste0(dataDir,"trimmed/fastQC/")
alignDir <- paste0(dataDir,"align/")
rsemDir <-  paste0(dataDir,"align/rsem/")
rsem.tmp <- paste0(dataDir,"align/rsem/tmp/")
logsDir <- paste(dataDir,"logs/",sep="")
scriptsDir <-  paste(projectDir, PI,scientist,project,"commands/",sep="")
rnaseqcDir <- paste(rsemDir,"rnaseqc/",sep="")
resultsDir <- paste0(projectDir, PI,scientist,project,"results/")

DESeq.Dir <- paste0(resultsDir,"DESeq/")
GSEA.Dir <- paste0( resultsDir, "GSEA/")

## Create a series of directories (if they don't exist) ##

if (!file.exists( paste0 ( projectDir, PI))) { dir.create(paste0 ( projectDir, PI))}
if (!file.exists( paste0 ( projectDir, PI,scientist))) { dir.create(paste0 ( projectDir, PI,scientist))}
if (!file.exists( paste0 ( projectDir, PI,scientist, project))) { dir.create(paste0 ( projectDir, PI,scientist,project))}

if (!file.exists( paste0 ( projectDir, PI,scientist, project, "/commands/"))) { dir.create(paste0 ( projectDir, PI,scientist,project,"/commands/"))}
if (!file.exists( paste0 ( projectDir, PI,scientist, project, "/data/"))) { dir.create(paste0 ( projectDir, PI,scientist,project,"/data/"))}
if (!file.exists( paste0 ( projectDir, PI,scientist, project, "/docs/"))) { dir.create(paste0 ( projectDir, PI,scientist,project,"/docs/"))}
if (!file.exists( resultsDir)) { dir.create( resultsDir)}
if (!file.exists( paste0 ( projectDir, PI,scientist, project, "/scripts/"))) { dir.create(paste0 ( projectDir, PI,scientist,project,"/scripts/"))}
if (!file.exists( paste0 ( projectDir, PI,scientist, project, "/test/"))) { dir.create(paste0 ( projectDir, PI,scientist,project,"/test/"))}

### Results from DESeq, etc are now stored in results directory

if (!file.exists( paste0(resultsDir, "picard/"))) { dir.create(paste0(resultsDir, "picard/"))}
if (!file.exists( paste0(resultsDir, "rnaseqc/"))) { dir.create(paste0(resultsDir, "rnaseqc/"))}
if (!file.exists( paste0 ( resultsDir, "DESeq/"))) { dir.create( paste0 ( resultsDir, "DESeq/"))}
if (!file.exists( paste0 ( resultsDir, "GSEA/"))) {dir.create ( paste0 ( resultsDir, "GSEA/"))}

### directories under data

if (!file.exists( paste0 ( dataDir, "fastq/"))) { dir.create( paste0 ( dataDir, "fastq/"))}
if (!file.exists( paste0 ( dataDir, "trimmed/"))) { dir.create( paste0 ( dataDir, "trimmed/"))}
if (!file.exists( paste0 ( dataDir, "trimmed/fastQC/"))) { dir.create( paste0 ( dataDir, "trimmed/fastQC/"))}
if (!file.exists( paste0 ( dataDir, "align/"))) { dir.create( paste0 ( dataDir, "align/"))}
if (!file.exists( paste0 ( dataDir, "align/rsem/"))) { dir.create( paste0 ( dataDir, "align/rsem/"))}
if (!file.exists( paste0 ( dataDir, "align/rsem/tmp/"))) { dir.create( paste0 ( dataDir, "align/rsem/tmp/"))}
if (!file.exists( paste0 ( dataDir, "align/rsem/rnaseqc/"))) { dir.create( paste0 ( dataDir, "align/rsem/rnaseqc/"))}
if (!file.exists( paste0 ( dataDir, "fastqscreen/"))) { dir.create( paste0 ( dataDir, "fastqscreen/"))}
if (!file.exists( paste0 ( dataDir, "logs/"))) { dir.create( paste0 ( dataDir, "logs/"))}

for (sample in shortsamplename){
  if (!file.exists( paste0( rsemDir, sample ))){ dir.create(paste0( rsemDir, sample ))}
  if (!file.exists( paste0( rnaseqcDir, sample ))){ dir.create(paste0( rnaseqcDir, sample ))}
  }

### Tools required

FastQC <- "/camp/apps/eb/software/FastQC/0.11.5-Java-1.8.0_92/fastqc"
Trimmomatic <- "/camp/apps/eb/software/Trimmomatic/0.36-Java-1.7.0_80/trimmomatic-0.36.jar"
flagstat <- "/camp/apps/eb/software/SAMtools/1.3.1-foss-2016b/bin/samtools flagstat"
rsem <- "/camp/apps/eb/software/RSEM/1.2.31-foss-2016b/bin/rsem-calculate-expression"
LibraryComplexity <- "/camp/apps/eb/software/picard/2.1.1-Java-1.8.0_112/picard.jar EstimateLibraryComplexity"
sort <- "/camp/apps/eb/software/SAMtools/1.3.1-foss-2016b/bin/samtools sort"
Index <- "/camp/apps/eb/software/SAMtools/1.3.1-foss-2016b/bin/samtools index"

Addreadgroup <- "/camp/apps/eb/software/picard/2.1.1-Java-1.8.0_112/picard.jar AddOrReplaceReadGroups"
BuildBamIndex <- "/camp/apps/eb/software/picard/2.1.1-Java-1.8.0_112/picard.jar BuildBamIndex"
MarkDuplicates <- "/camp/apps/eb/software/picard/2.1.1-Java-1.8.0_112/picard.jar MarkDuplicates"
ReorderSam <- "/camp/apps/eb/software/picard/2.1.1-Java-1.8.0_112/picard.jar ReorderSam"
SortSam <- "/camp/apps/eb/software/picard/2.1.1-Java-1.8.0_112/picard.jar SortSam"
rnaseqc <- "/camp/apps/eb/software/RNA-SeQC/1.1.8-Java-1.7.0_80/RNA-SeQC_v1.1.8.jar"

#### modules
mod.purge <- "module purge"
module.fastqc <- "module load FastQC/0.11.5-Java-1.8.0_92"
module.trimmom <- "module load Trimmomatic/0.36-Java-1.7.0_80"

module.rsem <- "module load RSEM/1.2.31-foss-2016b"
module.perl <- "module load Perl/5.24.0-foss-2016b"
module.star <-"module load STAR/2.5.2a-foss-2016b"
module.picard <- "module load picard/2.1.1-Java-1.8.0_112"
module.R <- "module load R/3.3.1-foss-2016b-bioc-3.3-libX11-1.6.3"
module.samtools <- "module load SAMtools/1.3.1-foss-2016b"

module.rnaseqc <- "module load RNA-SeQC/1.1.8-Java-1.7.0_80"
module.rseqc <- "module load RSeQC/2.6.4-foss-2016b-Python-2.7.12-R-3.3.1"

#i=1
###############################
#### READ ALINGMENT AND QC  ###
###############################
 for (sample in shortsamplename) {
   
   ###############################
   ### Concatenate fastq files ###
   ###############################
   ####### samples are now split between different lanes so they need to be join together ###
   # zcat.R1 <- paste( "zcat" , 
   #                   paste0(SeqDir, samples[grep(paste0(sample, "_"), samples)], "_R1_001.fastq.gz", collapse=" "),
   #                   paste0(SeqDir2, samples2[grep(paste0(sample, "_"), samples2)],"_R1_001.fastq.gz", collapse=" " ),
   #                   ">", paste0(fastqDir,sample, "_R1_001.fastq"))
   # 
   FQ1 <- paste0(SeqDir, samples[grep(paste0(sample, "_"), samples)], "_R1_001.fastq.gz")
   

   #########################
   #### Fastqscreen ########
   #########################
   
   module.fastqscreen <- "module load FastQ_Screen/0.12.1-foss-2018a-Perl-5.26.1"
   mod.bowtie<- "module load Bowtie2/2.2.9-foss-2016b"
   fastqscreenDir <- paste0("/camp/stp/babs/working/sopenam/projects/",PI,scientist,project,"/data/fastqscreen/")
   
   contam.screen1.cmd <- paste("fastq_screen", 
                               "--conf /camp/stp/babs/working/sopenam/Genomes/fastq_screen.conf",
                               "--aligner bowtie2",
                               "--outdir", fastqscreenDir,
                               FQ1,
                               ">",paste0(logsDir,sample,".00.fastqscreen1 2>&1"),
                               sep=" ")

   ####################
   ### Trimmomatic  ###
   ####################
   ### Make sure to use the correct Phred value and ADAPTORS, modify script below accordingly to data
   
   F1p <- paste(trimDir,sample, "_R1.paired.fastq", sep="")
   F1u <- paste(trimDir,sample, "_R1.unpaired.fastq", sep="")
 
   Trimming.cmd <- paste("java -Xmx8g -jar", Trimmomatic, "SE","-phred33",
                         FQ1, F1p, 
                         "ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36", 
                         sep=" ",
                         ">",paste(logsDir,sample,".02.trimming 2>&1",sep=''))
   
   ### FastQC after Trimming
   fastQCafterTrimming.cmd <- paste(FastQC,F1p, "--outdir",trim.fastQCDir, sep=" ", 
                                    ">",paste(logsDir,sample,".03.fQC_trimmed 2>&1",sep=''))
   
   
   ### Aling with RSEM/STAR
   ## RSEM paramaters
   ## library type this needs to be changed depending if it's unstranded (0.5), first strand(0), second strand(1)
   if (!file.exists(rsem.tmp)){ dir.create(rsem.tmp)}
   if (!file.exists( paste0( rsemDir, sample ))){ dir.create(paste0( rsemDir, sample ))}
   
   rsem.cmd <- paste(rsem, "--temporary-folder",
                     paste0(rsem.tmp,sample,"/"),
                     #"--paired-end",
                     "--star",
                     "--output-genome-bam", 
                     "--star-output-genome-bam",
                     "--time --num-threads 6",
                     #"--calc-ci --ci-memory 10240 --estimate-rspd",
                    #"--star-gzipped-read-file",
                     rsem.library.type,
                     F1p,
                     star.idx,
                     paste0(rsemDir,sample,"/",sample), 
                     sep=" ",
                     ">",paste0(logsDir,sample,".04.rsem 2>&1"))
   
   
   ###sort and index STAR.genome.bam
   ###  sorting
   I <- paste0(rsemDir,sample,"/",sample,".STAR.genome.bam")
   O <- paste0(rsemDir,sample,"/",sample,".STAR.genome.sorted.bam")
   sort.cmd <- paste(sort,"-o",O, I, sep = " ",
                     ">",paste(logsDir,sample,".05.sort 2>&1",sep=''))
   

   ### Index
   index.cmd <- paste(Index,O, sep = " ")
   
   ### Add read groups and sort STAR.genome.bam
   
   I <- paste("I=",rsemDir,sample,"/",sample, ".STAR.genome.sorted.bam",sep="")
   O <- paste("O=",rsemDir,sample,"/",sample, ".STAR.genome.rg.sorted.bam",sep="")
   S <- "SORT_ORDER=coordinate"
   RGID <- paste0("RGID=",sample)
   RGLB <- paste0("RGLB=",sample)
   RGPL <- "RGPL=illumina"
   RGSM <- paste0("RGSM=",sample)
   
   Addreadgroup.STAR.cmd <- paste("java -Xmx8g -jar", Addreadgroup, I, O, S,
                                  RGID, RGLB, RGPL, "RGPU=unit1", RGSM, "MAX_RECORDS_IN_RAM=2000000 VALIDATION_STRINGENCY=LENIENT", sep=" ",
                                  ">",paste(logsDir,sample,".06.AddRG_sort 2>&1",sep=''))
   
   
 
   ### Mark duplicates
   
   Input <-  paste0("INPUT=",paste0(rsemDir,sample,"/",sample, ".STAR.genome.rg.sorted.bam"))
   Output <- paste0("OUTPUT=",rsemDir,sample,"/",sample,".STAR.genome.rg.sorted.dedup.bam")
   Metrics <- paste0("METRICS_FILE=",resultsDir,"picard/",sample,".dedup.metrics.txt")
   MarkDuplicates.cmd <- paste("java -Xmx8g -jar", 
                               MarkDuplicates, 
                               Input, 
                               Output, 
                               Metrics,
                               "MAX_RECORDS_IN_RAM=2000000", 
                               "VALIDATION_STRINGENCY=LENIENT", 
                               "REMOVE_DUPLICATES=false",
                               "ASSUME_SORTED=true",
                               sep=" ",
                               ">",paste0(logsDir,sample,".07.markdups 2>&1"))
   
    
   ### Index .STAR.genome.rg.sorted.dedup.bam
   index.STAR.genome.rg.sorted.dedup.cmd <- paste(Index,paste0(rsemDir,sample,"/",sample,".STAR.genome.rg.sorted.dedup.bam"), sep = " ")
   
   ### Estimate library Complexity
   In <- paste("INPUT=",rsemDir,sample,"/",sample,".STAR.genome.rg.sorted.dedup.bam",sep="")
   Out <- paste("OUTPUT=",resultsDir,"picard/",sample,".complexity",sep="")
   LibraryComplexity.cmd <- paste("java -Xmx8g -jar", 
                                  LibraryComplexity, 
                                  "VALIDATION_STRINGENCY=SILENT", 
                                  In, 
                                  Out,
                                  sep=" ",
                                  ">",paste(logsDir,sample,".08.LibraryComplexity 2>&1",sep=''))
   
   ###CollectAlignmentSummaryMetrics
   R <- paste0("R=",ref,".fa")
   I <- paste0("I=",rsemDir,sample,"/",sample,".STAR.genome.rg.sorted.dedup.bam")
   O <- paste0("O=",resultsDir,"picard/",sample,".CollectAlignmentSummaryMetrics")
   CollectAlingSumMetrics.cmd <- paste("java -Xmx8g -jar ${EBROOTPICARD}/picard.jar CollectAlignmentSummaryMetrics",
                                       R,
                                       I,
                                       O,
                                       ">",paste0(logsDir,sample,".09.CollectAlignmentSummaryMetrics 2>&1"),
                                       sep=" ")
   
   ### CollectRnaSeqMetrics
   I <- paste0("I=",rsemDir,sample,"/",sample,".STAR.genome.rg.sorted.dedup.bam")
   O <- paste0("O=",resultsDir,"picard/",sample,".output.RNA_Metrics")
   
   CollectRnaSeqMetrics.cmd <- paste("java -Xmx8g -jar /camp/apps/eb/software/picard/2.1.1-Java-1.8.0_112/picard.jar CollectRnaSeqMetrics",
                                     I,
                                     O,
                                     paste0("REF_FLAT=",refflat),
                                     paste0("STRAND_SPECIFICITY=",s),
                                     paste0("RIBOSOMAL_INTERVALS=",riboInt),
                                     ">",paste0(logsDir,sample,".10.CollectRnaSeqMetrics 2>&1"),
                                     sep=" ")
   
   ### Flagstat
   bam <- paste(rsemDir,sample,"/",sample, ".STAR.genome.rg.sorted.dedup.bam",sep="")
   flagstat.cmd <- paste(flagstat, bam, ">", paste(rsemDir,sample,"/",sample,"_flagstat.txt",sep=""))
   
   ### RNA-SeQC using STAR.genome.bam
   if (!file.exists(paste0(rnaseqcDir, sample))) { dir.create( paste0(rnaseqcDir, sample))}
   ss <- paste0("'",sample,"|",rsemDir,sample,"/",sample,".STAR.genome.rg.sorted.dedup.bam","|","bam_file'")
   t <-  gtf_rnaseqc
   rnaseqc_rsemDir <- paste0(rsemDir,"rnaseqc/")
   
   o <- paste0(rnaseqc_rsemDir,sample)
   RNASeQC.cmd <- paste("java -Xmx8g -jar",
                        rnaseqc,
                        "-singleEnd",
                        "-o", o, 
                        "-r",paste0(ref,".fa"),
                        "-rRNA", rRNA,
                        "-s",ss, 
                        "-t", t,
                        "TMP_DIR=", rsem.tmp,
                        "-gatkFlags '-S SILENT -U ALLOW_SEQ_DICT_INCOMPATIBILITY'",
                        ">",paste(logsDir,sample,".16b.RNASeQC 2>&1",sep=''),
                        sep=" ")
   
   ### RSeQC to check Junction Saturation
   junct.sat <- paste("junction_saturation.py",
                      "-i", bam,
                      "-r", bed,
                      "-o",  paste0(rsemDir,sample,"/",sample,"_JunctionSaturation_metrics.txt"),
                      ">",paste0(logsDir,sample,".17.RSeQC.JuncSat 2>&1"),
                      sep=" ")
   
   junct.anno <- paste("junction_annotation.py",
                       "-i", bam,
                       "-r", bed,
                       "-o",  paste0(rsemDir,sample,"/",sample,"_JunctionAnno_metrics.txt"),
                       ">",paste0(logsDir,sample,".18.RSeQC.JuncAnno 2>&1"),
                       sep=" ")
   
   ######## Write out to file
   
   line1<- "#!/bin/sh"
   
   sink(paste(scriptsDir,sample,".RSEM.STAR.SE.sh",sep=""))
   cat("#!/bin/sh");cat("\n");
   cat(mod.purge); cat("\n");
   cat(module.fastqc);cat("\n");

  # cat(zcat.R1); cat("\n");

   cat(module.fastqscreen); cat("\n");
   cat(contam.screen1.cmd); cat("\n");

   cat(module.trimmom);cat("\n");
   cat(Trimming.cmd);cat("\n");
   cat('sleep 10'); cat("\n");
   cat(fastQCafterTrimming.cmd);cat("\n");
   cat('sleep 10');cat("\n");
   
   cat(module.rsem);cat("\n");
   cat(module.perl);cat("\n");
   cat(module.star);cat("\n");
   cat(rsem.cmd);cat("\n");
   
   cat(module.samtools);cat("\n");
   cat(sort.cmd);cat("\n");
   cat('sleep 10');cat("\n");
   cat(index.cmd);cat("\n\n");
   
   cat(module.picard);cat("\n");
   cat(module.R);cat("\n");
   cat(Addreadgroup.STAR.cmd);cat("\n");
   cat(MarkDuplicates.cmd);cat("\n\n");
   cat('sleep 10');cat("\n");
   cat(index.STAR.genome.rg.sorted.dedup.cmd);cat("\n\n");
   
   cat(LibraryComplexity.cmd);cat("\n\n");
   cat('sleep 10');cat("\n");
   cat(CollectAlingSumMetrics.cmd);cat("\n\n");
   cat('sleep 10');cat("\n");
   cat(CollectRnaSeqMetrics.cmd);cat("\n\n");
   cat('sleep 10');cat("\n\n");
   
   cat(module.samtools);cat("\n");
   cat(flagstat.cmd);cat("\n");
   cat('sleep 10');cat("\n\n");
   
   cat(module.rnaseqc);cat("\n");
   cat(RNASeQC.cmd);cat("\n\n");
   
   cat(module.rseqc);cat("\n");
   cat(junct.sat);cat("\n\n");
   cat(junct.anno);cat("\n\n");
   
   
   ### clean up, remove intermediate large files
   ## remove Fastq files
  # cat( "rm ", paste0( fastqDir, sample, "_R1.fastq"));cat("\n");
  # cat( "rm ", paste0( fastqDir, sample, "_R2.fastq"));cat("\n");
   
   ## remove Trimmed files
   cat("rm ",paste0( trimDir, sample, "_R1.paired.fastq"));cat ("\n");

   cat("rm ",paste0( trimDir, sample, "_R1.unpaired.fastq"));cat("\n");

   ## remove files from RSEM/STAR
   
   cat("rm ",paste0( rsemDir, sample, "/", sample, ".transcript.bam"));cat("\n");
   cat("rm ",paste0( rsemDir, sample, "/", sample, ".genome.bam"));cat("\n");
   
   cat("rm ",paste0( rsemDir, sample, "/", sample, ".STAR.genome.bam"));cat("\n");
   cat("rm ",paste0( rsemDir, sample, "/", sample, ".STAR.genome.rg.sorted.bam"));cat("\n");
   cat("rm ",paste0( rsemDir, sample, "/", sample, ".STAR.genome.sorted.bam"));cat("\n");
   cat("rm ",paste0( rsemDir, sample, "/", sample, ".STAR.genome.sorted.bam.bai"));cat("\n");
   #cat("rm ",paste0( rsemDir, sample, "/", sample, ".STAR.genome.rg.sorted.dedup.bam"));cat("\n");
   #cat("rm ",paste0( rsemDir, sample, "/", sample, ".STAR.genome.rg.sorted.dedup.bam.bai"));cat("\n");
   
   sink();
   
   
   #create a bash script to submit jobs
   sink(paste0( projectDir, PI, scientist, project, "scripts/submit_Rsemjobs.sh"), append = TRUE)
   cat("#!/bin/sh" ); cat ("\n\n");  
   cat("sbatch -c 6 --mem-per-cpu=6G --time=08:00:00"); cat(" "); 
   cat ( paste(scriptsDir,sample,".RSEM.STAR.SE.sh",sep=""));cat("\n"); 
   sink()
#   i=i+1
 }

system(paste( "bash", paste0( projectDir, PI, scientist, project, "scripts/submit_Rsemjobs.sh"), sep=" "))

### NOT SURE THIS WORKS ###
#submit jobs to cluster
#jobs <- system(paste( "bash", paste0( projectDir, PI, scientist, project, "scripts/submit_Rsemjobs.sh"), sep=" "))
               
##get the job id for the batch
#jobid <- strsplit( jobs, " ")[[1]][4]



#### Once the rnaseqc table has been created 
#### create a multiqc report

sink( paste0 ( scriptsDir,"Multiqc.sh" ))
cat("#!/bin/sh" ); cat ("\n\n");  
cat("module purge"); cat("\n\n");
cat( "#Create a summary file in results"); cat("\n");

cat( paste("touch",paste0(resultsDir, "rnaseqc/metrics.tsv") )); cat("\n");

cat("samples='", shortsamplename, "'", sep="\n"); cat("\n");

cat( "for sample in $samples"); cat("\n");
cat("do"); cat("\n");
cat(paste( "awk" , '-F \"\\t\" ', "'NR==2 { print  $1,$18,$32,$40,$21,$43,$5,$7,$29,$19,$31,$24,$10,$44,$3,$4,$47} ' OFS=\"\\t\"",
           paste0(rnaseqcDir, "$sample/metrics.tsv"), " |  tee -a", paste0(resultsDir, "rnaseqc/metrics.tsv")) ); cat("\n");


cat("done"); cat("\n");

cat("### add header to top of file ");cat("\n");
cat(paste( "awk" , 
           '-F \"\\t\" ',
           "'NR==1 { print  $1,$18,$32,$40,$21,$43,$5,$7,$29,$19,$31,$24,$10,$44,$3,$4,$47} ' OFS=\"\\t\"",
           paste0(rnaseqcDir, "$sample/metrics.tsv"), 
           " |  cat -", 
           paste0(resultsDir, "rnaseqc/metrics.tsv"), "> temp && mv temp ",paste0(resultsDir, "rnaseqc/metrics.tsv") )); cat("\n\n");

cat("### Run multiQC ####");cat("\n\n");

cat("module load MultiQC/1.6-Python-2.7.15-foss-2018a"); cat("\n");
cat("cp /camp/stp/babs/working/sopenam/scripts_template/multiqc16/multiqc_config.yaml");cat(" ");cat(resultsDir);cat("\n");
cat("cd "); cat(resultsDir); cat("\n");
cat("multiqc");cat(" ");
cat(trim.fastQCDir);cat(" "); 
cat(paste0(resultsDir, "rnaseqc/"));cat(" ");
cat(paste0 (resultsDir, "picard/"));cat ("\n");
sink()

sink( paste0(scriptsDir, "submit_multiqc.sh"))
cat("#!/bin/sh" ); cat ("\n\n");  
cat("sbatch -c 6 --time=01:00:00 --mem-per-cpu=6G"); cat(" "); 
#cat( paste0( "-d afterok:", jobid )); cat(" "); 
cat ( paste0( scriptsDir,"Multiqc.sh"));cat("\n"); 
sink()

system( paste( "bash", paste0( scriptsDir, "submit_multiqc.sh"), sep=" " ))

