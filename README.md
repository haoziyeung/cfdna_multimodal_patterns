# cfdna_multimodal_patterns


**Discovery multi-modal signals captured by fragmentation patterns of cell-free DNA**

-------------------
- **Primary data processing of low-pass whole genome sequencing(LD-WGS) data** :
``` bash
#quality control, adapter trimming, quality filtering and per-read quality pruning
fastp --cut_by_quality3 -l 25 --correction -w 2 -V -s 24 -z 1 -j Sample.json -h Sample.html --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTC --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTG -o Sample.Clean.R1.fq.gz -O Sample.Clean.R2.fq.gz -i Sample.Raw.R1.fastq.gz -I Sample.Raw.R2.fastq.gz

#speedseq align converts paired-end FASTQ sequences to a duplicate-marked, sorted, indexed BAM file that can be processed with other downstream modules.
speedseq align -t 30 -M 13 -o Sample -R "@RG:\tID:Sample\tSM:Sample\tPL:ILLUMINA\tCN:BerryOncology\tLB:Sample\tPM:NovaSeq" /disk/bundle/hs37d5/hs37d5.fa Sample.Clean.R1.gz Sample.Clean.R2.fq.gz
```
- **Multi-modal fragmentomic features derived from LD-WGS assay** :
``` bash
#Nucleosome Footprint: Inspired by a former study (Peter Ulz et al,. Nature Genetics 2016), we conducted NF(Nucleosome) score to distinguish between expressed and silent genes on the basis of plasma read coverage characteristics
samtools view -H Sample.bam | grep SQ | sed 's/SN://;s/LN://' | perl -lane'print "$F[1]\t1\t$F[2]"' | grep -v chrM > genome.file

wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz

zcat refGene.txt.gz | perl -lane'next if $F[2]=~/_/;$tss=$F[3] eq "+" ? $F[4] : $F[5];print join "\t",($F[2],$tss-2000,$tss-200,"$F[-4]:$F[1]:TSS1");print join "\t",($F[2],$tss-200,$tss+100,"$F[-4]:$F[1]:NDR");print join "\t",($F[2],$tss+100,$tss+2000,"$F[-4]:$F[1]:TSS2")' | grep -v chrMT | sort -k1V,1 -k2n,2 > hg19_RefSeq_TSS_NDR.bed

bedtools coverage -sorted -mean -g genome.file -a hg19_RefSeq_TSS_NDR.bed -b Sample.bam |sort -k1V,1 -k4 | gzip -c > Sample.gz

perl TSSNDR.pl Sample.gz > Sample.NF.score

#Fragment: Fragment score was built to calculate the length of insertion fragment and ratio of short/long fragment in different regions. 
python fragment.feature.py Sample.bam Sample HCC > Sample.Fragment.score

#Motif(a.k.a fragment endpoints): Fragment endpoint coordinates were extracted from BAM files with the SAMtools API. We identified 256 different types of 4-mer 5’end motif and calculated their percentages without considering chromosome Y and unidentifiable bases.
samtools view $1 chr{1..22} | awk -F '\t' '{if($9>0 && $9<170 && $5>15){print $3"\t"$4-1"\t"$4+3}}' | bedtools getfasta -fi /disk/bundle/hs37d5/hs37d5.fa -bed - -tab | cut -f 2 | sed 's/[a-z]/\u&/g' |grep -v 'N' | perl -pe 's/^\s+//g' | perl -pe 's/ +/\t/g' | pigz > Sample.motif.gz

perl motif.freq.pl Sample.motif.gz motif.256.list > Sample.Motif.score
```
- **Development of a HIFI model for early detection of hepatocellular carcinoma based on longitudinal  cell-free DNA signatures** :   `Sample.NF.score`, `Sample.Fragment.score` and `Sample.Motif.score` can be used as raw data frame(vector) with other downstream modules. The Wilcoxon rank-sum test was used to compare two datasets, HCC vs non-HCC (cirrhosis & healthy controls). Least Absolute Shrinkage and Selection Operator (LASSO) methods were applied to further reduce the number of markers in training set.  A machine learning method (support vector machine) was implemented for individual genomic feature-based model construction. Standardize features by removing the mean and scaling to unit variance. A machine learning method (support vector machine) was implemented for individual genomic feature-based model construction, basing on three parameters: (1) C: Penalty coefficient; (2) Kernel function; and (3) Gamma. For training dataset, 5-fold cross-validation was employed to figure out the best combination of the parameters. Three genomic features (Fragment, Motif and NF, as described above) were ensembled to establish a `HIFI` model to distinguish cancer individuals from non-cancer individuals.  The documentation includes more detailed https://github.com/scikit-learn/scikit-learn  about making a standard machine learning pipeline with scikit-learn.

**NOTE**: All scripts and binaries are provided as is, without any warrenty and for use at your own risk. This is not the release of a software package. We are only providing this information and code in addition to a description of methods for making it easier to reproduce our analyses. We are not providing any support for these scripts.
