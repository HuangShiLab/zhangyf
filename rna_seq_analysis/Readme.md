# Introduction

## Central dogma of molecular biology
<img width="964" alt="image" src="https://github.com/HuangShiLab/zhangyf/assets/170502144/78d81876-3272-4032-b3f8-efdea4e42a8b">

We harbor a collection of DNA, RNA, protein, and metabolites from both human and microbial cells. DNA can provide instructions for proteins. During the expression of a protein-coding gene, the central dogma of molecular biology states that genetic information flows only in one direction, from DNA, to RNA, to protein, ultimately affecting the metabolite levels. 

<img width="950" alt="image" src="https://github.com/HuangShiLab/zhangyf/assets/170502144/aaf64e44-5a01-442b-80dd-f19b27cca63b">

While performing microbiome data analysis, we collect microbiome samples and extract DNA for sequencing to uncover the potential capabilities of microbes by asking "Who are they?" and "What can they do?". However, the actual gene expression is reflected by transcription. Transcriptomics analysis involves isolation and purification of RNA from samples, followed by conversion from RNA to cDNA with reverse transcriptase. The cDNA is loaded in a sequencer and differential expressions are analysed.

## Transcriptome data analysis workflow
<img width="965" alt="image" src="https://github.com/HuangShiLab/zhangyf/assets/170502144/0aa4f270-751f-453e-846a-6cf8015b5d64">
The main workflow of transcriptome data analysis contains five steps.

### 1. Quanlity control
1.1 Sequencing data quality analysis
Software: FastQC v0.11.9
Inputdata: Raw data (fastq.gz file)
Commands:
```
output_dir=${outdir}/fastqc
mkdir -p ${output_dir}
fastqc -o "$output_dir" "$input_file1" "$input_file2"
```

1.2 Data Filtering
Software: Cutadapt version 4.8
Method description:
(1) remove the adapter sequences
(2) remove the 5' or 3' end bases that contains N’s or of quality values below 20
(3) remove reads that are less than 75 bp after trimming
Inputdata: Raw data (fastq.gz file)
Commands:
```
output_dir=${outdir}/cutadapt
mkdir -p ${output_dir}
cutadapt -q 20 -e 0.1 -O 1 --minimum-length=75 --max-n=0.1 \
  -a AGATCGGAAGAG -A AGATCGGAAGAG \
  -o ${output_dir}/trimmed_${sample_name}_R1.fastq.gz -p ${output_dir}/trimmed_${sample_name}_R2.fastq.gz \
  ${input_file1} ${input_file2}
```

Notes: Adapter sequences:  
•	Illumina Universal Adapter — AGATCGGAAGAG  
•	Illumina Small RNA 3' Adapter — TGGAATTCTCGG  
•	Illumina Small RNA 5' Adapter — GATCGTCGGACT  
•	Nextera Transposase Sequence — CTGTCTCTTATA  

2. 




