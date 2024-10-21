#!/bin/bash

# RNASeq data analysis pipeline
sample_name=$1
input_file1=$2
input_file2=$3
outdir=$4
bowtie2_index=$5

# fastqc: ~
output_dir=${outdir}/fastqc
mkdir -p ${output_dir}
fastqc ${input_file1} ${input_file2} -o ${output_dir}

# cutadapt: ~
output_dir=${outdir}/cutadapt
mkdir -p ${output_dir}
cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG -q 20 -e 0.1 -O 1 --minimum-length=75 --max-n=0.1 -o ${output_dir}/trimmed_${sample_name}_R1.fastq.gz -p ${output_dir}/trimmed_${sample_name}_R2.fastq.gz ${input_file1} ${input_file2}

# bowtie2: ~
input_dir=${outdir}/cutadapt
output_dir=${outdir}/bowtie2
mkdir -p ${output_dir}
bowtie2 -x ${bowtie2_index} -1 ${input_file1} -2 ${input_file2} -S ${output_dir}/${sample_name}.sam -p 10

# samtools:
input_dir=${outdir}/bowtie2
output_dir=${outdir}/bowtie2
samtools view -bS ${input_dir}/${sample_name}.sam | samtools sort -o ${output_dir}/${sample_name}_sorted.bam

# stringtie:
input_dir=${outdir}/bowtie2
output_dir=${outdir}/stringtie
mkdir -p ${output_dir}
stringtie ${input_dir}/${sample_name}_sorted.bam -p 8 -G /lustre1/g/aos_shihuang/data/RNA_Seq_Data/data_analysis/database/sequence.gff3 -o ${output_dir}/${sample_name}_assembled_transcripts.gtf
Rscript /lustre1/g/aos_shihuang/data/RNA_Seq_Data/data_analysis/htseq_reference_id_scripts/filter_assembled_transcripts_gtf_file.R -i ${output_dir}/${sample_name}_assembled_transcripts.gtf -o  ${output_dir}/filtered_${sample_name}_assembled_transcripts.gtf

# FPKM
input_dir=${outdir}/stringtie
output_dir=${outdir}/FPKM
mkdir -p ${output_dir}
python /lustre1/g/aos_shihuang/data/RNA_Seq_Data/data_analysis/FPKM_scripts/FPKM.py ${input_dir}/filtered_${sample_name}_assembled_transcripts.gtf ${output_dir}/${sample_name}_FPKM.txt

# TPM
input_dir=${outdir}/stringtie
output_dir=${outdir}/TPM
mkdir -p ${output_dir}
python /lustre1/g/aos_shihuang/ek/SaeRSRNASeq/TPM.py ${input_dir}/filtered_${sample_name}_assembled_transcripts.gtf ${output_dir}/${sample_name}_TPM.txt
