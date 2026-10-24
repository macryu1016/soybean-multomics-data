#!/bin/bash
# ============================================================
# Soybean salt stress multi-omics integration pipeline
# Author: Wenjiying
# Date: 2026-10
# ============================================================

set -euo pipefail

# ---------- 全局路径 ----------
BASE=/data/wenjiying/bean
CHIP=$BASE/chip
QC=$CHIP/QC
OUT=$BASE/tables
mkdir -p $QC $OUT

# ============================================================
# PART 1: ChIP-seq
# ============================================================
echo "=== PART 1: ChIP-seq ==="

# 1.1 比对（已在 run_chip_v2.pbs 中完成）
# fastp → bowtie2 → samtools sort/index
# 输出：$CHIP/alignment/*.bam

# 1.2 MACS2 peak calling（BAMPE + input）
cd $CHIP
mkdir -p peaks
for rep in 1 2 3; do
  case $rep in
    1) SRR=SRR9622336 ;;
    2) SRR=SRR9622337 ;;
    3) SRR=SRR9622338 ;;
  esac
  macs2 callpeak -t alignment/${SRR}.sorted.bam \
    -c alignment/Control_Input.bam \
    -f BAMPE -g 1.0e9 -n Control_IP_rep${rep} -q 0.05 \
    --outdir peaks
done

for rep in 1 2 3; do
  case $rep in
    1) SRR=SRR9622341 ;;
    2) SRR=SRR9622342 ;;
    3) SRR=SRR9622343 ;;
  esac
  macs2 callpeak -t alignment/${SRR}.sorted.bam \
    -c alignment/Salt_Input.bam \
    -f BAMPE -g 1.0e9 -n Salt_IP_rep${rep} -q 0.05 \
    --outdir peaks
done

# 1.3 QC：library size
for bam in alignment/*.bam; do
  echo -e "$(basename $bam)\t$(samtools view -c $bam)"
done > $QC/library_size.txt

# 1.4 QC：mapping rate
for log in alignment/*_bowtie2.log; do
  sample=$(basename $log _bowtie2.log)
  rate=$(grep "overall alignment rate" $log | awk '{print $1}')
  echo -e "$sample\t$rate"
done > $QC/mapping_rate.txt

# 1.5 QC：duplication rate（samtools markdup）
mkdir -p $QC/markdup
for bam in alignment/*.bam; do
  sample=$(basename $bam .bam)
  samtools sort -n -@ 4 -o $QC/markdup/${sample}.namesort.bam $bam
  samtools fixmate -m $QC/markdup/${sample}.namesort.bam $QC/markdup/${sample}.fixmate.bam
  samtools sort -@ 4 -o $QC/markdup/${sample}.positionsort.bam $QC/markdup/${sample}.fixmate.bam
  samtools markdup $QC/markdup/${sample}.positionsort.bam $QC/markdup/${sample}.markdup.bam
  total=$(samtools view -c $QC/markdup/${sample}.markdup.bam)
  dup=$(samtools view -f 1024 -c $QC/markdup/${sample}.markdup.bam)
  echo -e "$sample\t$total\t$dup\t$(echo "scale=4; $dup/$total" | bc)"
done > $QC/duplication.txt

# 1.6 QC：FRiP
for rep in 1 2 3; do
  for cond in Control Salt; do
    sample="${cond}_IP_rep${rep}"
    case "${cond}_${rep}" in
      Control_1) bam="alignment/SRR9622336.sorted.bam" ;;
      Control_2) bam="alignment/SRR9622337.sorted.bam" ;;
      Control_3) bam="alignment/SRR9622338.sorted.bam" ;;
      Salt_1) bam="alignment/SRR9622341.sorted.bam" ;;
      Salt_2) bam="alignment/SRR9622342.sorted.bam" ;;
      Salt_3) bam="alignment/SRR9622343.sorted.bam" ;;
    esac
    peak="peaks/${sample}_peaks.narrowPeak"
    bedtools bamtobed -i $bam > $QC/${sample}.bed
    total=$(wc -l < $QC/${sample}.bed)
    inpeak=$(bedtools intersect -a $QC/${sample}.bed -b $peak | wc -l)
    echo -e "${sample}\t${total}\t${inpeak}\t$(echo "scale=4; $inpeak/$total" | bc)"
  done
done > $QC/FRiP.txt

# 1.7 降采样 + MACS2 重跑
mkdir -p alignment_downsampled peaks_downsampled
for bam in alignment/*.bam; do
  sample=$(basename $bam .bam)
  total=$(samtools view -c $bam)
  frac=$(awk -v m=7024620 -v t=$total 'BEGIN{printf "%.4f", m/t}')
  keep=${frac#0.}
  samtools view -b -s 42.$keep -o alignment_downsampled/${sample}.bam $bam
  samtools index alignment_downsampled/${sample}.bam
done

for rep in 1 2 3; do
  case $rep in
    1) SRR=SRR9622336 ;;
    2) SRR=SRR9622337 ;;
    3) SRR=SRR9622338 ;;
  esac
  macs2 callpeak -t alignment_downsampled/${SRR}.sorted.bam \
    -c alignment_downsampled/Control_Input.bam \
    -f BAMPE -g 1.0e9 -n Control_IP_rep${rep}_downsampled -q 0.05 \
    --outdir peaks_downsampled
done

for rep in 1 2 3; do
  case $rep in
    1) SRR=SRR9622341 ;;
    2) SRR=SRR9622342 ;;
    3) SRR=SRR9622343 ;;
  esac
  macs2 callpeak -t alignment_downsampled/${SRR}.sorted.bam \
    -c alignment_downsampled/Salt_Input.bam \
    -f BAMPE -g 1.0e9 -n Salt_IP_rep${rep}_downsampled -q 0.05 \
    --outdir peaks_downsampled
done

for f in peaks_downsampled/*.narrowPeak; do
  echo -e "$(basename $f)\t$(wc -l < $f)"
done > $QC/downsampled_peak_counts.txt

# 1.8 Peak 数 CI + DiffBind FDR（R）
Rscript $BASE/scripts/chip_peak_CI_FDR.R

# ============================================================
# PART 2: m6A DART-seq
# ============================================================
echo "=== PART 2: m6A DART-seq ==="

Rscript $BASE/scripts/m6a_FDR_QC.R

# ============================================================
# PART 3: RNA–protein correlation
# ============================================================
echo "=== PART 3: RNA–protein ==="

Rscript $BASE/scripts/rna_protein_correlation.R

# ============================================================
# PART 4: Phosphoproteomics
# ============================================================
echo "=== PART 4: Phosphoproteomics ==="

Rscript $BASE/scripts/phospho_permutation.R

# ============================================================
# PART 5: TF activity
# ============================================================
echo "=== PART 5: TF activity ==="

Rscript $BASE/scripts/tf_activity_CI_FDR.R

# ============================================================
# PART 6: Metabolomics
# ============================================================
echo "=== PART 6: Metabolomics ==="

Rscript $BASE/scripts/metabolite_FDR_CI.R

# ============================================================
# PART 7: Pathway joint enrichment
# ============================================================
echo "=== PART 7: Pathway ==="

Rscript $BASE/scripts/pathway_joint_Fisher.R

echo "=== Pipeline complete ==="
