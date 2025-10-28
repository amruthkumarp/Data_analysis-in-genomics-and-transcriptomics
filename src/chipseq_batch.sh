#!/bin/bash
# =============================================
#   Automated Multi-sample ChIP-seq Pipeline

#   Tools: BWA, SAMtools, MACS2, HOMER, BEDTools
# =============================================

set -e  # stop if any command fails

# User Configurable Parameters

THREADS=4
GENOMESIZE=4.6e6
REF="E.coli_BW25113.fasta"
GFF="E.coli_BW25113_annotation.gff"
BED="E.coli_BW25113_annotation.bed"
TREAT_DIR="treatment"
CTRL_DIR="control"
OUTDIR="chipseq_results"

# Create folders
mkdir -p ${OUTDIR}/{logs,merged}


#  Index reference

echo "=== Checking reference index ==="
if [ ! -f "${REF}.bwt" ]; then
    echo "Indexing reference genome..."
    bwa index ${REF}
else
    echo "Reference genome already indexed."
fi


# Loop through treatment FASTQs

echo "=== Starting ChIP-seq analysis ==="

for CHIP in ${TREAT_DIR}/*.fastq.gz; do
    SAMPLE=$(basename "$CHIP" .fastq.gz)
    
    # Identify base name (remove CHR, R1, etc. if present)
    BASE=$(echo "$SAMPLE" | sed 's/_CHR1_R1//; s/_R1//; s/_rep[0-9]*//')

    # Try to find matching control
    CTRL=$(ls ${CTRL_DIR}/*${BASE}*fastq.gz 2>/dev/null | head -n 1)

    if [ -z "$CTRL" ]; then
        echo "  No control found for $SAMPLE, skipping..."
        continue
    fi

    echo ">>> Processing sample: ${BASE}"
    echo "    Treatment: $CHIP"
    echo "    Control:   $CTRL"

    SAMPLE_OUT="${OUTDIR}/${BASE}"
    mkdir -p ${SAMPLE_OUT}/{alignments,peaks,annotation,motifs}

    
    #  Align Treatment Reads
    
    echo "Aligning treatment sample..."
    bwa aln -q 20 -t ${THREADS} ${REF} ${CHIP} > ${SAMPLE_OUT}/alignments/${BASE}_chip.sai
    bwa samse ${REF} ${SAMPLE_OUT}/alignments/${BASE}_chip.sai ${CHIP} > ${SAMPLE_OUT}/alignments/${BASE}_chip.sam
    samtools view -bSq 20 ${SAMPLE_OUT}/alignments/${BASE}_chip.sam -o ${SAMPLE_OUT}/alignments/${BASE}_chip.bam -@${THREADS}
    samtools sort ${SAMPLE_OUT}/alignments/${BASE}_chip.bam -o ${SAMPLE_OUT}/alignments/${BASE}_chip_sorted.bam -@${THREADS}
    samtools index ${SAMPLE_OUT}/alignments/${BASE}_chip_sorted.bam -@${THREADS}

    
    # Align Control Reads
  
    echo "Aligning control sample..."
    bwa aln -q 20 -t ${THREADS} ${REF} ${CTRL} > ${SAMPLE_OUT}/alignments/${BASE}_ctrl.sai
    bwa samse ${REF} ${SAMPLE_OUT}/alignments/${BASE}_ctrl.sai ${CTRL} > ${SAMPLE_OUT}/alignments/${BASE}_ctrl.sam
    samtools view -bSq 20 ${SAMPLE_OUT}/alignments/${BASE}_ctrl.sam -o ${SAMPLE_OUT}/alignments/${BASE}_ctrl.bam -@${THREADS}
    samtools sort ${SAMPLE_OUT}/alignments/${BASE}_ctrl.bam -o ${SAMPLE_OUT}/alignments/${BASE}_ctrl_sorted.bam -@${THREADS}
    samtools index ${SAMPLE_OUT}/alignments/${BASE}_ctrl_sorted.bam -@${THREADS}

    
    #  Peak Calling
    
    echo "Running MACS3 peak calling..."
    macs3 callpeak \
        -t ${SAMPLE_OUT}/alignments/${BASE}_chip_sorted.bam \
        -c ${SAMPLE_OUT}/alignments/${BASE}_ctrl_sorted.bam \
        -f BAM -g ${GENOMESIZE} -n ${BASE} \
        -B -p 0.001 --nomodel \
        --outdir ${SAMPLE_OUT}/peaks

    
    # Peak Annotation
    
    echo "Annotating peaks..."
    annotatePeaks.pl ${SAMPLE_OUT}/peaks/${BASE}_peaks.narrowPeak ${REF} -gff ${GFF} > ${SAMPLE_OUT}/annotation/${BASE}_annotated.tsv

 
    # Extract FASTA sequences
    
    echo "Extracting peak sequences..."
    bedtools getfasta -fi ${REF} -bed ${SAMPLE_OUT}/peaks/${BASE}_peaks.narrowPeak -fo ${SAMPLE_OUT}/annotation/${BASE}_peaks.fa

   
    #  Motif Analysis

    echo "Running motif analysis..."
    findMotifs.pl ${SAMPLE_OUT}/annotation/${BASE}_peaks.fa fasta ${SAMPLE_OUT}/motifs/ -len 8,10,12 -norevopp

    echo "Completed: ${BASE}"
done

echo " All samples processed successfully!"

