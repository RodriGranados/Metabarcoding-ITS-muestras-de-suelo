#!/bin/bash

# Directorios
RAW_DIR=2.reads_limpios
OUT_DIR=3.reads_limpios
REPORT_DIR=3.reporte_limpieza

# Loop sobre los R1
for R1 in ${RAW_DIR}/*_R1.fastq.gz
do
    # Obtener R2 correspondiente
    R2=${R1/_R1.fastq.gz/_R2.fastq.gz}

    # Verificar que exista el R2
    if [[ ! -f "$R2" ]]; then
        echo "⚠️  No se encontró el par para $R1"
        continue
    fi

    # Nombre base de la muestra
    SAMPLE=$(basename $R1 _R1.fastq.gz)

    echo "Procesando muestra: ${SAMPLE}"

    fastp \
        -i ${R1} \
        -I ${R2} \
        -o ${OUT_DIR}/${SAMPLE}_clean3_R1.fastq.gz \
        -O ${OUT_DIR}/${SAMPLE}_clean3_R2.fastq.gz \
        --adapter_sequence TGGAATTCTCGGGTGCCAAGG \
        --adapter_sequence_r2 TGGAATTCTCGGGTGCCAAGG \
        --trim_poly_g \
        --detect_adapter_for_pe \
        -h ${REPORT_DIR}/${SAMPLE}_fastp.html \
        -j ${REPORT_DIR}/${SAMPLE}_fastp.json \
        --thread 32
done

echo "✅ Limpieza con fastp finalizada"
