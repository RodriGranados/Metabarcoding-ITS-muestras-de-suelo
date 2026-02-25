# =========================================
# 1. Importación de lecturas crudas
# =========================================

# Importa archivos FASTQ paired-end usando un archivo manifest.
# Se especifica que las secuencias están en formato Phred33.
# Genera un artefacto .qza con las lecturas demultiplexadas.

qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path manifest.csv \
--output-path demux.qza \
--input-format PairedEndFastqManifestPhred33V2


# =========================================
# 2. Eliminación de primers
# =========================================

# Remueve secuencias de primers forward y reverse usando cutadapt.
# -p-front-f y -p-front-r: secuencias de los primers ITS
# --p-minimum-length 100: descarta lecturas <100 pb tras el recorte
# --verbose: muestra detalles del proceso

qiime cutadapt trim-paired \
--i-demultiplexed-sequences ../importar_qiime2/demux.qza \
--p-front-f GCATCGATGAAGAACGCAGC \
--p-front-r TCCTCCGCTTATTGATATGC \
--p-minimum-length 100 \
--o-trimmed-sequences trimmed.qza \
--verbose


# =========================================
# 3. Evaluación de calidad post-trimming
# =========================================

# Genera visualización interactiva para inspeccionar calidad
# y decidir parámetros de truncamiento para DADA2.

qiime demux summarize \
--i-data trimmed.qza \
--o-visualization trimmed.qzv


# =========================================
# 4. Denoising y generación de ASVs con DADA2
# =========================================

# Modelo de error y corrección con DADA2:
# -p-trunc-len: longitud de truncamiento basada en calidad
# -p-max-ee: máximo número de errores esperados permitidos
# -p-n-threads 32: uso paralelo
# Salidas:
#   table.qza → tabla de abundancias (ASVs)
#   rep-seqs.qza → secuencias representativas
#   denoising-stats.qza → estadísticas de filtrado

qiime dada2 denoise-paired \
--i-demultiplexed-seqs ../limpieza_primer_ITS_filtro_longitud/trimmed.qza \
--p-trim-left-f 0 \
--p-trim-left-r 0 \
--p-trunc-len-f 279 \
--p-trunc-len-r 279 \
--p-max-ee-f 2.8 \
--p-max-ee-r 2.8 \
--p-n-threads 32 \
--o-table table.qza \
--o-representative-sequences rep-seqs.qza \
--o-denoising-stats denoising-stats.qza


# =========================================
# 5. Revisión de estadísticas de denoising
# =========================================

# Permite evaluar pérdida de lecturas en cada paso:
# filtrado, denoising, merge y chimera removal.

qiime metadata tabulate \
--m-input-file denoising-stats.qza \
--o-visualization denoising-stats.qzv


# =========================================
# 6. Resumen de tabla de abundancias
# =========================================

# Visualiza número de ASVs por muestra,
# profundidad de secuenciación y distribución.

qiime feature-table summarize \
--i-table table.qza \
--o-visualization table.qzv


# Visualización de secuencias representativas

qiime feature-table tabulate-seqs \
--i-data rep-seqs.qza \
--o-visualization rep-seqs.qzv


# =========================================
# 7. Identificación y evaluación de singletons
# =========================================

# Extrae features con frecuencia total = 1
# (ASVs presentes solo una vez en todo el dataset)

qiime feature-table filter-features \
--i-table table.qza \
--p-min-frequency 1 \
--p-max-frequency 1 \
--o-filtered-table singletons.qza


# (Repetición desde otra ruta — verificar redundancia)

qiime feature-table filter-features \
--i-table ../3.denoising/denoising/table.qza \
--p-min-frequency 1 \
--p-max-frequency 1 \
--o-filtered-table singletons.qza


# Resumen de singletons

qiime feature-table summarize \
--i-table singletons.qza \
--o-visualization singletons.qzv


# =========================================
# 8. Eliminación de singletons
# =========================================

# Retiene únicamente ASVs con frecuencia ≥2
# Reduce ruido y posibles artefactos de secuenciación.

qiime feature-table filter-features \
--i-table ../3.denoising/denoising/table.qza \
--p-min-frequency 2 \
--o-filtered-table table-no-singletons.qza


qiime feature-table summarize \
--i-table table-no-singletons.qza \
--o-visualization table-no-singletons.qzv


# =========================================
# 9. Asignación taxonómica
# =========================================

# Clasificación taxonómica usando Naive Bayes pre-entrenado
# basado en base de datos UNITE para ITS.
# Método sklearn implementado en QIIME 2.

qiime feature-classifier classify-sklearn \
--i-classifier unite_ver10_dynamic_19.02.2025-Q2-2024.10.qza \
--i-reads ../3.denoising/denoising/rep-seqs.qza \
--o-classification taxonomy.qza


# =========================================
# 10. Visualización de composición taxonómica
# =========================================

# Genera gráfico de barras por nivel taxonómico
# integrando metadata de muestras.

qiime taxa barplot \
--i-table ../4.eliminar_singlentons/table-no-singletons.qza \
--i-taxonomy taxonomy.qza \
--m-metadata-file ../3.denoising/importar_qiime2/metadata_ITS2.tsv \
--o-visualization taxa-bar-plots.qzv
