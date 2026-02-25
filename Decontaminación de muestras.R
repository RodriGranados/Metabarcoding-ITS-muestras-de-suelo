#### DECONTAMINACION DE CONTROLES ####

## Cargar librerias

library(phyloseq)
library(decontam) ## Identificacion estadistica de contaminantes
library(tidyverse)
library(dplyr)
library(readr)
library(tibble)


## ARREGLAR TAX
TAX <- TAX[, !colnames(TAX) %in% c("SH", "Confidence")]

add_prefix <- function(x, prefix) {
  ifelse(is.na(x), paste0(prefix), paste0(prefix, x))
}

TAX$Kingdom <- add_prefix(TAX$Kingdom, "k__")
TAX$Phylum  <- add_prefix(TAX$Phylum,  "p__")
TAX$Class   <- add_prefix(TAX$Class,   "c__")
TAX$Order   <- add_prefix(TAX$Order,   "o__")
TAX$Family  <- add_prefix(TAX$Family,  "f__")
TAX$Genus   <- add_prefix(TAX$Genus,   "g__")
TAX$Species <- add_prefix(TAX$Species, "s__")


## Leer metadata de decontaminación
metadata_decontam <- read_tsv("metadata_decontaminacion.txt")

## Usar la columna SampleID como rownames (nombres de fila)
metadata_decontam <- metadata_decontam %>% column_to_rownames("SampleID")

## Construccion del objeto phyloseq
ps <- phyloseq(
  otu_table(as.matrix(ASV), taxa_are_rows = TRUE),  ## Tabla de abundancias
  tax_table(as.matrix(TAX)),                        ## Taxonomia
  sample_data(metadata_decontam)                    ## Metadata de decontaminacion
)

## Creacion de la variable logica para decontaminacion
sample_data(ps)$is.neg <- sample_data(ps)$Tipo == "Control"

## Identificacion de decontaminantes con decotam
## Clasifica como contaminantes aquellos ASVs mas frecuentes en los
## controles que en las muestras reales
contam <- isContaminant(
  ps,
  method = "prevalence",
  neg = "is.neg",
  threshold = 0.1
)

## Resumen de ASVs contaminantes
## FALSE --> ASVs eliminados (controles)
## TRUE --> ASVs retenidos (muestras reales)
table(contam$contaminant)

## Dataset limpio (muestras decontaminadas) para siguientes analisis
ps.clean <- prune_taxa(!contam$contaminant, ps)

## Dataset solo con contaminantes para inspeccion
ps_contam <- prune_taxa(contam$contaminant, ps)

## Visualizacion de los contaminantes
plot_bar(ps_contam, fill = "Phylum")

tax_table(ps)[contam$contaminant, ] %>% head()



# Extraer taxonomía de ASVs contaminantes
tax_contam <- as.data.frame(
  tax_table(ps)[contam$contaminant, ]
)

# Agregar el ID del ASV como columna
tax_contam$ASV <- rownames(tax_contam)

# Reordenar columnas (ASV primero)
tax_contam <- tax_contam %>%
  select(ASV, everything())

## Inspeccion taxonomica
head(tax_contam)

table(tax_contam$Phylum)
table(tax_contam$Family)

## Añadir estadisticos de decotam
tax_contam$prev_stat <- contam$prev[contam$contaminant]  ## Estadistico de prevalencia
tax_contam$p_prev <- contam$p.prev[contam$contaminant]   ## p-valor del metodo

## Inspeccion de la tabla con ASVs, taxonomica y estadisticas 
## de los contaminantes removidos
View(tax_contam)


# Prevalencia real en controles
prev_control <- apply(
  otu_table(ps)[, sample_data(ps)$is.neg == TRUE] > 0,
  1,
  mean
)

# Prevalencia real en muestras
prev_muestra <- apply(
  otu_table(ps)[, sample_data(ps)$is.neg == FALSE] > 0,
  1,
  mean
)

# Añadir a tabla de ASVs contaminantes
tax_contam$prev_control <- prev_control[rownames(tax_contam)]
tax_contam$prev_muestra <- prev_muestra[rownames(tax_contam)]

View(tax_contam)

## Graficos
library(ggplot2)

## Dsitribucion taxonomica de contaminantes
ggplot(tax_contam, aes(x = Phylum)) +
  geom_bar() +
  coord_flip() +
  labs(
    x = "Phylum",
    y = "Número de ASVs contaminantes",
    title = "Distribución taxonómica de ASVs contaminantes"
  ) +
  theme_minimal()

## Grafico de significancia estadistica (p_prev)
ggplot(tax_contam, aes(x = p_prev)) +
  geom_histogram(bins = 30) +
  labs(
    x = "p-value (prevalence test)",
    y = "Número de ASVs",
    title = "Distribución de significancia estadística (Decontam)"
  ) +
  theme_minimal()




tax_contam$log_p <- -log10(tax_contam$p_prev)

ggplot(tax_contam,
       aes(x = Phylum,
           y = log_p)) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
  geom_hline(yintercept = -log10(0.1),
             linetype = "dashed",
             color = "red") +
  coord_flip() +
  labs(
    x = "Phylum",
    y = expression(-log[10](p[prev])),
    title = "Significancia de ASVs contaminantes por taxonomía (Phylum)",
    subtitle = "Línea roja: umbral Decontam (p = 0.1)"
  ) +
  theme_minimal()

ggplot(tax_contam,
       aes(x = Genus,
           y = log_p,
           color = Phylum)) +
  geom_jitter(width = 0.25, size = 2, alpha = 0.8) +
  geom_hline(yintercept = -log10(0.1),
             linetype = "dashed") +
  coord_flip() +
  labs(
    x = "Género",
    y = expression(-log[10](p[prev])),
    title = "Significancia de ASVs contaminantes por género"
  ) +
  theme_minimal()


ggplot(tax_contam,
       aes(x = Phylum,
           y = log_p)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +
  geom_jitter(width = 0.2, size = 1.8, alpha = 0.7) +
  geom_hline(yintercept = -log10(0.1),
             linetype = "dashed",
             color = "red") +
  coord_flip() +
  labs(
    x = "Phylum",
    y = expression(-log[10](p[prev])),
    title = "Distribución de significancia de ASVs contaminantes por Phylum"
  ) +
  theme_minimal()



ggplot(tax_contam,
       aes(x = Phylum,
           y = prev_muestra,
           fill = Phylum)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.7) +
  coord_flip() +
  labs(
    x = "Phylum",
    y = "Prevalencia en muestras",
    title = "Presencia residual de ASVs contaminantes en muestras"
  ) +
  theme_minimal() +
  guides(fill = "none")



library(tidyr)

tax_long <- tax_contam %>%
  select(ASV, Phylum, prev_control, prev_muestra, p_prev) %>%
  pivot_longer(cols = c(prev_control, prev_muestra),
               names_to = "Tipo",
               values_to = "Prevalencia")

ggplot(tax_long,
       aes(x = Tipo,
           y = ASV,
           fill = Prevalencia)) +
  geom_tile() +
  facet_wrap(~ Phylum, scales = "free_y") +
  scale_fill_viridis_c() +
  labs(
    x = "",
    y = "ASVs contaminantes",
    title = "Prevalencia de ASVs contaminantes en controles vs muestras"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_blank())



df_long$Stage <- factor(df_long$Stage, levels = c("Before", "After"))

ggplot(df_long,
       aes(x = Richness,
           fill = Stage)) +
  geom_density(alpha = 0.5) +
  labs(
    title = "Distribución de riqueza antes y después de decontam",
    x = "Riqueza observada",
    y = "Densidad"
  ) +
  theme_minimal()




## COMPARACION LECTURAS ANTES Y DESPUES

# Lecturas totales
total_before <- sum(sample_sums(ps))
total_after  <- sum(sample_sums(ps.clean))

lost_reads <- total_before - total_after
lost_pct   <- (lost_reads / total_before) * 100

data.frame(
  Total_antes = total_before,
  Total_despues = total_after,
  Lecturas_perdidas = lost_reads,
  Porcentaje_perdido = round(lost_pct, 2)
)

# Abundancias por muestra
ab_before <- sample_sums(ps)
ab_after  <- sample_sums(ps.clean)

loss_table <- data.frame(
  SampleID = names(ab_before),
  Lecturas_antes = ab_before,
  Lecturas_despues = ab_after,
  Lecturas_perdidas = ab_before - ab_after,
  Porcentaje_perdido = round(((ab_before - ab_after) / ab_before) * 100, 2)
)

head(loss_table)


metadata_1 <- data.frame(sample_data(ps))

loss_table <- cbind(loss_table, metadata_1[loss_table$SampleID, ])


library(ggplot2)

ggplot(loss_table, aes(x = reorder(SampleID, Porcentaje_perdido),
                       y = Porcentaje_perdido)) +
  geom_col(fill = "firebrick") +
  coord_flip() +
  labs(
    x = "Muestra",
    y = "% de lecturas eliminadas",
    title = "Pérdida de lecturas tras Decontam"
  ) +
  theme_bw()


ggplot(loss_table, aes(x = Tipo, y = Porcentaje_perdido)) +
  geom_boxplot(fill = "steelblue", alpha = 0.7) +
  geom_jitter(width = 0.1, size = 2) +
  labs(
    x = "",
    y = "% de lecturas eliminadas",
    title = "Impacto de Decontam en controles vs muestras"
  ) +
  theme_classic()
