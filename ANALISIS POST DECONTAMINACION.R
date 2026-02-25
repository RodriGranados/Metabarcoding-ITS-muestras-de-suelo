### LLAMAR LIBRERIAS
library(microeco)
library(readxl)
library(ggnested)
library(magrittr)
library(rgexf)

### CARGAR TABLA DE ABUNDANCIA, TAXONOMIA Y METADATA
ASV_post_decom <- as.data.frame(read_excel("feature table post decom.xlsx"))
rownames(ASV_post_decom) <- ASV_post_decom$ASV
ASV_post_decom$ASV <- NULL


TAX <- as.data.frame(read_excel("TAXONOMY.xlsx"))
rownames(TAX) <- TAX$ASV
TAX$ASV <- NULL

METADATA <- as.data.frame(read_excel("metadata_ingles.xlsx"))
rownames(METADATA) <- METADATA$...1
METADATA$...1 <- NULL

### CREAR OBJETO MICROTABLE
mt <- microtable$new(
  otu_table = ASV_post_decom,
  sample_table = METADATA,
  tax_table = TAX,
  phylo_tree = NULL,
  rep_fasta = NULL,
  auto_tidy = TRUE
)


## CALCULAR ABUNDANCIA - Graficos de composicion (barras)
mt$cal_abund()

mt$taxa_abund$Kingdom
mt$taxa_abund$Phylum
mt$taxa_abund$Genus

t1 <- trans_abund$new(dataset = mt, taxrank = "Phylum", ntaxa = 15)
t1$plot_bar(others_color = "grey70", facet = "m.s.n.m", xtext_keep = TRUE, legend_text_italic = TRUE)

t2 <- trans_abund$new(dataset = mt, taxrank = "Family", ntaxa = 15)
t2$plot_bar(others_color = "grey70", facet = "m.s.n.m", xtext_keep = TRUE, legend_text_italic = TRUE)

t3 <- trans_abund$new(dataset = mt, taxrank = "Genus", ntaxa = 15)
t3$plot_bar(others_color = "grey70", facet = "Texture", xtext_keep = TRUE, legend_text_italic = TRUE)
t3$plot_bar(others_color = "grey70", facet = c("Texture", "m.s.n.m"), xtext_keep = FALSE, legend_text_italic = FALSE, barwidth = 0.9)

## Graficos alluvial

t1$plot_bar(others_color = "grey70", facet = "m.s.n.m", xtext_keep = FALSE, legend_text_italic = TRUE, use_alluvium = TRUE)
t2$plot_bar(others_color = "grey70", facet = "m.s.n.m", xtext_keep = FALSE, legend_text_italic = TRUE, use_alluvium = TRUE)
t3$plot_bar(others_color = "grey70", facet = "m.s.n.m", xtext_keep = FALSE, legend_text_italic = TRUE, use_alluvium = TRUE)

## Graficos boxplot

t2 <- trans_abund$new(dataset = mt, taxrank = "Family", ntaxa = 15)
t2$plot_box(xtext_angle=30)
t2$plot_box(xtext_angle=30, group = "Localities")

t3 <- trans_abund$new(dataset = mt, taxrank = "Species", ntaxa = 15)
t3$plot_box(xtext_angle=30)
t3$plot_box(xtext_angle=30, group = "Localities")
t3$plot_box(xtext_angle=30, group = "m.s.n.m")


## Grafico de heatmap

t2$plot_heatmap(xtext_angle=45)
t2$plot_heatmap(xtext_angle=45, facet = "Localities")
t2$plot_heatmap(xtext_angle=45, facet = "m.s.n.m")

t3$plot_heatmap(xtext_angle=45)
t3$plot_heatmap(xtext_angle=45, facet = "Localities")
t3$plot_heatmap(xtext_angle=45, facet = "m.s.n.m")


## Grafico anidado

t4 <- trans_abund$new(dataset = mt, taxrank = "Class", ntaxa = 15,
                      high_level="Phylum", delete_taxonomy_prefix=FALSE)
t4$plot_bar(xtext_angle=60, ggnested=TRUE, barwidth = 1, facet = "m.s.n.m")


## Barplot + Clusterpplot

t5 <- trans_abund$new(dataset = mt, taxrank = "Class", ntaxa = 15,
                      high_level="Phylum", delete_taxonomy_prefix=FALSE, groupmean = "m.s.n.m")
t5$plot_bar(xtext_angle = 45)

t5$plot_bar(coord_flip = T)
t5$plot_bar(coord_flip = T, clustering_plot = T)


### NORMALIZACION (RAREFACCION)

sort(colSums(mt$otu_table))

tn1 <- trans_norm$new(dataset=mt)

mt_rarefy <- tn1$norm(method="rarefy", sample.size = 10236)


library(vegan)

otu_mat <- as.data.frame(mt$otu_table)

# vegan necesita muestras en filas
rarecurve(t(otu_mat),
          step = 500,
          sample = min(colSums(otu_mat)),
          cex = 0.6)



min_depth <- min(colSums(otu_mat))

rarecurve(t(otu_mat),
          step = 500,
          sample = min_depth,
          label = FALSE)



library(vegan)

otu_mat <- as.data.frame(mt$otu_table)

# Profundidad mínima
min_depth <- min(colSums(otu_mat))

# Número de muestras
n_samples <- ncol(otu_mat)

# Generar colores distintos
colores <- rainbow(n_samples)

rarecurve(t(otu_mat),
          step = 500,
          sample = min_depth,
          col = colores,
          label = FALSE)



library(RColorBrewer)

colores <- colorRampPalette(brewer.pal(8, "Dark2"))(n_samples)

rarecurve(t(otu_mat),
          step = 500,
          sample = min_depth,
          col = colores,
          label = FALSE)



library(vegan)

otu_mat <- as.data.frame(mt$otu_table)

min_depth <- min(colSums(otu_mat))
n_samples <- ncol(otu_mat)

colores <- rainbow(n_samples)

rarecurve(t(otu_mat),
          step = 500,
          sample = min_depth,
          col = colores,
          label = FALSE,
          xlab = "Número de lecturas",
          ylab = "Riqueza observada (ASVs)")

rarecurve(t(otu_mat),
          step = 500,
          sample = min_depth,
          col = colores,
          label = FALSE,
          xlab = "Número de lecturas",
          ylab = "Riqueza observada (ASVs)",
          cex.lab = 1.2,
          cex.axis = 1.1)

## DIVERSIDAD ALFA

mt_rarefy$cal_alphadiv()
mt_rarefy$alpha_diversity

t6 <-trans_alpha$new(dataset = mt_rarefy, group="Texture")

t6$cal_diff(method="KW")
t6$res_diff

t6$cal_diff(method="wilcox")
t6$res_diff

t6$cal_diff(method="t.test")
t6$res_diff

## Comparacion con 2 categorias

t7 <- trans_alpha$new(dataset = mt_rarefy, group = "Localities")

t7$cal_diff(method="anova", formula="Localities+Texture") ## Mas diferencias significativas (**)

t7$res_diff


## Comparacion de medidas de alpha diversidad segun GROUP
t8 <- trans_alpha$new(dataset = mt_rarefy, group = "Texture")

t8$cal_diff(method="anova")
t8$res_diff

t8$plot_alpha(measure = "Shannon", add = "jitter")

t8$plot_alpha(measure = "Simpson")

t8$plot_alpha(measure = "InvSimpson")

t8$plot_alpha(measure = "Pielou")


## DIVERSIDAD BETA

mt_rarefy$cal_betadiv()

mt_rarefy$beta_diversity$bray %>% View()

t9 <- trans_beta$new(dataset = mt_rarefy, group = "m.s.n.m", measure = "bray")
t9$cal_ordination(method="PCoA")
t9$plot_ordination(plot_color = "m.s.n.m", plot_type = c("point", "ellipse"))


t10 <- trans_beta$new(dataset = mt_rarefy, group = "Texture", measure = "bray")
t10$cal_ordination(method="PCoA")
t10$plot_ordination(plot_color = "Texture", plot_type = c("point", "ellipse"))


t11 <- trans_beta$new(dataset = mt_rarefy, group = "Localities", measure = "bray")
t11$cal_ordination(method="PCoA")
t11$plot_ordination(plot_color = "Localities", plot_type = c("point", "ellipse"))


## Distancias intragrupos

t10$cal_group_distance(within_group = T)
t10$cal_group_distance_diff(method="t.test")

t10$plot_group_distance()

## Distancias entre grupos

t10$cal_group_distance(within_group = F)
t10$cal_group_distance_diff(method="t.test")

t10$plot_group_distance()


## Clustering

t10$plot_clustering()

t10$plot_clustering(group = "Texture")


### TRANS DIFF

t12 <- trans_diff$new(dataset = mt, method = "lefse", group="m.s.n.m",
                     alpha=0.05, lefse_subgruop=NULL)

t12$res_diff %>% View()

t12$plot_diff_bar(threshold = 3)

t12$plot_diff_abund()


### CLADOGRAMAS

t12$plot_diff_cladogram(use_taxa_num = 100,
                       use_feature_num = 25,
                       clade_label_level=5,
                       group_order = c("700-900","500-700","300-500"))


### USAR FUNGAL TRAITS
mt$tidy_dataset()

# Crear un objeto trans_func a partir del microtable
# Este objeto se usa para mapear ASVs/OTUs a funciones ecológicas
t13 <- trans_func$new(mt)

# Asignar funciones ecológicas usando la base de datos FungalTraits
t13$cal_func(fungi_database = "FungalTraits")

# Calcular la redundancia funcional ponderada por abundancia
# Esto considera no solo la presencia del rasgo,
# sino también cuán abundantes son los taxa que lo poseen
t13$cal_func_FR(abundance_weighted = TRUE)


# Clonar el objeto mt para evitar modificar el original
tmp_mt <- clone(mt)

# Reemplazar manualmente la tabla de abundancia por la tabla funcional
tmp_mt$taxa_abund$func <- as.data.frame(t(t13$res_func_FR), check.names = FALSE)

# Crear un objeto trans_diff para pruebas estadísticas
t13 <- trans_diff$new(dataset = tmp_mt, method = "anova", group = "m.s.n.m", taxa_level = "func")
# Graficar diferencias de abundancia funcional
t13$plot_diff_abund(add_sig = T, simplify_names = FALSE) + ggplot2::ylab("Relative abundance (%)")


# Crear un objeto trans_network para construir la red
# cal_cor = "WGCNA" → método robusto para detectar módulos
# taxa_level = "OTU" → la red se construye a nivel de ASV/OTU
# filter_thres → elimina taxa extremadamente raros
# cor_method = "spearman" → correlación no paramétrica
t13 <- trans_network$new(dataset = mt, cal_cor = "WGCNA", taxa_level = "OTU", filter_thres = 0.0001, cor_method = "spearman")

# Calcular la red de correlación
# COR_p_thres → umbral de significancia
# COR_cut → fuerza mínima de correlación
t13$cal_network(COR_p_thres = 0.05, COR_cut = 0.6)

# Identificar módulos (clusters) dentro de la red
t13$cal_module()

# Convertir la información de módulos a un microtable
# Cada módulo se trata como una "comunidad"
meco_module <- t13$trans_comm(use_col = "module")

# Crear un nuevo objeto trans_func basado en los módulos
t14 <- trans_func$new(meco_module)

# Asignar funciones ecológicas a los módulos
t14$cal_func(fungi_database = "FungalTraits")

# Calcular redundancia funcional NO ponderada por abundancia
t14$cal_func_FR(abundance_weighted = FALSE)

# Graficar la redundancia funcional por módulo
t14$plot_func_FR(order_x = paste0("M", 1:10))





#### RED PARA CADA GRUPO DE LA METADATA

# Filtrar muestras con >= 10,000 lecturas
keep_samples <- colSums(mt$otu_table) >= 10000
mt_filt <- mt$clone()
mt_filt$otu_table <- mt$otu_table[, keep_samples]
mt_filt$sample_table <- mt$sample_table[keep_samples, ]
mt_filt$tidy_dataset()

mt_func <- mt_filt$clone()


### RENOMBRAR CON TAXONOMICA LOS ASVs
tax <- mt_func$tax_table

# Crear etiqueta con el nivel más específico disponible
tax$Tax_label <- apply(tax, 1, function(x) {
  x <- x[!is.na(x) & x != ""]
  tail(x, 1)
})


add_taxonomy_to_graph <- function(g, tax_table) {
  
  V(g)$Taxonomy <- tax_table[V(g)$name, "Tax_label"]
  
  # Si alguna ASV no tiene taxonomía
  V(g)$Taxonomy[is.na(V(g)$Taxonomy)] <- "Unclassified"
  
  return(g)
}


add_taxonomy_to_graph <- function(g, tax_table) {
  
  V(g)$Taxonomy <- tax_table[V(g)$name, "Tax_label"]
  V(g)$Taxonomy[is.na(V(g)$Taxonomy)] <- "Unclassified"
  
  return(g)
}


groups <- unique(mt_func$sample_table$Localities)#### Cambiar segun la meatada (Localities, m.s.n.m, Texture)

print(groups)
print(length(groups))


for (gname in groups) {
  
  message("Procesando grupo: ", gname)
  
  # Subconjunto de muestras
  samples_g <- rownames(
    mt_func$sample_table[mt_func$sample_table$Localities == gname, ]
  )
  
  # Subconjunto del microtable
  mt_g <- mt_func$clone()
  mt_g$otu_table <- mt_func$otu_table[, samples_g]
  mt_g$sample_table <- mt_func$sample_table[samples_g, ]
  mt_g$tidy_dataset()
  
  # Construcción de la red
  net_g <- trans_network$new(
    dataset    = mt_g,
    taxa_level = "OTU",
    cor_method = "spearman",
    filter_thres = 0.001
  )
  
  net_g$cal_network(
    COR_p_thres = 0.01,
    p_adjust_method = "fdr",
    COR_optimization = TRUE
  )
  
  net_g$cal_module(method = "cluster_fast_greedy")
  
  # 👉 AÑADIR TAXONOMÍA AQUÍ (CLAVE)
  g <- net_g$res_network
  g <- add_taxonomy_to_graph(g, tax)
  
  # Exportar a Gephi
  write.gexf(
    igraph.to.gexf(g),
    output = paste0("redes/Red_funcional_", gname, "_tax.gexf")
  )
}

meco_module_g <- net_g$trans_comm(use_col = "module")
tf_module <- trans_func$new(meco_module_g)
tf_module$cal_func(fungi_database = "FungalTraits")
tf_module$cal_func_FR(abundance_weighted = TRUE)
tf_module$res_func_FR


apply(tf_module$res_func_FR, 2, function(x) names(which.max(x)))


##### NUEVO ENFOQUE, RED GLOBAL


net_global <- trans_network$new(
  dataset = mt,
  cal_cor = "WGCNA",
  taxa_level = "OTU",
  filter_thres = 0.0001,
  cor_method = "spearman"
)

net_global$cal_network(COR_p_thres = 0.05, COR_cut = 0.7) #### COR_cut=0.7
net_global$cal_module()

module_assign <- net_global$res_module




g_global <- net_global$res_network

module_assign <- data.frame(
  module = V(g_global)$module,
  row.names = V(g_global)$name
)

head(module_assign)











add_module_to_graph <- function(g, module_table) {
  
  nodes <- V(g)$name
  
  modules <- module_table[nodes, "module"]
  
  modules[is.na(modules)] <- "Unassigned"
  
  V(g)$Module <- modules
  
  return(g)
}

module_colors <- c(
  M1 = "#C77CFF",
  M2 = "#7CAE00",
  M3 = "#4D4D4D",
  M4 = "#00BFC4",
  M5 = "#F8766D",
  M6 = "#E58700",
  M7 = "#00BA38",
  M8 = "#FF61C3",
  M9 = "#999999",
  M10 = "#A3A500",
  M11 = "#00A9FF",
  Unassigned = "#000000"
)



groups <- unique(mt$sample_table$Texture)


sum(V(g_global)$name %in% rownames(module_assign))


for (gname in groups) {
  
  message("Procesando grupo: ", gname)
  
  samples_g <- rownames(
    mt$sample_table[mt$sample_table$Texture == gname, ]
  )
  
  mt_g <- mt$clone()
  mt_g$otu_table <- mt$otu_table[, samples_g]
  mt_g$sample_table <- mt$sample_table[samples_g, ]
  mt_g$tidy_dataset()
  
  net_g <- trans_network$new(
    dataset = mt_g,
    taxa_level = "OTU",
    cor_method = "spearman",
    filter_thres = 0.001
  )
  
  net_g$cal_network(
    COR_p_thres = 0.01,
    p_adjust_method = "fdr"
  )
  
  g <- net_g$res_network
  
  # Añadir módulo global
  g <- add_module_to_graph(g, module_assign)
  
  # Asignar color fijo por módulo
  modules <- V(g)$Module
  colors_hex <- module_colors[modules]
  colors_hex[is.na(colors_hex)] <- "#000000"
  
  V(g)$color <- colors_hex
  
  # Reemplazar el objeto de red
  net_g$res_network <- g
  
  # Exportar correctamente
  net_g$res_network <- g
  
  net_g$save_network(
    filepath = paste0("redes1/Red_", gname, "_mod_global.gexf")
  )
}



tf_global <- trans_func$new(mt)
tf_global$cal_func(fungi_database = "FungalTraits")
head(tf_global$res_func)

module_assign


func_table <- tf_global$res_func

func_module <- merge(
  module_assign,
  func_table,
  by = "row.names",
  all.x = TRUE
)

colnames(func_module)[1] <- "OTU"

library(dplyr)

module_roles <- func_module %>%
  group_by(module) %>%
  summarise(across(where(is.numeric),
                   ~ sum(.x > 0, na.rm = TRUE)))


mat <- as.data.frame(module_roles)
rownames(mat) <- mat$module
mat$module <- NULL

# quitar columnas vacías
mat <- mat[, colSums(mat) > 0]

pheatmap(mat, scale = "none")



colnames(tf_global$res_func)

func_table <- tf_global$res_func

# SAPROTROPH
sap_cols <- grep("saprotroph|decomposer", colnames(func_table), value = TRUE)

# PLANT PATHOTROPH
plant_path_cols <- grep("plant_pathogen|root_pathogen|wood_pathogen|leaf/fruit/seed_pathogen", 
                        colnames(func_table), value = TRUE)

# ANIMAL PATHOTROPH
animal_path_cols <- grep("animal_parasite|vertebrate_parasite|arthropod_parasite|nematophagous", 
                         colnames(func_table), value = TRUE)

# SYMBIOTROPH (micorrizas + líquenes + unspecified_symbiotroph)
symbio_cols <- grep("mycorrhizal|lichenized|unspecified_symbiotroph|moss_symbiont|termite_symbiont", 
                    colnames(func_table), value = TRUE)

# ENDOPHYTE
endo_cols <- grep("endophyte", colnames(func_table), value = TRUE)



func_major <- data.frame(
  Saprotroph = ifelse(rowSums(func_table[, sap_cols], na.rm = TRUE) > 0, 1, 0),
  
  Plant_Pathotroph = ifelse(rowSums(func_table[, plant_path_cols], na.rm = TRUE) > 0, 1, 0),
  
  Animal_Pathotroph = ifelse(rowSums(func_table[, animal_path_cols], na.rm = TRUE) > 0, 1, 0),
  
  Symbiotroph = ifelse(rowSums(func_table[, symbio_cols], na.rm = TRUE) > 0, 1, 0),
  
  Endophyte = ifelse(rowSums(func_table[, endo_cols], na.rm = TRUE) > 0, 1, 0)
)

rownames(func_major) <- rownames(func_table)


func_module <- merge(
  module_assign,
  func_major,
  by = "row.names",
  all.x = TRUE
)

colnames(func_module)[1] <- "OTU"

library(dplyr)

module_prop <- func_module %>%
  group_by(module) %>%
  summarise(
    Saprotroph       = mean(Saprotroph, na.rm = TRUE),
    Plant_Pathotroph = mean(Plant_Pathotroph, na.rm = TRUE),
    Animal_Pathotroph= mean(Animal_Pathotroph, na.rm = TRUE),
    Symbiotroph      = mean(Symbiotroph, na.rm = TRUE),
    Endophyte        = mean(Endophyte, na.rm = TRUE),
    n_OTUs           = n()
  )


mat_prop <- as.data.frame(module_prop)
rownames(mat_prop) <- mat_prop$module
mat_prop$module <- NULL
mat_prop$n_OTUs <- NULL

library(pheatmap)

pheatmap(
  mat_prop,
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = NA
)





func_table <- tf_global$res_func
cn <- colnames(func_table)

strategy_patterns <- list(
  
  Saprotroph = "saprotroph|decomposer",
  
  Plant_pathogen = "plant_pathogen|root_pathogen|wood_pathogen|leaf/fruit/seed_pathogen",
  
  Animal_parasite = "animal_parasite|vertebrate_parasite|arthropod_parasite",
  
  Mycoparasite = "mycoparasite",
  
  Mycorrhizal = "mycorrhizal|ectomycorrhizal|arbuscular_mycorrhizal|ericoid_mycorrhizal",
  
  Endophyte = "endophyte",
  
  Lichenized = "lichenized",
  
  Symbiotroph_unspecified = "unspecified_symbiotroph",
  
  Wood_decay = "wood_saprotroph|Decay_substrate\\|wood",
  
  Soil_decomposer = "soil_saprotroph|Decay_substrate\\|soil",
  
  Litter_decomposer = "litter_saprotroph",
  
  Dung_decomposer = "dung_saprotroph|Decay_substrate\\|dung",
  
  Root_associated = "root-associated",
  
  Leaf_seed_decomposer = "leaf/fruit/seed",
  
  Arthropod_associated = "arthropod-associated",
  
  Vertebrate_associated = "vertebrate-associated",
  
  Nematophagous = "nematophagous",
  
  White_rot = "white_rot",
  
  Brown_rot = "brown_rot",
  
  Soft_rot = "soft_rot",
  
  Filamentous = "filamentous_mycelium",
  
  Yeast = "yeast",
  
  Zoosporic = "zoosporic"
)

func_strategies <- sapply(strategy_patterns, function(pattern) {
  
  cols <- grep(pattern, cn, value = TRUE)
  
  if (length(cols) == 0) {
    return(rep(0, nrow(func_table)))
  }
  
  as.numeric(rowSums(func_table[, cols, drop = FALSE], na.rm = TRUE) > 0)
})

func_strategies <- as.data.frame(func_strategies)
rownames(func_strategies) <- rownames(func_table)

func_module <- merge(
  module_assign,
  func_strategies,
  by = "row.names",
  all.x = TRUE
)

colnames(func_module)[1] <- "OTU"

library(dplyr)

module_prop20 <- func_module %>%
  group_by(module) %>%
  summarise(across(
    where(is.numeric),
    ~ mean(.x, na.rm = TRUE)
  ))

mat20 <- as.data.frame(module_prop20)
rownames(mat20) <- mat20$module
mat20$module <- NULL

pheatmap(
  mat20,
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = NA
)


mat20_t <- t(mat20)

pheatmap(
  mat20_t,
  scale = "none",
  cluster_rows = TRUE,     # cluster estrategias
  cluster_cols = TRUE,     # cluster módulos
  border_color = "black",  # activa bordes
  cellwidth = 20,
  cellheight = 12,
  fontsize = 10,
  angle_col = 0,
  fontsize_col = 9
)


for (gname in groups) {
  
  message("Procesando grupo: ", gname)
  
  samples_g <- rownames(
    mt$sample_table[mt$sample_table$Texture == gname, ]
  )
  
  mt_g <- mt$clone()
  mt_g$otu_table <- mt$otu_table[, samples_g]
  mt_g$sample_table <- mt$sample_table[samples_g, ]
  mt_g$tidy_dataset()
  
  # Construir red SIN detectar módulos
  net_g <- trans_network$new(
    dataset = mt_g,
    taxa_level = "OTU",
    cor_method = "spearman",
    filter_thres = 0.001
  )
  
  net_g$cal_network(
    COR_p_thres = 0.01,
    p_adjust_method = "fdr"
  )
  
  g <- net_g$res_network
  
  # 🔥 AÑADIR MÓDULO GLOBAL
  g <- add_module_to_graph(g, module_assign)
  
  # Exportar
  write.gexf(
    igraph.to.gexf(g),
    output = paste0("redes1/Red_", gname, "_mod_global.gexf")
  )
}


table(V(g)$Module)


