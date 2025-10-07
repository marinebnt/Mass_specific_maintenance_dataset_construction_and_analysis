
####  COMPARE THE c_m DATASETS : create a csv   ####


original <- read.csv("01-Dataset_construction/Outputs/phylosem_before25-09-16/phylosem_before_removingElasmo/TLstdmeca/output_SEMpsemFINALtot.csv")
nofundulus <- read.csv("01-Dataset_construction/Outputs/phylosem_nofundulus/phylosem/output_SEMpsemFINALtot.csv")
nofunoela <- read.csv("01-Dataset_construction/Outputs/phylosem_nofundulus_noelasmo/phylosem/output_SEMpsemFINALtot.csv")
killen <- read.csv("01-Dataset_construction/Outputs/phylosem_killen_nofundulus_noelasmo/phylosem/output_SEMpsemFINALtot.csv")
newest <- read.csv("01-Dataset_construction/Outputs/phylosem_clarke2024/phylosem/output_SEMpsemFINALtot.csv")
alldataset_ela <- read.csv("01-Dataset_construction/Outputs/phylosem_alldatasets_ela/output_SEMpsemFINALtot.csv")
alldataset <- read.csv("01-Dataset_construction/Outputs/phylosem_alldatasets/output_SEMpsemFINALtot.csv")

df <- data.frame("label"=original$label, "c_m"=original$c_m)
df <- dplyr::left_join(df, data.frame("label"=nofundulus$label, "c_mnofundulus"=nofundulus$c_m))
df <- dplyr::left_join(df, data.frame("label"=nofunoela$label, "c_mnofunoela"=nofunoela$c_m))
df <- dplyr::left_join(df, data.frame("label"=killen$label, "c_mkillen"=killen$c_m))
df <- dplyr::left_join(df, data.frame("label"=newest$label, "c_mclarke"=newest$c_m))
df <- dplyr::left_join(df, data.frame("label"=alldataset_ela$label, "c_malldataset_ela"=alldataset_ela$c_m))
df <- dplyr::left_join(df, data.frame("label"=alldataset$label, "c_malldataset"=alldataset$c_m))

write.csv(df, "compare_c_m.csv")


#### COMPARE THE c_m : produce plots and save them ####

library(ggplot2)

for (name in colnames(df)[-c(1:2)]){
  plot <- ggplot(df, aes(x=c_m, y=get(name)))+
    ggplot2::geom_point()+
    geom_abline(slope=1, intercept=0, col="red")+
    geom_smooth()+
    labs(title=paste0("Comparing c_m original and ", name),
        y = name)+
    geom_text(aes(label = label), check_overlap = TRUE)
  
  ggsave(paste0("c_m_comparison_", name, ".jpg"), plot)
}

  

#### COMPARE THE path coefficients : produce plots and save them ####

original <- read.csv("01-Dataset_construction/Outputs/phylosem_before25-09-16/phylosem_before_removingElasmo/TLstdmeca/coefnostd_std_SEM.csv")
nofundulus <- read.csv("01-Dataset_construction/Outputs/phylosem_nofundulus/phylosem/coefnostd_std_SEM.csv")
nofunoela <- read.csv("01-Dataset_construction/Outputs/phylosem_nofundulus_noelasmo/phylosem/coefnostd_std_SEM.csv")
killen <- read.csv("01-Dataset_construction/Outputs/phylosem_killen_nofundulus_noelasmo/phylosem/coefnostd_std_SEM.csv")
clarke <- read.csv("01-Dataset_construction/Outputs/phylosem_clarke2024/phylosem/coefnostd_std_SEM.csv")
alldataset_ela <- read.csv("01-Dataset_construction/Outputs/phylosem_alldatasets_ela/coefnostd_std_TLstdmeca.csv")
alldataset <- read.csv("01-Dataset_construction/Outputs/phylosem_alldatasets/coefnostd_std_TLstdmeca.csv")

df <- data.frame("label"=original$Path, "PC"=original$EstimateStd)
df <- dplyr::left_join(df, data.frame("label"=nofundulus$Path, "PC_nofundulus"=nofundulus$EstimateStd))
df <- dplyr::left_join(df, data.frame("label"=nofunoela$Path, "PC_nofunoela"=nofunoela$EstimateStd))
df <- dplyr::left_join(df, data.frame("label"=killen$Path, "PC_killen"=killen$EstimateStd))
df <- dplyr::left_join(df, data.frame("label"=clarke$Path, "PC_newest"=clarke$EstimateStd))
df <- dplyr::left_join(df, data.frame("label"=alldataset_ela$Path, "PC_alldataset_ela"=alldataset_ela$EstimateStd))
df <- dplyr::left_join(df, data.frame("label"=alldataset$Path, "PC_alldataset"=alldataset$EstimateStd))


library(ggplot2)

for (name in colnames(df)[-c(1:2)]){
  plot <- ggplot(df, aes(x=log(PC), y=log(get(name))))+
    ggplot2::geom_point()+
    geom_abline(slope=1, intercept=0, col="red")+
    # geom_smooth()+
    labs(title=paste0("Comparing path coefficients original and ", name),
         y = paste0("log(",name, ")"))+
    geom_text_repel(aes(label = label), max.overlaps = 20, col="grey55") #, check_overlap = TRUE)
  
  ggsave(paste0("PC_", name, ".jpg"), plot, height=9, width=6)
}




#### COMPARE THE SEMs - graphs plots  ####

library(igraph)
library(ggraph)
library(ggplot2)
library(dplyr)
library(patchwork)
library(ggnetwork)
library(ggrepel)



# --- YOUR DATA: make sure columns are exactly these names ----
original <- read.csv("01-Dataset_construction/Outputs/phylosem_before_removingElasmo/TLstdmeca/coefnostd_std_SEM.csv")
nofundulus <- read.csv("01-Dataset_construction/Outputs/phylosem_NOFUNDULUS/TLstdmeca/coefnostd_std_TLstdmeca.csv")
nofunoela <- read.csv("01-Dataset_construction/Outputs/phylosem_NOFUNDULUS_NOELASMO/TLstdmeca/coefnostd_std_TLstdmeca.csv")
killen <- read.csv("01-Dataset_construction/Outputs/phylosem_WKILLEN_NOELASMO_NOFUNDULUS/TLstdmeca/coefnostd_std_TLstdmeca.csv")
newest <- read.csv("01-Dataset_construction/Outputs/phylosem/TLstdmeca/phylosem/coefnostd_std_SEM_all_datasets.csv")
gravel <- read.csv("01-Dataset_construction/Outputs/phylosem_GRAVEL/Elasmo/coefnostd_std_TLstdmeca.csv")
sem_compare <- list("no_fundulus"=nofundulus, "nofund_noela"=nofunoela, "killen"=killen, "newest"=newest, "gravel"=gravel)

original <- original[-grep("V", original$Parameter), ]
original$lhs <-  stringr::str_remove(pattern = " ", unlist(stringr::str_split(original$Path, pattern = "->"))[seq(1,2*length(original$Path), by=2)])
original$rhs <- stringr::str_remove(pattern = " ", unlist(stringr::str_split(original$Path, pattern = "->"))[seq(2,2*length(original$Path), by=2)]) 
original_base <- original %>%
  transmute(from = lhs, to = rhs,
            est_std_original   = EstimateStd)

k=0
for (sem in sem_compare){
  k=k+1
  dataset <- names(sem_compare)[k]
  
  sem     <- sem[-grep("V", sem$Parameter), ]
  sem$lhs <- stringr::str_remove(pattern = " ", unlist(stringr::str_split(sem$Path, pattern = "->"))[seq(1,2*length(sem$Path), by=2)])
  sem$rhs <- stringr::str_remove(pattern = " ", unlist(stringr::str_split(sem$Path, pattern = "->"))[seq(2,2*length(sem$Path), by=2)]) 
  
  # --- Create an edges table containing both estimate types as edge attributes ---
  edges_base <- sem %>%
    transmute(from = lhs, to = rhs,
              est_std = EstimateStd)
  
  #full join
  edges_base <- full_join(edges_base, original_base)
  
  # add dummy edges
  # Identify disconnected components
  g_tmp <- graph_from_data_frame(edges_base, directed = TRUE)
  components <- components(g_tmp)
  
  # Add dummy edges between components
  comp_ids <- unique(components$membership)
  dummy_edges <- data.frame()
  
  for (i in seq_along(comp_ids)[-1]) {
    from_node <- names(components$membership)[components$membership == comp_ids[i]][1]
    to_node <- names(components$membership)[components$membership == comp_ids[i - 1]][1]
    
    dummy_edges <- rbind(dummy_edges, data.frame(
      from = from_node,
      to   = to_node,
      est_std = 0,
      est_std_original = 0
    ))
  }
  
  # Combine with real edges
  edges_base <- rbind(edges_base, dummy_edges)
  
  # Build graph with both attributes present on edges
  g <- graph_from_data_frame(edges_base, directed = TRUE)
  
  # A single layout for consistent node positions across plots
  set.seed(42)
  layout <- create_layout(g, layout="fr")#, layout = "fr")   # force-directed layout sugiyama kk
  
  # Swap x and y to rotate 90 degrees
  layout$x_old <- layout$x
  layout$x <- layout$y
  layout$y <- layout$x_old
  layout$x_old <- NULL
  layout$x <- layout$x * 1  # scale horizontal spacing
  layout$y <- layout$y * 1  # scale vertical spacing
  layout_df <- as.data.frame(layout)[, c("name","x","y")]
  
  # Helper: compute explicit label x/y/angle for each edge
  compute_edge_labels <- function(edges_df, layout_df, est_col, offset = 0.08){
    edges_df %>%
      dplyr::left_join(layout_df,   by = c("from" = "name")) %>%
      rename(x_from = x, y_from = y) %>%
      dplyr::left_join(layout_df,   by = c("to"   = "name")) %>%
      rename(x_to = x, y_to = y) %>%
      mutate(
        x_mid = (x_from + x_to) / 2,
        y_mid = (y_from + y_to) / 2,
        dx = x_to - x_from,
        dy = y_to - y_from,
        len = sqrt(dx^2 + dy^2),
        # normalized perpendicular vector for offset
        ux = ifelse(len == 0, 0, dx / len),
        uy = ifelse(len == 0, 0, dy / len),
        perp_x = -uy,
        perp_y = ux,
        # offset distance -- tweak this number for more/less separation from the link
        off = offset,
        x_label = x_mid + perp_x * off,
        y_label = y_mid + perp_y * off,
        # angle for text (degrees), keep text upright
        angle = atan2(dy, dx) * 180 / pi,
        angle = ifelse(angle > 90, angle - 180,
                       ifelse(angle < -90, angle + 180, angle)),
        label = sprintf("%.3f", .data[[est_col]])
      )
  }
  
  # Compute label frames for both coefficient types
  edge_labels_unstd <- compute_edge_labels(edges_base, layout_df, "est_std", offset = 0.08)
  edge_labels_std   <- compute_edge_labels(edges_base, layout_df, "est_std_original",   offset = 0.08)
  
  # ---- Plot function re-using the same layout (g contains both edge attrs) ----
  make_plot <- function(layout_obj, coef_attr, edge_labels_df, title) {
    ggraph(layout_obj) +
      # edge line: map to the correct edge attribute name (est_unstd or est_std)
      geom_edge_link(aes_string(width = paste0("abs(", coef_attr, ")"),
                                color = coef_attr),
                     arrow = arrow(length = unit(4, "mm")),
                     end_cap = circle(6, "mm"),
                     lineend = "round") +
      scale_edge_width(range = c(0.5, 2)) +
      scale_edge_colour_gradient2(low = "red", mid = "grey80", high = "blue", midpoint = 0, name = "coef") +
      geom_node_point(size = 6) +
      # explicit edge labels (midpoint + perpendicular offset)
      geom_text_repel(data = edge_labels_df, 
                      aes(x = x_label, y = y_label, label = label, angle = angle),
                      size = 3, inherit.aes = FALSE) +
      geom_node_text(aes(label = name), vjust = 1.8, size = 3, repel=F) +
      theme_void() +
      ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5))
  }
  
  # Create both plots (they will share the same node layout)
  p_unstd <- make_plot(layout_obj = layout, coef_attr = "est_std", edge_labels_df = edge_labels_unstd,
                       title =  paste0("SEM ", dataset))
  p_std   <- make_plot(layout, "est_std_original",   edge_labels_std,   paste0("SEM original"))
  
  # Combine side-by-side and save
  combined <- p_unstd + p_std + plot_layout(ncol = 1, nrow=2, widths = c(1,1))
  
  ggsave(paste0("SEM", dataset, ".jpg"), combined, width = 14, height = 10, dpi = 300)
  print(combined)
  
}  
  
  


#### correlation matrix ###########

# =============================
# Trait correlation analysis in R
# =============================

# Load libraries
library(tidyverse)    # data manipulation
library(corrplot)     # nice correlation plots
library(Hmisc)        # rcorr() for correlations + p-values
library(PerformanceAnalytics) # optional: chart.Correlation

# -----------------------------
# 1. Load and prepare your data
# -----------------------------
# Example: species_data is your dataframe
# Rows = species, Columns = traits
# Make sure only numeric trait columns are selected

# Replace with your dataset
# species_data <- read.csv("your_species_traits.csv")

# Example subset of numeric trait columns
trait_data <- dataphylo[-which(colnames(dataphylo) %in% c("Species", "na", "Genus", "Family", "Order", "Class"))] 

# Ensure numeric
trait_data <- trait_data %>% mutate(across(everything(), as.numeric))

# --------------------------------
# 2. Correlation matrix + p-values
# --------------------------------
cor_results <- rcorr(as.matrix(trait_data), type = "pearson")

# Extract correlation matrix
cor_matrix <- cor_results$r
# Extract p-values
p_matrix <- cor_results$P

# -----------------------------
# 3. Visualize correlations
# -----------------------------
# Corrplot with significance
corrplot(cor_matrix, method = "color", 
         type = "upper", 
         tl.col = "black", 
         tl.srt = 45, 
         p.mat = p_matrix, sig.level = 0.05, insig = "blank")

# Optional: simple pairs plot with histograms
chart.Correlation(trait_data, histogram=TRUE, pch=19)

# -----------------------------
# 4. Export results (optional)
# -----------------------------
write.csv(cor_matrix, "trait_correlations.csv")
write.csv(p_matrix, "trait_pvalues.csv")

  
