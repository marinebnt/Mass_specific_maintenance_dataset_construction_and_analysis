# ============================
# PCA + heatmap visualization
# ============================

AXESTOREPRESENT <- c(1,2)

# side note : I have compared the PCA and prcomp methods to produce the PCA, and I prefere the results produced by the PCA method. 
# PCA method has a first dimension with much higher values of the traits. 



contour_plot <- function(data_contour){

  # 1. RUN PCA

  pca <- runpca(method = method_plot, data_pca = data_contour)
  pca_coords <- pca[[1]]
  var_coords <- pca[[2]]
  var_expl <- pca[[3]]
  datatoadd <- data_to_add_func(data_contour, pca)
  # ----------------------------
  # 2. Visualization
  # ----------------------------
  interp_metab <- with(pca_coords,
                       akima::interp(x = Dim.1, y = Dim.2, z = metabolism,
                                     duplicate = "mean", linear = TRUE, nx = 200, ny = 200)
  )

  # Convert interpolation to dataframe
  interp_df <- as.data.frame(interp2xyz(interp_metab, data.frame = TRUE))
  # interp_df[is.na(interp_df)] <- 0

  ggplot(interp_df, aes(x = x, y = y)) +
    geom_contour_filled(bins = 15, aes(z = z)) +
    scale_fill_manual(
      values = colorRampPalette(c("blue","lightblue","yellow","red"))(15)
    ) +
    theme_minimal() +
    labs(x = paste0("PC1 (", round(var_expl[1], digits = 1), "%)"), y = paste0("PC2 (", round(var_expl[2], digits = 1), "%)"), fill = "RMR0 contpurs") +   # <- legend title
    theme(legend.position = "right")+
    new_scale_fill()+
    geom_point(
      data = datatoadd,
      aes(
        x = x, 
        y = y,
        shape = Archetypes,       
        fill = Archetypes         
      ),
      colour = "black",          
      size = 6,
      stroke = 1
    ) +
    scale_shape_manual(values = datatoadd$pchvec) +
    scale_fill_manual(values = datatoadd$colAA)+
    geom_segment(data = var_coords,
                 aes(x = 0, y = 0, xend = Dim.1*5, yend = Dim.2*5),
                 arrow = arrow(length = unit(0.25, "cm")),
                 color = "black", size = 0.8) +
    geom_text_repel(data = var_coords,
                    aes(x = Dim.1*5, y = Dim.2*5, label = varname),
                    color = "black", vjust = -0.5, size = 4, fontface = "bold")
}


runpca <- function(method, data_pca){
  # stop((method %in% c("prcomp", "PCA")))

  # Standardize automatically
  if (method == "PCA") {pca_res <- PCA(data_pca, scale.unit=T, graph = FALSE)}
  if (method == "prcomp") {pca_res <- prcomp(data_pca, scale. = T, center=T)}
  
  # Extract coordinates of individuals
  if (method == "PCA") {pca_coords <- as.data.frame(pca_res$ind$coord)
  pca_coords$Dim.1 <- -pca_coords$Dim.1}
  if (method == "prcomp") {
    pca_coords <- as.data.frame(pca_res$x)
    colnames(pca_coords) <- paste0("Dim.", 1:ncol(pca_coords))}
  
  # Add metabolism to dataframe
  pca_coords$metabolism <- metab
  pca_coords$species <- rownames(species_data)
  
  # Variable loadings (for arrows)
  if(method == "PCA") {var_coords <- as.data.frame(pca_res$var$coord)
  var_coords$Dim.1 <- -var_coords$Dim.1}
  if (method == "prcomp") {
    var_coords <- as.data.frame(pca_res$rotation)
    colnames(var_coords) <- paste0("Dim.", 1:ncol(var_coords))}
  var_coords$varname <- rownames(var_coords)
  
  if (method == "PCA") {var_expl <- pca_res$eig[1:2,2]}
  if (method == "prcomp") {
    variances <- pca_res$sdev^2
    var_expl <- variances / sum(variances) * 100}
  
  return(list(pca_coords, var_coords, var_expl))
  
}


preparedataforplot_heatmap <- function(numbPCA1, numbPCA2, dataacp, AA, PCA){
  dataacp <- dataacp
  points_pca <-PCA[[1]]
  dataacp$Dim.1 <- points_pca[, numbPCA1]  # indexing the first column
  dataacp$Dim.2 <- points_pca[, numbPCA2]  # indexing the second column
  rotation_data <- PCA[[2]][, c(numbPCA1, numbPCA2)]
  colnames(rotation_data) <- c(paste0("Dim.", AXESTOREPRESENT[1]), paste0("Dim.", AXESTOREPRESENT[2]))
  data_center <- as.matrix(t(data.frame(rep(0,nrow(rotation_data)))))
  colnames(data_center) <- rownames(rotation_data)
  rownames(data_center) <- NULL
  arch_corrected <- AA$archetypes[, PCA[[2]]$varname]
  archetypepointsPCA <- as.matrix(arch_corrected) %*% as.matrix(PCA[[2]][, 1:2])
  # eigenval <- get_eigenvalue(PCA)[c(numbPCA1, numbPCA2),"variance.percent"]
  return(list(rotate=rotation_data, AArotate=archetypepointsPCA, dataacp))
}



data_to_add_func <- function(data_for_plot, pca){
  colAA <- c()
  labelAA <- c()
  pchvec <- c()
  iopp <- which.max(as.data.frame(AAtot$archetypes)$K) 
  colAA[iopp] = "#F8766D"
  labelAA[iopp] = "Fast"
  pchvec[iopp] = c(24)
  ieq <- which.min(as.data.frame(AAtot$archetypes)$K) 
  colAA[ieq] = "#00BFC4"
  labelAA[ieq] = "Slow"
  pchvec[ieq] = c(22)
  iper <- c(1:3)[-c(ieq, iopp)]
  colAA[iper] = "#7CAE00"
  labelAA[iper] ="Intermediate"
  pchvec[iper] = c(21)
  
  listforplot <- preparedataforplot_heatmap(numbPCA1=1, numbPCA2=2, dataacp=data_for_plot, AA=AAtot, PCA=pca)
  rotation = listforplot[[1]]
  matrixAAinPCA = listforplot[[2]]
  dataacpPLOT = listforplot[[3]]
  dataacpPLOT = dataacpPLOT[, !duplicated(t(dataacpPLOT))]
  datatoadd <- matrixAAinPCA[,c(1,2)]
  rownames(datatoadd)<- NULL
  labels <- labelAA[1:nrow(AAtot$archetypes)]
  datatoadd <- data_frame(x=datatoadd[,1], y=datatoadd[,2], z=labels, pchvec=pchvec[1:nrow(AAtot$archetypes)])
  datatoadd$Archetypes <- labels
  datatoadd$Archetypes <- as.factor(labels)
  datatoadd$colAA <- colAA[1:nrow(AAtot$archetypes)]
  order_table <- order(datatoadd$Archetypes)
  order_table <- c(which((factor(labelAA)) == "Fast"),which((factor(labelAA)) == "Slow"),which((factor(labelAA)) == "Intermediate"))[1:nrow(AAtot$archetypes)]
  datatoadd <- datatoadd[order_table,]
  order_table <- c(which(levels(factor(labelAA)) == "Fast"),which(levels(factor(labelAA)) == "Slow"),which(levels(factor(labelAA)) == "Intermediate"))
  if(length(labels)>2){datatoadd$Archetypes <- factor(datatoadd$Archetypes, levels = levels(datatoadd$Archetypes)[order_table])}
  
  return(datatoadd)
}



for (method_plot in c("PCA", "prcomp")){
  

  
  ##########
  ######## CONTOURPLOT
  ##########
  
  species_data <- dataphylo
  
  # ----------------------------
  # 1. Load and prepare data
  # ----------------------------
  # Replace with your own data_contour
  # rows = species, columns = traits (including metabolism)
  # species_data <- read.csv("fish_traits.csv")
  
  # Example: separate focal trait (metabolism) from other traits
  repo <- stringr::str_remove(dir, "phylosem")
  
  
  traits_all <- species_data %>%
    select(-c(c_m, na, Genus, Family, Order, Class, Species, SpecCode, Habitat, osmose))   # all traits except metabolism
  
  traits_THV <- traits_all %>%
    select(c(Fecundity, Age.mat, Age.max, K, Mortality))   # all traits except metabolism
  
  traits_morpho <- traits_all %>%
    select(c(Fecundity, Min_caudalpeduncle_depth, Temperature, 
             Max_body_depth, Max_body_width, habitatbenthopelagic, 
             habitatdemersal, habitatpelagic, Lower_jaw_length, Age.mat))   # all traits except metabolism
  
  
  metab <- species_data$c_m
  
  
  k=0
  
  # PLOT CONTOUR PLOT AND SAVE
  # for (data in list(traits_all, traits_morpho, traits_THV)){
  #   k <- k+1
  #   assign(x = paste0("plot", k), contour_plot(data))
  # }
  
  # ggsave(filename = paste0(path_plots, "/", method_plot, "contourplot", repo, ".jpeg"),
  #        print(ggarrange(plot1, plot2, plot3, nrow = 2, ncol=2)),
  #        height = 9, width = 11)  
  # ggsave(filename = paste0(path_plots, "/", method_plot, "contourplotTHV", repo, ".jpeg"),
  #        print(contour_plot(traits_THV)),
  #        height = 7, width = 9)

  
  
  
  
  
  
  ##########
  ######## HEATMAP
  ##########
  
  # 1. Prepare the matrix
  # remove non-numeric columns (e.g. species names)
  traits_all <- species_data %>%
    select(c(Fecundity, Age.mat, Age.max, K, Mortality))   # THV  
  traits_c_m <- species_data %>%
    select(c(Fecundity, Age.mat, Age.max, K, Mortality, c_m))   # THV
  metab <- species_data$c_m
  # optional: scale traits (so they’re comparable)
  # trait_mat_scaled <- scale(traits_all) # dataphylo already scaled
  
  # 2. Basic pca + heatmap 
  
  pca <- runpca(method = method_plot, data_pca = traits_all)
  pca_coords <- pca[[1]]
  var_coords <- pca[[2]]
  var_expl <- pca[[3]]
  
  datatoadd <- data_to_add_func(traits_all, pca)
  
  df <- pca_coords %>%
    select(Dim.1, Dim.2, metabolism)
  

  pmin <- ggplot(df, aes(x = Dim.1, y = Dim.2)) +
    stat_summary_2d(bins = 25, fun = min, aes(z = metabolism)) +
    scale_fill_gradient(low = "blue", high = "red") +
    theme_minimal() +
    labs(x = paste0("PC1 (", round(var_expl[1], digits = 1), "%)"), y = paste0("PC2 (", round(var_expl[2], digits = 1), "%)"), z = "Metabolism", title="Min RMR0 values per area over PCA")+
    geom_segment(data = var_coords,
                 aes(x = 0, y = 0, xend = Dim.1*5, yend = Dim.2*5),
                 arrow = arrow(length = unit(0.25, "cm")),
                 color = "black", size = 0.8) +
    geom_text_repel(data = var_coords,
                    aes(x = Dim.1*5, y = Dim.2*5, label = varname),
                    color = "black", vjust = -0.5, size = 5, fontface = "bold")
  pmin <- pmin+ new_scale_fill()+
    geom_point(
      data = datatoadd,
      aes(
        x = x, 
        y = y,
        shape = Archetypes,       
        fill = Archetypes         
      ),
      colour = "black",          
      size = 6,
      stroke = 1
    ) +
    scale_shape_manual(values = datatoadd$pchvec) +
    scale_fill_manual(values = datatoadd$colAA)
  
  
  rm(max)
  rm(mean)
  rm(median)
  pmax <- ggplot(df, aes(x = Dim.1, y = Dim.2)) +
    stat_summary_2d(bins = 25, fun = max, aes(z = metabolism)) +
    scale_fill_gradient(low = "blue", high = "red") +
    theme_minimal() +
    labs(x = paste0("PC1 (", round(var_expl[1], digits = 1), "%)"), y = paste0("PC2 (", round(var_expl[2], digits = 1), "%)"), z = "Metabolism", title="Max RMR0 values per area over PCA")+
    geom_segment(data = var_coords,
                 aes(x = 0, y = 0, xend = Dim.1*5, yend = Dim.2*5),
                 arrow = arrow(length = unit(0.25, "cm")),
                 color = "black", size = 0.8) +
    geom_text_repel(data = var_coords,
                    aes(x = Dim.1*5, y = Dim.2*5, label = varname),
                    color = "black", vjust = -0.5, size = 5, fontface = "bold")  
  pmax <- pmax+ new_scale_fill()+
    geom_point(
      data = datatoadd,
      aes(
        x = x, 
        y = y,
        shape = Archetypes,       
        fill = Archetypes         
      ),
      colour = "black",          
      size = 7,
      stroke = 1
    ) +
    scale_shape_manual(values = datatoadd$pchvec) +
    scale_fill_manual(values = datatoadd$colAA)
  
  pmean <- ggplot(df, aes(x = Dim.1, y = Dim.2)) +
    stat_summary_2d(bins = 25, fun = mean, aes(z = metabolism)) +
    scale_fill_gradient(low = "blue", high = "red") +
    theme_minimal() +
    labs(x = paste0("PC1 (", round(var_expl[1], digits = 1), "%)"), y = paste0("PC2 (", round(var_expl[2], digits = 1), "%)"), z = "Metabolism", title="Mean RMR0 values per area over PCA")+
    geom_segment(data = var_coords,
                 aes(x = 0, y = 0, xend = Dim.1*5, yend = Dim.2*5),
                 arrow = arrow(length = unit(0.25, "cm")),
                 color = "black", size = 0.8) +
    geom_text_repel(data = var_coords,
                    aes(x = Dim.1*5, y = Dim.2*5, label = varname),
                    color = "black", vjust = -0.5, size = 5, fontface = "bold")
  pmean <- pmean+ new_scale_fill()+
    geom_point(
      data = datatoadd,
      aes(
        x = x, 
        y = y,
        shape = Archetypes,       
        fill = Archetypes         
      ),
      colour = "black",          
      size = 6,
      stroke = 1
    ) +
    scale_shape_manual(values = datatoadd$pchvec) +
    scale_fill_manual(values = datatoadd$colAA)
  pmedian <- ggplot(df, aes(x = Dim.1, y = Dim.2)) +
    stat_summary_2d(bins = 25, fun = median, aes(z = metabolism), alpha=0.8) +
    scale_fill_gradient(low = "blue", high = "red") +
    theme_bw() +
    labs(x = paste0("PC1: Life-cycle speed (", round(var_expl[1], digits = 1), "%)"), 
         y = paste0("PC2: Reproductive strategy (", round(var_expl[2], digits = 1), "%)"), 
         fill = "Standardized\nRMR0")+#, title="Median RMR0 values per area over PCA")+
    geom_segment(data = var_coords,
                 aes(x = 0, y = 0, xend = Dim.1*5, yend = Dim.2*5),
                 arrow = arrow(length = unit(0.25, "cm")),
                 color = "black", size = 0.8) +
    geom_text_repel(data = var_coords,
                    aes(x = Dim.1*5, y = Dim.2*5, label = varname),
                    color = "black", vjust = -0.5, size = 5, fontface = "bold")
  pmedian <- pmedian+ new_scale_fill()+
    geom_point(
      data = datatoadd,
      aes(
        x = x, 
        y = y,
        shape = Archetypes,       
        fill = Archetypes         
      ),
      colour = "black",          
      size = 8.5,
      stroke = 1
    ) +
    scale_shape_manual(values = datatoadd$pchvec) +
    scale_fill_manual(values = datatoadd$colAA)
  
  ggsave(filename = paste0(path_plots, "/", method_plot, "heatmap", repo, ".jpeg"),
         ggarrange(pmin, pmax, pmean, pmedian, ncol=2, nrow=2),
         height = 9, width = 11)
  ggsave(filename = paste0(path_plots, "/", method_plot, "MEDIANheatmap", repo, ".jpeg"),
         print(pmedian),
         height = 6, width = 9)

  
  
  
  
  ##########
  ######## Regression 
  ##########  
  
  
  # Open 3D plot
  
  # install.packages("reticulate")
  # reticulate::py_install("kaleido")
  
  # Fit regression (linear here, but could be GAM too)
  fit1 <- lm(metabolism ~ Dim.1 + Dim.2, data = pca_coords)
  fit2 <- lm(metabolism ~ Dim.1 : Dim.2, data = pca_coords)
  fit3 <- lm(metabolism ~ Dim.1 * Dim.2, data = pca_coords)
  anova(fit1, fit3)
  anova(fit2, fit3) 
  fit <- fit3
  
  # Create prediction grid
  x_seq <- seq(min(pca_coords$Dim.1), max(pca_coords$Dim.1), length = 50)
  y_seq <- seq(min(pca_coords$Dim.2), max(pca_coords$Dim.2), length = 50)
  grid <- expand.grid(Dim.1 = x_seq, Dim.2 = y_seq)
  grid$metabolism_pred <- predict(fit, newdata = grid)
  
  # Reshape for surface plotting
  z_matrix <- matrix(grid$metabolism_pred, nrow = length(x_seq), ncol = length(y_seq))
  
  # 3D scatter + regression surface
  preg <- plot_ly() %>%
    add_markers(data = pca_coords,
                x = ~Dim.1, y = ~Dim.2, z = ~metabolism,
                marker = list(color = 'red', size = 3),
                name = "Observed") %>%
    add_surface(x = x_seq, y = y_seq, z = z_matrix,
                colorscale = list(c(0,1), c("blue","red")),
                opacity = 0.7,
                name = "Regression surface") %>%
    layout(scene = list(
      xaxis = list(title = "PC1"),
      yaxis = list(title = "PC2"),
      zaxis = list(title = "Metabolism")
    ))
  
  
  saveWidget(preg, paste0(path_plots, "/", method_plot, "regression_", repo,".html"), selfcontained = TRUE)
  
  
  ##########
  ######## Project regression 
  ##########  
  
  # prepare input dataset
  metabolism_predicted <- predict(fit)
  df_reg <- cbind(pca_coords$Dim.1, pca_coords$Dim.2, metabolism_predicted)
  colnames(df_reg)[1:2] <- c("Dim.1", "Dim.2") 

  pmedian_predicted  <- ggplot(df_reg, aes(x = Dim.1, y = Dim.2)) +
    stat_summary_2d(bins = 25, fun = median, aes(z = metabolism_predicted), alpha=0.8) +
    scale_fill_gradient(low = "blue", high = "red") +
    theme_bw() +
    labs(x = paste0("PC1: Life-cycle speed (", round(var_expl[1], digits = 1), "%)"), 
         y = paste0("PC2: Reproductive strategy (", round(var_expl[2], digits = 1), "%)"), 
         fill = "Standardized\nRMR0")+ 
    #, title="Median RMR0 values per area over PCA")+
    geom_segment(data = var_coords,
                 aes(x = 0, y = 0, xend = Dim.1*5, yend = Dim.2*5),
                 arrow = arrow(length = unit(0.25, "cm")),
                 color = "black", size = 0.8) +
    geom_text_repel(data = var_coords,
                    aes(x = Dim.1*5, y = Dim.2*5, label = varname),
                    color = "black", vjust = -0.5, size = 5, fontface = "bold")
  pmedian_predicted  <- pmedian_predicted + new_scale_fill()+
    geom_point(
      data = datatoadd,
      aes(
        x = x, 
        y = y,
        shape = Archetypes,       
        fill = Archetypes         
      ),
      colour = "black",          
      size = 8.5,
      stroke = 1
    ) +
    scale_shape_manual(values = datatoadd$pchvec) +
    scale_fill_manual(values = datatoadd$colAA)
    
  
  
  ggsave(filename = paste0(path_plots, "/", method_plot, "heatmapMEDIANPREDICTED", repo, ".jpeg"),
         ggarrange(pmedian, pmedian_predicted, ncol=2, nrow=1, common.legend = T, legend = "right", labels = c("A", "B")),
         height = 6, width = 12)  
  ggsave(filename = paste0(path_plots, "/", method_plot, "heatmapMEDIANPREDICTED", repo, ".pdf"),
         ggarrange(pmedian, pmedian_predicted, ncol=2, nrow=1, common.legend = T, legend = "right", labels = c("A", "B")),
         height = 6, width = 12)
  
  
  #####
  # save summary of the linear model_h 
  #####
  
  
  model_h <- fit
  library(dplyr) # for pipe (%>%) command
  library(pixiedust)
  
  dust(model_h) %>% 
    sprinkle(cols = c("estimate", "std.error", "statistic"), round = 2) %>%
    sprinkle(cols = "p.value", fn = quote(pvalString(value))) %>% 
    sprinkle_colnames("Term", "Coefficient", "SE", "T-statistic", 
                      "P-value")

}



## plot correlation circle of c_m
# 
# res.pca <- PCA(traits_c_m, graph = FALSE)
# plot(res.pca, choix = "var")
# 



# 
# 
# 
# 
# 
# 
#   
#   p <- ggplot(pca_coords, aes(x = Dim.1, y = Dim.2, color = metabolism)) +
#     geom_point(size = 4, alpha = 0.8) +
#     geom_density_2d(color = "grey40", alpha = 0.5) +
#     scale_color_viridis_c(option = "plasma") +   # heatmap gradient
#     theme_minimal() +
#     labs(
#       x = paste0("PC1 (", round(pca_res$eig[1,2], 1), "%)"),
#       y = paste0("PC2 (", round(pca_res$eig[2,2], 1), "%)"),
#       color = "Metabolism"
#     ) +
#     theme(
#       text = element_text(size = 14),
#       legend.position = "right"
#     )
#   
#   # Print plot
#   print(p)
#   
#   ggplot(pca_coords, aes(x = Dim.1, y = Dim.2)) +
#     geom
#   # filled contour plot for metabolism
#   stat_summary_2d(aes(z = metabolism, fill = ..value..), bins = 15) +
#     scale_fill_viridis_c(option = "plasma") +
#     # species points on top
#     # geom_point(aes(color = metabolism), size = 3) +
#     scale_color_viridis_c(option = "plasma") +
#     theme_minimal() +
#     labs(
#       x = paste0("PC1 (", round(pca_res$eig[1,2], 1), "%)"),
#       y = paste0("PC2 (", round(pca_res$eig[2,2], 1), "%)"),
#       fill = "Metabolism",
#       color = "Metabolism"
#     ) +
#     theme(
#       text = element_text(size = 14),
#       legend.position = "right"
#     )
#   
#   
#   
#   
#   
#   # ----------------------------
#   # 2. PCA on all other traits
#   # ----------------------------
#   pca_res <- PCA(traits_THV, scale.unit = TRUE, graph = FALSE)
#   
#   # Species coordinates
#   pca_coords <- as.data.frame(pca_res$ind$coord)
#   pca_coords$metabolism <- metab
#   pca_coords$species <- rownames(species_data)
#   
#   # Variable loadings (for arrows)
#   var_coords <- as.data.frame(pca_res$var$coord)
#   var_coords$varname <- rownames(var_coords)
#   
#   # ----------------------------
#   # 3. Interpolate metabolism values in PCA space
#   # ----------------------------
#   interp_metab <- with(pca_coords,
#                        akima::interp(x = Dim.1, y = Dim.2, z = metabolism,
#                                      duplicate = "mean", linear = TRUE, nx = 200, ny = 200)
#   )
#   
#   # Convert to dataframe for ggplot
#   interp_df <- interp2xyz(interp_metab, data.frame = TRUE)
#   
#   # ----------------------------
#   # 4. Plot PCA biplot with contours
#   # ----------------------------
#   ggplot() +
#     # filled contour background
#     geom_contour_filled(data = interp_df,
#                         aes(x = x, y = y, z = z), bins = 15, alpha=0.7) +
#     # scale_fill_viridis_d(option = "plasma") +
#     scale_fill_manual(values = colorRampPalette(c("blue", "yellow"))(15))+
#     # add species points
#     # geom_point(data = pca_coords,
#     #            aes(x = Dim.1, y = Dim.2),
#     #            color = "black", size = 2, alpha = 0.7) +
#     
#     # add variable arrows (loadings)
#     geom_segment(data = var_coords,
#                  aes(x = 0, y = 0, xend = Dim.1*5, yend = Dim.2*5),
#                  arrow = arrow(length = unit(0.25, "cm")),
#                  color = "black", size = 0.8) +
#     
#     # add variable names
#     geom_text_repel(data = var_coords,
#                     aes(x = Dim.1*5, y = Dim.2*5, label = varname),
#                     color = "black", vjust = -0.5, size = 4) +
#     
#     theme_minimal() +
#     labs(
#       x = paste0("PC1 (", round(pca_res$eig[1,2], 1), "%)"),
#       y = paste0("PC2 (", round(pca_res$eig[2,2], 1), "%)"),
#       fill = "Metabolism"
#     ) +
#     theme(
#       text = element_text(size = 14),
#       legend.position = "right"
#     ) 
#   
#   
#   