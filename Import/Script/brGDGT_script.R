#### Library ####
library(reshape2)
library(ggplot2)
library(patchwork)
library(randomForest)
library(dplyr)
library(palaeoSig)
library(caret)
library(cluster)
library(FactoMineR)
library(factoextra)
library(e1071)
library(smotefamily)
library(tibble)
library(rgdal)
library(ggpmisc)
library(dismo) # BRT
library(ggnewscale)
library(ggforce)
library(RColorBrewer)

#### Graphical settings ####
Fill.scale <- scale_fill_gradientn(colors = c("grey95", "grey20"),
                                   name = "Count", guide = "none")

AI.class.scale <- scale_color_manual(name = "Aridity classes",
                                     values = c("1_Hyper-arid" = "#8c510a","2_Arid" = "#bf812e",
                                                "3_Semi-arid" = "#dfc27e","4_Dry sub-humid" = "#f5e9bf",
                                                "5_Humid" = "#80cec1", "Lacustrine" = "#75AADB",
                                                "Soil" = "#EB9122", "Moss" = "#9ABD55"), 
                                     labels = c("1_Hyper-arid" = "Hyper-arid", "2_Arid" = "Arid", "3_Semi-arid" = "Semi-arid", "Dry sub-humid", "Humid", "Lacustrine" = "Lacustrine", "Soil" = "Soil"))
Shape.scale <- scale_shape_manual(name = "Sample type", values = c("Soil" = 16, "Lacustrine" = 4))

Log10.scale <- scale_x_log10(
  breaks = scales::trans_breaks("log10", function(x) 10^x),
  labels = scales::trans_format("log10", scales::math_format(10^.x)))
Bin.size = 30; Annot.x = 0.95
R.size <- 4; V.step <- 0.085

myPalette <- grDevices::colorRampPalette(brewer.pal(11, "Spectral"))
sc <- scale_colour_gradientn(colours = myPalette(10))
sc2 <- scale_colour_manual(values = c("K-wet" = "3E8A70", "K-lake" = "#1EAEAE", "K-arid" = "#8c510a"))

#### Functions  ####

PCA.bioclim <- function(MP, transp_OK, Site.name, Type.samples, Ellipse, Shape = NULL, Show.centroid, Show.arrow, Show.site.lab = F,
                        Csv.sep, Scale.PCA, Groupes, Cluster.core, Cluster.core.lab, return.pick, Contrib, Size.MPS = 3, Alpha.MPS = 0.5,
                        Save.path, Manu.lim.x, Manu.lim.y, Dot.size, Dot.opac, Ellipse.opa, Density.contour, Doz.size.leg = NULL,
                        Opa.range, Reverse.dim, Show.annot, Show.Plotly, PCA.site, Marg.density.plot, Leg.nrow = NULL, GDGT = F,
                        Helinger.trans = F,
                        Symbol.path = NULL, Symbol.pos = NULL, Legend.position, Density.type, Save.plot, H, W, Num.facet, Legend.size){
  #### Settings ####
  library(vegan)
  library(ggrepel)
  library("FactoMineR") # FactoMineR pour la PCA (Le et al. 2008)
  library("factoextra")
  if(missing(Csv.sep)){Csv.sep = "\t"}
  if(missing(Site.name)){Site.name = "Site1"}
  if(missing(Cluster.core)){Cluster.core = NULL}
  if(missing(Num.facet)){Num.facet = NULL}
  if(missing(Ellipse.opa)){Ellipse.opa = 0.4}
  if(missing(Cluster.core.lab)){Cluster.core.lab = "Biome (Dinerstein et al., 2017)"}
  if(missing(PCA.site)){PCA.site = F}
  if(missing(Marg.density.plot)){Marg.density.plot = F}
  if(missing(Reverse.dim)){Reverse.dim = F}
  if(missing(Show.arrow)){Show.arrow = T}
  if(missing(Show.annot)){Show.annot = T}
  if(missing(Ellipse)){Ellipse = F}
  if(missing(Contrib)){Contrib = T}
  if(missing(return.pick)){return.pick = F}
  if(missing(Show.centroid)){Show.centroid = F}
  if(missing(transp_OK)){transp_OK = T}
  if(missing(Legend.position)){Legend.position = "right"}
  if(Show.centroid == F){Centroide = "quali"}
  if(Show.centroid == T){Centroide = NULL}
  if(missing(Legend.size)){Legend.size = 1}
  if(missing(Scale.PCA)){Scale.PCA = 1}
  if(missing(Save.path)){Save.path = NULL}
  if(missing(Groupes)){Groupes = NULL}
  if(missing(Manu.lim.x)){Manu.lim.x = NULL}
  if(missing(Manu.lim.y)){Manu.lim.y = NULL}
  if(missing(Type.samples)){Type.samples = NULL}
  if(missing(Show.Plotly)){Show.Plotly = F}
  if(missing(Save.plot)){Save.plot = NULL}
  if(missing(Density.contour)){Density.contour = F}
  if(missing(Density.type)){Density.type = "polygon"}
  if(missing(Dot.size)){Dot.size = NULL}
  if(missing(Dot.opac)){Dot.opac = NULL}
  if(missing(Opa.range)){Opa.range = c(0.01,0.1)}
  if(missing(W)){W = NULL}
  if(missing(H)){H = NULL}
  
  #### Save plots ####
  if(is.null(Save.plot) == F){
    Path.to.create <- gsub("(.*/).*\\.pdf.*","\\1", Save.plot)
    dir.create(file.path(Path.to.create), showWarnings = FALSE)
    if(is.null(W) == F & is.null(H) == F){
      pdf(file = Save.plot, width = W*0.01041666666667, height = H*0.01041666666667)}
    else{pdf(file = Save.plot)}}
  
  #### Data pulishing ####
  print("Lets PCA bioclimatic !")
  Remove.name <- c("Biome", "Biome.no", "Ecosystem", "Latitude", "Longitude", "Bioclim",
                   "Aridity", "Aridity2", "Type", "AP_NAP", "PFT", "GrowthForm", "Country")
  Keep.xdata <- MP[intersect(names(MP), Remove.name)]
  MP <- MP[setdiff(names(MP), Remove.name)]
  names(MP) <- gsub("TRY_", "", names(MP))
  
  #### Pour les GDGTs ####
  if(GDGT == T){
    # if(Remove.7Me == T){MP <- MP[setdiff(names(MP), names(MP)[grepl("7Me", names(MP))]),]}
    names(MP) <- gsub("f.", "", names(MP))
    names(MP) <- gsub("_5Me", "", names(MP))
    names(MP) <- gsub("_6Me", "\\'", names(MP))
    names(MP) <- gsub("_7Me", "\\''", names(MP))
  }
  #### Groupes ####
  if(is.null(Groupes) == F){
    Groupes <- melt(Groupes)
    Groupes <- Groupes$L1[match(names(MP), Groupes$value)]
    Groupes[is.na(Groupes)] <- "Unknown"
    Groupes <- factor(Groupes)
  }
  
  #### Transforming the data + PCA calcul ####
  if(nlevels(as.factor(is.na(MP))) >= 2){
    library(missMDA)
    nb <- estim_ncpPCA(MP, ncp.max = 5) ## Time consuming, nb = 2
    MP.comp <- imputePCA(MP, ncp = nb[[1]])
    MP <- MP.comp$completeObs
  }
  
  Sites.to.rm <- names(which(sapply(MP, function(x)all(is.na(x)))))
  if(length(Sites.to.rm) > 0){MP <- MP[!names(MP) %in% Sites.to.rm]}
  
  if(Helinger.trans == T){MP <- vegan::decostand(MP, method = "hellinger")}
  
  MP.pca <- PCA(MP, graph = FALSE, scale.unit = transp_OK)
  PCA <- data.frame(MP.pca$ind$coord)
  PCA <- PCA[,1:2]
  names(PCA) <- c("PC1", "PC2")
  
  if(is.null(Keep.xdata$Type) == F){
    PCA <- cbind(PCA, Type = Keep.xdata$Type, Biome = Keep.xdata[[Cluster.core]])
    PCA <- PCA[PCA$Type == "MPS",]
    if(is.null(Shape) == T){Shape <- 21}
    Point.surf <- geom_point(inherit.aes = F, PCA, mapping = aes(x = PC1, y = PC2, fill = Biome), shape = Shape, colour = "grey10", alpha = Alpha.MPS, na.rm = T, size = Size.MPS)
    Opac <- 0.3
    Ptaille <- 2
  }
  else{
    Point.surf <- NULL
    Opac <- 0.7
    Ptaille <- 1.5 
  }
  if(is.null(Dot.size) == F){Ptaille <- Dot.size}
  if(is.null(Dot.opac) == F){Opac <- Dot.opac}
  
  #### Color settings ####
  if(is.null(Groupes) == F){
    my_orange = c("Water" = "darkblue",
                  "Altitude" = "brown",
                  "Temperature" = "darkred")
    
    Scale.color.vec <- scale_color_manual(values = my_orange, name = "Proxies")
    Show.leg.arrow = T
    
  }
  else{
    Groupes <- "Unique"
    Show.leg.arrow = F
    my_orange <- data.frame(Unique = "royalblue")
    Scale.color.vec <- scale_color_manual(values = my_orange, guide = "none")
  }
  
  if(is.null(Cluster.core) == F){
    if(is.numeric(Keep.xdata[[Cluster.core]]) == T){
      my_orange2 = brewer.pal(n = 11, "RdYlBu")[Keep.col2[-c(3,4,5,7,8,9)]] 
      orange_palette2 = colorRampPalette(my_orange2)
      my_orange2 = rev(orange_palette2(length(seq(min(Keep.xdata[[Cluster.core]]), max(Keep.xdata[[Cluster.core]]), by = 200))))
      Scale.fill <- scale_fill_gradientn(colours = my_orange2, guide = "colourbar", 
                                         name = Cluster.core.lab,
                                         breaks = seq(round(min(Keep.xdata[[Cluster.core]]),digits = -3), max(Keep.xdata[[Cluster.core]]), by = 1000),
                                         na.value = "white")}
    else{
      values.bi = c("Deserts & Xeric Shrublands" = "#C88282",
                    "Temperate Grasslands, Savannas & Shrublands" = "#ECED8A",
                    "Montane Grasslands & Shrublands" = "#D0C3A7",
                    "Temperate Conifer Forests" = "#6B9A88",
                    "Temperate Broadleaf & Mixed Forests" = "#3E8A70",
                    "N/A" = "#FFEAAF",
                    "Tundra" = "#A9D1C2",
                    "Boreal Forests/Taiga" = "#8FB8E6",
                    "Tropical & Subtropical Coniferous Forests" = "#99CA81",
                    "Mangroves" = "#FE01C4",
                    "Flooded Grasslands & Savannas" = "#BEE7FF",
                    "Tropical & Subtropical Moist Broadleaf Forests" = "#38A700",
                    "Plant_height" = "royalblue",
                    "Leaf_thickness" = "darkorange",
                    "Photosynthesis_pathway" = "purple",
                    "TUSDB sites" = "#323232",
                    "Woodyness" = "#323232",
                    "Leaf_size" = "darkred",
                    "Variable" = "#aa373aff",
                    "Algal" = "#4666E9",
                    "NAP" = "#b5ab32ff",
                    "Herb" = "#b5ab32ff",
                    "Shrub" = "#aa373aff",
                    "Other" = "grey90",
                    "Unknown" = "grey90",
                    "AP" = "#0f6b31ff",
                    "Tree" = "#0f6b31ff",
                    "1_Hyper-arid" = "#8c510a", "2_Arid" = "#bf812e", "3_Semi-arid" = "#dfc27e", "4_Dry sub-humid" = "#f5e9bf", "5_Humid" = "#80cec1",
                    "1. Hyper-arid" = "#8c510a", "2. Arid" = "#bf812e", "3. Semi-arid" = "#dfc27e", "4. Dry sub-humid" = "#f5e9bf", "5. Humid" = "#80cec1",
                    Mongolia = "#3e96bdff", Chine = "#f02a26", Uzbekistan = "#6fb440", Armenia = "#54a697", "China, Tibet" = "#8c510a", "Northern Iran" = "#176E5B",
                    Tajikistan = "#e4af08", Russia = "#0035a9", Azerbaijan = "#094227", China = "#bb0202", "ACA lakes" =  "purple",
                    "Ch'ol cold desert-steppes" = "#7916C4", 
                    "Tugai riparian forest" = "#BB0268", 
                    "Ch'ol warm deserts" = "#bb0202", 
                    "Adyr desert-steppes" = "#ff5400", 
                    "Adyr steppes" = "#e6c607", 
                    "Tau riparian forest" = "#2C9740", 
                    "Tau thermophilous woodlands" = "#85682D", 
                    "Tau juniper steppe-forest" = "#176E5B",
                    "Tau steppes" = "#bab133",
                    "Alau cryophilous steppe-forest" = "#54a697",
                    "Alau meadows" = "#197CDA",
                    "Tugai riparian forests" = "#7916C4", 
                    "Ch'ol deserts" = "#bb0202", 
                    "Adyr pseudosteppes" = "#ff5400", 
                    "Adyr steppes" = "#ECED8A", 
                    "Tau xeric shrublands" = "#85682D", 
                    "Tau open woodlands" = "#1e8736",
                    "Tau steppes" = "#b0a62e",
                    "Alau cryophilous open woodlands" = "#54a697",
                    "Alau mesic grasslands" = "#197CDA" 
      )
      
      
      
      values.bi <- values.bi[which(names(values.bi) %in% unique(Keep.xdata[[Cluster.core]]))]
      Scale.fill <- scale_fill_manual(values = values.bi, name = Cluster.core.lab, drop = T)
      Scale.color <- scale_color_manual(values = values.bi, name = Cluster.core.lab, drop = T, guide = "none")
    }
    Col.select <- Keep.xdata[[Cluster.core]]
    Fill.select <- Keep.xdata[[Cluster.core]]
  }
  else{
    Scale.fill <- scale_fill_manual(values = "brown", name = "Species density")
    Scale.color <- scale_color_manual(values = "brown", name = "Species density")
    Opac = 0.2
    Col.select = "brown"
    Fill.select = "brown"
  }
  
  if(Marg.density.plot == F){
    My_title <- paste(Num.facet, Site.name, Type.samples, sep = " ")
  }
  else{My_title <- NULL}
  
  if(is.null(Leg.nrow) == F){
    Guide.color <- guides(fill = guide_legend(nrow = Leg.nrow), color = guide_legend(nrow = Leg.nrow))
    # Guide.color <- guide_legend(override.aes = list(nrow = Leg.nrow), title.position = "top")
    
    # guides(color = guide_legend(override.aes = list(size = Doz.size.leg, alpha = 1), title.position = "top"))+
    # print(Leg.nrow)
  }
  else{Guide.color <- NULL}
  
  if(is.null(Doz.size.leg)){Doz.size.leg <- Ptaille+2}
  
  #### Density contours #### 
  if(Density.contour == T){
    Data.contour <- data.frame(MP.pca$ind$coord)
    Data.contour$Col.select <- Col.select
    Data.contour <- unique(Data.contour)
    if(Density.type == "contour"){
      Density.color.line <- "grey20"
      Scale.opa <- NULL
      Point.surf <- NULL
      Dot.up <- "point"
    }
    if(Density.type == "polygon"){
      Density.color.line <- NA
      Scale.opa <- scale_alpha_continuous(range = Opa.range)
      if(is.null(Shape) == T){Shape <- 16}
      Point.surf <- geom_point(inherit.aes = F, Data.contour, mapping = aes(x = Dim.1, y = Dim.2, colour = Col.select), shape = Shape, alpha = Opac, na.rm = T, size = Ptaille)
      Dot.up <- "none"
    }
    Density.contour <- stat_density_2d(data = Data.contour, mapping = aes(x = Dim.1, y = Dim.2, alpha = ..level..),
                                       geom = Density.type, colour = Density.color.line, show.legend = F,
                                       bins = 6)
    
  }
  else{
    Dot.up <- "point"
    Scale.opa <- NULL
    Density.contour <- NULL}
  #### Labels ####
  if(Show.site.lab == T){
    PCA$Site <- row.names(PCA)
    Site.lab <- geom_text(inherit.aes = F, PCA, mapping = aes(x = PC1, y = PC2, label = Site), alpha = 0.5, na.rm = T, size = Size.MPS)
  }
  else{Site.lab <- NULL}
  
  
  #### Reverse dim ####
  A = round(MP.pca$eig[1,2], digits = 0)
  B = round(MP.pca$eig[2,2], digits = 0)
  
  if(Reverse.dim == T){
    if(Show.annot == T){
      Note.n <- annotate("text", x = min(min(MP.pca$ind$coord[,2]), min(MP.pca$var$coord[,2])), y = min(MP.pca$ind$coord[,1]), label = paste("n = ", nrow(MP), sep = ""), size = 5, hjust = 0)}
    else{Note.n <- NULL}
    Axes <- c(2,1)
    X.arrow <- c(2)
    Y.arrow <- c(1)
    Xlab <- labs(x = substitute(paste("PCA"[2], ~ "(", B, " %)", sep = " " )),
                 y = substitute(paste("PCA"[1], ~ "(", A, " %)", sep = " " )))
  }
  else{
    if(Show.annot == T){
      # Note.n <- annotate("text", x = min(min(MP.pca$ind$coord[,1]), min(MP.pca$var$coord[,1])),  y = min(min(MP.pca$ind$coord[,2]), min(MP.pca$var$coord[,2])), label = paste("n = ", nrow(MP), sep = ""), size = 5, hjust = 0, vjust = 0.5)}
      Note.n <- annotate("text", x = Inf, y = -Inf, label = paste("n = ", nrow(MP), sep = ""), size = 4, hjust = 1.1, vjust = -0.5)}
    else{Note.n <- NULL}
    
    Axes = c(1,2)
    X.arrow <- c(1)
    Y.arrow <- c(2)
    Xlab <- labs(x = substitute(paste("PCA"[1], ~ "(", A, " %)", sep = " " )),
                 y = substitute(paste("PCA"[2], ~ "(", B, " %)", sep = " " )))
  }
  if(is.null(Manu.lim.x)==F){
    Lim.x <- xlim(Manu.lim.x)
    if(Show.annot == T){Note.n <- annotate("text", x = Manu.lim.x[1], y = min(MP.pca$ind$coord[,2]), label = paste("n = ", nrow(MP), sep = ""), size = 5, hjust = 0, vjust = 0.5)}
    else{Note.n <- NULL}  
  }
  else{Lim.x <- NULL}
  if(is.null(Manu.lim.y)==F){Lim.y <- ylim(Manu.lim.y)}
  else{Lim.y <- NULL}
  
  #### Ajout vectors ####
  Pouet <- data.frame(X0 = 0,
                      Y0 = 0,
                      X1 = MP.pca$var$coord[,X.arrow]*Scale.PCA,
                      Y1 = MP.pca$var$coord[,Y.arrow]*Scale.PCA,
                      Lab = row.names(MP.pca$var$coord),
                      Contrib = MP.pca$var$contrib[,c(1)] + MP.pca$var$contrib[,c(2)] , 
                      Groupes = Groupes)
  Pouet$Contrib <- Pouet$Contrib/sum(Pouet$Contrib)
  
  Arrow.lab <- geom_text_repel(data = Pouet, aes(x = X1, y = Y1, label = Lab, color = Groupes), min.segment.length = 10, force = 20, show.legend = F)
  if(Contrib == T & Show.arrow == T){
    Arrow <- geom_segment(data = Pouet, aes(x = X0, y = Y0, xend = X1, yend = Y1, color = Groupes, size = Contrib), arrow = arrow(length=unit(0.2,"cm")), show.legend = Show.leg.arrow)}
  if(Contrib == F & Show.arrow == T){
    Arrow <- geom_segment(data = Pouet, aes(x = X0, y = Y0, xend = X1, yend = Y1, color = Groupes, size = Contrib), arrow = arrow(length=unit(0.2,"cm")), show.legend = Show.leg.arrow)}
  if(Show.arrow == F){
    Arrow <- NULL
    Arrow.lab <- NULL
  }
  Scale.size <- scale_size_continuous(range = c(0.2,1), guide = "none")
  
  #### Ajout symbole ####
  if(is.null(Symbol.path) == F){
    if(is.null(Symbol.pos)== T){Symbol.pos <- c(.9, .9, .16)}
    if(grepl("\\.png", Symbol.path)){
      library(png)
      library(grid)
      img <- readPNG(Symbol.path)
      g <- rasterGrob(x = Symbol.pos[1], y = Symbol.pos[2], width = Symbol.pos[3], height = Symbol.pos[3], img, interpolate = T)
    }
    
    if(grepl("\\.xml", Symbol.path)){
      library(grImport)
      img <- readPicture(Symbol.path)
      g <- pictureGrob(x = Symbol.pos[1], y = Symbol.pos[2], width = Symbol.pos[3], height = Symbol.pos[3], img)
    }
    
    Logo <- annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)
  }
  else{Logo <- NULL}
  
  #### PLOT  ####
  p <- fviz_pca_biplot(MP.pca, na.rm = T,
                       geom.ind = Dot.up, axes = Axes,
                       pointshape = 16,
                       pointsize = Ptaille,
                       fill.ind = Fill.select,
                       alpha.ind = Opac, #alpha.var ="contrib",
                       col.ind = Col.select,
                       repel = T, arrowsize = 0.5,
                       invisible = Centroide, # enlève ou ajoute le centroïde
                       addEllipses = Ellipse, ellipse.level = 0.75, ellipse.type = "norm", # convex or norm
                       ellipse.alpha = Ellipse.opa, 
                       title = My_title,
                       # col.var = Groupes #a enleve parfois
                       col.var = NA
  ) +
    guides(color = guide_legend(override.aes = list(size = Doz.size.leg, alpha = 1), title.position = "top"))+
    Guide.color +
    Xlab+ #Ylab+ 
    Lim.x+ Lim.y+
    # Axes+
    Scale.opa + 
    Density.contour +
    Scale.fill +
    Scale.color +
    Point.surf + Site.lab +
    new_scale_color() + 
    Arrow + Arrow.lab + 
    Scale.size + Scale.color.vec +
    Note.n + Logo + 
    guides(fill = guide_legend(title.position = "top")) +
    guides(color = guide_legend(title.position = "top")) +
    theme(plot.background = element_blank(),
          legend.position = Legend.position,
          legend.title = element_text(size = (Legend.size+2)),
          legend.text = element_text(size = Legend.size),
          legend.key = element_blank(),
          plot.margin=unit(c(0,0,0,0),"cm"),
          panel.grid = element_blank(),
          panel.border = element_rect(NA, "black", linewidth = 1),
          # panel.grid = element_line(linetype = "dashed"),
          axis.line = element_blank())
  
  #### Add margin density ####
  if(Marg.density.plot == T){
    MMarg <- data.frame(MP.pca$ind$coord)
    MMarg <- cbind(MMarg, Biome = Col.select)
    x_limits <- ggplot_build(p)$layout$panel_scales_x[[1]]$range$range
    y_limits <- ggplot_build(p)$layout$panel_scales_y[[1]]$range$range
    List.of.NA <- which(MMarg$Dim.1 < 1e-12 & MMarg$Dim.1 > -1e-12 & MMarg$Dim.2 < 1e-12 & MMarg$Dim.2 > -1e-12)
    if(length(List.of.NA) > 0){
      print("Remove NA from density.")
      MMarg <- MMarg[-List.of.NA,]}
    
    #### Density plots up ####
    plot_top <- ggplot(MMarg, aes(x = Dim.1, fill = Biome)) + 
      geom_density(alpha = 0.6, size = 0.1) + Scale.fill + 
      ggtitle(paste(Num.facet, Site.name, Type.samples, sep = " "))+ 
      scale_x_continuous(limits = x_limits, expand = c(0,0))+
      #### Theme ####
    theme(
      axis.line = element_blank(),
      axis.text.x = element_blank(), axis.text.y = element_text(hjust = 1, size = 6),
      axis.ticks.x.bottom = element_blank(),
      axis.title = element_blank(),
      axis.line.y = element_line(colour = "grey"),
      axis.ticks.y = element_line(colour = "grey"),
      #panel.border = element_rect(fill = NA, colour = "grey"),
      legend.title = element_text(),
      legend.key = element_blank(),
      legend.justification = c("center"),               # left, top, right, bottom
      legend.text = element_text(size = 8),
      panel.background = element_blank(),
      panel.spacing = unit(0.7, "lines"),
      legend.position = "none",
      strip.text.x = element_text(size = 12, angle = 0, face = "bold"),
      strip.placement = "outside",
      strip.background = element_rect(color = "white", fill = "white"),
      plot.margin=unit(c(0,0,0,0),"cm")
    )
    
    #### Density plots right ####
    plot_right <- ggplot(MMarg, aes(x = Dim.2, fill = Biome)) + 
      geom_density(alpha = 0.6, size = 0.1) + Scale.fill +
      scale_x_continuous(limits = y_limits, expand = c(0,0))+
      coord_flip() + 
      #### Theme ####
    theme(
      axis.line.y = element_blank(),
      axis.text.y = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
      axis.ticks.y = element_blank(),
      axis.title = element_blank(),
      axis.line.x = element_line(colour = "grey"),
      axis.ticks.x = element_line(colour = "grey"),
      #panel.border = element_rect(fill = NA),
      legend.title = element_text(),
      legend.key = element_blank(),
      legend.justification = c("center"),               # left, top, right, bottom
      legend.position = "none",
      legend.text = element_text(size = 8),
      panel.background = element_blank(),
      panel.spacing = unit(0.7, "lines"),
      strip.text.x = element_text(size = 12, angle = 0, face = "bold"),
      strip.placement = "outside",
      strip.background = element_rect(color = "white", fill = "white"),
      plot.margin=unit(c(0,0,0,0),"cm")
    )
    
    
    #### Export ####
    layout <- "AAAAAAA#
             CCCCCCCB
             CCCCCCCB
             CCCCCCCB
             CCCCCCCB
             CCCCCCCB
             "
    p <- plot_top + plot_right + p + plot_layout(design = layout)
  }
  
  print(p)
  
  #### Save html ####
  if(Show.Plotly == T){
    library(plotly)
    library(htmlwidgets)
    
    #### 3D PCA chart ####
    Pouet <- data.frame(MP.pca$ind$coord)
    Pouet$Col <- Col.select
    fig <- plot_ly(Pouet, x = ~Dim.1, y = ~Dim.2, z = ~Dim.3, color = ~Col, 
                   colors = values.bi, text = ~paste('Site: ', row.names(Pouet), "\n Biome:", Pouet$Col, sep = "")) #%>%
    options(warn = - 1)
    # add_markers(size = 60)
    # 
    # fig <- fig %>%
    #   layout(
    #     # title = tit,
    #     xaxis = list(title = "PCA1"),
    #     scene = list(bgcolor = "#e5ecf6")
    #   )
    
    #### Export ####
    Save.plot.html <- gsub("pdf", "html", Save.plot)
    Keep.name <- gsub(".*\\/", "", Save.plot.html)
    Path.root <- paste(gsub(Keep.name, "", Save.plot.html), "HTML_files/", sep = "")
    if(file.exists(Path.root) == F){dir.create(Path.root)}
    Save.plot.html <- paste(Path.root, Keep.name, sep = "")
    # # p1_ly <- ggplotly(p)
    # # p1_ly <- p1_ly %>% layout()
    options(warn = - 1)
    saveWidget(fig, file = Save.plot.html)
    options(warn = - 1)
    # saveWidget(p1_ly, file = Save.plot.html)
  }
  
  #### Export data ####
  if(is.null(Save.path) == F){
    Site.name <- gsub(" ","_",Site.name)
    Save.path.Site <- gsub("\\.csv", "_PCA_site.csv", Save.path)
    Save.path.Taxon <- gsub("\\.csv", "_PCA_clim_param.csv", Save.path)
    write.table(data.frame(MP.pca$ind$coord), file = Save.path.Site, col.names = T, sep = ",")
    write.table(data.frame(MP.pca$var$coord), file = Save.path.Taxon, col.names = T, sep = ",")}
  if(is.null(Save.plot) == F){
    dev.off()}
  
  if(return.pick == T){return(p)}
  if(return.pick == F){
    if(PCA.site == T){
      Mat.ret <- data.frame(MP.pca$ind$coord)
      # print(Mat.ret)
      return(Mat.ret)
    }
    if(Contrib == T){return(MP.pca$var$contrib)}
  }
}

PCA.pollen.surf <- function(MP, transp_OK, Helinger.trans = F, Alpha.dot = NULL, Site.name, Type.samples, Short.name, GDGT, Annot, Nb.contrib = NULL, Leg.size = 1,
                            Remove.7Me, Display.legends, Simple.title, Show.text, Symbol.loc, Symbol.path, Vector.show = NULL, Display.plot = T, Result.vector = F,
                            Cluster.path, Cluster.groups, Sort.eco, Csv.sep, Scale.PCA, Color.choice = NULL, Dot.size = 1.5, Reverse.x = F, Reverse.y = F,
                            Vector.scale = 1, cluster.from.PCA, Leg.loc,Leg.loc2, Save.path, Save.plot = NULL, Manu.lim, Title.inside = F){
  #### Settings ####
  library(vegan)
  Clustering = T
  if(missing(cluster.from.PCA)){cluster.from.PCA = F}
  if(missing(Show.text)){Show.text = F}
  if(missing(Cluster.path)){Clustering = F}
  if(missing(Short.name)){Short.name = T}
  if(missing(Remove.7Me)){Remove.7Me = F}
  if(missing(Display.legends)){Display.legends = T}
  if(missing(Simple.title)){Simple.title = F}
  if(missing(Cluster.groups)){Cluster.groups = "Ecosystems"}
  if(missing(Csv.sep)){Csv.sep = "\t"}
  if(missing(Scale.PCA)){Scale.PCA = 1}
  if(missing(Annot)){Annot = NULL}
  if(missing(Symbol.path)){Symbol.path = NULL}
  if(missing(Symbol.loc)){Symbol.loc = NULL}
  if(missing(Save.path)){Save.path = NULL}
  if(missing(Sort.eco)){Sort.eco = NULL}
  if(missing(Manu.lim)){Manu.lim = NULL}
  if(missing(Type.samples)){Type.samples = NULL}
  if(missing(Leg.loc)){Leg.loc = "bottomleft"}
  if(missing(Leg.loc2)){Leg.loc2 = "bottomright"}
  if(missing(Site.name)){Site.name = ""}
  if(missing(GDGT)){GDGT = NULL}
  
  #### Save plots ####
  if(is.null(Save.plot) == F){
    Path.to.create <- gsub("(.*/).*\\.pdf.*","\\1", Save.plot)
    dir.create(file.path(Path.to.create), showWarnings = FALSE)
    if(is.null(W) == F & is.null(H) == F){
      pdf(file = Save.plot, width = W*0.01041666666667, height = H*0.01041666666667)}
    else{pdf(file = Save.plot)}}
  
  #### Pour les GDGTs ####
  if(is.null(GDGT) == F){
    if(Remove.7Me == T){MP <- MP[setdiff(row.names(MP), row.names(MP)[grepl("7Me", row.names(MP))]),]}
    row.names(MP) <- gsub("f.", "", row.names(MP))
    row.names(MP) <- gsub("_5Me", "", row.names(MP))
    row.names(MP) <- gsub("_6Me", "\\'", row.names(MP))
    row.names(MP) <- gsub("_7Me", "\\''", row.names(MP))
  }
  
  #### Pour le pollen ####
  if(is.null(GDGT) == T){
    row.names(MP) <- gsub("aceae",".",row.names(MP))
    row.names(MP) <- gsub(".undiff","",row.names(MP))
    row.names(MP) <- gsub(".indet","",row.names(MP))
    row.names(MP) <- gsub(".*Po.$","Poa.",row.names(MP))
    row.names(MP) <- gsub(".spp.","",row.names(MP))
    row.names(MP) <- gsub(".type","",row.names(MP))
  }
  
  MP <- data.frame(t(MP), check.names = F)
  
  #### Transforming the data ####
  if(transp_OK == FALSE){
    MP.pca <- rda(MP, scale = F)
    MTitle = paste("PCA", Type.samples, "-", Site.name)}
  else{
    if(Helinger.trans == T){
      MP <- vegan::decostand(MP, method = "hellinger")
      MP.pca <- rda(MP, scale = F)}
    else{
      MP.pca <- rda(MP, scale = T)}
    
    MTitle = paste("PCA", Type.samples, "log-trans, ", Site.name)}
  if(Simple.title == T){MTitle = paste("PCA", Type.samples)}
  if(is.null(Annot) == F){MTitle <- paste(Annot, MTitle)}
  if(Title.inside == T){MTitle <- NULL}
  
  #### Test clustering from PCA ####
  if(cluster.from.PCA == T){
    # print(names(MP.pca))
    # View(MP.pca)
    library(FactoMineR)
    library(factoextra)
    # Compute PCA with ncp = 3
    res.pca <- PCA(MP, ncp = 5, graph = FALSE, scale.unit = F)
    # Compute hierarchical clustering on principal components
    res.hcpc <- HCPC(res.pca, graph = FALSE)
    # fit <- kmeans(pf$ind$coord[,1:4],2)
    # print(res.hcpc$data.clust)
    
    
    pdend <- fviz_dend(res.hcpc, 
                       cex = 0.7,                     # Label size
                       palette = "jco",               # Color palette see ?ggpubr::ggpar
                       rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
                       rect_border = "jco",           # Rectangle color
                       labels_track_height = 0.8      # Augment the room for labels
    )
    
    pclu <- fviz_cluster(res.hcpc,
                         repel = TRUE,            # Avoid label overlapping
                         show.clust.cent = TRUE, # Show cluster centers
                         palette = "jco",         # Color palette see ?ggpubr::ggpar
                         ggtheme = theme_minimal(),
                         main = "Factor map"
    )
    print(pdend)
    print(pclu)}
  
  #### Graphic parameters ####
  site.score  <- scores(MP.pca, choices=c(1,2), display = "wa", scaling = Scale.PCA)         # valeur PCA age
  taxa.score <- scores(MP.pca, choices=c(1,2), display = "species", scaling = Scale.PCA)
  
  par(mgp = c(1.7,0.6,0), mar=c(3,3,2,0.3)+0.1)
  
  PC1.MP <- round(MP.pca[["CA"]][["eig"]][["PC1"]]/MP.pca[["tot.chi"]]*100, digits = 0)
  PC2.MP <- round(MP.pca[["CA"]][["eig"]][["PC2"]]/MP.pca[["tot.chi"]]*100, digits = 0)
  xmin <- min(c(min(site.score[,1]),min(taxa.score[,1]*Vector.scale)))                             # lim axe PC1
  xmax <- max(c(max(site.score[,1]),max(taxa.score[,1]*Vector.scale)))
  ymin <- min(c(min(site.score[,2]),min(taxa.score[,2]*Vector.scale)))                             # lim axe PC2
  ymax <- max(c(max(site.score[,2]),max(taxa.score[,2]*Vector.scale)))
  
  xmin <- xmin + .05*xmin
  xmax <- xmax + .05*xmin
  ymin <- ymin + .05*ymin
  ymax <- ymax + .05*ymax
  
  if(is.null(Manu.lim) == F){
    xmin <- Manu.lim[1]
    xmax <- Manu.lim[2]
    ymin <- Manu.lim[3]
    ymax <- Manu.lim[4]
  }
  
  #### Reverse axis ####
  if(Reverse.x == T){
    taxa.score[,1] <- -taxa.score[,1]
    site.score[,1] <- -site.score[,1]}
  if(Reverse.y == T){
    taxa.score[,2] <- -taxa.score[,2]
    site.score[,2] <- -site.score[,2]}
  
  #### Contribution ####
  if(is.null(Vector.show) == F & is.null(Nb.contrib) == T){
    taxa.score <- taxa.score[row.names(taxa.score) %in% Vector.show,]
    MP.pca$CA$v <- MP.pca$CA$v[row.names(MP.pca$CA$v) %in% Vector.show,]
  }
  
  if(is.null(Nb.contrib) == F){
    Contribution <- sqrt(abs(taxa.score[,1])^2 + abs(taxa.score[,2])^2)
    Contribution <- Contribution[order(Contribution, decreasing = T)]
    Contribution <- names(Contribution[1:Nb.contrib])
    if(is.null(Vector.show) == F){
      Contribution <- unique(c(Contribution, Vector.show))
      print(paste("**** We will merge the", Nb.contrib, "most contributing taxa + the", length(Vector.show), "manually selected taxa.****" ))}
    
    taxa.score <- taxa.score[row.names(taxa.score) %in% Contribution,]
    MP.pca$CA$v <- MP.pca$CA$v[row.names(MP.pca$CA$v) %in% Contribution,]
  }
  
  
  #### Plot ####
  if(Display.plot == T | is.null(Save.plot) == F){
    #### Vector plot ####
    biplot(MP.pca,                                              # ajout des vecteurs taxa
           display = "species",                                  # 'species', vecteurs taxa
           type = "n",                                           # 't', txt aux vecteurs, 'n', vecteurs vides
           cex = 2,
           main = MTitle,
           xlim = c(xmin,xmax), 
           ylim = c(ymin,xmax),
           xlab = bquote(PCA[1] ~ "(" ~ .(format(PC1.MP, digits = 2)) ~ "%)"),
           ylab = bquote(PCA[2] ~ "(" ~ .(format(PC2.MP, digits = 2)) ~ "% )"))
    
    #### Color points setting ####
    if(Clustering == T){
      if(is.character(Cluster.path) == T){Cluster <- data.frame(read.csv(Cluster.path, sep = Csv.sep ,dec=".",header=T,row.names=1), stringsAsFactors = T)}
      if(is.data.frame(Cluster.path) == T){Cluster <- Cluster.path}
      Inter <- intersect(row.names(Cluster), row.names(site.score))
      Cluster <- Cluster[which(row.names(Cluster) %in% Inter),]
      M.eco <- subset(Cluster, select = Cluster.groups)
      M.eco <- data.frame(M.eco[][row.names(site.score),])
      row.names(M.eco) <- row.names(site.score)
      Fact.eco <- as.factor(M.eco[[1]])
      if(is.null(Sort.eco) == F){Fact.eco <- ordered(Fact.eco, levels = Sort.eco)}
      num_colors <- nlevels(Fact.eco)
      
      if(is.null(Color.choice) == F){diamond_color_colors <- Color.choice}
      else{
        colfunc    <- colorRampPalette(c("firebrick3", "darkorange", "goldenrod1", "#38A700", "darkgreen", "dodgerblue3", "grey10"))
        diamond_color_colors <- colfunc(num_colors)
      }
      Col.eco <- diamond_color_colors[Fact.eco]
    }
    else{Col.eco = 1}
    
    if(is.null(Alpha.dot) == F){
      Col.eco <- adjustcolor(Col.eco, alpha.f = Alpha.dot)
    }
    
    points(site.score, pch = 19, cex = Dot.size, col = Col.eco)
    # text(MP.pca, scaling = Scale.PCA, display = "species", cex = .95, col = "#6e2115ff")
    text(taxa.score*Vector.scale, colnames(MP), cex=.95, pos = 3, col = "#6e2115ff")                # textes taxa
    
    arrows(0,0, x1 = taxa.score[,1]*Vector.scale, y1 = taxa.score[,2]*Vector.scale, length = 0.05, lty = 1, col="#6e2115ff")
    if(Title.inside == T){
      usr <- par("usr")   # save old user/default/system coordinates
      par(usr = c(0, 1, 0, 1)) # new relative user coordinates
      text(0.01, 0.95, Annot, adj = 0, cex = 1.7)  # if that's what you want
      par(usr = usr) # restore original user coordinates
    }
    #### Add plots ####
    if(Short.name == T){Tsite <- sub("M","",row.names(MP))}
    else{Tsite <- row.names(MP)}
    if(Show.text == T){text(site.score, Tsite, cex=.6, pos= 3)}  # textes sites
    
    #### Legend #####
    if(Clustering == T & Display.legends == T){
      legend(Leg.loc,
             legend = levels(Fact.eco),
             col = diamond_color_colors,
             pch = 19, cex = Leg.size, pt.cex = Leg.size,
             y.intersp = 0.85,	          # espace entre y
             x.intersp = 0.6,           # espace entre x
             bty = "n" )
    }
    
    #### Legend 2 ####
    legend(Leg.loc2,
           legend = paste("n = ", nrow(MP)),
           pch = NA,
           y.intersp = 0.7,	                                     # espace entre y
           x.intersp = 0.4,                                      # espace entre x
           bty = "n" )
    
    #### Add symbol ####
    if(is.null(Symbol.path) == F){
      if(is.null(Symbol.loc)== T){Symbol.loc <- c(0.83,0.83,0.99,0.99)}
      if(grepl("\\.png", Symbol.path)){
        library(png)
        library(grid)
        pic <- readPNG(Symbol.path)
        usr <- par("usr")
        par(usr = c(0, 1, 0, 1))
        rasterImage(pic,Symbol.loc[1], Symbol.loc[2], Symbol.loc[3], Symbol.loc[4])
        par(usr = usr)
      }
      
      if(grepl("\\.xml", Symbol.path)){
        library(grImport)
        pic <- readPicture(Symbol.path)
        usr <- par("usr")
        # print(Symbol.loc[1])
        par(usr = c(0, 1, 0, 1))
        picture(pic,Symbol.loc[1], Symbol.loc[2], Symbol.loc[3], Symbol.loc[4])
        par(usr = usr)
      }
    }
    
  }
  
  if(is.null(Save.plot) == F){dev.off()}
  
  #### Export data ####
  if(is.null(Save.path) == F){
    Site.name <- gsub(" ","_",Site.name)
    Save.path.Site <- gsub("\\.csv", "_PCA_site.csv", Save.path)
    Save.path.Taxon <- gsub("\\.csv", "_PCA_taxa.csv", Save.path)
    Save.path.Cluster <- gsub("\\.csv", "_PCA_HCPC.csv", Save.path)
    
    write.table(t(site.score), file = Save.path.Site, col.names = FALSE, sep = ",")
    write.table(t(taxa.score), file = Save.path.Taxon, col.names = TRUE, sep = ",")
    
    if(cluster.from.PCA == T){
      
      write.table(res.hcpc$data.clust, file = Save.path.Cluster, col.names = FALSE, sep = ",")
      
    }
  }
  if(Result.vector == T){return(taxa.score)}
  else{return(site.score)}
}

# This function plot the boxplot for each br-GDGT fractional abundance
# If you add Mtype, it make the groups for moss, soils and sediment...
# Graph présent chez Ding et al. 2015, Hopmans et al. 2004
# Si ajout iso.GDGT, Hopmans et al. 2004 
GDGT.histo.plot.surf.core <- function(Mcore, Msurf, Mtype, Select.type, Show.Plotly, Reorder.group, Global.box = F, Annot.size = 6,
                                      Keep.br, Leg.nb.lines, Leg.iso = F, Box.linewidth = .5, Leg.size = 13,
                                      Remove.8Me, Remove.7Me, Color.choice, Iso.GDGT, Leg.pos, Return.plot = F, Boxplot.title = NULL, Leg.box = F,
                                      Zoom1.comp = NULL, Zoom2.comp = NULL, Zoom1.Ymax, Zoom2.Ymax, Insert1.loc = NULL, Insert2.loc = NULL,
                                      Name.untype, Ymax, Save.path, W, H, Dot.pop, Overlap.OK){
  #### Initialization values ####
  if(missing(Mtype)){Mtype = NULL}
  if(missing(Msurf)){Msurf = NULL}
  if(missing(Mcore)){Mcore = NULL}
  if(missing(Reorder.group)){Reorder.group = NULL}
  if(missing(Show.Plotly)){Show.Plotly = F}
  if(missing(Save.path)){Save.path = NULL}
  if(missing(W)){W = NULL}
  if(missing(H)){H = NULL}
  if(missing(Keep.br)){Keep.br = NULL}
  if(missing(Leg.nb.lines)){Leg.nb.lines = NULL}
  if(missing(Dot.pop)){Dot.pop = NULL}
  if(missing(Overlap.OK)){Overlap.OK = F}
  if(missing(Name.untype)){Name.untype = "Other"}
  if(missing(Select.type)){Select.type = "Sample.type"}
  if(missing(Color.choice)){Color.choice = NULL}
  if(missing(Remove.8Me)){Remove.8Me = F}
  if(missing(Remove.7Me)){Remove.7Me = F}
  if(missing(Iso.GDGT)){Iso.GDGT = F}
  if(missing(Leg.pos)){Leg.pos = c(0.31,.86)}
  
  #### Extract br-GDGT ####
  MBRcore <- Mcore[,grep("^f.I", colnames(Mcore))]
  MBRsurf <- Msurf[,grep("^f.I", colnames(Msurf))]
  
  #### Remove 8 Me ####
  if(Remove.8Me == T){
    MBRcore <- MBRcore[,!grepl("_8Me", colnames(MBRcore))]
    MBRsurf <- MBRsurf[,!grepl("_8Me", colnames(MBRsurf))]
  }
  
  #### Remove 7 Me ####
  if(Remove.7Me == T){
    MBRcore <- MBRcore[,!grepl("_7Me", colnames(MBRcore))]
    MBRsurf <- MBRsurf[,!grepl("_7Me", colnames(MBRsurf))]
  }
  
  if(is.null(MBRsurf) == F){MBRsurf <- MBRsurf[setdiff(names(MBRsurf), names(MBRsurf)[which(colSums(MBRsurf) == 0)])]}
  if(is.null(MBRcore) == F){MBRcore <- MBRcore[setdiff(names(MBRcore), names(MBRcore)[which(colSums(MBRcore) == 0)])]}
  
  #### Test if the names are egual ####
  if(is.null(MBRcore) == F & is.null(MBRsurf) == F){
    '%nin%' <- Negate('%in%')
    test1 <- length(names(MBRcore)[names(MBRcore) %nin% names(MBRsurf)])
    test2 <- length(names(MBRsurf)[names(MBRsurf) %nin% names(MBRsurf)])
    
    if(test1 > 0 & test2 > 0){
      print("Names are different.")
      break
    }
    if(test1 == 0 & test2 == 0){
      print("Let's go !")
      Mplot <- rbind(MBRcore, MBRsurf)   
    }
    else{print("Some molecules are missing !")
      print(setdiff(names(MBRsurf), names(MBRcore)))
      print(setdiff(names(MBRcore), names(MBRsurf)))
    }
  }
  else{
    if(is.null(MBRcore) == F){Mplot = MBRcore[,grep("^f.I", colnames(MBRcore))]}
    if(is.null(MBRsurf) == F){Mplot = MBRsurf[,grep("^f.I", colnames(MBRsurf))]}
  }
  
  #### Extract iso-GDGT ####
  if(Iso.GDGT == T){
    #### Recuperation des iso GDGT dans M tot ####
    if(is.null(Mcore) == F){
      Mcore.iso <- cbind(Mcore[,grep("^GDGT", colnames(Mcore))], Crenarch = Mcore$Crenarch, Crenarch.p = Mcore$Crenarch.p)
      if(is.null(Mcore.iso$GDGT0.Crenar) == F){Mcore.iso <- subset(Mcore.iso, select = -c(GDGT0.Crenar))} 
      Mcore.iso <- Mcore.iso/rowSums(Mcore.iso)*100
      
      #### Stats iso-GDGTs ####
    }
    
    if(is.null(Msurf) == F){
      Msurf.iso <- cbind(Msurf[,grep("^GDGT", colnames(Msurf))], Crenarch = Msurf$Crenarch, Crenarch.p = Msurf$Crenarch.p)
      if(is.null(Msurf.iso$GDGT0.Crenar) == F){Msurf.iso <- Msurf.iso[!grepl("GDGT0.Crenar", names(Msurf.iso))]}
      Msurf.iso <- Msurf.iso/rowSums(Msurf.iso)*100
    }
    
    if(is.null(Mcore) == F & is.null(Msurf) == F){
      Ymax.iso = max(max(Mcore.iso), max(Msurf.iso))
      L.iso = rbind(Msurf.iso, Mcore.iso)}
    if(is.null(Mcore) == T & is.null(Msurf) == F){
      Ymax.iso = max(Msurf.iso)
      L.iso = Msurf.iso}
    if(is.null(Mcore) == F & is.null(Msurf) == T){
      Ymax.iso = max(Mcore.iso)
      L.iso = Mcore.iso}
    #### Mise en forme des labels des axes ####
    if(is.null(Mtype) == F){
      Mtype <- Mtype[Select.type]
      names(Mtype) <- "Sample.type"
      A.iso <- merge(L.iso, Mtype, all.x = T, by = 0, sort = F)
      row.names(A.iso)<- A.iso$Row.names
      A.iso <- A.iso[,-1]
      levels(A.iso$Sample.type)[length(levels(A.iso$Sample.type)) + 1] <- Name.untype
      A.iso$Sample.type[is.na(A.iso$Sample.type)] <- Name.untype
    }
    else{ 
      A.iso <- L.iso 
      A.iso["Sample.type"]="Samples"}
    
    
    
    Model.lab <- gsub("GDGT", "GDGT-", names(A.iso))
    # Model.lab <- gsub("0", "0]", Model.lab)
    # Model.lab <- gsub("1", "1]", Model.lab)
    # Model.lab <- gsub("2", "2]", Model.lab)
    # Model.lab <- gsub("3", "3]", Model.lab)
    # Model.lab <- gsub("4", "4]", Model.lab)
    #Model.lab <- Model.lab[-1]                  # on enleve Row.names()
    Model.lab <- Model.lab[-length(Model.lab)]     # on enleve Samples.type
    Model.lab <- Model.lab[-length(Model.lab)]     # on enleve Crenarch.p et plus tard on rajoute Crenarch'
    
    #### If little pop to dot activated ####
    if(is.null(Dot.pop) == F){ # Echantillons choisis comme points
      Little.pop.iso <- A.iso[A.iso$Sample.type %in% Dot.pop,]
      for(i in 1:length(Dot.pop)){
        New.lab.iso <- unique(A.iso$Sample.type)[grepl(Dot.pop[i], unique(A.iso$Sample.type))]
        Little.pop.iso$Sample.type <- gsub(Dot.pop[i], New.lab.iso, Little.pop.iso$Sample.type)
        if(Overlap.OK == T){A.iso[A.iso$Sample.type %in% New.lab.iso, names(A.iso[,-length(names(A.iso))])] <- 100}}
      
      Not.in.dot.iso <- setdiff(unique(as.character(A.iso$Sample.type)), Dot.pop)
      New.row.0.iso <- data.frame(Sample.type = Not.in.dot.iso)
      Little.pop.iso <- rbind.fill(Little.pop.iso, New.row.0.iso)
      Little.pop.iso[is.na(Little.pop.iso)] <- 100
      A.little.iso <- melt(Little.pop.iso, id ='Sample.type')
      Add.points.size.lim.iso <- geom_dotplot(data = A.little.iso, aes(x=variable, y=value, fill=Sample.type),
                                              binaxis='y', stackdir='center', 
                                              position = position_dodge(.85),
                                              show.legend = FALSE,
                                              binwidth = 1)
      A.iso$Sample.type <- as.character(A.iso$Sample.type)
      A.iso <- melt(A.iso, id ='Sample.type')}
    
    else{                   # Pas échatillon comme points
      Add.points.size.lim.iso <- NULL
      A.iso$Sample.type <- as.character(A.iso$Sample.type)
      A.iso <- melt(A.iso, id ='Sample.type')
    }
  }
  else{A.iso = NULL}
  
  
  #### Label type ####
  Mplot <- 100*Mplot
  AZER <- gsub("f.", "", names(Mplot))
  AZER <- gsub("_5Me", "", AZER)
  AZER <- gsub("_6Me", "\\'", AZER)
  AZER <- gsub("_7Me", "\\''", AZER)
  AZER <- gsub("_8Me", "\\'''", AZER)
  names(Mplot)  <- AZER
  if(missing(Ymax)){Ymax = max(Mplot)}
  
  #### Merging datas ####
  if(is.null(Mtype) == F){
    if(Iso.GDGT == F){
      Mtype <- Mtype[Select.type]
      names(Mtype) <- "Sample.type"
    }
    A <- merge(Mplot, Mtype, all.x = T, by = 0, sort = F)
    row.names(A)<- A$Row.names
    A <- A[,-1]
    levels(A$Sample.type)[length(levels(A$Sample.type)) + 1] <- Name.untype
    A[is.na(A)] <- Name.untype
  }
  
  else{ 
    A = Mplot 
    A["Sample.type"]="Samples"}
  
  #### Counting n of each samples ####
  M.type <- as.character(A$Sample.type)
  N.type <- unique(M.type)
  NT <- sapply(N.type, function(x) length(M.type[M.type == x]))
  A.save <- A
  for(i in 1:length(NT)){A$Sample.type <- gsub(names(NT)[i], paste(names(NT)[i], ", n = ", NT[i], sep = ""), A$Sample.type)}
  A$Sample.type <- as.factor(A$Sample.type)
  
  if(Iso.GDGT == T & Leg.iso == T){
    for(i in 1:length(NT)){A.iso$Sample.type <- gsub(names(NT)[i], paste(names(NT)[i], ", n = ", NT[i], sep = ""), A.iso$Sample.type)}
    A.iso$Sample.type <- as.factor(A.iso$Sample.type)
  }
  
  #### Little sample size br-GDGT ####
  if(is.null(Dot.pop) == F){ # Echantillons choisis comme points
    Little.pop <- A.save[A.save$Sample.type %in% Dot.pop,]
    for(i in 1:length(Dot.pop)){
      New.lab <- unique(A$Sample.type)[grepl(Dot.pop[i], unique(A$Sample.type))]
      Little.pop$Sample.type <- gsub(Dot.pop[i], New.lab, Little.pop$Sample.type)
      if(Overlap.OK == T){A[A$Sample.type %in% New.lab, names(A[,-length(names(A))])] <- 100}}
    
    Not.in.dot <- setdiff(levels(A$Sample.type), Dot.pop)
    New.row.0 <- data.frame(Sample.type = Not.in.dot)
    Little.pop <- rbind.fill(Little.pop, New.row.0)
    Little.pop[is.na(Little.pop)] <- 100
    
    B.little <- melt(Little.pop, id ='Sample.type')
    
    if(is.null(Keep.br) == F){
      Keep.br.little <- levels(B.little$variable)[Keep.br] 
      B.little <- B.little[which(B.little$variable %in% Keep.br.little),]}
    
    Add.points.size.lim <- geom_dotplot(data = B.little, aes(x=variable, y=value, fill=Sample.type),
                                        binaxis='y', stackdir='center', 
                                        position = position_dodge(.85),
                                        show.legend = FALSE,
                                        binwidth = .5)
    B <- melt(A, id ='Sample.type')
  }
  
  else{                   # Pas échatillon comme points
    Add.points.size.lim <- NULL
    B <- melt(A, id ='Sample.type')
  }
  
  #### Color settings ####
  if(is.null(Color.choice) == T){
    if(length(NT) == 1){Color.vec <- c("grey")}
    if(length(NT) == 2){Color.vec <- c("#74D43B", "#1EAEAE")}
    if(length(NT) == 3){Color.vec <- c("#74D43B", "#1EAEAE", "#F3C643")}
    if(length(NT) == 4){Color.vec <- c("#74D43B", "#1EAEAE", "#F3C643", "grey")}
    if(length(NT) == 5){Color.vec <- c("#74D43B", "#1EAEAE", "grey", "#F3C643",  "darkorange")}
    if(length(NT) == 6){Color.vec <- c("#74D43B", "#1EAEAE", "grey", "#F3C643",  "darkorange", "indianred2")}
    if(length(NT) == 7){Color.vec <- c("#74D43B", "#1EAEAE", "grey", "#F3C643",  "darkorange", "indianred2", "darksalmon")}
    if(length(NT) == 8){Color.vec <- c("#74D43B", "#1EAEAE", "grey", "#F3C643",  "darkorange", "indianred2", "darksalmon", "darkred")}
    if(length(NT) > 8){Color.vec <- colorRampPalette(rev(brewer.pal(11, "Spectral")))(NT)}}
  else{Color.vec <- Color.choice}
  
  #### Other graphical settings ####
  if(Global.box == T){My_border <- element_rect(fill = NA)}
  else{My_border <- element_blank()}
  
  if(is.null(Boxplot.title) == F){My_title <- ggtitle(Boxplot.title)}
  else{My_title <- NULL}
  
  if(Leg.box == T){My_legbox <- element_rect(fill = "white", color = "grey")}
  else{My_legbox <- NULL}
  
  #### Change order of boxplot ####
  if(is.null(Reorder.group) == F){
    if(Iso.GDGT == T){
      Old.order.name <- unique(A.iso$Sample.type)[order(unique(A.iso$Sample.type))]
      A.iso$Sample.type <- factor(A.iso$Sample.type, levels = Old.order.name[Reorder.group], ordered = T)}
    
    B$Sample.type <- factor(B$Sample.type, levels = levels(B$Sample.type)[Reorder.group], ordered = T)
    Color.vec <- Color.vec[Reorder.group]}
  
  #### brGDGTs last settings ####
  if(Iso.GDGT == F){A.lab = ylab("Fractional Abundance (%)")}
  if(Iso.GDGT == T){A.lab = ylab("")}
  
  if(is.null(Keep.br) == F){
    Keep.br <- levels(B$variable)[Keep.br] 
    B <- B[which(B$variable %in% Keep.br),]}
  
  if(is.null(Leg.nb.lines) == F){
    Leg.guide <- guides(fill = guide_legend(nrow = Leg.nb.lines))
  }
  else{Leg.guide <- NULL}  
  
  if(Iso.GDGT == T & Leg.iso == T){
    Leg.pos.iso <- Leg.pos
    Leg.pos <- "none"
  }
  else{Leg.pos.iso <- "none"}
  
  #### Add ligne en tirets ####
  Xligne <- c(
    max(grep("III", levels(B$variable)))+0.5,
    max(grep("II", levels(B$variable)))+0.5)
  
  Xtext <- c(Xligne[1]/2-0.25, (Xligne[1]+Xligne[2])/2-0.25, (Xligne[2]+nlevels(B$variable))/2)
  My_lab <- c("Hexa.", "Penta.", "Tetra.")
  
  #### Insert LR ####
  get_inset <- function(B, Zoom1.Ymax, My_annot){
    pinsert <- ggplot(B, aes(x=variable, y=value, fill = Sample.type)) + 
      coord_cartesian(ylim = c(0, Zoom1.Ymax)) +
      geom_boxplot(outlier.shape = NA, position = position_dodge(.85), linewidth = Box.linewidth) +
      Add.points.size.lim + Leg.guide +
      scale_fill_manual(name = "Sample type", values = Color.vec
                        , labels = c("0" = "Foo", "1" = "Bar"))+
      annotate("text", x = 0.8, y = Zoom1.Ymax - 0.1*Zoom1.Ymax, label = My_annot, size = 5, fontface = 2)+
      theme(legend.position = "none", panel.background = element_blank(),
            legend.key = element_blank(), legend.background = element_blank(),
            axis.title = element_blank(),
            axis.text = element_text(size = 12), panel.grid = element_blank(),
            plot.margin = unit(x = c(1, 2, 2, 0),units="mm"), # Bas / Droite / Haut / Gauche
            axis.line = element_blank(),
            legend.text = element_text(size = 13),
            panel.border = element_rect(fill = NA),
            legend.title = element_text(size = 14),
            plot.background = element_blank())+
      A.lab
    return(pinsert)
  }
  
  if(is.null(Zoom1.comp) == F){
    B.insert <- B[which(B$variable %in% levels(B$variable)[Zoom1.comp]),]
    if(missing(Zoom1.Ymax)){Zoom1.Ymax = max(B.insert)}
    if(missing(Insert1.loc)){Insert1.loc = c(3, Xligne[1]-0.25, 10, 40)}
    
    My_insert1 <- get_inset(B.insert, Zoom1.Ymax, "B1")
    My_insert1 <- annotation_custom(ggplotGrob(My_insert1), xmin = Insert1.loc[1], xmax = Insert1.loc[2], ymin = Insert1.loc[3], ymax = Insert1.loc[4])
    Box1 <- geom_rect(xmin = min(Zoom1.comp)-0.5, xmax = Insert1.loc[2], ymin = -1, ymax = Zoom1.Ymax*2, fill = NA,
                      alpha = 1, color = "grey30", linewidth = .3, linetype = 2)
  }
  else{My_insert1 = NULL; Box1 <- NULL}
  
  if(is.null(Zoom2.comp) == F){
    B.insert <- B[which(B$variable %in% levels(B$variable)[Zoom2.comp]),]
    if(missing(Zoom2.Ymax)){Zoom2.Ymax = max(B.insert)}
    if(missing(Insert2.loc)){Insert2.loc = c(11, Xligne[2]-0.25, 10, 40)}
    
    My_insert2 <- get_inset(B.insert, Zoom2.Ymax, "B2")
    My_insert2 <- annotation_custom(ggplotGrob(My_insert2), xmin = Insert2.loc[1], xmax = Insert2.loc[2], ymin = Insert2.loc[3], ymax = Insert2.loc[4])
    Box2 <- geom_rect(xmin = min(Zoom2.comp)-0.5, xmax = Insert2.loc[2], ymin = -1, ymax = Zoom2.Ymax*2, fill = NA,
                      alpha = 1, color = "grey30", linewidth = .3, linetype = 2)
  }
  else{My_insert2 = NULL; Box2 <- NULL}
  
  #### Plot brGDGT ####
  p2 <- ggplot(B, aes(x=variable, y=value, fill = Sample.type)) + 
    coord_cartesian(ylim = c(0, Ymax)) + 
    xlab("brGDGTs") +
    geom_boxplot(outlier.shape = NA, position = position_dodge(.85), linewidth = Box.linewidth) +
    Add.points.size.lim + Leg.guide + My_title +
    geom_vline(xintercept = c(Xligne), lty="dotted")+
    annotate("text", x = Xtext, y = Ymax, label = My_lab, size = Annot.size, hjust = 0.25)+
    scale_fill_manual(name = "Sample type", values = Color.vec
                      , labels = c("0" = "Foo", "1" = "Bar"))+
    theme(legend.position = Leg.pos, panel.background = element_blank(),
          legend.key = element_blank(), legend.background = My_legbox,
          axis.title = element_text(size = 14), panel.border = My_border,
          axis.text = element_text(size = 12), panel.grid = element_blank(),
          plot.margin = unit(x = c(1, 2, 2, 0),units="mm"), # Bas / Droite / Haut / Gauche
          axis.line = element_line(colour = "black"), legend.text = element_text(size = Leg.size),
          legend.title = element_text(size = Leg.size*1.3),
          plot.background = element_blank())+
    A.lab + My_insert1 + My_insert2 + Box1 + Box2
  
  plot_build <- ggplot_build(p2)
  boxplot_stats <- plot_build$data[[1]]
  boxplot_stats <- round(boxplot_stats[c(13,3:5)], digits = 0)
  
  #### Plot isoGDGT ####
  if(Iso.GDGT == T){
    p1 <- ggplot(A.iso, aes(x=variable, y=value, fill=Sample.type)) + 
      coord_cartesian(ylim = c(0, Ymax.iso)) + 
      ylab("Fractional Abundance (%)")+
      xlab("isoGDGTs") +
      geom_boxplot(outlier.shape = NA, position = position_dodge(.85), linewidth = Box.linewidth) +
      scale_x_discrete(labels = c(parse(text = Model.lab), "Crenach\'"))+
      Add.points.size.lim.iso +
      scale_fill_manual(name = "Sample type", values = Color.vec
                        , labels = c("0" = "Foo", "1" = "Bar"))+
      theme(axis.text.x = element_text(angle = 25, hjust = 1),
            legend.position = Leg.pos.iso, panel.grid = element_blank(),
            panel.background = element_blank(), legend.background = element_blank(),
            axis.title = element_text(size = 14),
            axis.text = element_text(size = 12),
            axis.line = element_line(colour = "black"),
            legend.text = element_text(size = 13),
            legend.title = element_text(size = 14),
            plot.margin = unit(x = c(0, 0, 0, 0),units="mm"), 
            plot.background = element_blank()
      )
    
    # p3 <- plot_grid(p1, p2, ncol = 2, axis="tb", align = "h", rel_widths = c(3/10, 7/10), labels = c("A","B"), label_size = 18, label_x = -0.01) # 2 graphs
    p3 <- p1 + p2 + plot_layout(ncol = 2, widths = c(3/10, 7/10)) + plot_annotation(tag_levels = 'A') &
      theme(plot.tag = element_text(size = 18, face = "bold", hjust = 0.5, vjust = -3.5), strip.background = element_blank(),
            plot.background = element_blank(), plot.margin = unit(x = c(0, 0, 0, 0),units="mm"))
    
  }
  
  #### Save html ####
  if(Show.Plotly == T){
    library(plotly)
    library(htmlwidgets)
    Save.plot.html <- gsub("pdf", "html", Save.path)
    Keep.name <- gsub(".*\\/", "", Save.plot.html)
    Path.root <- paste(gsub(Keep.name, "", Save.plot.html), "HTML_files/", sep = "")
    if(file.exists(Path.root) == F){dir.create(Path.root)}
    Save.plot.html <- paste(Path.root, Keep.name, sep = "")
    p1_ly <- ggplotly(p2)
    p1_ly <- p1_ly %>% layout(boxmode = "group", boxpoints = F, legend = list(font = list(size = 12)))
    options(warn = - 1) 
    saveWidget(p1_ly, file = Save.plot.html)
  }
  #### Save plots ####
  if(is.null(Save.path) == F){
    if(is.null(W) == F & is.null(H) == F){
      ggsave(Save.path, width = W*0.026458333, height = H*0.026458333, units = "cm")}
    else{ggsave(Save.path)}}
  else(
    if(Iso.GDGT == T){return(p3)}
    else{return(p2)})
  # if(is.null(p2) == F){View(ggplot_build(p1)$data[[1]])}
  
  if(Return.plot == T){
    if(Iso.GDGT == T){return(p3)}
    else{return(p2)}
  }
}

Diag.ternaire.methylation <- function(MGDGT, Mcol, My_colors, Show.Dearing, Show.Naafs.peat, Show.Plotly, Return.plot = F,
                                      Full.labels = T, Legend.pos = "right", Annot = NULL, Show.arrows = T,
                                      Show.lake = F, Add.facet = F, Show.soil = F, Show.peat = F, Remove.ACA = F,
                                      Alpha.dot, Size.dot, Export.to.chart.studio = F, Save.path, W, H){
  #### Initialization values ####
  library(ggplot2)
  library(ggtern)            # permet de faire des graphiques ternaires
  if(missing(MGDGT)){warning("Import a matrix of GDGT")}
  if(missing(Show.Naafs.peat)){Show.Naafs.peat = F}
  if(missing(Show.Plotly)){Show.Plotly = F}
  if(missing(Show.Dearing)){Show.Dearing = F}
  if(missing(Save.path)){Save.path = NULL}
  if(missing(Alpha.dot)){Alpha.dot = 0.6}
  if(missing(Size.dot)){Size.dot = 2.5}
  if(missing(W)){W = NULL}
  if(missing(H)){H = NULL}
  if(missing(Mcol)){Mcol = NULL}
  if(missing(My_colors)){My_colors = NULL}

  #### Import WDB data ####
  Msurf.mean.WDB <- readRDS("Import/Training/MbrGDGT_WDB.Rds")
  Meco.WDB <- readRDS("Import/Training/Meco_WDB.Rds")

  #### Extract br-GDGT ####
  Cal.Mtern <- function(M){
    M <- M[,grep("^f.I", colnames(M))]
    III <- grep("III", names(M))
    III.II <- grep("II", names(M))
    I.II.III <- grep("I", names(M))
    II <- setdiff(III.II, III)
    I <- setdiff(I.II.III, III.II)
    Mtern <- M[0]
    Mtern[["P.hexa"]] <- rowSums(M[,III])
    Mtern[["P.penta"]] <- rowSums(M[,II])
    Mtern[["P.tetra"]] <- rowSums(M[,I])
    Mtern <- data.frame(apply(Mtern, 2, function(x) x/rowSums(Mtern)))
    return(Mtern)}
  Mtern <- Cal.Mtern(MGDGT)
  if(is.null(Mcol)==T){Mtern[["Group"]] <- "Surface"}
  else{Mtern[["Group"]] <- Mcol}

  #### Add Other databases ####
  if(Show.Dearing == T){
    Msurf.Dearing <- data.frame(read.csv(file="Import/World_DB/GDGT/Dearing_Crampton-Flood-etal_2019.csv",sep="\t",dec=".",header=T, row.names = 1))        # GDGT indexe Ayrag
    Mtern.dear <- Cal.Mtern(Msurf.Dearing)
    Mtern.dear$Group <- paste(Msurf.Dearing[["Samp.type"]], "Dearing")
    Mtern <- rbind(Mtern.dear, Mtern)
  }
  if(Show.Naafs.peat == T){
    Msurf.Naaf.peat <- data.frame(read.csv(file="Import/World_DB/GDGT/GDGT_Naaf_peat_2017.csv",sep="\t",dec=".",header=T, row.names=1))        # GDGT indexe Ayrag
    Mtern.naaf <- Cal.Mtern(Msurf.Naaf.peat)
    Mtern.naaf$Group <- "Peat Naafs"
    Mtern <- rbind(Mtern.naaf, Mtern)
  }
  if(Remove.ACA == T){
    Msurf.mean.WDB <- Msurf.mean.WDB[Meco.WDB$DB != "ACA",]
    Meco.WDB <- Meco.WDB[Meco.WDB$DB != "ACA",]
  }
  if(Show.lake == T){
    Msurf.Raberg <- Msurf.mean.WDB[Meco.WDB$Sample.type == "Lacustrine",]
    Mtern.Raberg <- Cal.Mtern(Msurf.Raberg)
    Mtern.Raberg$Group <- "Lacustrine WDB"
    Mtern <- rbind(Mtern.Raberg, Mtern)

  }
  if(Show.soil == T){
    Msurf.Raberg <- Msurf.mean.WDB[Meco.WDB$Sample.type == "Soil",]
    Mtern.Raberg <- Cal.Mtern(Msurf.Raberg)
    Mtern.Raberg$Group <- "Soil WDB"
    Mtern <- rbind(Mtern.Raberg, Mtern)
  }

  #### Color & shape settings ####
  Lab.list <- c(
    Lake_Ayrag = "grey30", Surface = "grey30", WZDB = "grey30",
    Lakes = "royalblue", Lacustrine = "royalblue",
    'Lacustrine WDB' = "royalblue", "Soil WDB" = "darkorange",
    "Fibrous soil" = "#bf812e", "K-warm/arid" = "#8c510a", "K-cold/wet" = "#80cec1", "K-lacustrine" = "#0073C2",
    "Topcore silk" = "#8c510a",
    "1_Fresh" = "royalblue","2_Brakish" = "#e4af08","3_Salted" = "#8c510a",
    "1_Acid" = "royalblue","2_Neutral" = "#e4af08","3_Alkalin" = "#8c510a",
    "1_Hyper-arid" = "#8c510a","2_Arid" = "#bf812e", "3_Semi-arid" = "#dfc27e","4_Dry sub-humid" = "#f5e9bf", "5_Humid" = "#80cec1",
    "1. Hyper-arid" = "#8c510a","2. Arid" = "#bf812e", "3. Semi-arid" = "#dfc27e","4. Dry sub-humid" = "#f5e9bf", "5. Humid" = "#80cec1",
    "Dry-Cold" = "darkblue","Wet-Cold" = "darkblue","Dry-Warm" = "darkorange","Wet-Warm" = "darkorange", "Lake Raberg" = "royalblue",
    Soil = "darkorange", "Soil Dearing" = "darkorange", "Soil ACA" = "#c3510aff",
    Peat = "#74D43B", "peat Dearing" = "#74D43B", "Peat Naafs" = "#74D43B",
    Moss = "#74D43B", "Moss ACA" = "#1d5140ff",
    Azerbaijan = "darkorange", Mongolia = "royalblue", Chine = "red", Uzbekistan = "green",
    Tajikistan = "pink", Fazilman = "darkblue", Tuya = "grey30", Ogshagil = "#1d5140ff", Russia = "#74D43B",
    L1 = "#515a14ff", L2 =  "#abc837f6", L3 =  "#70ae77f6", L4 = "#e8a147ff", L5 = "#313695", L6 =  "#577AB7", L7 = "#b65928ff", L8 = "#7e332fff"
  )

  if(is.null(My_colors) == F){
    names(My_colors) <- unique(Mcol)
    Lab.list <- c(Lab.list, My_colors)}

  Shape.list <- c(
    Lakes = 16, Surface = 16, Lacustrine = 16,
    "Topcore silk" = 16, Soil = 4, "Soil ACA" = 16, "Fibrous soil" = 4, "Soil Dearing" = 4,
    "Dry-Cold" = 4, "Wet-Cold" = 19, # "Dry-Warm" = 4, "Wet-Warm" = 19,
    "Lacustrine WDB" = 4, "Soil WDB" = 4, WZDB = 2,
    "1_Hyper-arid" = 16,"2_Arid" = 16, "3_Semi-arid" = 16,"4_Dry sub-humid" = 16, "5_Humid" = 16,
    "1. Hyper-arid" = 16,"2. Arid" = 16, "3. Semi-arid" = 16,"4. Dry sub-humid" = 16, "5. Humid" = 16,
    "Lake Raberg" = 3, Peat = 3, "peat Dearing" = 3, "Peat Naafs" = 3, Moss = 3, "Moss ACA" = 16
  )

  if(length(setdiff(unique(Mtern$Group), names(Shape.list))) > 0){
    Shape.list2 <- rep(19, length(unique(Mcol)))
    names(Shape.list2) <- unique(Mcol)
    Shape.list <- c(Shape.list, Shape.list2)
    Shape.list <- Shape.list[names(Shape.list)%in% unique(Mtern$Group)]
  }

  #### Labels settings ####
  if(Full.labels == T){My_labs <- labs(x = "Hexamethylated \n brGDGTs", y = "Pentamethylated \n brGDGTs", z = "Tetramethylated \n brGDGTs")}
  if(Full.labels == F){My_labs <- labs(x = "Hexa.", y = "Penta.", z = "Tetra.")}

  if(is.null(Annot) == F){My_annot <- labs(title = Annot)}
  else{My_annot <- NULL}

  if(Show.arrows == T){Arrows <- NULL}
  else{Arrows <- ggtern::theme_noarrows()}

  My_title <- element_text(hjust = 0, vjust = -6.5, face = "bold", size = 18)

  #### Plot ####
  if(Add.facet == F){
    p1 <- ggtern::ggtern(data = Mtern, aes(P.hexa, P.penta, P.tetra, color = Group, shape = Group)) +
      ggtern::theme_rgbw() + My_annot +
      Arrows +
      theme(legend.key = element_blank(), legend.position = Legend.pos, plot.margin=unit(c(0,0,0,0),"cm"),
            plot.background = element_blank(), panel.background = element_blank(),
            legend.background = element_blank(), strip.background = element_blank(),
            plot.title = My_title,
            legend.key.height = unit(1, "line"))+
      geom_point(alpha = Alpha.dot, size = Size.dot, na.rm = T) +
      scale_color_manual(name = "Samples type", values = Lab.list) +
      scale_shape_manual(name = "Samples type", values = Shape.list)+
      My_labs}

  #### Add facets ####
  if(Add.facet == T){
    p1 <- ggtern::ggtern(data = Mtern, aes(P.hexa, P.penta, P.tetra, color = Group, shape = Group)) +
      My_annot +
      ggtern::theme_nogrid()+
      ggtern::theme_ticksinside()+
      ggtern::theme_nolabels()+ ggtern::theme_notitles()+
      theme(legend.key = element_blank(), legend.position = "none", plot.margin=unit(c(0,0,0,0),"cm"),
            plot.background = element_blank(), panel.background = element_blank(),
            legend.background = element_blank(), strip.background = element_blank(),
            plot.title = My_title, panel.grid = element_blank(),
            legend.key.height = unit(1, "line"))+
      geom_point(size = 0.1) +
      facet_wrap(vars(Group))+
      ggtern::stat_density_tern(geom = 'polygon', bins = 10, aes(fill = Group, alpha = ..level..), color = NA, linewidth = .005) +
      scale_color_manual(name = "Samples type", values = Lab.list) +
      scale_fill_manual(name = "Samples type", values = Lab.list) +
      My_labs
  }
  #### Save html ####
  if(Show.Plotly == T){
    library(plotly)
    library(htmlwidgets)

    m <- list(
      l = 150,
      r = 0,
      b = 0,
      t = 0,
      pad = 0
    )

    Save.plot.html <- gsub("pdf", "html", Save.path)
    Keep.name <- gsub(".*\\/", "", Save.plot.html)
    Path.root <- paste(gsub(Keep.name, "", Save.plot.html), "HTML_files/", sep = "")
    if(file.exists(Path.root) == F){dir.create(Path.root)}
    Save.plot.html <- paste(Path.root, Keep.name, sep = "")
    Mtern$label <- row.names(Mtern)

    axis <- function(title) {list(title = title,
                                  titlefont = list(size = 18),tickfont = list(size = 13),
                                  tickcolor = 'rgba(0,0,0,0)',ticklen = 5)}

    fig <- Mtern %>% plot_ly()
    fig <- fig %>% add_trace(type = 'scatterternary',mode = 'markers',
                             a = ~P.penta,b = ~P.hexa,c = ~P.tetra,
                             text = ~label,
                             color = ~Group, colors = Lab.list,
                             symbol = ~Group, symbols = Shape.list,
                             marker = list(size = 10, line = list('width' = 2)))

    fig <- fig %>% layout(
      # title = "Simple Ternary Plot with Markers",
      margin = m,
      ternary = list(sum = 100,
                     aaxis = axis('Pentamethylated brGDGTs'),
                     baxis = axis('Hexamethylated brGDGTs'),
                     caxis = axis('Tetramethylated brGDGTs')))

    saveWidget(fig, file = Save.plot.html)

    if(Export.to.chart.studio == T){
      Sys.setenv("plotly_username"="Lucas_Dugerdil")
      Sys.setenv("plotly_api_key"="fBxVNd3kb1T6msCwEJCL")
      Name.chartplot <- gsub("\\.html", "", Keep.name)
      api_create(fig, filename = Name.chartplot)
    }
  }
  #### Save plots ####
  if(is.null(Save.path) == F){
    if(is.null(W) == F & is.null(H) == F){
      ggsave(Save.path, width = W*0.026458333, height = H*0.026458333, units = "cm")}
    else{ggsave(Save.path)}}

  detach("package:ggtern", unload=TRUE)
  
  if(Return.plot == T){return(p1)}
  else(return(Mtern))
}


# RDA sur la veget avec DB abiotique
# Soit Climat, soil...
# Param necessaire : MP, MClim
# Param option : transp_OK, Site.name, Csv.sep, Scale.PCA
RDA.pollen.surf <- function(MP, MClim, Choose.clim, Cluster.path = NULL, Alpha.dot = NULL, Type.samples, Shape.groups = NULL, Sort.eco, Color.choice = NULL, Result.vector = F,
                            Cluster.groups, Display.legends, transp_OK, Manu.lim, GDGT, Annot, Vector.show = NULL, Nb.contrib = NULL, Leg.size = 1, Sort.shape = NULL, Result.clim = F,
                            Remove.7Me, Leg.loc, Helinger.trans = F, Display.plot = T, Leg.loc2, Simple.title, Show.text, Symbol.loc, Symbol.path, VIF = F, VIF.seuil = NULL, Dot.size = 1.5,
                            Display.VIF = T, ggplot.display = F, Marg.density.plot = F, return.VIF = F, Vector.env.scale = 1, Vector.scale = 1, Site.name, Csv.sep, Scale.taxa, Scale.sites, Save.path, Save.plot = NULL, Title.inside = F, return.pick = F){
  #### Settings ####
  library(vegan)
  # library(missMDA)
  if(missing(Type.samples)){Type.samples = NULL}
  if(missing(Csv.sep)){Csv.sep = "\t"}
  if(missing(Annot)){Annot = NULL}
  if(missing(Sort.eco)){Sort.eco = NULL}
  if(missing(Show.text)){Show.text = F}
  if(missing(Remove.7Me)){Remove.7Me = F}
  if(missing(Simple.title)){Simple.title = F}
  if(missing(Cluster.groups)){
    if(is.null(Cluster.path) == F){
      Cluster.groups = c("Ecosystems", "Biomes", "Vegetation", "Aridity2", "Aridity")
      Cluster.groups <- names(Cluster.path)[names(Cluster.path)  %in% Cluster.groups][1]
    }
  }
  if(missing(Display.legends)){Display.legends = T}
  if(missing(Choose.clim)){print("Select the climat variable to perform the RDA.")}
  if(missing(Scale.taxa)){Scale.taxa = 1}
  if(missing(Scale.sites)){Scale.sites = 2}
  if(missing(Save.path)){Save.path = NULL}
  if(missing(Manu.lim)){Manu.lim = NULL}
  if(missing(GDGT)){GDGT = NULL}
  if(missing(Site.name)){Site.name = ""}
  if(missing(Leg.loc)){Leg.loc = "bottomleft"}
  if(missing(Leg.loc2)){Leg.loc2 = "bottomright"}
  if(missing(Symbol.path)){Symbol.path = NULL}
  if(missing(Symbol.loc)){Symbol.loc = NULL}
  
  #### Save plots ####
  if(is.null(Save.plot) == F){
    Path.to.create <- gsub("(.*/).*\\.pdf.*","\\1", Save.plot)
    dir.create(file.path(Path.to.create), showWarnings = FALSE)
    if(is.null(W) == F & is.null(H) == F){
      pdf(file = Save.plot, width = W*0.01041666666667, height = H*0.01041666666667)}
    else{pdf(file = Save.plot)}}
  
  #### Pour les GDGTs ####
  if(is.null(GDGT) == F){
    if(Remove.7Me == T){MP <- MP[setdiff(row.names(MP), row.names(MP)[grepl("7Me", row.names(MP))]),]}
    row.names(MP) <- gsub("f.", "", row.names(MP))
    row.names(MP) <- gsub("_5Me", "", row.names(MP))
    row.names(MP) <- gsub("_6Me", "\\'", row.names(MP))
    row.names(MP) <- gsub("_7Me", "\\''", row.names(MP))
  }
  
  #### Pour le pollen ####
  if(is.null(GDGT) == T){
    row.names(MP) <- gsub("aceae",".",row.names(MP))
    row.names(MP) <- gsub(".undiff","",row.names(MP))
    row.names(MP) <- gsub(".indet","",row.names(MP))
    row.names(MP) <- gsub(".*Po.$","Poa.",row.names(MP))
    row.names(MP) <- gsub(".spp.","",row.names(MP))
    row.names(MP) <- gsub(".type","",row.names(MP))
  }
  
  MP <- data.frame(t(MP), check.names = F)
  
  #### Import Climat + verif match DB ####
  if(is.character(MClim) == T){Clim <- data.frame(read.csv(MClim, sep = Csv.sep, dec = ".", header=T, row.names=1), stringsAsFactors = T)}
  if(is.data.frame(MClim) == T){Clim <- MClim}
  C <- colnames(Clim)
  Inter <- intersect(row.names(Clim), row.names(MP))
  Clim <- Clim[which(row.names(Clim) %in% Inter),]
  MP <- MP[which(row.names(MP) %in% Inter),]
  Clim <- data.frame(Clim[][row.names(MP),])
  row.names(Clim) <- row.names(MP)
  colnames(Clim) <- C
  Clim.rda <- subset(Clim, select = Choose.clim)
  
  #### Complete missing informations ####
  # if(nlevels(as.factor(is.na(Clim.rda))) >= 2){
  #   library(missMDA)
  #   nb <- estim_ncpPCA(Clim.rda, ncp.max = 5) ## Time consuming, nb = 2
  #   Clim.rda.comp <- imputePCA(Clim.rda, ncp = nb[[1]])
  #   Clim.rda <- Clim.rda.comp$completeObs
  # }
  
  #### Transforming the data ####
  if(transp_OK == F){MP.pca <- rda(MP, Clim.rda, scale = T)
  MTitle = paste("RDA", Type.samples, "/ climate -", Site.name)}
  else{
    if(Helinger.trans == T){
      MP <- vegan::decostand(MP, method = "hellinger")
      MP.pca <- rda(MP, scale(Clim.rda), scale = F)
    }
    else{MP.pca <- rda(log1p(MP), scale(Clim.rda), scale = T)}
    
    MTitle = paste("RDA", Type.samples, "/ climate - log-trans, ", Site.name)}
  if(Simple.title == T){MTitle = paste("RDA", Type.samples, "/ climate")}
  if(is.null(Annot) == F){MTitle <- paste(Annot, MTitle)}
  if(Title.inside == T){MTitle <- NULL}
  
  #### VIF test ####
  if(VIF == T){
    my_VIF <- vif.cca(MP.pca)
    
    if(is.null(VIF.seuil) == F){
      if(Display.VIF == T){print("**** VIF results ****"); print(round(my_VIF, digits = 1))}
      my_VIF <- my_VIF[my_VIF < VIF.seuil]
      Clim.rda <- Clim.rda[which(names(Clim.rda) %in% names(my_VIF))]
      
      if(transp_OK == F){MP.pca <- rda(MP, Clim.rda, scale = T)}
      else{
        if(Helinger.trans == T){
          MP <- vegan::decostand(MP, method = "hellinger")
          MP.pca <- rda(MP, scale(Clim.rda), scale = F)
        }
        else{MP.pca <- rda(log1p(MP), scale(Clim.rda), scale = T)}
        # MP.pca <- rda(log1p(MP), Clim.rda, scale = T)
      }
    }
    if(Display.VIF == T){print("**** VIF results **** If you want to remove useless vectors, use VIF.seuil."); print(round(my_VIF, digits = 1))}
    
  }
  #### Calcul parameters ####
  C_PC.rda <- coef(MP.pca)    # récupère les coefficients canoniques (equivalents R en multivar (donc il faut faire au 2 pour avoir R2)) pour chaque variable clim
  PClim.sc <- scores(MP.pca, choices=1:2, scaling = Scale.taxa, display = "sp")
  PClim.sc.site <- scores(MP.pca, choices=1:2, scaling = Scale.sites, display = "sites")
  
  clim.score <- MP.pca[["CCA"]][["biplot"]]
  R2 <- RsquareAdj(MP.pca)$r.squared
  Tot.eig <- sum(MP.pca[["CCA"]]$eig) + sum(MP.pca[["CA"]]$eig)
  RDA1 <- MP.pca[["CCA"]]$eig[1]/Tot.eig*100
  RDA2 <- MP.pca[["CCA"]]$eig[2]/Tot.eig*100
  
  #### Graphical parameters ####
  xmin <- min(c(min(PClim.sc.site[,1]),2*min(PClim.sc[,1])))                             # lim axe PC1
  xmax <- max(c(max(PClim.sc.site[,1]),2*max(PClim.sc[,1])))
  ymin <- min(c(min(PClim.sc.site[,2]),2*min(PClim.sc[,2])))                             # lim axe PC2
  ymax <- max(c(max(PClim.sc.site[,2]),2*max(PClim.sc[,2])))
  
  xmin <- xmin + .08*xmin
  xmax <- xmax + .08*xmin
  ymin <- ymin + .08*ymin
  ymax <- ymax + .08*ymax
  
  if(is.null(Manu.lim) == F){
    xmin <- Manu.lim[1]
    xmax <- Manu.lim[2]
    ymin <- Manu.lim[3]
    ymax <- Manu.lim[4]
  }
  #par(mar=c(5,4,4,2)+0.1)#,xpd=TRUE)
  par(mgp = c(1.7,0.6,0), mar=c(3,3,2,0.3)+0.1)
  
  #### Contribution ####
  if(is.null(Vector.show) == F & is.null(Nb.contrib) == T){
    PClim.sc <- PClim.sc[row.names(PClim.sc) %in% Vector.show,]
    MP.pca$CCA$v <- MP.pca$CCA$v[row.names(MP.pca$CCA$v) %in% Vector.show,]
    MP.pca$CA$v <- MP.pca$CA$v[row.names(MP.pca$CA$v) %in% Vector.show,]
    MP <- MP[colnames(MP) %in% Vector.show]}
  
  if(is.null(Nb.contrib) == F){
    Contribution <- sqrt(abs(scores(MP.pca, display = "sp", scale = 0)[,1])^2 + abs(scores(MP.pca, display = "sp", scale = 0)[,2])^2)
    Contribution <- Contribution[order(Contribution, decreasing = T)]
    Contribution <- names(Contribution[1:Nb.contrib])
    if(is.null(Vector.show) == F){
      Contribution <- unique(c(Contribution, Vector.show))
      print(paste("**** We will merge the", Nb.contrib, "most contributing taxa + the", length(Vector.show), "manually selected taxa.****" ))}
    
    PClim.sc <- PClim.sc[row.names(PClim.sc) %in% Contribution,]
    MP.pca$CCA$v <- MP.pca$CCA$v[row.names(MP.pca$CCA$v) %in% Contribution,]
    MP.pca$CA$v <- MP.pca$CA$v[row.names(MP.pca$CA$v) %in% Contribution,]
    MP <- MP[colnames(MP) %in% Contribution]}
  
  #### Clustering colors ####
  if(is.null(Cluster.path) == F){
    if(is.character(Cluster.path) == T){Cluster <- data.frame(read.csv(Cluster.path, sep = Csv.sep ,dec=".",header=T,row.names=1), stringsAsFactors = T)}
    if(is.data.frame(Cluster.path) == T){Cluster <- Cluster.path}
    Inter <- intersect(row.names(Cluster), row.names(PClim.sc.site))
    Cluster <- Cluster[which(row.names(Cluster) %in% Inter),]
    M.eco <- subset(Cluster, select = Cluster.groups)
    M.eco <- data.frame(M.eco[][row.names(PClim.sc.site),])
    row.names(M.eco) <- row.names(PClim.sc.site)
    Fact.eco <- as.factor(M.eco[[1]])
    if(is.null(Sort.eco) == F){Fact.eco <- ordered(Fact.eco, levels = Sort.eco)}
    num_colors <- nlevels(Fact.eco)
    if(is.null(Color.choice) == F){diamond_color_colors <- Color.choice}
    else{
      colfunc    <- colorRampPalette(c("firebrick3", "darkorange", "goldenrod1", "#38A700", "darkgreen", "dodgerblue3", "grey10"))
      diamond_color_colors <- colfunc(num_colors)
    }
    Col.eco <- diamond_color_colors[Fact.eco]
  }
  else{Col.eco = 1}
  
  if(is.null(Alpha.dot) == F){
    Col.eco <- adjustcolor(Col.eco, alpha.f = Alpha.dot)
  }
  
  #### Color scales ggplots ####
  if(is.null(Cluster.path) == F){
    
    values.bi = c("Deserts & Xeric Shrublands" = "#C88282",
                  "Temperate Grasslands, Savannas & Shrublands" = "#ECED8A",
                  "Montane Grasslands & Shrublands" = "#D0C3A7",
                  "Temperate Conifer Forests" = "#6B9A88",
                  "Temperate Broadleaf & Mixed Forests" = "#3E8A70",
                  "N/A" = "#FFEAAF",
                  "Tundra" = "#A9D1C2",
                  "Boreal Forests/Taiga" = "#8FB8E6",
                  "Tropical & Subtropical Coniferous Forests" = "#99CA81",
                  "Mangroves" = "#FE01C4",
                  "Flooded Grasslands & Savannas" = "#BEE7FF",
                  "Tropical & Subtropical Moist Broadleaf Forests" = "#38A700",
                  "Plant_height" = "royalblue",
                  "Leaf_thickness" = "darkorange",
                  "Photosynthesis_pathway" = "purple",
                  "TUSDB sites" = "#323232",
                  "Woodyness" = "#323232",
                  "Leaf_size" = "darkred",
                  "Variable" = "#aa373aff",
                  "Algal" = "#4666E9",
                  "NAP" = "#b5ab32ff",
                  "Herb" = "#b5ab32ff",
                  "Shrub" = "#aa373aff",
                  "Other" = "grey90",
                  "Unknown" = "grey90",
                  "AP" = "#0f6b31ff",
                  "Tree" = "#0f6b31ff",
                  "1_Hyper-arid" = "#8c510a", "2_Arid" = "#bf812e", "3_Semi-arid" = "#dfc27e", "4_Dry sub-humid" = "#f5e9bf", "5_Humid" = "#80cec1",
                  "1. Hyper-arid" = "#8c510a", "2. Arid" = "#bf812e", "3. Semi-arid" = "#dfc27e", "4. Dry sub-humid" = "#f5e9bf", "5. Humid" = "#80cec1",
                  Mongolia = "#3e96bdff", Chine = "#f02a26", Uzbekistan = "#6fb440", Armenia = "#54a697", "China, Tibet" = "#8c510a", "Northern Iran" = "#176E5B",
                  Tajikistan = "#e4af08", Russia = "#0035a9", Azerbaijan = "#094227", China = "#bb0202", "ACA lakes" =  "purple",
                  "Ch'ol cold desert-steppes" = "#7916C4", 
                  "Tugai riparian forest" = "#BB0268", 
                  "Ch'ol warm deserts" = "#bb0202", 
                  "Adyr desert-steppes" = "#ff5400", 
                  "Adyr steppes" = "#e6c607", 
                  "Tau riparian forest" = "#2C9740", 
                  "Tau thermophilous woodlands" = "#85682D", 
                  "Tau juniper steppe-forest" = "#176E5B",
                  "Tau steppes" = "#bab133",
                  "Alau cryophilous steppe-forest" = "#54a697",
                  "Alau meadows" = "#197CDA",
                  "Tugai riparian forests" = "#7916C4", 
                  "Ch'ol deserts" = "#bb0202", 
                  "Adyr pseudosteppes" = "#ff5400", 
                  "Adyr steppes" = "#ECED8A", 
                  "Tau xeric shrublands" = "#85682D", 
                  "Tau open woodlands" = "#1e8736",
                  "Tau steppes" = "#b0a62e",
                  "Alau cryophilous open woodlands" = "#54a697",
                  "Alau mesic grasslands" = "#197CDA" 
    )
    
    
    
    # values.bi <- values.bi[which(names(values.bi) %in% unique(Keep.xdata[[Cluster.core]]))]
    Scale.fill <- scale_fill_manual(values = values.bi, drop = T)
    Scale.color <- scale_color_manual(values = values.bi, drop = T, guide = "none")
  }
  else{Scale.fill <- NULL; Scale.color <- NULL}
  #### Plot ggplot ####
  if(Display.plot == T & ggplot.display == T){
    #### Settings ####
    Text.size = 4
    Mgg <- setNames(data.frame(PClim.sc.site), c("RDA1", "RDA2"))
    Mgg$My_color <- M.eco[[1]]
    LAB <- paste("list(italic(R^2) == ", round(R2, 2), ", n == ", nrow(MP), ")")
    annotations <- data.frame(
      annotateText = LAB,
      xpos = c(Inf),
      ypos =  c(-Inf),
      hjustvar = c(1.1) ,
      vjustvar = c(-0.5))
    My_annot <- geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), parse = T, size = Text.size)
    
    if(Marg.density.plot == F){Title.global <- ggtitle(MTitle)}
    else{Title.global <- NULL}
    print(Vector.scale)
    #### Main plot ####
    p <- ggplot(data = Mgg, aes(x = RDA1, y = RDA2)) +
      geom_point(aes(color = My_color), size = Dot.size, alpha = Alpha.dot) +
      geom_segment(data = PClim.sc, aes(x = 0, y = 0, xend = RDA1*Vector.scale, yend = RDA2*Vector.scale), arrow = arrow(length = unit(0.2, "cm")), color = "royalblue") +
      geom_text(data = PClim.sc, aes(x = RDA1*Vector.scale, y = RDA2*Vector.scale, label = rownames(PClim.sc)), color = "royalblue", vjust = -0.5, size = Text.size) +
      geom_segment(data = clim.score, aes(x = 0, y = 0, xend = RDA1*Vector.env.scale, yend = RDA2*Vector.env.scale), arrow = arrow(length = unit(0.4, "cm")), color = "darkred", linewidth = 1) +
      geom_text(data = clim.score, aes(x = RDA1*Vector.env.scale, y = RDA2*Vector.env.scale, label = rownames(clim.score)), color = "darkred", vjust = -0.5, size = Text.size*1.3) +
      My_annot +
      geom_vline(xintercept = 0, lty = "dashed")+
      geom_hline(yintercept = 0, lty = "dashed")+
      xlab(bquote(RDA[1] ~ "(" ~ .(format(RDA1, digits = 2)) ~ "%)")) +
      ylab(bquote(RDA[2] ~ "(" ~ .(format(RDA2, digits = 2)) ~ "% )")) +
      Title.global+ Scale.fill + Scale.color +
      theme_classic()+
      theme(axis.line = element_blank(), legend.background = element_blank(),
            plot.background = element_blank(), panel.background = element_blank(),
            strip.placement = "outside", legend.position = Leg.loc,
            panel.border = element_rect(NA, "black", linewidth = 1))
    
    #### Add margin density ####
    if(Marg.density.plot == T){
      x_limits <- ggplot_build(p)$layout$panel_scales_x[[1]]$range$range
      y_limits <- ggplot_build(p)$layout$panel_scales_y[[1]]$range$range
      List.of.NA <- which(Mgg$RDA1 < 1e-12 & Mgg$RDA1 > -1e-12 & Mgg$RDA2 < 1e-12 & Mgg$RDA2 > -1e-12)
      if(length(List.of.NA) > 0){
        print("Remove NA from density.")
        Mgg <- Mgg[-List.of.NA,]}
      
      #### Density plots up ####
      plot_top <- ggplot(Mgg, aes(x = RDA1, fill = My_color)) + 
        geom_density(alpha = 0.6, size = 0.1) + Scale.fill + 
        scale_x_continuous(limits = x_limits, expand = c(0,0))+
        ggtitle(MTitle)+ 
        #### Theme ####
      theme(
        axis.line = element_blank(),
        axis.text.x = element_blank(), axis.text.y = element_text(hjust = 1, size = 6),
        axis.ticks.x.bottom = element_blank(),
        axis.title = element_blank(),
        axis.line.y = element_line(colour = "grey"),
        axis.ticks.y = element_line(colour = "grey"),
        #panel.border = element_rect(fill = NA, colour = "grey"),
        legend.title = element_text(),
        legend.key = element_blank(),
        legend.justification = c("center"),               # left, top, right, bottom
        legend.text = element_text(size = 8),
        panel.background = element_blank(),
        panel.spacing = unit(0.7, "lines"),
        legend.position = "none",
        strip.text.x = element_text(size = 12, angle = 0, face = "bold"),
        strip.placement = "outside",
        strip.background = element_rect(color = "white", fill = "white"),
        plot.margin=unit(c(0,0,0,0),"cm")
      )
      
      #### Density plots right ####
      plot_right <- ggplot(Mgg, aes(x = RDA2, fill = My_color)) + 
        geom_density(alpha = 0.6, size = 0.1) + Scale.fill +
        scale_x_continuous(limits = y_limits, expand = c(0,0))+
        coord_flip() + 
        #### Theme ####
      theme(
        axis.line.y = element_blank(),
        axis.text.y = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        axis.line.x = element_line(colour = "grey"),
        axis.ticks.x = element_line(colour = "grey"),
        legend.title = element_text(),
        legend.key = element_blank(),
        legend.justification = c("center"),               # left, top, right, bottom
        legend.position = "none",
        legend.text = element_text(size = 8),
        panel.background = element_blank(),
        panel.spacing = unit(0.7, "lines"),
        strip.text.x = element_text(size = 12, angle = 0, face = "bold"),
        strip.placement = "outside",
        strip.background = element_rect(color = "white", fill = "white"),
        plot.margin=unit(c(0,0,0,0),"cm")
      )
      
      
      #### Export ####
      layout <- "AAAAAAA#
                 CCCCCCCB
                 CCCCCCCB
                 CCCCCCCB
                 CCCCCCCB
                 CCCCCCCB
                 "
      p <- plot_top + plot_right + p + plot_layout(design = layout) & theme(plot.margin = unit(c(0,0,0,0),"cm"))
    }
    Display.plot <- F
  }
  
  #### Plot (R base) ####
  if(Display.plot == T | is.null(Save.plot) == F){
    #### Vector plot ####
    plot(MP.pca,
         type = "n",
         scaling = Scale.sites, 
         main = MTitle,
         xlim = c(xmin,xmax),
         ylim = c(ymin,xmax),
         xlab = bquote(RDA[1] ~ "(" ~ .(format(RDA1, digits = 2)) ~ "%)"),
         ylab = bquote(RDA[2] ~ "(" ~ .(format(RDA2, digits = 2)) ~ "% )")
    )
    
    
    
    #### Shapes ####
    if(is.null(Shape.groups) == F){
      M.shape <- subset(Cluster, select = Shape.groups)
      M.shape <- data.frame(M.shape[][row.names(PClim.sc.site),])
      row.names(M.shape) <- row.names(PClim.sc.site)
      Fact.shape <- as.factor(M.shape[[1]])
      if(is.null(Sort.shape) == F){Fact.shape <- ordered(Fact.shape, levels = Sort.shape)}
      num_shape <- nlevels(Fact.shape)
      Shape.list <- c(19,18,20,12,11,10,9,18, 14)
      Shape.list <- Shape.list[1:num_shape]
      names(Shape.list) <- levels(Fact.shape)
      Shape.eco <- as.vector(Shape.list[Fact.shape])
    }
    else{Shape.eco = 19}
    
    #### Add plots ####
    if(Show.text == T){text(PClim.sc.site, sub("M","", row.names(MP)), cex=.5, pos = 3)}            # textes sites
    points(PClim.sc.site, pch = Shape.eco, cex = Dot.size, col = Col.eco)
    arrows(0,0,x1 = PClim.sc[,1]*Vector.scale, y1 = PClim.sc[,2]*Vector.scale, length = 0.05, lty = 1, col="#6e2115ff")
    arrows(0,0,x1 = clim.score[,1]*Vector.env.scale, y1 = clim.score[,2]*Vector.env.scale, length = 0.08, lty = 1, lwd = 2, col="#1c4871ff")
    text(PClim.sc*Vector.scale, colnames(MP), cex=.95, pos = 3, col = "#6e2115ff")                # textes taxa
    text(clim.score*Vector.env.scale, row.names(clim.score), cex=1.3, pos = 3, col = "#1c4871ff")      # textes clim
    if(Title.inside == T){
      usr <- par("usr")   # save old user/default/system coordinates
      par(usr = c(0, 1, 0, 1)) # new relative user coordinates
      text(0.01, 0.95, Annot, adj = 0, cex = 1.7)  # if that's what you want
      par(usr = usr) # restore original user coordinates
    }  
    
    #### Legend #####
    if(is.null(Cluster.path) == F & Display.legends == T){
      legend(Leg.loc,
             legend = levels(Fact.eco),
             col = diamond_color_colors,
             pch = 19, cex = 0.9,
             y.intersp = 0.85,	          # espace entre y
             x.intersp = 0.6,           # espace entre x
             bty = "n" )
    }
    
    #### Legend 2 ####
    legend(Leg.loc2,
           legend = bquote(italic(R)^2  == ~ .(format(R2, digits = 2)) ~ ", n = " ~ .(nrow(MP))),
           pch = NA,
           y.intersp = 0.7,	                                     # espace entre y
           x.intersp = 0.4,                                      # espace entre x
           bty = "n" )
    
    #### Add symbol ####
    if(is.null(Symbol.path) == F){
      if(is.null(Symbol.loc)== T){Symbol.loc <- c(0.83,0.83,0.99,0.99)}
      if(grepl("\\.png", Symbol.path)){
        library(png)
        library(grid)
        pic <- readPNG(Symbol.path)
        usr <- par("usr")
        par(usr = c(0, 1, 0, 1))
        rasterImage(pic,Symbol.loc[1], Symbol.loc[2], Symbol.loc[3], Symbol.loc[4])
        par(usr = usr)
      }
      
      if(grepl("\\.xml", Symbol.path)){
        library(grImport)
        pic <- readPicture(Symbol.path)
        usr <- par("usr")
        par(usr = c(0, 1, 0, 1))
        picture(pic,Symbol.loc[1], Symbol.loc[2], Symbol.loc[3], Symbol.loc[4])
        par(usr = usr)
      }
    }
    
  }
  #### Export data ####
  if(is.null(Save.path) == F){
    Site.name <- gsub(" ","_",Site.name)
    Save.path.Site <- gsub("\\.csv", "_RDA_site.csv", Save.path)
    Save.path.Taxon <- gsub("\\.csv", "_RDA_taxa.csv", Save.path)
    write.table(t(PClim.sc.site), file = Save.path.Site, col.names = FALSE, sep = ",")
    write.table(t(PClim.sc), file = Save.path.Taxon, col.names = TRUE, sep = ",")
  }
  
  if(is.null(Save.plot) == F){
    if(ggplot.display == T){
      ggsave(p, file = Save.plot, width = W*0.026458333, height = H*0.026458333, units = "cm")
    }
    else{dev.off()}
  }
  
  if(return.VIF == T){return(round(my_VIF, digits = 1))}
  else{
    if(Result.vector == T){return(PClim.sc)}
    if(return.pick == T){return(p)}
    if(Result.clim == T){return(clim.score)}
    else{return(MP.pca)}
    
  }
}


PCoI.vegetation <- function(MV, MP, Mclim.MP = NULL, Mclim.MV = NULL, Meco = NULL, H, W,  Dot.opac = 0.6, Dot.size = 2.5, Legend.pos = "none", PCA.display = T,
                            Stats.pos = NULL, Show.outliers = F, Show.site.lab = F, Eco.col = NULL, Show.grp.effectif = F, PCoI.Only.vectors = F,
                            Symbol.path = NULL, Symbol.pos = NULL, Title = NULL,  Eco.lab = NULL, Manu.lim.x = NULL, Manu.lim.y = NULL,
                            Helinger.trans = F, Scale.PCA = 1, GDGT = F, Choose.clim = NULL, Text.size = 2,
                            Save.plot, Reverse.dim = F, Show.errors = T, Return.plot = F){
  #### Settings ####
  library(ggrepel) # nom des lignes à côté
  library(ggnewscale) # function new_scale_color
  library(missMDA)
  if(missing(H)){H = 400}
  if(missing(W)){W = 800}
  if(missing(Save.plot)){Save.plot = NULL}
  
  #### Fullfill NA's ####
  if(is.null(Mclim.MV) == T | is.null(Mclim.MP) == T){
    MP <- data.frame(t(data.frame(MP)))
    if(nlevels(as.factor(is.na(MP))) >= 2){
      print("NA's fullfilled for the PCA for MP.")
      nb <- estim_ncpPCA(MP, ncp.max=5) ## Time consuming, nb = 2
      MP.comp <- imputePCA(MP, ncp = nb[[1]])
      MP <- data.frame(t(data.frame(MP.comp$completeObs)))
    }
    else{MP <- data.frame(t(data.frame(MP)))}
    
    MV <- data.frame(t(data.frame(MV)))
    if(nlevels(as.factor(is.na(MV))) >= 2){
      print("NA's fullfilled for the PCA for MV.")
      nb <- estim_ncpPCA(MV, ncp.max=5) ## Time consuming, nb = 2
      MV.comp <- imputePCA(MV, ncp = nb[[1]])
      MV <- data.frame(t(data.frame(MV.comp$completeObs)))
    }
    else{MV <- data.frame(t(data.frame(MV)))}
  }
  
  #### PCA calculations ####
  if(PCoI.Only.vectors == F){
    MP <- MP[,which(names(MP) %in% names(MV))]
    MV <- MV[,which(names(MV) %in% names(MP))]
    MP <- MP[,match(names(MV), names(MP))]
  }
  
  if(PCA.display == T){
    par(mfrow=c(1,2))
    if(GDGT == T){
      PCA.MV <- PCA.pollen.surf(MV, Csv.sep =",", transp_OK = T, Scale.PCA = Scale.PCA, Helinger.trans = Helinger.trans, Site.name = "MV", Display.plot = T, Result.vector = T, GDGT = GDGT)
      PCA.MP <- PCA.pollen.surf(MP, Csv.sep =",", transp_OK = T, Scale.PCA = Scale.PCA, Helinger.trans = Helinger.trans, Site.name = "MP", Display.plot = T, Result.vector = T, GDGT = GDGT)}
    else{
      PCA.MV <- PCA.vegetation(MV, Csv.sep =",", transp_OK = T, Scale.PCA = Scale.PCA, Site.name = "MV", PCA.display = T)
      PCA.MP <- PCA.vegetation(MP, Csv.sep =",", transp_OK = T, Scale.PCA = Scale.PCA, Site.name = "MP", PCA.display = T)}
    
    if(is.null(Mclim.MV) == F & is.null(Mclim.MP) == F){
      PCA.MV <- RDA.pollen.surf(MV, MClim = Mclim.MV, Csv.sep =",", transp_OK = T, Helinger.trans = Helinger.trans, Scale.sites = Scale.PCA, Scale.taxa = Scale.PCA, Site.name = "MV", Display.plot = T, Result.vector = T, GDGT = GDGT, Choose.clim = Choose.clim)
      PCA.MP <- RDA.pollen.surf(MP, MClim = Mclim.MP, Csv.sep =",", transp_OK = T, Helinger.trans = Helinger.trans, Scale.sites = Scale.PCA, Scale.taxa = Scale.PCA, Site.name = "MP", Display.plot = T, Result.vector = T, GDGT = GDGT, Choose.clim = Choose.clim)
    }
  }
  else{
    if(GDGT == T){
      PCA.MV <- PCA.pollen.surf(MV, Csv.sep =",", transp_OK = T, Scale.PCA = Scale.PCA, Helinger.trans = Helinger.trans, Site.name = "MV", Display.plot = F, Result.vector = T, GDGT = GDGT)
      PCA.MP <- PCA.pollen.surf(MP, Csv.sep =",", transp_OK = T, Scale.PCA = Scale.PCA, Helinger.trans = Helinger.trans, Site.name = "MP", Display.plot = F, Result.vector = T, GDGT = GDGT)}
    else{
      PCA.MV <- PCA.vegetation(MV, Csv.sep =",", transp_OK = T, Scale.PCA = Scale.PCA, Site.name = "MV", PCA.display = F)
      PCA.MP <- PCA.vegetation(MP, Csv.sep =",", transp_OK = T, Scale.PCA = Scale.PCA, Site.name = "MP", PCA.display = F)}
    
    if(is.null(Mclim.MV) == F & is.null(Mclim.MP) == F){
      PCA.MV <- RDA.pollen.surf(MV, MClim = Mclim.MV, Csv.sep =",", transp_OK = T, Helinger.trans = Helinger.trans, Scale.sites = Scale.PCA, Scale.taxa = Scale.PCA, Site.name = "MV", Display.plot = F, Result.vector = T, GDGT = GDGT, Choose.clim = Choose.clim)
      PCA.MP <- RDA.pollen.surf(MP, MClim = Mclim.MP, Csv.sep =",", transp_OK = T, Helinger.trans = Helinger.trans, Scale.sites = Scale.PCA, Scale.taxa = Scale.PCA, Site.name = "MP", Display.plot = F, Result.vector = T, GDGT = GDGT, Choose.clim = Choose.clim)
      PCA.MV.clim <- RDA.pollen.surf(MV, MClim = Mclim.MV, Csv.sep =",", transp_OK = T, Helinger.trans = Helinger.trans, Scale.sites = Scale.PCA, Scale.taxa = Scale.PCA, Site.name = "MV", Display.plot = F, Result.clim = T, GDGT = GDGT, Choose.clim = Choose.clim)
      PCA.MP.clim <- RDA.pollen.surf(MP, MClim = Mclim.MP, Csv.sep =",", transp_OK = T, Helinger.trans = Helinger.trans, Scale.sites = Scale.PCA, Scale.taxa = Scale.PCA, Site.name = "MP", Display.plot = F, Result.clim = T, GDGT = GDGT, Choose.clim = Choose.clim)
    }
  }
  
  #### Procrustes calculations (full PCA) ####
  pro <- procrustes(X = PCA.MV, Y = PCA.MP, symmetric = F)
  PROTEST <- protest(PCA.MV, PCA.MP)
  ctest <- data.frame(rda1 = pro$Yrot[,1],
                      rda2 = pro$Yrot[,2],
                      xrda1 = pro$X[,1],
                      xrda2 = pro$X[,2])
  
  if(is.null(Mclim.MV) == F & is.null(Mclim.MP) == F){
    pro.clim <- procrustes(X = PCA.MV.clim, Y = PCA.MP.clim, symmetric = F)
    ctest.clim <- data.frame(rda1 = pro.clim$Yrot[,1],
                             rda2 = pro.clim$Yrot[,2],
                             xrda1 = pro.clim$X[,1],
                             xrda2 = pro.clim$X[,2])
    ctest$Ecosystem <- "Vectors"
    ctest.clim$Ecosystem <- "Climat vectors"
    ctest <- rbind(ctest, ctest.clim)
  }
  
  if(is.null(Meco) == F){
    Meco <- Meco[names(MV),]
    ctest$Ecosystem <- as.factor(Meco[which(row.names(Meco) %in% row.names(ctest)), "Ecosystem"])
  }
  else{
    if(PCoI.Only.vectors == T & is.null(Mclim.MV) == T & is.null(Mclim.MP) == T){ctest$Ecosystem <- "Vectors"}
    if(PCoI.Only.vectors == F & is.null(Mclim.MV) == T & is.null(Mclim.MP) == T){ctest$Ecosystem <- "Sites"}}
  
  if(PCoI.Only.vectors == T){My_shape <- 21}
  else{My_shape <- 19}
  
  #### Color ecosystems ####
  values.bi = c("Light taiga" = "#1874CD",
                "Sites" = "#658D94",
                "Vectors" = "royalblue", "Climat vectors" = "darkred",
                "Dark taiga" = "#658D94",
                "Steppe-desert" = "#DD5925",
                "Desert-steppe" = "#DD5925",
                "Desert" = "#CD2626",
                "Steppe" = "#EE8D25",
                "Alpine meadow" = "#FFC125",
                "Steppe" = "#1874CD",
                "Desert" = "firebrick3", 
                "Desert-steppe" = "darkorange", 
                "Mountain steppes meadows"= "darkorange", 
                "Forest-steppes"= "#6789CE", "Juniperus woodland"= "dodgerblue3",
                "Forest" = "darkgreen", "Riparian forest" = "darkgreen", 
                "Steppe" = "goldenrod1",
                "Forest-steppe" = "#38A700", 
                "Woodland" = "#9FBB67",
                "Alpine steppe" = "dodgerblue3",
                "Anthropic" = "grey10",
                "Steppe-forest" = "#B2A75C",
                "Forest-steppe" = "#B2A75C", 
                "Ch'ol cold desert-steppes" = "#7916C4",
                "Tugai riparian forest" = "#BB0268",
                "Ch'ol warm deserts" = "#bb0202",
                "Adyr desert-steppes" = "#ff5400",
                "Adyr steppes" = "#e6c607",
                "Tau riparian forest" = "#2C9740",
                "Tau thermophilous woodlands" = "#85682D",
                "Tau juniper steppe-forest" = "#176E5B",
                "Tau steppes" = "#bab133",
                "Alau cryophilous steppe-forest" = "#54a697",
                "Alau meadows" = "#197CDA",
                "Tugai riparian forests" = "#7916C4",
                "Ch'ol deserts" = "#bb0202",
                "Adyr pseudosteppes" = "#ff5400",
                "Adyr steppes" = "#ECED8A",
                "Tau xeric shrublands" = "#85682D",
                "Tau open woodlands" = "#1e8736",
                "Tau steppes" = "#b0a62e",
                "Alau cryophilous open woodlands" = "#54a697",
                "Alau mesic grasslands" = "#197CDA" 
  )
  
  ctest$Lab <- row.names(ctest)
  
  if(is.null(Eco.col) == F & is.null(Meco) == F){
    Eco.col <- Eco.col[unique(Meco$Ecosystem)]
    Scale.color <- scale_color_manual(values = Eco.col)
    Scale.fill <- scale_fill_manual(values = Eco.col)
    
    if(is.null(Eco.lab) ==F){
      Scale.color <- scale_color_manual(values = Eco.col, labels = Eco.lab, name = "Vegetation-types")
      Scale.fill <- scale_fill_manual(values = Eco.col, labels = Eco.lab, name = "Vegetation-types")
    }
  }
  else{
    Scale.color <- scale_color_manual(values = values.bi, guide = "none")
    Scale.fill <- scale_fill_manual(values = values.bi, guide = "none")
  }
  
  if(Show.outliers == T){Col.out = "grey60"}
  if(Show.outliers == F){Col.out = NA}
  
  #### Other param graph + stats ####
  if(Reverse.dim == T){
    My_labs <- labs(x = substitute(paste("PCoI"[2])),y = substitute(paste("PCoI"[1])))
    names(ctest) <- names(ctest)[c(2,1,4,3,5,6)]
  }
  else{My_labs <- labs(x = substitute(paste("PCoI"[1])),y = substitute(paste("PCoI"[2])))}
  
  if(is.null(Title) == F){Mytit <- ggtitle(Title)}
  else{Mytit <- NULL}
  
  if(is.null(Stats.pos) == F){
    Note.n <- annotate("text", x = Stats.pos[1], y = Stats.pos[2],
                       label = paste("list(italic(m[12]^2) == ", round(pro$ss, digits = 2), ", n == ", ncol(MP), ")", sep = ""), size = 3, color = "grey10", hjust = 0, vjust = 1, parse = T)
    
    Note.n2 <- annotate("text", x = Stats.pos[1], y = Stats.pos[2]-.08,
                        label = paste("list(PROTEST~italic(r) == ", round(PROTEST$t0, digits = 2),
                                      ", italic(p) > ", round(PROTEST$signif, digits = 3), ")", sep = ""), size = 3, color = "grey10", hjust = 0, vjust = 1, parse = T)
    
  }
  else{Note.n <- NULL}
  
  if(Show.site.lab == T){Site.lab <- geom_label_repel(aes(x=xrda1, y=xrda2, label = Lab, color = Ecosystem), size = Text.size)}
  else{Site.lab <- NULL}
  
  if(is.null(Manu.lim.y)==F){Lim.y <- ylim(Manu.lim.y)}
  else{Lim.y <- NULL}
  if(is.null(Manu.lim.x)==F){Lim.x <- xlim(Manu.lim.x)}
  else{Lim.x <- NULL}
  
  #### Save plots ####
  if(is.null(Save.plot) == F){
    Path.to.create <- gsub("(.*/).*\\.pdf.*","\\1", Save.plot)
    dir.create(file.path(Path.to.create), showWarnings = FALSE)
    if(is.null(W) == F & is.null(H) == F){
      pdf(file = Save.plot, width = W*0.01041666666667, height = H*0.01041666666667)}
    else{pdf(file = Save.plot)}}
  
  #### Ajout symbole ####
  if(is.null(Symbol.path) == F){
    if(is.null(Symbol.pos)== T){Symbol.pos <- c(.9, .9, .16)}
    if(grepl("\\.png", Symbol.path)){
      library(png)
      library(grid)
      img <- readPNG(Symbol.path)
      g <- rasterGrob(x = Symbol.pos[1], y = Symbol.pos[2], width = Symbol.pos[3], height = Symbol.pos[3], img, interpolate = T)
    }
    
    if(grepl("\\.xml", Symbol.path)){
      library(grImport)
      img <- readPicture(Symbol.path)
      g <- pictureGrob(x = Symbol.pos[1], y = Symbol.pos[2], width = Symbol.pos[3], height = Symbol.pos[3], img)
    }
    
    Logo <- annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)
  }
  else{Logo <- NULL}
  
  #### Plot 1 ####
  p <- ggplot(ctest) +
    geom_vline(xintercept = 0, lty = "dashed")+
    geom_hline(yintercept = 0, lty = "dashed")+
    geom_segment(aes(x=rda1, y=rda2, xend=xrda1, yend=xrda2, color = Ecosystem), alpha = 0.6, arrow = arrow(length=unit(0.25,"cm")))+
    geom_point(aes(x=rda1, y=rda2, color = Ecosystem), size = Dot.size, alpha = Dot.opac, shape = My_shape) +
    My_labs + Mytit + Site.lab +
    Scale.color+Scale.fill+
    new_scale_color() + 
    Logo + Note.n + Note.n2 + Lim.y + Lim.x + 
    theme_classic()+
    theme(axis.line = element_blank(), legend.background = element_blank(),
          strip.background = element_blank(), plot.background = element_blank(), panel.background = element_blank(), 
          strip.placement = "outside", legend.position = Legend.pos,
          panel.border = element_rect(NA, "black", linewidth = 1)
    )
  
  #### Procrustes errors ####
  ctest$Procrustes_error <- sqrt((ctest$rda1-ctest$xrda1)^2+(ctest$rda2-ctest$xrda2)^2)
  ctest <- ctest[order(ctest$Procrustes_error,decreasing = T),]
  row.names(ctest) <- seq(1:nrow(ctest))
  ctest$ID <- as.numeric(row.names(ctest))
  
  if(Show.grp.effectif == T){
    ctest <- group_by(ctest, Ecosystem)
    ctest <- dplyr::mutate(ctest, N = n())
    ctest$Xmax <-  max(ctest$Procrustes_error)
    ctest$Lab <- paste(ctest$Ecosystem, ", n = ", ctest$N, sep = "")
  }
  else{ctest$Lab <- ctest$Ecosystem}
  
  p2 <- ggplot(ctest, aes(x = ID, y = Procrustes_error)) +
    geom_bar(aes(fill = Ecosystem), position = 'dodge', stat='identity')+
    geom_hline(yintercept = mean(ctest$Procrustes_error), color="grey30", lwd = 0.2)+
    geom_hline(yintercept = quantile(ctest$Procrustes_error)[c(2,4)], color="grey30", lwd = 0.2, linetype = "longdash")+
    coord_flip()+
    Scale.color+Scale.fill+ 
    theme_classic()+
    theme(axis.line.x = element_blank(), legend.position = Legend.pos, axis.ticks.x = element_blank(), axis.title = element_blank(), axis.text.x = element_blank())
  
  #### Boxplot errors #### 
  p3 <- ggplot(ctest) +
    geom_boxplot(aes(Procrustes_error, y = Lab, fill = Ecosystem), outlier.colour = Col.out)+
    geom_vline(xintercept = mean(ctest$Procrustes_error), color="grey30", lwd = 0.2)+
    geom_vline(xintercept = quantile(ctest$Procrustes_error)[c(2,4)], color="grey30", lwd = 0.2, linetype = "longdash")+
    Scale.color+Scale.fill+ 
    xlab("Procrustes error")+ 
    theme_classic()+
    theme(axis.title.y = element_blank(), legend.position = Legend.pos)
  
  #### Save and export ####
  if(Show.errors == T){p <- p|(p2 / p3 + plot_layout(widths = c(3/4, 1/4), guides = "collect")) + plot_layout(heights = c(3/5, 2/5))}
  
  if(is.null(Save.plot) == F){
    print(p)
    dev.off()}
  
  if(Return.plot == T){return(p)}
  else{return(ctest)}
}

LR.FA.clim <- function(M, Mclim, Mtype, RegLig, Add.R2, yLim, xLim, Facet = F, Compare.all = F, Stat.lm = F, Ylab = NULL, Xlab = NULL,
                       Bootstraps = F, T_test = F, Legend.position = "bottom", Segment.size = .2, R2.pos = "bottomright", R2.size = 2.7, z.size = 3, 
                       Subset.type, Select.type, Show.Plotly, Return.plot = F, Save.csv = NULL, Save.plot, W, H){
  #### Initialization values ####
  if(missing(RegLig)){RegLig = NULL}
  if(missing(Mtype)){Mtype = NULL}
  if(missing(Select.type)){Select.type = NULL}
  if(missing(Subset.type)){Subset.type = NULL}
  if(missing(Show.Plotly)){Show.Plotly = F}
  if(missing(Add.R2)){Add.R2 = F}
  if(missing(yLim)){yLim = NULL}
  if(missing(xLim)){xLim = NULL}
  if(missing(Save.plot)){Save.plot = NULL}
  if(missing(W)){W = NULL}
  if(missing(H)){H = NULL}
  
  #### Extract data ####
  M <- cbind(M, Mclim)
  Param.lab <- names(Mclim)
  names(M)[ncol(M)] <- "Param.clim"
  
  if(is.null(Mtype) == F){
    if(is.null(Select.type) == T){Select.type <- names(Mtype)[2]}
    M <- cbind(M, Type = Mtype[[Select.type]])
    
    if(is.null(Subset.type) == F){
      M <- M[M$Type %in% Subset.type,]}
  }
  else{M$Type <- "Samples"}
  
  M_MAAT <- melt(M, id = c("Param.clim", "Type"))
  M_MAAT$Facet <- M_MAAT$Type
  
  #### Ordering factors ####
  if(is.factor(M_MAAT$Type) == T){
    Keep.factor <- T
    My_levels <- levels(M_MAAT$Type)
  }
  else{Keep.factor <- F}
  
  #### Compare to all ####
  if(Compare.all == T){
    M_MAAT_all <- M_MAAT
    if(any(grepl("[[:digit:]]", unique(M_MAAT_all$Type)))){My_all <- paste(length(unique(M_MAAT_all$Type)[!is.na(unique(M_MAAT_all$Type))])+1, "All", sep = "_")}
    else{My_all <- "All"}
    
    if(Keep.factor == T){My_levels <- c(My_levels, My_all)}
    
    M_MAAT_all$Type <- My_all  
    M_MAAT_all$Facet <- My_all
    M_MAAT <- rbind(M_MAAT_all, M_MAAT)
    for(i in 1:length(M_MAAT$Type)){
      M_MAAT_all$Facet <- unique(M_MAAT$Type)[i]
      M_MAAT <- rbind(M_MAAT_all, M_MAAT)
    }
    
    M_MAAT_all <- M_MAAT_all[!is.na(M_MAAT_all$Type),]
    M_MAAT_all <- M_MAAT_all[!is.na(M_MAAT_all$Facet),]
    M_MAAT_all <- M_MAAT_all[!duplicated(M_MAAT_all),]
  }
  M_MAAT <- M_MAAT[!is.na(M_MAAT$Type),]
  M_MAAT <- M_MAAT[!is.na(M_MAAT$Facet),]
  M_MAAT <- M_MAAT[!duplicated(M_MAAT),]
  
  #### Graph param ####
  DF.line <- data.frame(Slope = 1, Int = 0, variable = unique(M_MAAT$variable))
  if(is.null(yLim) == F){
    if(typeof(yLim) == "double"){
      yLim <- ylim(yLim)
      facet_scale <- "fixed"
    }
    if(typeof(yLim) == "character"){
      facet_scale <- yLim
      yLim <- NULL
    }
  }
  else{
    facet_scale <- "fixed"
    yLim <- NULL}
  if(is.null(xLim) == F){
    if(typeof(xLim) == "double"){
      xLim <- xlim(xLim)
      facet_scale <- "fixed"
    }
    if(typeof(xLim) == "character"){
      facet_scale <- xLim
      xLim <- NULL
    }
  }
  else{
    facet_scale <- "fixed"
    xLim <- NULL}
  
  if(is.null(Ylab) == F){Ylab <- ylab(Ylab)}
  
  if(is.null(Xlab) == F){Xlab <- xlab(Xlab)}
  if(is.null(Xlab) == T){Xlab <- xlab(Param.lab)}
  
  
  #### Color & shape settings ####
  Lab.list <- c(Lake_Ayrag = "grey30", Surface = "grey30", "4_All" = "grey70", "5_All" = "grey70", "All" = "grey70",
                Lakes = "royalblue", Lacustrine = "royalblue", "Wet-Cold" = "royalblue",
                "Fibrous soil" = "#bf812e",
                "Topcore silk" = "#8c510a", "Dry-Warm" = "#8c510a",
                "1_Hyper-arid" = "#8c510a", "2_Arid" = "#bf812e", "3_Semi-arid" = "#dfc27e", "4_Dry sub-humid" = "#f5e9bf", "5_Humid" = "#80cec1",
                "Hyper-arid" = "#8c510a", "Arid" = "#bf812e", "Semi-arid" = "#dfc27e", "Dry sub-humid" = "#f5e9bf", "Humid" = "#80cec1",
                "No-data" = "grey30", "Fresh" = "royalblue","Hyposaline" = "#A9D1C2","Saline" = "darkorange", "Hypersaline" = "darkred",
                "0_No-data" = "grey30", "1_Fresh" = "royalblue","2_Hyposaline" = "#A9D1C2","3_Saline" = "darkorange", "4_Hypersaline" = "darkred",
                "0_No-data" = "grey30", "1_Fresh" = "royalblue","2_Brakish" = "#e4af08","3_Low-salted" = "darkorange", "4_Hyper-salted" = "#8c510a",
                "0_No-data" = "grey30", "1_Fresh" = "royalblue","2_Brakish" = "#e4af08","3_Salted" = "#8c510a",
                "1_Acid" = "royalblue","2_Neutral" = "#e4af08","3_Alkalin" = "#8c510a",
                "Acid" = "royalblue","Neutral" = "#e4af08","Alkalin" = "#8c510a",
                Low_S = "royalblue", High_S = "darkorange",
                Low_IR = "royalblue", High_IR = "darkorange",
                Mongolia = "#3e96bdff", Chine = "#f02a26", Uzbekistan = "#6fb440", 
                Tajikistan = "#e4af08", Russia = "#0035a9", Azerbaijan = "#094227",
                Soil = "darkorange", "Soil Dearing" = "darkorange", "Soil ACA" = "#c3510aff",
                Peat = "#74D43B", "peat Dearing" = "#74D43B", "Peat Naafs" = "#74D43B",
                Moss = "#74D43B", "Moss ACA" = "#1d5140ff",
                "Deserts & Xeric Shrublands" = "#C88282",
                "Temperate Grasslands, Savannas & Shrublands" = "#ECED8A",
                "Montane Grasslands & Shrublands" = "#D0C3A7",
                "Temperate Conifer Forests" = "#6B9A88",
                "Temperate Broadleaf & Mixed Forests" = "#3E8A70",
                "N/A" = "#FFEAAF",
                "Tundra" = "#A9D1C2",
                "Boreal Forests/Taiga" = "#8FB8E6")
  Lab.list <- Lab.list[names(Lab.list)%in% unique(M_MAAT$Type)]
  
  if(Compare.all == T){
    Shape.list <- rep(1.5, length(Lab.list))
    names(Shape.list) <- names(Lab.list)
    Shape.list[names(Shape.list) == My_all] <- .3
    My_shapes <- scale_size_manual(name = "point", values = Shape.list, guide = "none")
    My_dots <- geom_point(aes(color = Type, size = Type))
  }
  else{My_shapes <- NULL; My_dots <- geom_point(aes(color = Type))}
  #### R2 position ####
  if(R2.pos == "bottomleft"){
    R2.y = "bottom"
    R2.x = "left"}
  if(R2.pos == "bottomright"){
    R2.y = "bottom"
    R2.x = "right"}
  if(R2.pos == "topright"){
    R2.y = "top"
    R2.x = "right"}
  if(R2.pos == "topleft"){
    R2.y = "top"
    R2.x = "left"}
  if(R2.pos == "none"){
    R2.y = "none"
    R2.x = "none"}
  
  #### Select Reg Lin type ####
  if(is.null(RegLig) == T){
    RegLin <- NULL
    Add.r2 <- NULL
  }
  else{
    if(RegLig == "Global"){RegLin <- geom_smooth(method = "lm", se = F, size = .8, color = "grey20", linetype = "longdash", formula = 'y ~ x')}
    if(RegLig == "Local"){RegLin <- geom_smooth(method = "lm", se = F, size = .8, aes(color = Type), formula = 'y ~ x')}
    
    if(Add.R2 == T){
      library(ggpmisc)
      if(RegLig == "Global"){
        Add.r2 <- stat_poly_eq(label.y = R2.y, label.x = R2.x, size = R2.size, small.r = F, 
                               # aes(label = sprintf("%s*\", \"*%s" , after_stat(rr.label), after_stat(p.value.label)))
                               aes(label = after_stat(rr.label))
        )}
      
      if(RegLig == "Local"){
        Add.r2 <- stat_poly_eq(label.y = R2.y, label.x = R2.x, size = R2.size, small.r = F,
                               # aes(label = sprintf("%s*\", \"*%s" , after_stat(rr.label), after_stat(p.value.label)))
                               aes(color = Type, label = sprintf("%s*\", \"*%s" , after_stat(adj.rr.label), after_stat(n.label)))
                               # aes(color = Type, label = after_stat(adj.rr.label))
        )}
    }
    else{Add.r2 <- NULL}
  }
  
  #### My_facet #####
  if(Facet == F){
    My_facet <- facet_wrap(vars(variable), scales = facet_scale)
    My_striplab <- element_text()
    
  }
  else{
    if(Keep.factor == T){M_MAAT$Facet <- factor(M_MAAT$Facet, My_levels, ordered = T)}
    My_facet <- facet_grid(variable~Facet, scales = facet_scale)
    My_striplab <- element_blank()
  }
  #### Plot ####
  p <- ggplot(data = M_MAAT, aes(x = Param.clim, y = value))+
    # geom_abline(inherit.aes = F, data = DF.line, aes(slope = Slope, intercept = Int), size = .8, color = "grey60", linetype = "dashed")+
    My_facet + Ylab + Xlab +
    My_dots +
    RegLin + Add.r2 + yLim + xLim +
    My_shapes + 
    scale_color_manual(name = Select.type, values = Lab.list) +
    theme(
      axis.line = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(fill = NA),
      plot.background = element_blank(), 
      legend.position = Legend.position,
      panel.grid = element_blank(), 
      legend.key = element_blank(),
      strip.text.x = element_text(size = 12),
      strip.text.y = My_striplab,
      plot.margin = unit(c(0,0,0,0), 'pt'),
      strip.background = element_blank(),
    )
  
  #### Stats lm ####
  if(Stat.lm == T){
    #### Cleaning data ####
    print("**** Let's add the t-test graphs !****")
    MGA <- M
    My_res <- setNames(data.frame(matrix(ncol = 7, nrow = 0)), c('Intercept', 'Coef', 'SE', 'SE_inter', 'n', 'R_adj', 'RMSE'))
    Keep.param <- names(MGA)[1]
    names(MGA)[1] <- "value"
    
    MGA$Type <- as.character(MGA$Type)
    if(Compare.all == T){
      MGA2 <- MGA
      MGA2$Type <- My_all
      MGA <- rbind(MGA, MGA2)}
    
    #### p-values to stars ####
    p_value_to_stars <- function(p) {
      ifelse(p <= 0.001, "****", 
             ifelse(p <= 0.01, "***", 
                    ifelse(p <= 0.05, "**", 
                           ifelse(p <= 0.1, "*", ""))))}
    
    #### Calcul LR ####
    MGA <- MGA[!is.na(MGA$Type),]
    for(i in unique(MGA$Type)){
      MGA.i <- MGA[MGA$Type == i,]
      LM <- lm(Param.clim ~ value, MGA.i)
      Res <- summary(LM)
      My_res <- rbind(My_res, c(Intercept = LM$coefficients[[1]], Coef = LM$coefficients[[2]], SE = Res$coefficients[2,2], SE_inter = Res$coefficients[1,2], n = nrow(MGA.i), R_adj = round(Res$adj.r.squared, digits = 2), RMSE = round(Res$sigma, digits = 2)))
    }
    names(My_res) <- c('Intercept', 'Coef', 'SE', 'SE_inter', 'n', 'R_adj', 'RMSE')
    row.names(My_res) <- unique(MGA$Type)
    if(Keep.factor == T){My_res <- My_res[match(row.names(My_res), My_levels),]}
    My_res$Formula <- paste(Param.lab, "_{", gsub("ine", ".", gsub(".*_", "", row.names(My_res))), "} = ", My_res$Intercept, " + ", My_res$Coef, "\times MBT'_{5Me}, (n = ", My_res$n, ", R^{2}_{adj} = ", My_res$R_adj, ", RMSE = ", My_res$RMSE, ")", sep = "")
    
    #### Calcul z ####
    compare.coeff <- function(b1,se1,b2,se2){return((b1-b2)/sqrt(se1^2+se2^2))}
    NCom <- combn(1:nrow(My_res), 2)
    
    My_res3 <- setNames(data.frame(matrix(ncol = 6, nrow = 0)), c('Begin', 'End', 'z', 'p_val', 'Cor1', 'Cor2'))
    if(Bootstraps == T){My_res3$Bootstrap <- numeric(0)}
    if(T_test == T){
      My_res3$t <- numeric(0)
      My_res3$p_val_t <- numeric(0)
    }
    
    for(i in 1:ncol(NCom)){
      #### Slopes comparaison ####
      z <-  compare.coeff(My_res$Coef[NCom[1,i]],My_res$SE[NCom[1,i]],My_res$Coef[NCom[2,i]],My_res$SE[NCom[2,i]])
      p_value <-  2*pnorm(-abs(z))
      Res.i <- c(Begin = NCom[1,i], End = NCom[2,i], z = round(z, digits = 1), p_val = p_value, Cor1 = row.names(My_res)[NCom[1,i]], Cor2 = row.names(My_res)[NCom[2,i]])
      
      #### Intercept comparation settings ####
      if(Bootstraps == T | T_test == T){
        data1 <- MGA[MGA$Type == row.names(My_res)[NCom[1,i]],]
        data2 <- MGA[MGA$Type == row.names(My_res)[NCom[2,i]],]
        intercept1 <- My_res$Intercept[NCom[1,i]]
        intercept2 <- My_res$Intercept[NCom[2,i]]
        se1 <- My_res$SE_inter[NCom[1,i]]
        se2 <- My_res$SE_inter[NCom[2,i]]
      }
      
      #### Intercept comparaison  (bootstrap) ####
      if(Bootstraps == T){
        set.seed(123)  # Set seed for reproducibility
        n_bootstrap <- 10
        bootstrap_diff <- numeric(n_bootstrap)
        
        for (j in 1:n_bootstrap) {
          data1_resample <- data1[sample(1:nrow(data1), replace = T), ]
          data2_resample <- data2[sample(1:nrow(data2), replace = T), ]
          model1_resample <- lm(Param.clim ~ value, data = data1_resample)
          model2_resample <- lm(Param.clim ~ value, data = data2_resample)
          bootstrap_diff[j] <- coef(model2_resample)[1] - coef(model1_resample)[1]
        }
        observed_diff <- intercept2 - intercept1
        p_value_bootstrap <- mean(abs(bootstrap_diff) >= abs(observed_diff))
        Res.i <- c(Res.i, Bootstrap = p_value_bootstrap)
      }
      
      #### Intercept comparaison  (t-test) ####
      if(T_test == T){
        model1 <- lm(Param.clim ~ value, data = data1)
        model2 <- lm(Param.clim ~ value, data = data2)
        intercept1 <- coef(model1)[1]
        intercept2 <- coef(model2)[1]
        se1 <- summary(model1)$coefficients[1, 2]
        se2 <- summary(model2)$coefficients[1, 2]
        
        intercept_diff <- intercept2 - intercept1
        se_diff <- sqrt(se1^2 + se2^2)
        t_stat <- intercept_diff / se_diff
        df <- min(length(data1$Param.clim) - 2, length(data2$value) - 2)
        p_value_t_test <- 2 * pt(-abs(t_stat), df)
        # p_value_t_test <- 2*pnorm(-abs(t_stat))
        Res.i <- c(Res.i, t = round(t_stat, digits = 1), p_val_t = round(p_value_t_test, digits = 3))
      }
      My_res3 <- rbind(My_res3, Res.i)
      
    }
    Clean.names <- c('Begin', 'End', 'z', 'p_val', 'Cor1', 'Cor2')
    if(Bootstraps == T){Clean.names <-  c(Clean.names, 'Bootstrap')}
    if(T_test == T){Clean.names <-  c(Clean.names, 't', 'p_val_t')}
    names(My_res3) <- Clean.names
    
    #### Table graph t-test building ####
    My_res3$Begin <- as.numeric(My_res3$Begin)
    My_res3$End <- as.numeric(My_res3$End)
    My_res3$p_val <- paste("z(a) = ", My_res3$z, sapply(My_res3$p_val, p_value_to_stars), sep = "")
    if(T_test == T){My_res3$p_val <- paste(My_res3$p_val, ", z(b) = ", My_res3$t, sapply(My_res3$p_val_t, p_value_to_stars), sep = "")}
    My_res3$xlab <- (My_res3$Begin + My_res3$End)/2
    My_res3$yline <- abs(My_res3$Begin - My_res3$End)
    My_res3$Begin[My_res3$yline == 4] <- My_res3$Begin[My_res3$yline == 4] -0.3
    My_res3$End[My_res3$yline == 4] <- My_res3$End[My_res3$yline == 4] +0.3
    My_res3$Begin[My_res3$yline == 3] <- My_res3$Begin[My_res3$yline == 3] -0.1
    My_res3$End[My_res3$yline == 3] <- My_res3$End[My_res3$yline == 3] +0.1
    My_res3$Begin[My_res3$yline == 1] <- My_res3$xlab[My_res3$yline == 1] -0.4
    My_res3$End[My_res3$yline == 1] <- My_res3$xlab[My_res3$yline == 1] +0.4
    My_res3$yline[My_res3$yline == 1] <- 0.3
    My_res3$yline[My_res3$yline == 2][c(2)] <- 1.2
    My_res3$yline[My_res3$yline == 3][c(1)] <-2.9 
    My_res3$yline[My_res3$yline == 3] <-3.8 
    My_res3$yline[My_res3$yline == 4] <- 4.6
    Ylim.max <- nrow(My_res)
    
    #### Table graph t-test building (adaptation to 5 groups) ####
    if(length(unique(MGA$Type)) >= 6){
      print("Add some space up.")  
      My_res3$yline[My_res3$yline == 5] <- 6.9
      My_res3$yline[My_res3$yline == 2][c(3)] <- 2.9
      My_res3$yline[My_res3$yline == 4.6][c(1)] <- 6.1
      My_res3$yline[My_res3$yline == 3.8][c(1)] <- 5.3
      My_res3$yline[My_res3$yline == 4.6][c(1)] <- 4.5
      My_res3$yline[My_res3$yline == 3.8][c(1)] <- 3.7
      Ylim.max <- Ylim.max + 1.2
      
      My_res3$End[My_res3$yline == 2.9][c(1)] <- 4
      My_res3$Begin[My_res3$yline == 2.9][c(2)] <- 4.1
      My_res3$Begin[My_res3$yline == 6.9] <- 0.7
      My_res3$Begin[My_res3$yline == 6.1] <- 0.85
      
    }
    
    My_res3$yline2 <- My_res3$yline - 0.3
    My_res3$Begin[My_res3$yline == 2] <- My_res3$Begin[My_res3$yline == 2] +0.05
    My_res3$End[My_res3$yline == 2] <- My_res3$End[My_res3$yline == 2] -0.05
    print(My_res3)
    
    #### Plot scale ####
    pScale <- ggplot(data = My_res3)+
      geom_segment(aes(x = Begin, xend = End, y = yline, yend = yline), linewidth = Segment.size)+ 
      geom_segment(aes(x = Begin, xend = Begin, y = yline, yend = yline2), linewidth = Segment.size)+ 
      geom_segment(aes(x = End, xend = End, y = yline, yend = yline2), linewidth = Segment.size)+ 
      scale_x_continuous(limits = c(0.5,nrow(My_res)+0.5),expand = c(0, 0)) +
      scale_y_continuous(limits = c(-0.1,Ylim.max),expand = c(0, 0)) +
      geom_label(aes(x = xlab, y = yline, label = p_val), size = z.size, vjust = 0.5, hjust = 0.5, show.legend = F, fill = "white")+
      
      theme(axis.title = element_blank(), axis.ticks = element_blank(), plot.background = element_blank(), panel.grid = element_blank(), panel.background = element_blank(),
            axis.text = element_blank(), plot.margin = unit(c(0,0,0,0), 'pt'),  strip.clip = "off")
    
    p <- pScale / p + plot_layout(heights = c(0.35,0.65))
    
    #### Save table ####
    write.table(My_res, Save.csv, sep = ",", dec = ".")
  }
  
  #### Save html ####
  if(Show.Plotly == T){
    library(plotly)
    library(htmlwidgets)
    # print(names(Mtern))
    Save.plot.html <- gsub("pdf", "html", Save.plot)
    Keep.name <- gsub(".*\\/", "", Save.plot.html)
    Path.root <- paste(gsub(Keep.name, "", Save.plot.html), "HTML_files/", sep = "")
    if(file.exists(Path.root) == F){dir.create(Path.root)}
    Save.plot.html <- paste(Path.root, Keep.name, sep = "")
    p1_ly <- ggplotly(p)
    p1_ly <- p1_ly %>% layout(boxmode = "group", boxpoints = F)
    options(warn = - 1) 
    saveWidget(p1_ly, file = Save.plot.html)}
  
  #### Save plots ####
  if(is.null(Save.plot) == F){
    if(is.null(W) == F & is.null(H) == F){
      ggsave(Save.plot, width = W*0.026458333, height = H*0.026458333, units = "cm")}
    else{ggsave(Save.plot)}}
  
  if(Return.plot == T){return(p)}
  else{
    if(Stat.lm == T){
      return(list(My_res, My_res3))}
    else{return(M_MAAT)}
  }
  
}