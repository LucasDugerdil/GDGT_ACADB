#### Import training data and script ####
source("Import/Script/brGDGT_script.R")
Meco <- readRDS("Import/Training/Meco_ACADB.Rds")
Mclim <- readRDS("Import/Training/Mclim_ACADB.Rds")
Msurf.mean <- readRDS("Import/Training/MbrGDGT_ACADB.Rds")
Meco.WDB <- readRDS("Import/Training/Meco_WDB.Rds")
Mclim.WDB <- readRDS("Import/Training/Mclim_WDB.Rds")
Msurf.mean.WDB <- readRDS("Import/Training/MbrGDGT_WDB.Rds")

#### Cleaning data ####
M.br.GDGT <- Msurf.mean[,grep("^f.I", colnames(Msurf.mean))]
M.br.GDGT <- M.br.GDGT[,!grepl("_8Me", colnames(M.br.GDGT))]
M.br.GDGT <- M.br.GDGT[,!grepl("_7Me", colnames(M.br.GDGT))]

Meco$Salinity_classe <- gsub(".*_", "", Meco$Salinity_classe)
Meco$Salinity_classe <- factor(Meco$Salinity_classe, c('Fresh', 'Hyposaline', 'Saline', "Hypersaline"))
Meco$Sample.type[Meco$Sample.type == "Moss"] <- "Soil"
Meco$Acidity <- gsub(".*_", "", Meco$Acidity)
Meco$Acidity <- factor(Meco$Acidity, c('Acid', 'Neutral', 'Alkalin'))

Mfull.WDB <- cbind(Msurf.mean.WDB, Meco.WDB, subset(Mclim.WDB, select = -c(Latitude, Longitude)))

Mfull <- cbind(Msurf.mean, Meco, subset(Mclim, select = -c(Latitude, Longitude)))

#### Select figures / tables to plot ####
Table_1 = T; Table_3 = T; Table_S1 = T
Fig_1 = T; Fig_3 = T; Fig_4 = T; Fig_5 = T; Fig_S1 = T; Fig_S2 = T; Fig_6 = T; Fig_7 = T; Fig_8 = T; Fig_S3 = T; Fig_9 = T; Fig_10 = T; Fig_S4 = T; Fig_S5 = T; Fig_S6 = T
# Fig_1 = F; Fig_3 = F; Fig_4 = F; Fig_5 = F; Fig_S1 = F; Fig_S2 = F; Fig_6 = F; Fig_7 = F; Fig_8 = F; Fig_S3 = F; Fig_9 = F; Fig_10 = F; Fig_S4 = F; Fig_S5 = F; Fig_S6 = F

#### Table 1 (ACADB settings ) #### 
if(Table_1 == T){
  #### Merge DB ####
  Keep.param <- c("MAAT", "MAF", "MAP", "AI", "Altitude", "Country")
  Keep.param2 <- c("MAAT", "MAF", "MAP", "AI", "Altitude", "DB")
  Mclim.stats <- cbind(Mclim, Meco)
  Mclim.stats <- Mclim.stats[Keep.param]
  
  Mclim.stats.WDB <- cbind(Msurf.mean.WDB, Meco.WDB, subset(Mclim.WDB, select = -c(Latitude, Longitude)))
  Mclim.stats.WDB <- Mclim.stats.WDB[Mclim.stats.WDB$DB == "World DB",]
  Mclim.stats.WDB <- Mclim.stats.WDB[Keep.param2]
  names(Mclim.stats.WDB)[ncol(Mclim.stats.WDB)] <- "Country"
  
  Mclim.stats <- rbind(Mclim.stats, Mclim.stats.WDB)
  
  #### Stats (mean + SD) calculation ####
  Ncountry <- dplyr::count(Mclim.stats, Mclim.stats$Country)
  Mclim.stats <- na.omit(Mclim.stats)
  Mclim.stats.mean <- aggregate(Mclim.stats[-c(ncol(Mclim.stats))], by = Mclim.stats["Country"], FUN = mean, na.action = na.omit)
  Mclim.stats.mean <- rbind(Mclim.stats.mean, c(NA, colMeans(Mclim.stats.mean[Mclim.stats.mean$Country != "World DB", -c(1)], na.rm = T)))
  Mclim.stats.sd <- aggregate(Mclim.stats[-c(ncol(Mclim.stats))], by = Mclim.stats["Country"], FUN = sd)
  Mclim.stats.sd <- rbind(Mclim.stats.sd, c(NA, colMeans(Mclim.stats.sd[Mclim.stats.sd$Country != "World DB", -c(1)], na.rm = T)))
  
  Mclim.stats.sd <- melt(Mclim.stats.sd, "Country")
  names(Mclim.stats.sd) <- c("Country", "Param", "sd")
  Mclim.stats.mean <- melt(Mclim.stats.mean, "Country")
  names(Mclim.stats.mean) <- c("Country", "Param", "mean")
  Mclim.stats.mean <- full_join(Mclim.stats.mean, Mclim.stats.sd, by = join_by(Country, Param))
  Mclim.stats.mean$sd <- signif(Mclim.stats.mean$sd, digits = 2)
  Mclim.stats.mean$mean <- signif(Mclim.stats.mean$mean, digits = 2)
  Mclim.stats.mean$clean <- paste(Mclim.stats.mean$mean, Mclim.stats.mean$sd, sep = " $\\pm$ ")
  Mclim.stats.mean <- suppressMessages(reshape2::dcast(Mclim.stats.mean, Country ~ Param))
  Mclim.stats.mean[is.na(Mclim.stats.mean)] <- "Total ACA"
  
  #### Mise-en-forme Table stats #### 
  Pub2 <- unique(Meco[c("Country", "Reference")])
  Pub2 <- rbind(Pub2, data.frame(Country = "World DB", Reference = "Raberg et al. (2022)$^{(a)}$"))
  Pub2 <- aggregate(Reference~Country, data = Pub2, paste0, collapse = ", ")
  row.names(Pub2) <- Pub2$Country
  Pub2$Bibtek <- paste("\\citet{", gsub("\\)", "}", gsub(" et al. \\(", "_", Pub2$Reference)), sep = "")
  Pub2$Bibtek <- gsub("\\\\citet\\{This study", "This study", Pub2$Bibtek)
  Pub2$Bibtek <- gsub("\\}/", ", ", Pub2$Bibtek)
  Pub2$Bibtek <- gsub("\\},", ", ", Pub2$Bibtek)
  Pub2$Bibtek <- gsub("-", "_", Pub2$Bibtek)
  Pub2$Bibtek <- gsub("study, ", "study, \\\\citet{", Pub2$Bibtek)
  Mclim.stats.mean <- cbind(Mclim.stats.mean, N = c(Ncountry$n, sum(Ncountry$n[-length(Ncountry$n)])), Original_publication = Pub2[Mclim.stats.mean$Country, "Reference"]) 
  Mclim.stats.mean <- Mclim.stats.mean[c(1,7,2:6,8)]
  Save.path <- "Tables/Table_1.csv"
  write.table(Mclim.stats.mean, Save.path, row.names = F)
  Mclim.stats.mean$Original_publication <- Pub2[Mclim.stats.mean$Country, "Bibtek"]
  library(xtable)
  LateX.caption <- "Presentation of the different dataset used in this study with their associated average climate parameters, data description (covered countries, dataset size, date of acquisition) and original publications."
  Save.path.tex <- gsub("\\.csv", "\\.tex", Save.path)
  Mclim.stats.mean$N[which(Mclim.stats.mean$Original_publication == "This study")] <- paste(Mclim.stats.mean$N[which(Mclim.stats.mean$Original_publication == "This study")], "*", sep = "")
  #print(Mclim.stats.mean)  
  Tlatex <- xtable(Mclim.stats.mean, caption = LateX.caption, type = "latex", label = "Table_ACA")
  print(Tlatex, file = Save.path.tex, booktabs = T, include.rownames = F, comment = F,
        caption.placement = "top", sanitize.text.function = function(x){x}, hline.after = )
}

#### Fig. 1 ####
if(Fig_1 == T){
  Mclim.bio <- cbind(Mclim, Meco)
  BioCl.pca.ACA.gdgt <- PCA.bioclim(Mclim.bio[c(1:9,44,68)], transp_OK = T, Scale.PCA = 4,
                                    Cluster.core = "Aridity", Shape = 21, Legend.size = 7, Dot.size = 4,
                                    Site.name = "", Num.facet = "", return.pick = F, Legend.position = "none",
                                    Save.plot = "Figures/Fig_1.pdf", H = 370, W = 370)
}

#### Fig. 3 f[br-GDGT] FA distribution ####
if(Fig_3 == T){
  Meco$Aridity2 <- gsub("_", ". ", Meco$Aridity)
  Boxplot.ACADB <- GDGT.histo.plot.surf.core(Msurf = Msurf.mean,
                                             Mtype = Meco, Select.type = "Aridity2", #Leg.pos = c(0.65,0.82), 
                                             Leg.iso = T, Return.plot = T, Leg.pos = "left",
                                             Iso.GDGT = F, Remove.8Me = T, Remove.7Me = F, Overlap.OK = T, Show.Plotly = F,
                                             Color.choice = c("#8c510a", "#bf812e", "#dfc27e", "#f5e9bf", "#80cec1"),
                                             Zoom1.comp = c(4:8), Box.linewidth = 0.1, Zoom1.Ymax = 4, Zoom2.comp = c(14:16), Zoom2.Ymax = 2.5, Insert2.loc = c(11,16.5,20,42), Insert1.loc = c(2.5,8.5,20,42),
                                             Ymax = 45, W = 1100, H = 400,Save.path = "Figures/Fig_3.pdf")
  
  
}

#### Fig. 4 Ternary diagram methylation ####
if(Fig_4 == T){
  Meco$Aridity2 <- gsub("_", ". ", Meco$Aridity)
  Mtern.ACA.soil <- Diag.ternaire.methylation(MGDGT = Msurf.mean[Meco$Sample.type == "Soil",],
                                              Mcol = Meco[Meco$Sample.type == "Soil","Aridity2"], Return.plot = T, Full.labels = F,
                                              Annot = "A1", Show.arrows = F, Add.facet = F,
                                              Show.Dearing = F, Show.Naafs.peat = F,  Show.soil = T, Show.Plotly = F, Show.lake = F, Remove.ACA = T,
                                              W = 1000, H = 1000, Alpha.dot = 0.8, Size.dot = 2.5, Export.to.chart.studio = F)
  
  Mtern.ACA.soil.facet <- Diag.ternaire.methylation(MGDGT = Msurf.mean[Meco$Sample.type == "Soil",],
                                                    Mcol = Meco[Meco$Sample.type == "Soil","Aridity2"], Return.plot = T, Full.labels = F,
                                                    Annot = "A2", Show.arrows = F, Add.facet = T,
                                                    Show.Dearing = F, Show.Naafs.peat = F,  Show.soil = T, Show.Plotly = F, Show.lake = F, Remove.ACA = T,
                                                    W = 1000, H = 1000, Alpha.dot = 0.8, Size.dot = 2.5, Export.to.chart.studio = F)
  
  Mtern.ACA.lacustrine <- Diag.ternaire.methylation(MGDGT = Msurf.mean[Meco$Sample.type == "Lacustrine",], 
                                                    Mcol = Meco[Meco$Sample.type == "Lacustrine","Aridity2"], Return.plot = T, Full.labels = F, Add.facet = F,
                                                    Annot = "B1", Show.Dearing = F, Show.Naafs.peat = F, Show.Plotly = F, Show.lake = T, Show.arrows = F,
                                                    W = 1000, H = 1000, Alpha.dot = 0.8, Size.dot = 2.5, Export.to.chart.studio = F, Remove.ACA = T)
  
  Mtern.ACA.lacustrine.facet <- Diag.ternaire.methylation(MGDGT = Msurf.mean[Meco$Sample.type == "Lacustrine",], 
                                                          Mcol = Meco[Meco$Sample.type == "Lacustrine","Aridity2"], Return.plot = T, Full.labels = F, Add.facet = T,
                                                          Annot = "B2", Show.Dearing = F, Show.Naafs.peat = F, Show.Plotly = F, Show.lake = T, Show.arrows = F,
                                                          W = 1000, H = 1000, Alpha.dot = 0.8, Size.dot = 2.5, Export.to.chart.studio = F, Remove.ACA = T)
  
  
  W = 1300; H = 1000; Save.path = "Figures/Fig_4.pdf"
  Fig.4 <- (((Mtern.ACA.soil + Mtern.ACA.soil.facet) / (Mtern.ACA.lacustrine + Mtern.ACA.lacustrine.facet)) + plot_layout(guides = "collect"))&
    theme(axis.title = element_blank(), plot.background = element_blank(), legend.background = element_blank(), plot.caption = element_blank(),
          plot.margin = unit(c(0,0,0,0),"cm"), panel.background = element_blank(), plot.subtitle = element_blank())
  ggsave(Save.path, Fig.4, width = W*0.026458333, height = H*0.026458333, units = "cm")
}

#### Fig. 5 PCA / RDA br-GGT vs. climat ####
if(Fig_5 == T){
  #### Import function / datas ####
  Meco$Aridity2 <- gsub("_", "\\. ", Meco$Aridity)
  Cluster.path = Meco ; Cluster.groups = "Aridity2"; Myshapes = "Sample.type"
  Ordin.aridity2 <- c("1. Hyper-arid", "2. Arid", "3. Semi-arid", "4. Dry sub-humid", "5. Humid")
  
  Meco.RDA <- Meco[complete.cases(Meco[c("Salinity", "pH")]),]
  Meco.RDA$Sample.type[Meco.RDA$Sample.type == "Moss"] <- "Soil"
  
  M.br.GDGT.lac <- M.br.GDGT[Meco.RDA$Sample.type == "Lacustrine",]
  M.br.GDGT.soil <- M.br.GDGT[Meco.RDA$Sample.type == "Soil",]
  Meco.RDA.soil <- Meco.RDA[Meco.RDA$Sample.type == "Soil",]
  Meco.RDA.lac <- Meco.RDA[Meco.RDA$Sample.type == "Lacustrine",]
  
  Meco3 <- data.frame(cbind(subset(Meco, select = c(pH, Salinity)), Mclim))
  Meco3 <- Meco3[!is.na(Meco3$Salinity),]
  
  M1 <- left_join(rownames_to_column(M.br.GDGT.soil), rownames_to_column(Meco))
  M2 <- left_join(rownames_to_column(M.br.GDGT.lac), rownames_to_column(Meco))
  
  #### Procrutes (vecteurs PCA) ####
  PCoI <- PCoI.vegetation(MV = data.frame(t(M.br.GDGT.soil)), MP = data.frame(t(M.br.GDGT.lac)),
                          Show.errors = F, Show.outliers = F, Show.site.lab = T, PCoI.Only.vectors = T,
                          PCA.display = F, Helinger.trans = T, Scale.PCA = 3, GDGT = T, Return.plot = T,
                          Stats.pos = c(0.02,0.6), Text.size = 3.5)
  
  #### Procrutes (vecteurs RDA) ####
  PCoI.2 <- PCoI.vegetation(MV = data.frame(t(M.br.GDGT.soil[intersect(row.names(M.br.GDGT.soil), row.names(Meco.RDA.soil)),]), check.names = F), Mclim.MV = Meco3[intersect(row.names(Meco3), row.names(Meco)),],
                            MP = data.frame(t(M.br.GDGT.lac[intersect(row.names(M.br.GDGT.lac), row.names(Meco.RDA.lac)),]), check.names = F), Mclim.MP = Meco3[intersect(row.names(Meco3), row.names(Meco)),],
                            Choose.clim = c("MAAT", "MAF", "AI", "Salinity", "pH", "MPWAQ"),
                            Show.errors = F, Show.outliers = F, Show.site.lab = T, PCoI.Only.vectors = T,
                            PCA.display = F, Helinger.trans = T, Scale.PCA = 2, GDGT = T, Return.plot = T,
                            Stats.pos = c(0.02,0.8), Text.size = 3.5)
  
  #### PCA soil ####
  PCA.soil <- PCA.bioclim(M1[c(grep("^f.", names(M1)), 46)], Dot.opac = 0.7, Dot.size = 3, Cluster.core = "Aridity2",
                          Ellipse = F, Cluster.core.lab = "",  Legend.size = 11, GDGT = T, Helinger.trans = T,
                          Density.contour = F, Opa.range = c(0.1,.4), Density.type = "polygon", 
                          transp_OK = F, Scale.PCA = 6, return.pick = T, Num.facet = "(A)", Show.annot = T, Show.site.lab = F,
                          Site.name = "PCA brGDGT (soil)", Legend.position = "none", Show.centroid = F, Marg.density.plot = T)
  #### PCA soil ####
  PCA.lac <- PCA.bioclim(M2[c(grep("^f.", names(M2)), 46)], Dot.opac = 0.7, Dot.size = 3, Cluster.core = "Aridity2",
                         Ellipse = F, Cluster.core.lab = "", Legend.size = 6, Shape = 18,
                         Density.contour = F, Opa.range = c(0.1,.4), Density.type = "polygon", Show.annot = T, Show.site.lab = F,
                         transp_OK = F, Scale.PCA = 6, return.pick = T, Num.facet = "(B)", Reverse.dim = F, GDGT = T, Helinger.trans = T,
                         Site.name = "PCA brGDGT (lacustrine)", Legend.position = "none", Show.centroid = F, Marg.density.plot = T)
  
  #### RDA ####
  VIF.full.soil <- RDA.pollen.surf(data.frame(t(M.br.GDGT.soil[intersect(row.names(M.br.GDGT.soil), row.names(Meco.RDA)),]), check.names = F), MClim = Meco3[intersect(row.names(Meco3), row.names(Meco)),],
                                   Choose.clim = c("MAAT", "MAF", "AI", "Salinity", "pH", "MPWAQ", "Altitude", "MPCOQ", "MTCOQ", "MTWAQ", "MAP"), Display.plot = F,
                                   Remove.7Me = T, Csv.sep =",", transp_OK = T, Helinger.trans = T, VIF = T, Display.VIF = F, return.VIF = T)
  
  RDA.soil <- RDA.pollen.surf(data.frame(t(M.br.GDGT.soil[intersect(row.names(M.br.GDGT.soil), row.names(Meco.RDA.soil)),]), check.names = F), MClim = Meco3[intersect(row.names(Meco3), row.names(Meco)),],
                              Choose.clim = c("MAAT", "MAF", "AI", "Salinity", "pH", "MPWAQ"),
                              Cluster.path = Meco.RDA.soil, Sort.shape = c("Soil", "Lacustrine"),
                              Cluster.groups = Cluster.groups, Shape.groups = Myshapes, ggplot.display = T, Marg.density.plot = T,
                              Display.legends = F, Remove.7Me = T, Simple.title = T, Show.text = F,
                              Csv.sep =",", transp_OK = T, Helinger.trans = T, Scale.sites = 2, Scale.taxa = 2,
                              Manu.lim = c(-1.05,1.05,-1.05,1.05),
                              Dot.size = 3, Alpha.dot = .7, Vector.env.scale = 0.9, Vector.scale = 1.5,
                              Annot = "(D)",  Sort.eco = Ordin.aridity2, VIF = T, return.VIF = F, Display.VIF = F,  return.pick = T,
                              Color.choice = c("#8c510a", "#bf812e", "#dfc27e", "#f5e9bf","#80cec1"),
                              Save.path = "Results/MGDGT_mong.csv",
                              Type.samples = "brGDGT (soil)", GDGT = T)
  
  VIF.full.lac <- RDA.pollen.surf(data.frame(t(M.br.GDGT.lac[intersect(row.names(M.br.GDGT.lac), row.names(Meco.RDA.lac)),]), check.names = F), MClim = Meco3[intersect(row.names(Meco3), row.names(Meco)),],
                                  Choose.clim = c("MAAT", "MAF", "AI", "Salinity", "pH", "MPWAQ", "Altitude", "MPCOQ", "MTCOQ", "MTWAQ", "MAP"), Display.plot = F,
                                  Remove.7Me = T, Csv.sep =",", transp_OK = T, Helinger.trans = T, VIF = T, Display.VIF = F, return.VIF = T)
  
  RDA.lac <- RDA.pollen.surf(data.frame(t(M.br.GDGT.lac[intersect(row.names(M.br.GDGT.lac), row.names(Meco.RDA.lac)),]), check.names = F), MClim = Meco3[intersect(row.names(Meco3), row.names(Meco)),],
                             Choose.clim = c("MAAT", "MAF", "AI", "Salinity", "pH", "MPWAQ"),
                             Cluster.path = Meco.RDA.lac, 
                             Cluster.groups = Cluster.groups, Shape.groups = Myshapes, ggplot.display = T, Marg.density.plot = T,
                             Display.legends = F, Remove.7Me = T, Simple.title = T, Show.text = F, 
                             Csv.sep =",", transp_OK = T, Helinger.trans = T, Scale.sites = 2, Scale.taxa = 2,
                             Manu.lim = c(-1.05,1.05,-1.05,1.05),
                             Dot.size = 3, Alpha.dot = .7, Vector.env.scale = 0.9, Vector.scale = 1.5,
                             Annot = "(E)",  Sort.eco = Ordin.aridity2, VIF = T, Display.VIF = F, return.VIF = F, return.pick = T,
                             Color.choice = c("#8c510a", "#bf812e", "#dfc27e", "#f5e9bf","#80cec1"),
                             Save.path = "Results/MGDGT_mong.csv",
                             Type.samples = "brGDGT (lacustrine)", GDGT = T)
  
  #### Merge all plots ####
  layout <- 
    "
            AAAAAAAAA#DDDDDDDDD#OOOOOOOOO
            CCCCCCCCCBFFFFFFFFFEGGGGGGGGG
            CCCCCCCCCBFFFFFFFFFEGGGGGGGGG
            CCCCCCCCCBFFFFFFFFFEGGGGGGGGG
            CCCCCCCCCBFFFFFFFFFEGGGGGGGGG
            CCCCCCCCCBFFFFFFFFFEGGGGGGGGG
            CCCCCCCCCBFFFFFFFFFEGGGGGGGGG
            CCCCCCCCCBFFFFFFFFFEGGGGGGGGG
            HHHHHHHHH#KKKKKKKKK#PPPPPPPPP
            JJJJJJJJJIMMMMMMMMMLNNNNNNNNN
            JJJJJJJJJIMMMMMMMMMLNNNNNNNNN
            JJJJJJJJJIMMMMMMMMMLNNNNNNNNN
            JJJJJJJJJIMMMMMMMMMLNNNNNNNNN
            JJJJJJJJJIMMMMMMMMMLNNNNNNNNN
            JJJJJJJJJIMMMMMMMMMLNNNNNNNNN
            JJJJJJJJJIMMMMMMMMMLNNNNNNNNN
             "
  pca.full <- PCA.soil[[1]] + PCA.soil[[2]] + PCA.soil[[3]] + PCA.lac[[1]] + PCA.lac[[2]] + PCA.lac[[3]] + PCoI + 
    RDA.soil[[1]] + RDA.soil[[2]] + RDA.soil[[3]] + RDA.lac[[1]] + RDA.lac[[2]] + RDA.lac[[3]] + PCoI.2 + 
    plot_spacer() + plot_spacer() +
    plot_layout(design = layout) & theme(plot.margin = unit(c(0,0,0,0),"cm"))
  
  W = 1300; H = 1000; ggsave(filename = "Figures/Fig_5.pdf", pca.full, width = W*0.026458333, height = H*0.026458333, units = "cm")
  
}

#### Fig. S.1 (sensitivity analysis pH vs. CBT') ####
if(Fig_S1 == T){
  Mph <- Mfull[Mfull$Sample.type == "Soil", c("CBTp", "pH")]
  Mph$Classes <- NA
  E <- 0.01
  Seq.pH <- seq(4.6, 9.5, E)
  Res <- Seq.pH; Res2 <- Seq.pH; Res3 <- Seq.pH; Res4 <- Seq.pH; Res5 <- Seq.pH
  
  for(i in 1:length(Seq.pH)){
    Mph$Classes <- ifelse(Mph$pH<Seq.pH[i], "K1", "K2")
    Res[i] <- summary(lm(CBTp~pH*Classes, Mph))$r.squared
    
    Mph.K1 <- Mph[which(Mph$Classes == "K1"),]
    Res2[i] <- summary(lm(CBTp~pH, Mph.K1))$adj.r.squared
    
    Mph.K2 <- Mph[which(Mph$Classes == "K2"),]
    Res3[i] <- summary(lm(CBTp~pH, Mph.K2))$adj.r.squared
    Res4[i] <- nrow(Mph.K1)
    Res5[i] <- nrow(Mph.K2)
    
  }
  Vline <- which(Res == max(Res, na.rm = T))
  Res <- data.frame(cbind(Res, Seq.pH, Res2, Res3))
  Vline <- mean(Seq.pH[Vline])
  
  ppH <- ggplot(Res, aes(y = Res, x = Seq.pH))+
    geom_vline(xintercept = Vline, linetype = "dashed", linewidth = .6)+ xlab("pH threshold") + ylab(expression(Multiple~R^2))+
    annotate("text", x = 7.7, y = 0.25, label = paste("pH =", round(Vline, digits = 1)))+
    geom_line(aes(y = Res2, x = Seq.pH), color = "royalblue")+ geom_point(aes(y = Res2, x = Seq.pH), , color = "royalblue")+
    geom_line(aes(y = Res3, x = Seq.pH), color = "darkorange")+ geom_point(aes(y = Res3, x = Seq.pH), , color = "darkorange")+
    geom_line()+ geom_point()+ 
    theme(panel.background = element_rect(fill = NA, color = 'black'), axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(),
          panel.grid = element_blank())      
  
  pn <- ggplot(Res)+ xlab("CI threshold") +
    geom_vline(xintercept = Vline, linetype = "dashed", linewidth = .6)+ xlab("pH threshold") + ylab(expression(N[soil]))+
    geom_line(aes(y = Res4, x = Seq.pH), color = "darkorange")+ geom_point(aes(y = Res4, x = Seq.pH), , color = "darkorange", alpha = .5)+
    geom_line(aes(y = Res5, x = Seq.pH), color = "royalblue")+ geom_point(aes(y = Res5, x = Seq.pH), , color = "royalblue", alpha = .5)+
    theme(panel.background = element_rect(fill = NA, color = 'black'), panel.grid = element_blank())
  
  ppH <- ppH / pn
  
  H = 600; W = 600; Save.plot ="Figures/Fig_S1.pdf" 
  ggsave(filename = Save.plot, ppH, width = W*0.026458333, height = H*0.026458333, units = "cm")
}

#### Fig. S.2 (salinity effect on FA) #####
if(Fig_S2 == T){
  library(ggpubr)
  M7 <- Msurf.mean[grepl("_7Me", names(Msurf.mean)) & grepl("^f.", names(Msurf.mean))]
  M6 <- Msurf.mean[grepl("_6Me", names(Msurf.mean)) & grepl("^f.", names(Msurf.mean))]
  M5 <- Msurf.mean[grepl("_5Me", names(Msurf.mean)) & grepl("^f.", names(Msurf.mean))]
  M5 <- cbind(M5, Msurf.mean[intersect(names(Msurf.mean), c("f.Ia", "f.Ib", "f.Ic"))])
  Mfull$Sum5Me<- rowSums(M5)
  Mfull$Sum6Me<- rowSums(M6)
  Mfull$Sum7Me<- rowSums(M7)
  Mfull$Sum6_7Me<- rowSums(M7+M6)
  Ms <- Mfull[c("Sum5Me","Sum6Me", "Sum7Me", "Sum6_7Me", "Salinity", "AI", "Salinity_classe", "Sample.type", "Aridity")]
  Ms <- melt(Ms, id = c("Salinity", "AI", "Salinity_classe", "Aridity", "Sample.type"))
  names(Ms)[c(ncol(Ms), ncol(Ms)-1)] <- c("Index.val", "Index.lab")
  Ms <- melt(Ms, id = c("Salinity_classe", "Aridity", "Index.val", "Index.lab", "Sample.type"))
  names(Ms)[c(ncol(Ms), ncol(Ms)-1)] <- c("Param.val", "Param.lab")
  Ms$Index.lab <- gsub("_", "+", Ms$Index.lab)
  Ms$Index.lab <- paste(gsub("Sum", "sum(", Ms$Index.lab), ")", sep = "")
  Ms$Index.lab <- factor(Ms$Index.lab, levels = unique(Ms$Index.lab), ordered = T)
  P1 <- ggplot(Ms, aes(y = Index.val, x = Param.val, color = Sample.type))+
    facet_wrap(Param.lab~Index.lab, scale = "free_x", nrow = 2)+
    ylim(c(0,1.2))+ AI.class.scale +
    geom_point(alpha = .2) + geom_smooth(method = "lm", se = F, formula = 'y ~ x', linewidth = 0.7)+ 
    stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, cor.coef.name = c("r"))+
    theme_bw()+ ylab(expression(Sigma[brGDGT]~("%")))+
    xlab("Bioclimate Parameters")+
    theme(legend.position = "bottom", panel.grid = element_blank(), 
          title = element_text(size = 13))
  
  PFULL<-P1
  H = 700; W = 1000; Save.plot ="Figures/Fig_S2.pdf" 
  ggsave(filename = Save.plot, PFULL, width = W*0.026458333, height = H*0.026458333, units = "cm")
}

#### Fig. 6 (MBT'5Me and MBT'6Me) ####
if(Fig_6 == T){
  #### Plots ####
  P1 <- ggplot(Mfull, aes(y = MBTp5Me, color = Sample.type, shape = Sample.type, x = pH))+
    geom_hex(data = Mfull.WDB, mapping = aes(y = MBTp5Me, x = pH, fill = ..count..), bins = Bin.size, color = NA)+ Fill.scale+
    geom_point(aes(color = Aridity)) + geom_smooth(method = "lm", formula = 'y ~ x', se = F, linewidth = 0.7) + stat_poly_eq(size = R.size, vstep = V.step, label.y = "bottom")+ylim(c(0,1)) + ggtitle(expression(paste(MBT,"'"[5~Me])))+ annotate("text", x = min(Mfull.WDB$pH, na.rm = T), y = Annot.x, label = "(A)", hjust = 0, size = 5)
  
  P2 <- ggplot(Mfull, aes(y = MBTp6Me, color = Sample.type, shape = Sample.type, x = pH))+
    geom_hex(data = Mfull.WDB, mapping = aes(y = MBTp6Me, x = pH, fill = ..count..), bins = Bin.size, color = NA)+ Fill.scale+
    geom_point(aes(color = Aridity)) + geom_smooth(method = "lm", formula = 'y ~ x', se = F, linewidth = 0.7) + stat_poly_eq(size = R.size, vstep = V.step, label.y = "bottom")+ylim(c(0,1)) + ggtitle(expression(paste(MBT,"'"[6~Me])))+ annotate("text", x = min(Mfull.WDB$pH, na.rm = T), y = Annot.x, label = "(B)", hjust = 0, size = 5)
  
  P3 <- ggplot(Mfull, aes(y = MBTp5Me, color = Sample.type, shape = Sample.type, x = AI))+
    Log10.scale+ annotation_logticks(sides = 'b')+
    geom_hex(data = Mfull.WDB, mapping = aes(y = MBTp5Me, x = AI, fill = ..count..), bins = Bin.size, color = NA)+ Fill.scale+
    geom_point(aes(color = Aridity)) + geom_smooth(method = "lm", formula = 'y ~ x', se = F, linewidth = 0.7) + stat_poly_eq(size = R.size, vstep = V.step, label.y = "bottom", label.x = "left")+ylim(c(0,1))+ annotate("text", x = min(Mfull$AI, na.rm = T), y = Annot.x, label = "(C)", hjust = 0, size = 5)
  
  P4 <- ggplot(Mfull, aes(y = MBTp6Me, color = Sample.type, shape = Sample.type, x = AI))+
    Log10.scale+ annotation_logticks(sides = 'b')+
    geom_hex(data = Mfull.WDB, mapping = aes(y = MBTp6Me, x = AI, fill = ..count..), bins = Bin.size, color = NA)+ Fill.scale+
    geom_point(aes(color = Aridity)) + geom_smooth(method = "lm", formula = 'y ~ x', se = F, linewidth = 0.7) + stat_poly_eq(size = R.size, vstep = V.step, label.y = "top", label.x = "right")+ylim(c(0,1))+ annotate("text", x =  min(Mfull$AI, na.rm = T), y = Annot.x, label = "(D)", hjust = 0, size = 5)
  
  P5 <- ggplot(Mfull, aes(y = MBTp5Me, color = Sample.type, shape = Sample.type, x = MAAT))+
    geom_hex(data = Mfull.WDB, mapping = aes(y = MBTp5Me, x = MAAT, fill = ..count..), bins = Bin.size, color = NA)+ Fill.scale+ xlab("MAAT (°C)")+
    geom_point(aes(color = Aridity)) + geom_smooth(method = "lm", formula = 'y ~ x', se = F, linewidth = 0.7) + stat_poly_eq(size = R.size, vstep = V.step, label.y = "bottom", label.x = "right")+ylim(c(0,1))+ annotate("text", x = min(Mfull.WDB$MAAT, na.rm = T), y = Annot.x, label = "(E)", hjust = 0, size = 5)
  
  P6 <- ggplot(Mfull, aes(y = MBTp6Me, color = Sample.type, shape = Sample.type, x = MAAT))+
    geom_hex(data = Mfull.WDB, mapping = aes(y = MBTp6Me, x = MAAT, fill = ..count..), bins = Bin.size, color = NA)+ Fill.scale+ xlab("MAAT (°C)")+
    geom_point(aes(color = Aridity)) + geom_smooth(method = "lm", formula = 'y ~ x', se = F, linewidth = 0.7) + stat_poly_eq(size = R.size, vstep = V.step, label.y = "top", label.x = "right")+ylim(c(0,1))+ annotate("text", x = min(Mfull.WDB$MAAT, na.rm = T), y = Annot.x, label = "(F)", hjust = 0, size = 5)
  
  P7 <- ggplot(Mfull, aes(y = MBTp5Me, color = Sample.type, shape = Sample.type, x = Salinity))+
    Log10.scale+ annotation_logticks(sides = 'b')+
    geom_hex(data = Mfull.WDB, mapping = aes(y = MBTp5Me, x = Salinity, fill = ..count..), bins = Bin.size, color = NA)+ Fill.scale+ xlab(expression(paste(Salinity~(mg.L^-1))))+
    geom_point(aes(color = Aridity)) + geom_smooth(method = "lm", formula = 'y ~ x', se = F, linewidth = 0.7)+ stat_poly_eq(size = R.size, vstep = V.step, label.y = "top", label.x = "right") +ylim(c(0,1))+ annotate("text", x = min(Mfull.WDB$Salinity, na.rm = T), y = Annot.x, label = "(G)", hjust = 0, size = 5)
  
  P8 <- ggplot(Mfull, aes(y = MBTp6Me, color = Sample.type, shape = Sample.type, x = Salinity))+
    Log10.scale+ annotation_logticks(sides = 'b')+
    geom_hex(data = Mfull.WDB, mapping = aes(y = MBTp6Me, x = Salinity, fill = ..count..), bins = Bin.size, color = NA)+ Fill.scale+ xlab(expression(paste(Salinity~(mg.L^-1))))+
    geom_point(aes(color = Aridity)) + geom_smooth(method = "lm", formula = 'y ~ x', se = F, linewidth = 0.7)+ stat_poly_eq(size = R.size, vstep = V.step, label.y = "top", label.x = "right")+ylim(c(0,1))+ annotate("text", x = min(Mfull.WDB$Salinity, na.rm = T), y = Annot.x, label = "(H)", hjust = 0, size = 5)
  
  #### Scales ####
  PFULL <- ((P1 + P2) / (P3 + P4) / (P5 + P6) / (P7 + P8)) + plot_layout(guides = "collect") & theme_bw() & 
    theme(legend.position = "bottom", panel.grid = element_blank(), legend.title.position = "top",
          axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5, face = "bold")) & AI.class.scale & Shape.scale &
    guides(color = guide_legend(nrow = 3, override.aes = list(size = 3)), shape = guide_legend(nrow = 2, override.aes = list(size = 3)))
  
  #### Export ####
  H = 1200; W = 500; Save.plot ="Figures/Fig_6.pdf" 
  ggsave(filename = Save.plot, PFULL, width = W*0.026458333, height = H*0.026458333, units = "cm")
}

#### Fig. 7 (pH, CBT, IR) ####
if(Fig_7 == T){
  #### Plots ####
  P1 <- ggplot(Mfull, aes(y = CBTp, color = Sample.type, shape = Sample.type, x = pH))+
    geom_vline(xintercept = 7, lty = "dashed")+
    geom_hex(data = Mfull.WDB, mapping = aes(y = CBTp, x = pH, fill = ..count..), bins = Bin.size, color = NA)+ Fill.scale+ 
    geom_point(aes(color = Aridity)) + geom_smooth(method = "lm", formula = 'y ~ x', se = F, linewidth = 0.7) + stat_poly_eq(size = 3, vstep = 0.07, label.y = "bottom", label.x = "left") + annotate("text", x = Inf, y = Inf, label = "(A)", hjust = 1.1, vjust = 1.3, size = 5) + ylab("CBT'")
  
  P2 <- ggplot(Mfull, aes(y = CBT5Me, color = Sample.type, shape = Sample.type, x = pH))+
    geom_vline(xintercept = 7, lty = "dashed")+
    geom_hex(data = Mfull.WDB, mapping = aes(y = CBT5Me, x = pH, fill = ..count..), bins = Bin.size, color = NA)+ Fill.scale+ ylim(c(-2,3))+
    geom_point(aes(color = Aridity)) + geom_smooth(method = "lm", formula = 'y ~ x', se = F, linewidth = 0.7) + stat_poly_eq(size = 3, vstep = 0.07, label.y = "bottom") + annotate("text", x = Inf, y = Inf, label = "(B)", hjust = 1.1, vjust = 1.3, size = 5)+ylab(expression(paste(CBT,"'"[5~Me])))
  
  P3 <- ggplot(Mfull, aes(y = IR6Me, color = Sample.type, shape = Sample.type, x = pH))+
    geom_hex(data = Mfull.WDB, mapping = aes(y = IR6Me, x = pH, fill = ..count..), bins = Bin.size, color = NA)+ Fill.scale+
    geom_point(aes(color = Aridity)) + geom_smooth(method = "lm", formula = 'y ~ x', se = F, linewidth = 0.7) + stat_poly_eq(size = 3, vstep = 0.07, label.y = "top")+ylim(c(0,1))+ annotate("text", x = Inf, y = Inf, label = "(C)", hjust = 1.1, vjust = 1.3, size = 5)+ylab(expression(IR[6~Me]))
  
  P4 <- ggplot(Mfull, aes(y = IRp6_7Me, color = Sample.type, shape = Sample.type, x = pH))+
    geom_hex(data = Mfull.WDB, mapping = aes(y = IRp6_7Me, x = pH, fill = ..count..), bins = Bin.size, color = NA)+ Fill.scale+
    geom_point(aes(color = Aridity)) + geom_smooth(method = "lm", formula = 'y ~ x', se = F, linewidth = 0.7) + stat_poly_eq(size = 3, vstep = 0.07, label.y = "top")+ylim(c(0,1))+ annotate("text", x = Inf, y = Inf, label = "(D)", hjust = 1.1, vjust = 1.3, size = 5)+ylab(expression(paste(IR, "'"[6+7~Me])))
  
  P5 <- ggplot(Mfull, aes(y = CBTp, color = Sample.type, shape = Sample.type, x = AI)) +
    Log10.scale+ annotation_logticks(sides = 'b')+
    geom_hex(data = Mfull.WDB, mapping = aes(y = CBTp, x = AI, fill = ..count..), bins = Bin.size, color = NA)+ Fill.scale+
    geom_point(aes(color = Aridity)) + geom_smooth(method = "lm", formula = 'y ~ x', se = F, linewidth = 0.7) + stat_poly_eq(size = 3, vstep = 0.07, label.y = "top",  label.x = "left")+ylim(c(-2,3))+ annotate("text", x = Inf, y = Inf, label = "(E)", hjust = 1.1, vjust = 1.3, size = 5)+ ylab("CBT'")
  
  P6 <- ggplot(Mfull, aes(y = CBTp, color = Sample.type, shape = Sample.type, x = Salinity)) +
    Log10.scale+ annotation_logticks(sides = 'b')+
    geom_hex(data = Mfull.WDB, mapping = aes(y = CBTp, x = Salinity, fill = ..count..), bins = Bin.size, color = NA)+ Fill.scale+
    geom_point(aes(color = Aridity)) + geom_smooth(method = "lm", formula = 'y ~ x', se = F, linewidth = 0.7) + stat_poly_eq(size = 3, vstep = 0.07, label.y = "top")+ylim(c(-2,3))+ annotate("text", x = Inf, y = Inf, label = "(F)", hjust = 1.1, vjust = 1.3, size = 5)+ ylab("CBT'")
  
  P7 <- ggplot(Mfull, aes(y = CBTp, color = Sample.type, shape = Sample.type, x = IR6Me)) +
    geom_hex(data = Mfull.WDB, mapping = aes(y = CBTp, x = IR6Me, fill = ..count..), bins = Bin.size, color = NA)+ Fill.scale+
    geom_point(aes(color = Aridity)) + geom_smooth(method = "lm", formula = 'y ~ x', se = F, linewidth = 0.7)+ stat_poly_eq(size = 3, vstep = 0.07,  label.y = "bottom", label.x = "left")+ylim(c(-2,3))+ annotate("text", x = Inf, y = Inf, label = "(G)", hjust = 1.1, vjust = 1.3, size = 5) +xlab(expression(IR[6~Me]))+ ylab("CBT'")
  
  P8 <- ggplot(Mfull, aes(y = CBTp, color = Sample.type, shape = Sample.type, x = IRp6_7Me)) +
    geom_hex(data = Mfull.WDB, mapping = aes(y = CBTp, x = IRp6_7Me, fill = ..count..), bins = Bin.size, color = NA)+ Fill.scale+
    geom_point(aes(color = Aridity)) + geom_smooth(method = "lm", formula = 'y ~ x', se = F, linewidth = 0.7)+ stat_poly_eq(size = 3, vstep = 0.07, label.y = "bottom")+ylim(c(-2,3))+ annotate("text", x = Inf, y = Inf, label = "(H)", hjust = 1.1, vjust = 1.3, size = 5) +xlab(expression(paste(IR, "'"[6+7~Me]))) + ylab("CBT'")
  
  #### Scales ####
  PFULL <- (P1 | P2 | P3 | P4) / (P5 | P6 | P7 | P8) + plot_layout(guides = "collect") & theme_bw() & 
    theme(legend.position = "bottom", panel.grid = element_blank(), 
          title = element_text(size = 13)) & AI.class.scale & Shape.scale &
    guides(color = guide_legend(nrow = 1, override.aes = list(size = 3)), shape = guide_legend(nrow = 1, override.aes = list(size = 3)))
  
  #### Export ####
  H = 690*.8; W = 1300*.8; Save.plot ="Figures/Fig_7.pdf" 
  ggsave(filename = Save.plot, PFULL, width = W*0.026458333, height = H*0.026458333, units = "cm")
}

#### Fig. 8 (Salinity) ####
if(Fig_8 == T){
  #### Plots ####
  P1 <- ggplot(Mfull, aes(y = IR6Me, color = Sample.type, shape = Sample.type, x = Salinity))+
    Log10.scale+ annotation_logticks(sides = 'b')+
    geom_hex(data = Mfull.WDB, mapping = aes(y = IR6Me, x = Salinity, fill = ..count..), bins = Bin.size, color = NA)+ Fill.scale+ xlab(expression(paste(Salinity~(mg.L^-1))))+
    geom_point(aes(color = Aridity)) + geom_smooth(method = "lm", formula = 'y ~ x', se = F, linewidth = 0.7) + stat_poly_eq(size = 3, vstep = 0.07, label.y = "bottom", label.x = "right")+ylim(c(0,1))+ annotate("text", x = max(Mfull.WDB$Salinity, na.rm = T), y = Annot.x, label = "(A)", hjust = 1, size = 5)+ylab(expression(IR[6~Me]))
  
  P2 <- ggplot(Mfull, aes(y = IR7Me, color = Sample.type, shape = Sample.type, x = Salinity))+
    Log10.scale+ annotation_logticks(sides = 'b')+
    geom_hex(data = Mfull.WDB, mapping = aes(y = IR7Me, x = Salinity, fill = ..count..), bins = Bin.size, color = NA)+ Fill.scale+ xlab(expression(paste(Salinity~(mg.L^-1))))+
    geom_point(aes(color = Aridity)) + geom_smooth(method = "lm", formula = 'y ~ x', se = F, linewidth = 0.7) + stat_poly_eq(size = 3, vstep = 0.07)+ylim(c(0,1))+ annotate("text", x = max(Mfull.WDB$Salinity, na.rm = T), y = Annot.x, label = "(B)", hjust = 1, size = 5)+ylab(expression(IR[7~Me]))
  
  P3 <- ggplot(Mfull, aes(y = IR6_7Me, color = Sample.type, shape = Sample.type, x = Salinity))+
    Log10.scale+ annotation_logticks(sides = 'b')+
    geom_hex(data = Mfull.WDB, mapping = aes(y = IR6_7Me, x = Salinity, fill = ..count..), bins = Bin.size, color = NA)+ Fill.scale+ xlab(expression(paste(Salinity~(mg.L^-1))))+
    geom_point(aes(color = Aridity)) + geom_smooth(method = "lm", formula = 'y ~ x', se = F, linewidth = 0.7) + stat_poly_eq(size = 3, vstep = 0.07, label.y = "top")+ylim(c(0,1))+ annotate("text", x = max(Mfull.WDB$Salinity, na.rm = T), y = Annot.x, label = "(C)", hjust = 1, size = 5)+ylab(expression(IR[6+7~Me]))
  
  P4 <- ggplot(Mfull, aes(y = IRp6_7Me, color = Sample.type, shape = Sample.type, x = Salinity))+
    Log10.scale+ annotation_logticks(sides = 'b')+
    geom_hex(data = Mfull.WDB, mapping = aes(y = IRp6_7Me, x = Salinity, fill = ..count..), bins = Bin.size, color = NA)+ Fill.scale+ xlab(expression(paste(Salinity~(mg.L^-1))))+
    geom_point(aes(color = Aridity)) + geom_smooth(method = "lm", formula = 'y ~ x', se = F, linewidth = 0.7) + stat_poly_eq(size = 3, vstep = 0.07, label.y = "top")+ylim(c(0,1))+ annotate("text", x = max(Mfull.WDB$Salinity, na.rm = T), y = Annot.x, label = "(D)", hjust = 1, size = 5)+ylab(expression(paste(IR, "'"[6+7~Me])))
  
  P5 <- ggplot(Mfull, aes(y = IR6_7Me, color = Sample.type, shape = Sample.type, x = AI)) +
    Log10.scale+ annotation_logticks(sides = 'b')+
    geom_hex(data = Mfull.WDB, mapping = aes(y = IR6_7Me, x = AI, fill = ..count..), bins = Bin.size, color = NA)+ Fill.scale+
    geom_point(aes(color = Aridity)) + geom_smooth(method = "lm", formula = 'y ~ x', se = F, linewidth = 0.7) + stat_poly_eq(size = 3, vstep = 0.07, label.y = "top")+ylim(c(0,1))+ annotate("text", x = max(Mfull.WDB$AI, na.rm = T), y = Annot.x, label = "(E)", hjust = 1, size = 5)+ylab(expression(IR[6+7~Me]))
  
  P6 <- ggplot(Mfull, aes(y = IRp6_7Me, color = Sample.type, shape = Sample.type, x = AI)) +
    Log10.scale+ annotation_logticks(sides = 'b')+
    geom_hex(data = Mfull.WDB, mapping = aes(y = IRp6_7Me, x = AI, fill = ..count..), bins = Bin.size, color = NA)+ Fill.scale+
    geom_point(aes(color = Aridity)) + geom_smooth(method = "lm", formula = 'y ~ x', se = F, linewidth = 0.7) + stat_poly_eq(size = 3, vstep = 0.07, label.y = "top")+ylim(c(0,1))+ annotate("text", x = max(Mfull.WDB$AI, na.rm = T), y = Annot.x, label = "(F)", hjust = 1, size = 5)+ylab(expression(paste(IR, "'"[6+7~Me])))
  
  P7 <- ggplot(Mfull, aes(y = MBTp5Me, color = Sample.type, shape = Sample.type, x = IRp6_7Me)) +
    geom_hex(data = Mfull.WDB, mapping = aes(y = MBTp5Me, x = IRp6_7Me, fill = ..count..), bins = Bin.size, color = NA)+ Fill.scale+
    geom_point(aes(color = Aridity)) + geom_smooth(method = "lm", formula = 'y ~ x', se = F, linewidth = 0.7)+ stat_poly_eq(size = 3, vstep = 0.07,  label.y = "bottom", label.x = "right")+ylim(c(0,1))+ annotate("text", x = max(Mfull$IRp6_7Me, na.rm = T), y = Annot.x, label = "(G)", hjust = 1, size = 5)+xlab(expression(paste(IR, "'"[6+7~Me])))+ ylab(expression(paste(MBT,"'"[5~Me])))
  
  P8 <- ggplot(Mfull, aes(y = MBTp6Me, color = Sample.type, shape = Sample.type, x = IRp6_7Me)) +
    geom_hex(data = Mfull.WDB, mapping = aes(y = MBTp6Me, x = IRp6_7Me, fill = ..count..), bins = Bin.size, color = NA)+ Fill.scale+
    geom_point(aes(color = Aridity)) + geom_smooth(method = "lm", formula = 'y ~ x', se = F, linewidth = 0.7)+ stat_poly_eq(size = 3, vstep = 0.07, label.y = "bottom")+ylim(c(0,1))+ annotate("text", x = max(Mfull$IRp6_7Me, na.rm = T), y = Annot.x, label = "(H)", hjust = 1, size = 5)+xlab(expression(paste(IR, "'"[6+7~Me])))+ ylab(expression(paste(MBT,"'"[6~Me])))
  
  #### Scales ####
  PFULL <- (P1 | P2 | P3 | P4) / (P5 | P6 | P7 | P8) + plot_layout(guides = "collect") & theme_bw() & 
    theme(legend.position = "bottom", panel.grid = element_blank(), 
          title = element_text(size = 13)) & AI.class.scale & Shape.scale &
    guides(color = guide_legend(nrow = 1, override.aes = list(size = 3)), shape = guide_legend(nrow = 1, override.aes = list(size = 3)))
  
  #### Export ####
  H = 690*.8; W = 1300*.8; Save.plot ="Figures/Fig_8.pdf" 
  ggsave(filename = Save.plot, PFULL, width = W*0.026458333, height = H*0.026458333, units = "cm")
  
}

#### Table 3 (ANOVA, MANOVA) ####
if(Table_3 == T){
  #### P-values as stars ####
  p_value_to_stars <- function(p) {
    ifelse(p <= 0.001, "***", 
           ifelse(p <= 0.01, "** ", 
                  ifelse(p <= 0.05, "*  ", 
                         ifelse(p <= 0.1, "   ", "   "))))}
  p_value_to_test <- function(p) {
    ifelse(p <= 0.001, "**", 
           ifelse(p <= 0.01, "*", ""))}
  
  #### Cleaning data and test ####
  Mfull.mano <- Mfull
  Nb.evironmt <- 4
  Mfull.mano <- Mfull[complete.cases(Mfull[c("Salinity", "pH")]),]
  print(paste("Database size for MANOVA:", nrow(Mfull.mano)))
  List.indices <- c("MBTp5Me", "MBTp6Me", "IR6Me", "IR6_7Me", "IRp6_7Me", "CBTp", "CBT5Me")
  List.FA <- c("f.IIIa_5Me", "f.IIIa_6Me", "f.IIIb_5Me", "f.IIIb_6Me", "f.IIIc_5Me", "f.IIIc_6Me", "f.IIa_5Me", "f.IIa_6Me", "f.IIb_5Me", "f.IIb_6Me", "f.IIc_5Me", "f.IIc_6Me", "f.Ia", "f.Ib", "f.Ic")
  List.param <- c("\\textbf{F-statistics}", "f(IIIa)", "f(IIIa')", "f(IIIb)", "f(IIIb')", "f(IIIc)", "f(IIIc')", "f(IIa)", "f(IIa')", "f(IIb)", "f(IIb')", "f(IIc)", "f(IIc')", "f(Ia)", "f(Ib)", "f(Ic)", "\\textbf{F-statistics}", "$\\mathrm{MBT'_{5Me}}$", "$\\mathrm{MBT'_{6Me}}$", "$\\mathrm{IR_{6Me}}$", "$\\mathrm{IR_{6+7Me}}$", "$\\mathrm{IR'_{6+7Me}}$", "CBT'", "$\\mathrm{CBT'_{5Me}}$")
  Group.names <- c("Model", "Compound", "\\textbf{Alkalinity}", "\\textbf{Aridity}", "\\textbf{Salinity}", "\\textbf{Sample type}")
  
  #### 1. Assumption of Multivariate Normality (using Mardia's test) ####
  library(MVN)
  result <- mvn(data = Mfull.mano[, List.indices], mvnTest = "mardia")
  result <- mvn(data = Mfull.mano[, List.FA], mvnTest = "mardia")
  
  #### 2. Assumption of Homogeneity of Variance-Covariance (Equality of Covariance Matrices) using Levene test ####
  library(car)
  apply_levene_test <- function(df, group_var) {
    df[[group_var]] <- as.factor(df[[group_var]])
    
    results <- lapply(df[, sapply(df, is.numeric)], function(x) {
      levene_result <- leveneTest(x ~ df[[group_var]], data = df)
      return(levene_result$`Pr(>F)`[1])  # Extract p-value
    })
    
    result_df <- data.frame(
      Variable = names(results),
      P_Value = unlist(results)
    )
    return(result_df)
  }
  
  
  levene_results1 <- cbind(Acidity  = apply_levene_test(Mfull.mano[, c(List.FA, "Acidity")], "Acidity")[,1:2],
                           Aridity = apply_levene_test(Mfull.mano[, c(List.FA, "Aridity")], "Aridity")[,2],
                           Salinity_classe = apply_levene_test(Mfull.mano[, c(List.FA, "Salinity_classe")], "Salinity_classe")[,2],
                           Sample.type = apply_levene_test(Mfull.mano[, c(List.FA, "Sample.type")], "Sample.type")[,2]
  )
  levene_results2 <- cbind(Acidity  = apply_levene_test(Mfull.mano[, c(List.indices, "Acidity")], "Acidity")[,1:2],
                           Aridity = apply_levene_test(Mfull.mano[, c(List.indices, "Aridity")], "Aridity")[,2],
                           Salinity_classe = apply_levene_test(Mfull.mano[, c(List.indices, "Salinity_classe")], "Salinity_classe")[,2],
                           Sample.type = apply_levene_test(Mfull.mano[, c(List.indices, "Sample.type")], "Sample.type")[,2]
  )
  
  levene_pval1 <- sapply(levene_results1, p_value_to_test)
  levene_pval2 <- sapply(levene_results2, p_value_to_test)
  levene_pval1 <- rbind(levene_pval1, levene_pval2)
  row.names(levene_pval1) <- List.param[-c(1,17)]
  row.names(levene_pval1) <- List.param[-c(1,17)]
  levene_pval1 <- data.frame(levene_pval1[,-c(1)])
  names(levene_pval1) <- Group.names[3:length(Group.names)]
  
  Save.path.tex <- "Tables/Table_S5.tex"; library(xtable)
  LateX.caption <- "Levene's test results for the two MANOVA models used to evaluate the grouping factors effect accross brGFGT FA and indices. The violation of the assumption of homogeneity of variance-covariance is marked by ** (p-values < 0.001) or * (p-values < 0.01)"
  Tlatex <- xtable(levene_pval1, caption = LateX.caption, type = "latex", label = "SI_Table_Lavene")
  print(Tlatex, file = Save.path.tex, booktabs = T, include.rownames = T, comment = F,
        caption.placement = "top", sanitize.text.function = function(x){x},
        hline.after = c(-1,0,15,nrow(levene_pval1)))
  
  #### Model ####
  model1 <- manova(cbind(f.IIIa_5Me, f.IIIa_6Me, f.IIIb_5Me, f.IIIb_6Me, f.IIIc_5Me, f.IIIc_6Me, f.IIa_5Me, f.IIa_6Me, f.IIb_5Me, f.IIb_6Me, f.IIc_5Me, f.IIc_6Me, f.Ia, f.Ib, f.Ic) ~ Acidity*Aridity*Salinity_classe*Sample.type, data = Mfull.mano)
  MANOVA.res1 <- summary(model1)
  AOV.res1 <- summary.aov(model1)
  
  model2 <- manova(cbind(MBTp5Me, MBTp6Me, IR6Me, IR6_7Me, IRp6_7Me, CBTp, CBT5Me) ~ Acidity*Aridity*Salinity_classe*Sample.type, data = Mfull.mano)
  MANOVA.res2 <- summary(model2)
  AOV.res2 <- summary.aov(model2)
  
  #### Table results ####
  Tab.result <- setNames(data.frame(matrix(data = NA, nrow = 24, ncol = Nb.evironmt)), gsub(' ','', row.names(MANOVA.res1$stats)[1:Nb.evironmt]))
  Tab.result[1,] <- paste(round(MANOVA.res1$stats[1:Nb.evironmt,3], digits = 1), sapply(MANOVA.res1$stats[1:Nb.evironmt,6], p_value_to_stars), sep = "")
  Tab.result[2,] <- paste(round(AOV.res1$` Response f.IIIa_5Me`$`F value`[1:Nb.evironmt], digits = 1), sapply(AOV.res1$` Response f.IIIa_5Me`$`Pr(>F)`[1:Nb.evironmt], p_value_to_stars), sep = "")
  Tab.result[3,] <- paste(round(AOV.res1$` Response f.IIIa_6Me`$`F value`[1:Nb.evironmt], digits = 1), sapply(AOV.res1$` Response f.IIIa_6Me`$`Pr(>F)`[1:Nb.evironmt], p_value_to_stars), sep = "")
  Tab.result[4,] <- paste(round(AOV.res1$` Response f.IIIc_5Me`$`F value`[1:Nb.evironmt], digits = 1), sapply(AOV.res1$` Response f.IIIc_5Me`$`Pr(>F)`[1:Nb.evironmt], p_value_to_stars), sep = "")
  Tab.result[5,] <- paste(round(AOV.res1$` Response f.IIIc_6Me`$`F value`[1:Nb.evironmt], digits = 1), sapply(AOV.res1$` Response f.IIIc_6Me`$`Pr(>F)`[1:Nb.evironmt], p_value_to_stars), sep = "")
  Tab.result[6,] <- paste(round(AOV.res1$` Response f.IIIb_5Me`$`F value`[1:Nb.evironmt], digits = 1), sapply(AOV.res1$` Response f.IIIb_5Me`$`Pr(>F)`[1:Nb.evironmt], p_value_to_stars), sep = "")
  Tab.result[7,] <- paste(round(AOV.res1$` Response f.IIIb_6Me`$`F value`[1:Nb.evironmt], digits = 1), sapply(AOV.res1$` Response f.IIIb_6Me`$`Pr(>F)`[1:Nb.evironmt], p_value_to_stars), sep = "")
  Tab.result[8,] <- paste(round(AOV.res1$` Response f.IIa_5Me`$`F value`[1:Nb.evironmt], digits = 1), sapply(AOV.res1$` Response f.IIa_5Me`$`Pr(>F)`[1:Nb.evironmt], p_value_to_stars), sep = "")
  Tab.result[9,] <- paste(round(AOV.res1$` Response f.IIa_6Me`$`F value`[1:Nb.evironmt], digits = 1), sapply(AOV.res1$` Response f.IIa_6Me`$`Pr(>F)`[1:Nb.evironmt], p_value_to_stars), sep = "")
  Tab.result[10,] <- paste(round(AOV.res1$` Response f.IIb_5Me`$`F value`[1:Nb.evironmt], digits = 1), sapply(AOV.res1$` Response f.IIb_5Me`$`Pr(>F)`[1:Nb.evironmt], p_value_to_stars), sep = "")
  Tab.result[11,] <- paste(round(AOV.res1$` Response f.IIb_6Me`$`F value`[1:Nb.evironmt], digits = 1), sapply(AOV.res1$` Response f.IIb_6Me`$`Pr(>F)`[1:Nb.evironmt], p_value_to_stars), sep = "")
  Tab.result[12,] <- paste(round(AOV.res1$` Response f.IIc_5Me`$`F value`[1:Nb.evironmt], digits = 1), sapply(AOV.res1$` Response f.IIc_5Me`$`Pr(>F)`[1:Nb.evironmt], p_value_to_stars), sep = "")
  Tab.result[13,] <- paste(round(AOV.res1$` Response f.IIc_6Me`$`F value`[1:Nb.evironmt], digits = 1), sapply(AOV.res1$` Response f.IIc_6Me`$`Pr(>F)`[1:Nb.evironmt], p_value_to_stars), sep = "")
  Tab.result[14,] <- paste(round(AOV.res1$` Response f.Ia`$`F value`[1:Nb.evironmt], digits = 1), sapply(AOV.res1$` Response f.Ia`$`Pr(>F)`[1:Nb.evironmt], p_value_to_stars), sep = "")
  Tab.result[15,] <- paste(round(AOV.res1$` Response f.Ib`$`F value`[1:Nb.evironmt], digits = 1), sapply(AOV.res1$` Response f.Ib`$`Pr(>F)`[1:Nb.evironmt], p_value_to_stars), sep = "")
  Tab.result[16,] <- paste(round(AOV.res1$` Response f.Ic`$`F value`[1:Nb.evironmt], digits = 1), sapply(AOV.res1$` Response f.Ic`$`Pr(>F)`[1:Nb.evironmt], p_value_to_stars), sep = "")
  
  Tab.result[17,] <- paste(round(MANOVA.res2$stats[1:Nb.evironmt,3], digits = 1), sapply(MANOVA.res2$stats[1:Nb.evironmt,6], p_value_to_stars), sep = "")
  Tab.result[18,] <- paste(round(AOV.res2$` Response MBTp5Me`$`F value`[1:Nb.evironmt], digits = 1), sapply(AOV.res2$` Response MBTp5Me`$`Pr(>F)`[1:Nb.evironmt], p_value_to_stars), sep = "")
  Tab.result[19,] <- paste(round(AOV.res2$` Response MBTp6Me`$`F value`[1:Nb.evironmt], digits = 1), sapply(AOV.res2$` Response MBTp6Me`$`Pr(>F)`[1:Nb.evironmt], p_value_to_stars), sep = "")
  Tab.result[20,] <- paste(round(AOV.res2$` Response IR6Me`$`F value`[1:Nb.evironmt], digits = 1), sapply(AOV.res2$` Response IR6Me`$`Pr(>F)`[1:Nb.evironmt], p_value_to_stars), sep = "")
  Tab.result[21,] <- paste(round(AOV.res2$` Response IR6_7Me`$`F value`[1:Nb.evironmt], digits = 1), sapply(AOV.res2$` Response IR6_7Me`$`Pr(>F)`[1:Nb.evironmt], p_value_to_stars), sep = "")
  Tab.result[22,] <- paste(round(AOV.res2$` Response IRp6_7Me`$`F value`[1:Nb.evironmt], digits = 1), sapply(AOV.res2$` Response IRp6_7Me`$`Pr(>F)`[1:Nb.evironmt], p_value_to_stars), sep = "")
  
  Tab.result[23,] <- paste(round(AOV.res2$` Response CBTp`$`F value`[1:Nb.evironmt], digits = 1), sapply(AOV.res2$` Response CBTp`$`Pr(>F)`[1:Nb.evironmt], p_value_to_stars), sep = "")
  Tab.result[24,] <- paste(round(AOV.res2$` Response CBT5Me`$`F value`[1:Nb.evironmt], digits = 1), sapply(AOV.res2$` Response CBT5Me`$`Pr(>F)`[1:Nb.evironmt], p_value_to_stars), sep = "")
  
  #### Export table in LateX ####
  Tab.result <- cbind(c("\\textbf{Model 1 (FA)}", rep("", 15), "\\textbf{Model 2 (Indices)}", rep("", 7)), List.param, Tab.result)
  names(Tab.result) <- Group.names
  
  Save.path.tex <- "Tables/Table_3.tex"; library(xtable)
  LateX.caption <- "Statistical results of the two MANOVA carried to test (1) the brGDGT FA responses to MAAT, MAF, Salinity, AI and the sample types and (2) the $\\mathrm{MBT'_{5Me}}$, $\\mathrm{MBT'_{6Me}}$ and $\\mathrm{IR'_{6+7Me}}$ response to the same environmental parameters."
  Tlatex <- xtable(Tab.result, caption = LateX.caption, type = "latex", label = "Table_MANOVA")
  print(Tlatex, file = Save.path.tex, booktabs = T, include.rownames = F, comment = F,
        caption.placement = "top", sanitize.text.function = function(x){x},
        hline.after = c(-1,0,16,nrow(Tab.result)))
  
}

#### Fig. S.3 (MANOVA histogram) ####
if(Fig_S3 == T){
  #### Preparation ####
  List.indices <- c("MBTp5Me", "MBTp6Me", "IR6Me", "IR6_7Me", "IRp6_7Me")
  List.confounding <- c("Acidity", "Aridity", "Salinity_classe", "Sample.type")
  Mfull.hist <- Mfull[complete.cases(Mfull[c("Salinity", "pH")]),]
  Mfull.hist <- Mfull.hist[c(List.indices,List.confounding)]
  Mfull.hist <- reshape2::melt(Mfull.hist, id = List.confounding)
  Mfull.hist <- reshape2::melt(Mfull.hist, id = c("variable","value"))
  names(Mfull.hist) <- c("Indices", "value", "Grouping.factors", "Grouping.factors.val")
  
  List.indices2 <- c("CBTp", "CBT5Me")
  Mfull.hist2 <- Mfull[complete.cases(Mfull[c("Salinity", "pH")]),]
  Mfull.hist2 <- Mfull.hist2[c(List.indices2,List.confounding)]
  Mfull.hist2 <- reshape2::melt(Mfull.hist2, id = List.confounding)
  Mfull.hist2 <- reshape2::melt(Mfull.hist2, id = c("variable","value"))
  names(Mfull.hist2) <- c("Indices", "value", "Grouping.factors", "Grouping.factors.val")
  
  Clean.labels <- c("Acidity" = "Alkalinity",
                    "Aridity" = "Aridity",
                    "Salinity_classe" = "Salinity",
                    "Sample.type" = "Sample type")
  
  #### Plot ####
  phist <- ggplot(data = Mfull.hist, aes(x = Indices, y = value))+
    geom_boxplot(aes(fill = Grouping.factors.val), outliers = F)+
    facet_col(vars(Grouping.factors), labeller = labeller(Grouping.factors = Clean.labels))+
    scale_x_discrete(name = "brGDGT-based ratios",
                     labels = c(expression(paste(MBT,"'"[5~Me])),
                                expression(paste(MBT,"'"[6~Me])), 
                                expression(IR["6Me"]), 
                                expression(IR[6+7~Me]),
                                expression(paste(IR,"'"[6+7~Me]))))+
    scale_y_continuous(name = "Ratio values")
  
  phist2 <- ggplot(data = Mfull.hist2, aes(x = Indices, y = value))+
    geom_boxplot(aes(fill = Grouping.factors.val), outliers = F)+
    facet_col(vars(Grouping.factors), labeller = labeller(Grouping.factors = Clean.labels))+
    scale_x_discrete(name = "brGDGT-based ratios",
                     labels = c(expression(paste(CBT,"'")),
                                expression(paste(CBT,"'"[5~Me]))))+
    scale_y_continuous(name = NULL)
  
  phist <- phist + phist2 + plot_layout(guides = "collect", widths = c(5/7, 2/7))&
    theme(panel.background = element_rect(fill = NA, colour = 'black', linewidth = .5), panel.border = element_blank(), strip.background = element_blank(), strip.text = element_text(size = 13),
          panel.grid = element_blank(), axis.line.x = element_blank())
  
  #### Export ####
  H = 1000; W = 800; Save.plot ="Figures/Fig_S3.pdf" 
  ggsave(filename = Save.plot, phist, width = W*0.026458333, height = H*0.026458333, units = "cm")
  
}

#### Fig. 9 (Aridity delta MBT) ####
if(Fig_9 == T){
  #### Preparation ####
  library(ggforce)
  List.indices <- c("MBTp5Me", "MBTp6Me")
  List.confounding <- c("Aridity2", "Sample.type")
  Mfull$Aridity2 <- gsub("_", ". ", Mfull$Aridity)
  Mdelta <- Mfull[complete.cases(Mfull[c("Salinity", "pH")]),]
  Mdelta <- Mdelta[c(List.indices,List.confounding)]
  Mdelta <- reshape2::melt(Mdelta, id = List.confounding)
  Mdelta <- reshape2::melt(Mdelta, id = c("variable","value","Sample.type"))
  names(Mdelta) <- c("Indices", "value", "Sample.type", "Grouping.factors", "Grouping.factors.val")
  
  Clean.labels <- c("Acidity" = "Alkalinity",
                    "Aridity2" = "Aridity classes",
                    "Salinity_classe" = "Salinity",
                    "Sample.type" = "Sample type")
  
  #### Plot ####
  phist <- ggplot(data = Mdelta, aes(x = Indices, y = value))+
    geom_boxplot(aes(fill = Grouping.factors.val, colour = Sample.type), outliers = F)+
    scale_colour_manual(values = c("Soil" = "black", "Lacustrine" = "darkred"))+
    facet_col(vars(Grouping.factors), labeller = labeller(Grouping.factors = Clean.labels))+
    scale_x_discrete(name = "brGDGT-based ratios",
                     labels = c(expression(paste(MBT,"'"[5~Me])),
                                expression(paste(MBT,"'"[6~Me]))))+
    scale_y_continuous(name = "Ratio values")+
    scale_fill_manual(name = "Aridity classes",
                      values = c("1. Hyper-arid" = "#8c510a","2. Arid" = "#bf812e",
                                 "3. Semi-arid" = "#dfc27e","4. Dry sub-humid" = "#f5e9bf",
                                 "5. Humid" = "#80cec1"), labels = c("Hyper-arid", "Arid", "Semi-arid", "Dry sub-humid", "Humid"))+
    theme(axis.title = element_blank())
  
  Mfull$DeltaMBT <- Mfull$MBTp5Me - Mfull$MBTp6Me 
  P1 <- ggplot(Mfull, aes(y = MBTp5Me, x = AI, color = Sample.type))+ geom_point(alpha = .2)+ scale_y_continuous(name = expression(paste("MBT","'"[5~Me]))) +geom_smooth(method = "lm", formula = 'y ~ x', se = F)+ stat_poly_eq(size = R.size, vstep = V.step, label.y = "top", label.x = "right")+ scale_color_manual(values = c("Lacustrine" = "#75AADB","Soil" = "#EB9122"))
  P2 <- ggplot(Mfull, aes(y = MBTp6Me, x = AI, color = Sample.type))+ geom_point(alpha = .2)+ scale_y_continuous(name = expression(paste("MBT","'"[6~Me]))) +geom_smooth(method = "lm", formula = 'y ~ x', se = F)+ stat_poly_eq(size = R.size, vstep = V.step, label.y = "top", label.x = "left")+ scale_color_manual(values = c("Lacustrine" = "#75AADB","Soil" = "#EB9122"))
  P3 <- ggplot(Mfull, aes(y = DeltaMBT, x = AI, color = Sample.type))+geom_point(alpha = .2)+ scale_y_continuous(name = expression(paste(Delta,"(MBT","'"[5~Me], ", MBT","'"[6~Me], ")"))) +geom_smooth(method = "lm", formula = 'y ~ x', se = F)+ stat_poly_eq(size = R.size, vstep = V.step, label.y = "top", label.x = "right")+ scale_color_manual(values = c("Lacustrine" = "#75AADB","Soil" = "#EB9122"))
  
  layout <- "AAAA
               BBCC
               DDEE"
  phist <- phist + P1 + P2 + P3 + guide_area() + plot_annotation(tag_levels = 'A', tag_prefix = '(', tag_suffix = ")") +
    plot_layout(design = layout, guides = "collect")&
    theme(panel.background = element_rect(fill = NA, colour = 'black', linewidth = .5), panel.border = element_blank(), strip.background = element_blank(), strip.text = element_text(size = 13),
          panel.grid = element_blank(),
          axis.line.x = element_blank())
  
  #### Plot p2 ####
  Nclass <- 14
  Mfull$Temp.class <- cut_number(n = Nclass, Mfull$MAAT)
  Mfull <- Mfull[Mfull$Sample.type == "Soil",]
  levels(Mfull$Temp.class) <- paste("T", seq(1, Nclass), sep = "")
  A <- summary(lm(data = Mfull, AI ~ DeltaMBT*Temp.class))
  B <- summary(lm(data = Mfull, AI ~ DeltaMBT))
  C <- count(Mfull, Mfull$Temp.class)
  annotations <- data.frame(
    annotateText = paste("list(italic(R[multi.]^2) == ", round(A$r.squared, 2), 
                         ", n[k-MAAT] == ", Nclass, 
                         ", mu[n-k] == ", round(mean(C$n),0), 
                         ")"),
    xpos = c(-Inf), ypos =  c(-Inf), hjustvar = c(-0.05), vjustvar = c(-0.5))
  My_annot <- geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText), parse = T, size = 4)
  
  P4 <- ggplot(Mfull, aes(y = DeltaMBT, x = AI))+
    geom_point(aes(color = Temp.class), alpha = .2)+
    geom_smooth(aes(color = Temp.class), method = "lm", formula = 'y ~ x', se = F, alpha = .5)+
    geom_smooth(method = "lm", se = F, formula = 'y ~ x', color = "#EB9122", linewidth = 1.5)+
    My_annot +
    scale_y_continuous(name = expression(paste(Delta,"(MBT","'"[5~Me], ", MBT","'"[6~Me], ")")))+
    scale_color_manual(values = colorRampPalette(c("darkblue", "grey80", "darkred"))(Nclass))+
    theme(legend.position = "none",
          panel.grid = element_blank(),
          panel.background = element_rect(fill = NA, color = "black")) 
  
  P5 <- ggplot(Mfull, aes(y = MAAT, x = 1))+ 
    geom_boxplot(aes(color = Temp.class), alpha = .2, outliers = F)+ 
    scale_color_manual(values = colorRampPalette(c("darkblue", "grey80", "darkred"))(Nclass))+
    scale_x_continuous(name = "MAAT", position = "top")+
    theme(legend.position = "none", axis.ticks.x = element_blank(), axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill = NA, color = "black")) 
  
  P4 <- P4 + P5 + plot_layout(widths = c(.7,.3))
  
  phist <- (phist) / P4 + plot_layout(heights = c(.75,.25))
  
  #### Export ####
  H = 800; W = 500; Save.plot ="Figures/Fig_9.pdf" 
  ggsave(filename = Save.plot, phist, width = W*0.026458333, height = H*0.026458333, units = "cm")
}

#### Fig. 10 (Calibration for Salinity classes) ####
if(Fig_10 == T){
  LR.FA.clim(Msurf.mean["MBTp5Me"], Mclim["MAAT"], 
             Mtype = Meco, RegLig = "Local", 
             yLim = c(0,1), Ylab = expression(paste(MBT,"'"[5~Me])), Xlab = "MAAT (°C)",
             Facet = T, Bootstraps = F, Compare.all = T, Legend.position = "none", Stat.lm = T, xLim = c(-10,20), 
             Show.Plotly = F, Return.plot = T, T_test = T, Add.R2 = T, Select.type = c("Salinity_classe"),
             Save.csv = "Results/Coef_lr_MAAT_MBT5ME_salinity.csv",
             W = 1000*.9, H = 500*.9, Save.plot = "Figures/Fig_10.pdf")
}

#### Fig. S.4 (Calibration for pH classes) ####
if(Fig_S4 == T){
  LR.FA.clim(Msurf.mean["MBTp5Me"], Mclim["MAAT"], Return.plot = T,
             Mtype = Meco, RegLig = "Local", # Global / Local
             yLim = c(0,1), Facet = T, Compare.all = T, Stat.lm = T, T_test = T, Bootstraps = F,
             Ylab = expression(paste(MBT,"'"[5~Me])), Xlab = "MAAT (°C)", Legend.position = "none",
             Show.Plotly = F, Add.R2 = T, Select.type = "Acidity", xLim = c(-10,20),
             Save.csv = "Results/Coef_lr_MAAT_MBT5ME_acidity.csv",
             W = 1000*.9, H = 500*.9, Save.plot = "Figures/Fig_S4.pdf")}

#### Fig. S.5 (Calibration for Aridity classes) ####
if(Fig_S5 == T){
  Meco$Aridity <- gsub(".*_", "", Meco$Aridity)
  Meco$Aridity <- factor(Meco$Aridity, c('Hyper-arid', 'Arid', 'Semi-arid', 'Dry sub-humid', 'Humid'))
  LR.FA.clim(Msurf.mean["MBTp5Me"], Mclim["MAAT"], Return.plot = T,
             Mtype = Meco, RegLig = "Local", 
             yLim = c(0,1), Facet = T, Compare.all = T, Stat.lm = T, T_test = T, Bootstraps = F,
             Ylab = expression(paste(MBT,"'"[5~Me])), Xlab = "MAAT (°C)", Legend.position = "none",
             Show.Plotly = F, Add.R2 = T, Select.type = "Aridity", xLim = c(-10,20),
             Save.csv = "Results/Coef_lr_MAAT_MBT5ME_aridity.csv",
             W = 1500*.9, H = 600*.9, Save.plot = "Figures/Fig_S5.pdf")}

#### Fig. S.6 (Calibration for Sample type classes) ####
if(Fig_S6 == T){
  LR.FA.clim(Msurf.mean["MBTp5Me"], Mclim["MAAT"], Return.plot = T,
             Mtype = Meco, RegLig = "Local", 
             yLim = c(0,1), Facet = T, Compare.all = T, Stat.lm = T, T_test = T, Bootstraps = F,
             Ylab = expression(paste(MBT,"'"[5~Me])), Xlab = "MAAT (°C)", Legend.position = "none",
             Show.Plotly = F, Add.R2 = T, Select.type = "Sample.type", xLim = c(-10,20),
             Save.csv = "Results/Coef_lr_MAAT_MBT5ME_sample_type.csv",
             W = 1000*.9, H = 500*.9, Save.plot = "Figures/Fig_S6.pdf")}


#### Table S.1 (meta data of the ACADB new sites ) ####
if(Table_S1 == T){
  MSI <- Meco[Meco$Reference == "This study",c(1,4,3,8,19,21)]
  MSI <- cbind(Mclim[Meco$Reference == "This study",c(1,2,45)], MSI)
  MSI$Sample.type[MSI$Sample.type == "Moss"] <- "Soil"
  MSI$Altitude <- as.character(MSI$Altitude)
  MSI$Ecosystem <- gsub("Alau", "\\\\textit{Alau}", MSI$Ecosystem)
  MSI$Ecosystem <- gsub("Tau", "\\\\textit{Tau}", MSI$Ecosystem)
  MSI$Ecosystem <- gsub("Adyr", "\\\\textit{Adyr}", MSI$Ecosystem)
  MSI$Ecosystem <- gsub("Chol", "\\\\textit{Chol}", MSI$Ecosystem)
  MSI$Ecosystem <- gsub("Reaumuria soongorica", "\\\\textit{Reaumuria soongorica}", MSI$Ecosystem)
  MSI$Ecosystem <- gsub("Stipa breviflora, S.bungeana", "\\\\textit{Stipa breviflora, S.bungeana}", MSI$Ecosystem)
  MSI$Ecosystem <- gsub("Dasipbora fruticosa", "\\\\textit{Dasiphora fruticosa}", MSI$Ecosystem)
  MSI$Ecosystem <- gsub("Achnatherum splendens", "\\\\textit{Achnatherum splendens}", MSI$Ecosystem)
  MSI$Ecosystem <- gsub("Artemisia arenaria", "\\\\textit{Artemisia arenaria}", MSI$Ecosystem)
  MSI$Ecosystem <- gsub("Kobresia pygmaea", "\\\\textit{Kobresia pygmaea}", MSI$Ecosystem)
  MSI$Ecosystem <- gsub("Kobresia spp.", "\\\\textit{Kobresia} spp.", MSI$Ecosystem)
  MSI$Ecosystem <- gsub("Carex", "\\\\textit{Carex}", MSI$Ecosystem)
  MSI$Ecosystem <- gsub("Salix gilashanica", "\\\\textit{Salix gilashanica}", MSI$Ecosystem)
  MSI$Ecosystem <- gsub("Saussurea", "\\\\textit{Saussurea}", MSI$Ecosystem)
  MSI$Ecosystem <- gsub("Stipa purpurea", "\\\\textit{Stipa purpurea}", MSI$Ecosystem)
  MSI$Ecosystem <- gsub("Stipa breviflora", "\\\\textit{Stipa breviflora}", MSI$Ecosystem)
  
  Save.path.tex <- "Tables/Table_S1.tex"
  LateX.caption <- "Geographical and biological presentation of the new surface sites of the ACADB presented in this study."
  Tlatex <- xtable(MSI, caption = LateX.caption, type = "latex", label = "Table_SI_metadata")
  print(Tlatex, file = Save.path.tex, booktabs = T, include.rownames = T, comment = F,
        caption.placement = "top", sanitize.text.function = function(x){x})
}