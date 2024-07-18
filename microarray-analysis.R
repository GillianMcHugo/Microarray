# microarray analysis
## analysis of microarray data including qc, pca, de and functional enrichment analyses
### install required packages
if (!require(BiocManager)) install.packages("BiocManager")
if (!require(affy)) BiocManager::install("affy")
if (!require(arrayQualityMetrics)) BiocManager::install("arrayQualityMetrics")
if (!require(Biobase)) BiocManager::install("Biobase")
if (!require(devtools)) install.packages("devtools")
if (!require(farms)) devtools::install_github("bioinf-jku/farms")
if (!require(ComplexUpset)) install.packages("ComplexUpset")
if (!require(dplyr)) install.packages("dplyr")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(gprofiler2)) install.packages("gprofiler2")
if (!require(limma)) install.packages("limma")
if (!require(preprocessCore)) BiocManager::install("preprocessCore",
                                                   configure.args = c(preprocessCore = "--disable-threading"),
                                                   force = T,
                                                   update = T,
                                                   type = "source")
if (!require(stringr)) install.packages("stringr")

### list files
files <- list.files("01_data-sources/03_merged")

### extract sample information
samples <- stringr::str_split_fixed(files,
                                    "_",
                                    4) %>%
  as.data.frame() %>%
  dplyr::transmute("number" = dplyr::row_number(),
                   "population" = V1,
                   "id" = stringr::str_c(V1,
                                         V2,
                                         sep = "_"),
                   "tissue" = V3,
                   "dpi" = stringr::str_replace(V4,
                                                ".CEL",
                                                ""),
                   "set" = stringr::str_c(population,
                                          tissue,
                                          dpi,
                                          sep = "_"),
                   "name" = stringr::str_replace(files,
                                                 ".CEL",
                                                 ""),
                   "file" = files) %>%
  readr::write_csv("01_data-sources/samples.csv")

### read cel files
cel <- affy::ReadAffy(celfile.path = "01_data-sources/03_merged",
                      filenames = samples$file,
                      phenoData = Biobase::AnnotatedDataFrame(samples),
                      sampleNames = samples$name)

### function to draw and save expression boxplots
expression_boxplot <- function(input, filename){
  # function to generate colour palette
  smooth_rainbow <- khroma::colour("smooth rainbow")
  # generate colour palettes
  pal <- smooth_rainbow(31,
                        range = c(0.1,
                                  0.9))
  viridis <- viridis::viridis(41, direction = -1)
  # set colour palette
  cols <-  c(pal[28],
             pal[17],
             pal[9],
             pal[21])
  # set order of samples
  samples_order <- samples %>% 
    dplyr::arrange(tissue,
                   dpi) %>%
    dplyr::pull(name)
  # draw plot
  plot <- Biobase::exprs(input) %>%
    reshape2::melt(value.name = "intensity",
                   varnames = c("probe",
                                "sample")) %>%
    dplyr::mutate(population = stringr::str_split_i(sample,
                                                    "_",
                                                    1),
                  tissue = stringr::str_split_i(sample,
                                                "_",
                                                3),
                  dpi = stringr::str_split_i(sample,
                                             "_",
                                             4)) %>%
    ggplot2::ggplot(ggplot2::aes(x = interaction(factor(sample,
                                                        levels = samples_order),
                                                 factor(dpi)),
                                 y = intensity,
                                 fill = tissue)) +
    ggplot2::geom_boxplot(outlier.colour = "grey70") +
    ggh4x::facet_wrap2(~ factor(population,
                                levels = c("NDAM",
                                           "BORA")),
                       dir = "v",
                       scales = "free_x",
                       strip = ggh4x::strip_themed(background_y = ggh4x::elem_list_rect(fill = c(pal[5],
                                                                                                 pal[26])),
                                                   by_layer_y = F),
                       strip.position = "right") +
    ggplot2::scale_x_discrete(guide = ggh4x::guide_axis_nested()) +
    ggplot2::scale_y_continuous(expand = c(0,
                                           0),
                                limits = c(4,
                                           16)) +
    ggplot2::scale_fill_manual(name = "Tissue",
                               values = cols) +
    ggplot2::labs(x = "Days Post Infection",
                  y = "Log<sub>2</sub> Expression Intensity") +
    ggplot2::theme_light() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.title.y = ggtext::element_markdown(),
                   ggh4x.axis.nestline = ggplot2::element_line(colour = c(viridis[1],
                                                                          viridis[13],
                                                                          viridis[15],
                                                                          viridis[16],
                                                                          viridis[19],
                                                                          viridis[22],
                                                                          viridis[26],
                                                                          viridis[27],
                                                                          viridis[30],
                                                                          viridis[33],
                                                                          viridis[35],
                                                                          viridis[36]),
                                                               linewidth = 0.5),
                   ggh4x.axis.nesttext.x = ggplot2::element_text(),
                   legend.position = "top",
                   panel.border = ggplot2::element_rect(linewidth = ggplot2::rel(0.5)))
  # save plot
  png(paste0("02_quality-control/02_figures/",
             filename,
             ".png"),
      height = 6,
      res = 300,
      units = "in",
      width = 9)
  print(plot)
  dev.off()
  pdf(paste0("02_quality-control/02_figures/",
             filename,
             ".pdf"),
      height = 6,
      width = 9)
  print(plot)
  dev.off()
  # return plot
  plot
}

### apply expression boxplot function to raw data
raw_boxplot <- expression_boxplot(input = affy::expresso(cel,
                                                         bg.correct = F,
                                                         normalize = F,
                                                         pmcorrect.method = "pmonly",
                                                         summary.method = "farms",
                                                         summary.param = list(weight = 0.5,
                                                                              mu = 0,
                                                                              weighted.mean = F,
                                                                              robust = T,
                                                                              correction = 0,
                                                                              laplacian = F,
                                                                              centering = "median",
                                                                              spuriousCorrelation = 0)),
                                  filename = "raw-boxplot")

### perform quality control
arrayQualityMetrics::arrayQualityMetrics(do.logtransform = T,
                                         expressionset = cel,
                                         outdir = "02_quality-control/01_output")

### filter samples that failed 2 or more tests
samples_qc <- dplyr::filter(samples,
                            !name %in% c("BORA_090_LI_00",
                                         "BORA_100_LI_00",
                                         "BORA_171_SP_00",
                                         "BORA_174_LN_00",
                                         "NDAM_171_LI_26",
                                         "NDAM_171_LN_35",
                                         "NDAM_173_LN_21",
                                         "NDAM_178_LI_00",
                                         "NDAM_179_LI_26")) %>%
  readr::write_csv("02_quality-control/samples-qc.csv")

### read CEL files of filtered samples
cel_qc <- affy::ReadAffy(celfile.path = "01_data-sources/03_merged",
                         filenames = samples_qc$file,
                         phenoData = Biobase::AnnotatedDataFrame(samples_qc),
                         sampleNames = samples_qc$name)

### normalise the arrays
norm <- farms::qFarms(cel_qc)

### apply expression boxplot function to normalised data
normalised_boxplot <- expression_boxplot(input = norm,
                                         filename = "normalised-boxplot")

### combine boxplots
combined_boxplots <- raw_boxplot +
  normalised_boxplot +
  patchwork::plot_annotation(tag_levels = 'A') +
  patchwork::plot_layout(axes = "collect",
                         guides = "collect",
                         ncol = 1) &
  ggplot2::theme(legend.position = "top")

### save combined boxplots plot
png("02_quality-control/02_figures/combined-boxplots.png",
    height = 12,
    res = 300,
    units = "in",
    width = 9)
print(combined_boxplots)
dev.off()
pdf("02_quality-control/02_figures/combined-boxplots.pdf",
    height = 12,
    width = 9)
print(combined_boxplots)
dev.off()

### identify informative probes
ini <- farms::INIcalls(norm)

#### check percentage of informative probes
farms::summary(ini)

### extract expression set of informative probes
eset <- farms::getI_Eset(ini) 

### principal component analysis
#### function to perform and plot principal component analysis
microarray_pca <- function(pc2 = 2, pc1 = 1){
  # perform pca
  mds <- limma::plotMDS(eset,
                        dim.plot = c(pc1,
                                     pc2),
                        gene.selection = "common",
                        top = nrow(eset))
  # read and plot eigenvalues
  eigen <- as.data.frame(mds$eigen.values) %>%
    dplyr::rename(X1 = "mds$eigen.values") %>%
    dplyr::slice(1:10) %>%
    dplyr::mutate(prop = X1/sum(X1),
                  pc = row_number()) %>%
    ggplot2::ggplot(ggplot2::aes(pc,
                                 prop,
                                 fill = factor(pc))) +
    ggplot2::geom_bar(position = "dodge",
                      stat = "identity") +
    ggplot2::scale_x_continuous(breaks = c(1:10),
                                expand = c(0,
                                           0),
                                limits = c(0.55,
                                           10.45)) +
    ggplot2::scale_y_continuous(expand = c(0,
                                           0),
                                limits = c(0,
                                           0.85)) +
    ggplot2::scale_fill_manual(values = c(sapply(1:10,
                                                 function(i)
                                                   dplyr::case_when(pc1 == i ~ "#000000",
                                                                    pc2 == i ~ "#000000",
                                                                    .default = "#666666")))) +
    ggplot2::labs(x = "Principal Component",
                  y = "Proportion of Variance") +
    ggplot2::coord_equal(ratio = (10.45-0.55)/0.85) +
    ggplot2::theme_light() +
    ggplot2::theme(legend.position = "none")
  # set limits
  lims <- c("NDAM",
            "BORA")
  # function to generesponse colour palette
  smooth_rainbow <- khroma::colour("smooth rainbow")
  # generesponse colour palette 
  pal <- smooth_rainbow(31,
                        range = c(0.1,
                                  0.9))
  # set colour palette
  cols <-  c(pal[28],
             pal[17],
             pal[9],
             pal[21])
  # read and plot pca data
  pca <- samples_qc %>%
    dplyr::mutate(pca_x = mds$x,
                  pca_y = mds$y) %>%
    ggplot2::ggplot(ggplot2::aes(pca_x,
                                 pca_y,
                                 col = tissue,
                                 fill = tissue,
                                 pch = population)) +
    ggplot2::geom_point(size = 2.5,
                        stroke = 0.5) +
    ggplot2::scale_colour_manual(name = "Tissue",
                                 values = cols) +
    ggplot2::scale_fill_manual(name = "Tissue",
                               values = ggplot2::alpha(cols,
                                                       0.5)) +
    ggplot2::scale_shape_manual(limits = lims,
                                name = "Population",
                                values = c(21,
                                           25)) +
    ggplot2::labs(x = paste0("Principal Component ",
                             pc1,
                             " (",
                             stringr::str_pad(round(eigen$data$prop[pc1]*100,
                                                    2),
                                              4,
                                              pad = "0",
                                              side = "right"),
                             "%)"),
                  y = paste0("Principal Component ",
                             pc2,
                             " (",
                             stringr::str_pad(round(eigen$data$prop[pc2]*100,
                                                    2),
                                              4,
                                              pad = "0",
                                              side = "right"),
                             "%)")) +
    ggplot2::theme_light()
  # set layout
  layout <- "AAB
             AAC"
  # combine plots
  plot <- pca +
    patchwork::guide_area() +
    eigen +
    patchwork::plot_annotation(tag_levels = 'A') +
    patchwork::plot_layout(design = layout,
                           guides = "collect")
  # save plot
  png(paste0("03_principal-component-analysis/pc-",
             stringr::str_pad(pc1,
                              2,
                              side = "left",
                              pad = "0"),
             "-",
             stringr::str_pad(pc2,
                              2,
                              side = "left",
                              pad = "0"),
             ".png"),
      height = 6,
      res = 300,
      units = "in",
      width = 9)
  print(plot)
  dev.off()
  pdf(paste0("03_principal-component-analysis/pc-",
             stringr::str_pad(pc1,
                              2,
                              side = "left",
                              pad = "0"),
             "-",
             stringr::str_pad(pc2,
                              2,
                              side = "left",
                              pad = "0"),
             ".pdf"),
      height = 6,
      width = 9)
  print(plot)
  dev.off()
}

#### apply function
microarray_pca()
microarray_pca(3)
microarray_pca(4)
microarray_pca(5)

#### function to perform and plot principal component analysis coloured by days post infection
microarray_pca_dpi <- function(pc2 = 2, pc1 = 1){
  # perform pca
  mds <- limma::plotMDS(eset,
                        dim.plot = c(pc1,
                                     pc2),
                        gene.selection = "common",
                        top = nrow(eset))
  # read and plot eigenvalues
  eigen <- as.data.frame(mds$eigen.values) %>%
    dplyr::rename(X1 = "mds$eigen.values") %>%
    dplyr::slice(1:10) %>%
    dplyr::mutate(prop = X1/sum(X1),
                  pc = row_number()) %>%
    ggplot2::ggplot(ggplot2::aes(pc,
                                 prop,
                                 fill = factor(pc))) +
    ggplot2::geom_bar(position = "dodge",
                      stat="identity") +
    ggplot2::scale_x_continuous(breaks = c(1:10),
                                expand = c(0,
                                           0),
                                limits = c(0.55,
                                           10.45)) +
    ggplot2::scale_y_continuous(expand = c(0,
                                           0),
                                limits = c(0,
                                           0.85)) +
    ggplot2::scale_fill_manual(values = c(sapply(1:10,
                                                 function(i)
                                                   dplyr::case_when(pc1 == i ~ "#000000",
                                                                    pc2 == i ~ "#000000",
                                                                    .default = "#666666")))) +
    ggplot2::labs(x = "Principal Component",
                  y = "Proportion of Variance") +
    ggplot2::coord_equal(ratio = (10.45-0.55)/0.85) +
    ggplot2::theme_light() +
    ggplot2::theme(legend.position = "none")
  # set limits
  lims <- c("NDAM",
            "BORA")
  # function to generesponse colour palette
  smooth_rainbow <- khroma::colour("smooth rainbow")
  # generesponse colour palette 
  pal <- smooth_rainbow(31,
                        range = c(0.1,
                                  0.9))
  # set colour palette
  cols <-  c(pal[28],
             pal[17],
             pal[9],
             pal[21])
  # read and plot pca data
  pca <- samples_qc %>%
    dplyr::mutate(dpi = as.numeric(dpi),
                  pca_x = mds$x,
                  pca_y = mds$y) %>%
    ggplot2::ggplot(ggplot2::aes(pca_x,
                                 pca_y,
                                 col = tissue,
                                 fill = dpi,
                                 pch = population)) +
    ggplot2::geom_point(size = 2.5,
                        stroke = 0.5) +
    ggplot2::scale_colour_manual(name = "Tissue",
                                 values = cols) +
    ggplot2::scale_fill_viridis_c(alpha = 0.5,
                                  begin = 1-((1/41)*36),
                                  direction = -1,
                                  name = "DPI") +
    ggplot2::scale_shape_manual(limits = lims,
                                name = "Population",
                                values = c(21,
                                           25)) +
    ggplot2::labs(x = paste0("Principal Component ",
                             pc1,
                             " (",
                             stringr::str_pad(round(eigen$data$prop[pc1]*100,
                                                    2),
                                              4,
                                              pad = "0",
                                              side = "right"),
                             "%)"),
                  y = paste0("Principal Component ",
                             pc2,
                             " (",
                             stringr::str_pad(round(eigen$data$prop[pc2]*100,
                                                    2),
                                              4,
                                              pad = "0",
                                              side = "right"),
                             "%)")) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(shape = 21)),
                    fill = ggplot2::guide_colourbar(draw.ulim = F,
                                                    draw.llim = F,
                                                    order = 1)) +
    ggplot2::theme_light()
  # set layout
  layout <- "AAB
             AAC"
  # combine plots
  plot <- pca +
    patchwork::guide_area() +
    eigen +
    patchwork::plot_annotation(tag_levels = 'A') +
    patchwork::plot_layout(design = layout,
                           guides = "collect") &
    ggplot2::theme(legend.box = "horizontal")
  # save plot
  png(paste0("03_principal-component-analysis/pc-",
             stringr::str_pad(pc1,
                              2,
                              side = "left",
                              pad = "0"),
             "-",
             stringr::str_pad(pc2,
                              2,
                              side = "left",
                              pad = "0"),
             "-dpi.png"),
      height = 6,
      res = 300,
      units = "in",
      width = 9)
  print(plot)
  dev.off()
  pdf(paste0("03_principal-component-analysis/pc-",
             stringr::str_pad(pc1,
                              2,
                              side = "left",
                              pad = "0"),
             "-",
             stringr::str_pad(pc2,
                              2,
                              side = "left",
                              pad = "0"),
             "-dpi.pdf"),
      height = 6,
      width = 9)
  print(plot)
  dev.off()
}

#### apply function
microarray_pca_dpi()
microarray_pca_dpi(3)
microarray_pca_dpi(4)
microarray_pca_dpi(5)

### set design matrix
design <- model.matrix(~ 0 + factor(samples_qc$set,
                                    levels = unique(samples_qc$set)))

### rename the columns of the design matrix
colnames(design) <- unique(samples_qc$set)

### identify correlation between individuals with multiple samples in the normalised data
correlation <- limma::duplicateCorrelation(norm,
                                           design = design,
                                           ndups = 1,
                                           block = samples_qc$id)

### fit a linear model to the normalised data
fit <- limma::lmFit(norm,
                    design = design,
                    ndups = 1, 
                    cor = correlation$consensus,
                    block = samples_qc$id)

### make contrasts
contrasts <- limma::makeContrasts("RESP_BL_14" = "(NDAM_BL_14-NDAM_BL_00)-(BORA_BL_14-BORA_BL_00)",
                                  "RESP_BL_25" = "(NDAM_BL_25-NDAM_BL_00)-(BORA_BL_25-BORA_BL_00)",
                                  "RESP_BL_34" = "(NDAM_BL_34-NDAM_BL_00)-(BORA_BL_34-BORA_BL_00)",
                                  "RESP_LI_12" = "(NDAM_LI_12-NDAM_LI_00)-(BORA_LI_12-BORA_LI_00)",
                                  "RESP_LI_15" = "(NDAM_LI_15-NDAM_LI_00)-(BORA_LI_15-BORA_LI_00)",
                                  "RESP_LI_18" = "(NDAM_LI_18-NDAM_LI_00)-(BORA_LI_18-BORA_LI_00)",
                                  "RESP_LI_21" = "(NDAM_LI_21-NDAM_LI_00)-(BORA_LI_21-BORA_LI_00)",
                                  "RESP_LI_26" = "(NDAM_LI_26-NDAM_LI_00)-(BORA_LI_26-BORA_LI_00)",
                                  "RESP_LI_29" = "(NDAM_LI_29-NDAM_LI_00)-(BORA_LI_29-BORA_LI_00)",
                                  "RESP_LI_32" = "(NDAM_LI_32-NDAM_LI_00)-(BORA_LI_32-BORA_LI_00)",
                                  "RESP_LI_35" = "(NDAM_LI_35-NDAM_LI_00)-(BORA_LI_35-BORA_LI_00)",
                                  "RESP_LN_21" = "(NDAM_LN_21-NDAM_LN_00)-(BORA_LN_21-BORA_LN_00)",
                                  "RESP_LN_35" = "(NDAM_LN_35-NDAM_LN_00)-(BORA_LN_35-BORA_LN_00)",
                                  "RESP_SP_21" = "(NDAM_SP_21-NDAM_SP_00)-(BORA_SP_21-BORA_SP_00)",
                                  "RESP_SP_35" = "(NDAM_SP_35-NDAM_SP_00)-(BORA_SP_35-BORA_SP_00)",
                                  "DIRE_BL_00" = "NDAM_BL_00-BORA_BL_00",
                                  "DIRE_BL_14" = "NDAM_BL_14-BORA_BL_14",
                                  "DIRE_BL_25" = "NDAM_BL_25-BORA_BL_25",
                                  "DIRE_BL_34" = "NDAM_BL_34-BORA_BL_34",
                                  "DIRE_LI_00" = "NDAM_LI_00-BORA_LI_00",
                                  "DIRE_LI_12" = "NDAM_LI_12-BORA_LI_12",
                                  "DIRE_LI_15" = "NDAM_LI_15-BORA_LI_15",
                                  "DIRE_LI_18" = "NDAM_LI_18-BORA_LI_18",
                                  "DIRE_LI_21" = "NDAM_LI_21-BORA_LI_21",
                                  "DIRE_LI_26" = "NDAM_LI_26-BORA_LI_26",
                                  "DIRE_LI_29" = "NDAM_LI_29-BORA_LI_29",
                                  "DIRE_LI_32" = "NDAM_LI_32-BORA_LI_32",
                                  "DIRE_LI_35" = "NDAM_LI_35-BORA_LI_35",
                                  "DIRE_LN_00" = "NDAM_LN_00-BORA_LN_00",
                                  "DIRE_LN_21" = "NDAM_LN_21-BORA_LN_21",
                                  "DIRE_LN_35" = "NDAM_LN_35-BORA_LN_35",
                                  "DIRE_SP_00" = "NDAM_SP_00-BORA_SP_00",
                                  "DIRE_SP_21" = "NDAM_SP_21-BORA_SP_21",
                                  "DIRE_SP_35" = "NDAM_SP_35-BORA_SP_35",
                                  "NDAM_BL_14" = "NDAM_BL_14-NDAM_BL_00",
                                  "NDAM_BL_25" = "NDAM_BL_25-NDAM_BL_00",
                                  "NDAM_BL_34" = "NDAM_BL_34-NDAM_BL_00",
                                  "NDAM_LI_12" = "NDAM_LI_12-NDAM_LI_00",
                                  "NDAM_LI_15" = "NDAM_LI_15-NDAM_LI_00",
                                  "NDAM_LI_18" = "NDAM_LI_18-NDAM_LI_00",
                                  "NDAM_LI_21" = "NDAM_LI_21-NDAM_LI_00",
                                  "NDAM_LI_26" = "NDAM_LI_26-NDAM_LI_00",
                                  "NDAM_LI_29" = "NDAM_LI_29-NDAM_LI_00",
                                  "NDAM_LI_32" = "NDAM_LI_32-NDAM_LI_00",
                                  "NDAM_LI_35" = "NDAM_LI_35-NDAM_LI_00",
                                  "NDAM_LN_21" = "NDAM_LN_21-NDAM_LN_00",
                                  "NDAM_LN_35" = "NDAM_LN_35-NDAM_LN_00",
                                  "NDAM_SP_21" = "NDAM_SP_21-NDAM_SP_00",
                                  "NDAM_SP_35" = "NDAM_SP_35-NDAM_SP_00",
                                  "BORA_BL_14" = "BORA_BL_14-BORA_BL_00",
                                  "BORA_BL_25" = "BORA_BL_25-BORA_BL_00",
                                  "BORA_BL_34" = "BORA_BL_34-BORA_BL_00",
                                  "BORA_LI_12" = "BORA_LI_12-BORA_LI_00",
                                  "BORA_LI_15" = "BORA_LI_15-BORA_LI_00",
                                  "BORA_LI_18" = "BORA_LI_18-BORA_LI_00",
                                  "BORA_LI_21" = "BORA_LI_21-BORA_LI_00",
                                  "BORA_LI_26" = "BORA_LI_26-BORA_LI_00",
                                  "BORA_LI_29" = "BORA_LI_29-BORA_LI_00",
                                  "BORA_LI_32" = "BORA_LI_32-BORA_LI_00",
                                  "BORA_LI_35" = "BORA_LI_35-BORA_LI_00",
                                  "BORA_LN_21" = "BORA_LN_21-BORA_LN_00",
                                  "BORA_LN_35" = "BORA_LN_35-BORA_LN_00",
                                  "BORA_SP_21" = "BORA_SP_21-BORA_SP_00",
                                  "BORA_SP_35" = "BORA_SP_35-BORA_SP_00",
                                  levels = design)

### set order of contrasts
contrast_ids <- c("DIRE_BL_00",
                  "DIRE_LI_00",
                  "DIRE_LN_00",
                  "DIRE_SP_00",
                  "RESP_LI_12",
                  "DIRE_LI_12",
                  "NDAM_LI_12",
                  "BORA_LI_12",
                  "RESP_BL_14",
                  "DIRE_BL_14",
                  "NDAM_BL_14",
                  "BORA_BL_14",
                  "RESP_LI_15",
                  "DIRE_LI_15",
                  "NDAM_LI_15",
                  "BORA_LI_15",
                  "RESP_LI_18",
                  "DIRE_LI_18",
                  "NDAM_LI_18",
                  "BORA_LI_18",
                  "RESP_LI_21",
                  "DIRE_LI_21",
                  "NDAM_LI_21",
                  "BORA_LI_21",
                  "RESP_LN_21",
                  "DIRE_LN_21",
                  "NDAM_LN_21",
                  "BORA_LN_21",
                  "RESP_SP_21",
                  "DIRE_SP_21",
                  "NDAM_SP_21",
                  "BORA_SP_21",
                  "RESP_BL_25",
                  "DIRE_BL_25",
                  "NDAM_BL_25",
                  "BORA_BL_25",
                  "RESP_LI_26",
                  "DIRE_LI_26",
                  "NDAM_LI_26",
                  "BORA_LI_26",
                  "RESP_LI_29",
                  "DIRE_LI_29",
                  "NDAM_LI_29",
                  "BORA_LI_29",
                  "RESP_LI_32",
                  "DIRE_LI_32",
                  "NDAM_LI_32",
                  "BORA_LI_32",
                  "RESP_BL_34",
                  "DIRE_BL_34",
                  "NDAM_BL_34",
                  "BORA_BL_34",
                  "RESP_LI_35",
                  "DIRE_LI_35",
                  "NDAM_LI_35",
                  "BORA_LI_35",
                  "RESP_LN_35",
                  "DIRE_LN_35",
                  "NDAM_LN_35",
                  "BORA_LN_35",
                  "RESP_SP_35",
                  "DIRE_SP_35",
                  "NDAM_SP_35",
                  "BORA_SP_35")

### fit contrasts Apply contrasts to the linear model farms_fit
contrasts_fit <- limma::contrasts.fit(fit,
                                      contrasts)

### filter the fitted contrasts using the informative probes expression set
contrasts_fit_informative <- contrasts_fit[rownames(contrasts_fit) %in% rownames(exprs(eset)),] 

### calculate empirical bayes statistics
ebayes <- limma::eBayes(contrasts_fit_informative)

### apply global multiple correction
results <- limma::decideTests(ebayes,
                              method = "global")

### save results
limma::write.fit(ebayes,
                 results = results,
                 file = "04_differential-expression-analysis/01_output/de-results.csv",
                 adjust = "BH",
                 sep = ",",
                 method = "global",
                 F.adjust = "BH")

### read results
de_results <- readr::read_csv("04_differential-expression-analysis/01_output/de-results.csv") %>%
  dplyr::rename(probe = "...1")

### convert probes to genes
genes <- gprofiler2::gconvert(query = de_results$probe,
                              organism = "btaurus") %>%
  dplyr::rename(probe = input,
                ensembl = target,
                symbol = name)

### combine de results with genes
de_genes <- dplyr::full_join(de_results,
                             genes) %>%
  readr::write_csv("04_differential-expression-analysis/01_output/de-genes.csv")

### extract all de genes
all_de_genes <- data.frame(set = NA, symbol = NA)
for (i in 1:length(contrast_ids)){
  all_increased_genes <- de_genes %>%
    tidyr::drop_na(symbol) %>%
    dplyr::filter(!!rlang::sym(paste0("Results.",
                                      contrast_ids[i])) ==  1) %>%
    dplyr::distinct(symbol,
                    .keep_all = T) %>%
    dplyr::select(symbol) %>%
    dplyr::mutate(set = paste0(contrast_ids[i],
                               "_increased"))
  all_decreased_genes <- de_genes %>%
    tidyr::drop_na(symbol) %>%
    dplyr::filter(!!rlang::sym(paste0("Results.",
                                      contrast_ids[i])) ==  -1) %>%
    dplyr::distinct(symbol,
                    .keep_all = T) %>%
    dplyr::select(symbol) %>%
    dplyr::mutate(set = paste0(contrast_ids[i],
                               "_decreased"))
  all_de_genes <- dplyr::full_join(all_de_genes,
                                   all_increased_genes)
  all_de_genes <- dplyr::full_join(all_de_genes,
                                   all_decreased_genes)
}

### extract all de gene symbols
all_de_gene_symbols <- data.frame(set = NA, symbol = NA)
for (i in 1:length(contrast_ids)){
  all_increased_gene_symbols <- de_genes %>%
    tidyr::drop_na(symbol) %>%
    dplyr::filter(!stringr::str_detect(symbol, "^ENSBTAG")) %>%
    dplyr::filter(!!rlang::sym(paste0("Results.",
                                      contrast_ids[i])) ==  1) %>%
    dplyr::distinct(symbol,
                    .keep_all = T) %>%
    dplyr::select(symbol) %>%
    dplyr::mutate(set = paste0(contrast_ids[i],
                               "_increased"))
  all_decreased_gene_symbols <- de_genes %>%
    tidyr::drop_na(symbol) %>%
    dplyr::filter(!stringr::str_detect(symbol, "^ENSBTAG")) %>%
    dplyr::filter(!!rlang::sym(paste0("Results.",
                                      contrast_ids[i])) ==  -1) %>%
    dplyr::distinct(symbol,
                    .keep_all = T) %>%
    dplyr::select(symbol) %>%
    dplyr::mutate(set = paste0(contrast_ids[i],
                               "_decreased"))
  all_de_gene_symbols <- dplyr::full_join(all_de_gene_symbols,
                                          all_increased_gene_symbols)
  all_de_gene_symbols <- dplyr::full_join(all_de_gene_symbols,
                                          all_decreased_gene_symbols)
}

### count number of duplicated genes and gene symbols
all_de_gene_counts <- all_de_genes %>%
  dplyr::select(symbol) %>%
  tidyr::drop_na() %>%
  group_by(symbol) %>%
  summarize(n_contrasts_de = n()) %>% 
  readr::write_csv("04_differential-expression-analysis/01_output/all-de-gene-counts.csv")

all_de_gene_symbol_counts <- all_de_gene_symbols %>%
  dplyr::select(symbol) %>%
  tidyr::drop_na() %>%
  group_by(symbol) %>%
  summarize(n_contrasts_de = n()) %>% 
  readr::write_csv("04_differential-expression-analysis/01_output/all-de-gene-symbol-counts.csv")

### extract top 10 de genes
top_de_genes <- data.frame(set = NA, symbol = NA)
for (i in 1:length(contrast_ids)){
  top_increased_genes <- de_genes %>%
    tidyr::drop_na(symbol) %>%
    dplyr::filter(!!rlang::sym(paste0("Results.",
                                      contrast_ids[i])) ==  1) %>%
    dplyr::arrange(!!rlang::sym(paste0("P.value.adj.",
                                       contrast_ids[i]))) %>%
    dplyr::distinct(symbol,
                    .keep_all = T) %>%
    dplyr::slice_min(n = 10,
                     order_by = !!rlang::sym(paste0("P.value.adj.",
                                                    contrast_ids[i]))) %>%
    dplyr::select(symbol) %>%
    dplyr::mutate(set = paste0(contrast_ids[i],
                               "_increased"))
  top_decreased_genes <- de_genes %>%
    tidyr::drop_na(symbol) %>%
    dplyr::filter(!!rlang::sym(paste0("Results.",
                                      contrast_ids[i])) ==  -1) %>%
    dplyr::distinct(symbol,
                    .keep_all = T) %>%
    dplyr::slice_min(n = 10,
                     order_by = !!rlang::sym(paste0("P.value.adj.",
                                                    contrast_ids[i]))) %>%
    dplyr::select(symbol) %>%
    dplyr::mutate(set = paste0(contrast_ids[i],
                               "_decreased"))
  top_de_genes <- dplyr::full_join(top_de_genes,
                                   top_increased_genes)
  top_de_genes <- dplyr::full_join(top_de_genes,
                                   top_decreased_genes)
}
top_de_genes %>% 
  dplyr::group_by(set) %>%
  dplyr::summarise(genes = toString(symbol)) %>%
  tidyr::drop_na() %>%
  dplyr::mutate(contrast_type = stringr::str_split_i(set,
                                                     "_",
                                                     1),
                contrast_type = factor(contrast_type,
                                       levels = c("RESP",
                                                  "DIRE",
                                                  "NDAM",
                                                  "BORA")),
                expression = stringr::str_split_i(set,
                                                  "_",
                                                  4),
                tissue = stringr::str_split_i(set,
                                              "_",
                                              2),
                dpi = stringr::str_split_i(set,
                                           "_",
                                           3)) %>%
  dplyr::arrange(contrast_type, desc(expression), tissue, dpi) %>%
  dplyr::relocate(genes,
                  .after = last_col()) %>% 
  readr::write_csv("04_differential-expression-analysis/01_output/top-de-genes.csv")

### extract top 10 de gene symbols
top_de_gene_symbols <- data.frame(set = NA, symbol = NA)
for (i in 1:length(contrast_ids)){
  top_increased_gene_symbols <- de_genes %>%
    tidyr::drop_na(symbol) %>%
    dplyr::filter(!stringr::str_detect(symbol, "^ENSBTAG")) %>%
    dplyr::filter(!!rlang::sym(paste0("Results.",
                                      contrast_ids[i])) ==  1) %>%
    dplyr::arrange(!!rlang::sym(paste0("P.value.adj.",
                                       contrast_ids[i]))) %>%
    dplyr::distinct(symbol,
                    .keep_all = T) %>%
    dplyr::slice_min(n = 10,
                     order_by = !!rlang::sym(paste0("P.value.adj.",
                                                    contrast_ids[i]))) %>%
    dplyr::select(symbol) %>%
    dplyr::mutate(set = paste0(contrast_ids[i],
                               "_increased"))
  top_decreased_gene_symbols <- de_genes %>%
    tidyr::drop_na(symbol) %>%
    dplyr::filter(!stringr::str_detect(symbol, "^ENSBTAG")) %>%
    dplyr::filter(!!rlang::sym(paste0("Results.",
                                      contrast_ids[i])) ==  -1) %>%
    dplyr::distinct(symbol,
                    .keep_all = T) %>%
    dplyr::slice_min(n = 10,
                     order_by = !!rlang::sym(paste0("P.value.adj.",
                                                    contrast_ids[i]))) %>%
    dplyr::select(symbol) %>%
    dplyr::mutate(set = paste0(contrast_ids[i],
                               "_decreased"))
  top_de_gene_symbols <- dplyr::full_join(top_de_gene_symbols,
                                          top_increased_gene_symbols)
  top_de_gene_symbols <- dplyr::full_join(top_de_gene_symbols,
                                          top_decreased_gene_symbols)
}
top_de_gene_symbols %>% 
  dplyr::group_by(set) %>%
  dplyr::summarise(genes = toString(symbol)) %>%
  tidyr::drop_na() %>%
  dplyr::mutate(contrast_type = stringr::str_split_i(set,
                                                     "_",
                                                     1),
                contrast_type = factor(contrast_type,
                                       levels = c("RESP",
                                                  "DIRE",
                                                  "NDAM",
                                                  "BORA")),
                expression = stringr::str_split_i(set,
                                                  "_",
                                                  4),
                tissue = stringr::str_split_i(set,
                                              "_",
                                              2),
                dpi = stringr::str_split_i(set,
                                           "_",
                                           3)) %>%
  dplyr::arrange(contrast_type,
                 desc(expression),
                 tissue,
                 dpi) %>%
  dplyr::relocate(genes,
                  .after = last_col()) %>% 
  readr::write_csv("04_differential-expression-analysis/01_output/top-de-gene-symbols.csv")

### count number of unique genes and gene symbols in top 10
top_de_gene_counts <- top_de_genes %>%
  dplyr::select(symbol) %>%
  tidyr::drop_na() %>%
  dplyr::group_by(symbol) %>%
  dplyr::summarize(n_contrasts_de = n()) %>% 
  readr::write_csv("04_differential-expression-analysis/01_output/top-de-gene-counts.csv")

top_de_gene_symbol_counts <- top_de_gene_symbols %>%
  dplyr::select(symbol) %>%
  tidyr::drop_na() %>%
  dplyr::group_by(symbol) %>%
  dplyr::summarize(n_contrasts_de = n()) %>% 
  readr::write_csv("04_differential-expression-analysis/01_output/top-de-gene-symbol-counts.csv")

top_de_gene_symbols %>%
  dplyr::select(symbol) %>%
  tidyr::drop_na() %>%
  dplyr::count()

top_de_gene_symbol_counts %>%
  dplyr::count()

top_de_gene_symbol_counts %>%
  dplyr::filter(n_contrasts_de == 1) %>%
  dplyr::count()

top_de_gene_symbol_counts %>%
  dplyr::filter(n_contrasts_de != 1) %>%
  dplyr::count()

top_de_gene_symbol_counts %>%
  dplyr::filter(n_contrasts_de != 1) %>%
  dplyr::summarise(mean(n_contrasts_de))

### extract unique gene symbols
top_unique_de_gene_symbols <- top_de_gene_symbol_counts %>%
  dplyr::filter(n_contrasts_de == 1) %>%
  dplyr::left_join(top_de_gene_symbols) %>%
  dplyr::mutate(contrast_type = stringr::str_split_i(set,
                                                     "_",
                                                     1)) %>%
  readr::write_csv("04_differential-expression-analysis/01_output/top-unique-de-gene-symbols.csv")

top_unique_de_gene_symbols %>%
  dplyr::filter(contrast_type == "RESP") %>%
  dplyr::count()

top_unique_de_gene_symbols %>%
  dplyr::filter(contrast_type == "DIRE") %>%
  dplyr::count()

top_unique_de_gene_symbols %>%
  dplyr::filter(contrast_type == "NDAM") %>%
  dplyr::count()

top_unique_de_gene_symbols %>%
  dplyr::filter(contrast_type == "BORA") %>%
  dplyr::count()

### extract top 10 de gene symbols for resp contrasts
top_de_gene_symbol_counts_RESP <- top_de_gene_symbols %>%
  tidyr::drop_na() %>%
  dplyr::mutate(contrast_type = stringr::str_split_i(set,
                                                     "_",
                                                     1),
                expression = stringr::str_split_i(set,
                                                  "_",
                                                  4)) %>%
  dplyr::filter(contrast_type == "RESP") %>%
  dplyr::group_by(symbol) %>%
  dplyr::summarize(n_contrasts_de = n()) %>% 
  readr::write_csv("04_differential-expression-analysis/01_output/top-de-gene-symbol-counts-RESP.csv")

### prepare numbers of de genes for bar plots
de_numbers <- as.data.frame(summary(results)) %>%
  dplyr::rename(de = Var1,
                contrast_id = Var2,
                number = Freq) %>%
  dplyr::mutate(contrast_id = factor(contrast_id,
                                     levels = contrast_ids),
                contrast_type = stringr::str_split_i(contrast_id,
                                                     "_",
                                                     1),
                tissue = stringr::str_split_i(contrast_id,
                                              "_",
                                              2),
                dpi = stringr::str_split_i(contrast_id,
                                           "_",
                                           3),
                dpi = dplyr::case_when(dpi == "00" ~ "0",
                                       .default = as.character(dpi))) %>%
  readr::write_csv("04_differential-expression-analysis/01_output/de-numbers.csv")
de_numbers <- readr::read_csv("04_differential-expression-analysis/01_output/de-numbers.csv")

### bar plot of all contrasts
#### function to generate colour palette
smooth_rainbow <- khroma::colour("smooth rainbow")

#### generate colour palettes 
pal <- smooth_rainbow(31,
                      range = c(0.1,
                                0.9))

full_pal <- smooth_rainbow(34)

viridis <- viridis::viridis(41, direction = -1)

#### set colour palettes
cols <-  c(pal[28],
           pal[17],
           pal[9],
           pal[21])

fill_cols <- c(full_pal[34],
               "#666666",
               pal[5],
               pal[26])

#### set limits
lims <- c("RESP",
          "DIRE",
          "NDAM",
          "BORA")

#### draw plot
bar_plot <- ggplot2::ggplot(de_numbers,
                            ggplot2::aes(interaction(factor(contrast_id),
                                                     factor(dpi)),
                                         colour = tissue,
                                         fill = contrast_type)) +
  ggplot2::geom_col(ggplot2::aes(y = number),
                    data = subset(de_numbers,
                                  de == "Up"),
                    linewidth = 1,
                    width = 0.8) +
  ggplot2::geom_col(ggplot2::aes(y = -number),
                    data = subset(de_numbers,
                                  de == "Down"),
                    linewidth = 1,
                    width = 0.8) +
  ggplot2::annotate("text",
                    label = c("Increased Expression",
                              "Decreased Expression"),
                    hjust = -0.05,
                    x = -Inf,
                    y = c(2750,
                          -2750)) +
  ggplot2::geom_hline(colour = "grey70",
                      linewidth = rel(0.25),
                      yintercept = 0, ) +
  ggplot2::scale_colour_manual(name = "Tissue",
                               values = cols) +
  ggplot2::scale_fill_manual(limits = lims,
                             name = "Contrast",
                             values = fill_cols) +
  ggplot2::scale_x_discrete(guide = ggh4x::guide_axis_nested()) +
  ggplot2::scale_y_continuous(expand = c(0, 0),
                              limits = c(-3000, 3000)) +
  ggplot2::labs(x = "Days Post-Infection",
                y = "Number of Differentially Expressed Genes") +
  ggplot2::guides(colour = ggplot2::guide_legend(order = 2,
                                                 override.aes = list(fill = "white")),
                  fill = ggplot2::guide_legend(order = 1)) +
  ggplot2::theme_light() +
  ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                 axis.ticks.x = ggplot2::element_blank(),
                 ggh4x.axis.nesttext.x = ggplot2::element_text(),
                 ggh4x.axis.nestline = ggplot2::element_line(colour = c(viridis[1],
                                                                        viridis[13],
                                                                        viridis[15],
                                                                        viridis[16],
                                                                        viridis[19],
                                                                        viridis[22],
                                                                        viridis[26],
                                                                        viridis[27],
                                                                        viridis[30],
                                                                        viridis[33],
                                                                        viridis[35],
                                                                        viridis[36]),
                                                             linewidth = 0.5),
                 legend.position = "top",
                 panel.border = element_rect(linewidth = ggplot2::rel(0.5)))

#### save plot
png("04_differential-expression-analysis/02_figures/01_bar-plots/all-bar.png",
    height = 6,
    res = 300,
    units = "in",
    width = 9)
print(bar_plot)
dev.off()
pdf("04_differential-expression-analysis/02_figures/01_bar-plots/all-bar.pdf",
    height = 6,
    width = 9)
print(bar_plot)
dev.off()

### bar plots for each contrast
#### function to draw bar plots for each contrast
bar_plot <- function(contrast, limit){
  # function to generate colour palette
  smooth_rainbow <- khroma::colour("smooth rainbow")
  # generate colour palettes 
  pal <- smooth_rainbow(31,
                        range = c(0.1,
                                  0.9))
  viridis <- viridis::viridis(41, direction = -1)
  # set colour palettes
  cols <-  c(pal[28],
             pal[17],
             pal[9],
             pal[21])
  dpi_cols <- c(viridis[1],
                viridis[13],
                viridis[15],
                viridis[16],
                viridis[19],
                viridis[22],
                viridis[26],
                viridis[27],
                viridis[30],
                viridis[33],
                viridis[35],
                viridis[36])
  # draw plot
  plot <- ggplot2::ggplot(subset(de_numbers,
                                 contrast_type == contrast),
                          ggplot2::aes(interaction(factor(contrast_id),
                                                   factor(dpi)),
                                       colour = tissue,
                                       fill = tissue)) +
    ggplot2::geom_col(ggplot2::aes(y = number),
                      data = subset(de_numbers,
                                    contrast_type == contrast & de == "Up"),
                      linewidth = 1) +
    ggplot2::geom_col(aes(y = -number),
                      data = subset(de_numbers,
                                    contrast_type == contrast & de == "Down"),
                      linewidth = 1) +
    ggplot2::annotate("text",
                      label = c("Increased Expression",
                                "Decreased Expression"),
                      hjust = -0.05,
                      x = -Inf,
                      y = c(limit*0.917,
                            -limit*0.917)) +
    ggplot2::geom_hline(colour = "grey70",
                        linewidth = ggplot2::rel(0.25),
                        yintercept = 0, ) +
    ggplot2::scale_colour_manual(name = "Tissue",
                                 values = cols) +
    ggplot2::scale_fill_manual(name = "Tissue",
                               values = cols) +
    ggplot2::scale_x_discrete(guide = ggh4x::guide_axis_nested()) +
    ggplot2::scale_y_continuous(expand = c(0,
                                           0),
                                limits = c(-limit,
                                           limit)) +
    ggplot2::labs(x = "Days Post-Infection",
                  y = "Number of Differentially Expressed Genes") +
    ggplot2::theme_light() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   ggh4x.axis.nesttext.x = ggplot2::element_text(),
                   ggh4x.axis.nestline = ggplot2::element_line(colour = dpi_cols[(13-length(unique(subset(de_numbers,
                                                                                                          contrast_type == contrast)$dpi))):12],
                                                               linewidth = 0.5),
                   legend.position = "top",
                   panel.border = element_rect(linewidth = ggplot2::rel(0.5)))
  # save plot
  png(paste0("04_differential-expression-analysis/02_figures/01_bar-plots/",
             contrast,
             "-bar.png"),
      height = 6,
      res = 300,
      units = "in",
      width = 9)
  print(plot)
  dev.off()
  pdf(paste0("04_differential-expression-analysis/02_figures/01_bar-plots/",
             contrast,
             "-bar.pdf"),
      height = 6,
      width = 9)
  print(plot)
  dev.off()
}

#### apply function
bar_plot("RESP",
         limit = 500)

bar_plot("DIRE",
         limit = 1100)

bar_plot("NDAM",
         limit = 3000)

bar_plot("BORA",
         limit = 3000)

### upset plots to show numbers of genes in common
#### prepare data for upset plots
upset <- readr::read_csv("04_differential-expression-analysis/01_output/de-results.csv") %>%
  dplyr::select(dplyr::starts_with("Results.")) %>%
  dplyr::rename_with(~ stringr::str_remove_all(.,
                                               pattern = "Results.")) %>%
  dplyr::mutate(dplyr::across(dplyr::everything(), ~ dplyr::case_match(.x,
                                                                       c(1, -1) ~ TRUE,
                                                                       0 ~ FALSE)))

upset_data <- as.data.frame(contrast_ids) %>%
  dplyr::transmute(set = contrast_ids,
                   intersection = contrast_ids,
                   contrast_type = stringr::str_split_i(contrast_ids,
                                                        "_",
                                                        1),
                   tissue = stringr::str_split_i(contrast_ids,
                                                 "_",
                                                 2),
                   dpi  = stringr::str_split_i(contrast_ids,
                                               "_",
                                               3))

#### set labels
all_labels <- c("RESP_BL_14" = "RESP BL 14",
                "RESP_BL_25" = "RESP BL 25",
                "RESP_BL_34" = "RESP BL 34",
                "RESP_LI_12" = "RESP LI 12",
                "RESP_LI_15" = "RESP LI 15",
                "RESP_LI_18" = "RESP LI 18",
                "RESP_LI_21" = "RESP LI 21",
                "RESP_LI_26" = "RESP LI 26",
                "RESP_LI_29" = "RESP LI 29",
                "RESP_LI_32" = "RESP LI 32",
                "RESP_LI_35" = "RESP LI 35",
                "RESP_LN_21" = "RESP LN 21",
                "RESP_LN_35" = "RESP LN 35",
                "RESP_SP_21" = "RESP SP 21",
                "RESP_SP_35" = "RESP SP 35",
                "DIRE_BL_00" = "DIRE BL 00",
                "DIRE_BL_14" = "DIRE BL 14",
                "DIRE_BL_25" = "DIRE BL 25",
                "DIRE_BL_34" = "DIRE BL 34",
                "DIRE_LI_00" = "DIRE LI 00",
                "DIRE_LI_12" = "DIRE LI 12",
                "DIRE_LI_15" = "DIRE LI 15",
                "DIRE_LI_18" = "DIRE LI 18",
                "DIRE_LI_21" = "DIRE LI 21",
                "DIRE_LI_26" = "DIRE LI 26",
                "DIRE_LI_29" = "DIRE LI 29",
                "DIRE_LI_32" = "DIRE LI 32",
                "DIRE_LI_35" = "DIRE LI 35",
                "DIRE_LN_00" = "DIRE LN 00",
                "DIRE_LN_21" = "DIRE LN 21",
                "DIRE_LN_35" = "DIRE LN 35",
                "DIRE_SP_00" = "DIRE SP 00",
                "DIRE_SP_21" = "DIRE SP 21",
                "DIRE_SP_35" = "DIRE SP 35",
                "NDAM_BL_14" = "NDAM BL 14",
                "NDAM_BL_25" = "NDAM BL 25",
                "NDAM_BL_34" = "NDAM BL 34",
                "NDAM_LI_12" = "NDAM LI 12",
                "NDAM_LI_15" = "NDAM LI 15",
                "NDAM_LI_18" = "NDAM LI 18",
                "NDAM_LI_21" = "NDAM LI 21",
                "NDAM_LI_26" = "NDAM LI 26",
                "NDAM_LI_29" = "NDAM LI 29",
                "NDAM_LI_32" = "NDAM LI 32",
                "NDAM_LI_35" = "NDAM LI 35",
                "NDAM_LN_21" = "NDAM LN 21",
                "NDAM_LN_35" = "NDAM LN 35",
                "NDAM_SP_21" = "NDAM SP 21",
                "NDAM_SP_35" = "NDAM SP 35",
                "BORA_BL_14" = "BORA BL 14",
                "BORA_BL_25" = "BORA BL 25",
                "BORA_BL_34" = "BORA BL 34",
                "BORA_LI_12" = "BORA LI 12",
                "BORA_LI_15" = "BORA LI 15",
                "BORA_LI_18" = "BORA LI 18",
                "BORA_LI_21" = "BORA LI 21",
                "BORA_LI_26" = "BORA LI 26",
                "BORA_LI_29" = "BORA LI 29",
                "BORA_LI_32" = "BORA LI 32",
                "BORA_LI_35" = "BORA LI 35",
                "BORA_LN_21" = "BORA LN 21",
                "BORA_LN_35" = "BORA LN 35",
                "BORA_SP_21" = "BORA SP 21",
                "BORA_SP_35" = "BORA SP 35")

#### set queries
all_queries <- list(ComplexUpset::upset_query(set = "DIRE_BL_34",
                                              color = pal[28],
                                              fill = "#666666",
                                              only_components = "overall_sizes"),
                    ComplexUpset::upset_query(intersect = "DIRE_BL_34",
                                              color = pal[28],
                                              fill = "#666666",
                                              only_components = "Intersection size"),
                    ComplexUpset::upset_query(set = "DIRE_LI_00",
                                              color = pal[17],
                                              fill = "#666666",
                                              only_components = "overall_sizes"),
                    ComplexUpset::upset_query(intersect = "DIRE_LI_00",
                                              color = pal[17],
                                              fill = "#666666",
                                              only_components = "Intersection size"),
                    ComplexUpset::upset_query(set = "NDAM_BL_25",
                                              color = pal[28],
                                              fill = pal[5],
                                              only_components = "overall_sizes"),
                    ComplexUpset::upset_query(set = "NDAM_BL_34",
                                              color = pal[28],
                                              fill = pal[5],
                                              only_components = "overall_sizes"),
                    ComplexUpset::upset_query(set = "NDAM_LI_12",
                                              color = pal[17],
                                              fill = pal[5],
                                              only_components = "overall_sizes"),
                    ComplexUpset::upset_query(intersect = "NDAM_LI_12",
                                              color = pal[17],
                                              fill = pal[5],
                                              only_components = "Intersection size"),
                    ComplexUpset::upset_query(set = "NDAM_LI_15",
                                              color = pal[17],
                                              fill = pal[5],
                                              only_components = "overall_sizes"),
                    ComplexUpset::upset_query(set = "NDAM_LI_18",
                                              color = pal[17],
                                              fill = pal[5],
                                              only_components = "overall_sizes"),
                    ComplexUpset::upset_query(set = "NDAM_LI_21",
                                              color = pal[17],
                                              fill = pal[5],
                                              only_components = "overall_sizes"),
                    ComplexUpset::upset_query(intersect = "NDAM_LI_21",
                                              color = pal[17],
                                              fill = pal[5],
                                              only_components = "Intersection size"),
                    ComplexUpset::upset_query(set = "NDAM_LI_26",
                                              color = pal[17],
                                              fill = pal[5],
                                              only_components = "overall_sizes"),
                    ComplexUpset::upset_query(set = "NDAM_LI_29",
                                              color = pal[17],
                                              fill = pal[5],
                                              only_components = "overall_sizes"),
                    ComplexUpset::upset_query(set = "NDAM_LI_32",
                                              color = pal[17],
                                              fill = pal[5],
                                              only_components = "overall_sizes"),
                    ComplexUpset::upset_query(set = "NDAM_LI_35",
                                              color = pal[17],
                                              fill = pal[5],
                                              only_components = "overall_sizes"),
                    ComplexUpset::upset_query(intersect = "NDAM_LI_35",
                                              color = pal[17],
                                              fill = pal[5],
                                              only_components = "Intersection size"),
                    ComplexUpset::upset_query(set = "NDAM_LN_21",
                                              color = pal[9],
                                              fill = pal[5],
                                              only_components = "overall_sizes"),
                    ComplexUpset::upset_query(set = "NDAM_LN_35",
                                              color = pal[9],
                                              fill = pal[5],
                                              only_components = "overall_sizes"),
                    ComplexUpset::upset_query(intersect = "NDAM_LN_35",
                                              color = pal[9],
                                              fill = pal[5],
                                              only_components = "Intersection size"),
                    ComplexUpset::upset_query(set = "NDAM_SP_21",
                                              color = pal[21],
                                              fill = pal[5],
                                              only_components = "overall_sizes"),
                    ComplexUpset::upset_query(set = "NDAM_SP_35",
                                              color = pal[21],
                                              fill = pal[5],
                                              only_components = "overall_sizes"),
                    ComplexUpset::upset_query(intersect = c("NDAM_LN_21",
                                                            "NDAM_LN_35"),
                                              color = pal[9],
                                              fill = pal[5],
                                              only_components = "Intersection size"),
                    ComplexUpset::upset_query(set = "BORA_BL_25",
                                              color = pal[28],
                                              fill = pal[26],
                                              only_components = "overall_sizes"),
                    ComplexUpset::upset_query(set = "BORA_BL_34",
                                              color = pal[28],
                                              fill = pal[26],
                                              only_components = "overall_sizes"),
                    ComplexUpset::upset_query(set = "BORA_LI_15",
                                              color = pal[17],
                                              fill = pal[26],
                                              only_components = "overall_sizes"),
                    ComplexUpset::upset_query(set = "BORA_LI_18",
                                              color = pal[17],
                                              fill = pal[26],
                                              only_components = "overall_sizes"),
                    ComplexUpset::upset_query(set = "BORA_LI_21",
                                              color = pal[17],
                                              fill = pal[26],
                                              only_components = "overall_sizes"),
                    ComplexUpset::upset_query(intersect = "BORA_LI_21",
                                              color = pal[17],
                                              fill = pal[26],
                                              only_components = "Intersection size"),
                    ComplexUpset::upset_query(set = "BORA_LI_26",
                                              color = pal[17],
                                              fill = pal[26],
                                              only_components = "overall_sizes"),
                    ComplexUpset::upset_query(set = "BORA_LI_29",
                                              color = pal[17],
                                              fill = pal[26],
                                              only_components = "overall_sizes"),
                    ComplexUpset::upset_query(set = "BORA_LI_32",
                                              color = pal[17],
                                              fill = pal[26],
                                              only_components = "overall_sizes"),
                    ComplexUpset::upset_query(set = "BORA_LI_35",
                                              color = pal[17],
                                              fill = pal[26],
                                              only_components = "overall_sizes"),
                    ComplexUpset::upset_query(intersect = "BORA_LI_35",
                                              color = pal[17],
                                              fill = pal[26],
                                              only_components = "Intersection size"),
                    ComplexUpset::upset_query(set = "BORA_LN_21",
                                              color = pal[9],
                                              fill = pal[26],
                                              only_components = "overall_sizes"),
                    ComplexUpset::upset_query(set = "BORA_LN_35",
                                              color = pal[9],
                                              fill = pal[26],
                                              only_components = "overall_sizes"),
                    ComplexUpset::upset_query(set = "BORA_SP_21",
                                              color = pal[21],
                                              fill = pal[26],
                                              only_components = "overall_sizes"),
                    ComplexUpset::upset_query(set = "BORA_SP_35",
                                              color = pal[21],
                                              fill = pal[26],
                                              only_components = "overall_sizes"),
                    ComplexUpset::upset_query(intersect = c("BORA_LI_21",
                                                            "BORA_LI_35"),
                                              color = pal[17],
                                              fill = pal[26],
                                              only_components = "Intersection size"),
                    ComplexUpset::upset_query(intersect = c("NDAM_BL_25",
                                                            "NDAM_BL_34",
                                                            "BORA_BL_25",
                                                            "BORA_BL_34"),
                                              color = pal[28],
                                              only_components = "Intersection size"),
                    ComplexUpset::upset_query(intersect = c("NDAM_LI_21",
                                                            "NDAM_LI_35",
                                                            "BORA_LI_21",
                                                            "BORA_LI_35"),
                                              color = pal[17],
                                              only_components = "Intersection size"),
                    ComplexUpset::upset_query(intersect = c("NDAM_LI_21",
                                                            "BORA_LI_21",
                                                            "BORA_LI_35"),
                                              color = pal[17],
                                              only_components = "Intersection size"),
                    ComplexUpset::upset_query(intersect = c("NDAM_SP_21",
                                                            "NDAM_SP_35",
                                                            "BORA_SP_21",
                                                            "BORA_SP_35"),
                                              color = pal[21],
                                              only_components = "Intersection size"),
                    ComplexUpset::upset_query(intersect = c("NDAM_LI_15",
                                                            "NDAM_LI_18",
                                                            "NDAM_LI_21",
                                                            "NDAM_LI_26",
                                                            "NDAM_LI_29",
                                                            "NDAM_LI_32",
                                                            "NDAM_LI_35",
                                                            "BORA_LI_15",
                                                            "BORA_LI_18",
                                                            "BORA_LI_21",
                                                            "BORA_LI_26",
                                                            "BORA_LI_29",
                                                            "BORA_LI_32",
                                                            "BORA_LI_35"),
                                              color = pal[17],
                                              only_components = "Intersection size"),
                    ComplexUpset::upset_query(intersect = c("NDAM_LI_18",
                                                            "NDAM_LI_21",
                                                            "NDAM_LI_26",
                                                            "NDAM_LI_29",
                                                            "NDAM_LI_32",
                                                            "NDAM_LI_35",
                                                            "BORA_LI_18",
                                                            "BORA_LI_21",
                                                            "BORA_LI_26",
                                                            "BORA_LI_29",
                                                            "BORA_LI_32",
                                                            "BORA_LI_35"),
                                              color = pal[17],
                                              only_components = "Intersection size"))

#### select columns
upset_cols <- colnames(upset)

#### draw plot
upset_plot_all <- ComplexUpset::upset(upset,
                                      intersect = upset_cols,
                                      encode_sets = F,
                                      guide = "over",
                                      min_degree = 1,
                                      n_intersections = 20,
                                      name = "Intersections",
                                      height_ratio = 1.5,
                                      width_ratio = 0.25,
                                      labeller = ggplot2::as_labeller(all_labels),
                                      queries = all_queries,
                                      base_annotations = list('Intersection size' = ComplexUpset::intersection_size(mapping = ggplot2::aes(fill = "#000000"),
                                                                                                                    linewidth = 1,
                                                                                                                    width = 0.8,
                                                                                                                    text = list(size = ggplot2::rel(3))) +
                                                                ggplot2::scale_y_continuous(expand = c(0,
                                                                                                       0),
                                                                                            limits = c(0,
                                                                                                       125)) +
                                                                ggplot2::scale_fill_identity() +
                                                                ggplot2::ylab("Number of Differentially\nExpressed Genes in Common") ),
                                      matrix = ComplexUpset::intersection_matrix(geom = ggplot2::geom_point(size = 1.5)),
                                      set_sizes = ComplexUpset::upset_set_size(geom = ggplot2::geom_bar(linewidth = 1,
                                                                                                        width = 0.75)) +
                                        scale_y_reverse(expand = c(0,
                                                                   0),
                                                        limits = c(5500,
                                                                   0)) +
                                        ylab("Total Number of\nDifferentially Expressed Genes"),
                                      stripes = ComplexUpset::upset_stripes(geom = ggplot2::geom_segment(alpha = 0.5,
                                                                                                         linewidth = 3.5),
                                                                            mapping = ggplot2::aes(colour = tissue),
                                                                            colors = c("BL" = pal[28],
                                                                                       "LI" = pal[17],
                                                                                       "LN" = pal[9],
                                                                                       "SP" = pal[21]),
                                                                            data = upset_data),
                                      themes = list('Intersection size' = ggplot2::theme_light() +
                                                      ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                                                                     axis.ticks.x = ggplot2::element_blank(),
                                                                     axis.title.x = ggplot2::element_blank()),
                                                    intersections_matrix = ggplot2::theme_light() +
                                                      ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                                                                     axis.ticks = ggplot2::element_blank(),
                                                                     axis.title.y = ggplot2::element_blank(),
                                                                     axis.text.y = ggplot2::element_text(family = "mono",
                                                                                                         margin = ggplot2::margin(l = -40)),
                                                                     legend.position = "none"),
                                                    overall_sizes = ggplot2::theme_light() +
                                                      ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                                                                     axis.ticks.y = ggplot2::element_blank(),
                                                                     axis.title.y = ggplot2::element_blank(),
                                                                     plot.margin = ggplot2::margin(r = 8))))

#### draw legend plot
legend_plot <- ggplot2::ggplot(upset_data,
                               ggplot2::aes(x = set,
                                            y = dpi,
                                            colour = tissue,
                                            fill = contrast_type)) +
  ggplot2::geom_col(linewidth = 1) +
  ggplot2::scale_colour_manual(name = "Tissue",
                               values = cols) +
  ggplot2::scale_fill_manual(limits = lims,
                             name = "Contrast",
                             values = fill_cols) +
  ggplot2::guides(colour = ggplot2::guide_legend(order = 2,
                                                 override.aes = list(fill = ggplot2::alpha(cols, 0.5),
                                                                     stroke = 2),),
                  fill = ggplot2::guide_legend(order = 1)) +
  ggplot2::theme_light()

#### set layout  
layout <- "AB
           CE"

#### combine plots
plot <- upset_plot_all +
  legend_plot +
  patchwork::plot_layout(design = layout,
                         guides = "collect") &
  ggplot2::theme(legend.box = "horizontal")

#### save plot
png("04_differential-expression-analysis/02_figures/02_upset-plots/all-upset.png",
    height = 6,
    res = 300,
    units = "in",
    width = 9)
print(plot)
dev.off()
pdf("04_differential-expression-analysis/02_figures/02_upset-plots/all-upset.pdf",
    height = 6,
    width = 9)
print(plot)
dev.off()

#### function to draw upset plots for each contrast type
upset_plot <- function(contrast_type, total_limit, common_limit, stripe_width = 4){
  # set labels
  labels <- c("RESP_BL_14" = "BL 14",
              "RESP_BL_25" = "BL 25",
              "RESP_BL_34" = "BL 34",
              "RESP_LI_12" = "LI 12",
              "RESP_LI_15" = "LI 15",
              "RESP_LI_18" = "LI 18",
              "RESP_LI_21" = "LI 21",
              "RESP_LI_26" = "LI 26",
              "RESP_LI_29" = "LI 29",
              "RESP_LI_32" = "LI 32",
              "RESP_LI_35" = "LI 35",
              "RESP_LN_21" = "LN 21",
              "RESP_LN_35" = "LN 35",
              "RESP_SP_21" = "SP 21",
              "RESP_SP_35" = "SP 35",
              "DIRE_BL_00" = "BL 00",
              "DIRE_BL_14" = "BL 14",
              "DIRE_BL_25" = "BL 25",
              "DIRE_BL_34" = "BL 34",
              "DIRE_LI_00" = "LI 00",
              "DIRE_LI_12" = "LI 12",
              "DIRE_LI_15" = "LI 15",
              "DIRE_LI_18" = "LI 18",
              "DIRE_LI_21" = "LI 21",
              "DIRE_LI_26" = "LI 26",
              "DIRE_LI_29" = "LI 29",
              "DIRE_LI_32" = "LI 32",
              "DIRE_LI_35" = "LI 35",
              "DIRE_LN_00" = "LN 00",
              "DIRE_LN_21" = "LN 21",
              "DIRE_LN_35" = "LN 35",
              "DIRE_SP_00" = "SP 00",
              "DIRE_SP_21" = "SP 21",
              "DIRE_SP_35" = "SP 35",
              "NDAM_BL_14" = "BL 14",
              "NDAM_BL_25" = "BL 25",
              "NDAM_BL_34" = "BL 34",
              "NDAM_LI_12" = "LI 12",
              "NDAM_LI_15" = "LI 15",
              "NDAM_LI_18" = "LI 18",
              "NDAM_LI_21" = "LI 21",
              "NDAM_LI_26" = "LI 26",
              "NDAM_LI_29" = "LI 29",
              "NDAM_LI_32" = "LI 32",
              "NDAM_LI_35" = "LI 35",
              "NDAM_LN_21" = "LN 21",
              "NDAM_LN_35" = "LN 35",
              "NDAM_SP_21" = "SP 21",
              "NDAM_SP_35" = "SP 35",
              "BORA_BL_14" = "BL 14",
              "BORA_BL_25" = "BL 25",
              "BORA_BL_34" = "BL 34",
              "BORA_LI_12" = "LI 12",
              "BORA_LI_15" = "LI 15",
              "BORA_LI_18" = "LI 18",
              "BORA_LI_21" = "LI 21",
              "BORA_LI_26" = "LI 26",
              "BORA_LI_29" = "LI 29",
              "BORA_LI_32" = "LI 32",
              "BORA_LI_35" = "LI 35",
              "BORA_LN_21" = "LN 21",
              "BORA_LN_35" = "LN 35",
              "BORA_SP_21" = "SP 21",
              "BORA_SP_35" = "SP 35")
  # set queries
  RESP_queries <- list(ComplexUpset::upset_query(set = "RESP_BL_14",
                                                 color = pal[28],
                                                 fill = viridis[15],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(intersect = "RESP_BL_14",
                                                 color = pal[28],
                                                 fill = viridis[15],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(set = "RESP_BL_25",
                                                 color = pal[28],
                                                 fill = viridis[26],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(intersect = "RESP_BL_25",
                                                 color = pal[28],
                                                 fill = viridis[26],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(set = "RESP_BL_34",
                                                 color = pal[28],
                                                 fill = viridis[35],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(intersect = "RESP_BL_34",
                                                 color = pal[28],
                                                 fill = viridis[35],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(set = "RESP_LI_12",
                                                 color = pal[17],
                                                 fill = viridis[13],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(intersect = "RESP_LI_12",
                                                 color = pal[17],
                                                 fill = viridis[13],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(set = "RESP_LI_15",
                                                 color = pal[17],
                                                 fill = viridis[16],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(intersect = "RESP_LI_15",
                                                 color = pal[17],
                                                 fill = viridis[16],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(set = "RESP_LI_18",
                                                 color = pal[17],
                                                 fill = viridis[19],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(intersect = "RESP_LI_18",
                                                 color = pal[17],
                                                 fill = viridis[19],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(set = "RESP_LI_21",
                                                 color = pal[17],
                                                 fill = viridis[22],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(intersect = "RESP_LI_21",
                                                 color = pal[17],
                                                 fill = viridis[22],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(set = "RESP_LI_29",
                                                 color = pal[17],
                                                 fill = viridis[30],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(intersect = "RESP_LI_29",
                                                 color = pal[17],
                                                 fill = viridis[30],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(set = "RESP_LI_32",
                                                 color = pal[17],
                                                 fill = viridis[33],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(intersect = "RESP_LI_32",
                                                 color = pal[17],
                                                 fill = viridis[33],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(set = "RESP_LI_35",
                                                 color = pal[17],
                                                 fill = viridis[36],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(intersect = "RESP_LI_35",
                                                 color = pal[17],
                                                 fill = viridis[36],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(set = "RESP_LN_21",
                                                 color = pal[9],
                                                 fill = viridis[22],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(intersect = "RESP_LN_21",
                                                 color = pal[9],
                                                 fill = viridis[22],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(set = "RESP_LN_35",
                                                 color = pal[9],
                                                 fill = viridis[36],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(intersect = "RESP_LN_35",
                                                 color = pal[9],
                                                 fill = viridis[36],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(set = "RESP_SP_21",
                                                 color = pal[21],
                                                 fill = viridis[22],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(intersect = "RESP_SP_21",
                                                 color = pal[21],
                                                 fill = viridis[22],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(set = "RESP_SP_35",
                                                 color = pal[21],
                                                 fill = viridis[36],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(intersect = "RESP_SP_35",
                                                 color = pal[21],
                                                 fill = viridis[36],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(intersect = c("RESP_BL_14",
                                                               "RESP_BL_34"),
                                                 color = pal[28],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(intersect = c("RESP_BL_25",
                                                               "RESP_BL_34"),
                                                 color = pal[28],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(intersect = c("RESP_LI_29",
                                                               "RESP_LI_35"),
                                                 color = pal[17],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(intersect = c("RESP_LN_21",
                                                               "RESP_LN_35"),
                                                 color = pal[9],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(intersect = c("RESP_SP_21",
                                                               "RESP_SP_35"),
                                                 color = pal[21],
                                                 only_components = "Intersection size"))
  DIRE_queries <- list(ComplexUpset::upset_query(set = "DIRE_BL_00",
                                                 color = pal[28],
                                                 fill = viridis[1],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(intersect = "DIRE_BL_00",
                                                 color = pal[28],
                                                 fill = viridis[1],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(set = "DIRE_BL_14",
                                                 color = pal[28],
                                                 fill = viridis[15],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(intersect = "DIRE_BL_14",
                                                 color = pal[28],
                                                 fill = viridis[15],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(set = "DIRE_BL_25",
                                                 color = pal[28],
                                                 fill = viridis[26],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(intersect = "DIRE_BL_25",
                                                 color = pal[28],
                                                 fill = viridis[26],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(set = "DIRE_BL_34",
                                                 color = pal[28],
                                                 fill = viridis[35],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(intersect = "DIRE_BL_34",
                                                 color = pal[28],
                                                 fill = viridis[35],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(set = "DIRE_LI_00",
                                                 color = pal[17],
                                                 fill = viridis[1],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(intersect = "DIRE_LI_00",
                                                 color = pal[17],
                                                 fill = viridis[1],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(set = "DIRE_LI_12",
                                                 color = pal[17],
                                                 fill = viridis[13],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(intersect = "DIRE_LI_12",
                                                 color = pal[17],
                                                 fill = viridis[13],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(set = "DIRE_LI_15",
                                                 color = pal[17],
                                                 fill = viridis[16],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(intersect = "DIRE_LI_15",
                                                 color = pal[17],
                                                 fill = viridis[16],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(set = "DIRE_LI_29",
                                                 color = pal[17],
                                                 fill = viridis[30],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(intersect = "DIRE_LI_29",
                                                 color = pal[17],
                                                 fill = viridis[30],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(set = "DIRE_LI_32",
                                                 color = pal[17],
                                                 fill = viridis[33],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(intersect = "DIRE_LI_32",
                                                 color = pal[17],
                                                 fill = viridis[33],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(set = "DIRE_LI_35",
                                                 color = pal[17],
                                                 fill = viridis[36],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(intersect = "DIRE_LI_35",
                                                 color = pal[17],
                                                 fill = viridis[36],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(set = "DIRE_LN_00",
                                                 color = pal[9],
                                                 fill = viridis[1],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(intersect = "DIRE_LN_00",
                                                 color = pal[9],
                                                 fill = viridis[1],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(set = "DIRE_LN_21",
                                                 color = pal[9],
                                                 fill = viridis[22],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(intersect = "DIRE_LN_21",
                                                 color = pal[9],
                                                 fill = viridis[22],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(set = "DIRE_LN_35",
                                                 color = pal[9],
                                                 fill = viridis[36],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(intersect = "DIRE_LN_35",
                                                 color = pal[9],
                                                 fill = viridis[36],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(set = "DIRE_SP_00",
                                                 color = pal[21],
                                                 fill = viridis[1],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(intersect = "DIRE_SP_00",
                                                 color = pal[21],
                                                 fill = viridis[1],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(set = "DIRE_SP_35",
                                                 color = pal[21],
                                                 fill = viridis[36],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(intersect = "DIRE_SP_35",
                                                 color = pal[21],
                                                 fill = viridis[36],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(intersect = c("DIRE_BL_00",
                                                               "DIRE_BL_25"),
                                                 color = pal[28],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(intersect = c("DIRE_BL_00",
                                                               "DIRE_BL_34"),
                                                 color = pal[28],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(intersect = c("DIRE_BL_00",
                                                               "DIRE_BL_25",
                                                               "DIRE_BL_34"),
                                                 color = pal[28],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(intersect = c("DIRE_BL_25",
                                                               "DIRE_BL_34"),
                                                 color = pal[28],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(intersect = c("DIRE_LI_00",
                                                               "DIRE_LI_12"),
                                                 color = pal[17],
                                                 only_components = "Intersection size"))
  NDAM_queries <- list(ComplexUpset::upset_query(set = "NDAM_BL_14",
                                                 color = pal[28],
                                                 fill = viridis[15],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(set = "NDAM_BL_25",
                                                 color = pal[28],
                                                 fill = viridis[26],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(set = "NDAM_BL_34",
                                                 color = pal[28],
                                                 fill = viridis[35],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(intersect = "NDAM_BL_34",
                                                 color = pal[28],
                                                 fill = viridis[35],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(set = "NDAM_LI_15",
                                                 color = pal[17],
                                                 fill = viridis[16],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(set = "NDAM_LI_18",
                                                 color = pal[17],
                                                 fill = viridis[19],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(set = "NDAM_LI_21",
                                                 color = pal[17],
                                                 fill = viridis[22],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(intersect = "NDAM_LI_21",
                                                 color = pal[17],
                                                 fill = viridis[22],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(set = "NDAM_LI_26",
                                                 color = pal[17],
                                                 fill = viridis[27],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(set = "NDAM_LI_29",
                                                 color = pal[17],
                                                 fill = viridis[30],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(set = "NDAM_LI_32",
                                                 color = pal[17],
                                                 fill = viridis[33],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(set = "NDAM_LI_35",
                                                 color = pal[17],
                                                 fill = viridis[36],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(intersect = "NDAM_LI_35",
                                                 color = pal[17],
                                                 fill = viridis[36],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(set = "NDAM_LN_21",
                                                 color = pal[9],
                                                 fill = viridis[22],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(set = "NDAM_LN_35",
                                                 color = pal[9],
                                                 fill = viridis[36],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(intersect = "NDAM_LN_35",
                                                 color = pal[9],
                                                 fill = viridis[36],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(set = "NDAM_SP_21",
                                                 color = pal[21],
                                                 fill = viridis[22],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(intersect = "NDAM_SP_21",
                                                 color = pal[21],
                                                 fill = viridis[22],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(set = "NDAM_SP_35",
                                                 color = pal[21],
                                                 fill = viridis[36],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(intersect = c("NDAM_BL_14",
                                                               "NDAM_BL_25",
                                                               "NDAM_BL_34"),
                                                 color = pal[28],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(intersect = c("NDAM_BL_25",
                                                               "NDAM_BL_34"),
                                                 color = pal[28],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(intersect = c("NDAM_LI_15",
                                                               "NDAM_LI_18",
                                                               "NDAM_LI_21",
                                                               "NDAM_LI_26",
                                                               "NDAM_LI_29",
                                                               "NDAM_LI_32",
                                                               "NDAM_LI_35"),
                                                 color = pal[17],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(intersect = c("NDAM_LI_18",
                                                               "NDAM_LI_21",
                                                               "NDAM_LI_26",
                                                               "NDAM_LI_29",
                                                               "NDAM_LI_32",
                                                               "NDAM_LI_35"),
                                                 color = pal[17],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(intersect = c("NDAM_LI_21",
                                                               "NDAM_LI_35"),
                                                 color = pal[17],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(intersect = c("NDAM_LN_21",
                                                               "NDAM_LN_35"),
                                                 color = pal[9],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(intersect = c("NDAM_SP_21",
                                                               "NDAM_SP_35"),
                                                 color = pal[21],
                                                 only_components = "Intersection size"))
  BORA_queries <- list(ComplexUpset::upset_query(set = "BORA_BL_25",
                                                 color = pal[28],
                                                 fill = viridis[26],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(intersect = "BORA_BL_25",
                                                 color = pal[28],
                                                 fill = viridis[26],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(set = "BORA_BL_34",
                                                 color = pal[28],
                                                 fill = viridis[35],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(intersect = "BORA_BL_34",
                                                 color = pal[28],
                                                 fill = viridis[35],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(set = "BORA_LI_15",
                                                 color = pal[17],
                                                 fill = viridis[16],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(set = "BORA_LI_18",
                                                 color = pal[17],
                                                 fill = viridis[19],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(set = "BORA_LI_21",
                                                 color = pal[17],
                                                 fill = viridis[22],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(intersect = "BORA_LI_21",
                                                 color = pal[17],
                                                 fill = viridis[22],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(set = "BORA_LI_26",
                                                 color = pal[17],
                                                 fill = viridis[27],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(set = "BORA_LI_29",
                                                 color = pal[17],
                                                 fill = viridis[30],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(intersect = "BORA_LI_29",
                                                 color = pal[17],
                                                 fill = viridis[30],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(set = "BORA_LI_32",
                                                 color = pal[17],
                                                 fill = viridis[33],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(set = "BORA_LI_35",
                                                 color = pal[17],
                                                 fill = viridis[36],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(intersect = "BORA_LI_35",
                                                 color = pal[17],
                                                 fill = viridis[36],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(set = "BORA_LN_21",
                                                 color = pal[9],
                                                 fill = viridis[22],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(set = "BORA_LN_35",
                                                 color = pal[9],
                                                 fill = viridis[36],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(intersect = "BORA_LN_35",
                                                 color = pal[9],
                                                 fill = viridis[36],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(set = "BORA_SP_21",
                                                 color = pal[21],
                                                 fill = viridis[22],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(intersect = "BORA_SP_21",
                                                 color = pal[21],
                                                 fill = viridis[22],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(set = "BORA_SP_35",
                                                 color = pal[21],
                                                 fill = viridis[36],
                                                 only_components = "overall_sizes"),
                       ComplexUpset::upset_query(intersect = "BORA_SP_35",
                                                 color = pal[21],
                                                 fill = viridis[36],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(intersect = c("BORA_BL_25",
                                                               "BORA_BL_34"),
                                                 color = pal[28],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(intersect = c("BORA_LI_15",
                                                               "BORA_LI_18",
                                                               "BORA_LI_21",
                                                               "BORA_LI_26",
                                                               "BORA_LI_29",
                                                               "BORA_LI_32",
                                                               "BORA_LI_35"),
                                                 color = pal[17],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(intersect = c("BORA_LI_18",
                                                               "BORA_LI_21",
                                                               "BORA_LI_26",
                                                               "BORA_LI_29",
                                                               "BORA_LI_32",
                                                               "BORA_LI_35"),
                                                 color = pal[17],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(intersect = c("BORA_LI_21",
                                                               "BORA_LI_26",
                                                               "BORA_LI_29",
                                                               "BORA_LI_32",
                                                               "BORA_LI_35"),
                                                 color = pal[17],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(intersect = c("BORA_LI_21",
                                                               "BORA_LI_29",
                                                               "BORA_LI_32",
                                                               "BORA_LI_35"),
                                                 color = pal[17],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(intersect = c("BORA_LI_21",
                                                               "BORA_LI_29",
                                                               "BORA_LI_35"),
                                                 color = pal[17],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(intersect = c("BORA_LI_21",
                                                               "BORA_LI_35"),
                                                 color = pal[17],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(intersect = c("BORA_LN_21",
                                                               "BORA_LN_35"),
                                                 color = pal[9],
                                                 only_components = "Intersection size"),
                       ComplexUpset::upset_query(intersect = c("BORA_SP_21",
                                                               "BORA_SP_35"),
                                                 color = pal[21],
                                                 only_components = "Intersection size"))
  # select columns
  upset_cols <- upset %>%
    select(starts_with(contrast_type)) %>% 
    colnames()
  # draw plot
  upset_plot <- ComplexUpset::upset(upset,
                                    intersect = upset_cols,
                                    encode_sets = F,
                                    guide = "over",
                                    min_degree = 1,
                                    n_intersections = 20,
                                    name = "Intersections",
                                    width_ratio = 0.25,
                                    labeller = ggplot2::as_labeller(labels),
                                    queries = get(paste0(contrast_type,
                                                         "_queries")),
                                    base_annotations = list('Intersection size' = ComplexUpset::intersection_size(mapping = ggplot2::aes(fill = "#000000"),
                                                                                                                  linewidth = 1,
                                                                                                                  width = 0.8,
                                                                                                                  text = list(size = ggplot2::rel(3))) +
                                                              ggplot2::scale_y_continuous(expand = c(0,
                                                                                                     0),
                                                                                          limits = c(0,
                                                                                                     common_limit)) +
                                                              ggplot2::scale_fill_identity() +
                                                              ggplot2::ylab("Number of Differentially\nExpressed Genes in Common")),
                                    matrix = ComplexUpset::intersection_matrix(geom = ggplot2::geom_point(size = 1.5)),
                                    set_sizes = ComplexUpset::upset_set_size(geom = ggplot2::geom_bar(linewidth = 1,
                                                                                                      width = 0.75)) +
                                      ggplot2::scale_y_reverse(expand = c(0,
                                                                          0),
                                                               limits = c(total_limit,
                                                                          0)) +
                                      ggplot2::ylab("Total Number of\nDifferentially Expressed Genes"),
                                    stripes = ComplexUpset::upset_stripes(geom = ggplot2::geom_segment(alpha = 0.5,
                                                                                                       linewidth = stripe_width),
                                                                          mapping = ggplot2::aes(colour = tissue),
                                                                          colors = c("BL" = pal[28],
                                                                                     "LI" = pal[17],
                                                                                     "LN" = pal[9],
                                                                                     "SP" = pal[21]),
                                                                          data = upset_data),
                                    themes = list('Intersection size' = ggplot2::theme_light() +
                                                    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                                                                   axis.ticks.x = ggplot2::element_blank(),
                                                                   axis.title.x = ggplot2::element_blank(),
                                                                   axis.title.y = ggplot2::element_text(margin = ggplot2::margin(l = -40))),
                                                  intersections_matrix = ggplot2::theme_light() +
                                                    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                                                                   axis.ticks = ggplot2::element_blank(),
                                                                   axis.title.y = ggplot2::element_blank(),
                                                                   axis.text.y = ggplot2::element_text(family = "mono",
                                                                                                       margin = ggplot2::margin(l = -40)),
                                                                   legend.position = "none"),
                                                  overall_sizes = ggplot2::theme_light() +
                                                    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                                                                   axis.ticks.y = ggplot2::element_blank(),
                                                                   axis.title.y = ggplot2::element_blank(),
                                                                   plot.margin = ggplot2::margin(r = 8))))
  # draw legend plot
  legend_plot <- ggplot2::ggplot(upset_data,
                                 ggplot2::aes(x = set,
                                              y = dpi,
                                              colour = tissue,
                                              fill = as.numeric(dpi))) +
    ggplot2::geom_col(linewidth = 1) +
    ggplot2::scale_colour_manual(name = "Tissue",
                                 values = cols) +
    ggplot2::scale_fill_viridis_c(begin = 1-((1/41)*36),
                                  direction = -1,
                                  name = "DPI") +
    ggplot2::guides(colour = ggplot2::guide_legend(order = 2,
                                                   override.aes = list(fill = ggplot2::alpha(cols,
                                                                                             0.5),
                                                                       stroke = 2),),
                    fill = ggplot2::guide_colourbar(draw.ulim = F,
                                                    order = 1)) +
    ggplot2::theme_light()
  # set layout  
  layout <- "AB
           CE"
  # combine plots
  plot <- upset_plot +
    legend_plot +
    patchwork::plot_layout(design = layout,
                           guides = "collect") &
    ggplot2::theme(legend.box = "horizontal")
  #save plot
  png(paste0("04_differential-expression-analysis/02_figures/02_upset-plots/",
             contrast_type,
             "-upset.png"),
      height = 6,
      res = 300,
      units = "in",
      width = 9)
  print(plot)
  dev.off()
  pdf(paste0("04_differential-expression-analysis/02_figures/02_upset-plots/",
             contrast_type,
             "-upset.pdf"),
      height = 6,
      width = 9)
  print(plot)
  dev.off()
}

#### apply function
upset_plot("RESP",
           total_limit = 800,
           common_limit = 400)

upset_plot("DIRE",
           total_limit = 2000,
           common_limit = 650,
           stripe_width = 3.75)

upset_plot("NDAM",
           total_limit = 5000,
           common_limit = 750)

upset_plot("BORA",
           total_limit = 5500,
           common_limit = 500,
           stripe_width = 4.25)

### volcano plots
#### function to draw and save volcano plots
volcano_plot <- function(contrast){
  # function to generate colour palette
  smooth_rainbow <- khroma::colour("smooth rainbow")
  # generate colour palettes 
  pal <- smooth_rainbow(31,
                        range = c(0.1,
                                  0.9))
  # set colour palettes
  cols <-  c(pal[26],
             "#666666",
             pal[5])
  # set labels
  labs <- c("Decreased Expression",
            "Not Significant",
            "Increased Expression")
  # select top 10 genes
  top_10_up <- de_genes %>%
    tidyr::drop_na(symbol) %>%
    dplyr::filter(!stringr::str_detect(symbol,
                                       "^ENSBTAG")) %>%
    dplyr::filter(!!rlang::sym(paste0("Results.",
                                      contrast)) == 1) %>%
    dplyr::arrange(!!rlang::sym(paste0("P.value.adj.",
                                       contrast))) %>%
    dplyr::distinct(symbol,
                    .keep_all = T) %>%
    dplyr::slice_min(n = 10,
                     order_by = !!rlang::sym(paste0("P.value.adj.",
                                                    contrast)))
  top_10_down <- de_genes %>%
    tidyr::drop_na(symbol) %>%
    dplyr::filter(!stringr::str_detect(symbol,
                                       "^ENSBTAG")) %>%
    dplyr::filter(!!rlang::sym(paste0("Results.",
                                      contrast)) == -1) %>%
    dplyr::arrange(!!rlang::sym(paste0("P.value.adj.",
                                       contrast))) %>%
    dplyr::distinct(symbol,
                    .keep_all = T) %>%
    dplyr::slice_min(n = 10,
                     order_by = !!rlang::sym(paste0("P.value.adj.",
                                                    contrast)))
  top_20 <- dplyr::full_join(top_10_up,
                             top_10_down)
  #set x max
  x_max <- de_genes %>%
    dplyr::slice_max(n = 1,
                     with_ties = F,
                     order_by = abs(!!rlang::sym(paste0("Coef.",
                                                        contrast)))) %>%
    dplyr::pull(!!rlang::sym(paste0("Coef.",
                                    contrast))) %>%
    abs()
  # draw plot
  plot <- de_genes %>%
    ggplot2::ggplot(ggplot2::aes(x = .data[[paste0("Coef.",
                                                   contrast)]],
                                 y = -log(.data[[paste0("P.value.adj.",
                                                        contrast)]],
                                          base = 10),
                                 colour = as.factor(.data[[paste0("Results.",
                                                                  contrast)]]),
                                 fill = as.factor(.data[[paste0("Results.",
                                                                contrast)]]))) +
    ggplot2::geom_point(shape = 21) +
    ggplot2::geom_hline(colour = "grey70",
                        linetype = "dashed",
                        yintercept = -log(0.05,
                                          base = 10)) +
    ggplot2::labs(x = "Log<sub>2</sub> Fold Change",
                  y = "-Log<sub>10</sub>*P*<sub>adj.</sub>") +
    ggplot2::scale_colour_manual(labels = labs,
                                 values = ggplot2::alpha(cols, c(1,
                                                                 0.1,
                                                                 1))) +
    ggplot2::scale_fill_manual(labels = labs,
                               values = ggplot2::alpha(cols, c(0.5,
                                                               0.05,
                                                               0.5))) +
    ggplot2::scale_x_continuous(limits = c(-x_max,
                                           x_max)) +
    ggplot2::theme_light() +
    ggplot2::theme(axis.title.x = ggtext::element_markdown(),
                   axis.title.y = ggtext::element_markdown(),
                   legend.position = "top",
                   legend.title = ggplot2::element_blank())
  # save plot
  png(paste0("04_differential-expression-analysis/02_figures/03_volcano-plots/",
             contrast,
             "-volcano.png"),
      height = 6,
      res = 300,
      units = "in",
      width = 9)
  print(plot)
  dev.off()
  pdf(paste0("04_differential-expression-analysis/02_figures/03_volcano-plots/",
             contrast,
             "-volcano.pdf"),
      height = 6,
      width = 9)
  print(plot)
  dev.off()
  # add labels 
  plot_labels <- plot +
    ggrepel::geom_text_repel(data = top_20,
                             ggplot2::aes(label = symbol),
                             bg.color = "white",
                             bg.r = 0.1,
                             box.padding = 0.5,
                             fontface = "italic",
                             force = 20,
                             force_pull = 0,
                             lineheight = 0.7,
                             max.iter = 1000000000,
                             max.overlaps = Inf,
                             max.time = 10,
                             min.segment.length = 0,
                             segment.size = 0.25,
                             show.legend = F,
                             size = ggplot2::rel(3))
  # save plot
  png(paste0("04_differential-expression-analysis/02_figures/03_volcano-plots/",
             contrast,
             "-volcano-labels.png"),
      height = 6,
      res = 300,
      units = "in",
      width = 9)
  print(plot_labels)
  dev.off()
  pdf(paste0("04_differential-expression-analysis/02_figures/03_volcano-plots/",
             contrast,
             "-volcano-labels.pdf"),
      height = 6,
      width = 9)
  print(plot_labels)
  dev.off()
}

#### apply function
lapply(contrast_ids,
       volcano_plot)

### functional enrichment
#### function to perform functional enrichment
functional_enrichment <- function(contrast){
  input <- de_genes %>%
    dplyr::select(probe,
                  AveExpr,
                  paste0("Coef.",
                         contrast),
                  paste0("t.",
                         contrast),
                  paste0("P.value.",
                         contrast),
                  paste0("P.value.adj.",
                         contrast),
                  F,
                  F.p.value,
                  F.p.value.adj,
                  paste0("Results.",
                         contrast),
                  ensembl,
                  symbol,
                  description) %>%
    readr::write_csv(paste0("05_functional-enrichment-analysis/01_input/",
                            stringr::str_replace_all(contrast,
                                                     "_",
                                                     "-"),
                            ".csv"))
  # set background
  background <- input %>% 
    dplyr::pull(ensembl)
  # set query genes
  a_up <- input %>%
    dplyr::filter(!!rlang::sym(paste0("Results.",
                                      contrast)) == 1) %>% 
    dplyr::select(ensembl)
  b_down <- input %>%
    dplyr::filter(!!rlang::sym(paste0("Results.",
                                      contrast)) == -1) %>% 
    dplyr::select(ensembl)
  # gprofiler functional analysis
  go <- gprofiler2::gost(correction_method = "gSCS",
                         custom_bg = background,
                         domain_scope = "custom_annotated",
                         evcodes = T,
                         highlight = T,
                         organism = "btaurus",
                         sources = c("GO"),
                         query = list("a_up" = a_up,
                                      "b_down" = b_down))
  # select top ten highlighted terms for each query
  try(top_terms <- subset(go[["result"]],
                          highlighted == T) %>%
        dplyr::group_by(query) %>% 
        dplyr::arrange(p_value) %>%
        dplyr::slice(1:10) %>%
        dplyr::mutate(top_term = T) %>%
        dplyr::ungroup())
  # add top term variable and save gprofiler results
  try(go[["result"]] <- dplyr::full_join(go[["result"]],
                                         top_terms) %>%
        dplyr::mutate(intersection_ratio = intersection_size/term_size,
                      plot_p_value = dplyr::case_when(-log10(p_value) > 16 ~ 16,
                                                      .default = -log10(p_value)),
                      source_label = factor(source,
                                            labels = c("GO:BP" = "'GO:Biological Process'",
                                                       "GO:CC" = "'GO:Cellular<br>Component'",
                                                       "GO:MF" = "'GO:Molecular Function'"),
                                            levels = c("GO:BP",
                                                       "GO:CC",
                                                       "GO:MF")),
                      term_name_wrap = stringr::str_wrap(term_name,
                                                         width = 10),
                      query_label = factor(query,
                                           labels = c("a_up" = "'Increased expression'",
                                                      "b_down" = "'Decreased expression'"),
                                           levels = c("a_up",
                                                      "b_down"))) %>%
        readr::write_csv(paste0("05_functional-enrichment-analysis/02_output/01_gprofiler/",
                                stringr::str_replace_all(contrast,
                                                         "_",
                                                         "-"),
                                "-gprofiler.csv")))
  # get numbers of terms
  bp_terms <- go[["meta"]][["result_metadata"]][["GO:BP"]][["number_of_terms"]]
  cc_terms <- go[["meta"]][["result_metadata"]][["GO:CC"]][["number_of_terms"]]
  mf_terms <- go[["meta"]][["result_metadata"]][["GO:MF"]][["number_of_terms"]]
  all_terms <- sum(bp_terms,
                   cc_terms,
                   mf_terms)
  # plot gprofiler results
  plot <- NULL
  try(plot <- ggplot2::ggplot(go[["result"]],
                              ggplot2::aes(x = source_order,
                                           y = plot_p_value)) +
        ggplot2::geom_point(ggplot2::aes(colour = query,
                                         fill = query,
                                         size = intersection_ratio),
                            shape = 21,
                            stroke = 0.5) +
        ggh4x::facet_grid2(ggplot2::vars(query_label),
                           ggplot2::vars(source_label),
                           drop = F,
                           labeller = ggplot2::label_parsed,
                           scales = "free_x",
                           strip = ggh4x::strip_vanilla(clip = "off"),
                           switch = "x") +
        ggplot2::geom_point(data = subset(go[["result"]],
                                          top_term == T),
                            ggplot2::aes(size = intersection_ratio),
                            shape = 21,
                            stroke = 0.5) +
        ggplot2::continuous_scale(aesthetics = c("size",
                                                 "point.size"), 
                                  breaks = c(0.25,
                                             0.5,
                                             0.75,
                                             1),
                                  labels = c(0.25,
                                             0.5,
                                             0.75,
                                             1),
                                  limits = c(0,
                                             1),
                                  scale_name = "size",
                                  palette = scales::area_pal()) +
        ggrepel::geom_text_repel(data = subset(go[["result"]],
                                               top_term == T),
                                 ggplot2::aes(label = term_name_wrap,
                                              point.size = intersection_ratio),
                                 bg.color = "white",
                                 bg.r = 0.1,
                                 box.padding = 0.5,
                                 force = 20,
                                 force_pull = 0,
                                 lineheight = 0.7,
                                 max.iter = 1000000000,
                                 max.overlaps = Inf,
                                 max.time = 10,
                                 min.segment.length = 0,
                                 segment.size = 0.25,
                                 size = ggplot2::rel(3)) +
        ggplot2::labs(size = "Intersection ratio",
                      y = "-Log<sub>10</sub>*P*<sub>adj.</sub>") +
        ggplot2::guides(colour = "none",
                        fill = "none",
                        point.size = "none",
                        size = ggplot2::guide_legend(nrow = 1,
                                                     order = 1)) +
        ggplot2::scale_colour_manual(values = c("a_up" = pal[5],
                                                "b_down" = pal[26])) +
        ggplot2::scale_fill_manual(values = ggplot2::alpha(c("a_up" = pal[5],
                                                             "b_down" = pal[26]),
                                                           0.8)) +
        ggh4x::scale_x_facet(COL == 1,
                             breaks = c(0,
                                        cc_terms,
                                        cc_terms*2,
                                        cc_terms*3,
                                        cc_terms*4,
                                        cc_terms*5,
                                        cc_terms*6,
                                        cc_terms*7),
                             minor_breaks = NULL,
                             limits = c(0,
                                        bp_terms)) +
        ggh4x::scale_x_facet(COL == 2,
                             minor_breaks = NULL,
                             breaks = c(0,
                                        cc_terms),
                             limits = c(0,
                                        cc_terms)) +
        ggh4x::scale_x_facet(COL == 3,
                             minor_breaks = NULL,
                             breaks = c(0,
                                        cc_terms,
                                        cc_terms*2),
                             limits = c(0,
                                        mf_terms)) +
        ggplot2::scale_y_continuous(breaks = c(0,
                                               4,
                                               8,
                                               12,
                                               16),
                                    expand = c(0,
                                               0),
                                    labels = c(0,
                                               4,
                                               8,
                                               12,
                                               expression("">=16)),
                                    limits = c(0,
                                               17)) +
        ggh4x::force_panelsizes(cols = c(bp_terms/all_terms,
                                         cc_terms/all_terms,
                                         mf_terms/all_terms),
                                respect = F) +
        ggplot2::theme_light() +
        ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                       axis.ticks.x = ggplot2::element_blank(),
                       axis.title.x = ggplot2::element_blank(),
                       axis.title.y = ggplot2::ggtext::element_markdown(),
                       legend.margin = ggplot2::margin(0,
                                                       0,
                                                       0,
                                                       0),
                       legend.position = "top",
                       strip.background = ggplot2::element_blank(),
                       strip.text = ggtext::element_markdown(colour = "black",
                                                             size = ggplot2::rel(1))))
  # save gprofiler plot
  png(paste0("05_functional-enrichment-analysis/03_figures/01_gprofiler/",
             stringr::str_replace_all(contrast,
                                      "_",
                                      "-"),
             "-gprofiler.png"),
      height = 6,
      res = 300,
      units = "in",
      width = 9)
  print(plot)
  dev.off()
  pdf(paste0("05_functional-enrichment-analysis/03_figures/01_gprofiler/",
             stringr::str_replace_all(contrast,
                                      "_",
                                      "-"),
             "-gprofiler.pdf"),
      height = 6,
      width = 9)
  print(plot)
  dev.off()
  # format to gem
  gem <- go$result[,c("query",
                      "term_id",
                      "term_name",
                      "p_value",
                      "intersection")]
  try(colnames(gem) <- c("query",
                         "GO.ID",
                         "Description",
                         "p.Val",
                         "Genes"))
  gem$FDR <- gem$p.Val
  gem$Phenotype = "+1"
  # write separate files for queries
  try(gem %>% dplyr::group_by(query) %>%
        dplyr::group_walk(~ write.table(data.frame(.x[,c("GO.ID",
                                                         "Description",
                                                         "p.Val",
                                                         "FDR",
                                                         "Phenotype",
                                                         "Genes")]),
                                        file = paste0("05_functional-enrichment-analysis/02_output/02_gem-files/01_all/",
                                                      stringr::str_replace_all(contrast,
                                                                               "_",
                                                                               "-"),
                                                      "-",
                                                      stringr::str_replace_all(unique(.y$query),
                                                                               "_",
                                                                               "-"),
                                                      "-gem.txt"),
                                        sep = "\t",
                                        quote = F,
                                        row.names = F)))
}

#### apply function
lapply(contrast_ids,
       functional_enrichment)

#### remove empty pdf files
file.info(dir("05_functional-enrichment-analysis/03_figures/01_gprofiler")) %>%
  tibble::rownames_to_column(var = "name") %>%
  dplyr::transmute(name = stringr::str_split_i(name,
                                               "\\.",
                                               1)) %>%
  dplyr::group_by(name) %>%
  dplyr::filter(n()==1) %>%
  dplyr::transmute(name = stringr::str_c("05_functional-enrichment-analysis/03_figures/01_gprofiler/",
                                         name,
                                         ".pdf")) %>%
  dplyr::pull() %>%
  file.remove()

#### combine gem files
combine_gem_files <- function(contrast_type, tissue, direction){
  list.files("05_functional-enrichment-analysis/02_output/02_gem-files/01_all/",
             full.names = T,
             pattern = paste0(contrast_type,
                              "-",
                              tissue,
                              ".*",
                              dplyr::case_match(direction,
                                                "up" ~ "a-up",
                                                "down" ~ "b-down"))) %>%
    purrr::map_dfr(readr::read_delim) %>% 
    readr::write_delim(paste0("05_functional-enrichment-analysis/02_output/02_gem-files/02_combined/",
                              contrast_type,
                              "-",
                              direction,
                              "/",
                              contrast_type,
                              "-",
                              tissue,
                              "-",
                              direction,
                              "-gem.txt"),
                       delim = "\t")
}

contrast_types <- c("RESP",
                    "DIRE",
                    "NDAM",
                    "BORA")

#### apply function
lapply(contrast_types,
       combine_gem_files,
       tissue = "BL",
       direction = "up")
lapply(contrast_types,
       combine_gem_files,
       tissue = "BL",
       direction = "down")

lapply(contrast_types,
       combine_gem_files,
       tissue = "LI",
       direction = "up")
lapply(contrast_types,
       combine_gem_files,
       tissue = "LI",
       direction = "down")

lapply(contrast_types,
       combine_gem_files,
       tissue = "LN",
       direction = "up")
lapply(contrast_types,
       combine_gem_files,
       tissue = "LN",
       direction = "down")

lapply(contrast_types,
       combine_gem_files,
       tissue = "SP",
       direction = "up")
lapply(contrast_types,
       combine_gem_files,
       tissue = "SP",
       direction = "down")

#### draw enrichmentmap legend plot
enrichmentmap_legend_plot <- de_numbers %>%
  dplyr::filter(de != "NotSig") %>%
  ggplot2::ggplot(ggplot2::aes(x = contrast_id,
                               y = dpi,
                               colour = de,
                               fill = tissue)) +
  ggplot2::geom_point(size = 10,
                      shape = 21) +
  ggplot2::scale_colour_manual(name = "Expression",
                               labels = c("Increased",
                                          "Decreased"),
                               limits = c("up",
                                          "down"),
                               values = c(pal[5],
                                          pal[26])) +
  ggplot2::scale_fill_manual(name = "Tissue",
                             values = cols) +
  ggplot2::guides(colour = ggplot2::guide_legend(order = 1,
                                                 override.aes = list(fill = ggplot2::alpha(c(pal[5],
                                                                                             pal[26]),
                                                                                           0.25))),
                  fill = guide_legend(order = 2,
                                      override.aes = list(colour = cols))) +
  ggplot2::theme_light() +
  ggplot2::theme(legend.position = "top",
                 legend.text = ggplot2::element_text(size = 20),
                 legend.title = ggplot2::element_text(size = 24))

#### save plot
pdf("05_functional-enrichment-analysis/03_figures/02_enrichmentmap/legend.pdf",
    height = 6,
    width = 12)
print(enrichmentmap_legend_plot)
dev.off()