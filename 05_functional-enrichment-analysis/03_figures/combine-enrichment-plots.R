# combine enrichment plots
## combine results of enrichment analyses of microarray data
### install required packages
if (!require(dplyr)) install.packages("dplyr")
if (!require(khroma)) install.packages("khroma")
if (!require(magick)) install.packages("magick")
if (!require(magrittr)) install.packages("magrittr")
if (!require(stringr)) install.packages("stringr")
if (!require(viridis)) install.packages("viridis")

### function to combine gprofiler enrichment plots
combine_gprofiler_plots <- function(tissue, contrast_type){
  # get contrasts
  contrasts <- list.files("05_functional-enrichment-analysis/03_figures/01_gprofiler/",
                          pattern = paste0(contrast_type,
                                           "-",
                                           tissue,
                                           ".*",
                                           "-gprofiler.png"))
  contrast_info <- as.data.frame(contrasts) %>%
    dplyr::mutate(dpi = stringr::str_split_i(contrasts,
                                             "-",
                                             3),
                  dpi = dplyr::case_when(dpi == "00" ~ "0",
                                         .default = as.character(dpi)))
  # function to generate colour palette
  smooth_rainbow <- khroma::colour("smooth rainbow")
  # generate colour palette
  pal <- smooth_rainbow(31,
                        range = c(0.1,
                                  0.9))
  viridis <- viridis::viridis(41,
                              direction = -1)
  # read, crop and add space to legend
  legend <- magick::image_read(paste0("05_functional-enrichment-analysis/03_figures/01_gprofiler/",
                                      contrasts[1])) %>%
    magick::image_crop("2700x120")
  legend_space <- magick::image_append(c(magick::image_blank(20,
                                                             120,
                                                             color = "white"),
                                         legend))
  # add legend to stack of plots
  plots <- c(legend_space)
  # read, crop and annotate plots
  for (i in rev(1:length(contrasts))){
    assign(paste0(contrasts[i],
                  "_plot"),
           magick::image_read(paste0("05_functional-enrichment-analysis/03_figures/01_gprofiler/",
                                     contrasts[i])) %>%
             magick::image_crop("2620x1535+80+120"))
    assign(paste0(contrasts[i],
                  "_plot_label"),
           magick::image_append(c(magick::image_blank(100,
                                                      1535,
                                                      color = "white"),
                                  get(paste0(contrasts[i],
                                             "_plot")))) %>%
             magick::image_annotate(paste0(contrast_info$dpi[i],
                                           " Days Post-Infection"),
                                    degrees = 270,
                                    gravity = "west",
                                    location = "+45+323",
                                    size = 60) %>%
             magick::image_composite(magick::image_blank(5,
                                                         1450,
                                                         color = viridis[[as.numeric(contrast_info$dpi[i]) + 1]]),
                                     gravity = "west",
                                     offset = "+100+23") %>%
             magick::image_composite(magick::image_blank(5,
                                                         739,
                                                         color = pal[5]),
                                     gravity = "east",
                                     offset = "+110-419") %>%
             magick::image_composite(magick::image_blank(5,
                                                         739,
                                                         color = pal[26]),
                                     gravity = "east",
                                     offset = "+110+385"))
    plots <- magick::image_append(c(plots,
                                    get(paste0(contrasts[i],
                                               "_plot_label"))),
                                  stack = T)
  }
  # read, crop and add space to x-axis title
  x <- magick::image_read(paste0("05_functional-enrichment-analysis/03_figures/01_gprofiler/",
                                 contrasts[1])) %>%
    magick::image_crop("2700x145+0+1655")
  x_space <- magick::image_append(c(magick::image_blank(20,
                                                        120,
                                                        color = "white"),
                                    x))
  # add x-axis title
  plots_x <- magick::image_append(c(plots,
                                    x_space),
                                  stack = T)
  # read, crop and add space to y-axis title
  y <- magick::image_read(paste0("05_functional-enrichment-analysis/03_figures/01_gprofiler/",
                                 contrasts[1])) %>%
    magick::image_crop("80x1800")
  dim <- magick::image_info(plots_x)
  y_space <- magick::image_composite(magick::image_blank(80,
                                                         dim$height,
                                                         color = "white"),
                                     y,
                                     gravity = "west")
  # add y-axis title
  plot <- magick::image_append(c(y_space,
                                 plots_x))
  # save plot
  magick::image_write(plot,
                      paste0("05_functional-enrichment-analysis/03_figures/01_gprofiler/",
                             contrast_type,
                             "-",
                             tissue,
                             "-gprofiler.png"))
}

### set contrasts
DIRE_BL_contrasts <- c("DIRE_BL_34")
RESP_BL_contrasts <- c("RESP_BL_14",
                       "RESP_BL_25",
                       "RESP_BL_34")
NDAM_BL_contrasts <- c("NDAM_BL_14",
                       "NDAM_BL_25",
                       "NDAM_BL_34")
BORA_BL_contrasts <- c("BORA_BL_14",
                       "BORA_BL_25",
                       "BORA_BL_34")
DIRE_LI_contrasts <- c("DIRE_LI_00",
                       "DIRE_LI_12",
                       "DIRE_LI_21",
                       "DIRE_LI_26",
                       "DIRE_LI_29",
                       "DIRE_LI_32",
                       "DIRE_LI_35")
RESP_LI_contrasts <- c("RESP_LI_26",
                       "RESP_LI_32",
                       "RESP_LI_35")
NDAM_LI_contrasts <- c("NDAM_LI_12",
                       "NDAM_LI_15",
                       "NDAM_LI_18",
                       "NDAM_LI_21",
                       "NDAM_LI_26",
                       "NDAM_LI_29",
                       "NDAM_LI_32",
                       "NDAM_LI_35")
BORA_LI_contrasts <- c("BORA_LI_12",
                       "BORA_LI_15",
                       "BORA_LI_18",
                       "BORA_LI_21",
                       "BORA_LI_26",
                       "BORA_LI_29",
                       "BORA_LI_32",
                       "BORA_LI_35")
DIRE_LN_contrasts <- c("DIRE_LN_00",
                       "DIRE_LN_35")
RESP_LN_contrasts <- c("RESP_LN_21",
                       "RESP_LN_35")
NDAM_LN_contrasts <- c("NDAM_LN_21",
                       "NDAM_LN_35")
BORA_LN_contrasts <- c("BORA_LN_21",
                       "BORA_LN_35")
DIRE_SP_contrasts <- c("DIRE_SP_00")
RESP_SP_contrasts <- c("RESP_SP_21",
                       "RESP_SP_35")
NDAM_SP_contrasts <- c("NDAM_SP_21",
                       "NDAM_SP_35")
BORA_SP_contrasts <- c("BORA_SP_21",
                       "BORA_SP_35")

### set tissues
tissues <- c("BL",
             "LI",
             "LN",
             "SP")

### apply function
mapply(combine_gprofiler_plots,
       tissues,
       contrast_type = "DIRE")

mapply(combine_gprofiler_plots,
       tissues,
       contrast_type = "RESP")

mapply(combine_gprofiler_plots,
       tissues,
       contrast_type = "NDAM")

mapply(combine_gprofiler_plots,
       tissues,
       contrast_type = "BORA")

### function to combine enrichmentmap plots
combine_enrichmentmap_plots <- function(contrast_type){
  # read and crop legend
  legend <- magick::image_read_pdf("05_functional-enrichment-analysis/03_figures/02_enrichmentmap/legend.pdf") %>%
    magick::image_crop("3508x200+100")
  # read plots
  up <- magick::image_read_pdf(paste0("05_functional-enrichment-analysis/03_figures/02_enrichmentmap/",
                                      contrast_type,
                                      "-up.pdf")) %>%
    magick::image_crop("3508x2250")
  down <- magick::image_read_pdf(paste0("05_functional-enrichment-analysis/03_figures/02_enrichmentmap/",
                                        contrast_type,
                                        "-down.pdf")) %>%
    magick::image_crop("3508x2250")
  # add x-axis title
  plot <- magick::image_append(c(legend,
                                 up,
                                 down),
                               stack = T)
  # save plot
  magick::image_write(plot,
                      paste0("05_functional-enrichment-analysis/03_figures/02_enrichmentmap/",
                             contrast_type,
                             "-enrichmentmap.png"))
}

### set contrast types
contrast_types <- c("RESP",
                    "DIRE",
                    "NDAM",
                    "BORA")

### apply function
mapply(combine_enrichmentmap_plots,
       contrast_types)