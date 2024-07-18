# functional enrichment analysis
## functional enrichment analysis using enrichmentmap
### download gmt file from the data sources section of gprofiler https://biit.cs.ut.ee/gprofiler
wget https://biit.cs.ut.ee/gprofiler//static/gprofiler_full_btaurus.ENSG.gmt -O 01_data-sources/04_gprofiler/gprofiler_full_btaurus.ENSG.gmt

### open cytoscape
### click apps > enrichmentmap > scan a folder for enrichment data
### select a folder in 05_functional-enrichment-analysis\02_output\02_gem-files\02_combined and click open
### click common files included in all data sets
### click the 3 dots next to the blank gmt file line
### select gprofiler_full_btaurus.ENSG.gmt and click open
### click build
### click the enrichmentmap tab, in the style section use the dropdown menu to set the chart data to colour by data set
### click the options button beside the data sets heading in the filter section
### click change colours
### click on the colour assigned to the edges
### click rgb and set the colour code to 666666
### repeat the process for each data set gem file, setting 
#### blood results to D62221
#### liver results to 80B973
#### lymph node results to 4F74C1
#### spleen results to D4B13F
### click the style tab and set node border width to 0.0, node image/chart border width to 0.0 in options, node colour transparency to 0, and edge transparency set to 50
### click the network tab, left click the network and click rename network, edit the network name to match the folder of gmt files and click ok
### click apps > autoannotate > new annotation set
### tick the box to layout network to prevent cluster overlap and click create annotations
### in the autoannotate display set the fill and border colour t0
#### 8E509A for up and
#### E35E2C for down
### set the min font size to 12 and font scale to 40%
### tick the box to activate word wrap and set the wrap length to 15
### click layout > yfiles remove overlaps
### click the autoannotate tab and click on each cluster name to highlight it and move it as necessary
### left click the cluster to rename if required
### click file > export > network to image

### combine enrichment plots
Rscript 05_functional-enrichment-analysis/03_figures/combine-enrichment-plots.R