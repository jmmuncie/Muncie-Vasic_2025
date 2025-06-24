###----------------Description-------------------------------------

#This script is written to create tables for all Mef2c-dependent DARs in the
#cell types of interest


###---------------Load Libraries and Set Threads------------------------------

#To install ArchR release_1.0.2 branch:
devtools::install_github("GreenleafLab/ArchR", ref="release_1.0.2", repos = BiocManager::repositories())

library(ArchR); library(Seurat); library(dplyr)
library(GenomeInfoDb); library(ensembldb)
library(ggplot2); library(patchwork)
set.seed(1234)
library(BSgenome) 
library(hexbin)

#Set ArchR Threads to 4:
addArchRThreads(threads = 4, force = FALSE)


###---------------E7.75 DARs--------------------------------------------

#Load ArchR proj and set working directory
proj_Mef2c_v13_E775_subset <- loadArchRProject("~/Mef2c_ArchR_working/Mef2c_v13_E775_subset")
setwd("~/Mef2c_ArchR_working")

markerTest_CPs <- getMarkerFeatures(
  ArchRProj = proj_Mef2c_v13_E775_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "CellTypeByGenotype",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "CPs_x_KO",
  bgdGroups = "CPs_x_WT"
)

markerTest_CPs_list_upinWT <- getMarkers(markerTest_CPs, cutOff = "FDR <= 0.15 & Log2FC <= -0.5", returnGR = TRUE)
markerTest_CPs_list_upinKO <- getMarkers(markerTest_CPs, cutOff = "FDR <= 0.15 & Log2FC >= 0.5", returnGR = TRUE)

#Let's convert the upinWT list to a GRanges object, sorted by Log2FC
markerpeaksgr_E775_CPs_WT <- markerTest_CPs_list_upinWT$`CPs_x_KO`
order_E775_CPs_WT <- order(markerpeaksgr_E775_CPs_WT$Log2FC, start(markerpeaksgr_E775_CPs_WT))
markerpeaks_sorted_E775_CPs_WT <- markerpeaksgr_E775_CPs_WT[order_E775_CPs_WT]

#Let's convert the upinKO list to a GRanges object, sorted by Log2FC
markerpeaksgr_E775_CPs_KO <- markerTest_CPs_list_upinKO$`CPs_x_KO`
order_E775_CPs_KO <- order(markerpeaksgr_E775_CPs_KO$Log2FC, start(markerpeaksgr_E775_CPs_KO))
markerpeaks_sorted_E775_CPs_KO <- markerpeaksgr_E775_CPs_KO[order_E775_CPs_KO]

#Intersect these peaks with the full peak set to get nearestGene and peakType info
#Subset matching peaks
pgr_E775 <- getPeakSet(proj_Mef2c_v13_E775_subset_CPspool)
sub_pgr_E775_upWT_idx <- findOverlaps(markerpeaks_sorted_E775_CPs_WT, pgr_E775)
sub_pgr_E775_upWT <- pgr_E775[sub_pgr_E775_upWT_idx@to]
sub_pgr_E775_upKO_idx <- findOverlaps(markerpeaks_sorted_E775_CPs_KO, pgr_E775)
sub_pgr_E775_upKO <- pgr_E775[sub_pgr_E775_upKO_idx@to]

#Create data frame for all WT DARs 
WT_DARs_E775_CPs_df <- data.frame(
  mm10PeakRange = paste0(markerpeaks_sorted_E775_CPs_WT@seqnames,':',markerpeaks_sorted_E775_CPs_WT@ranges@start,'-',markerpeaks_sorted_E775_CPs_WT@ranges@start + 500),
  peakType = sub_pgr_E775_upWT$peakType,
  peakLog2FC = markerpeaks_sorted_E775_CPs_WT$Log2FC,
  peakFDR = markerpeaks_sorted_E775_CPs_WT$FDR,
  nearestTSS =  sub_pgr_E775_upWT$nearestTSS, 
  distToTSS = sub_pgr_E775_upWT$distToTSS
)

#Create data frame for all KO DARs 
KO_DARs_E775_CPs_df <- data.frame(
  mm10PeakRange = paste0(markerpeaks_sorted_E775_CPs_KO@seqnames,':',markerpeaks_sorted_E775_CPs_KO@ranges@start,'-',markerpeaks_sorted_E775_CPs_KO@ranges@start + 500),
  peakType = sub_pgr_E775_upKO$peakType,
  peakLog2FC = markerpeaks_sorted_E775_CPs_KO$Log2FC,
  peakFDR = markerpeaks_sorted_E775_CPs_KO$FDR,
  nearestTSS =  sub_pgr_E775_upKO$nearestTSS, 
  distToTSS = sub_pgr_E775_upKO$distToTSS
)

#Combine dataframes
E775_CPs_DARs <- rbind(WT_DARs_E775_CPs_df, KO_DARs_E775_CPs_df)

# Save dataframe to a CSV file
write.csv(E775_CPs_DARs, "~/Desktop/Working_Directory/E775_CPs_DARs.csv", row.names = FALSE)


###---------------E85 DARs--------------------------------------------

#Load ArchR proj and set working directory
proj_Mef2c_v13_E85_subset <- loadArchRProject("~/Mef2c_ArchR_working/Mef2c_v13_E85_subset")
setwd("~/Mef2c_ArchR_working")


#CMs_IFT
markerTest_E85_CMs_IFT <- getMarkerFeatures(
  ArchRProj = proj_Mef2c_v13_E85_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "CellTypeByGenotype",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "CMs_IFT_x_KO",
  bgdGroups = "CMs_IFT_x_WT"
)

markerTest_E85_CMs_IFT_list_upinWT <- getMarkers(markerTest_E85_CMs_IFT, cutOff = "FDR <= 0.15 & Log2FC <= -0.5", returnGR = TRUE)
markerTest_E85_CMs_IFT_list_upinKO <- getMarkers(markerTest_E85_CMs_IFT, cutOff = "FDR <= 0.15 & Log2FC >= 0.5", returnGR = TRUE)

#Let's convert the upinWT list to a GRanges object, sorted by Log2FC
markerpeaksgr_E85_IFT_WT <- markerTest_E85_CMs_IFT_list_upinWT$`CMs_IFT_x_KO`
order_E85_IFT_WT <- order(markerpeaksgr_E85_IFT_WT$Log2FC, start(markerpeaksgr_E85_IFT_WT))
markerpeaks_sorted_E85_IFT_WT <- markerpeaksgr_E85_IFT_WT[order_E85_IFT_WT]

#Let's convert the upinKO list to a GRanges object, sorted by Log2FC
markerpeaksgr_E85_IFT_KO <- markerTest_E85_CMs_IFT_list_upinKO$`CMs_IFT_x_KO`
order_E85_IFT_KO <- order(markerpeaksgr_E85_IFT_KO$Log2FC, start(markerpeaksgr_E85_IFT_KO))
markerpeaks_sorted_E85_IFT_KO <- markerpeaksgr_E85_IFT_KO[order_E85_IFT_KO]

#Intersect these peaks with the full peak set to get nearestGene and peakType info
#Subset matching peaks
pgr_E85 <- getPeakSet(proj_Mef2c_v13_E85_subset)
sub_pgr_E85_IFT_upWT_idx <- findOverlaps(markerpeaks_sorted_E85_IFT_WT, pgr_E85)
sub_pgr_E85_IFT_upWT <- pgr_E85[sub_pgr_E85_IFT_upWT_idx@to]
sub_pgr_E85_IFT_upKO_idx <- findOverlaps(markerpeaks_sorted_E85_IFT_KO, pgr_E85)
sub_pgr_E85_IFT_upKO <- pgr_E85[sub_pgr_E85_IFT_upKO_idx@to]

#Create data frame for all WT DARs 
WT_DARs_E85_CMs_IFT_df <- data.frame(
  mm10PeakRange = paste0(markerpeaks_sorted_E85_IFT_WT@seqnames,':',markerpeaks_sorted_E85_IFT_WT@ranges@start,'-',markerpeaks_sorted_E85_IFT_WT@ranges@start + 500),
  peakType = sub_pgr_E85_IFT_upWT$peakType,
  peakLog2FC = markerpeaks_sorted_E85_IFT_WT$Log2FC,
  peakFDR = markerpeaks_sorted_E85_IFT_WT$FDR,
  nearestTSS =  sub_pgr_E85_IFT_upWT$nearestTSS, 
  distToTSS = sub_pgr_E85_IFT_upWT$distToTSS
)

#Create data frame for all KO DARs 
KO_DARs_E85_CMs_IFT_df <- data.frame(
  mm10PeakRange = paste0(markerpeaks_sorted_E85_IFT_KO@seqnames,':',markerpeaks_sorted_E85_IFT_KO@ranges@start,'-',markerpeaks_sorted_E85_IFT_KO@ranges@start + 500),
  peakType = sub_pgr_E85_IFT_upKO$peakType,
  peakLog2FC = markerpeaks_sorted_E85_IFT_KO$Log2FC,
  peakFDR = markerpeaks_sorted_E85_IFT_KO$FDR,
  nearestTSS =  sub_pgr_E85_IFT_upKO$nearestTSS, 
  distToTSS = sub_pgr_E85_IFT_upKO$distToTSS
)

#Combine dataframes
E85_CMs_IFT_DARs <- rbind(WT_DARs_E85_CMs_IFT_df, KO_DARs_E85_CMs_IFT_df)

# Save dataframe to a CSV file
write.csv(E85_CMs_IFT_DARs, "~/Desktop/Working_Directory/E85_CMs_IFT_DARs.csv", row.names = FALSE)


#CMs_V
markerTest_E85_CMs_V <- getMarkerFeatures(
  ArchRProj = proj_Mef2c_v13_E85_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "CellTypeByGenotype",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "CMs_V_x_KO",
  bgdGroups = "CMs_V_x_WT"
)

markerTest_E85_CMs_V_list_upinWT <- getMarkers(markerTest_E85_CMs_V, cutOff = "FDR <= 0.15 & Log2FC <= -0.5", returnGR = TRUE)
markerTest_E85_CMs_V_list_upinKO <- getMarkers(markerTest_E85_CMs_V, cutOff = "FDR <= 0.15 & Log2FC >= 0.5", returnGR = TRUE)

#Let's convert the upinWT list to a GRanges object, sorted by Log2FC
markerpeaksgr_E85_V_WT <- markerTest_E85_CMs_V_list_upinWT$`CMs_V_x_KO`
order_E85_V_WT <- order(markerpeaksgr_E85_V_WT$Log2FC, start(markerpeaksgr_E85_V_WT))
markerpeaks_sorted_E85_V_WT <- markerpeaksgr_E85_V_WT[order_E85_V_WT]

#Let's convert the upinKO list to a GRanges object, sorted by Log2FC
markerpeaksgr_E85_V_KO <- markerTest_E85_CMs_V_list_upinKO$`CMs_V_x_KO`
order_E85_V_KO <- order(markerpeaksgr_E85_V_KO$Log2FC, start(markerpeaksgr_E85_V_KO))
markerpeaks_sorted_E85_V_KO <- markerpeaksgr_E85_V_KO[order_E85_V_KO]

#Intersect these peaks with the full peak set to get nearestGene and peakType info
#Subset matching peaks
pgr_E85 <- getPeakSet(proj_Mef2c_v13_E85_subset)
sub_pgr_E85_V_upWT_idx <- findOverlaps(markerpeaks_sorted_E85_V_WT, pgr_E85)
sub_pgr_E85_V_upWT <- pgr_E85[sub_pgr_E85_V_upWT_idx@to]
sub_pgr_E85_V_upKO_idx <- findOverlaps(markerpeaks_sorted_E85_V_KO, pgr_E85)
sub_pgr_E85_V_upKO <- pgr_E85[sub_pgr_E85_V_upKO_idx@to]

#Create data frame for all WT DARs 
WT_DARs_E85_CMs_V_df <- data.frame(
  mm10PeakRange = paste0(markerpeaks_sorted_E85_V_WT@seqnames,':',markerpeaks_sorted_E85_V_WT@ranges@start,'-',markerpeaks_sorted_E85_V_WT@ranges@start + 500),
  peakType = sub_pgr_E85_V_upWT$peakType,
  peakLog2FC = markerpeaks_sorted_E85_V_WT$Log2FC,
  peakFDR = markerpeaks_sorted_E85_V_WT$FDR,
  nearestTSS =  sub_pgr_E85_V_upWT$nearestTSS, 
  distToTSS = sub_pgr_E85_V_upWT$distToTSS
)

#Create data frame for all KO DARs 
KO_DARs_E85_CMs_V_df <- data.frame(
  mm10PeakRange = paste0(markerpeaks_sorted_E85_V_KO@seqnames,':',markerpeaks_sorted_E85_V_KO@ranges@start,'-',markerpeaks_sorted_E85_V_KO@ranges@start + 500),
  peakType = sub_pgr_E85_V_upKO$peakType,
  peakLog2FC = markerpeaks_sorted_E85_V_KO$Log2FC,
  peakFDR = markerpeaks_sorted_E85_V_KO$FDR,
  nearestTSS =  sub_pgr_E85_V_upKO$nearestTSS, 
  distToTSS = sub_pgr_E85_V_upKO$distToTSS
)

#Combine dataframes
E85_CMs_V_DARs <- rbind(WT_DARs_E85_CMs_V_df, KO_DARs_E85_CMs_V_df)

# Save dataframe to a CSV file
write.csv(E85_CMs_V_DARs, "~/Desktop/Working_Directory/E85_CMs_V_DARs.csv", row.names = FALSE)


#CMs_OFT
#Note: due to low cell numbers, OFT CMs from E85 and E9 were combined for DAR
#testing, leading to a single set of DARs for these two timepoints

#Load combined OFT project
proj_Mef2c_v13_E85_E9_subset <- loadArchRProject("~/Mef2c_ArchR_working/Mef2c_v13_E85_E9_subset")
setwd("~/Mef2c_ArchR_working")

#MakerTest
markerTest_CMs_OFT <- getMarkerFeatures(
  ArchRProj = proj_Mef2c_v13_E85_E9_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "CellTypeByGenotype",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "CMs_OFT_x_KO",
  bgdGroups = "CMs_OFT_x_WT"
)

markerTest_CMs_OFT_list_upinWT <- getMarkers(markerTest_CMs_OFT, cutOff = "FDR <= 0.15 & Log2FC <= -0.5", returnGR = TRUE)
markerTest_CMs_OFT_list_upinKO <- getMarkers(markerTest_CMs_OFT, cutOff = "FDR <= 0.15 & Log2FC >= 0.5", returnGR = TRUE)

#Let's convert the upinWT list to a GRanges object, sorted by Log2FC
markerpeaksgr_OFT_WT <- markerTest_CMs_OFT_list_upinWT$`CMs_OFT_x_KO`
order_OFT_WT <- order(markerpeaksgr_OFT_WT$Log2FC, start(markerpeaksgr_OFT_WT))
markerpeaks_sorted_OFT_WT <- markerpeaksgr_OFT_WT[order_OFT_WT]

#Let's convert the upinKO list to a GRanges object, sorted by Log2FC
markerpeaksgr_OFT_KO <- markerTest_CMs_OFT_list_upinKO$`CMs_OFT_x_KO`
order_OFT_KO <- order(markerpeaksgr_OFT_KO$Log2FC, start(markerpeaksgr_OFT_KO))
markerpeaks_sorted_OFT_KO <- markerpeaksgr_OFT_KO[order_OFT_KO]

#Intersect these peaks with the full peak set to get nearestGene and peakType info
#Subset matching peaks
pgr_OFT <- getPeakSet(proj_Mef2c_v13_E85_E9_subset)
sub_pgr_OFT_upWT_idx <- findOverlaps(markerpeaks_sorted_OFT_WT, pgr_OFT)
sub_pgr_OFT_upWT <- pgr_OFT[sub_pgr_OFT_upWT_idx@to]
sub_pgr_OFT_upKO_idx <- findOverlaps(markerpeaks_sorted_OFT_KO, pgr_OFT)
sub_pgr_OFT_upKO <- pgr_OFT[sub_pgr_OFT_upKO_idx@to]

#Create data frame for all WT DARs 
WT_DARs_OFT_df <- data.frame(
  mm10PeakRange = paste0(markerpeaks_sorted_OFT_WT@seqnames,':',markerpeaks_sorted_OFT_WT@ranges@start,'-',markerpeaks_sorted_OFT_WT@ranges@start + 500),
  peakType = sub_pgr_OFT_upWT$peakType,
  peakLog2FC = markerpeaks_sorted_OFT_WT$Log2FC,
  peakFDR = markerpeaks_sorted_OFT_WT$FDR,
  nearestTSS =  sub_pgr_OFT_upWT$nearestTSS, 
  distToTSS = sub_pgr_OFT_upWT$distToTSS
)

#Create data frame for all KO DARs 
KO_DARs_OFT_df <- data.frame(
  mm10PeakRange = paste0(markerpeaks_sorted_OFT_KO@seqnames,':',markerpeaks_sorted_OFT_KO@ranges@start,'-',markerpeaks_sorted_OFT_KO@ranges@start + 500),
  peakType = sub_pgr_OFT_upKO$peakType,
  peakLog2FC = markerpeaks_sorted_OFT_KO$Log2FC,
  peakFDR = markerpeaks_sorted_OFT_KO$FDR,
  nearestTSS =  sub_pgr_OFT_upKO$nearestTSS, 
  distToTSS = sub_pgr_OFT_upKO$distToTSS
)

#Combine dataframes
E85andE9_CMs_OFT_DARs <- rbind(WT_DARs_OFT_df, KO_DARs_OFT_df)

# Save dataframe to a CSV file
write.csv(E85andE9_CMs_OFT_DARs, "~/Desktop/Working_Directory/E85andE9_CMs_OFT_DARs.csv", row.names = FALSE)


###---------------E9 DARs--------------------------------------------

#Load ArchR proj and set working directory
proj_Mef2c_v13_E9_subset <- loadArchRProject("~/Mef2c_ArchR_working/Mef2c_v13_E9_subset")
setwd("~/Mef2c_ArchR_working")


#CMs_IFT
markerTest_E9_CMs_IFT <- getMarkerFeatures(
  ArchRProj = proj_Mef2c_v13_E9_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "CellTypeByGenotype",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "CMs_IFT_x_KO",
  bgdGroups = "CMs_IFT_x_WT"
)

markerTest_E9_CMs_IFT_list_upinWT <- getMarkers(markerTest_E9_CMs_IFT, cutOff = "FDR <= 0.15 & Log2FC <= -0.5", returnGR = TRUE)
markerTest_E9_CMs_IFT_list_upinKO <- getMarkers(markerTest_E9_CMs_IFT, cutOff = "FDR <= 0.15 & Log2FC >= 0.5", returnGR = TRUE)

#Let's convert the upinWT list to a GRanges object, sorted by Log2FC
markerpeaksgr_E9_IFT_WT <- markerTest_E9_CMs_IFT_list_upinWT$`CMs_IFT_x_KO`
order_E9_IFT_WT <- order(markerpeaksgr_E9_IFT_WT$Log2FC, start(markerpeaksgr_E9_IFT_WT))
markerpeaks_sorted_E9_IFT_WT <- markerpeaksgr_E9_IFT_WT[order_E9_IFT_WT]

#Let's convert the upinKO list to a GRanges object, sorted by Log2FC
markerpeaksgr_E9_IFT_KO <- markerTest_E9_CMs_IFT_list_upinKO$`CMs_IFT_x_KO`
order_E9_IFT_KO <- order(markerpeaksgr_E9_IFT_KO$Log2FC, start(markerpeaksgr_E9_IFT_KO))
markerpeaks_sorted_E9_IFT_KO <- markerpeaksgr_E9_IFT_KO[order_E9_IFT_KO]

#Intersect these peaks with the full peak set to get nearestGene and peakType info
#Subset matching peaks
pgr_E9 <- getPeakSet(proj_Mef2c_v13_E9_subset)
sub_pgr_E9_IFT_upWT_idx <- findOverlaps(markerpeaks_sorted_E9_IFT_WT, pgr_E9)
sub_pgr_E9_IFT_upWT <- pgr_E9[sub_pgr_E9_IFT_upWT_idx@to]
sub_pgr_E9_IFT_upKO_idx <- findOverlaps(markerpeaks_sorted_E9_IFT_KO, pgr_E9)
sub_pgr_E9_IFT_upKO <- pgr_E9[sub_pgr_E9_IFT_upKO_idx@to]

#Create data frame for all WT DARs 
WT_DARs_E9_CMs_IFT_df <- data.frame(
  mm10PeakRange = paste0(markerpeaks_sorted_E9_IFT_WT@seqnames,':',markerpeaks_sorted_E9_IFT_WT@ranges@start,'-',markerpeaks_sorted_E9_IFT_WT@ranges@start + 500),
  peakType = sub_pgr_E9_IFT_upWT$peakType,
  peakLog2FC = markerpeaks_sorted_E9_IFT_WT$Log2FC,
  peakFDR = markerpeaks_sorted_E9_IFT_WT$FDR,
  nearestTSS =  sub_pgr_E9_IFT_upWT$nearestTSS, 
  distToTSS = sub_pgr_E9_IFT_upWT$distToTSS
)

#Create data frame for all KO DARs 
KO_DARs_E9_CMs_IFT_df <- data.frame(
  mm10PeakRange = paste0(markerpeaks_sorted_E9_IFT_KO@seqnames,':',markerpeaks_sorted_E9_IFT_KO@ranges@start,'-',markerpeaks_sorted_E9_IFT_KO@ranges@start + 500),
  peakType = sub_pgr_E9_IFT_upKO$peakType,
  peakLog2FC = markerpeaks_sorted_E9_IFT_KO$Log2FC,
  peakFDR = markerpeaks_sorted_E9_IFT_KO$FDR,
  nearestTSS =  sub_pgr_E9_IFT_upKO$nearestTSS, 
  distToTSS = sub_pgr_E9_IFT_upKO$distToTSS
)

#Combine dataframes
E9_CMs_IFT_DARs <- rbind(WT_DARs_E9_CMs_IFT_df, KO_DARs_E9_CMs_IFT_df)

# Save dataframe to a CSV file
write.csv(E9_CMs_IFT_DARs, "~/Desktop/Working_Directory/E9_CMs_IFT_DARs.csv", row.names = FALSE)


#CMs_V
markerTest_E9_CMs_V <- getMarkerFeatures(
  ArchRProj = proj_Mef2c_v13_E9_subset, 
  useMatrix = "PeakMatrix",
  groupBy = "CellTypeByGenotype",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "CMs_V_x_KO",
  bgdGroups = "CMs_V_x_WT"
)

markerTest_E9_CMs_V_list_upinWT <- getMarkers(markerTest_E9_CMs_V, cutOff = "FDR <= 0.15 & Log2FC <= -0.5", returnGR = TRUE)
markerTest_E9_CMs_V_list_upinKO <- getMarkers(markerTest_E9_CMs_V, cutOff = "FDR <= 0.15 & Log2FC >= 0.5", returnGR = TRUE)

#Let's convert the upinWT list to a GRanges object, sorted by Log2FC
markerpeaksgr_E9_V_WT <- markerTest_E9_CMs_V_list_upinWT$`CMs_V_x_KO`
order_E9_V_WT <- order(markerpeaksgr_E9_V_WT$Log2FC, start(markerpeaksgr_E9_V_WT))
markerpeaks_sorted_E9_V_WT <- markerpeaksgr_E9_V_WT[order_E9_V_WT]

#Let's convert the upinKO list to a GRanges object, sorted by Log2FC
markerpeaksgr_E9_V_KO <- markerTest_E9_CMs_V_list_upinKO$`CMs_V_x_KO`
order_E9_V_KO <- order(markerpeaksgr_E9_V_KO$Log2FC, start(markerpeaksgr_E9_V_KO))
markerpeaks_sorted_E9_V_KO <- markerpeaksgr_E9_V_KO[order_E9_V_KO]

#Intersect these peaks with the full peak set to get nearestGene and peakType info
#Subset matching peaks
pgr_E9 <- getPeakSet(proj_Mef2c_v13_E9_subset)
sub_pgr_E9_V_upWT_idx <- findOverlaps(markerpeaks_sorted_E9_V_WT, pgr_E9)
sub_pgr_E9_V_upWT <- pgr_E9[sub_pgr_E9_V_upWT_idx@to]
sub_pgr_E9_V_upKO_idx <- findOverlaps(markerpeaks_sorted_E9_V_KO, pgr_E9)
sub_pgr_E9_V_upKO <- pgr_E9[sub_pgr_E9_V_upKO_idx@to]

#Create data frame for all WT DARs 
WT_DARs_E9_CMs_V_df <- data.frame(
  mm10PeakRange = paste0(markerpeaks_sorted_E9_V_WT@seqnames,':',markerpeaks_sorted_E9_V_WT@ranges@start,'-',markerpeaks_sorted_E9_V_WT@ranges@start + 500),
  peakType = sub_pgr_E9_V_upWT$peakType,
  peakLog2FC = markerpeaks_sorted_E9_V_WT$Log2FC,
  peakFDR = markerpeaks_sorted_E9_V_WT$FDR,
  nearestTSS =  sub_pgr_E9_V_upWT$nearestTSS, 
  distToTSS = sub_pgr_E9_V_upWT$distToTSS
)

#Create data frame for all KO DARs 
KO_DARs_E9_CMs_V_df <- data.frame(
  mm10PeakRange = paste0(markerpeaks_sorted_E9_V_KO@seqnames,':',markerpeaks_sorted_E9_V_KO@ranges@start,'-',markerpeaks_sorted_E9_V_KO@ranges@start + 500),
  peakType = sub_pgr_E9_V_upKO$peakType,
  peakLog2FC = markerpeaks_sorted_E9_V_KO$Log2FC,
  peakFDR = markerpeaks_sorted_E9_V_KO$FDR,
  nearestTSS =  sub_pgr_E9_V_upKO$nearestTSS, 
  distToTSS = sub_pgr_E9_V_upKO$distToTSS
)

#Combine dataframes
E9_CMs_V_DARs <- rbind(WT_DARs_E9_CMs_V_df, KO_DARs_E9_CMs_V_df)

# Save dataframe to a CSV file
write.csv(E9_CMs_V_DARs, "~/Desktop/Working_Directory/E9_CMs_V_DARs.csv", row.names = FALSE)

