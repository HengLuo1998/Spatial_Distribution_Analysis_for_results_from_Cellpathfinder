library(tidyverse)
library(ggpubr)
library(ggplot2)
library(readxl)
library(viridis)
library(reshape)
library(dplyr)
library(nabor)
theme_set(theme_bw(base_size = 16))

cell_stats_fra = read.csv("~/Sample/SampleE.csv") # Read epithelial1 table
cell_stats_epi = read.csv("~/Samplei/SampleE.csv") # Read epithelial2 table
cell_stats_str = read.csv("~/Samplei/SampleS.csv") # Read stromal table
cell_stats_cd45 = read.csv("~/Samplei/SampleCD45.csv") # Read CD45RO+ T cell table
cell_stats_cd8 = read.csv("~/Samplei/SampleCD8.csv") # Read CD8+ T cell table
cell_stats_cd68 = read.csv("~/Samplei/SampleCD68.csv") # Read CD68+ Macrophage table

str_status = cell_stats_str$Str # Extract the stromal status information from the stromal table
cd45_status = cell_stats_cd45$CD45 # Extract the cd45 status information from the CD45RO+ T cell table
cd8_status = cell_stats_cd8$CD8 # Extract the cd8 status information from the CD8+ T cell table
cd68_status = cell_stats_cd68$Mp # Extract the macrophage status information from the CD68+ macrophage table
fra_status_mean = cell_stats_fra$X.Nucleus1..MeanIntensity.CH4 # Extract the FRA(CH4) mean intensity information from the epithelial1 table
fra_status_total = cell_stats_fra$X.Nucleus1..TotalIntensity.CH4 # Extract the FRA(CH4) total intensity information from the epithelial1 table

cell_stats_epi$Str = str_status # Add the stromal status information into the epithelial table 
cell_stats_epi$cd45 = cd45_status # Add the cd45 status information into the epithelial table 
cell_stats_epi$cd8 = cd8_status # Add the cd8 status information into the epithelial table 
cell_stats_epi$cd68 = cd68_status # Add the cd68 status information into the epithelial table 
cell_stats_epi$X.Nucleus1..TotalIntensity.CH4 = fra_status_total # Add the FRA(CH4) total intensity infiormation into the epithelial table 
cell_stats_epi$X.Nucleus1..MeanIntensity.CH4 = fra_status_mean # Add the FRA(CH4) mean intensity infiormation into the epithelial table 
cell_stats_full = cell_stats_epi

#Pre-process for T cluster
cell_stats_full$CD45_T = 0 # Non T cell
cell_stats_full$CD45_T[cell_stats_full$cd45 == 1 & cell_stats_full$cd8 == 0] = 1 # CD45RO+ T cell
cell_stats_full$CD8_T = 0 # CD8- cell
cell_stats_full$CD8_T[cell_stats_full$cd45 == 0 & cell_stats_full$cd8 == 1] = 1 # CD8+ T cell
cell_stats_full$CD8_T[cell_stats_full$cd45 == 1 & cell_stats_full$cd8 == 1] = 2 # CD45RO+ CD8+ T cell
cell_stats_full$CD68_Mp = 0 # CD68- cell
cell_stats_full$CD68_Mp[cell_stats_full$cd68 == 1 & cell_stats_full$cd45 == 0 & cell_stats_full$cd8 == 0] = 1 # CD68+ Macrophage
cell_stats_full$Imm = 0 # Non-T cell
cell_stats_full$Imm[cell_stats_full$CD45_T == 1 | cell_stats_full$CD8_T == 1 | cell_stats_full$CD8_T == 2 | cell_stats_full$CD68_Mp == 1] = 1 # Immune cell
cell_stats_full$ImmDetail = "Non Immune cell"
cell_stats_full$ImmDetail[cell_stats_full$CD45_T == 1] = "CD45RO+ T cell"
cell_stats_full$ImmDetail[cell_stats_full$CD8_T == 1] = "CD8+ T cell"
cell_stats_full$ImmDetail[cell_stats_full$CD8_T == 2] = "CD45RO+ CD8+ T cell"
cell_stats_full$ImmDetail[cell_stats_full$CD68_Mp == 1] = "CD68+ Macrophage"


#Establish column for signalling pathway and oncogenic trait level
cell_stats_full$PARP = 0
cell_stats_full$FRA1 = 0
cell_stats_full$YAP1 = 0

imm = cell_stats_full%>% # Select all the nuclei with an T state of 1 (1:Immune,0:Not T)
  filter(cell_stats_full$ImmDetail != "Non Immune cell")
T_cells = cell_stats_full%>%
  filter(cell_stats_full$ImmDetail == "CD45RO+ T cell" | cell_stats_full$ImmDetail == "CD8+ T cell" | cell_stats_full$ImmDetail == "CD45RO+ CD8+ T cell")
CD8_T_cells = cell_stats_full%>%
  filter(cell_stats_full$ImmDetail == "CD8+ T cell" | cell_stats_full$ImmDetail == "CD45RO+ CD8+ T cell")
Macrophage = cell_stats_full%>%
  filter(cell_stats_full$ImmDetail == "CD68+ Macrophage") 
epi = cell_stats_full %>% # Select all the nuclei with an epithelial state of 1 (1:Epithelial,0:Not epithelial)
  filter(cell_stats_full$Epi == 1 & cell_stats_full$Imm == 0)
str = cell_stats_full %>% # Select all the nuclei with a stromal state of 1 (1:Stromal,0:Not Stromal)
  filter(cell_stats_full$Str == 1 & cell_stats_full$Imm == 0)

#Exclude all the Str/Epi/Imm double/triple positive nuclei
epi_trim = epi %>% 
  filter(epi$Str == 0 & epi$Imm == 0) 

str_trim = str %>%
  filter(str$Epi == 0 & str$Imm == 0)


q_PARP = quantile(epi_trim$X.Nucleus1..MeanIntensity.CH3)
threshold_PARP = q_PARP[4] + 1.5*IQR(epi_trim$X.Nucleus1..MeanIntensity.CH3)
q_FRA1 = quantile(epi_trim$X.Nucleus1..MeanIntensity.CH4)
threshold_FRA1 = q_FRA1[4] + 1.5*IQR(epi_trim$X.Nucleus1..MeanIntensity.CH4)
q_YAP1 = quantile(epi_trim$X.Nucleus1..MeanIntensity.CH5)
threshold_YAP1 = q_YAP1[4] + 1.5*IQR(epi_trim$X.Nucleus1..MeanIntensity.CH5)
epi_trim$PARP[epi_trim$X.Nucleus1..MeanIntensity.CH3 >= threshold_PARP] = 1
epi_trim$FRA1[epi_trim$X.Nucleus1..MeanIntensity.CH4 >= threshold_FRA1] = 1
epi_trim$YAP1[epi_trim$X.Nucleus1..MeanIntensity.CH5 >= threshold_YAP1] = 1


################################################################################
#Updated version of script in searching epithelial cells close to the stromal-epithelial interface (d0 = 20, 30 and 40)
epi234 = epi_trim
epi234$Str_proximity20 = 0 #Total number of stromal nuclei that is within 20um to the individual epithelial nucleus
epi234$Str_proximity30 = 0 #Total number of stromal nuclei that is between 20um to 30um to the individual epithelial nucleus
epi234$Str_proximity40 = 0 #Total number of stromal nuclei that is between 30um to 40um to the individual epithelial nucleus
epi234$Imm_proximity20 = 0 #Total number of T nuclei that is within 20um to the individual epithelial nucleus
epi234$Imm_proximity30 = 0 #Total number of T nuclei that is between 20um to 30um to the individual epithelial nucleus
epi234$Imm_proximity40 = 0 #Total number of T nuclei that is between 30um to 40um to the individual epithelial nucleus
epi234$T_proximity20 = 0 
epi234$T_proximity30 = 0 
epi234$T_proximity40 = 0 
epi234$CD8_proximity20 = 0 
epi234$CD8_proximity30 = 0 
epi234$CD8_proximity40 = 0 
epi234$CD68Mp_proximity20 = 0 
epi234$CD68Mp_proximity30 = 0 
epi234$CD68Mp_proximity40 = 0 
# build kd-tree for all points
# TODO: please check, should we compare against str_trim or str? In the old code you used str... please change if I took the wrong one here
str_points <- cbind(str_trim$X.Nucleus1..X.CoordinateInWell, str_trim$X.Nucleus1..Y.CoordinateInWell)
imm_points <- cbind(imm$X.Nucleus1..X.CoordinateInWell, imm$X.Nucleus1..Y.CoordinateInWell)
T_points <- cbind(T_cells$X.Nucleus1..X.CoordinateInWell, T_cells$X.Nucleus1..Y.CoordinateInWell)
CD8_points <- cbind(CD8_T_cells$X.Nucleus1..X.CoordinateInWell, CD8_T_cells$X.Nucleus1..Y.CoordinateInWell)
CD68_points <- cbind(Macrophage$X.Nucleus1..X.CoordinateInWell, Macrophage$X.Nucleus1..Y.CoordinateInWell)

epi_points <- cbind(epi234$X.Nucleus1..X.CoordinateInWell, epi234$X.Nucleus1..Y.CoordinateInWell)
str_kdtree <- WKNND(str_points)
imm_kdtree <- WKNND(imm_points)
T_kdtree <- WKNND(T_points)
CD8_kdtree <- WKNND(CD8_points)
CD68_kdtree <- WKNND(CD68_points)

epi_query1 <- str_kdtree$query(epi_points, k=nrow(str_points), eps=0, radius=40) # ugly code: how to set k here to not query against all?
epi_query2 <- imm_kdtree$query(epi_points, k=nrow(imm_points), eps=0, radius=200) # ugly code: how to set k here to not query against all?
epi_query3 <- T_kdtree$query(epi_points, k=nrow(T_points), eps=0, radius=200)
epi_query4 <- CD8_kdtree$query(epi_points, k=nrow(CD8_points), eps=0, radius=200)
epi_query5 <- CD68_kdtree$query(epi_points, k=nrow(CD68_points), eps=0, radius=200)
# get number of cells within 20um / 30um / 40 um distance to Stromal or Immune cells
epi234$Str_proximity20 <- rowSums(epi_query1$nn.dists <= 20)
epi234$Str_proximity30 <- rowSums(epi_query1$nn.dists > 20 & epi_query1$nn.dists <= 30)
epi234$Str_proximity40 <- rowSums(epi_query1$nn.dists > 30 & epi_query1$nn.dists <= 40)

epi234$Imm_proximity20 <- rowSums(epi_query2$nn.dists <= 100)
epi234$Imm_proximity30 <- rowSums(epi_query2$nn.dists > 100 & epi_query2$nn.dists <= 150)
epi234$Imm_proximity40 <- rowSums(epi_query2$nn.dists > 150 & epi_query2$nn.dists <= 200)

epi234$T_proximity20 = rowSums(epi_query3$nn.dists <= 100) 
epi234$T_proximity30 = rowSums(epi_query3$nn.dists > 100 & epi_query3$nn.dists <= 150) 
epi234$T_proximity40 = rowSums(epi_query3$nn.dists > 150 & epi_query3$nn.dists <= 200) 

epi234$CD8_proximity20 = rowSums(epi_query4$nn.dists <= 100)
epi234$CD8_proximity30 = rowSums(epi_query4$nn.dists > 100 & epi_query4$nn.dists <= 150) 
epi234$CD8_proximity40 = rowSums(epi_query4$nn.dists > 150 & epi_query4$nn.dists <= 200) 

epi234$CD68Mp_proximity20 = rowSums(epi_query5$nn.dists <= 100) 
epi234$CD68Mp_proximity30 = rowSums(epi_query5$nn.dists > 100 & epi_query5$nn.dists <= 150) 
epi234$CD68Mp_proximity40 = rowSums(epi_query5$nn.dists > 150 & epi_query5$nn.dists <= 200) 

str_trim$Str_proximity20 = 0 #Create same column in stromal part for further row bind.
str_trim$Str_proximity30 = 0
str_trim$Str_proximity40 = 0
str_trim$Imm_proximity20 = 0 
str_trim$Imm_proximity30 = 0
str_trim$Imm_proximity40 = 0
str_trim$T_proximity20 = 0 
str_trim$T_proximity30 = 0 
str_trim$T_proximity40 = 0 
str_trim$CD8_proximity20 = 0 
str_trim$CD8_proximity30 = 0 
str_trim$CD8_proximity40 = 0 
str_trim$CD68Mp_proximity20 = 0 
str_trim$CD68Mp_proximity30 = 0 
str_trim$CD68Mp_proximity40 = 0 

imm$Str_proximity20 = 0 #Create same column in T part for further row bind.
imm$Str_proximity30 = 0
imm$Str_proximity40 = 0
imm$Imm_proximity20 = 0 
imm$Imm_proximity30 = 0
imm$Imm_proximity40 = 0
imm$T_proximity20 = 0 
imm$T_proximity30 = 0 
imm$T_proximity40 = 0 
imm$CD8_proximity20 = 0 
imm$CD8_proximity30 = 0 
imm$CD8_proximity40 = 0 
imm$CD68Mp_proximity20 = 0 
imm$CD68Mp_proximity30 = 0 
imm$CD68Mp_proximity40 = 0 

T_cells$Str_proximity20 = 0 #Create same column in T_cells part for further row bind.
T_cells$Str_proximity30 = 0
T_cells$Str_proximity40 = 0
T_cells$Imm_proximity20 = 0 
T_cells$Imm_proximity30 = 0
T_cells$Imm_proximity40 = 0
T_cells$T_proximity20 = 0 
T_cells$T_proximity30 = 0 
T_cells$T_proximity40 = 0 
T_cells$CD8_proximity20 = 0 
T_cells$CD8_proximity30 = 0 
T_cells$CD8_proximity40 = 0 
T_cells$CD68Mp_proximity20 = 0 
T_cells$CD68Mp_proximity30 = 0 
T_cells$CD68Mp_proximity40 = 0 

CD8_T_cells$Str_proximity20 = 0 #Create same column in CD8_T_cells part for further row bind.
CD8_T_cells$Str_proximity30 = 0
CD8_T_cells$Str_proximity40 = 0
CD8_T_cells$Imm_proximity20 = 0 
CD8_T_cells$Imm_proximity30 = 0
CD8_T_cells$Imm_proximity40 = 0
CD8_T_cells$T_proximity20 = 0 
CD8_T_cells$T_proximity30 = 0 
CD8_T_cells$T_proximity40 = 0 
CD8_T_cells$CD8_proximity20 = 0 
CD8_T_cells$CD8_proximity30 = 0 
CD8_T_cells$CD8_proximity40 = 0 
CD8_T_cells$CD68Mp_proximity20 = 0 
CD8_T_cells$CD68Mp_proximity30 = 0 
CD8_T_cells$CD68Mp_proximity40 = 0 

Macrophage$Str_proximity20 = 0 #Create same column in Macrophage part for further row bind.
Macrophage$Str_proximity30 = 0
Macrophage$Str_proximity40 = 0
Macrophage$Imm_proximity20 = 0 
Macrophage$Imm_proximity30 = 0
Macrophage$Imm_proximity40 = 0
Macrophage$T_proximity20 = 0 
Macrophage$T_proximity30 = 0 
Macrophage$T_proximity40 = 0 
Macrophage$CD8_proximity20 = 0 
Macrophage$CD8_proximity30 = 0 
Macrophage$CD8_proximity40 = 0 
Macrophage$CD68Mp_proximity20 = 0 
Macrophage$CD68Mp_proximity30 = 0 
Macrophage$CD68Mp_proximity40 = 0 

################################################################################
#Original table with epithelial, stromal and T cluster#####################
################################################################################
full234 = rbind(epi234,str_trim,imm) #Combine the stromal part, T part and processed epithelial part to form a full table 
full234$Str_proximity_status234 = "40 - ∞" #Initiate a new column in the full table for distance annotation
full234$Str_proximity_status234[full234$Str_proximity20 > 0] = "0 - 20"
full234$Str_proximity_status234[full234$Str_proximity20 == 0 & full234$Str_proximity30 > 0] = "20 - 30"
full234$Str_proximity_status234[full234$Str_proximity20 == 0 & full234$Str_proximity30 == 0 & full234$Str_proximity40 > 0] = "30 - 40"
full234$Str_proximity_status234[full234$Str == 1 & full234$Imm == 0] = "Stromal"
full234$Str_proximity_status234[full234$Imm == 1] = "Immune"

full234$Imm_proximity_status234 = "200 - ∞" #Initiate a new column in the full table for distance annotation
full234$Imm_proximity_status234[full234$Imm_proximity20 > 0] = "0 - 100"
full234$Imm_proximity_status234[full234$Imm_proximity20 == 0 & full234$Imm_proximity30 > 0] = "100 - 150"
full234$Imm_proximity_status234[full234$Imm_proximity20 == 0 & full234$Imm_proximity30 == 0 & full234$Imm_proximity40 > 0] = "150 - 200"
full234$Imm_proximity_status234[full234$Str == 1 & full234$Imm == 0] = "Stromal"
full234$Imm_proximity_status234[full234$Imm == 1] = "Immune"

cd8 = cell_stats_full%>%
  filter(cell_stats_full$cd8 == 1 & cell_stats_full$cd45 == 0 & cell_stats_full$cd68 == 0)
cd8_cd45 = cell_stats_full%>%
  filter(cell_stats_full$cd8 == 1 & cell_stats_full$cd45 == 1 & cell_stats_full$cd68 == 0)
cd8_cd45_cd68 = cell_stats_full%>%
  filter(cell_stats_full$cd8 == 1 & cell_stats_full$cd45 == 1 & cell_stats_full$cd68 == 1)
cd45 = cell_stats_full%>%
  filter(cell_stats_full$cd8 == 0 & cell_stats_full$cd45 == 1 & cell_stats_full$cd68 == 0)
cd45_cd68 = cell_stats_full%>%
  filter(cell_stats_full$cd8 == 0 & cell_stats_full$cd45 == 1 & cell_stats_full$cd68 == 1)
cd68 = cell_stats_full%>%
  filter(cell_stats_full$cd8 == 0 & cell_stats_full$cd45 == 0 & cell_stats_full$cd68 == 1)
cd68_cd8 = cell_stats_full%>%
  filter(cell_stats_full$cd8 == 1 & cell_stats_full$cd45 == 0 & cell_stats_full$cd68 == 1)
nSTR = nrow(str_trim)
nEPI = nrow(epi_trim)
nTOTAL = nrow(cell_stats_full)
nunkown = nTOTAL - nrow(full234)
nCD45 = nrow(cd45)
nCD8 = nrow(cd8)
nCD68 = nrow(cd68)
nCD45_CD8 = nrow(cd8_cd45)
nCD8_CD68 = nrow(cd68_cd8)
nCD45_CD68 = nrow(cd45_cd68)
nCD45_CD8_CD68 = nrow(cd8_cd45_cd68)
vfra_max = max(epi_trim$X.Nucleus1..MeanIntensity.CH4)
vfra_min = min(epi_trim$X.Nucleus1..MeanIntensity.CH4)
vyap_max = max(epi_trim$X.Nucleus1..MeanIntensity.CH5)
vyap_min = min(epi_trim$X.Nucleus1..MeanIntensity.CH5)
vparp_max = max(epi_trim$X.Nucleus1..MeanIntensity.CH3)
vparp_min = min(epi_trim$X.Nucleus1..MeanIntensity.CH3)

CELL_STATS_TABLE = data.frame(Epithelial = c(nEPI),
                              Stromal = c(nSTR),
                              CD45 = c(nCD45),
                              CD8 = c(nCD8),
                              CD68 = c(nCD68),
                              CD45_CD8 = c(nCD45_CD8),
                              CD8_CD68 = c(nCD8_CD68),
                              CD45_CD68 = c(nCD45_CD68),
                              CD45_CD8_CD68 = c(nCD45_CD8_CD68),
                              Unassigned = c(nunkown),
                              Total = c(nTOTAL),
                              fra_max = c(vfra_max),
                              fra_min = c(vfra_min),
                              yap_max = c(vyap_max),
                              yap_min = c(vyap_min),
                              parp_max = c(vparp_max),
                              parp_min = c(vparp_min))
CELL_STATS_TABLE

################################################################################
#New table with only T cells and epithelial cells###############################
################################################################################
T234 = rbind(epi234,T_cells) #Combine the stromal part, T part and processed epithelial part to form a full table 
T234$T_proximity_status234 = "200 - ∞" #Initiate a new column in the T table for distance annotation
T234$T_proximity_status234[T234$T_proximity20 > 0] = "0 - 100"
T234$T_proximity_status234[T234$T_proximity20 == 0 & T234$T_proximity30 > 0] = "100 - 150"
T234$T_proximity_status234[T234$T_proximity20 == 0 & T234$T_proximity30 == 0 & T234$T_proximity40 > 0] = "150 - 200"
T234$T_proximity_status234[T234$ImmDetail == "CD45RO+ T cell" | T234$ImmDetail == "CD8+ T cell" | T234$ImmDetail == "CD45RO+ CD8+ T cell"] = "T cells"
################################################################################
#New table with only CD8+ T cells and epithelial cells##########################
################################################################################
CD8_T234 = rbind(epi234,CD8_T_cells) #Combine the stromal part, T part and processed epithelial part to form a full table 
CD8_T234$CD8_proximity_status234 = "200 - ∞" #Initiate a new column in the CD8_T table for distance annotation
CD8_T234$CD8_proximity_status234[CD8_T234$CD8_proximity20 > 0] = "0 - 100"
CD8_T234$CD8_proximity_status234[CD8_T234$CD8_proximity20 == 0 & CD8_T234$CD8_proximity30 > 0] = "100 - 150"
CD8_T234$CD8_proximity_status234[CD8_T234$CD8_proximity20 == 0 & CD8_T234$CD8_proximity30 == 0 & CD8_T234$CD8_proximity40 > 0] = "150 - 200"
CD8_T234$CD8_proximity_status234[CD8_T234$ImmDetail == "CD8+ T cell" | CD8_T234$ImmDetail == "CD45RO+ CD8+ T cell"] = "CD8+ T cells"
################################################################################
#New table with only CD68+ Macrophage and epithelial cells######################
################################################################################
CD68Mp234 = rbind(epi234,Macrophage) #Combine the stromal part, T part and processed epithelial part to form a full table 
CD68Mp234$CD68Mp_proximity_status234 = "200 - ∞" #Initiate a new column in the CD68Mp table for distance annotation
CD68Mp234$CD68Mp_proximity_status234[CD68Mp234$CD68Mp_proximity20 > 0] = "0 - 100"
CD68Mp234$CD68Mp_proximity_status234[CD68Mp234$CD68Mp_proximity20 == 0 & CD68Mp234$CD68Mp_proximity30 > 0] = "100 - 150"
CD68Mp234$CD68Mp_proximity_status234[CD68Mp234$CD68Mp_proximity20 == 0 & CD68Mp234$CD68Mp_proximity30 == 0 & CD68Mp234$CD68Mp_proximity40 > 0] = "150 - 200"
CD68Mp234$CD68Mp_proximity_status234[CD68Mp234$ImmDetail == "CD68+ Macrophage"] = "CD68+ Macrophage"

################################################################################
###Start of distance analysis against stromal cluster###########################
################################################################################

################################################################################
#Use scatterplot to reconstruct the sptial information
#General info
Spatial_ESI = ggplot(full234, aes(x= X.Nucleus1..X.CoordinateInWell, y= X.Nucleus1..Y.CoordinateInWell, color = Str_proximity_status234)) + 
  geom_point() +
  scale_color_manual(values=c("red", "purple", "blue", "black", "green", "grey")) + 
  ggtitle("Spatial distribution of nuclei relative to stroma in Sample specimen") +
  ylab("Y - Coordinates") +
  xlab("X - Coordinates") +
  theme(legend.direction = "vertical", legend.box = "vertical")

Spatial_ESI
#PARP-high distribution
Spatial_ES_PARP = ggplot(full234, aes(x= X.Nucleus1..X.CoordinateInWell, y= X.Nucleus1..Y.CoordinateInWell, color = Str_proximity_status234)) + 
  geom_point(aes(size = factor(PARP))) +
  scale_color_manual(values=c("red", "purple", "blue", "black", "green", "grey")) + 
  theme(legend.direction = "vertical", legend.box = "vertical")

Spatial_ES_PARP
#FRA1-high distribution
Spatial_ES_FRA1 = ggplot(full234, aes(x= X.Nucleus1..X.CoordinateInWell, y= X.Nucleus1..Y.CoordinateInWell, color = Str_proximity_status234)) + 
  geom_point(aes(size = factor(FRA1))) +
  scale_color_manual(values=c("red", "purple", "blue", "black", "green", "grey")) + 
  theme(legend.direction = "vertical", legend.box = "vertical")

Spatial_ES_FRA1
#YAP1-high distribution
Spatial_ES_YAP1 = ggplot(full234, aes(x= X.Nucleus1..X.CoordinateInWell, y= X.Nucleus1..Y.CoordinateInWell, color = Str_proximity_status234)) + 
  geom_point(aes(size = factor(YAP1))) +
  scale_color_manual(values=c("red", "purple", "blue", "black", "green", "grey")) + 
  theme(legend.direction = "vertical", legend.box = "vertical")

Spatial_ES_YAP1

Vio_ES_PARP = ggplot(full234, aes(x = Str_proximity_status234, y = X.Nucleus1..MeanIntensity.CH3)) +
  geom_violin()+
  geom_boxplot(width = 0.1)+
  ggtitle("PARP Mean Intensity in epithelial cells with different distance to stroma") +
  ylab("PARP mean intensity") +
  xlab("Distance to stroma (µm)") 

Vio_ES_PARP

Vio_ES_FRA = ggplot(full234, aes(x = Str_proximity_status234, y = X.Nucleus1..MeanIntensity.CH4)) +
  geom_violin()+
  geom_boxplot(width = 0.1)+
  ggtitle("FRA1 Mean Intensity in epithelial cells with different distance to stroma") +
  ylab("FRA1 mean intensity") +
  xlab("Distance to stroma (µm)") 

Vio_ES_FRA

Vio_ES_YAP = ggplot(full234, aes(x = Str_proximity_status234, y = X.Nucleus1..MeanIntensity.CH5)) +
  geom_violin()+
  geom_boxplot(width = 0.1)+
  ggtitle("YAP1 Mean Intensity in epithelial cells with different distance to stroma") +
  ylab("YAP1 mean intensity") +
  xlab("Distance to stroma (µm)") 

Vio_ES_YAP

################################################################################
#Hypergeometric statistics for fra-high nuclei
full234_simple_parp = data.frame(PARP = full234$PARP, proximity_str = full234$Str_proximity_status234)
np20 = nrow(full234_simple_parp[full234_simple_parp$proximity_str == "0 - 20",])
np30 = nrow(full234_simple_parp[full234_simple_parp$proximity_str == "20 - 30",])
np40 = nrow(full234_simple_parp[full234_simple_parp$proximity_str == "30 - 40",])
npd = nrow(full234_simple_parp[full234_simple_parp$proximity_str == "40 - ∞",])
np20_parphigh = nrow(full234_simple_parp[full234_simple_parp$proximity_str == "0 - 20" & full234_simple_parp$PARP == 1,])
np30_parphigh = nrow(full234_simple_parp[full234_simple_parp$proximity_str == "20 - 30" & full234_simple_parp$PARP == 1,])
np40_parphigh = nrow(full234_simple_parp[full234_simple_parp$proximity_str == "30 - 40" & full234_simple_parp$PARP == 1,])
npd_parphigh = nrow(full234_simple_parp[full234_simple_parp$proximity_str == "40 - ∞" & full234_simple_parp$PARP == 1,])
nps_parphigh = nrow(full234_simple_parp[full234_simple_parp$proximity_str == "Stromal" & full234_simple_parp$PARP == 1,])
nps = nrow(full234) - np20 - np30 - np40 - npd
np_parphigh = nrow(full234_simple_parp[full234_simple_parp$PARP == 1,])
np_parplow = np20 + np30 + np40 + npd - np_parphigh

p_parphigh_20 = dhyper(np20_parphigh, np_parphigh, np_parplow, np20)
p_parphigh_30 = dhyper(np30_parphigh, np_parphigh, np_parplow, np30)
p_parphigh_40 = dhyper(np40_parphigh, np_parphigh, np_parplow, np40)
p_parphigh_distal = dhyper(npd_parphigh, np_parphigh, np_parplow, npd)

dhyper_ES_PARP = data.frame(Proximity = c("0 - 20", "20 - 30", "30 - 40", "40 - ∞"),
                            number_of_parphigh_drawn = c(np20_parphigh,np30_parphigh,np40_parphigh,npd_parphigh),
                            number_of_parphigh = c(np_parphigh,np_parphigh,np_parphigh,np_parphigh),
                            number_of_parplow = c(np_parplow,np_parplow,np_parplow,np_parplow),
                            number_of_drawn = c(np20,np30,np40,npd),
                            probability = c(p_parphigh_20, p_parphigh_30, p_parphigh_40, p_parphigh_distal))
esp20 = c(np20_parphigh,(np20 - np20_parphigh))
esp30 = c(np30_parphigh,(np30 - np30_parphigh))
esp40 = c(np40_parphigh,(np40 - np40_parphigh))
espd = c(npd_parphigh,(npd - npd_parphigh))
esp = data.frame(esp20, esp30, esp40, espd)
espc = chisq.test(esp)



#Visualize the idea that the probability of seeing parp-high epithelial nuclei increases as the distance to the stroma decreases.
table_ES_parp_distance = data.frame(PARP_low = c(np20 - np20_parphigh, np30 - np30_parphigh, np40 - np40_parphigh, npd - npd_parphigh), PARP_high = c(np20_parphigh, np30_parphigh, np40_parphigh, npd_parphigh), row.names = c("0 - 20","20 - 30","30 - 40","40 - ∞"))

table_ES_parp_distance$Distance_to_stroma_µm = row.names(table_ES_parp_distance)
mtable_ES_parp_distance = melt(table_ES_parp_distance, id.vars = "Distance_to_stroma_µm")
bar_ES_parp_distance = ggplot(mtable_ES_parp_distance, aes(Distance_to_stroma_µm, value, fill = variable)) +
  geom_bar(position = "fill", stat = "identity")+
  ggtitle("Relationship between the probability of seeing a PARP-1 high epithelial nuclei and the distance of the nuclei to the stroma")
bar_ES_parp_distance

per_ES_parphigh_distance = data.frame (PARP_high_fraction = c(np20_parphigh / np20, np30_parphigh / np30, np40_parphigh / np40, npd_parphigh / npd),Distance_to_stroma_µm = c("0 - 20","20 - 30","30 - 40","40 - ∞") )
bar_per_ES_parphigh_distance = ggplot(per_ES_parphigh_distance, aes(x = Distance_to_stroma_µm, y = PARP_high_fraction)) +
  geom_bar(aes(fill = Distance_to_stroma_µm), position = "dodge", stat = "identity")+
  scale_fill_manual(values=c("red", "purple", "blue", "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
bar_per_ES_parphigh_distance

ES_distance = data.frame (Number_of_nuclei = c(np20, np30, np40, npd),Distance_to_stroma_µm = c("0 - 20","20 - 30","30 - 40","40 - ∞") )
bar_ES_distance = ggplot(ES_distance, aes(x = Distance_to_stroma_µm, y = Number_of_nuclei)) +
  geom_bar(aes(fill = Distance_to_stroma_µm), position = "dodge", stat = "identity")+
  scale_fill_manual(values=c("red", "purple", "blue", "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle("Number of epithelial nuclei with different distances to the stroma")
bar_ES_distance

ES_parphigh_distance = data.frame (Number_of_nuclei = c(np20_parphigh, np30_parphigh, np40_parphigh, npd_parphigh),Distance_to_stroma_µm = c("0 - 20","20 - 30","30 - 40","40 - ∞") )
bar_ES_parphigh_distance = ggplot(ES_parphigh_distance, aes(x = Distance_to_stroma_µm, y = Number_of_nuclei)) +
  geom_bar(aes(fill = Distance_to_stroma_µm), position = "dodge", stat = "identity")+
  scale_fill_manual(values=c("red", "purple", "blue", "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle("Number of parp-high epithelial nuclei with different distances to the stroma")
bar_ES_parphigh_distance


################################################################################
#Hypergeometric statistics for fra-high nuclei
full234_simple_fra = data.frame(FRA1 = full234$FRA1, proximity_str = full234$Str_proximity_status234)
nf20 = nrow(full234_simple_fra[full234_simple_fra$proximity_str == "0 - 20",])
nf30 = nrow(full234_simple_fra[full234_simple_fra$proximity_str == "20 - 30",])
nf40 = nrow(full234_simple_fra[full234_simple_fra$proximity_str == "30 - 40",])
nfd = nrow(full234_simple_fra[full234_simple_fra$proximity_str == "40 - ∞",])
nf20_frahigh = nrow(full234_simple_fra[full234_simple_fra$proximity_str == "0 - 20" & full234_simple_fra$FRA1 == 1,])
nf30_frahigh = nrow(full234_simple_fra[full234_simple_fra$proximity_str == "20 - 30" & full234_simple_fra$FRA1 == 1,])
nf40_frahigh = nrow(full234_simple_fra[full234_simple_fra$proximity_str == "30 - 40" & full234_simple_fra$FRA1 == 1,])
nfd_frahigh = nrow(full234_simple_fra[full234_simple_fra$proximity_str == "40 - ∞" & full234_simple_fra$FRA1 == 1,])
nfs_frahigh = nrow(full234_simple_fra[full234_simple_fra$proximity_str == "Stromal" & full234_simple_fra$FRA1 == 1,])
nfs = nrow(full234) - nf20 - nf30 - nf40 - nfd
nf_frahigh = nf20_frahigh + nf30_frahigh + nf40_frahigh + nfd_frahigh
nf_fralow = nf20 + nf30 + nf40 + nfd - nf_frahigh  

p_frahigh_20 = dhyper(nf20_frahigh, nf_frahigh, nf_fralow, nf20)
p_frahigh_30 = dhyper(nf30_frahigh, nf_frahigh, nf_fralow, nf30)
p_frahigh_40 = dhyper(nf40_frahigh, nf_frahigh, nf_fralow, nf40)
p_frahigh_distal = dhyper(nfd_frahigh, nf_frahigh, nf_fralow, nfd)

dhyper_ES_FRA = data.frame(Proximity = c("0 - 20", "20 - 30", "30 - 40", "40 - ∞"),
                           number_of_frahigh_drawn = c(nf20_frahigh,nf30_frahigh,nf40_frahigh,nfd_frahigh),
                           number_of_frahigh = c(nf_frahigh,nf_frahigh,nf_frahigh,nf_frahigh),
                           number_of_fralow = c(nf_fralow,nf_fralow,nf_fralow,nf_fralow),
                           number_of_drawn = c(nf20,nf30,nf40,nfd),
                           probability = c(p_frahigh_20, p_frahigh_30, p_frahigh_40, p_frahigh_distal))
dhyper_ES_FRA

esf20 = c(nf20_frahigh,(nf20 - nf20_frahigh))
esf30 = c(nf30_frahigh,(nf30 - nf30_frahigh))
esf40 = c(nf40_frahigh,(nf40 - nf40_frahigh))
esfd = c(nfd_frahigh,(nfd - nfd_frahigh))
esf = data.frame(esf20, esf30, esf40, esfd)
esfc = chisq.test(esf)

#Visualize the idea that the probability of seeing fra-high epithelial nuclei increases as the distance to the stroma decreases.
table_ES_fra1_distance = data.frame(Fra_low = c(nf20 - nf20_frahigh, nf30 - nf30_frahigh, nf40 - nf40_frahigh, nfd - nfd_frahigh), Fra_high = c(nf20_frahigh, nf30_frahigh, nf40_frahigh, nfd_frahigh), row.names = c("0 - 20","20 - 30","30 - 40","40 - ∞"))

table_ES_fra1_distance$Distance_to_stroma_µm = row.names(table_ES_fra1_distance)
mtable_ES_fra1_distance = melt(table_ES_fra1_distance, id.vars = "Distance_to_stroma_µm")
bar_ES_fra1_distance = ggplot(mtable_ES_fra1_distance, aes(Distance_to_stroma_µm, value, fill = variable)) +
  geom_bar(position = "fill", stat = "identity")+
  ggtitle("Relationship between the probability of seeing a Fra-1 high epithelial nuclei and the distance of the nuclei to the stroma")
bar_ES_fra1_distance

per_ES_frahigh_distance = data.frame (Fra_high_fraction = c(nf20_frahigh / nf20, nf30_frahigh / nf30, nf40_frahigh / nf40, nfd_frahigh / nfd),Distance_to_stroma_µm = c("0 - 20","20 - 30","30 - 40","40 - ∞") )
bar_per_ES_frahigh_distance = ggplot(per_ES_frahigh_distance, aes(x = Distance_to_stroma_µm, y = Fra_high_fraction)) +
  geom_bar(aes(fill = Distance_to_stroma_µm), position = "dodge", stat = "identity")+
  scale_fill_manual(values=c("red", "purple", "blue", "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
bar_per_ES_frahigh_distance

ES_distance = data.frame (Number_of_nuclei = c(nf20, nf30, nf40, nfd),Distance_to_stroma_µm = c("0 - 20","20 - 30","30 - 40","40 - ∞") )
bar_ES_distance = ggplot(ES_distance, aes(x = Distance_to_stroma_µm, y = Number_of_nuclei)) +
  geom_bar(aes(fill = Distance_to_stroma_µm), position = "dodge", stat = "identity")+
  scale_fill_manual(values=c("red", "purple", "blue", "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle("Number of epithelial nuclei with different distances to the stroma")
bar_ES_distance

ES_frahigh_distance = data.frame (Number_of_nuclei = c(nf20_frahigh, nf30_frahigh, nf40_frahigh, nfd_frahigh),Distance_to_stroma_µm = c("0 - 20","20 - 30","30 - 40","40 - ∞") )
bar_ES_frahigh_distance = ggplot(ES_frahigh_distance, aes(x = Distance_to_stroma_µm, y = Number_of_nuclei)) +
  geom_bar(aes(fill = Distance_to_stroma_µm), position = "dodge", stat = "identity")+
  scale_fill_manual(values=c("red", "purple", "blue", "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle("Number of fra-high epithelial nuclei with different distances to the stroma")
bar_ES_frahigh_distance

################################################################################
#Hypergeometric statistics for yap-high nuclei
full234_simple_yap = data.frame(YAP1 = full234$YAP1, proximity_str = full234$Str_proximity_status234)
ny20 = nrow(full234_simple_yap[full234_simple_yap$proximity_str == "0 - 20",])
ny30 = nrow(full234_simple_yap[full234_simple_yap$proximity_str == "20 - 30",])
ny40 = nrow(full234_simple_yap[full234_simple_yap$proximity_str == "30 - 40",])
nyd = nrow(full234_simple_yap[full234_simple_yap$proximity_str == "40 - ∞",])
ny20_yaphigh = nrow(full234_simple_yap[full234_simple_yap$proximity_str == "0 - 20" & full234_simple_yap$YAP1 == 1,])
ny30_yaphigh = nrow(full234_simple_yap[full234_simple_yap$proximity_str == "20 - 30" & full234_simple_yap$YAP1 == 1,])
ny40_yaphigh = nrow(full234_simple_yap[full234_simple_yap$proximity_str == "30 - 40" & full234_simple_yap$YAP1 == 1,])
nyd_yaphigh = nrow(full234_simple_yap[full234_simple_yap$proximity_str == "40 - ∞" & full234_simple_yap$YAP1 == 1,])
nys_yaphigh = nrow(full234_simple_yap[full234_simple_yap$proximity_str == "Stromal" & full234_simple_yap$YAP1 == 1,])
nys = nrow(full234) - ny20 - ny30 - ny40 - nyd
ny_yaphigh = ny20_yaphigh + ny30_yaphigh + ny40_yaphigh + nyd_yaphigh
ny_yaplow = ny20 + ny30 + ny40 + nyd - ny_yaphigh  

p_yaphigh_20 = dhyper(ny20_yaphigh, ny_yaphigh, ny_yaplow, ny20)
p_yaphigh_30 = dhyper(ny30_yaphigh, ny_yaphigh, ny_yaplow, ny30)
p_yaphigh_40 = dhyper(ny40_yaphigh, ny_yaphigh, ny_yaplow, ny40)
p_yaphigh_distal = dhyper(nyd_yaphigh, ny_yaphigh, ny_yaplow, nyd)

dhyper_ES_YAP1 = data.frame(Proximity = c("0 - 20", "20 - 30", "30 - 40", "40 - ∞"),
                            number_of_yaphigh_drawn = c(ny20_yaphigh,ny30_yaphigh,ny40_yaphigh,nyd_yaphigh),
                            number_of_yaphigh = c(ny_yaphigh,ny_yaphigh,ny_yaphigh,ny_yaphigh),
                            number_of_yaplow = c(ny_yaplow,ny_yaplow,ny_yaplow,ny_yaplow),
                            number_of_drawn = c(ny20,ny30,ny40,nyd),
                            probability = c(p_yaphigh_20, p_yaphigh_30, p_yaphigh_40, p_yaphigh_distal))

esy20 = c(ny20_yaphigh,(ny20 - ny20_yaphigh))
esy30 = c(ny30_yaphigh,(ny30 - ny30_yaphigh))
esy40 = c(ny40_yaphigh,(ny40 - ny40_yaphigh))
esyd = c(nyd_yaphigh,(nyd - nyd_yaphigh))
esy = data.frame(esy20, esy30, esy40, esyd)
esyc = chisq.test(esy)


#Visualize the idea that the probability of seeing yap-high epithelial nuclei increases as the distance to the stroma decreases.
table_ES_yap_distance = data.frame(YAP1_low = c(ny20 - ny20_yaphigh, ny30 - ny30_yaphigh, ny40 - ny40_yaphigh, nyd - nyd_yaphigh), YAP1_high = c(ny20_yaphigh, ny30_yaphigh, ny40_yaphigh, nyd_yaphigh), row.names = c("0 - 20","20 - 30","30 - 40","40 - ∞"))

table_ES_yap_distance$Distance_to_stroma_µm = row.names(table_ES_yap_distance)
mtable_ES_yap_distance = melt(table_ES_yap_distance, id.vars = "Distance_to_stroma_µm")
bar_ES_yap_distance = ggplot(mtable_ES_yap_distance, aes(Distance_to_stroma_µm, value, fill = variable)) +
  geom_bar(position = "fill", stat = "identity")+
  ggtitle("Relationship between the probability of seeing a YAP1-1 high epithelial nuclei and the distance of the nuclei to the stroma")
bar_ES_yap_distance

per_ES_yaphigh_distance = data.frame (YAP1_high_fraction = c(ny20_yaphigh / ny20, ny30_yaphigh / ny30, ny40_yaphigh / ny40, nyd_yaphigh / nyd),Distance_to_stroma_µm = c("0 - 20","20 - 30","30 - 40","40 - ∞") )
bar_per_ES_yaphigh_distance = ggplot(per_ES_yaphigh_distance, aes(x = Distance_to_stroma_µm, y = YAP1_high_fraction)) +
  geom_bar(aes(fill = Distance_to_stroma_µm), position = "dodge", stat = "identity")+
  scale_fill_manual(values=c("red", "purple", "blue", "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
bar_per_ES_yaphigh_distance

ES_distance = data.frame (Number_of_nuclei = c(ny20, ny30, ny40, nyd),Distance_to_stroma_µm = c("0 - 20","20 - 30","30 - 40","40 - ∞") )
bar_ES_distance = ggplot(ES_distance, aes(x = Distance_to_stroma_µm, y = Number_of_nuclei)) +
  geom_bar(aes(fill = Distance_to_stroma_µm), position = "dodge", stat = "identity")+
  scale_fill_manual(values=c("red", "purple", "blue", "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
bar_ES_distance

ES_yaphigh_distance = data.frame (Number_of_nuclei = c(ny20_yaphigh, ny30_yaphigh, ny40_yaphigh, nyd_yaphigh),Distance_to_stroma_µm = c("0 - 20","20 - 30","30 - 40","40 - ∞") )
bar_ES_yaphigh_distance = ggplot(ES_yaphigh_distance, aes(x = Distance_to_stroma_µm, y = Number_of_nuclei)) +
  geom_bar(aes(fill = Distance_to_stroma_µm), position = "dodge", stat = "identity")+
  scale_fill_manual(values=c("red", "purple", "blue", "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle("Number of yap-high epithelial nuclei with different distances to the stroma")
bar_ES_yaphigh_distance

################################################################################
###End of distance analysis against stromal cluster#############################
################################################################################


################################################################################
###Start of distance analysis against immune cluster############################
################################################################################

################################################################################
#Use scatterplot to reconstruct the sptial information
#General info
Spatial_EIS = ggplot(full234, aes(x= X.Nucleus1..X.CoordinateInWell, y= X.Nucleus1..Y.CoordinateInWell, color = Imm_proximity_status234)) + 
  geom_point() +
  scale_color_manual(values=c("red", "purple", "blue", "black", "green", "grey")) + 
  ggtitle("Spatial distribution of nuclei relative to immune cells in Sample specimen") +
  ylab("Y - Coordinates") +
  xlab("X - Coordinates") +
  theme(legend.direction = "vertical", legend.box = "vertical")

Spatial_EIS

Spatial_ET = ggplot(T234, aes(x= X.Nucleus1..X.CoordinateInWell, y= X.Nucleus1..Y.CoordinateInWell, color = T_proximity_status234)) + 
  geom_point() +
  scale_color_manual(values=c("red", "purple", "blue", "black", "green", "grey")) + 
  ggtitle("Spatial distribution of nuclei relative to immune cells in Sample specimen") +
  ylab("Y - Coordinates") +
  xlab("X - Coordinates") +
  theme(legend.direction = "vertical", legend.box = "vertical")

Spatial_ET

Spatial_E8 = ggplot(CD8_T234, aes(x= X.Nucleus1..X.CoordinateInWell, y= X.Nucleus1..Y.CoordinateInWell, color = CD8_proximity_status234)) + 
  geom_point() +
  scale_color_manual(values=c("red", "purple", "blue", "black", "green", "grey")) + 
  theme(legend.direction = "vertical", legend.box = "vertical")

Spatial_E8

Spatial_EM = ggplot(CD68Mp234, aes(x= X.Nucleus1..X.CoordinateInWell, y= X.Nucleus1..Y.CoordinateInWell, color = CD68Mp_proximity_status234)) + 
  geom_point() +
  scale_color_manual(values=c("red", "purple", "blue", "black", "green", "grey")) + 
  theme(legend.direction = "vertical", legend.box = "vertical")

Spatial_EM
#PARP-high distribution
Spatial_EI_PARP = ggplot(full234, aes(x= X.Nucleus1..X.CoordinateInWell, y= X.Nucleus1..Y.CoordinateInWell, color = Imm_proximity_status234)) + 
  geom_point(aes(size = factor(PARP))) +
  scale_color_manual(values=c("red", "purple", "blue", "black", "green", "grey")) + 
  ggtitle("Spatial distribution of PARP-high nuclei relative to immune cells in Sample specimen") +
  ylab("Y - Coordinates") +
  xlab("X - Coordinates") +
  theme(legend.direction = "vertical", legend.box = "vertical")

Spatial_EI_PARP
#FRA1-high distribution
Spatial_EI_FRA1 = ggplot(full234, aes(x= X.Nucleus1..X.CoordinateInWell, y= X.Nucleus1..Y.CoordinateInWell, color = Imm_proximity_status234)) + 
  geom_point(aes(size = factor(FRA1))) +
  scale_color_manual(values=c("red", "purple", "blue", "black", "green", "grey")) + 
  ggtitle("Spatial distribution of FRA1-high nuclei relative to immune cells in Sample specimen") +
  ylab("Y - Coordinates") +
  xlab("X - Coordinates") +
  theme(legend.direction = "vertical", legend.box = "vertical")

Spatial_EI_FRA1
#YAP1-high distribution
Spatial_EI_YAP1 = ggplot(full234, aes(x= X.Nucleus1..X.CoordinateInWell, y= X.Nucleus1..Y.CoordinateInWell, color = Imm_proximity_status234)) + 
  geom_point(aes(size = factor(YAP1))) +
  scale_color_manual(values=c("red", "purple", "blue", "black", "green", "grey")) + 
  ggtitle("Spatial distribution of YAP1-high nuclei relative to immune cells in Sample specimen") +
  ylab("Y - Coordinates") +
  xlab("X - Coordinates") +
  theme(legend.direction = "vertical", legend.box = "vertical")

Spatial_EI_YAP1
#PARP-high distribution
Spatial_ET_PARP = ggplot(T234, aes(x= X.Nucleus1..X.CoordinateInWell, y= X.Nucleus1..Y.CoordinateInWell, color = T_proximity_status234)) + 
  geom_point(aes(size = factor(PARP))) +
  scale_color_manual(values=c("red", "purple", "blue", "black", "green", "grey")) + 
  ggtitle("Spatial distribution of PARP-high nuclei relative to T cells in Sample specimen") +
  ylab("Y - Coordinates") +
  xlab("X - Coordinates") +
  theme(legend.direction = "vertical", legend.box = "vertical")

Spatial_ET_PARP
#FRA1-high distribution
Spatial_ET_FRA1 = ggplot(T234, aes(x= X.Nucleus1..X.CoordinateInWell, y= X.Nucleus1..Y.CoordinateInWell, color = T_proximity_status234)) + 
  geom_point(aes(size = factor(FRA1))) +
  scale_color_manual(values=c("red", "purple", "blue", "black", "green", "grey")) + 
  ggtitle("Spatial distribution of FRA1-high nuclei relative to T cells in Sample specimen") +
  ylab("Y - Coordinates") +
  xlab("X - Coordinates") +
  theme(legend.direction = "vertical", legend.box = "vertical")

Spatial_ET_FRA1
#YAP1-high distribution
Spatial_ET_YAP1 = ggplot(T234, aes(x= X.Nucleus1..X.CoordinateInWell, y= X.Nucleus1..Y.CoordinateInWell, color = T_proximity_status234)) + 
  geom_point(aes(size = factor(YAP1))) +
  scale_color_manual(values=c("red", "purple", "blue", "black", "green", "grey")) + 
  ggtitle("Spatial distribution of YAP1-high nuclei relative to T cells in Sample specimen") +
  ylab("Y - Coordinates") +
  xlab("X - Coordinates") +
  theme(legend.direction = "vertical", legend.box = "vertical")

Spatial_ET_YAP1
#PARP-high distribution
Spatial_E8_PARP = ggplot(CD8_T234, aes(x= X.Nucleus1..X.CoordinateInWell, y= X.Nucleus1..Y.CoordinateInWell, color = CD8_proximity_status234)) + 
  geom_point(aes(size = factor(PARP))) +
  scale_color_manual(values=c("red", "purple", "blue", "black", "green", "grey")) + 
  theme(legend.direction = "vertical", legend.box = "vertical")

Spatial_E8_PARP
#FRA1-high distribution
Spatial_E8_FRA1 = ggplot(CD8_T234, aes(x= X.Nucleus1..X.CoordinateInWell, y= X.Nucleus1..Y.CoordinateInWell, color = CD8_proximity_status234)) + 
  geom_point(aes(size = factor(FRA1))) +
  scale_color_manual(values=c("red", "purple", "blue", "black", "green", "grey")) + 
  theme(legend.direction = "vertical", legend.box = "vertical")

Spatial_E8_FRA1
#YAP1-high distribution
Spatial_E8_YAP1 = ggplot(CD8_T234, aes(x= X.Nucleus1..X.CoordinateInWell, y= X.Nucleus1..Y.CoordinateInWell, color = CD8_proximity_status234)) + 
  geom_point(aes(size = factor(YAP1))) +
  scale_color_manual(values=c("red", "purple", "blue", "black", "green", "grey")) + 
  theme(legend.direction = "vertical", legend.box = "vertical")

Spatial_E8_YAP1
#PARP-high distribution
Spatial_EM_PARP = ggplot(CD68Mp234, aes(x= X.Nucleus1..X.CoordinateInWell, y= X.Nucleus1..Y.CoordinateInWell, color = CD68Mp_proximity_status234)) + 
  geom_point(aes(size = factor(PARP))) +
  scale_color_manual(values=c("red", "purple", "blue", "black", "green", "grey")) + 
  theme(legend.direction = "vertical", legend.box = "vertical")

Spatial_EM_PARP
#FRA1-high distribution
Spatial_EM_FRA1 = ggplot(CD68Mp234, aes(x= X.Nucleus1..X.CoordinateInWell, y= X.Nucleus1..Y.CoordinateInWell, color = CD68Mp_proximity_status234)) + 
  geom_point(aes(size = factor(FRA1))) +
  scale_color_manual(values=c("red", "purple", "blue", "black", "green", "grey")) + 
  theme(legend.direction = "vertical", legend.box = "vertical")

Spatial_EM_FRA1
#YAP1-high distribution
Spatial_EM_YAP1 = ggplot(CD68Mp234, aes(x= X.Nucleus1..X.CoordinateInWell, y= X.Nucleus1..Y.CoordinateInWell, color = CD68Mp_proximity_status234)) + 
  geom_point(aes(size = factor(YAP1))) +
  scale_color_manual(values=c("red", "purple", "blue", "black", "green", "grey")) + 
  theme(legend.direction = "vertical", legend.box = "vertical")

Spatial_EM_YAP1
################################################################################
Vio_EI_PARP = ggplot(full234, aes(x = Imm_proximity_status234, y = X.Nucleus1..MeanIntensity.CH3)) +
  geom_violin()+
  geom_boxplot(width = 0.1)+
  ggtitle("PARP Mean Intensity in epithelial cells with different distance to immune cells") +
  ylab("PARP mean intensity") +
  xlab("Distance to stroma (µm)") 

Vio_EI_PARP

Vio_EI_FRA = ggplot(full234, aes(x = Imm_proximity_status234, y = X.Nucleus1..MeanIntensity.CH4)) +
  geom_violin()+
  geom_boxplot(width = 0.1)+
  ggtitle("FRA1 Mean Intensity in epithelial cells with different distance to immune cells") +
  ylab("FRA1 mean intensity") +
  xlab("Distance to stroma (µm)") 

Vio_EI_FRA

Vio_EI_YAP = ggplot(full234, aes(x = Imm_proximity_status234, y = X.Nucleus1..MeanIntensity.CH5)) +
  geom_violin()+
  geom_boxplot(width = 0.1)+
  ggtitle("YAP1 Mean Intensity in epithelial cells with different distance to immune cells") +
  ylab("YAP1 mean intensity") +
  xlab("Distance to stroma (µm)") 

Vio_EI_YAP

Vio_ET_PARP = ggplot(T234, aes(x = T_proximity_status234, y = X.Nucleus1..MeanIntensity.CH3)) +
  geom_violin()+
  geom_boxplot(width = 0.1)+
  ggtitle("PARP Mean Intensity in epithelial cells with different distance to T cells") +
  ylab("PARP mean intensity") +
  xlab("Distance to stroma (µm)") 

Vio_ET_PARP

Vio_ET_FRA = ggplot(T234, aes(x = T_proximity_status234, y = X.Nucleus1..MeanIntensity.CH4)) +
  geom_violin()+
  geom_boxplot(width = 0.1)+
  ggtitle("FRA1 Mean Intensity in epithelial cells with different distance to T cells") +
  ylab("FRA1 mean intensity") +
  xlab("Distance to stroma (µm)") 

Vio_ET_FRA

Vio_ET_YAP = ggplot(T234, aes(x = T_proximity_status234, y = X.Nucleus1..MeanIntensity.CH5)) +
  geom_violin()+
  geom_boxplot(width = 0.1)+
  ggtitle("YAP1 Mean Intensity in epithelial cells with different distance to T cells") +
  ylab("YAP1 mean intensity") +
  xlab("Distance to stroma (µm)") 

Vio_ET_YAP

Vio_E8_PARP = ggplot(CD8_T234, aes(x = CD8_proximity_status234, y = X.Nucleus1..MeanIntensity.CH3)) +
  geom_violin()+
  geom_boxplot(width = 0.1)+
  ggtitle("PARP Mean Intensity in epithelial cells with different distance to CD8+ T cells") +
  ylab("PARP mean intensity") +
  xlab("Distance to stroma (µm)") 

Vio_E8_PARP

Vio_E8_FRA = ggplot(CD8_T234, aes(x = CD8_proximity_status234, y = X.Nucleus1..MeanIntensity.CH4)) +
  geom_violin()+
  geom_boxplot(width = 0.1)+
  ggtitle("FRA1 Mean Intensity in epithelial cells with different distance to CD8+ T cells") +
  ylab("FRA1 mean intensity") +
  xlab("Distance to stroma (µm)") 

Vio_E8_FRA

Vio_E8_YAP = ggplot(CD8_T234, aes(x = CD8_proximity_status234, y = X.Nucleus1..MeanIntensity.CH5)) +
  geom_violin()+
  geom_boxplot(width = 0.1)+
  ggtitle("YAP1 Mean Intensity in epithelial cells with different distance to CD8+ T cells") +
  ylab("YAP1 mean intensity") +
  xlab("Distance to stroma (µm)") 

Vio_E8_YAP

Vio_EM_PARP = ggplot(CD68Mp234, aes(x = CD68Mp_proximity_status234, y = X.Nucleus1..MeanIntensity.CH3)) +
  geom_violin()+
  geom_boxplot(width = 0.1)+
  ggtitle("PARP Mean Intensity in epithelial cells with different distance to CD68+ macrophage") +
  ylab("PARP mean intensity") +
  xlab("Distance to stroma (µm)") 

Vio_EM_PARP

Vio_EM_FRA = ggplot(CD68Mp234, aes(x = CD68Mp_proximity_status234, y = X.Nucleus1..MeanIntensity.CH4)) +
  geom_violin()+
  geom_boxplot(width = 0.1)+
  ggtitle("FRA1 Mean Intensity in epithelial cells with different distance to CD68+ macrophage") +
  ylab("FRA1 mean intensity") +
  xlab("Distance to stroma (µm)") 

Vio_EM_FRA

Vio_EM_YAP = ggplot(CD68Mp234, aes(x = CD68Mp_proximity_status234, y = X.Nucleus1..MeanIntensity.CH5)) +
  geom_violin()+
  geom_boxplot(width = 0.1)+
  ggtitle("YAP1 Mean Intensity in epithelial cells with different distance to CD68+ macrophage") +
  ylab("YAP1 mean intensity") +
  xlab("Distance to stroma (µm)") 

Vio_EM_YAP
################################################################################
################################################################################
#Hypergeometric statistics for parp-high nuclei imm
full234_simple_parpi = data.frame(PARP = full234$PARP, proximity_imm = full234$Imm_proximity_status234)
npi20 = nrow(full234_simple_parpi[full234_simple_parpi$proximity_imm == "0 - 100",])
npi30 = nrow(full234_simple_parpi[full234_simple_parpi$proximity_imm == "100 - 150",])
npi40 = nrow(full234_simple_parpi[full234_simple_parpi$proximity_imm == "150 - 200",])
npid = nrow(full234_simple_parpi[full234_simple_parpi$proximity_imm == "200 - ∞",])
npi20_parphigh = nrow(full234_simple_parpi[full234_simple_parpi$proximity_imm == "0 - 100" & full234_simple_parpi$PARP == 1,])
npi30_parphigh = nrow(full234_simple_parpi[full234_simple_parpi$proximity_imm == "100 - 150" & full234_simple_parpi$PARP == 1,])
npi40_parphigh = nrow(full234_simple_parpi[full234_simple_parpi$proximity_imm == "150 - 200" & full234_simple_parpi$PARP == 1,])
npid_parphigh = nrow(full234_simple_parpi[full234_simple_parpi$proximity_imm == "200 - ∞" & full234_simple_parpi$PARP == 1,])
npis_parphigh = nrow(full234_simple_parpi[full234_simple_parpi$proximity_imm == "Stromal" & full234_simple_parpi$PARP == 1,])
npis = nrow(full234) - npi20 - npi30 - npi40 - npid
npi_parphigh = npi20_parphigh + npi30_parphigh + npi40_parphigh + npid_parphigh
npi_parplow = npi20 + npi30 + npi40 + npid - npi_parphigh  

pi_parphigh_20 = dhyper(npi20_parphigh, npi_parphigh, npi_parplow, npi20)
pi_parphigh_30 = dhyper(npi30_parphigh, npi_parphigh, npi_parplow, npi30)
pi_parphigh_40 = dhyper(npi40_parphigh, npi_parphigh, npi_parplow, npi40)
pi_parphigh_distal = dhyper(npid_parphigh, npi_parphigh, npi_parplow, npid)

dhyper_EI_PARP = data.frame(Proximity = c("0 - 100", "100 - 150", "150 - 200", "200 - ∞"),
                            number_of_parphigh_drawn = c(npi20_parphigh,npi30_parphigh,npi40_parphigh,npid_parphigh),
                            number_of_parphigh = c(npi_parphigh,npi_parphigh,npi_parphigh,npi_parphigh),
                            number_of_parplow = c(npi_parplow,npi_parplow,npi_parplow,npi_parplow),
                            number_of_drawn = c(npi20,npi30,npi40,npid),
                            probability = c(pi_parphigh_20, pi_parphigh_30, pi_parphigh_40, pi_parphigh_distal))

#Visualize the idea that the probability of seeing parp-high epithelial nuclei increases as the distance to the immune cells decreases.
table_EI_parp_distance = data.frame(PARP_low = c(npi20 - npi20_parphigh, npi30 - npi30_parphigh, npi40 - npi40_parphigh, npid - npid_parphigh), PARP_high = c(npi20_parphigh, npi30_parphigh, npi40_parphigh, npid_parphigh), row.names = c("0 - 100","100 - 150","150 - 200","200 - ∞"))

table_EI_parp_distance$Distance_to_immune_cells_µm = row.names(table_EI_parp_distance)
mtable_EI_parp_distance = melt(table_EI_parp_distance, id.vars = "Distance_to_immune_cells_µm")
bar_EI_parp_distance = ggplot(mtable_EI_parp_distance, aes(Distance_to_immune_cells_µm, value, fill = variable)) +
  geom_bar(position = "fill", stat = "identity")+
  ggtitle("Relationship between the probability of seeing a PARP-1 high epithelial nuclei and the distance of the nuclei to the immune cells")
bar_EI_parp_distance

per_EI_parphigh_distance = data.frame (PARP_high_fraction = c(npi20_parphigh / npi20, npi30_parphigh / npi30, npi40_parphigh / npi40, npid_parphigh / npid),Distance_to_immune_cells_µm = c("0 - 100","100 - 150","150 - 200","200 - ∞"))
bar_per_EI_parphigh_distance = ggplot(per_EI_parphigh_distance, aes(x = Distance_to_immune_cells_µm, y = PARP_high_fraction)) +
  geom_bar(aes(fill = Distance_to_immune_cells_µm), position = "dodge", stat = "identity")+
  scale_fill_manual(values=c("red", "purple", "blue", "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle("Percentage of PARP-high epithelial nuclei in the epithelial nuclei with different distances to the immune cells")
bar_per_EI_parphigh_distance

EI_distance = data.frame (Number_of_nuclei = c(npi20, npi30, npi40, npid),Distance_to_immune_cells_µm = c("0 - 100","100 - 150","150 - 200","200 - ∞"))
bar_EI_distance = ggplot(EI_distance, aes(x = Distance_to_immune_cells_µm, y = Number_of_nuclei)) +
  geom_bar(aes(fill = Distance_to_immune_cells_µm), position = "dodge", stat = "identity")+
  scale_fill_manual(values=c("red", "purple", "blue", "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle("Number of epithelial nuclei with different distances to the immune cells")
bar_EI_distance

EI_parphigh_distance = data.frame (Number_of_nuclei = c(npi20_parphigh, npi30_parphigh, npi40_parphigh, npid_parphigh),Distance_to_immune_cells_µm = c("0 - 100","100 - 150","150 - 200","200 - ∞"))
bar_EI_parphigh_distance = ggplot(EI_parphigh_distance, aes(x = Distance_to_immune_cells_µm, y = Number_of_nuclei)) +
  geom_bar(aes(fill = Distance_to_immune_cells_µm), position = "dodge", stat = "identity")+
  scale_fill_manual(values=c("red", "purple", "blue", "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle("Number of parp-high epithelial nuclei with different distances to the immune cells")
bar_EI_parphigh_distance


################################################################################
#Hypergeometric statistics for fra-high nuclei imm
full234_simple_frai = data.frame(FRA1 = full234$FRA1, proximity_imm = full234$Imm_proximity_status234)
nfi20 = nrow(full234_simple_frai[full234_simple_frai$proximity_imm == "0 - 100",])
nfi30 = nrow(full234_simple_frai[full234_simple_frai$proximity_imm == "100 - 150",])
nfi40 = nrow(full234_simple_frai[full234_simple_frai$proximity_imm == "150 - 200",])
nfid = nrow(full234_simple_frai[full234_simple_frai$proximity_imm == "200 - ∞",])
nfi20_frahigh = nrow(full234_simple_frai[full234_simple_frai$proximity_imm == "0 - 100" & full234_simple_frai$FRA1 == 1,])
nfi30_frahigh = nrow(full234_simple_frai[full234_simple_frai$proximity_imm == "100 - 150" & full234_simple_frai$FRA1 == 1,])
nfi40_frahigh = nrow(full234_simple_frai[full234_simple_frai$proximity_imm == "150 - 200" & full234_simple_frai$FRA1 == 1,])
nfid_frahigh = nrow(full234_simple_frai[full234_simple_frai$proximity_imm == "200 - ∞" & full234_simple_frai$FRA1 == 1,])
nfis_frahigh = nrow(full234_simple_frai[full234_simple_frai$proximity_imm == "Stromal" & full234_simple_frai$FRA1 == 1,])
nfis = nrow(full234) - nfi20 - nfi30 - nfi40 - nfid
nfi_frahigh = nfi20_frahigh + nfi30_frahigh + nfi40_frahigh + nfid_frahigh
nfi_fralow = nfi20 + nfi30 + nfi40 + nfid - nfi_frahigh  

pi_frahigh_20 = dhyper(nfi20_frahigh, nfi_frahigh, nfi_fralow, nfi20)
pi_frahigh_30 = dhyper(nfi30_frahigh, nfi_frahigh, nfi_fralow, nfi30)
pi_frahigh_40 = dhyper(nfi40_frahigh, nfi_frahigh, nfi_fralow, nfi40)
pi_frahigh_distal = dhyper(nfid_frahigh, nfi_frahigh, nfi_fralow, nfid)

dhyper_EI_FRA = data.frame(Proximity = c("0 - 100", "100 - 150", "150 - 200", "200 - ∞"),
                           number_of_frahigh_drawn = c(nfi20_frahigh,nfi30_frahigh,nfi40_frahigh,nfid_frahigh),
                           number_of_frahigh = c(nfi_frahigh,nfi_frahigh,nfi_frahigh,nfi_frahigh),
                           number_of_fralow = c(nfi_fralow,nfi_fralow,nfi_fralow,nfi_fralow),
                           number_of_drawn = c(nfi20,nfi30,nfi40,nfid),
                           probability = c(pi_frahigh_20, pi_frahigh_30, pi_frahigh_40, pi_frahigh_distal))

#Visualize the idea that the probability of seeing fra-high epithelial nuclei increases as the distance to the immune cells decreases.
table_EI_fra1_distance = data.frame(Fra_low = c(nfi20 - nfi20_frahigh, nfi30 - nfi30_frahigh, nfi40 - nfi40_frahigh, nfid - nfid_frahigh), Fra_high = c(nfi20_frahigh, nfi30_frahigh, nfi40_frahigh, nfid_frahigh), row.names = c("0 - 100","100 - 150","150 - 200","200 - ∞"))

table_EI_fra1_distance$Distance_to_immune_cells_µm = row.names(table_EI_fra1_distance)
mtable_EI_fra1_distance = melt(table_EI_fra1_distance, id.vars = "Distance_to_immune_cells_µm")
bar_EI_fra1_distance = ggplot(mtable_EI_fra1_distance, aes(Distance_to_immune_cells_µm, value, fill = variable)) +
  geom_bar(position = "fill", stat = "identity")+
  ggtitle("Relationship between the probability of seeing a Fra-1 high epithelial nuclei and the distance of the nuclei to the immune cells")
bar_EI_fra1_distance

per_EI_frahigh_distance = data.frame (Fra_high_fraction = c(nfi20_frahigh / nfi20, nfi30_frahigh / nfi30, nfi40_frahigh / nfi40, nfid_frahigh / nfid),Distance_to_immune_cells_µm = c("0 - 100","100 - 150","150 - 200","200 - ∞"))
bar_per_EI_frahigh_distance = ggplot(per_EI_frahigh_distance, aes(x = Distance_to_immune_cells_µm, y = Fra_high_fraction)) +
  geom_bar(aes(fill = Distance_to_immune_cells_µm), position = "dodge", stat = "identity")+
  scale_fill_manual(values=c("red", "purple", "blue", "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle("Percentage of Fra-high epithelial nuclei in the epithelial nuclei with different distances to the immune cells")
bar_per_EI_frahigh_distance

EI_distance = data.frame (Number_of_nuclei = c(nfi20, nfi30, nfi40, nfid),Distance_to_immune_cells_µm = c("0 - 100","100 - 150","150 - 200","200 - ∞"))
bar_EI_distance = ggplot(EI_distance, aes(x = Distance_to_immune_cells_µm, y = Number_of_nuclei)) +
  geom_bar(aes(fill = Distance_to_immune_cells_µm), position = "dodge", stat = "identity")+
  scale_fill_manual(values=c("red", "purple", "blue", "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle("Number of epithelial nuclei with different distances to the immune cells")
bar_EI_distance

EI_frahigh_distance = data.frame (Number_of_nuclei = c(nfi20_frahigh, nfi30_frahigh, nfi40_frahigh, nfid_frahigh),Distance_to_immune_cells_µm = c("0 - 100","100 - 150","150 - 200","200 - ∞"))
bar_EI_frahigh_distance = ggplot(EI_frahigh_distance, aes(x = Distance_to_immune_cells_µm, y = Number_of_nuclei)) +
  geom_bar(aes(fill = Distance_to_immune_cells_µm), position = "dodge", stat = "identity")+
  scale_fill_manual(values=c("red", "purple", "blue", "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle("Number of fra-high epithelial nuclei with different distances to the immune cells")
bar_EI_frahigh_distance

################################################################################
#Hypergeometric statistics for yap-high nuclei imm
full234_simple_yapi = data.frame(YAP1 = full234$YAP1, proximity_imm = full234$Imm_proximity_status234)
nyi20 = nrow(full234_simple_yapi[full234_simple_yapi$proximity_imm == "0 - 100",])
nyi30 = nrow(full234_simple_yapi[full234_simple_yapi$proximity_imm == "100 - 150",])
nyi40 = nrow(full234_simple_yapi[full234_simple_yapi$proximity_imm == "150 - 200",])
nyid = nrow(full234_simple_yapi[full234_simple_yapi$proximity_imm == "200 - ∞",])
nyi20_yaphigh = nrow(full234_simple_yapi[full234_simple_yapi$proximity_imm == "0 - 100" & full234_simple_yapi$YAP1 == 1,])
nyi30_yaphigh = nrow(full234_simple_yapi[full234_simple_yapi$proximity_imm == "100 - 150" & full234_simple_yapi$YAP1 == 1,])
nyi40_yaphigh = nrow(full234_simple_yapi[full234_simple_yapi$proximity_imm == "150 - 200" & full234_simple_yapi$YAP1 == 1,])
nyid_yaphigh = nrow(full234_simple_yapi[full234_simple_yapi$proximity_imm == "200 - ∞" & full234_simple_yapi$YAP1 == 1,])
nyis_yaphigh = nrow(full234_simple_yapi[full234_simple_yapi$proximity_imm == "Stromal" & full234_simple_yapi$YAP1 == 1,])
nyis = nrow(full234) - nyi20 - nyi30 - nyi40 - nyid
nyi_yaphigh = nyi20_yaphigh + nyi30_yaphigh + nyi40_yaphigh + nyid_yaphigh
nyi_yaplow = nyi20 + nyi30 + nyi40 + nyid - nyi_yaphigh

pi_yaphigh_20 = dhyper(nyi20_yaphigh, nyi_yaphigh, nyi_yaplow, nyi20)
pi_yaphigh_30 = dhyper(nyi30_yaphigh, nyi_yaphigh, nyi_yaplow, nyi30)
pi_yaphigh_40 = dhyper(nyi40_yaphigh, nyi_yaphigh, nyi_yaplow, nyi40)
pi_yaphigh_distal = dhyper(nyid_yaphigh, nyi_yaphigh, nyi_yaplow, nyid)

dhyper_EI_YAP1 = data.frame(Proximity = c("0 - 100", "100 - 150", "150 - 200", "200 - ∞"),
                            number_of_yaphigh_drawn = c(nyi20_yaphigh,nyi30_yaphigh,nyi40_yaphigh,nyid_yaphigh),
                            number_of_yaphigh = c(nyi_yaphigh,nyi_yaphigh,nyi_yaphigh,nyi_yaphigh),
                            number_of_yaplow = c(nyi_yaplow,nyi_yaplow,nyi_yaplow,nyi_yaplow),
                            number_of_drawn = c(nyi20,nyi30,nyi40,nyid),
                            probability = c(pi_yaphigh_20, pi_yaphigh_30, pi_yaphigh_40, pi_yaphigh_distal))

#Visualize the idea that the probability of seeing yap-high epithelial nuclei increases as the distance to the immune cells decreases.
table_EI_yap_distance = data.frame(YAP1_low = c(nyi20 - nyi20_yaphigh, nyi30 - nyi30_yaphigh, nyi40 - nyi40_yaphigh, nyid - nyid_yaphigh), YAP1_high = c(nyi20_yaphigh, nyi30_yaphigh, nyi40_yaphigh, nyid_yaphigh), row.names = c("0 - 100","100 - 150","150 - 200","200 - ∞"))

table_EI_yap_distance$Distance_to_immune_cells_µm = row.names(table_EI_yap_distance)
mtable_EI_yap_distance = melt(table_EI_yap_distance, id.vars = "Distance_to_immune_cells_µm")
bar_EI_yap_distance = ggplot(mtable_EI_yap_distance, aes(Distance_to_immune_cells_µm, value, fill = variable)) +
  geom_bar(position = "fill", stat = "identity")+
  ggtitle("Relationship between the probability of seeing a YAP1-1 high epithelial nuclei and the distance of the nuclei to the immune cells")
bar_EI_yap_distance

per_EI_yaphigh_distance = data.frame (YAP1_high_fraction = c(nyi20_yaphigh / nyi20, nyi30_yaphigh / nyi30, nyi40_yaphigh / nyi40, nyid_yaphigh / nyid),Distance_to_immune_cells_µm = c("0 - 100","100 - 150","150 - 200","200 - ∞"))
bar_per_EI_yaphigh_distance = ggplot(per_EI_yaphigh_distance, aes(x = Distance_to_immune_cells_µm, y = YAP1_high_fraction)) +
  geom_bar(aes(fill = Distance_to_immune_cells_µm), position = "dodge", stat = "identity")+
  scale_fill_manual(values=c("red", "purple", "blue", "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle("Percentage of YAP1-high epithelial nuclei in the epithelial nuclei with different distances to the immune cells")
bar_per_EI_yaphigh_distance

EI_distance = data.frame (Number_of_nuclei = c(nyi20, nyi30, nyi40, nyid),Distance_to_immune_cells_µm = c("0 - 100","100 - 150","150 - 200","200 - ∞"))
bar_EI_distance = ggplot(EI_distance, aes(x = Distance_to_immune_cells_µm, y = Number_of_nuclei)) +
  geom_bar(aes(fill = Distance_to_immune_cells_µm), position = "dodge", stat = "identity")+
  scale_fill_manual(values=c("red", "purple", "blue", "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle("Number of epithelial nuclei with different distances to the immune cells")
bar_EI_distance

EI_yaphigh_distance = data.frame (Number_of_nuclei = c(nyi20_yaphigh, nyi30_yaphigh, nyi40_yaphigh, nyid_yaphigh),Distance_to_immune_cells_µm = c("0 - 100","100 - 150","150 - 200","200 - ∞"))
bar_EI_yaphigh_distance = ggplot(EI_yaphigh_distance, aes(x = Distance_to_immune_cells_µm, y = Number_of_nuclei)) +
  geom_bar(aes(fill = Distance_to_immune_cells_µm), position = "dodge", stat = "identity")+
  scale_fill_manual(values=c("red", "purple", "blue", "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle("Number of yap-high epithelial nuclei with different distances to the immune cells")
bar_EI_yaphigh_distance
################################################################################
################################################################################
################################################################################
#Hypergeometric statistics for parp-high nuclei T
T234_simple_parp = data.frame(PARP = T234$PARP, proximity_T = T234$T_proximity_status234)
npt20 = nrow(T234_simple_parp[T234_simple_parp$proximity_T == "0 - 100",])
npt30 = nrow(T234_simple_parp[T234_simple_parp$proximity_T == "100 - 150",])
npt40 = nrow(T234_simple_parp[T234_simple_parp$proximity_T == "150 - 200",])
nptd = nrow(T234_simple_parp[T234_simple_parp$proximity_T == "200 - ∞",])
npt20_parphigh = nrow(T234_simple_parp[T234_simple_parp$proximity_T == "0 - 100" & T234_simple_parp$PARP == 1,])
npt30_parphigh = nrow(T234_simple_parp[T234_simple_parp$proximity_T == "100 - 150" & T234_simple_parp$PARP == 1,])
npt40_parphigh = nrow(T234_simple_parp[T234_simple_parp$proximity_T == "150 - 200" & T234_simple_parp$PARP == 1,])
nptd_parphigh = nrow(T234_simple_parp[T234_simple_parp$proximity_T == "200 - ∞" & T234_simple_parp$PARP == 1,])
npt_parphigh = npt20_parphigh + npt30_parphigh + npt40_parphigh + nptd_parphigh
npt_parplow = npt20 + npt30 + npt40 + nptd - npt_parphigh  

pt_parphigh_20 = dhyper(npt20_parphigh, npt_parphigh, npt_parplow, npt20)
pt_parphigh_30 = dhyper(npt30_parphigh, npt_parphigh, npt_parplow, npt30)
pt_parphigh_40 = dhyper(npt40_parphigh, npt_parphigh, npt_parplow, npt40)
pt_parphigh_distal = dhyper(nptd_parphigh, npt_parphigh, npt_parplow, nptd)

dhyper_ET_PARP = data.frame(Proximity = c("0 - 100", "100 - 150", "150 - 200", "200 - ∞"),
                            number_of_parphigh_drawn = c(npt20_parphigh,npt30_parphigh,npt40_parphigh,nptd_parphigh),
                            number_of_parphigh = c(npt_parphigh,npt_parphigh,npt_parphigh,npt_parphigh),
                            number_of_parplow = c(npt_parplow,npt_parplow,npt_parplow,npt_parplow),
                            number_of_drawn = c(npt20,npt30,npt40,nptd),
                            probability = c(pt_parphigh_20, pt_parphigh_30, pt_parphigh_40, pt_parphigh_distal))

#Visualize the idea that the probability of seeing parp-high epithelial nuclei increases as the distance to the T cells decreases.
table_ET_parp_distance = data.frame(PARP_low = c(npt20 - npt20_parphigh, npt30 - npt30_parphigh, npt40 - npt40_parphigh, nptd - nptd_parphigh), PARP_high = c(npt20_parphigh, npt30_parphigh, npt40_parphigh, nptd_parphigh), row.names = c("0 - 100","100 - 150","150 - 200","200 - ∞"))

table_ET_parp_distance$Distance_to_T_cells_µm = row.names(table_ET_parp_distance)
mtable_ET_parp_distance = melt(table_ET_parp_distance, id.vars = "Distance_to_T_cells_µm")
bar_ET_parp_distance = ggplot(mtable_ET_parp_distance, aes(Distance_to_T_cells_µm, value, fill = variable)) +
  geom_bar(position = "fill", stat = "identity")+
  ggtitle("Relationship between the probability of seeing a PARP-1 high epithelial nuclei and the distance of the nuclei to the T cells")
bar_ET_parp_distance

per_ET_parphigh_distance = data.frame (PARP_high_fraction = c(npt20_parphigh / npt20, npt30_parphigh / npt30, npt40_parphigh / npt40, nptd_parphigh / nptd),Distance_to_T_cells_µm = c("0 - 100","100 - 150","150 - 200","200 - ∞"))
bar_per_ET_parphigh_distance = ggplot(per_ET_parphigh_distance, aes(x = Distance_to_T_cells_µm, y = PARP_high_fraction)) +
  geom_bar(aes(fill = Distance_to_T_cells_µm), position = "dodge", stat = "identity")+
  scale_fill_manual(values=c("red", "purple", "blue", "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle("Percentage of PARP-high epithelial nuclei in the epithelial nuclei with different distances to the T cells")
bar_per_ET_parphigh_distance

ET_distance = data.frame (Number_of_nuclei = c(npt20, npt30, npt40, nptd),Distance_to_T_cells_µm = c("0 - 100","100 - 150","150 - 200","200 - ∞"))
bar_ET_distance = ggplot(ET_distance, aes(x = Distance_to_T_cells_µm, y = Number_of_nuclei)) +
  geom_bar(aes(fill = Distance_to_T_cells_µm), position = "dodge", stat = "identity")+
  scale_fill_manual(values=c("red", "purple", "blue", "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle("Number of epithelial nuclei with different distances to the T cells")
bar_ET_distance

ET_parphigh_distance = data.frame (Number_of_nuclei = c(npt20_parphigh, npt30_parphigh, npt40_parphigh, nptd_parphigh),Distance_to_T_cells_µm = c("0 - 100","100 - 150","150 - 200","200 - ∞"))
bar_ET_parphigh_distance = ggplot(ET_parphigh_distance, aes(x = Distance_to_T_cells_µm, y = Number_of_nuclei)) +
  geom_bar(aes(fill = Distance_to_T_cells_µm), position = "dodge", stat = "identity")+
  scale_fill_manual(values=c("red", "purple", "blue", "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle("Number of parp-high epithelial nuclei with different distances to the T cells")
bar_ET_parphigh_distance


################################################################################
#Hypergeometric statistics for fra-high nuclei T
T234_simple_fra = data.frame(FRA1 = T234$FRA1, proximity_T = T234$T_proximity_status234)
nft20 = nrow(T234_simple_fra[T234_simple_fra$proximity_T == "0 - 100",])
nft30 = nrow(T234_simple_fra[T234_simple_fra$proximity_T == "100 - 150",])
nft40 = nrow(T234_simple_fra[T234_simple_fra$proximity_T == "150 - 200",])
nftd = nrow(T234_simple_fra[T234_simple_fra$proximity_T == "200 - ∞",])
nft20_frahigh = nrow(T234_simple_fra[T234_simple_fra$proximity_T == "0 - 100" & T234_simple_fra$FRA1 == 1,])
nft30_frahigh = nrow(T234_simple_fra[T234_simple_fra$proximity_T == "100 - 150" & T234_simple_fra$FRA1 == 1,])
nft40_frahigh = nrow(T234_simple_fra[T234_simple_fra$proximity_T == "150 - 200" & T234_simple_fra$FRA1 == 1,])
nftd_frahigh = nrow(T234_simple_fra[T234_simple_fra$proximity_T == "200 - ∞" & T234_simple_fra$FRA1 == 1,])
nft_frahigh = nft20_frahigh + nft30_frahigh + nft40_frahigh + nftd_frahigh
nft_fralow = nft20 + nft30 + nft40 + nftd - nft_frahigh  

pt_frahigh_20 = dhyper(nft20_frahigh, nft_frahigh, nft_fralow, nft20)
pt_frahigh_30 = dhyper(nft30_frahigh, nft_frahigh, nft_fralow, nft30)
pt_frahigh_40 = dhyper(nft40_frahigh, nft_frahigh, nft_fralow, nft40)
pt_frahigh_distal = dhyper(nftd_frahigh, nft_frahigh, nft_fralow, nftd)

dhyper_ET_FRA = data.frame(Proximity = c("0 - 100", "100 - 150", "150 - 200", "200 - ∞"),
                           number_of_frahigh_drawn = c(nft20_frahigh,nft30_frahigh,nft40_frahigh,nftd_frahigh),
                           number_of_frahigh = c(nft_frahigh,nft_frahigh,nft_frahigh,nft_frahigh),
                           number_of_fralow = c(nft_fralow,nft_fralow,nft_fralow,nft_fralow),
                           number_of_drawn = c(nft20,nft30,nft40,nftd),
                           probability = c(pt_frahigh_20, pt_frahigh_30, pt_frahigh_40, pt_frahigh_distal))

#Visualize the idea that the probability of seeing fra-high epithelial nuclei increases as the distance to the T cells decreases.
table_ET_fra1_distance = data.frame(Fra_low = c(nft20 - nft20_frahigh, nft30 - nft30_frahigh, nft40 - nft40_frahigh, nftd - nftd_frahigh), Fra_high = c(nft20_frahigh, nft30_frahigh, nft40_frahigh, nftd_frahigh), row.names = c("0 - 100","100 - 150","150 - 200","200 - ∞"))

table_ET_fra1_distance$Distance_to_T_cells_µm = row.names(table_ET_fra1_distance)
mtable_ET_fra1_distance = melt(table_ET_fra1_distance, id.vars = "Distance_to_T_cells_µm")
bar_ET_fra1_distance = ggplot(mtable_ET_fra1_distance, aes(Distance_to_T_cells_µm, value, fill = variable)) +
  geom_bar(position = "fill", stat = "identity")+
  ggtitle("Relationship between the probability of seeing a Fra-1 high epithelial nuclei and the distance of the nuclei to the T cells")
bar_ET_fra1_distance

per_ET_frahigh_distance = data.frame (Fra_high_fraction = c(nft20_frahigh / nft20, nft30_frahigh / nft30, nft40_frahigh / nft40, nftd_frahigh / nftd),Distance_to_T_cells_µm = c("0 - 100","100 - 150","150 - 200","200 - ∞"))
bar_per_ET_frahigh_distance = ggplot(per_ET_frahigh_distance, aes(x = Distance_to_T_cells_µm, y = Fra_high_fraction)) +
  geom_bar(aes(fill = Distance_to_T_cells_µm), position = "dodge", stat = "identity")+
  scale_fill_manual(values=c("red", "purple", "blue", "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle("Percentage of Fra-high epithelial nuclei in the epithelial nuclei with different distances to the T cells")
bar_per_ET_frahigh_distance

ET_distance = data.frame (Number_of_nuclei = c(nft20, nft30, nft40, nftd),Distance_to_T_cells_µm = c("0 - 100","100 - 150","150 - 200","200 - ∞"))
bar_ET_distance = ggplot(ET_distance, aes(x = Distance_to_T_cells_µm, y = Number_of_nuclei)) +
  geom_bar(aes(fill = Distance_to_T_cells_µm), position = "dodge", stat = "identity")+
  scale_fill_manual(values=c("red", "purple", "blue", "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle("Number of epithelial nuclei with different distances to the T cells")
bar_ET_distance

ET_frahigh_distance = data.frame (Number_of_nuclei = c(nft20_frahigh, nft30_frahigh, nft40_frahigh, nftd_frahigh),Distance_to_T_cells_µm = c("0 - 100","100 - 150","150 - 200","200 - ∞"))
bar_ET_frahigh_distance = ggplot(ET_frahigh_distance, aes(x = Distance_to_T_cells_µm, y = Number_of_nuclei)) +
  geom_bar(aes(fill = Distance_to_T_cells_µm), position = "dodge", stat = "identity")+
  scale_fill_manual(values=c("red", "purple", "blue", "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle("Number of fra-high epithelial nuclei with different distances to the T cells")
bar_ET_frahigh_distance

################################################################################
#Hypergeometric statistics for yap-high nuclei T
T234_simple_yap = data.frame(YAP1 = T234$YAP1, proximity_T = T234$T_proximity_status234)
nyt20 = nrow(T234_simple_yap[T234_simple_yap$proximity_T == "0 - 100",])
nyt30 = nrow(T234_simple_yap[T234_simple_yap$proximity_T == "100 - 150",])
nyt40 = nrow(T234_simple_yap[T234_simple_yap$proximity_T == "150 - 200",])
nytd = nrow(T234_simple_yap[T234_simple_yap$proximity_T == "200 - ∞",])
nyt20_yaphigh = nrow(T234_simple_yap[T234_simple_yap$proximity_T == "0 - 100" & T234_simple_yap$YAP1 == 1,])
nyt30_yaphigh = nrow(T234_simple_yap[T234_simple_yap$proximity_T == "100 - 150" & T234_simple_yap$YAP1 == 1,])
nyt40_yaphigh = nrow(T234_simple_yap[T234_simple_yap$proximity_T == "150 - 200" & T234_simple_yap$YAP1 == 1,])
nytd_yaphigh = nrow(T234_simple_yap[T234_simple_yap$proximity_T == "200 - ∞" & T234_simple_yap$YAP1 == 1,])
nyt_yaphigh = nyt20_yaphigh + nyt30_yaphigh + nyt40_yaphigh + nytd_yaphigh
nyt_yaplow = nyt20 + nyt30 + nyt40 + nytd - nyt_yaphigh  

pt_yaphigh_20 = dhyper(nyt20_yaphigh, nyt_yaphigh, nyt_yaplow, nyt20)
pt_yaphigh_30 = dhyper(nyt30_yaphigh, nyt_yaphigh, nyt_yaplow, nyt30)
pt_yaphigh_40 = dhyper(nyt40_yaphigh, nyt_yaphigh, nyt_yaplow, nyt40)
pt_yaphigh_distal = dhyper(nytd_yaphigh, nyt_yaphigh, nyt_yaplow, nytd)

dhyper_ET_YAP1 = data.frame(Proximity = c("0 - 100", "100 - 150", "150 - 200", "200 - ∞"),
                            number_of_yaphigh_drawn = c(nyt20_yaphigh,nyt30_yaphigh,nyt40_yaphigh,nytd_yaphigh),
                            number_of_yaphigh = c(nyt_yaphigh,nyt_yaphigh,nyt_yaphigh,nyt_yaphigh),
                            number_of_yaplow = c(nyt_yaplow,nyt_yaplow,nyt_yaplow,nyt_yaplow),
                            number_of_drawn = c(nyt20,nyt30,nyt40,nytd),
                            probability = c(pt_yaphigh_20, pt_yaphigh_30, pt_yaphigh_40, pt_yaphigh_distal))

#Visualize the idea that the probability of seeing yap-high epithelial nuclei increases as the distance to the T cells decreases.
table_ET_yap_distance = data.frame(YAP1_low = c(nyt20 - nyt20_yaphigh, nyt30 - nyt30_yaphigh, nyt40 - nyt40_yaphigh, nytd - nytd_yaphigh), YAP1_high = c(nyt20_yaphigh, nyt30_yaphigh, nyt40_yaphigh, nytd_yaphigh), row.names = c("0 - 100","100 - 150","150 - 200","200 - ∞"))

table_ET_yap_distance$Distance_to_T_cells_µm = row.names(table_ET_yap_distance)
mtable_ET_yap_distance = melt(table_ET_yap_distance, id.vars = "Distance_to_T_cells_µm")
bar_ET_yap_distance = ggplot(mtable_ET_yap_distance, aes(Distance_to_T_cells_µm, value, fill = variable)) +
  geom_bar(position = "fill", stat = "identity")+
  ggtitle("Relationship between the probability of seeing a YAP1-1 high epithelial nuclei and the distance of the nuclei to the T cells")
bar_ET_yap_distance

per_ET_yaphigh_distance = data.frame (YAP1_high_fraction = c(nyt20_yaphigh / nyt20, nyt30_yaphigh / nyt30, nyt40_yaphigh / nyt40, nytd_yaphigh / nytd),Distance_to_T_cells_µm = c("0 - 100","100 - 150","150 - 200","200 - ∞"))
bar_per_ET_yaphigh_distance = ggplot(per_ET_yaphigh_distance, aes(x = Distance_to_T_cells_µm, y = YAP1_high_fraction)) +
  geom_bar(aes(fill = Distance_to_T_cells_µm), position = "dodge", stat = "identity")+
  scale_fill_manual(values=c("red", "purple", "blue", "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle("Percentage of YAP1-high epithelial nuclei in the epithelial nuclei with different distances to the T cells")
bar_per_ET_yaphigh_distance

ET_distance = data.frame (Number_of_nuclei = c(nyt20, nyt30, nyt40, nytd),Distance_to_T_cells_µm = c("0 - 100","100 - 150","150 - 200","200 - ∞"))
bar_ET_distance = ggplot(ET_distance, aes(x = Distance_to_T_cells_µm, y = Number_of_nuclei)) +
  geom_bar(aes(fill = Distance_to_T_cells_µm), position = "dodge", stat = "identity")+
  scale_fill_manual(values=c("red", "purple", "blue", "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle("Number of epithelial nuclei with different distances to the T cells")
bar_ET_distance

ET_yaphigh_distance = data.frame (Number_of_nuclei = c(nyt20_yaphigh, nyt30_yaphigh, nyt40_yaphigh, nytd_yaphigh),Distance_to_T_cells_µm = c("0 - 100","100 - 150","150 - 200","200 - ∞"))
bar_ET_yaphigh_distance = ggplot(ET_yaphigh_distance, aes(x = Distance_to_T_cells_µm, y = Number_of_nuclei)) +
  geom_bar(aes(fill = Distance_to_T_cells_µm), position = "dodge", stat = "identity")+
  scale_fill_manual(values=c("red", "purple", "blue", "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle("Number of yap-high epithelial nuclei with different distances to the T cells")
bar_ET_yaphigh_distance
################################################################################
#Hypergeometric statistics for parp-high nuclei CD8
CD8_T234_simple_parp = data.frame(PARP = CD8_T234$PARP, proximity_CD8 = CD8_T234$CD8_proximity_status234)
np820 = nrow(CD8_T234_simple_parp[CD8_T234_simple_parp$proximity_CD8 == "0 - 100",])
np830 = nrow(CD8_T234_simple_parp[CD8_T234_simple_parp$proximity_CD8 == "100 - 150",])
np840 = nrow(CD8_T234_simple_parp[CD8_T234_simple_parp$proximity_CD8 == "150 - 200",])
np8d = nrow(CD8_T234_simple_parp[CD8_T234_simple_parp$proximity_CD8 == "200 - ∞",])
np820_parphigh = nrow(CD8_T234_simple_parp[CD8_T234_simple_parp$proximity_CD8 == "0 - 100" & CD8_T234_simple_parp$PARP == 1,])
np830_parphigh = nrow(CD8_T234_simple_parp[CD8_T234_simple_parp$proximity_CD8 == "100 - 150" & CD8_T234_simple_parp$PARP == 1,])
np840_parphigh = nrow(CD8_T234_simple_parp[CD8_T234_simple_parp$proximity_CD8 == "150 - 200" & CD8_T234_simple_parp$PARP == 1,])
np8d_parphigh = nrow(CD8_T234_simple_parp[CD8_T234_simple_parp$proximity_CD8 == "200 - ∞" & CD8_T234_simple_parp$PARP == 1,])
np8_parphigh = np820_parphigh + np830_parphigh + np840_parphigh + np8d_parphigh
np8_parplow = np820 + np830 + np840 + np8d - np8_parphigh  

p8_parphigh_20 = dhyper(np820_parphigh, np8_parphigh, np8_parplow, np820)
p8_parphigh_30 = dhyper(np830_parphigh, np8_parphigh, np8_parplow, np830)
p8_parphigh_40 = dhyper(np840_parphigh, np8_parphigh, np8_parplow, np840)
p8_parphigh_distal = dhyper(np8d_parphigh, np8_parphigh, np8_parplow, np8d)

dhyper_E8_PARP = data.frame(Proximity = c("0 - 100", "100 - 150", "150 - 200", "200 - ∞"),
                            number_of_parphigh_drawn = c(np820_parphigh,np830_parphigh,np840_parphigh,np8d_parphigh),
                            number_of_parphigh = c(np8_parphigh,np8_parphigh,np8_parphigh,np8_parphigh),
                            number_of_parplow = c(np8_parplow,np8_parplow,np8_parplow,np8_parplow),
                            number_of_drawn = c(np820,np830,np840,np8d),
                            probability = c(p8_parphigh_20, p8_parphigh_30, p8_parphigh_40, p8_parphigh_distal))
e8p20 = c(np820_parphigh,(np820 - np820_parphigh))
e8p30 = c(np830_parphigh,(np830 - np830_parphigh))
e8p40 = c(np840_parphigh,(np840 - np840_parphigh))
e8pd = c(np8d_parphigh,(np8d - np8d_parphigh))
e8p = data.frame(e8p20, e8p30, e8p40, e8pd)
e8pc = chisq.test(e8p)


#Visualize the idea that the probability of seeing parp-high epithelial nuclei increases as the distance to the CD8+ T cells decreases.
table_E8_parp_distance = data.frame(PARP_low = c(np820 - np820_parphigh, np830 - np830_parphigh, np840 - np840_parphigh, np8d - np8d_parphigh), PARP_high = c(np820_parphigh, np830_parphigh, np840_parphigh, np8d_parphigh), row.names = c("0 - 100","100 - 150","150 - 200","200 - ∞"))

table_E8_parp_distance$Distance_to_CD8_T_cells_µm = row.names(table_E8_parp_distance)
mtable_E8_parp_distance = melt(table_E8_parp_distance, id.vars = "Distance_to_CD8_T_cells_µm")
bar_E8_parp_distance = ggplot(mtable_E8_parp_distance, aes(Distance_to_CD8_T_cells_µm, value, fill = variable)) +
  geom_bar(position = "fill", stat = "identity")+
  ggtitle("Relationship between the probability of seeing a PARP-1 high epithelial nuclei and the distance of the nuclei to the CD8+ T cells")
bar_E8_parp_distance

per_E8_parphigh_distance = data.frame (PARP_high_fraction = c(np820_parphigh / np820, np830_parphigh / np830, np840_parphigh / np840, np8d_parphigh / np8d),Distance_to_CD8_T_cells_µm = c("0 - 100","100 - 150","150 - 200","200 - ∞"))
bar_per_E8_parphigh_distance = ggplot(per_E8_parphigh_distance, aes(x = Distance_to_CD8_T_cells_µm, y = PARP_high_fraction)) +
  geom_bar(aes(fill = Distance_to_CD8_T_cells_µm), position = "dodge", stat = "identity")+
  scale_fill_manual(values=c("red", "purple", "blue", "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
bar_per_E8_parphigh_distance

E8_distance = data.frame (Number_of_nuclei = c(np820, np830, np840, np8d),Distance_to_CD8_T_cells_µm = c("0 - 100","100 - 150","150 - 200","200 - ∞"))
bar_E8_distance = ggplot(E8_distance, aes(x = Distance_to_CD8_T_cells_µm, y = Number_of_nuclei)) +
  geom_bar(aes(fill = Distance_to_CD8_T_cells_µm), position = "dodge", stat = "identity")+
  scale_fill_manual(values=c("red", "purple", "blue", "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle("Number of epithelial nuclei with different distances to the CD8+ T cells")
bar_E8_distance

E8_parphigh_distance = data.frame (Number_of_nuclei = c(np820_parphigh, np830_parphigh, np840_parphigh, np8d_parphigh),Distance_to_CD8_T_cells_µm = c("0 - 100","100 - 150","150 - 200","200 - ∞"))
bar_E8_parphigh_distance = ggplot(E8_parphigh_distance, aes(x = Distance_to_CD8_T_cells_µm, y = Number_of_nuclei)) +
  geom_bar(aes(fill = Distance_to_CD8_T_cells_µm), position = "dodge", stat = "identity")+
  scale_fill_manual(values=c("red", "purple", "blue", "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle("Number of parp-high epithelial nuclei with different distances to the CD8+ T cells")
bar_E8_parphigh_distance


################################################################################
#Hypergeometric statistics for fra-high nuclei CD8
CD8_T234_simple_fra = data.frame(FRA1 = CD8_T234$FRA1, proximity_CD8 = CD8_T234$CD8_proximity_status234)
nf820 = nrow(CD8_T234_simple_fra[CD8_T234_simple_fra$proximity_CD8 == "0 - 100",])
nf830 = nrow(CD8_T234_simple_fra[CD8_T234_simple_fra$proximity_CD8 == "100 - 150",])
nf840 = nrow(CD8_T234_simple_fra[CD8_T234_simple_fra$proximity_CD8 == "150 - 200",])
nf8d = nrow(CD8_T234_simple_fra[CD8_T234_simple_fra$proximity_CD8 == "200 - ∞",])
nf820_frahigh = nrow(CD8_T234_simple_fra[CD8_T234_simple_fra$proximity_CD8 == "0 - 100" & CD8_T234_simple_fra$FRA1 == 1,])
nf830_frahigh = nrow(CD8_T234_simple_fra[CD8_T234_simple_fra$proximity_CD8 == "100 - 150" & CD8_T234_simple_fra$FRA1 == 1,])
nf840_frahigh = nrow(CD8_T234_simple_fra[CD8_T234_simple_fra$proximity_CD8 == "150 - 200" & CD8_T234_simple_fra$FRA1 == 1,])
nf8d_frahigh = nrow(CD8_T234_simple_fra[CD8_T234_simple_fra$proximity_CD8 == "200 - ∞" & CD8_T234_simple_fra$FRA1 == 1,])
nf8_frahigh = nf820_frahigh + nf830_frahigh + nf840_frahigh + nf8d_frahigh
nf8_fralow = nf820 + nf830 + nf840 + nf8d - nf8_frahigh  

p8_frahigh_20 = dhyper(nf820_frahigh, nf8_frahigh, nf8_fralow, nf820)
p8_frahigh_30 = dhyper(nf830_frahigh, nf8_frahigh, nf8_fralow, nf830)
p8_frahigh_40 = dhyper(nf840_frahigh, nf8_frahigh, nf8_fralow, nf840)
p8_frahigh_distal = dhyper(nf8d_frahigh, nf8_frahigh, nf8_fralow, nf8d)

dhyper_E8_FRA = data.frame(Proximity = c("0 - 100", "100 - 150", "150 - 200", "200 - ∞"),
                           number_of_frahigh_drawn = c(nf820_frahigh,nf830_frahigh,nf840_frahigh,nf8d_frahigh),
                           number_of_frahigh = c(nf8_frahigh,nf8_frahigh,nf8_frahigh,nf8_frahigh),
                           number_of_fralow = c(nf8_fralow,nf8_fralow,nf8_fralow,nf8_fralow),
                           number_of_drawn = c(nf820,nf830,nf840,nf8d),
                           probability = c(p8_frahigh_20, p8_frahigh_30, p8_frahigh_40, p8_frahigh_distal))

e8f20 = c(nf820_frahigh,(nf820 - nf820_frahigh))
e8f30 = c(nf830_frahigh,(nf830 - nf830_frahigh))
e8f40 = c(nf840_frahigh,(nf840 - nf840_frahigh))
e8fd = c(nf8d_frahigh,(nf8d - nf8d_frahigh))
e8f = data.frame(e8f20, e8f30, e8f40, e8fd)
e8fc = chisq.test(e8f)


#Visualize the idea that the probability of seeing fra-high epithelial nuclei increases as the distance to the CD8+ T cells decreases.
table_E8_fra1_distance = data.frame(Fra_low = c(nf820 - nf820_frahigh, nf830 - nf830_frahigh, nf840 - nf840_frahigh, nf8d - nf8d_frahigh), Fra_high = c(nf820_frahigh, nf830_frahigh, nf840_frahigh, nf8d_frahigh), row.names = c("0 - 100","100 - 150","150 - 200","200 - ∞"))

table_E8_fra1_distance$Distance_to_CD8_T_cells_µm = row.names(table_E8_fra1_distance)
mtable_E8_fra1_distance = melt(table_E8_fra1_distance, id.vars = "Distance_to_CD8_T_cells_µm")
bar_E8_fra1_distance = ggplot(mtable_E8_fra1_distance, aes(Distance_to_CD8_T_cells_µm, value, fill = variable)) +
  geom_bar(position = "fill", stat = "identity")+
  ggtitle("Relationship between the probability of seeing a Fra-1 high epithelial nuclei and the distance of the nuclei to the CD8+ T cells")
bar_E8_fra1_distance

per_E8_frahigh_distance = data.frame (Fra_high_fraction = c(nf820_frahigh / nf820, nf830_frahigh / nf830, nf840_frahigh / nf840, nf8d_frahigh / nf8d),Distance_to_CD8_T_cells_µm = c("0 - 100","100 - 150","150 - 200","200 - ∞"))
bar_per_E8_frahigh_distance = ggplot(per_E8_frahigh_distance, aes(x = Distance_to_CD8_T_cells_µm, y = Fra_high_fraction)) +
  geom_bar(aes(fill = Distance_to_CD8_T_cells_µm), position = "dodge", stat = "identity")+
  scale_fill_manual(values=c("red", "purple", "blue", "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
bar_per_E8_frahigh_distance

E8_distance = data.frame (Number_of_nuclei = c(nf820, nf830, nf840, nf8d),Distance_to_CD8_T_cells_µm = c("0 - 100","100 - 150","150 - 200","200 - ∞"))
bar_E8_distance = ggplot(E8_distance, aes(x = Distance_to_CD8_T_cells_µm, y = Number_of_nuclei)) +
  geom_bar(aes(fill = Distance_to_CD8_T_cells_µm), position = "dodge", stat = "identity")+
  scale_fill_manual(values=c("red", "purple", "blue", "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle("Number of epithelial nuclei with different distances to the CD8+ T cells")
bar_E8_distance

E8_frahigh_distance = data.frame (Number_of_nuclei = c(nf820_frahigh, nf830_frahigh, nf840_frahigh, nf8d_frahigh),Distance_to_CD8_T_cells_µm = c("0 - 100","100 - 150","150 - 200","200 - ∞"))
bar_E8_frahigh_distance = ggplot(E8_frahigh_distance, aes(x = Distance_to_CD8_T_cells_µm, y = Number_of_nuclei)) +
  geom_bar(aes(fill = Distance_to_CD8_T_cells_µm), position = "dodge", stat = "identity")+
  scale_fill_manual(values=c("red", "purple", "blue", "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle("Number of fra-high epithelial nuclei with different distances to the CD8+ T cells")
bar_E8_frahigh_distance

################################################################################
#Hypergeometric statistics for yap-high nuclei CD8
CD8_T234_simple_yap = data.frame(YAP1 = CD8_T234$YAP1, proximity_CD8 = CD8_T234$CD8_proximity_status234)
ny820 = nrow(CD8_T234_simple_yap[CD8_T234_simple_yap$proximity_CD8 == "0 - 100",])
ny830 = nrow(CD8_T234_simple_yap[CD8_T234_simple_yap$proximity_CD8 == "100 - 150",])
ny840 = nrow(CD8_T234_simple_yap[CD8_T234_simple_yap$proximity_CD8 == "150 - 200",])
ny8d = nrow(CD8_T234_simple_yap[CD8_T234_simple_yap$proximity_CD8 == "200 - ∞",])
ny820_yaphigh = nrow(CD8_T234_simple_yap[CD8_T234_simple_yap$proximity_CD8 == "0 - 100" & CD8_T234_simple_yap$YAP1 == 1,])
ny830_yaphigh = nrow(CD8_T234_simple_yap[CD8_T234_simple_yap$proximity_CD8 == "100 - 150" & CD8_T234_simple_yap$YAP1 == 1,])
ny840_yaphigh = nrow(CD8_T234_simple_yap[CD8_T234_simple_yap$proximity_CD8 == "150 - 200" & CD8_T234_simple_yap$YAP1 == 1,])
ny8d_yaphigh = nrow(CD8_T234_simple_yap[CD8_T234_simple_yap$proximity_CD8 == "200 - ∞" & CD8_T234_simple_yap$YAP1 == 1,])
ny8_yaphigh = ny820_yaphigh + ny830_yaphigh + ny840_yaphigh + ny8d_yaphigh
ny8_yaplow = ny820 + ny830 + ny840 + ny8d - ny8_yaphigh  

p8_yaphigh_20 = dhyper(ny820_yaphigh, ny8_yaphigh, ny8_yaplow, ny820)
p8_yaphigh_30 = dhyper(ny830_yaphigh, ny8_yaphigh, ny8_yaplow, ny830)
p8_yaphigh_40 = dhyper(ny840_yaphigh, ny8_yaphigh, ny8_yaplow, ny840)
p8_yaphigh_distal = dhyper(ny8d_yaphigh, ny8_yaphigh, ny8_yaplow, ny8d)

dhyper_E8_YAP1 = data.frame(Proximity = c("0 - 100", "100 - 150", "150 - 200", "200 - ∞"),
                            number_of_yaphigh_drawn = c(ny820_yaphigh,ny830_yaphigh,ny840_yaphigh,ny8d_yaphigh),
                            number_of_yaphigh = c(ny8_yaphigh,ny8_yaphigh,ny8_yaphigh,ny8_yaphigh),
                            number_of_yaplow = c(ny8_yaplow,ny8_yaplow,ny8_yaplow,ny8_yaplow),
                            number_of_drawn = c(ny820,ny830,ny840,ny8d),
                            probability = c(p8_yaphigh_20, p8_yaphigh_30, p8_yaphigh_40, p8_yaphigh_distal))


e8y20 = c(ny820_yaphigh,(ny820 - ny820_yaphigh))
e8y30 = c(ny830_yaphigh,(ny830 - ny830_yaphigh))
e8y40 = c(ny840_yaphigh,(ny840 - ny840_yaphigh))
e8yd = c(ny8d_yaphigh,(ny8d - ny8d_yaphigh))
e8y = data.frame(e8y20, e8y30, e8y40, e8yd)
e8yc = chisq.test(e8y)

#Visualize the idea that the probability of seeing yap-high epithelial nuclei increases as the distance to the CD8+ T cells decreases.
table_E8_yap_distance = data.frame(YAP1_low = c(ny820 - ny820_yaphigh, ny830 - ny830_yaphigh, ny840 - ny840_yaphigh, ny8d - ny8d_yaphigh), YAP1_high = c(ny820_yaphigh, ny830_yaphigh, ny840_yaphigh, ny8d_yaphigh), row.names = c("0 - 100","100 - 150","150 - 200","200 - ∞"))

table_E8_yap_distance$Distance_to_CD8_T_cells_µm = row.names(table_E8_yap_distance)
mtable_E8_yap_distance = melt(table_E8_yap_distance, id.vars = "Distance_to_CD8_T_cells_µm")
bar_E8_yap_distance = ggplot(mtable_E8_yap_distance, aes(Distance_to_CD8_T_cells_µm, value, fill = variable)) +
  geom_bar(position = "fill", stat = "identity")+
  ggtitle("Relationship between the probability of seeing a YAP1-1 high epithelial nuclei and the distance of the nuclei to the CD8+ T cells")
bar_E8_yap_distance

per_E8_yaphigh_distance = data.frame (YAP1_high_fraction = c(ny820_yaphigh / ny820, ny830_yaphigh / ny830, ny840_yaphigh / ny840, ny8d_yaphigh / ny8d),Distance_to_CD8_T_cells_µm = c("0 - 100","100 - 150","150 - 200","200 - ∞"))
bar_per_E8_yaphigh_distance = ggplot(per_E8_yaphigh_distance, aes(x = Distance_to_CD8_T_cells_µm, y = YAP1_high_fraction)) +
  geom_bar(aes(fill = Distance_to_CD8_T_cells_µm), position = "dodge", stat = "identity")+
  scale_fill_manual(values=c("red", "purple", "blue", "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
bar_per_E8_yaphigh_distance

E8_distance = data.frame (Number_of_nuclei = c(ny820, ny830, ny840, ny8d),Distance_to_CD8_T_cells_µm = c("0 - 100","100 - 150","150 - 200","200 - ∞"))
bar_E8_distance = ggplot(E8_distance, aes(x = Distance_to_CD8_T_cells_µm, y = Number_of_nuclei)) +
  geom_bar(aes(fill = Distance_to_CD8_T_cells_µm), position = "dodge", stat = "identity")+
  scale_fill_manual(values=c("red", "purple", "blue", "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
bar_E8_distance

E8_yaphigh_distance = data.frame (Number_of_nuclei = c(ny820_yaphigh, ny830_yaphigh, ny840_yaphigh, ny8d_yaphigh),Distance_to_CD8_T_cells_µm = c("0 - 100","100 - 150","150 - 200","200 - ∞"))
bar_E8_yaphigh_distance = ggplot(E8_yaphigh_distance, aes(x = Distance_to_CD8_T_cells_µm, y = Number_of_nuclei)) +
  geom_bar(aes(fill = Distance_to_CD8_T_cells_µm), position = "dodge", stat = "identity")+
  scale_fill_manual(values=c("red", "purple", "blue", "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle("Number of yap-high epithelial nuclei with different distances to the CD8+ T cells")
bar_E8_yaphigh_distance
################################################################################
#Hypergeometric statistics for parp-high nuclei CD68Mp
CD68Mp_simple_parp = data.frame(PARP = CD68Mp234$PARP, proximity_CD68Mp = CD68Mp234$CD68Mp_proximity_status234)
npm20 = nrow(CD68Mp_simple_parp[CD68Mp_simple_parp$proximity_CD68Mp == "0 - 100",])
npm30 = nrow(CD68Mp_simple_parp[CD68Mp_simple_parp$proximity_CD68Mp == "100 - 150",])
npm40 = nrow(CD68Mp_simple_parp[CD68Mp_simple_parp$proximity_CD68Mp == "150 - 200",])
npmd = nrow(CD68Mp_simple_parp[CD68Mp_simple_parp$proximity_CD68Mp == "200 - ∞",])
npm20_parphigh = nrow(CD68Mp_simple_parp[CD68Mp_simple_parp$proximity_CD68Mp == "0 - 100" & CD68Mp_simple_parp$PARP == 1,])
npm30_parphigh = nrow(CD68Mp_simple_parp[CD68Mp_simple_parp$proximity_CD68Mp == "100 - 150" & CD68Mp_simple_parp$PARP == 1,])
npm40_parphigh = nrow(CD68Mp_simple_parp[CD68Mp_simple_parp$proximity_CD68Mp == "150 - 200" & CD68Mp_simple_parp$PARP == 1,])
npmd_parphigh = nrow(CD68Mp_simple_parp[CD68Mp_simple_parp$proximity_CD68Mp == "200 - ∞" & CD68Mp_simple_parp$PARP == 1,])
npm_parphigh = npm20_parphigh + npm30_parphigh + npm40_parphigh + npmd_parphigh
npm_parplow = npm20 + npm30 + npm40 + npmd - npm_parphigh  

pm_parphigh_20 = dhyper(npm20_parphigh, npm_parphigh, npm_parplow, npm20)
pm_parphigh_30 = dhyper(npm30_parphigh, npm_parphigh, npm_parplow, npm30)
pm_parphigh_40 = dhyper(npm40_parphigh, npm_parphigh, npm_parplow, npm40)
pm_parphigh_distal = dhyper(npmd_parphigh, npm_parphigh, npm_parplow, npmd)

dhyper_EM_PARP = data.frame(Proximity = c("0 - 100", "100 - 150", "150 - 200", "200 - ∞"),
                            number_of_parphigh_drawn = c(npm20_parphigh,npm30_parphigh,npm40_parphigh,npmd_parphigh),
                            number_of_parphigh = c(npm_parphigh,npm_parphigh,npm_parphigh,npm_parphigh),
                            number_of_parplow = c(npm_parplow,npm_parplow,npm_parplow,npm_parplow),
                            number_of_drawn = c(npm20,npm30,npm40,npmd),
                            probability = c(pm_parphigh_20, pm_parphigh_30, pm_parphigh_40, pm_parphigh_distal))

emp20 = c(npm20_parphigh,(npm20 - npm20_parphigh))
emp30 = c(npm30_parphigh,(npm30 - npm30_parphigh))
emp40 = c(npm40_parphigh,(npm40 - npm40_parphigh))
empd = c(npmd_parphigh,(npmd - npmd_parphigh))
emp = data.frame(emp20, emp30, emp40, empd)
empc = chisq.test(emp)

#Visualize the idea that the probability of seeing parp-high epithelial nuclei increases as the distance to the CD68+ macrophage decreases.
table_EM_parp_distance = data.frame(PARP_low = c(npm20 - npm20_parphigh, npm30 - npm30_parphigh, npm40 - npm40_parphigh, npmd - npmd_parphigh), PARP_high = c(npm20_parphigh, npm30_parphigh, npm40_parphigh, npmd_parphigh), row.names = c("0 - 100","100 - 150","150 - 200","200 - ∞"))

table_EM_parp_distance$Distance_to_macrophage_µm = row.names(table_EM_parp_distance)
mtable_EM_parp_distance = melt(table_EM_parp_distance, id.vars = "Distance_to_macrophage_µm")
bar_EM_parp_distance = ggplot(mtable_EM_parp_distance, aes(Distance_to_macrophage_µm, value, fill = variable)) +
  geom_bar(position = "fill", stat = "identity")+
  ggtitle("Relationship between the probability of seeing a PARP-1 high epithelial nuclei and the distance of the nuclei to the CD68+ macrophage")
bar_EM_parp_distance

per_EM_parphigh_distance = data.frame (PARP_high_fraction = c(npm20_parphigh / npm20, npm30_parphigh / npm30, npm40_parphigh / npm40, npmd_parphigh / npmd),Distance_to_macrophage_µm = c("0 - 100","100 - 150","150 - 200","200 - ∞"))
bar_per_EM_parphigh_distance = ggplot(per_EM_parphigh_distance, aes(x = Distance_to_macrophage_µm, y = PARP_high_fraction)) +
  geom_bar(aes(fill = Distance_to_macrophage_µm), position = "dodge", stat = "identity")+
  scale_fill_manual(values=c("red", "purple", "blue", "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
bar_per_EM_parphigh_distance

EM_distance = data.frame (Number_of_nuclei = c(npm20, npm30, npm40, npmd),Distance_to_macrophage_µm = c("0 - 100","100 - 150","150 - 200","200 - ∞"))
bar_EM_distance = ggplot(EM_distance, aes(x = Distance_to_macrophage_µm, y = Number_of_nuclei)) +
  geom_bar(aes(fill = Distance_to_macrophage_µm), position = "dodge", stat = "identity")+
  scale_fill_manual(values=c("red", "purple", "blue", "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle("Number of epithelial nuclei with different distances to the CD68+ macrophage")
bar_EM_distance

EM_parphigh_distance = data.frame (Number_of_nuclei = c(npm20_parphigh, npm30_parphigh, npm40_parphigh, npmd_parphigh),Distance_to_macrophage_µm = c("0 - 100","100 - 150","150 - 200","200 - ∞"))
bar_EM_parphigh_distance = ggplot(EM_parphigh_distance, aes(x = Distance_to_macrophage_µm, y = Number_of_nuclei)) +
  geom_bar(aes(fill = Distance_to_macrophage_µm), position = "dodge", stat = "identity")+
  scale_fill_manual(values=c("red", "purple", "blue", "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle("Number of parp-high epithelial nuclei with different distances to the CD68+ macrophage")
bar_EM_parphigh_distance


################################################################################
#Hypergeometric statistics for fra-high nuclei CD68Mp
CD68Mp_simple_fra = data.frame(FRA1 = CD68Mp234$FRA1, proximity_CD68Mp = CD68Mp234$CD68Mp_proximity_status234)
nfm20 = nrow(CD68Mp_simple_fra[CD68Mp_simple_fra$proximity_CD68Mp == "0 - 100",])
nfm30 = nrow(CD68Mp_simple_fra[CD68Mp_simple_fra$proximity_CD68Mp == "100 - 150",])
nfm40 = nrow(CD68Mp_simple_fra[CD68Mp_simple_fra$proximity_CD68Mp == "150 - 200",])
nfmd = nrow(CD68Mp_simple_fra[CD68Mp_simple_fra$proximity_CD68Mp == "200 - ∞",])
nfm20_frahigh = nrow(CD68Mp_simple_fra[CD68Mp_simple_fra$proximity_CD68Mp == "0 - 100" & CD68Mp_simple_fra$FRA1 == 1,])
nfm30_frahigh = nrow(CD68Mp_simple_fra[CD68Mp_simple_fra$proximity_CD68Mp == "100 - 150" & CD68Mp_simple_fra$FRA1 == 1,])
nfm40_frahigh = nrow(CD68Mp_simple_fra[CD68Mp_simple_fra$proximity_CD68Mp == "150 - 200" & CD68Mp_simple_fra$FRA1 == 1,])
nfmd_frahigh = nrow(CD68Mp_simple_fra[CD68Mp_simple_fra$proximity_CD68Mp == "200 - ∞" & CD68Mp_simple_fra$FRA1 == 1,])
nfm_frahigh = nfm20_frahigh + nfm30_frahigh + nfm40_frahigh + nfmd_frahigh
nfm_fralow = nfm20 + nfm30 + nfm40 + nfmd - nfm_frahigh  

pm_frahigh_20 = dhyper(nfm20_frahigh, nfm_frahigh, nfm_fralow, nfm20)
pm_frahigh_30 = dhyper(nfm30_frahigh, nfm_frahigh, nfm_fralow, nfm30)
pm_frahigh_40 = dhyper(nfm40_frahigh, nfm_frahigh, nfm_fralow, nfm40)
pm_frahigh_distal = dhyper(nfmd_frahigh, nfm_frahigh, nfm_fralow, nfmd)

dhyper_EM_FRA = data.frame(Proximity = c("0 - 100", "100 - 150", "150 - 200", "200 - ∞"),
                           number_of_frahigh_drawn = c(nfm20_frahigh,nfm30_frahigh,nfm40_frahigh,nfmd_frahigh),
                           number_of_frahigh = c(nfm_frahigh,nfm_frahigh,nfm_frahigh,nfm_frahigh),
                           number_of_fralow = c(nfm_fralow,nfm_fralow,nfm_fralow,nfm_fralow),
                           number_of_drawn = c(nfm20,nfm30,nfm40,nfmd),
                           probability = c(pm_frahigh_20, pm_frahigh_30, pm_frahigh_40, pm_frahigh_distal))

emf20 = c(nfm20_frahigh,(nfm20 - nfm20_frahigh))
emf30 = c(nfm30_frahigh,(nfm30 - nfm30_frahigh))
emf40 = c(nfm40_frahigh,(nfm40 - nfm40_frahigh))
emfd = c(nfmd_frahigh,(nfmd - nfmd_frahigh))
emf = data.frame(emf20, emf30, emf40, emfd)
emfc = chisq.test(emf)

#Visualize the idea that the probability of seeing fra-high epithelial nuclei increases as the distance to the CD68+ macrophage decreases.
table_EM_fra1_distance = data.frame(Fra_low = c(nfm20 - nfm20_frahigh, nfm30 - nfm30_frahigh, nfm40 - nfm40_frahigh, nfmd - nfmd_frahigh), Fra_high = c(nfm20_frahigh, nfm30_frahigh, nfm40_frahigh, nfmd_frahigh), row.names = c("0 - 100","100 - 150","150 - 200","200 - ∞"))

table_EM_fra1_distance$Distance_to_macrophage_µm = row.names(table_EM_fra1_distance)
mtable_EM_fra1_distance = melt(table_EM_fra1_distance, id.vars = "Distance_to_macrophage_µm")
bar_EM_fra1_distance = ggplot(mtable_EM_fra1_distance, aes(Distance_to_macrophage_µm, value, fill = variable)) +
  geom_bar(position = "fill", stat = "identity")+
  ggtitle("Relationship between the probability of seeing a Fra-1 high epithelial nuclei and the distance of the nuclei to the CD68+ macrophage")
bar_EM_fra1_distance

per_EM_frahigh_distance = data.frame (Fra_high_fraction = c(nfm20_frahigh / nfm20, nfm30_frahigh / nfm30, nfm40_frahigh / nfm40, nfmd_frahigh / nfmd),Distance_to_macrophage_µm = c("0 - 100","100 - 150","150 - 200","200 - ∞"))
bar_per_EM_frahigh_distance = ggplot(per_EM_frahigh_distance, aes(x = Distance_to_macrophage_µm, y = Fra_high_fraction)) +
  geom_bar(aes(fill = Distance_to_macrophage_µm), position = "dodge", stat = "identity")+
  scale_fill_manual(values=c("red", "purple", "blue", "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
bar_per_EM_frahigh_distance

EM_distance = data.frame (Number_of_nuclei = c(nfm20, nfm30, nfm40, nfmd),Distance_to_macrophage_µm = c("0 - 100","100 - 150","150 - 200","200 - ∞"))
bar_EM_distance = ggplot(EM_distance, aes(x = Distance_to_macrophage_µm, y = Number_of_nuclei)) +
  geom_bar(aes(fill = Distance_to_macrophage_µm), position = "dodge", stat = "identity")+
  scale_fill_manual(values=c("red", "purple", "blue", "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle("Number of epithelial nuclei with different distances to the CD68+ macrophage")
bar_EM_distance

EM_frahigh_distance = data.frame (Number_of_nuclei = c(nfm20_frahigh, nfm30_frahigh, nfm40_frahigh, nfmd_frahigh),Distance_to_macrophage_µm = c("0 - 100","100 - 150","150 - 200","200 - ∞"))
bar_EM_frahigh_distance = ggplot(EM_frahigh_distance, aes(x = Distance_to_macrophage_µm, y = Number_of_nuclei)) +
  geom_bar(aes(fill = Distance_to_macrophage_µm), position = "dodge", stat = "identity")+
  scale_fill_manual(values=c("red", "purple", "blue", "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle("Number of fra-high epithelial nuclei with different distances to the CD68+ macrophage")
bar_EM_frahigh_distance

################################################################################
#Hypergeometric statistics for yap-high nuclei CD68Mp
CD68Mp_simple_yap = data.frame(YAP1 = CD68Mp234$YAP1, proximity_CD68Mp = CD68Mp234$CD68Mp_proximity_status234)
nym20 = nrow(CD68Mp_simple_yap[CD68Mp_simple_yap$proximity_CD68Mp == "0 - 100",])
nym30 = nrow(CD68Mp_simple_yap[CD68Mp_simple_yap$proximity_CD68Mp == "100 - 150",])
nym40 = nrow(CD68Mp_simple_yap[CD68Mp_simple_yap$proximity_CD68Mp == "150 - 200",])
nymd = nrow(CD68Mp_simple_yap[CD68Mp_simple_yap$proximity_CD68Mp == "200 - ∞",])
nym20_yaphigh = nrow(CD68Mp_simple_yap[CD68Mp_simple_yap$proximity_CD68Mp == "0 - 100" & CD68Mp_simple_yap$YAP1 == 1,])
nym30_yaphigh = nrow(CD68Mp_simple_yap[CD68Mp_simple_yap$proximity_CD68Mp == "100 - 150" & CD68Mp_simple_yap$YAP1 == 1,])
nym40_yaphigh = nrow(CD68Mp_simple_yap[CD68Mp_simple_yap$proximity_CD68Mp == "150 - 200" & CD68Mp_simple_yap$YAP1 == 1,])
nymd_yaphigh = nrow(CD68Mp_simple_yap[CD68Mp_simple_yap$proximity_CD68Mp == "200 - ∞" & CD68Mp_simple_yap$YAP1 == 1,])
nym_yaphigh = nym20_yaphigh + nym30_yaphigh + nym40_yaphigh + nymd_yaphigh
nym_yaplow = nym20 + nym30 + nym40 + nymd - nym_yaphigh  

pm_yaphigh_20 = dhyper(nym20_yaphigh, nym_yaphigh, nym_yaplow, nym20)
pm_yaphigh_30 = dhyper(nym30_yaphigh, nym_yaphigh, nym_yaplow, nym30)
pm_yaphigh_40 = dhyper(nym40_yaphigh, nym_yaphigh, nym_yaplow, nym40)
pm_yaphigh_distal = dhyper(nymd_yaphigh, nym_yaphigh, nym_yaplow, nymd)

dhyper_EM_YAP1 = data.frame(Proximity = c("0 - 100", "100 - 150", "150 - 200", "200 - ∞"),
                            number_of_yaphigh_drawn = c(nym20_yaphigh,nym30_yaphigh,nym40_yaphigh,nymd_yaphigh),
                            number_of_yaphigh = c(nym_yaphigh,nym_yaphigh,nym_yaphigh,nym_yaphigh),
                            number_of_yaplow = c(nym_yaplow,nym_yaplow,nym_yaplow,nym_yaplow),
                            number_of_drawn = c(nym20,nym30,nym40,nymd),
                            probability = c(pm_yaphigh_20, pm_yaphigh_30, pm_yaphigh_40, pm_yaphigh_distal))

emy20 = c(nym20_yaphigh,(nym20 - nym20_yaphigh))
emy30 = c(nym30_yaphigh,(nym30 - nym30_yaphigh))
emy40 = c(nym40_yaphigh,(nym40 - nym40_yaphigh))
emyd = c(nymd_yaphigh,(nymd - nymd_yaphigh))
emy = data.frame(emy20, emy30, emy40, emyd)
emyc = chisq.test(emy)

#Visualize the idea that the probability of seeing yap-high epithelial nuclei increases as the distance to the CD68+ macrophage decreases.
table_EM_yap_distance = data.frame(YAP1_low = c(nym20 - nym20_yaphigh, nym30 - nym30_yaphigh, nym40 - nym40_yaphigh, nymd - nymd_yaphigh), YAP1_high = c(nym20_yaphigh, nym30_yaphigh, nym40_yaphigh, nymd_yaphigh), row.names = c("0 - 100","100 - 150","150 - 200","200 - ∞"))

table_EM_yap_distance$Distance_to_macrophage_µm = row.names(table_EM_yap_distance)
mtable_EM_yap_distance = melt(table_EM_yap_distance, id.vars = "Distance_to_macrophage_µm")
bar_EM_yap_distance = ggplot(mtable_EM_yap_distance, aes(Distance_to_macrophage_µm, value, fill = variable)) +
  geom_bar(position = "fill", stat = "identity")+
  ggtitle("Relationship between the probability of seeing a YAP1-1 high epithelial nuclei and the distance of the nuclei to the CD68+ macrophage")
bar_EM_yap_distance

per_EM_yaphigh_distance = data.frame (YAP1_high_fraction = c(nym20_yaphigh / nym20, nym30_yaphigh / nym30, nym40_yaphigh / nym40, nymd_yaphigh / nymd),Distance_to_macrophage_µm = c("0 - 100","100 - 150","150 - 200","200 - ∞"))
bar_per_EM_yaphigh_distance = ggplot(per_EM_yaphigh_distance, aes(x = Distance_to_macrophage_µm, y = YAP1_high_fraction)) +
  geom_bar(aes(fill = Distance_to_macrophage_µm), position = "dodge", stat = "identity")+
  scale_fill_manual(values=c("red", "purple", "blue", "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
bar_per_EM_yaphigh_distance

EM_distance = data.frame (Number_of_nuclei = c(nym20, nym30, nym40, nymd),Distance_to_macrophage_µm = c("0 - 100","100 - 150","150 - 200","200 - ∞"))
bar_EM_distance = ggplot(EM_distance, aes(x = Distance_to_macrophage_µm, y = Number_of_nuclei)) +
  geom_bar(aes(fill = Distance_to_macrophage_µm), position = "dodge", stat = "identity")+
  scale_fill_manual(values=c("red", "purple", "blue", "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
bar_EM_distance

EM_yaphigh_distance = data.frame (Number_of_nuclei = c(nym20_yaphigh, nym30_yaphigh, nym40_yaphigh, nymd_yaphigh),Distance_to_macrophage_µm = c("0 - 100","100 - 150","150 - 200","200 - ∞"))
bar_EM_yaphigh_distance = ggplot(EM_yaphigh_distance, aes(x = Distance_to_macrophage_µm, y = Number_of_nuclei)) +
  geom_bar(aes(fill = Distance_to_macrophage_µm), position = "dodge", stat = "identity")+
  scale_fill_manual(values=c("red", "purple", "blue", "black")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ggtitle("Number of yap-high epithelial nuclei with different distances to the CD68+ macrophage")
bar_EM_yaphigh_distance

################################################################################
################################################################################
png(filename = "/~/Analysis/Samplex/Samplex_Spatial_ESI.png", width = 2000, height = 1500, units = "px")
Spatial_ESI
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_Spatial_ES_PARP.png", width = 1000, height = 750, units = "px")
Spatial_ES_PARP
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_Spatial_ES_FRA1.png", width = 1000, height = 750, units = "px")
Spatial_ES_FRA1
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_Spatial_ES_YAP1.png", width = 1000, height = 750, units = "px")
Spatial_ES_YAP1
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_Spatial_EIS.png", width = 2000, height = 1500, units = "px")
Spatial_EIS
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_Spatial_EI_PARP.png", width = 2000, height = 1500, units = "px")
Spatial_EI_PARP
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_Spatial_EI_FRA1.png", width = 2000, height = 1500, units = "px")
Spatial_EI_FRA1
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_Spatial_EI_YAP1.png", width = 2000, height = 1500, units = "px")
Spatial_EI_YAP1
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_Spatial_ET.png", width = 2000, height = 1500, units = "px")
Spatial_ET
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_Spatial_ET_PARP.png", width = 2000, height = 1500, units = "px")
Spatial_ET_PARP
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_Spatial_ET_FRA1.png", width = 2000, height = 1500, units = "px")
Spatial_ET_FRA1
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_Spatial_ET_YAP1.png", width = 2000, height = 1500, units = "px")
Spatial_ET_YAP1
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_Spatial_E8.png", width = 2000, height = 1500, units = "px")
Spatial_E8
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_Spatial_E8_PARP.png", width = 1000, height = 750, units = "px")
Spatial_E8_PARP
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_Spatial_E8_FRA1.png", width = 1000, height = 750, units = "px")
Spatial_E8_FRA1
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_Spatial_E8_YAP1.png", width = 1000, height = 750, units = "px")
Spatial_E8_YAP1
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_Spatial_EM.png", width = 2000, height = 1500, units = "px")
Spatial_EM
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_Spatial_EM_PARP.png", width = 1000, height = 750, units = "px")
Spatial_EM_PARP
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_Spatial_EM_FRA1.png", width = 1000, height = 750, units = "px")
Spatial_EM_FRA1
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_Spatial_EM_YAP1.png", width = 1000, height = 750, units = "px")
Spatial_EM_YAP1
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_bar_ES_distance.png", width = 700, height = 350, units = "px")
bar_ES_distance
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_bar_ES_parp_distance.png", width = 1000, height = 1000, units = "px")
bar_ES_parp_distance
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_bar_ES_parphigh_distance.png", width = 1000, height = 1000, units = "px")
bar_ES_parphigh_distance
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_bar_per_ES_parphigh_distance.png", width = 700, height = 350, units = "px")
bar_per_ES_parphigh_distance
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_bar_ES_fra1_distance.png", width = 1000, height = 1000, units = "px")
bar_ES_fra1_distance
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_bar_ES_frahigh_distance.png", width = 1000, height = 1000, units = "px")
bar_ES_frahigh_distance
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_bar_per_ES_frahigh_distance.png", width = 700, height = 350, units = "px")
bar_per_ES_frahigh_distance
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_bar_ES_yap_distance.png", width = 1000, height = 1000, units = "px")
bar_ES_yap_distance
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_bar_ES_yaphigh_distance.png", width = 1000, height = 1000, units = "px")
bar_ES_yaphigh_distance
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_bar_per_ES_yaphigh_distance.png", width = 700, height = 350, units = "px")
bar_per_ES_yaphigh_distance
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_bar_EI_distance.png", width = 1000, height = 1000, units = "px")
bar_EI_distance
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_bar_EI_parp_distance.png", width = 1000, height = 1000, units = "px")
bar_EI_parp_distance
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_bar_EI_parphigh_distance.png", width = 1000, height = 1000, units = "px")
bar_EI_parphigh_distance
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_bar_per_EI_parphigh_distance.png", width = 1000, height = 1000, units = "px")
bar_per_EI_parphigh_distance
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_bar_EI_fra1_distance.png", width = 1000, height = 1000, units = "px")
bar_EI_fra1_distance
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_bar_EI_frahigh_distance.png", width = 1000, height = 1000, units = "px")
bar_EI_frahigh_distance
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_bar_per_EI_frahigh_distance.png", width = 1000, height = 1000, units = "px")
bar_per_EI_frahigh_distance
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_bar_EI_yap_distance.png", width = 1000, height = 1000, units = "px")
bar_EI_yap_distance
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_bar_EI_yaphigh_distance.png", width = 1000, height = 1000, units = "px")
bar_EI_yaphigh_distance
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_bar_per_EI_yaphigh_distance.png", width = 1000, height = 1000, units = "px")
bar_per_EI_yaphigh_distance
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_bar_ET_distance.png", width = 1000, height = 1000, units = "px")
bar_ET_distance
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_bar_ET_parp_distance.png", width = 1000, height = 1000, units = "px")
bar_ET_parp_distance
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_bar_ET_parphigh_distance.png", width = 1000, height = 1000, units = "px")
bar_ET_parphigh_distance
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_bar_per_ET_parphigh_distance.png", width = 1000, height = 1000, units = "px")
bar_per_ET_parphigh_distance
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_bar_ET_fra1_distance.png", width = 1000, height = 1000, units = "px")
bar_ET_fra1_distance
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_bar_ET_frahigh_distance.png", width = 1000, height = 1000, units = "px")
bar_ET_frahigh_distance
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_bar_per_ET_frahigh_distance.png", width = 1000, height = 1000, units = "px")
bar_per_ET_frahigh_distance
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_bar_ET_yap_distance.png", width = 1000, height = 1000, units = "px")
bar_ET_yap_distance
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_bar_ET_yaphigh_distance.png", width = 1000, height = 1000, units = "px")
bar_ET_yaphigh_distance
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_bar_per_ET_yaphigh_distance.png", width = 1000, height = 1000, units = "px")
bar_per_ET_yaphigh_distance
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_bar_E8_distance.png", width = 700, height = 350, units = "px")
bar_E8_distance
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_bar_E8_parp_distance.png", width = 1000, height = 1000, units = "px")
bar_E8_parp_distance
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_bar_E8_parphigh_distance.png", width = 1000, height = 1000, units = "px")
bar_E8_parphigh_distance
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_bar_per_E8_parphigh_distance.png", width = 700, height = 350, units = "px")
bar_per_E8_parphigh_distance
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_bar_E8_fra1_distance.png", width = 1000, height = 1000, units = "px")
bar_E8_fra1_distance
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_bar_E8_frahigh_distance.png", width = 1000, height = 1000, units = "px")
bar_E8_frahigh_distance
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_bar_per_E8_frahigh_distance.png", width = 700, height = 350, units = "px")
bar_per_E8_frahigh_distance
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_bar_E8_yap_distance.png", width = 1000, height = 1000, units = "px")
bar_E8_yap_distance
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_bar_E8_yaphigh_distance.png", width = 1000, height = 1000, units = "px")
bar_E8_yaphigh_distance
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_bar_per_E8_yaphigh_distance.png", width = 700, height = 350, units = "px")
bar_per_E8_yaphigh_distance
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_bar_EM_distance.png", width = 700, height = 350, units = "px")
bar_EM_distance
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_bar_EM_parp_distance.png", width = 1000, height = 1000, units = "px")
bar_EM_parp_distance
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_bar_EM_parphigh_distance.png", width = 1000, height = 1000, units = "px")
bar_EM_parphigh_distance
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_bar_per_EM_parphigh_distance.png", width = 700, height = 350, units = "px")
bar_per_EM_parphigh_distance
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_bar_EM_fra1_distance.png", width = 1000, height = 1000, units = "px")
bar_EM_fra1_distance
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_bar_EM_frahigh_distance.png", width = 1000, height = 1000, units = "px")
bar_EM_frahigh_distance
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_bar_per_EM_frahigh_distance.png", width = 700, height = 350, units = "px")
bar_per_EM_frahigh_distance
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_bar_EM_yap_distance.png", width = 1000, height = 1000, units = "px")
bar_EM_yap_distance
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_bar_EM_yaphigh_distance.png", width = 1000, height = 1000, units = "px")
bar_EM_yaphigh_distance
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_bar_per_EM_yaphigh_distance.png", width = 700, height = 350, units = "px")
bar_per_EM_yaphigh_distance
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_Vio_ES_PARP.png", width = 1000, height = 1000, units = "px")
Vio_ES_PARP
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_Vio_ES_FRA.png", width = 1000, height = 1000, units = "px")
Vio_ES_FRA
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_Vio_ES_YAP.png", width = 1000, height = 1000, units = "px")
Vio_ES_YAP
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_Vio_EI_PARP.png", width = 1000, height = 1000, units = "px")
Vio_EI_PARP
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_Vio_EI_FRA.png", width = 1000, height = 1000, units = "px")
Vio_EI_FRA
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_Vio_EI_YAP.png", width = 1000, height = 1000, units = "px")
Vio_EI_YAP
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_Vio_ET_PARP.png", width = 1000, height = 1000, units = "px")
Vio_ET_PARP
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_Vio_ET_FRA.png", width = 1000, height = 1000, units = "px")
Vio_ET_FRA
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_Vio_ET_YAP.png", width = 1000, height = 1000, units = "px")
Vio_ET_YAP
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_Vio_E8_PARP.png", width = 1000, height = 1000, units = "px")
Vio_E8_PARP
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_Vio_E8_FRA.png", width = 1000, height = 1000, units = "px")
Vio_E8_FRA
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_Vio_E8_YAP.png", width = 1000, height = 1000, units = "px")
Vio_E8_YAP
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_Vio_EM_PARP.png", width = 1000, height = 1000, units = "px")
Vio_EM_PARP
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_Vio_EM_FRA.png", width = 1000, height = 1000, units = "px")
Vio_EM_FRA
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_Vio_EM_YAP.png", width = 1000, height = 1000, units = "px")
Vio_EM_YAP
dev.off()


library(spatstat);library(splancs)
epifra = epi234%>%
  filter(epi234$FRA1 == 1)
epifralow = epi234%>%
  filter(epi234$FRA1 == 0)
epiyap = epi234%>%
  filter(epi234$YAP1 == 1)
epiyaplow = epi234%>%
  filter(epi234$YAP1 == 0)
epiparp = epi234%>%
  filter(epi234$PARP == 1)
epiparplow = epi234%>%
  filter(epi234$PARP == 0)

sf234 = rbind(epifra, str_trim)
sf234$cell_type = "Stromal"
sf234$cell_type[sf234$Epi == 1] = "Epithelial"
cell_type1f = factor(sf234$cell_type)
csf234 = ppp(x = sf234$X.Nucleus1..X.CoordinateInWell,y = sf234$X.Nucleus1..Y.CoordinateInWell,marks = cell_type1f, window = owin(c(-2000,2000), c(-2000,2000)))
Kcsf = envelope(csf234,Kcross, r = c(0:60))
plot(Kcsf)
sfl234 = rbind(epifralow, str_trim)
sfl234$cell_type = "Stromal"
sfl234$cell_type[sfl234$Epi == 1] = "Epithelial"
cell_type1fl = factor(sfl234$cell_type)
csfl234 = ppp(x = sfl234$X.Nucleus1..X.CoordinateInWell,y = sfl234$X.Nucleus1..Y.CoordinateInWell,marks = cell_type1fl, window = owin(c(-2000,2000), c(-2000,2000)))
Kcsfl = envelope(csfl234,Kcross, r = c(0:60))
plot(Kcsfl)
#hsf = hyperframe(ppp = csf234, group = cell_type1f)
#hsf
#sf_test = studpermu.test(hsf, ppp ~ group, nperm = 10)

Gcsf = Gcross(csf234, i = "Epithelial", j = "Stromal",r = c(0:60), correction = "km")
plot(Gcsf)

sy234 = rbind(epiyap, str_trim)
sy234$cell_type = "Stromal"
sy234$cell_type[sy234$Epi == 1] = "Epithelial"
cell_type1y = factor(sy234$cell_type)
csy234 = ppp(x = sy234$X.Nucleus1..X.CoordinateInWell,y = sy234$X.Nucleus1..Y.CoordinateInWell,marks = cell_type1y, window = owin(c(-2000,2000), c(-2000,2000)))
Kcsy = envelope(csy234,Kcross, r = c(0:60))
plot(Kcsy)
syl234 = rbind(epiyaplow, str_trim)
syl234$cell_type = "Stromal"
syl234$cell_type[syl234$Epi == 1] = "Epithelial"
cell_type1yl = factor(syl234$cell_type)
csyl234 = ppp(x = syl234$X.Nucleus1..X.CoordinateInWell,y = syl234$X.Nucleus1..Y.CoordinateInWell,marks = cell_type1yl, window = owin(c(-2000,2000), c(-2000,2000)))
Kcsyl = envelope(csyl234,Kcross, r = c(0:60))
plot(Kcsyl)
#plot(envelope(csyl234,Kcross,r = c(0,60)))
Gcsy = Gcross(csy234, i = "Epithelial", j = "Stromal", r = c(0:60), correction = "km")
plot(Gcsy)

sp234 = rbind(epiparp, str_trim)
sp234$cell_type = "Stromal"
sp234$cell_type[sp234$Epi == 1] = "Epithelial"
cell_type1p = factor(sp234$cell_type)
csp234 = ppp(x = sp234$X.Nucleus1..X.CoordinateInWell,y = sp234$X.Nucleus1..Y.CoordinateInWell,marks = cell_type1p, window = owin(c(-2000,2000), c(-2000,2000)))
Kcsp = envelope(csp234,Kcross, r = c(0:60))
plot(Kcsp)
spl234 = rbind(epiparplow, str_trim)
spl234$cell_type = "Stromal"
spl234$cell_type[spl234$Epi == 1] = "Epithelial"
cell_type1pl = factor(spl234$cell_type)
cspl234 = ppp(x = spl234$X.Nucleus1..X.CoordinateInWell,y = spl234$X.Nucleus1..Y.CoordinateInWell,marks = cell_type1pl, window = owin(c(-2000,2000), c(-2000,2000)))
Kcspl = envelope(cspl234,Kcross, r = c(0:60))
plot(Kcspl)
Gcsp = Gcross(csp234, i = "Epithelial", j = "Stromal", r = c(0:60), correction = "km")
plot(Gcsp)

Tf234 = rbind(epifra, T_cells)
Tf234$cell_type = "T cell"
Tf234$cell_type[Tf234$Epi == 1] = "Epithelial"
cell_type2f = factor(Tf234$cell_type)
cTf234 = ppp(x = Tf234$X.Nucleus1..X.CoordinateInWell,y = Tf234$X.Nucleus1..Y.CoordinateInWell,marks = cell_type2f, window = owin(c(-2000,2000), c(-2000,2000)))
KcTf = envelope(cTf234,Kcross, r = c(0:200))
plot(KcTf)
GcTf = Gcross(cTf234, i = "Epithelial", j = "T cell", r = c(0:200), correction = "km")
plot(GcTf)

Ty234 = rbind(epiyap, T_cells)
Ty234$cell_type = "T cell"
Ty234$cell_type[Ty234$Epi == 1] = "Epithelial"
cell_type2y = factor(Ty234$cell_type)
cTy234 = ppp(x = Ty234$X.Nucleus1..X.CoordinateInWell,y = Ty234$X.Nucleus1..Y.CoordinateInWell,marks = cell_type2y, window = owin(c(-2000,2000), c(-2000,2000)))
KcTy = envelope(cTy234,Kcross, r = c(0:200))
plot(KcTy)
GcTy = Gcross(cTy234, i = "Epithelial", j = "T cell", r = c(0:200), correction = "km")
plot(GcTy)

Tp234 = rbind(epiparp, T_cells)
Tp234$cell_type = "T cell"
Tp234$cell_type[Tp234$Epi == 1] = "Epithelial"
cell_type2p = factor(Tp234$cell_type)
cTp234 = ppp(x = Tp234$X.Nucleus1..X.CoordinateInWell,y = Tp234$X.Nucleus1..Y.CoordinateInWell,marks = cell_type2p, window = owin(c(-2000,2000), c(-2000,2000)))
KcTp = envelope(cTp234,Kcross, r = c(0:200))
plot(KcTp)
GcTp = Gcross(cTp234, i = "Epithelial", j = "T cell", r = c(0:200), correction = "km")
plot(GcTp)

CD8_Tf234 = rbind(epifra, CD8_T_cells)
CD8_Tf234$cell_type = "CD8+ T cell"
CD8_Tf234$cell_type[CD8_Tf234$Epi == 1] = "Epithelial"
cell_type3f = factor(CD8_Tf234$cell_type)
cCD8_Tf234 = ppp(x = CD8_Tf234$X.Nucleus1..X.CoordinateInWell,y = CD8_Tf234$X.Nucleus1..Y.CoordinateInWell,marks = cell_type3f, window = owin(c(-2000,2000), c(-2000,2000)))
KcCD8_Tf = envelope(cCD8_Tf234,Kcross, r = c(0:200))
plot(KcCD8_Tf)
CD8_Tfl234 = rbind(epifralow, CD8_T_cells)
CD8_Tfl234$cell_type = "CD8+ T cell"
CD8_Tfl234$cell_type[CD8_Tfl234$Epi == 1] = "Epithelial"
cell_type3fl = factor(CD8_Tfl234$cell_type)
cCD8_Tfl234 = ppp(x = CD8_Tfl234$X.Nucleus1..X.CoordinateInWell,y = CD8_Tfl234$X.Nucleus1..Y.CoordinateInWell,marks = cell_type3fl, window = owin(c(-2000,2000), c(-2000,2000)))
KcCD8_Tfl = envelope(cCD8_Tfl234,Kcross, r = c(0:200))
plot(KcCD8_Tfl)
GcCD8_Tf = Gcross(cCD8_Tf234, i = "Epithelial", j = "CD8+ T cell", r = c(0:200), correction = "km")
plot(GcCD8_Tf)

CD8_Ty234 = rbind(epiyap, CD8_T_cells)
CD8_Ty234$cell_type = "CD8+ T cell"
CD8_Ty234$cell_type[CD8_Ty234$Epi == 1] = "Epithelial"
cell_type3y = factor(CD8_Ty234$cell_type)
cCD8_Ty234 = ppp(x = CD8_Ty234$X.Nucleus1..X.CoordinateInWell,y = CD8_Ty234$X.Nucleus1..Y.CoordinateInWell,marks = cell_type3y, window = owin(c(-2000,2000), c(-2000,2000)))
KcCD8_Ty = envelope(cCD8_Ty234,Kcross, r = c(0:200))
plot(KcCD8_Ty)
CD8_Tyl234 = rbind(epiyaplow, CD8_T_cells)
CD8_Tyl234$cell_Type = "CD8+ T cell"
CD8_Tyl234$cell_Type[CD8_Tyl234$Epi == 1] = "Epithelial"
cell_Type3yl = factor(CD8_Tyl234$cell_Type)
cCD8_Tyl234 = ppp(x = CD8_Tyl234$X.Nucleus1..X.CoordinateInWell,y = CD8_Tyl234$X.Nucleus1..Y.CoordinateInWell,marks = cell_Type3yl, window = owin(c(-2000,2000), c(-2000,2000)))
KcCD8_Tyl = envelope(cCD8_Tyl234,Kcross, r = c(0:200))
plot(KcCD8_Tyl)
GcCD8_Ty = Gcross(cCD8_Ty234, i = "Epithelial", j = "CD8+ T cell", r = c(0:200), correction = "km")
plot(GcCD8_Ty)

CD8_Tp234 = rbind(epiparp, CD8_T_cells)
CD8_Tp234$cell_type = "CD8+ T cell"
CD8_Tp234$cell_type[CD8_Tp234$Epi == 1] = "Epithelial"
cell_type3p = factor(CD8_Tp234$cell_type)
cCD8_Tp234 = ppp(x = CD8_Tp234$X.Nucleus1..X.CoordinateInWell,y = CD8_Tp234$X.Nucleus1..Y.CoordinateInWell,marks = cell_type3p, window = owin(c(-2000,2000), c(-2000,2000)))
KcCD8_Tp = envelope(cCD8_Tp234,Kcross, r = c(0:200))
plot(KcCD8_Tp)
CD8_Tpl234 = rbind(epiparplow, CD8_T_cells)
CD8_Tpl234$cell_type = "CD8+ T cell"
CD8_Tpl234$cell_type[CD8_Tpl234$Epi == 1] = "Epithelial"
cell_type3pl = factor(CD8_Tpl234$cell_type)
cCD8_Tpl234 = ppp(x = CD8_Tpl234$X.Nucleus1..X.CoordinateInWell,y = CD8_Tpl234$X.Nucleus1..Y.CoordinateInWell,marks = cell_type3pl, window = owin(c(-2000,2000), c(-2000,2000)))
KcCD8_Tpl = envelope(cCD8_Tpl234,Kcross, r = c(0:200))
plot(KcCD8_Tpl)
GcCD8_Tp = Gcross(cCD8_Tp234, i = "Epithelial", j = "CD8+ T cell", r = c(0:200), correction = "km")
plot(GcCD8_Tp)

CD68Mpf234 = rbind(epifra, Macrophage)
CD68Mpf234$cell_type = "CD68+ Macrophage"
CD68Mpf234$cell_type[CD68Mpf234$Epi == 1] = "Epithelial"
cell_type4f = factor(CD68Mpf234$cell_type)
cCD68Mpf234 = ppp(x = CD68Mpf234$X.Nucleus1..X.CoordinateInWell,y = CD68Mpf234$X.Nucleus1..Y.CoordinateInWell,marks = cell_type4f, window = owin(c(-2000,2000), c(-2000,2000)))
KcCD68Mpf = envelope(cCD68Mpf234,Kcross, r = c(0:200))
plot(KcCD68Mpf)
CD68Mpfl234 = rbind(epifralow, Macrophage)
CD68Mpfl234$cell_type = "CD68+ Macrophage"
CD68Mpfl234$cell_type[CD68Mpfl234$Epi == 1] = "Epithelial"
cell_type4fl = factor(CD68Mpfl234$cell_type)
cCD68Mpfl234 = ppp(x = CD68Mpfl234$X.Nucleus1..X.CoordinateInWell,y = CD68Mpfl234$X.Nucleus1..Y.CoordinateInWell,marks = cell_type4fl, window = owin(c(-2000,2000), c(-2000,2000)))
KcCD68Mpfl = envelope(cCD68Mpfl234,Kcross, r = c(0:200))
plot(KcCD68Mpfl)
GcCD68Mpf = Gcross(cCD68Mpf234, i = "Epithelial", j = "CD68+ Macrophage", r = c(0:200), correction = "km")
plot(GcCD68Mpf)

CD68Mpy234 = rbind(epiyap, Macrophage)
CD68Mpy234$cell_type = "CD68+ Macrophage"
CD68Mpy234$cell_type[CD68Mpy234$Epi == 1] = "Epithelial"
cell_type4y = factor(CD68Mpy234$cell_type)
cCD68Mpy234 = ppp(x = CD68Mpy234$X.Nucleus1..X.CoordinateInWell,y = CD68Mpy234$X.Nucleus1..Y.CoordinateInWell,marks = cell_type4y, window = owin(c(-2000,2000), c(-2000,2000)))
KcCD68Mpy = envelope(cCD68Mpy234,Kcross, r = c(0:200))
plot(KcCD68Mpy)
CD68Mpyl234 = rbind(epiyaplow, Macrophage)
CD68Mpyl234$cell_type = "CD68+ Macrophage"
CD68Mpyl234$cell_type[CD68Mpyl234$Epi == 1] = "Epithelial"
cell_type4yl = factor(CD68Mpyl234$cell_type)
cCD68Mpyl234 = ppp(x = CD68Mpyl234$X.Nucleus1..X.CoordinateInWell,y = CD68Mpyl234$X.Nucleus1..Y.CoordinateInWell,marks = cell_type4yl, window = owin(c(-2000,2000), c(-2000,2000)))
KcCD68Mpyl = envelope(cCD68Mpyl234,Kcross, r = c(0:200))
plot(KcCD68Mpyl)
GcCD68Mpy = Gcross(cCD68Mpy234, i = "Epithelial", j = "CD68+ Macrophage", r = c(0:200), correction = "km")
plot(GcCD68Mpy)

CD68Mpp234 = rbind(epiparp, Macrophage)
CD68Mpp234$cell_type = "CD68+ Macrophage"
CD68Mpp234$cell_type[CD68Mpp234$Epi == 1] = "Epithelial"
cell_type4p = factor(CD68Mpp234$cell_type)
cCD68Mpp234 = ppp(x = CD68Mpp234$X.Nucleus1..X.CoordinateInWell,y = CD68Mpp234$X.Nucleus1..Y.CoordinateInWell,marks = cell_type4p, window = owin(c(-2000,2000), c(-2000,2000)))
KcCD68Mpp = envelope(cCD68Mpp234,Kcross, r = c(0:200))
plot(KcCD68Mpp)
CD68Mppl234 = rbind(epiparplow, Macrophage)
CD68Mppl234$cell_type = "CD68+ Macrophage"
CD68Mppl234$cell_type[CD68Mppl234$Epi == 1] = "Epithelial"
cell_type4pl = factor(CD68Mppl234$cell_type)
cCD68Mppl234 = ppp(x = CD68Mppl234$X.Nucleus1..X.CoordinateInWell,y = CD68Mppl234$X.Nucleus1..Y.CoordinateInWell,marks = cell_type4pl, window = owin(c(-2000,2000), c(-2000,2000)))
KcCD68Mppl = envelope(cCD68Mppl234,Kcross, r = c(0:200))
plot(KcCD68Mppl)
GcCD68Mpp = Gcross(cCD68Mpp234, i = "Epithelial", j = "CD68+ Macrophage", r = c(0:200), correction = "km")
plot(GcCD68Mpp)

png(filename = "/~/Analysis/Samplex/Samplex_Kcsf.png", width = 700, height = 600, units = "px")
plot(Kcsf)
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_Kcsy.png", width = 700, height = 600, units = "px")
plot(Kcsy)
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_Kcsp.png", width = 700, height = 600, units = "px")
plot(Kcsp)
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_Kcsfl.png", width = 700, height = 600, units = "px")
plot(Kcsfl)
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_Kcsyl.png", width = 700, height = 600, units = "px")
plot(Kcsyl)
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_Kcspl.png", width = 700, height = 600, units = "px")
plot(Kcspl)
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_Gcsf.png", width = 700, height = 600, units = "px")
plot(Gcsf)
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_Gcsy.png", width = 700, height = 600, units = "px")
plot(Gcsy)
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_Gcsp.png", width = 700, height = 600, units = "px")
plot(Gcsp)
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_KcTf.png", width = 700, height = 600, units = "px")
plot(KcTf)
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_KcTy.png", width = 700, height = 600, units = "px")
plot(KcTy)
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_KcTp.png", width = 700, height = 600, units = "px")
plot(KcTp)
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_GcTf.png", width = 700, height = 600, units = "px")
plot(GcTf)
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_GcTy.png", width = 700, height = 600, units = "px")
plot(GcTy)
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_GcTp.png", width = 700, height = 600, units = "px")
plot(GcTp)
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_KcCD8_Tf.png", width = 700, height = 600, units = "px")
plot(KcCD8_Tf)
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_KcCD8_Ty.png", width = 700, height = 600, units = "px")
plot(KcCD8_Ty)
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_KcCD8_Tp.png", width = 700, height = 600, units = "px")
plot(KcCD8_Tp)
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_KcCD8_Tfl.png", width = 700, height = 600, units = "px")
plot(KcCD8_Tfl)
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_KcCD8_Tyl.png", width = 700, height = 600, units = "px")
plot(KcCD8_Tyl)
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_KcCD8_Tpl.png", width = 700, height = 600, units = "px")
plot(KcCD8_Tpl)
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_GcCD8_Tf.png", width = 700, height = 600, units = "px")
plot(GcCD8_Tf)
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_GcCD8_Ty.png", width = 700, height = 600, units = "px")
plot(GcCD8_Ty)
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_GcCD8_Tp.png", width = 700, height = 600, units = "px")
plot(GcCD8_Tp)
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_KcCD68Mpf.png", width = 700, height = 600, units = "px")
plot(KcCD68Mpf)
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_KcCD68Mpy.png", width = 700, height = 600, units = "px")
plot(KcCD68Mpy)
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_KcCD68Mpp.png", width = 700, height = 600, units = "px")
plot(KcCD68Mpp)
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_KcCD68Mpfl.png", width = 700, height = 600, units = "px")
plot(KcCD68Mpfl)
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_KcCD68Mpyl.png", width = 700, height = 600, units = "px")
plot(KcCD68Mpyl)
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_KcCD68Mppl.png", width = 700, height = 600, units = "px")
plot(KcCD68Mppl)
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_GcCD68Mpf.png", width = 700, height = 600, units = "px")
plot(GcCD68Mpf)
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_GcCD68Mpy.png", width = 700, height = 600, units = "px")
plot(GcCD68Mpy)
dev.off()

png(filename = "/~/Analysis/Samplex/Samplex_GcCD68Mpp.png", width = 700, height = 600, units = "px")
plot(GcCD68Mpp)
dev.off()

write.csv(dhyper_ES_FRA, "/~/Analysis/Samplex/Samplex_dhyper_ES_FRA.csv", row.names=FALSE)
write.csv(dhyper_ES_PARP, "/~/Analysis/Samplex/Samplex_dhyper_ES_PARP.csv", row.names=FALSE)
write.csv(dhyper_ES_YAP1, "/~/Analysis/Samplex/Samplex_dhyper_ES_YAP1.csv", row.names=FALSE)
write.csv(dhyper_EI_FRA, "/~/Analysis/Samplex/Samplex_dhyper_EI_FRA.csv", row.names=FALSE)
write.csv(dhyper_EI_PARP, "/~/Analysis/Samplex/Samplex_dhyper_EI_PARP.csv", row.names =FALSE)
write.csv(dhyper_EI_YAP1, "/~/Analysis/Samplex/Samplex_dhyper_EI_YAP1.csv", row.names=FALSE)
write.csv(dhyper_ET_FRA, "/~/Analysis/Samplex/Samplex_dhyper_ET_FRA.csv", row.names=FALSE)
write.csv(dhyper_ET_PARP, "/~/Analysis/Samplex/Samplex_dhyper_ET_PARP.csv", row.names=FALSE)
write.csv(dhyper_ET_YAP1, "/~/Analysis/Samplex/Samplex_dhyper_ET_YAP1.csv", row.names=FALSE)
write.csv(dhyper_E8_FRA, "/~/Analysis/Samplex/Samplex_dhyper_E8_FRA.csv", row.names=FALSE)
write.csv(dhyper_E8_PARP, "/~/Analysis/Samplex/Samplex_dhyper_E8_PARP.csv", row.names=FALSE)
write.csv(dhyper_E8_YAP1, "/~/Analysis/Samplex/Samplex_dhyper_E8_YAP1.csv", row.names=FALSE)
write.csv(dhyper_EM_FRA, "/~/Analysis/Samplex/Samplex_dhyper_EM_FRA.csv", row.names=FALSE)
write.csv(dhyper_EM_PARP, "/~/Analysis/Samplex/Samplex_dhyper_EM_PARP.csv", row.names=FALSE)
write.csv(dhyper_EM_YAP1, "/~/Analysis/Samplex/Samplex_dhyper_EM_YAP1.csv", row.names=FALSE)
write.csv(CELL_STATS_TABLE, "/~/Analysis/Samplex/Samplex_CELL_STATS_TABLE.csv", row.names=FALSE)


gf1 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH4[epi_trim$X.Nucleus1..MeanIntensity.CH4 >= 0 & epi_trim$X.Nucleus1..MeanIntensity.CH4 <= 10]))
gf2 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH4[epi_trim$X.Nucleus1..MeanIntensity.CH4 > 10 & epi_trim$X.Nucleus1..MeanIntensity.CH4 <= 20]))
gf3 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH4[epi_trim$X.Nucleus1..MeanIntensity.CH4 > 20 & epi_trim$X.Nucleus1..MeanIntensity.CH4 <= 30]))
gf4 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH4[epi_trim$X.Nucleus1..MeanIntensity.CH4 > 30 & epi_trim$X.Nucleus1..MeanIntensity.CH4 <= 40]))
gf5 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH4[epi_trim$X.Nucleus1..MeanIntensity.CH4 > 40 & epi_trim$X.Nucleus1..MeanIntensity.CH4 <= 50]))
gf6 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH4[epi_trim$X.Nucleus1..MeanIntensity.CH4 > 50 & epi_trim$X.Nucleus1..MeanIntensity.CH4 <= 60]))
gf7 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH4[epi_trim$X.Nucleus1..MeanIntensity.CH4 > 60 & epi_trim$X.Nucleus1..MeanIntensity.CH4 <= 70]))
gf8 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH4[epi_trim$X.Nucleus1..MeanIntensity.CH4 > 70 & epi_trim$X.Nucleus1..MeanIntensity.CH4 <= 80]))
gf9 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH4[epi_trim$X.Nucleus1..MeanIntensity.CH4 > 80 & epi_trim$X.Nucleus1..MeanIntensity.CH4 <= 90]))
gf10 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH4[epi_trim$X.Nucleus1..MeanIntensity.CH4 > 90 & epi_trim$X.Nucleus1..MeanIntensity.CH4 <= 100]))
gf11= nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH4[epi_trim$X.Nucleus1..MeanIntensity.CH4 > 100 & epi_trim$X.Nucleus1..MeanIntensity.CH4 <= 110]))
gf12 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH4[epi_trim$X.Nucleus1..MeanIntensity.CH4 > 110 & epi_trim$X.Nucleus1..MeanIntensity.CH4 <= 120]))
gf13 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH4[epi_trim$X.Nucleus1..MeanIntensity.CH4 > 120 & epi_trim$X.Nucleus1..MeanIntensity.CH4 <= 130]))
gf14 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH4[epi_trim$X.Nucleus1..MeanIntensity.CH4 > 130 & epi_trim$X.Nucleus1..MeanIntensity.CH4 <= 140]))
gf15 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH4[epi_trim$X.Nucleus1..MeanIntensity.CH4 > 140 & epi_trim$X.Nucleus1..MeanIntensity.CH4 <= 150]))
gf16 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH4[epi_trim$X.Nucleus1..MeanIntensity.CH4 > 150 & epi_trim$X.Nucleus1..MeanIntensity.CH4 <= 160]))
gf17 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH4[epi_trim$X.Nucleus1..MeanIntensity.CH4 > 160 & epi_trim$X.Nucleus1..MeanIntensity.CH4 <= 170]))
gf18 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH4[epi_trim$X.Nucleus1..MeanIntensity.CH4 > 170 & epi_trim$X.Nucleus1..MeanIntensity.CH4 <= 180]))
gf19 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH4[epi_trim$X.Nucleus1..MeanIntensity.CH4 > 180 & epi_trim$X.Nucleus1..MeanIntensity.CH4 <= 190]))
gf20 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH4[epi_trim$X.Nucleus1..MeanIntensity.CH4 > 190 & epi_trim$X.Nucleus1..MeanIntensity.CH4 <= 200]))
gf21 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH4[epi_trim$X.Nucleus1..MeanIntensity.CH4 > 200 & epi_trim$X.Nucleus1..MeanIntensity.CH4 <= 210]))
gf22 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH4[epi_trim$X.Nucleus1..MeanIntensity.CH4 > 210 & epi_trim$X.Nucleus1..MeanIntensity.CH4 <= 220]))
gf23 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH4[epi_trim$X.Nucleus1..MeanIntensity.CH4 > 220 & epi_trim$X.Nucleus1..MeanIntensity.CH4 <= 230]))
gf24 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH4[epi_trim$X.Nucleus1..MeanIntensity.CH4 > 230 & epi_trim$X.Nucleus1..MeanIntensity.CH4 <= 240]))
gf25 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH4[epi_trim$X.Nucleus1..MeanIntensity.CH4 > 240 & epi_trim$X.Nucleus1..MeanIntensity.CH4 <= 250]))
gf26 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH4[epi_trim$X.Nucleus1..MeanIntensity.CH4 > 250]))

FRA_DATA = data.frame(f1 = gf1,
                      f2 = gf2,
                      f3 = gf3,
                      f4 = gf4,
                      f5 = gf5,
                      f6 = gf6,
                      f7 = gf7,
                      f8 = gf8,
                      f9 = gf9,
                      f10 = gf10,
                      f11 = gf11,
                      f12 = gf12,
                      f13 = gf13,
                      f14 = gf14,
                      f15 = gf15,
                      f16 = gf16,
                      f17 = gf17,
                      f18 = gf18,
                      f19 = gf19,
                      f20 = gf20,
                      f21 = gf21,
                      f22 = gf22,
                      f23 = gf23,
                      f24 = gf24,
                      f25 = gf25,
                      f26 = gf26)

gy1 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH5[epi_trim$X.Nucleus1..MeanIntensity.CH5 >= 0 & epi_trim$X.Nucleus1..MeanIntensity.CH5 <= 25]))
gy2 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH5[epi_trim$X.Nucleus1..MeanIntensity.CH5 > 25 & epi_trim$X.Nucleus1..MeanIntensity.CH5 <= 50]))
gy3 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH5[epi_trim$X.Nucleus1..MeanIntensity.CH5 > 50 & epi_trim$X.Nucleus1..MeanIntensity.CH5 <= 75]))
gy4 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH5[epi_trim$X.Nucleus1..MeanIntensity.CH5 > 75 & epi_trim$X.Nucleus1..MeanIntensity.CH5 <= 100]))
gy5 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH5[epi_trim$X.Nucleus1..MeanIntensity.CH5 > 100 & epi_trim$X.Nucleus1..MeanIntensity.CH5 <= 125]))
gy6 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH5[epi_trim$X.Nucleus1..MeanIntensity.CH5 > 125 & epi_trim$X.Nucleus1..MeanIntensity.CH5 <= 150]))
gy7 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH5[epi_trim$X.Nucleus1..MeanIntensity.CH5 > 150 & epi_trim$X.Nucleus1..MeanIntensity.CH5 <= 175]))
gy8 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH5[epi_trim$X.Nucleus1..MeanIntensity.CH5 > 175 & epi_trim$X.Nucleus1..MeanIntensity.CH5 <= 200]))
gy9 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH5[epi_trim$X.Nucleus1..MeanIntensity.CH5 > 200 & epi_trim$X.Nucleus1..MeanIntensity.CH5 <= 225]))
gy10 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH5[epi_trim$X.Nucleus1..MeanIntensity.CH5 > 225 & epi_trim$X.Nucleus1..MeanIntensity.CH5 <= 250]))
gy11 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH5[epi_trim$X.Nucleus1..MeanIntensity.CH5 > 250 & epi_trim$X.Nucleus1..MeanIntensity.CH5 <= 275]))
gy12 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH5[epi_trim$X.Nucleus1..MeanIntensity.CH5 > 275 & epi_trim$X.Nucleus1..MeanIntensity.CH5 <= 300]))
gy13 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH5[epi_trim$X.Nucleus1..MeanIntensity.CH5 > 300 & epi_trim$X.Nucleus1..MeanIntensity.CH5 <= 325]))
gy14 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH5[epi_trim$X.Nucleus1..MeanIntensity.CH5 > 325 & epi_trim$X.Nucleus1..MeanIntensity.CH5 <= 350]))
gy15 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH5[epi_trim$X.Nucleus1..MeanIntensity.CH5 > 350 & epi_trim$X.Nucleus1..MeanIntensity.CH5 <= 375]))
gy16 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH5[epi_trim$X.Nucleus1..MeanIntensity.CH5 > 375 & epi_trim$X.Nucleus1..MeanIntensity.CH5 <= 400]))
gy17 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH5[epi_trim$X.Nucleus1..MeanIntensity.CH5 > 400 & epi_trim$X.Nucleus1..MeanIntensity.CH5 <= 425]))
gy18 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH5[epi_trim$X.Nucleus1..MeanIntensity.CH5 > 425 & epi_trim$X.Nucleus1..MeanIntensity.CH5 <= 450]))
gy19 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH5[epi_trim$X.Nucleus1..MeanIntensity.CH5 > 450 & epi_trim$X.Nucleus1..MeanIntensity.CH5 <= 475]))
gy20 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH5[epi_trim$X.Nucleus1..MeanIntensity.CH5 > 475 & epi_trim$X.Nucleus1..MeanIntensity.CH5 <= 500]))
gy21 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH5[epi_trim$X.Nucleus1..MeanIntensity.CH5 > 500]))

YAP_DATA = data.frame(y1 = gy1,
                      y2 = gy2,
                      y3 = gy3,
                      y4 = gy4,
                      y5 = gy5,
                      y6 = gy6,
                      y7 = gy7,
                      y8 = gy8,
                      y9 = gy9,
                      y10 = gy10,
                      y11 = gy11,
                      y12 = gy12,
                      y13 = gy13,
                      y14 = gy14,
                      y15 = gy15,
                      y16 = gy16,
                      y17 = gy17,
                      y18 = gy18,
                      y19 = gy19,
                      y20 = gy20,
                      y21 = gy21)

gp1 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH3[epi_trim$X.Nucleus1..MeanIntensity.CH3 >= 0 & epi_trim$X.Nucleus1..MeanIntensity.CH3 <= 25]))
gp2 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH3[epi_trim$X.Nucleus1..MeanIntensity.CH3 > 25 & epi_trim$X.Nucleus1..MeanIntensity.CH3 <= 50]))
gp3 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH3[epi_trim$X.Nucleus1..MeanIntensity.CH3 > 50 & epi_trim$X.Nucleus1..MeanIntensity.CH3 <= 75]))
gp4 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH3[epi_trim$X.Nucleus1..MeanIntensity.CH3 > 75 & epi_trim$X.Nucleus1..MeanIntensity.CH3 <= 100]))
gp5 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH3[epi_trim$X.Nucleus1..MeanIntensity.CH3 > 100 & epi_trim$X.Nucleus1..MeanIntensity.CH3 <= 125]))
gp6 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH3[epi_trim$X.Nucleus1..MeanIntensity.CH3 > 125 & epi_trim$X.Nucleus1..MeanIntensity.CH3 <= 150]))
gp7 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH3[epi_trim$X.Nucleus1..MeanIntensity.CH3 > 150 & epi_trim$X.Nucleus1..MeanIntensity.CH3 <= 175]))
gp8 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH3[epi_trim$X.Nucleus1..MeanIntensity.CH3 > 175 & epi_trim$X.Nucleus1..MeanIntensity.CH3 <= 200]))
gp9 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH3[epi_trim$X.Nucleus1..MeanIntensity.CH3 > 200 & epi_trim$X.Nucleus1..MeanIntensity.CH3 <= 225]))
gp10 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH3[epi_trim$X.Nucleus1..MeanIntensity.CH3 > 225 & epi_trim$X.Nucleus1..MeanIntensity.CH3 <= 250]))
gp11 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH3[epi_trim$X.Nucleus1..MeanIntensity.CH3 > 250 & epi_trim$X.Nucleus1..MeanIntensity.CH3 <= 275]))
gp12 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH3[epi_trim$X.Nucleus1..MeanIntensity.CH3 > 275 & epi_trim$X.Nucleus1..MeanIntensity.CH3 <= 300]))
gp13 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH3[epi_trim$X.Nucleus1..MeanIntensity.CH3 > 300 & epi_trim$X.Nucleus1..MeanIntensity.CH3 <= 325]))
gp14 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH3[epi_trim$X.Nucleus1..MeanIntensity.CH3 > 325 & epi_trim$X.Nucleus1..MeanIntensity.CH3 <= 350]))
gp15 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH3[epi_trim$X.Nucleus1..MeanIntensity.CH3 > 350 & epi_trim$X.Nucleus1..MeanIntensity.CH3 <= 375]))
gp16 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH3[epi_trim$X.Nucleus1..MeanIntensity.CH3 > 375 & epi_trim$X.Nucleus1..MeanIntensity.CH3 <= 400]))
gp17 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH3[epi_trim$X.Nucleus1..MeanIntensity.CH3 > 400 & epi_trim$X.Nucleus1..MeanIntensity.CH3 <= 425]))
gp18 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH3[epi_trim$X.Nucleus1..MeanIntensity.CH3 > 425 & epi_trim$X.Nucleus1..MeanIntensity.CH3 <= 450]))
gp19 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH3[epi_trim$X.Nucleus1..MeanIntensity.CH3 > 450 & epi_trim$X.Nucleus1..MeanIntensity.CH3 <= 475]))
gp20 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH3[epi_trim$X.Nucleus1..MeanIntensity.CH3 > 475 & epi_trim$X.Nucleus1..MeanIntensity.CH3 <= 500]))
gp21 = nrow(as.data.frame(epi_trim$X.Nucleus1..MeanIntensity.CH3[epi_trim$X.Nucleus1..MeanIntensity.CH3 > 500]))

PARP_DATA = data.frame(p1 = gp1,
                       p2 = gp2,
                       p3 = gp3,
                       p4 = gp4,
                       p5 = gp5,
                       p6 = gp6,
                       p7 = gp7,
                       p8 = gp8,
                       p9 = gp9,
                       p10 = gp10,
                       p11 = gp11,
                       p12 = gp12,
                       p13 = gp13,
                       p14 = gp14,
                       p15 = gp15,
                       p16 = gp16,
                       p17 = gp17,
                       p18 = gp18,
                       p19 = gp19,
                       p20 = gp20,
                       p21 = gp21)

write.csv(FRA_DATA, "/~/Analysis/Samplex/Samplex_aFRA_DATA.csv", row.names=FALSE)
write.csv(YAP_DATA, "/~/Analysis/Samplex/Samplex_aYAP_DATA.csv", row.names=FALSE)
write.csv(PARP_DATA, "/~/Analysis/Samplex/Samplex_aPARP_DATA.csv", row.names=FALSE)

chi_sq_results = data.frame(ESP = c(espc$p.value),
                            ESF = c(esfc$p.value),
                            ESY = c(esyc$p.value),
                            E8P = c(e8pc$p.value),
                            E8F = c(e8fc$p.value),
                            E8Y = c(e8yc$p.value),
                            EMP = c(empc$p.value),
                            EMF = c(emfc$p.value),
                            EMY = c(emyc$p.value))

write.csv(chi_sq_results, "/~/Analysis/Samplex/Samplex_a_chi_sq_results.csv", row.names=FALSE)

