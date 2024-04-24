##Mapping of VPN synapses to DN model nodes for dendrograms in figure 2(B-F left)

library(nat)
library(neuprintr)
library(natverse)
library(tidyverse)
library(catmaid)
library(hemibrainr)
library(fafbseg)
library(reticulate)
library(elmr)
library(gargle)
library(nat.flybrains)
library(nat.jrcbrains)
library(nat.nblast)
library(shiny)
options(max.print = 10000)
library(here)

# Define file names
file_names <- c("DNp01_flywire_syn_labeled_um.xlsx",
                "DNp02_flywire_syn_labeled_um.xlsx",
                "DNp03_flywire_syn_labeled_um.xlsx",
                "DNp04_flywire_syn_labeled_um.xlsx",
                "DNp06_flywire_syn_labeled_um.xlsx")

# Read each Excel file into a list of data frames
data_list <- lapply(file_names, function(file) {
  read_excel(here(file))
})

# Assign names to the data frames in the list
names(data_list) <- c("DNp01_syn_fly", "DNp02_syn_fly", "DNp03_syn_fly", "DNp04_syn_fly", "DNp06_syn_fly")

# Now you can access each data frame using their respective names
DNp01_syn_fly <- data_list$DNp01_syn_fly
DNp02_syn_fly <- data_list$DNp02_syn_fly
DNp03_syn_fly <- data_list$DNp03_syn_fly
DNp04_syn_fly <- data_list$DNp04_syn_fly
DNp06_syn_fly <- data_list$DNp06_syn_fly

VPNs = c("LC4","LC6", "LC22", "LPLC1", "LPLC2", "LPLC4")
DNp01_syn_fly_updated <- DNp01_syn_fly[DNp01_syn_fly$type %in% VPNs,]
DNp01_syn_fly_updated <- DNp01_syn_fly_updated %>%
  arrange(type)

DNp02_syn_fly_updated <- DNp02_syn_fly[DNp02_syn_fly$type %in% VPNs,]
DNp02_syn_fly_updated <- DNp02_syn_fly_updated %>%
  arrange(type)

DNp03_syn_fly_updated <- DNp03_syn_fly[DNp03_syn_fly$type %in% VPNs,]
DNp03_syn_fly_updated <- DNp03_syn_fly_updated %>%
  arrange(type)

DNp04_syn_fly_updated <- DNp04_syn_fly[DNp04_syn_fly$type %in% VPNs,]
DNp04_syn_fly_updated <- DNp04_syn_fly_updated %>%
  arrange(type)

DNp06_syn_fly_updated <- DNp06_syn_fly[DNp06_syn_fly$type %in% VPNs,]
DNp06_syn_fly_updated <- DNp06_syn_fly_updated %>%
  arrange(type)

synapseLocations_DNp01 <- DNp01_syn_fly_updated[, c("post_x", "post_y", "post_z")]
synapseLocations_DNp02 <- DNp02_syn_fly_updated[, c("post_x", "post_y", "post_z")]
synapseLocations_DNp03 <- DNp03_syn_fly_updated[, c("post_x", "post_y", "post_z")]
synapseLocations_DNp04 <- DNp04_syn_fly_updated[, c("post_x", "post_y", "post_z")]
synapseLocations_DNp06 <- DNp06_syn_fly_updated[, c("post_x", "post_y", "post_z")]

#change model here to map
DNp01 <- "DNp01_um_model.swc"
DNp02 <- "DNp02_um_model.swc"
DNp03 <- "DNp03_um_model.swc"
DNp04 <- "DNp04_um_model.swc"
DNp06 <- "DNp06_um_model.swc"

DNp01_skel = read.neuron(file.path(getwd(), DNp01))
DNp02_skel = read.neuron(file.path(getwd(), DNp02))
DNp03_skel = read.neuron(file.path(getwd(), DNp03))
DNp04_skel = read.neuron(file.path(getwd(), DNp04))
DNp06_skel = read.neuron(file.path(getwd(), DNp06))

near = nabor::knn(query= xyzmatrix(synapseLocations_DNp01),data=nat::xyzmatrix(DNp01_skel$d),k=1)$nn.idx
LC4_DNp01 = near[0:544]
LPLC2_DNp01 = near[545:1122, ]

near = nabor::knn(query= xyzmatrix(synapseLocations_DNp02),data=nat::xyzmatrix(DNp02_skel$d),k=1)$nn.idx
LC4_DNp02 = near[0:792,]

near = nabor::knn(query= xyzmatrix(synapseLocations_DNp03),data=nat::xyzmatrix(DNp03_skel$d),k=1)$nn.idx
LC4_DNp03 = near[221:747]
LC22_DNp03 = near[0:220]
LPLC1_DNp03 = near[747:1849]
LPLC2_DNp03 = near[1850:1865]
LPLC4_DNp03 = near[1866:3089]

near = nabor::knn(query= xyzmatrix(synapseLocations_DNp04),data=nat::xyzmatrix(DNp04_skel$d),k=1)$nn.idx
LC4_DNp04 = near[0:1587]
LPLC1_DNp04 = near[1588:1598, ]
LPLC2_DNp04 = near[1599:2444, ]
plot3d(xyzmatrix(synapseLocations_DNp04))
plot3d(DNp04_AM_trace, WithNodes = FALSE)

near = nabor::knn(query= xyzmatrix(synapseLocations_DNp06),data=nat::xyzmatrix(DNp06_skel$d),k=1)$nn.idx
LC4_DNp06 = near[0:99]
LC6_DNp06 = near[100:151]
LPLC1_DNp06 = near[0:693]
LPLC2_DNp06 = near[694:1337]


