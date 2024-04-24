#Figure 2(B-F left) and Figure 4(A-B)

#Import libraries
library(catmaid)
library(neuprintr)
library(hemibrainr)
library(natverse)
library(plotly)
library(dplyr)
library(ggplot2)
library(elmr)
library(csv)
library(rgl)
library(readxl)

# Loading in the receptive field estimated COMs

LC4_receptive_field <- read_excel("LC4_receptivefield_COM.xlsx")
LC6_receptive_field <- read_excel("LC6_receptivefield_COM.xlsx")
LC22_receptive_field <- read_excel("LC22_receptivefield_COM.xlsx")
LPLC1_receptive_field <- read_excel("LPLC1_receptivefield_COM.xlsx")
LPLC2_receptive_field <- read_excel("LPLC2_receptivefield_COM.xlsx")
LPLC4_receptive_field <- read_excel("LPLC4_receptivefield_COM.xlsx")

# Loading in the synapse data
DNp01_syn_fly <- read_excel("DNp01_flywire_syn_labeled_um.xlsx")
DNp02_syn_fly <- read_excel("DNp02_flywire_syn_labeled_um.xlsx")
DNp03_syn_fly <- read_excel("DNp03_flywire_syn_labeled_um.xlsx")
DNp04_syn_fly <- read_excel("DNp04_flywire_syn_labeled_um.xlsx")
DNp06_syn_fly <- read_excel("DNp06_flywire_syn_labeled_um.xlsx")

##DN mesh data
choose_segmentation("flywire")
DNp01_mesh = read_cloudvolume_meshes("720575940622838154")
DNp02_mesh = read_cloudvolume_meshes('720575940619654053')
DNp03_mesh = read_cloudvolume_meshes("720575940627645514")
DNp04_mesh = read_cloudvolume_meshes("720575940604954289")
DNp06_mesh = read_cloudvolume_meshes("720575940622673860")

# Subsetting all synapses to only VPN synapses
VPNs = c("LC4", "LC6", "LC22", "LPLC1", "LPLC2", "LPLC4")
DNp01_syn_fly_points <- DNp01_syn_fly[DNp01_syn_fly$type %in% VPNs,]
DNp02_syn_fly_points <- DNp02_syn_fly[DNp02_syn_fly$type %in% VPNs,]
DNp03_syn_fly_points <- DNp03_syn_fly[DNp03_syn_fly$type %in% VPNs,]
DNp04_syn_fly_points <- DNp04_syn_fly[DNp04_syn_fly$type %in% VPNs,]
DNp06_syn_fly_points <- DNp06_syn_fly[DNp06_syn_fly$type %in% VPNs,]

# Getting the top 10 neurons in the anterior/posterior axis and the dorsal/ventral axis for LC4
DV_LC4 <- LC4_receptive_field[with(LC4_receptive_field,order(-DV_norm)),]
AP_LC4 <- LC4_receptive_field[with(LC4_receptive_field,order(-AP_norm)),]
Ventral_LC4 <- DV_LC4[1:10,]
Dorsal_LC4 <- DV_LC4[45:54,]
Anterior_LC4 <- AP_LC4[1:10,]
Posterior_LC4 <- AP_LC4[45:54,]
Ventral_ids = Ventral_LC4$updated_ids
Dorsal_ids = Dorsal_LC4$updated_ids

#Same as above but for LC6
DV_LC6 <- LC6_receptive_field[with(LC6_receptive_field,order(-DV_norm)),]
AP_LC6 <- LC6_receptive_field[with(LC6_receptive_field,order(-AP_norm)),]
Ventral_LC6 <- DV_LC6[1:10,]
Dorsal_LC6 <- DV_LC6[53:62,]
Anterior_LC6 <- AP_LC6[1:10,]
Posterior_LC6 <- AP_LC6[53:62,]

#Same as above but for LC22
DV_LC22 <- LC22_receptive_field[with(LC22_receptive_field,order(-DV_norm)),]
AP_LC22 <- LC22_receptive_field[with(LC22_receptive_field,order(-AP_norm)),]
Ventral_LC22 <- DV_LC22[1:10,]
Dorsal_LC22 <- DV_LC22[52:61,]
Anterior_LC22 <- AP_LC22[1:10,]
Posterior_LC22 <- AP_LC22[52:61,]


#Same as above but for LPLC1
DV_LPLC1 <- LPLC1_receptive_field[with(LPLC1_receptive_field,order(-DV_norm)),]
AP_LPLC1 <- LPLC1_receptive_field[with(LPLC1_receptive_field,order(-AP_norm)),]
Ventral_LPLC1 <- DV_LPLC1[1:10,]
Dorsal_LPLC1 <- DV_LPLC1[56:65,]
Anterior_LPLC1 <- AP_LPLC1[1:10,]
Posterior_LPLC1 <- AP_LPLC1[56:65,]

#Same as above but for LPLC2
DV_LPLC2 <- LPLC2_receptive_field[with(LPLC2_receptive_field,order(-DV_norm)),]
AP_LPLC2 <- LPLC2_receptive_field[with(LPLC2_receptive_field,order(-AP_norm)),]
Ventral_LPLC2 <- DV_LPLC2[1:5,]
Dorsal_LPLC2 <- DV_LPLC2[103:107,]
Anterior_LPLC2 <- AP_LPLC2[1:10,]
Posterior_LPLC2 <- AP_LPLC2[98:107,]

#Same as above but for LPLC4
DV_LPLC4 <- LPLC4_receptive_field[with(LPLC4_receptive_field,order(-DV_norm)),]
AP_LPLC4 <- LPLC4_receptive_field[with(LPLC4_receptive_field,order(-AP_norm)),]
Ventral_LPLC4 <- DV_LPLC4[1:10,]
Dorsal_LPLC4 <- DV_LPLC4[45:54,]
Anterior_LPLC4 <- AP_LPLC4[1:10,]
Posterior_LPLC4 <- AP_LPLC4[45:54,]

#Assign the skeletons to read, just assign the VPN population loaded above, however you can also read the mesh directly from the flywire website,
#note that some may not plot if the ID has changed since the skeletonization, so use the meshes down below in case of that, or go back and re-skeletonize
neu <- LPLC4

for (id in Anterior_LPLC4$updated_ids) {
  plot3d(neu[[id]], col = "red", WithNodes = FALSE, axes = FALSE, xlab = "", ylab = "", zlab = "")
}
for (id in Posterior_LPLC4$updated_ids) {
  plot3d(neu[[id]], col = "blue", WithNodes = FALSE, axes = FALSE, xlab = "", ylab = "", zlab = "")
}

#Utilizing the meshes for each population: Takes approximately 20-30 seconds per command.
#This is purely for visualization of populations in each axis.
#LC4
ant_LC4 = read_cloudvolume_meshes(Anterior_LC4$updated_ids)
post_LC4 = read_cloudvolume_meshes(Posterior_LC4$updated_ids)
dors_LC4 = read_cloudvolume_meshes(Dorsal_LC4$updated_ids)
vent_LC4 = read_cloudvolume_meshes(Ventral_LC4$updated_ids)

#LC6
ant_LC6 = read_cloudvolume_meshes(Anterior_LC6$updated_ids)
post_LC6 = read_cloudvolume_meshes(Posterior_LC6$updated_ids)
dors_LC6 = read_cloudvolume_meshes(Dorsal_LC6$updated_ids)
vent_LC6 = read_cloudvolume_meshes(Ventral_LC6$updated_ids)

#LC22
ant_LC22 = read_cloudvolume_meshes(Anterior_LC22$updated_ids)
post_LC22 = read_cloudvolume_meshes(Posterior_LC22$updated_ids)
dors_LC22 = read_cloudvolume_meshes(Dorsal_LC22$updated_ids)
vent_LC22 = read_cloudvolume_meshes(Ventral_LC22$updated_ids)

#LPLC1
ant_LPLC1 = read_cloudvolume_meshes(Anterior_LPLC1$updated_ids)
post_LPLC1 = read_cloudvolume_meshes(Posterior_LPLC1$updated_ids)
dors_LPLC1 = read_cloudvolume_meshes(Dorsal_LPLC1$updated_ids)
vent_LPLC1 = read_cloudvolume_meshes(Ventral_LPLC1$updated_ids)

#LPLC2
ant_LPLC2 = read_cloudvolume_meshes(Anterior_LPLC2$updated_ids)
post_LPLC2 = read_cloudvolume_meshes(Posterior_LPLC2$updated_ids)
dors_LPLC2 = read_cloudvolume_meshes(Dorsal_LPLC2$updated_ids)
vent_LPLC2 = read_cloudvolume_meshes(Ventral_LPLC2$updated_ids)

#LPLC4
ant_LPLC4 = read_cloudvolume_meshes(Anterior_LPLC4$updated_ids)
post_LPLC4 = read_cloudvolume_meshes(Posterior_LPLC4$updated_ids)
dors_LPLC4 = read_cloudvolume_meshes(Dorsal_LPLC4$updated_ids)
vent_LPLC4 = read_cloudvolume_meshes(Ventral_LPLC4$updated_ids)

#Plotting of the obtained meshes, just change the cell type, note this is only 10 from each directional extreme of the lobula

plot3d(ant_LPLC4, col = 'purple')
plot3d(post_LPLC4, col = 'cyan')
plot3d(dors_LC4, col = 'red')
plot3d(vent_LC4, col = 'blue')
plot3d(FAFB14)
#rgl.snapshot(filename = "change_filename_here.png",fmt = "png")


### Temp work area to plot the synapses to each DN by their axis ###

# Update the IDs from the receptive fields, here you can change which direction you want to look at (do one at a time to keep track)
ids <- Posterior_LPLC4$updated_ids
ids_updated = flywire_updateids(ids)

#LC4
Dorsal_LC4<- cbind(select(Dorsal_LC4, -updated_ids), updated_ids = ids_updated)
Anterior_LC4<- cbind(select(Anterior_LC4, -updated_ids), updated_ids = ids_updated)
Ventral_LC4<- cbind(select(Ventral_LC4, -updated_ids), updated_ids = ids_updated)
Posterior_LC4<- cbind(select(Posterior_LC4, -updated_ids), updated_ids = ids_updated)

#LC6
Dorsal_LC6<- cbind(select(Dorsal_LC6, -updated_ids), updated_ids = ids_updated)
Anterior_LC6<- cbind(select(Anterior_LC6, -updated_ids), updated_ids = ids_updated)
Ventral_LC6<- cbind(select(Ventral_LC6, -updated_ids), updated_ids = ids_updated)
Posterior_LC6<- cbind(select(Posterior_LC6, -updated_ids), updated_ids = ids_updated)

#LC22
Dorsal_LC22<- cbind(select(Dorsal_LC22, -updated_ids), updated_ids = ids_updated)
Anterior_LC22<- cbind(select(Anterior_LC22, -updated_ids), updated_ids = ids_updated)
Ventral_LC22<- cbind(select(Ventral_LC22, -updated_ids), updated_ids = ids_updated)
Posterior_LC22<- cbind(select(Posterior_LC22, -updated_ids), updated_ids = ids_updated)

#LPLC1
Dorsal_LPLC1<- cbind(select(Dorsal_LPLC1, -updated_ids), updated_ids = ids_updated)
Anterior_LPLC1<- cbind(select(Anterior_LPLC1, -updated_ids), updated_ids = ids_updated)
Ventral_LPLC1<- cbind(select(Ventral_LPLC1, -updated_ids), updated_ids = ids_updated)
Posterior_LPLC1<- cbind(select(Posterior_LPLC1, -updated_ids), updated_ids = ids_updated)

#LPLC2
Dorsal_LPLC2<- cbind(select(Dorsal_LPLC2, -updated_ids), updated_ids = ids_updated)
Anterior_LPLC2<- cbind(select(Anterior_LPLC2, -updated_ids), updated_ids = ids_updated)
Ventral_LPLC2<- cbind(select(Ventral_LPLC2, -updated_ids), updated_ids = ids_updated)
Posterior_LPLC2<- cbind(select(Posterior_LPLC2, -updated_ids), updated_ids = ids_updated)

#LPLC4
Dorsal_LPLC4<- cbind(select(Dorsal_LPLC4, -updated_ids), updated_ids = ids_updated)
Anterior_LPLC4<- cbind(select(Anterior_LPLC4, -updated_ids), updated_ids = ids_updated)
Ventral_LPLC4<- cbind(select(Ventral_LPLC4, -updated_ids), updated_ids = ids_updated)
Posterior_LPLC4<- cbind(select(Posterior_LPLC4, -updated_ids), updated_ids = ids_updated)

#Update the IDs from the synapse data
DNp01_ids <-DNp01_syn_fly_points$pre
DNp01_ids_updated = flywire_updateids(DNp01_ids)
DNp01_syn_fly_points<- cbind(select(DNp01_syn_fly_points, -pre), updated_ids = DNp01_ids_updated)

DNp02_ids <-DNp02_syn_fly_points$pre
DNp02_ids_updated = flywire_updateids(DNp02_ids)
DNp02_syn_fly_points<- cbind(select(DNp02_syn_fly_points, -pre), updated_ids = DNp02_ids_updated)

DNp03_ids <-DNp03_syn_fly_points$pre
DNp03_ids_updated = flywire_updateids(DNp03_ids)
DNp03_syn_fly_points<- cbind(select(DNp03_syn_fly_points, -pre), updated_ids = DNp03_ids_updated)

DNp04_ids <-DNp04_syn_fly_points$pre
DNp04_ids_updated = flywire_updateids(DNp04_ids)
DNp04_syn_fly_points<- cbind(select(DNp04_syn_fly_points, -pre), updated_ids = DNp04_ids_updated)

DNp06_ids <-DNp06_syn_fly_points$pre
DNp06_ids_updated = flywire_updateids(DNp06_ids)
DNp06_syn_fly_points<- cbind(select(DNp06_syn_fly_points, -pre), updated_ids = DNp06_ids_updated)

## Subsetting the VPN synapses for DN of interest the IDs from the receptive fields.
#DNp01 VPNS: LC4 and LPLC2
Dorsal_LC4_DNp01 <- DNp01_syn_fly_points[DNp01_syn_fly_points$updated_ids %in% Dorsal_LC4$updated_ids, ]
Ventral_LC4_DNp01 <- DNp01_syn_fly_points[DNp01_syn_fly_points$updated_ids %in% Ventral_LC4$updated_ids, ]
Anterior_LC4_DNp01 <- DNp01_syn_fly_points[DNp01_syn_fly_points$updated_ids %in% Anterior_LC4$updated_ids, ]
Posterior_LC4_DNp01 <- DNp01_syn_fly_points[DNp01_syn_fly_points$updated_ids %in% Posterior_LC4$updated_ids, ]

Dorsal_LPLC2_DNp01 <- DNp01_syn_fly_points[DNp01_syn_fly_points$updated_ids %in% Dorsal_LPLC2$updated_ids, ]
Ventral_LPLC2_DNp01 <- DNp01_syn_fly_points[DNp01_syn_fly_points$updated_ids %in% Ventral_LPLC2$updated_ids, ]
Anterior_LPLC2_DNp01 <- DNp01_syn_fly_points[DNp01_syn_fly_points$updated_ids %in% Anterior_LPLC2$updated_ids, ]
Posterior_LPLC2_DNp01 <- DNp01_syn_fly_points[DNp01_syn_fly_points$updated_ids %in% Posterior_LPLC2$updated_ids, ]

#DNp02 VPNS: LC4 
Dorsal_LC4_DNp02 <- DNp02_syn_fly_points[DNp02_syn_fly_points$updated_ids %in% Dorsal_LC4$updated_ids, ]
Ventral_LC4_DNp02 <- DNp02_syn_fly_points[DNp02_syn_fly_points$updated_ids %in% Ventral_LC4$updated_ids, ]
Anterior_LC4_DNp02 <- DNp02_syn_fly_points[DNp02_syn_fly_points$updated_ids %in% Anterior_LC4$updated_ids, ]
Posterior_LC4_DNp02 <- DNp02_syn_fly_points[DNp02_syn_fly_points$updated_ids %in% Posterior_LC4$updated_ids, ]

#DNp03 VPNS: LC4, LC22, LPLC1, and LPLC4
Dorsal_LC4_DNp03 <- DNp03_syn_fly_points[DNp03_syn_fly_points$updated_ids %in% Dorsal_LC4$updated_ids, ]
Ventral_LC4_DNp03 <- DNp03_syn_fly_points[DNp03_syn_fly_points$updated_ids %in% Ventral_LC4$updated_ids, ]
Anterior_LC4_DNp03 <- DNp03_syn_fly_points[DNp03_syn_fly_points$updated_ids %in% Anterior_LC4$updated_ids, ]
Posterior_LC4_DNp03 <- DNp03_syn_fly_points[DNp03_syn_fly_points$updated_ids %in% Posterior_LC4$updated_ids, ]

Dorsal_LC22_DNp03 <- DNp03_syn_fly_points[DNp03_syn_fly_points$updated_ids %in% Dorsal_LC22$updated_ids, ]
Ventral_LC22_DNp03 <- DNp03_syn_fly_points[DNp03_syn_fly_points$updated_ids %in% Ventral_LC22$updated_ids, ]
Anterior_LC22_DNp03 <- DNp03_syn_fly_points[DNp03_syn_fly_points$updated_ids %in% Anterior_LC22$updated_ids, ]
Posterior_LC22_DNp03 <- DNp03_syn_fly_points[DNp03_syn_fly_points$updated_ids %in% Posterior_LC22$updated_ids, ]

Dorsal_LPLC1_DNp03 <- DNp03_syn_fly_points[DNp03_syn_fly_points$updated_ids %in% Dorsal_LPLC1$updated_ids, ]
Ventral_LPLC1_DNp03 <- DNp03_syn_fly_points[DNp03_syn_fly_points$updated_ids %in% Ventral_LPLC1$updated_ids, ]
Anterior_LPLC1_DNp03 <- DNp03_syn_fly_points[DNp03_syn_fly_points$updated_ids %in% Anterior_LPLC1$updated_ids, ]
Posterior_LPLC1_DNp03 <- DNp03_syn_fly_points[DNp03_syn_fly_points$updated_ids %in% Posterior_LPLC1$updated_ids, ]

Dorsal_LPLC4_DNp03 <- DNp03_syn_fly_points[DNp03_syn_fly_points$updated_ids %in% Dorsal_LPLC4$updated_ids, ]
Ventral_LPLC4_DNp03 <- DNp03_syn_fly_points[DNp03_syn_fly_points$updated_ids %in% Ventral_LPLC4$updated_ids, ]
Anterior_LPLC4_DNp03 <- DNp03_syn_fly_points[DNp03_syn_fly_points$updated_ids %in% Anterior_LPLC4$updated_ids, ]
Posterior_LPLC4_DNp03 <- DNp03_syn_fly_points[DNp03_syn_fly_points$updated_ids %in% Posterior_LPLC4$updated_ids, ]

#DNp04 VPNS: LC4, and LPLC2
Dorsal_LC4_DNp04 <- DNp04_syn_fly_points[DNp04_syn_fly_points$updated_ids %in% Dorsal_LC4$updated_ids, ]
Ventral_LC4_DNp04 <- DNp04_syn_fly_points[DNp04_syn_fly_points$updated_ids %in% Ventral_LC4$updated_ids, ]
Anterior_LC4_DNp04 <- DNp04_syn_fly_points[DNp04_syn_fly_points$updated_ids %in% Anterior_LC4$updated_ids, ]
Posterior_LC4_DNp04 <- DNp04_syn_fly_points[DNp04_syn_fly_points$updated_ids %in% Posterior_LC4$updated_ids, ]

Dorsal_LPLC2_DNp04 <- DNp04_syn_fly_points[DNp04_syn_fly_points$updated_ids %in% Dorsal_LPLC2$updated_ids, ]
Ventral_LPLC2_DNp04 <- DNp04_syn_fly_points[DNp04_syn_fly_points$updated_ids %in% Ventral_LPLC2$updated_ids, ]
Anterior_LPLC2_DNp04 <- DNp04_syn_fly_points[DNp04_syn_fly_points$updated_ids %in% Anterior_LPLC2$updated_ids, ]
Posterior_LPLC2_DNp04 <- DNp04_syn_fly_points[DNp04_syn_fly_points$updated_ids %in% Posterior_LPLC2$updated_ids, ]

# DNp06 VPNS: LC4, LC6, LPLC1, and LPLC2
Dorsal_LC4_DNp06 <- DNp06_syn_fly_points[DNp06_syn_fly_points$updated_ids %in% Dorsal_LC4$updated_ids, ]
Ventral_LC4_DNp06 <- DNp06_syn_fly_points[DNp06_syn_fly_points$updated_ids %in% Ventral_LC4$updated_ids, ]
Anterior_LC4_DNp06 <- DNp06_syn_fly_points[DNp06_syn_fly_points$updated_ids %in% Anterior_LC4$updated_ids, ]
Posterior_LC4_DNp06 <- DNp06_syn_fly_points[DNp06_syn_fly_points$updated_ids %in% Posterior_LC4$updated_ids, ]

Dorsal_LC6_DNp06 <- DNp06_syn_fly_points[DNp06_syn_fly_points$updated_ids %in% Dorsal_LC6$updated_ids, ]
Ventral_LC6_DNp06 <- DNp06_syn_fly_points[DNp06_syn_fly_points$updated_ids %in% Ventral_LC6$updated_ids, ]
Anterior_LC6_DNp06 <- DNp06_syn_fly_points[DNp06_syn_fly_points$updated_ids %in% Anterior_LC6$updated_ids, ]
Posterior_LC6_DNp06 <- DNp06_syn_fly_points[DNp06_syn_fly_points$updated_ids %in% Posterior_LC6$updated_ids, ]

Dorsal_LPLC1_DNp06 <- DNp06_syn_fly_points[DNp06_syn_fly_points$updated_ids %in% Dorsal_LPLC1$updated_ids, ]
Ventral_LPLC1_DNp06 <- DNp06_syn_fly_points[DNp06_syn_fly_points$updated_ids %in% Ventral_LPLC1$updated_ids, ]
Anterior_LPLC1_DNp06 <- DNp06_syn_fly_points[DNp06_syn_fly_points$updated_ids %in% Anterior_LPLC1$updated_ids, ]
Posterior_LPLC1_DNp06 <- DNp06_syn_fly_points[DNp06_syn_fly_points$updated_ids %in% Posterior_LPLC1$updated_ids, ]

Dorsal_LPLC2_DNp06 <- DNp06_syn_fly_points[DNp06_syn_fly_points$updated_ids %in% Dorsal_LPLC2$updated_ids, ]
Ventral_LPLC2_DNp06 <- DNp06_syn_fly_points[DNp06_syn_fly_points$updated_ids %in% Ventral_LPLC2$updated_ids, ]
Anterior_LPLC2_DNp06 <- DNp06_syn_fly_points[DNp06_syn_fly_points$updated_ids %in% Anterior_LPLC2$updated_ids, ]
Posterior_LPLC2_DNp06 <- DNp06_syn_fly_points[DNp06_syn_fly_points$updated_ids %in% Posterior_LPLC2$updated_ids, ]

# Counting of synapses across the neurons 
unique(Posterior_LPLC2_DNp06$updated_ids)

Posterior_LPLC2_syns_to_DNp06 = aggregate(data.frame(synapses = Posterior_LPLC2_DNp06$updated_ids), list(LPLC2_ID = Posterior_LPLC2_DNp06$updated_ids), length)
Posterior_LPLC2_syn_count_to_DNp06 = length(Posterior_LPLC2_DNp06$updated_ids) # 59 synapses 

#(Figure 4A-B, but remaining VPN-DNs below as well)
#Plotting of select synapses with a given neuron

# LC4 > DNp01
points3d(Dorsal_LC4_DNp01$post_x, Dorsal_LC4_DNp01$post_y, Dorsal_LC4_DNp01$post_z, col = 'red', size= 16)
points3d(Ventral_LC4_DNp01$post_x, Ventral_LC4_DNp01$post_y, Ventral_LC4_DNp01$post_z, col = 'blue', size= 16)
points3d(Anterior_LC4_DNp01$post_x, Anterior_LC4_DNp01$post_y, Anterior_LC4_DNp01$post_z, col = 'purple', size= 16)
points3d(Posterior_LC4_DNp01$post_x, Posterior_LC4_DNp01$post_y, Posterior_LC4_DNp01$post_z, col = 'cyan', size= 16)
plot3d(DNp01_mesh_fly/1000, col ='black', alpha = 0.3)

# LPLC2 > DNp01
points3d(Dorsal_LPLC2_DNp01$post_x, Dorsal_LPLC2_DNp01$post_y, Dorsal_LPLC2_DNp01$post_z, col = 'red', size= 16)
points3d(Ventral_LPLC2_DNp01$post_x, Ventral_LPLC2_DNp01$post_y, Ventral_LPLC2_DNp01$post_z, col = 'blue', size= 16)
points3d(Anterior_LPLC2_DNp01$post_x, Anterior_LPLC2_DNp01$post_y, Anterior_LPLC2_DNp01$post_z, col = 'purple', size= 16)
points3d(Posterior_LPLC2_DNp01$post_x, Posterior_LPLC2_DNp01$post_y, Posterior_LPLC2_DNp01$post_z, col = 'cyan', size= 16)

#DNp01 VPN zoomed in dendrites 
view3d(fov=0,zoom=0.25,userMatrix=rotationMatrix(50/180*pi,1,0,0) %*% rotationMatrix(30/180*pi,0,0,1) %*% rotationMatrix(-55/180*pi,0,1,0))
rgl.snapshot(filename = "LC4_DV_synapses_to_DNp01.png",fmt = "png")

#DNp03 VPN zoomed in dendrites 
# LPLC4 > DNp03
points3d(Dorsal_LPLC4_DNp03$post_x, Dorsal_LPLC4_DNp03$post_y, Dorsal_LPLC4_DNp03$post_z, col = 'red', size= 16)
points3d(Ventral_LPLC4_DNp03$post_x, Ventral_LPLC4_DNp03$post_y, Ventral_LPLC4_DNp03$post_z, col = 'blue', size= 16)
points3d(Anterior_LPLC4_DNp03$post_x, Anterior_LPLC4_DNp03$post_y, Anterior_LPLC4_DNp03$post_z, col = 'purple', size= 16)
points3d(Posterior_LPLC4_DNp03$post_x, Posterior_LPLC4_DNp03$post_y, Posterior_LPLC4_DNp03$post_z, col = 'cyan', size= 16)
plot3d(DNp03_mesh_fly/1000, col ='black', alpha = 0.3)
view3d(fov=10,zoom=0.22,userMatrix=rotationMatrix(31/180*pi,1,0,0) %*% rotationMatrix(20/180*pi,0,0,1) %*% rotationMatrix(-18/180*pi,0,1,0))
view3d(zoom =.22)

rgl.snapshot(filename = "LPLC4_DV_synapses_to_DNp03.png",fmt = "png")

###### Below is the synapse count per 10 VPNs per axis and the plotting of VPN synapses to other DNs. 
Dorsal_LC4_syn_count_to_DNp01 = length(Dorsal_LC4_DNp01$updated_ids) # 86 synapses
Ventral_LC4_syn_count_to_DNp01 = length(Ventral_LC4_DNp01$updated_ids) # 85 synapses
Anterior_LC4_syn_count_to_DNp01 = length(Anterior_LC4_DNp01$updated_ids) # 121 synapses
Posterior_LC4_syn_count_to_DNp01 = length(Posterior_LC4_DNp01$updated_ids) # 105 synapses



Dorsal_LPLC2_syn_count_to_DNp01 = length(Dorsal_LPLC2_DNp01$updated_ids) # 116 synapses
Ventral_LPLC2_syn_count_to_DNp01 = length(Ventral_LPLC2_DNp01$updated_ids) # 16 synapses
Anterior_LPLC2_syn_count_to_DNp01 = length(Anterior_LPLC2_DNp01$updated_ids) # 39 synapses
Posterior_LPLC2_syn_count_to_DNp01 = length(Posterior_LPLC2_DNp01$updated_ids) # 32 synapses

# LC4 > DNp02
points3d(Dorsal_LC4_DNp02$post_x, Dorsal_LC4_DNp02$post_y, Dorsal_LC4_DNp02$post_z, col = 'red', size= 10)
points3d(Ventral_LC4_DNp02$post_x, Ventral_LC4_DNp02$post_y, Ventral_LC4_DNp02$post_z, col = 'blue', size= 10)
points3d(Anterior_LC4_DNp02$post_x, Anterior_LC4_DNp02$post_y, Anterior_LC4_DNp02$post_z, col = 'purple', size= 10)
points3d(Posterior_LC4_DNp02$post_x, Posterior_LC4_DNp02$post_y, Posterior_LC4_DNp02$post_z, col = 'cyan', size= 10)

Dorsal_LC4_syn_count_to_DNp02 = length(Dorsal_LC4_DNp02$updated_ids) # 121 synapses
Ventral_LC4_syn_count_to_DNp02 = length(Ventral_LC4_DNp02$updated_ids) # 163 synapses
Anterior_LC4_syn_count_to_DNp02 = length(Anterior_LC4_DNp02$updated_ids) # 221 synapses
Posterior_LC4_syn_count_to_DNp02 = length(Posterior_LC4_DNp02$updated_ids) # 47 synapses

# LC4 > DNp03
points3d(Dorsal_LC4_DNp03$post_x, Dorsal_LC4_DNp03$post_y, Dorsal_LC4_DNp03$post_z, col = 'red', size= 10)
points3d(Ventral_LC4_DNp03$post_x, Ventral_LC4_DNp03$post_y, Ventral_LC4_DNp03$post_z, col = 'blue', size= 10)
points3d(Anterior_LC4_DNp03$post_x, Anterior_LC4_DNp03$post_y, Anterior_LC4_DNp03$post_z, col = 'purple', size= 10)
points3d(Posterior_LC4_DNp03$post_x, Posterior_LC4_DNp03$post_y, Posterior_LC4_DNp03$post_z, col = 'cyan', size= 10)

Dorsal_LC4_syn_count_to_DNp03 = length(Dorsal_LC4_DNp03$updated_ids) # 82 synapses
Ventral_LC4_syn_count_to_DNp03 = length(Ventral_LC4_DNp03$updated_ids) # 111 synapses
Anterior_LC4_syn_count_to_DNp03 = length(Anterior_LC4_DNp03$updated_ids) # 102 synapses
Posterior_LC4_syn_count_to_DNp03 = length(Posterior_LC4_DNp03$updated_ids) # 116 synapses

# LC22 > DNp03
points3d(Dorsal_LC22_DNp03$post_x, Dorsal_LC22_DNp03$post_y, Dorsal_LC22_DNp03$post_z, col = 'red', size= 10)
points3d(Ventral_LC22_DNp03$post_x, Ventral_LC22_DNp03$post_y, Ventral_LC22_DNp03$post_z, col = 'blue', size= 10)
points3d(Anterior_LC22_DNp03$post_x, Anterior_LC22_DNp03$post_y, Anterior_LC22_DNp03$post_z, col = 'purple', size= 10)
points3d(Posterior_LC22_DNp03$post_x, Posterior_LC22_DNp03$post_y, Posterior_LC22_DNp03$post_z, col = 'cyan', size= 10)

Dorsal_LC22_syn_count_to_DNp03 = length(Dorsal_LC22_DNp03$updated_ids) # 5 synapses
Ventral_LC22_syn_count_to_DNp03 = length(Ventral_LC22_DNp03$updated_ids) # 31 synapses
Anterior_LC22_syn_count_to_DNp03 = length(Anterior_LC22_DNp03$updated_ids) # 94 synapses
Posterior_LC22_syn_count_to_DNp03 = length(Posterior_LC22_DNp03$updated_ids) # 6 synapses

# LPLC1 > DNp03
points3d(Dorsal_LPLC1_DNp03$post_x, Dorsal_LPLC1_DNp03$post_y, Dorsal_LPLC1_DNp03$post_z, col = 'red', size= 10)
points3d(Ventral_LPLC1_DNp03$post_x, Ventral_LPLC1_DNp03$post_y, Ventral_LPLC1_DNp03$post_z, col = 'blue', size= 10)
points3d(Anterior_LPLC1_DNp03$post_x, Anterior_LPLC1_DNp03$post_y, Anterior_LPLC1_DNp03$post_z, col = 'purple', size= 10)
points3d(Posterior_LPLC1_DNp03$post_x, Posterior_LPLC1_DNp03$post_y, Posterior_LPLC1_DNp03$post_z, col = 'cyan', size= 10)

Dorsal_LPLC1_syn_count_to_DNp03 = length(Dorsal_LPLC1_DNp03$updated_ids) # 95 synapses
Ventral_LPLC1_syn_count_to_DNp03 = length(Ventral_LPLC1_DNp03$updated_ids) # 304 synapses
Anterior_LPLC1_syn_count_to_DNp03 = length(Anterior_LPLC1_DNp03$updated_ids) # 90 synapses
Posterior_LPLC1_syn_count_to_DNp03 = length(Posterior_LPLC1_DNp03$updated_ids) # 183 synapses

Dorsal_LPLC4_syn_count_to_DNp03 = length(Dorsal_LPLC4_DNp03$updated_ids) # 170 synapses
Ventral_LPLC4_syn_count_to_DNp03 = length(Ventral_LPLC4_DNp03$updated_ids) # 197 synapses
Anterior_LPLC4_syn_count_to_DNp03 = length(Anterior_LPLC4_DNp03$updated_ids) # 248 synapses
Posterior_LPLC4_syn_count_to_DNp03 = length(Posterior_LPLC4_DNp03$updated_ids) # 188 synapses

# LC4 > DNp04
points3d(Dorsal_LC4_DNp04$post_x, Dorsal_LC4_DNp04$post_y, Dorsal_LC4_DNp04$post_z, col = 'red', size= 10)
points3d(Ventral_LC4_DNp04$post_x, Ventral_LC4_DNp04$post_y, Ventral_LC4_DNp04$post_z, col = 'blue', size= 10)
points3d(Anterior_LC4_DNp04$post_x, Anterior_LC4_DNp04$post_y, Anterior_LC4_DNp04$post_z, col = 'purple', size= 10)
points3d(Posterior_LC4_DNp04$post_x, Posterior_LC4_DNp04$post_y, Posterior_LC4_DNp04$post_z, col = 'cyan', size= 10)

Dorsal_LC4_syn_count_to_DNp04 = length(Dorsal_LC4_DNp04$updated_ids) # 166 synapses
Ventral_LC4_syn_count_to_DNp04 = length(Ventral_LC4_DNp04$updated_ids) # 480 synapses
Anterior_LC4_syn_count_to_DNp04 = length(Anterior_LC4_DNp04$updated_ids) # 353 synapses
Posterior_LC4_syn_count_to_DNp04 = length(Posterior_LC4_DNp04$updated_ids) # 294 synapses

# LPLC2 > DNp04
points3d(Dorsal_LPLC2_DNp04$post_x, Dorsal_LPLC2_DNp04$post_y, Dorsal_LPLC2_DNp04$post_z, col = 'red', size= 10)
points3d(Ventral_LPLC2_DNp04$post_x, Ventral_LPLC2_DNp04$post_y, Ventral_LPLC2_DNp04$post_z, col = 'blue', size= 10)
points3d(Anterior_LPLC2_DNp04$post_x, Anterior_LPLC2_DNp04$post_y, Anterior_LPLC2_DNp04$post_z, col = 'purple', size= 10)
points3d(Posterior_LPLC2_DNp04$post_x, Posterior_LPLC2_DNp04$post_y, Posterior_LPLC2_DNp04$post_z, col = 'cyan', size= 10)

Dorsal_LPLC2_syn_count_to_DNp04 = length(Dorsal_LPLC2_DNp04$updated_ids) # 88 synapses
Ventral_LPLC2_syn_count_to_DNp04 = length(Ventral_LPLC2_DNp04$updated_ids) # 77 synapses
Anterior_LPLC2_syn_count_to_DNp04 = length(Anterior_LPLC2_DNp04$updated_ids) # 111 synapses
Posterior_LPLC2_syn_count_to_DNp04 = length(Posterior_LPLC2_DNp04$updated_ids) # 49 synapses

# LC4 > DNp06
points3d(Dorsal_LC4_DNp06$post_x, Dorsal_LC4_DNp06$post_y, Dorsal_LC4_DNp06$post_z, col = 'red', size= 10)
points3d(Ventral_LC4_DNp06$post_x, Ventral_LC4_DNp06$post_y, Ventral_LC4_DNp06$post_z, col = 'blue', size= 10)
points3d(Anterior_LC4_DNp06$post_x, Anterior_LC4_DNp06$post_y, Anterior_LC4_DNp06$post_z, col = 'purple', size= 10)
points3d(Posterior_LC4_DNp06$post_x, Posterior_LC4_DNp06$post_y, Posterior_LC4_DNp06$post_z, col = 'cyan', size= 10)

Dorsal_LC4_syn_count_to_DNp06 = length(Dorsal_LC4_DNp06$updated_ids) # 11 synapses
Ventral_LC4_syn_count_to_DNp06 = length(Ventral_LC4_DNp06$updated_ids) # 23 synapses
Anterior_LC4_syn_count_to_DNp06 = length(Anterior_LC4_DNp06$updated_ids) # 19 synapses
Posterior_LC4_syn_count_to_DNp06 = length(Posterior_LC4_DNp06$updated_ids) # 24 synapses

# LC6 > DNp06
points3d(Dorsal_LC6_DNp06$post_x, Dorsal_LC6_DNp06$post_y, Dorsal_LC6_DNp06$post_z, col = 'red', size= 10)
points3d(Ventral_LC6_DNp06$post_x, Ventral_LC6_DNp06$post_y, Ventral_LC6_DNp06$post_z, col = 'blue', size= 10)
points3d(Anterior_LC6_DNp06$post_x, Anterior_LC6_DNp06$post_y, Anterior_LC6_DNp06$post_z, col = 'purple', size= 10)
points3d(Posterior_LC6_DNp06$post_x, Posterior_LC6_DNp06$post_y, Posterior_LC6_DNp06$post_z, col = 'cyan', size= 10)

Dorsal_LC6_syn_count_to_DNp06 = length(Dorsal_LC6_DNp06$updated_ids) # 8 synapses
Ventral_LC6_syn_count_to_DNp06 = length(Ventral_LC6_DNp06$updated_ids) # 8 synapses
Anterior_LC6_syn_count_to_DNp06 = length(Anterior_LC6_DNp06$updated_ids) # 1 synapses
Posterior_LC6_syn_count_to_DNp06 = length(Posterior_LC6_DNp06$updated_ids) # 14 synapses

# LPLC1 > DNp06
points3d(Dorsal_LPLC1_DNp06$post_x, Dorsal_LPLC1_DNp06$post_y, Dorsal_LPLC1_DNp06$post_z, col = 'red', size= 10)
points3d(Ventral_LPLC1_DNp06$post_x, Ventral_LPLC1_DNp06$post_y, Ventral_LPLC1_DNp06$post_z, col = 'blue', size= 10)
points3d(Anterior_LPLC1_DNp06$post_x, Anterior_LPLC1_DNp06$post_y, Anterior_LPLC1_DNp06$post_z, col = 'purple', size= 10)
points3d(Posterior_LPLC1_DNp06$post_x, Posterior_LPLC1_DNp06$post_y, Posterior_LPLC1_DNp06$post_z, col = 'cyan', size= 10)

Dorsal_LPLC1_syn_count_to_DNp06 = length(Dorsal_LPLC1_DNp06$updated_ids) # 47 synapses
Ventral_LPLC1_syn_count_to_DNp06 = length(Ventral_LPLC1_DNp06$updated_ids) # 215 synapses
Anterior_LPLC1_syn_count_to_DNp06 = length(Anterior_LPLC1_DNp06$updated_ids) # 50 synapses
Posterior_LPLC1_syn_count_to_DNp06 = length(Posterior_LPLC1_DNp06$updated_ids) # 128 synapses

# LPLC2 > DNp06
points3d(Dorsal_LPLC2_DNp06$post_x, Dorsal_LPLC2_DNp06$post_y, Dorsal_LPLC2_DNp06$post_z, col = 'red', size= 10)
points3d(Ventral_LPLC2_DNp06$post_x, Ventral_LPLC2_DNp06$post_y, Ventral_LPLC2_DNp06$post_z, col = 'blue', size= 10)
points3d(Anterior_LPLC2_DNp06$post_x, Anterior_LPLC2_DNp06$post_y, Anterior_LPLC2_DNp06$post_z, col = 'purple', size= 10)
points3d(Posterior_LPLC2_DNp06$post_x, Posterior_LPLC2_DNp06$post_y, Posterior_LPLC2_DNp06$post_z, col = 'cyan', size= 10)

Dorsal_LPLC2_syn_count_to_DNp06 = length(Dorsal_LPLC2_DNp06$updated_ids) # 72 synapses
Ventral_LPLC2_syn_count_to_DNp06 = length(Ventral_LPLC2_DNp06$updated_ids) # 130 synapses
Anterior_LPLC2_syn_count_to_DNp06 = length(Anterior_LPLC2_DNp06$updated_ids) # 71 synapses
Posterior_LPLC2_syn_count_to_DNp06 = length(Posterior_LPLC2_DNp06$updated_ids) # 59 synapses


#(Figure 2B-F Left)

#VPNs to DNs
VPNs_to_DNp01 <- rbind(LC4_receptive_field, LPLC2_receptive_field)
VPNs_to_DNp02 <- rbind(LC4_receptive_field)
VPNs_to_DNp03 <- rbind(LC4_receptive_field, LPLC1_receptive_field, LPLC2_receptive_field, LPLC4_receptive_field, LC22_receptive_field)
VPNs_to_DNp04 <- rbind(LC4_receptive_field, LPLC1_receptive_field, LPLC2_receptive_field)
VPNs_to_DNp06 <- rbind(LC4_receptive_field, LC6_receptive_field, LPLC1_receptive_field, LPLC2_receptive_field)
VPNs_to_DN <- rbind(LC4_receptive_field, LC6_receptive_field, LC22_receptive_field, LPLC1_receptive_field, LPLC2_receptive_field, LPLC4_receptive_field)

#Color coding individual synapses for plotting
VPNs = c("LC4", "LC6", "LC22", "LPLC1", "LPLC2", "LPLC4")

DNp01_syn_fly_points= DNp01_syn_fly %>% mutate(VPN_color = case_when(type == "LC4" ~ "blue", type == "LPLC1" ~ "red", type == "LPLC4" ~ "green", type == "LC22" ~ "yellow",type == "LPLC2" ~ "orange", type == "LC6" ~ "maroon1"))
DNp01_syn_fly_points <- DNp01_syn_fly_points[DNp01_syn_fly_points$type %in% VPNs,]

DNp02_syn_fly_points= DNp02_syn_fly %>% mutate(VPN_color = case_when(type == "LC4" ~ "blue", type == "LPLC1" ~ "red", type == "LPLC4" ~ "green", type == "LC22" ~ "yellow",type == "LPLC2" ~ "orange", type == "LC6" ~ "maroon1"))
DNp02_syn_fly_points <- DNp02_syn_fly_points[DNp02_syn_fly_points$type %in% VPNs,]

DNp03_syn_fly_points= DNp03_syn_fly %>% mutate(VPN_color = case_when(type == "LC4" ~ "blue", type == "LPLC1" ~ "red", type == "LPLC4" ~ "green", type == "LC22" ~ "yellow",type == "LPLC2" ~ "orange", type == "LC6" ~ "maroon1"))
DNp03_syn_fly_points <- DNp03_syn_fly_points[DNp03_syn_fly_points$type %in% VPNs,]

DNp04_syn_fly_points= DNp04_syn_fly %>% mutate(VPN_color = case_when(type == "LC4" ~ "blue", type == "LPLC1" ~ "red", type == "LPLC4" ~ "green", type == "LC22" ~ "yellow",type == "LPLC2" ~ "orange", type == "LC6" ~ "maroon1"))
DNp04_syn_fly_points <- DNp04_syn_fly_points[DNp04_syn_fly_points$type %in% VPNs,]

DNp06_syn_fly_points= DNp06_syn_fly %>% mutate(VPN_color = case_when(type == "LC4" ~ "blue", type == "LPLC1" ~ "red", type == "LPLC4" ~ "green", type == "LC22" ~ "yellow",type == "LPLC2" ~ "orange", type == "LC6" ~ "maroon1"))
DNp06_syn_fly_points <- DNp06_syn_fly_points[DNp06_syn_fly_points$type %in% VPNs,]

#Plotting of VPN synapses with mesh data
#Synapse location data has been converted to ums, hence dividing the mesh from nm units to um.
plot3d(DNp01_mesh/1000, col='black')
points3d(DNp01_syn_fly_points$post_x, DNp01_syn_fly_points$post_y, DNp01_syn_fly_points$post_z, col = DNp01_syn_fly_points$VPN_color, size= 6)

plot3d(DNp02_mesh/1000, col='black')
points3d(DNp02_syn_fly_points$post_x, DNp02_syn_fly_points$post_y, DNp02_syn_fly_points$post_z, col = DNp02_syn_fly_points$VPN_color, size= 6)

plot3d(DNp03_mesh/1000, col = 'black')
points3d(DNp03_syn_fly_points$post_x, DNp03_syn_fly_points$post_y, DNp03_syn_fly_points$post_z, col = DNp03_syn_fly_points$VPN_color, size= 6)

plot3d(DNp04_mesh/1000, col = 'black')
points3d(DNp04_syn_fly_points$post_x, DNp04_syn_fly_points$post_y, DNp04_syn_fly_points$post_z, col = DNp04_syn_fly_points$VPN_color, size= 6)

plot3d(DNp06_mesh/1000, col = 'black')
points3d(DNp06_syn_fly_points$post_x, DNp06_syn_fly_points$post_y, DNp06_syn_fly_points$post_z, col = DNp06_syn_fly_points$VPN_color, size= 6)
