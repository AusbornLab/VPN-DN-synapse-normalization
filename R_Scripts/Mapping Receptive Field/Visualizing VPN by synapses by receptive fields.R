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

LC4_receptive_field <- read_excel("datafiles/Receptive_field_files/LC4_receptive_field.xlsx")
LC6_receptive_field <- read_excel("datafiles/Receptive_field_files/LC6_receptive_field.xlsx")
LC22_receptive_field <- read_excel("datafiles/Receptive_field_files/LC22_receptive_field.xlsx")
LPLC1_receptive_field <- read_excel("datafiles/Receptive_field_files/LPLC1_receptive_field.xlsx")
LPLC2_receptive_field <- read_excel("datafiles/Receptive_field_files/LPLC2_receptive_field.xlsx")
LPLC4_receptive_field <- read_excel("datafiles/Receptive_field_files/LPLC4_receptive_field.xlsx")

# Loading in the synapse data
DNp01_syn_fly <- read_excel("datafiles/morphologyData/DNp01_morphData/DNp01_caveclient_syn_labeled_um_8_19_2025.xlsx")
DNp02_syn_fly <- read_excel("datafiles/morphologyData/DNp02_morphData/DNp02_caveclient_syn_labeled_um_8_19_2025.xlsx")
DNp03_syn_fly <- read_excel("datafiles/morphologyData/DNp03_morphData/DNp03_caveclient_syn_labeled_um_8_19_2025.xlsx")
DNp04_syn_fly <- read_excel("datafiles/morphologyData/DNp04_morphData/DNp04_caveclient_syn_labeled_um_8_19_2025.xlsx")
DNp06_syn_fly <- read_excel("datafiles/morphologyData/DNp06_morphData/DNp06_caveclient_syn_labeled_um_8_19_2025.xlsx")

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


### Ploting the synapses to each DN by their axis ###

## Subsetting the VPN synapses for DN of interest the IDs from the receptive fields.
#DNp01 VPNS: LC4 and LPLC2
Dorsal_LC4_DNp01 <- DNp01_syn_fly_points[DNp01_syn_fly_points$pre %in% Dorsal_LC4$updated_ids, ]
Ventral_LC4_DNp01 <- DNp01_syn_fly_points[DNp01_syn_fly_points$pre %in% Ventral_LC4$updated_ids, ]
Anterior_LC4_DNp01 <- DNp01_syn_fly_points[DNp01_syn_fly_points$pre %in% Anterior_LC4$updated_ids, ]
Posterior_LC4_DNp01 <- DNp01_syn_fly_points[DNp01_syn_fly_points$pre %in% Posterior_LC4$updated_ids, ]

Dorsal_LPLC2_DNp01 <- DNp01_syn_fly_points[DNp01_syn_fly_points$pre %in% Dorsal_LPLC2$updated_ids, ]
Ventral_LPLC2_DNp01 <- DNp01_syn_fly_points[DNp01_syn_fly_points$pre %in% Ventral_LPLC2$updated_ids, ]
Anterior_LPLC2_DNp01 <- DNp01_syn_fly_points[DNp01_syn_fly_points$pre %in% Anterior_LPLC2$updated_ids, ]
Posterior_LPLC2_DNp01 <- DNp01_syn_fly_points[DNp01_syn_fly_points$pre %in% Posterior_LPLC2$updated_ids, ]

#DNp02 VPNS: LC4 
Dorsal_LC4_DNp02 <- DNp02_syn_fly_points[DNp02_syn_fly_points$pre %in% Dorsal_LC4$updated_ids, ]
Ventral_LC4_DNp02 <- DNp02_syn_fly_points[DNp02_syn_fly_points$pre %in% Ventral_LC4$updated_ids, ]
Anterior_LC4_DNp02 <- DNp02_syn_fly_points[DNp02_syn_fly_points$pre %in% Anterior_LC4$updated_ids, ]
Posterior_LC4_DNp02 <- DNp02_syn_fly_points[DNp02_syn_fly_points$pre %in% Posterior_LC4$updated_ids, ]

#DNp03 VPNS: LC4, LC22, LPLC1, and LPLC4
Dorsal_LC4_DNp03 <- DNp03_syn_fly_points[DNp03_syn_fly_points$pre %in% Dorsal_LC4$updated_ids, ]
Ventral_LC4_DNp03 <- DNp03_syn_fly_points[DNp03_syn_fly_points$pre %in% Ventral_LC4$updated_ids, ]
Anterior_LC4_DNp03 <- DNp03_syn_fly_points[DNp03_syn_fly_points$pre %in% Anterior_LC4$updated_ids, ]
Posterior_LC4_DNp03 <- DNp03_syn_fly_points[DNp03_syn_fly_points$pre %in% Posterior_LC4$updated_ids, ]

Dorsal_LC22_DNp03 <- DNp03_syn_fly_points[DNp03_syn_fly_points$pre %in% Dorsal_LC22$updated_ids, ]
Ventral_LC22_DNp03 <- DNp03_syn_fly_points[DNp03_syn_fly_points$pre %in% Ventral_LC22$updated_ids, ]
Anterior_LC22_DNp03 <- DNp03_syn_fly_points[DNp03_syn_fly_points$pre %in% Anterior_LC22$updated_ids, ]
Posterior_LC22_DNp03 <- DNp03_syn_fly_points[DNp03_syn_fly_points$pre %in% Posterior_LC22$updated_ids, ]

Dorsal_LPLC1_DNp03 <- DNp03_syn_fly_points[DNp03_syn_fly_points$pre %in% Dorsal_LPLC1$updated_ids, ]
Ventral_LPLC1_DNp03 <- DNp03_syn_fly_points[DNp03_syn_fly_points$pre %in% Ventral_LPLC1$updated_ids, ]
Anterior_LPLC1_DNp03 <- DNp03_syn_fly_points[DNp03_syn_fly_points$pre %in% Anterior_LPLC1$updated_ids, ]
Posterior_LPLC1_DNp03 <- DNp03_syn_fly_points[DNp03_syn_fly_points$pre %in% Posterior_LPLC1$updated_ids, ]

Dorsal_LPLC4_DNp03 <- DNp03_syn_fly_points[DNp03_syn_fly_points$pre %in% Dorsal_LPLC4$updated_ids, ]
Ventral_LPLC4_DNp03 <- DNp03_syn_fly_points[DNp03_syn_fly_points$pre %in% Ventral_LPLC4$updated_ids, ]
Anterior_LPLC4_DNp03 <- DNp03_syn_fly_points[DNp03_syn_fly_points$pre %in% Anterior_LPLC4$updated_ids, ]
Posterior_LPLC4_DNp03 <- DNp03_syn_fly_points[DNp03_syn_fly_points$pre %in% Posterior_LPLC4$updated_ids, ]

#DNp04 VPNS: LC4, and LPLC2
Dorsal_LC4_DNp04 <- DNp04_syn_fly_points[DNp04_syn_fly_points$pre %in% Dorsal_LC4$updated_ids, ]
Ventral_LC4_DNp04 <- DNp04_syn_fly_points[DNp04_syn_fly_points$pre %in% Ventral_LC4$updated_ids, ]
Anterior_LC4_DNp04 <- DNp04_syn_fly_points[DNp04_syn_fly_points$pre %in% Anterior_LC4$updated_ids, ]
Posterior_LC4_DNp04 <- DNp04_syn_fly_points[DNp04_syn_fly_points$pre %in% Posterior_LC4$updated_ids, ]

Dorsal_LPLC2_DNp04 <- DNp04_syn_fly_points[DNp04_syn_fly_points$pre %in% Dorsal_LPLC2$updated_ids, ]
Ventral_LPLC2_DNp04 <- DNp04_syn_fly_points[DNp04_syn_fly_points$pre %in% Ventral_LPLC2$updated_ids, ]
Anterior_LPLC2_DNp04 <- DNp04_syn_fly_points[DNp04_syn_fly_points$pre %in% Anterior_LPLC2$updated_ids, ]
Posterior_LPLC2_DNp04 <- DNp04_syn_fly_points[DNp04_syn_fly_points$pre %in% Posterior_LPLC2$updated_ids, ]

# DNp06 VPNS: LC4, LC6, LPLC1, and LPLC2
Dorsal_LC4_DNp06 <- DNp06_syn_fly_points[DNp06_syn_fly_points$pre %in% Dorsal_LC4$updated_ids, ]
Ventral_LC4_DNp06 <- DNp06_syn_fly_points[DNp06_syn_fly_points$pre %in% Ventral_LC4$updated_ids, ]
Anterior_LC4_DNp06 <- DNp06_syn_fly_points[DNp06_syn_fly_points$pre %in% Anterior_LC4$updated_ids, ]
Posterior_LC4_DNp06 <- DNp06_syn_fly_points[DNp06_syn_fly_points$pre %in% Posterior_LC4$updated_ids, ]

Dorsal_LC6_DNp06 <- DNp06_syn_fly_points[DNp06_syn_fly_points$pre %in% Dorsal_LC6$updated_ids, ]
Ventral_LC6_DNp06 <- DNp06_syn_fly_points[DNp06_syn_fly_points$pre %in% Ventral_LC6$updated_ids, ]
Anterior_LC6_DNp06 <- DNp06_syn_fly_points[DNp06_syn_fly_points$pre %in% Anterior_LC6$updated_ids, ]
Posterior_LC6_DNp06 <- DNp06_syn_fly_points[DNp06_syn_fly_points$pre %in% Posterior_LC6$updated_ids, ]

Dorsal_LPLC1_DNp06 <- DNp06_syn_fly_points[DNp06_syn_fly_points$pre %in% Dorsal_LPLC1$updated_ids, ]
Ventral_LPLC1_DNp06 <- DNp06_syn_fly_points[DNp06_syn_fly_points$pre %in% Ventral_LPLC1$updated_ids, ]
Anterior_LPLC1_DNp06 <- DNp06_syn_fly_points[DNp06_syn_fly_points$pre %in% Anterior_LPLC1$updated_ids, ]
Posterior_LPLC1_DNp06 <- DNp06_syn_fly_points[DNp06_syn_fly_points$pre %in% Posterior_LPLC1$updated_ids, ]

Dorsal_LPLC2_DNp06 <- DNp06_syn_fly_points[DNp06_syn_fly_points$pre %in% Dorsal_LPLC2$updated_ids, ]
Ventral_LPLC2_DNp06 <- DNp06_syn_fly_points[DNp06_syn_fly_points$pre %in% Ventral_LPLC2$updated_ids, ]
Anterior_LPLC2_DNp06 <- DNp06_syn_fly_points[DNp06_syn_fly_points$pre %in% Anterior_LPLC2$updated_ids, ]
Posterior_LPLC2_DNp06 <- DNp06_syn_fly_points[DNp06_syn_fly_points$pre %in% Posterior_LPLC2$updated_ids, ]

# Counting of synapses across the neurons 
length(Posterior_LC4_DNp01$pre)

Posterior_LPLC2_syns_to_DNp06 = aggregate(data.frame(synapses = Posterior_LPLC2_DNp06$updated_ids), list(LPLC2_ID = Posterior_LPLC2_DNp06$updated_ids), length)
Posterior_LPLC2_syn_count_to_DNp06 = length(Posterior_LPLC2_DNp06$updated_ids) # 59 synapses 

#(Figure 4A-B, but remaining VPN-DNs below as well)
#Plotting of select synapses with a given neuron

# LC4 > DNp01
points3d(Dorsal_LC4_DNp01$post_x, Dorsal_LC4_DNp01$post_y, Dorsal_LC4_DNp01$post_z, col = 'red', size= 16)
points3d(Ventral_LC4_DNp01$post_x, Ventral_LC4_DNp01$post_y, Ventral_LC4_DNp01$post_z, col = 'blue', size= 16)
points3d(Anterior_LC4_DNp01$post_x, Anterior_LC4_DNp01$post_y, Anterior_LC4_DNp01$post_z, col = 'purple', size= 16)
points3d(Posterior_LC4_DNp01$post_x, Posterior_LC4_DNp01$post_y, Posterior_LC4_DNp01$post_z, col = 'cyan', size= 16)
plot3d(DNp01_mesh/1000, col ='black', alpha = 0.3)

# LPLC2 > DNp01
points3d(Dorsal_LPLC2_DNp01$post_x, Dorsal_LPLC2_DNp01$post_y, Dorsal_LPLC2_DNp01$post_z, col = 'red', size= 16)
points3d(Ventral_LPLC2_DNp01$post_x, Ventral_LPLC2_DNp01$post_y, Ventral_LPLC2_DNp01$post_z, col = 'blue', size= 16)
points3d(Anterior_LPLC2_DNp01$post_x, Anterior_LPLC2_DNp01$post_y, Anterior_LPLC2_DNp01$post_z, col = 'purple', size= 16)
points3d(Posterior_LPLC2_DNp01$post_x, Posterior_LPLC2_DNp01$post_y, Posterior_LPLC2_DNp01$post_z, col = 'cyan', size= 16)

#DNp01 VPN zoomed in dendrites 
view3d(fov=0,zoom=0.25,userMatrix=rotationMatrix(50/180*pi,1,0,0) %*% rotationMatrix(30/180*pi,0,0,1) %*% rotationMatrix(-55/180*pi,0,1,0))

# Place scalebar near bottom-left-front corner with adjusted position
scalebar_length <- 15
margin <- 0.8 # x and y margin
depth_fraction <- 0.3  # how far forward along Z (0 = back, 1 = front)

xrange <- diff(bbox[1:2])
yrange <- diff(bbox[3:4])
zrange <- diff(bbox[5:6])

# Starting point of the scalebar
bar_start <- c(
  bbox[1] + margin * xrange,
  bbox[3] + margin * yrange,
  bbox[5] + depth_fraction * zrange
)

# End point along x-axis
bar_end <- bar_start + c(scalebar_length, 0, 0)

# Draw the scalebar line
segments3d(rbind(bar_start, bar_end), col = "black", lwd = 4)

# Add text label centered above the line
text3d(
  x = mean(c(bar_start[1], bar_end[1])),
  y = bar_start[2] + 0.02 * yrange,  # slight vertical offset
  z = bar_start[3],
  texts = paste0(scalebar_length, "um"),
  col = "black",
  cex = 1
)
rgl.snapshot(filename = "LC4_AP_synapses_to_DNp01.png",fmt = "png")

#DNp03 VPN zoomed in dendrites 
# LPLC4 > DNp03
points3d(Dorsal_LPLC4_DNp03$post_x, Dorsal_LPLC4_DNp03$post_y, Dorsal_LPLC4_DNp03$post_z, col = 'red', size= 16)
points3d(Ventral_LPLC4_DNp03$post_x, Ventral_LPLC4_DNp03$post_y, Ventral_LPLC4_DNp03$post_z, col = 'blue', size= 16)
points3d(Anterior_LPLC4_DNp03$post_x, Anterior_LPLC4_DNp03$post_y, Anterior_LPLC4_DNp03$post_z, col = 'purple', size= 16)
points3d(Posterior_LPLC4_DNp03$post_x, Posterior_LPLC4_DNp03$post_y, Posterior_LPLC4_DNp03$post_z, col = 'cyan', size= 16)
plot3d(DNp03_mesh/1000, col ='black', alpha = 0.3)

scalebar_length <- 15
margin <- 0.45 # x and y margin
depth_fraction <- 0.3  # how far forward along Z (0 = back, 1 = front)

xrange <- diff(bbox[1:2])
yrange <- diff(bbox[3:4])
zrange <- diff(bbox[5:6])

# Starting point of the scalebar
bar_start <- c(
  bbox[1] + margin * xrange,
  bbox[3] + margin * yrange,
  bbox[5] + depth_fraction * zrange
)

# End point along x-axis
bar_end <- bar_start + c(scalebar_length, 0, 2)

# Draw the scalebar line
segments3d(rbind(bar_start, bar_end), col = "black", lwd = 4)

# Add text label centered above the line
text3d(
  x = mean(c(bar_start[1], bar_end[1])),
  y = bar_start[2] + 0.02 * yrange,  # slight vertical offset
  z = bar_start[3],
  texts = paste0(scalebar_length, "um"),
  col = "black",
  cex = 1
)
view3d(fov=10,zoom=0.25,userMatrix=rotationMatrix(31/180*pi,1,0,0) %*% rotationMatrix(20/180*pi,0,0,1) %*% rotationMatrix(-18/180*pi,0,1,0))
rgl.snapshot(filename = "LPLC4_DV_synapses_to_DNp03.png",fmt = "png")

###### Below is the synapse count per 10 VPNs per axis and the plotting of VPN synapses to other DNs. 
get_syn_count <- function(region = "Dorsal", vpn = "LC4", dn = "DNp01") {
  # Construct the full variable name
  var_name <- paste0(region, "_", vpn, "_", dn)
  
  # Access the variable and count synapses
  syn_count <- length(get(var_name)$pre)
  
  return(syn_count)
}

# Example usage:
get_syn_count("Dorsal", "LC4", "DNp01")
get_syn_count("Ventral", "LPLC4", "DNp03")



# LC4 > DNp02
points3d(Dorsal_LC4_DNp02$post_x, Dorsal_LC4_DNp02$post_y, Dorsal_LC4_DNp02$post_z, col = 'red', size= 10)
points3d(Ventral_LC4_DNp02$post_x, Ventral_LC4_DNp02$post_y, Ventral_LC4_DNp02$post_z, col = 'blue', size= 10)
points3d(Anterior_LC4_DNp02$post_x, Anterior_LC4_DNp02$post_y, Anterior_LC4_DNp02$post_z, col = 'purple', size= 10)
points3d(Posterior_LC4_DNp02$post_x, Posterior_LC4_DNp02$post_y, Posterior_LC4_DNp02$post_z, col = 'cyan', size= 10)


# LC4 > DNp03
points3d(Dorsal_LC4_DNp03$post_x, Dorsal_LC4_DNp03$post_y, Dorsal_LC4_DNp03$post_z, col = 'red', size= 10)
points3d(Ventral_LC4_DNp03$post_x, Ventral_LC4_DNp03$post_y, Ventral_LC4_DNp03$post_z, col = 'blue', size= 10)
points3d(Anterior_LC4_DNp03$post_x, Anterior_LC4_DNp03$post_y, Anterior_LC4_DNp03$post_z, col = 'purple', size= 10)
points3d(Posterior_LC4_DNp03$post_x, Posterior_LC4_DNp03$post_y, Posterior_LC4_DNp03$post_z, col = 'cyan', size= 10)


# LC22 > DNp03
points3d(Dorsal_LC22_DNp03$post_x, Dorsal_LC22_DNp03$post_y, Dorsal_LC22_DNp03$post_z, col = 'red', size= 10)
points3d(Ventral_LC22_DNp03$post_x, Ventral_LC22_DNp03$post_y, Ventral_LC22_DNp03$post_z, col = 'blue', size= 10)
points3d(Anterior_LC22_DNp03$post_x, Anterior_LC22_DNp03$post_y, Anterior_LC22_DNp03$post_z, col = 'purple', size= 10)
points3d(Posterior_LC22_DNp03$post_x, Posterior_LC22_DNp03$post_y, Posterior_LC22_DNp03$post_z, col = 'cyan', size= 10)

# LPLC1 > DNp03
points3d(Dorsal_LPLC1_DNp03$post_x, Dorsal_LPLC1_DNp03$post_y, Dorsal_LPLC1_DNp03$post_z, col = 'red', size= 10)
points3d(Ventral_LPLC1_DNp03$post_x, Ventral_LPLC1_DNp03$post_y, Ventral_LPLC1_DNp03$post_z, col = 'blue', size= 10)
points3d(Anterior_LPLC1_DNp03$post_x, Anterior_LPLC1_DNp03$post_y, Anterior_LPLC1_DNp03$post_z, col = 'purple', size= 10)
points3d(Posterior_LPLC1_DNp03$post_x, Posterior_LPLC1_DNp03$post_y, Posterior_LPLC1_DNp03$post_z, col = 'cyan', size= 10)

# LC4 > DNp04
points3d(Dorsal_LC4_DNp04$post_x, Dorsal_LC4_DNp04$post_y, Dorsal_LC4_DNp04$post_z, col = 'red', size= 10)
points3d(Ventral_LC4_DNp04$post_x, Ventral_LC4_DNp04$post_y, Ventral_LC4_DNp04$post_z, col = 'blue', size= 10)
points3d(Anterior_LC4_DNp04$post_x, Anterior_LC4_DNp04$post_y, Anterior_LC4_DNp04$post_z, col = 'purple', size= 10)
points3d(Posterior_LC4_DNp04$post_x, Posterior_LC4_DNp04$post_y, Posterior_LC4_DNp04$post_z, col = 'cyan', size= 10)

# LPLC2 > DNp04
points3d(Dorsal_LPLC2_DNp04$post_x, Dorsal_LPLC2_DNp04$post_y, Dorsal_LPLC2_DNp04$post_z, col = 'red', size= 10)
points3d(Ventral_LPLC2_DNp04$post_x, Ventral_LPLC2_DNp04$post_y, Ventral_LPLC2_DNp04$post_z, col = 'blue', size= 10)
points3d(Anterior_LPLC2_DNp04$post_x, Anterior_LPLC2_DNp04$post_y, Anterior_LPLC2_DNp04$post_z, col = 'purple', size= 10)
points3d(Posterior_LPLC2_DNp04$post_x, Posterior_LPLC2_DNp04$post_y, Posterior_LPLC2_DNp04$post_z, col = 'cyan', size= 10)


# LC4 > DNp06
points3d(Dorsal_LC4_DNp06$post_x, Dorsal_LC4_DNp06$post_y, Dorsal_LC4_DNp06$post_z, col = 'red', size= 10)
points3d(Ventral_LC4_DNp06$post_x, Ventral_LC4_DNp06$post_y, Ventral_LC4_DNp06$post_z, col = 'blue', size= 10)
points3d(Anterior_LC4_DNp06$post_x, Anterior_LC4_DNp06$post_y, Anterior_LC4_DNp06$post_z, col = 'purple', size= 10)
points3d(Posterior_LC4_DNp06$post_x, Posterior_LC4_DNp06$post_y, Posterior_LC4_DNp06$post_z, col = 'cyan', size= 10)

# LC6 > DNp06
points3d(Dorsal_LC6_DNp06$post_x, Dorsal_LC6_DNp06$post_y, Dorsal_LC6_DNp06$post_z, col = 'red', size= 10)
points3d(Ventral_LC6_DNp06$post_x, Ventral_LC6_DNp06$post_y, Ventral_LC6_DNp06$post_z, col = 'blue', size= 10)
points3d(Anterior_LC6_DNp06$post_x, Anterior_LC6_DNp06$post_y, Anterior_LC6_DNp06$post_z, col = 'purple', size= 10)
points3d(Posterior_LC6_DNp06$post_x, Posterior_LC6_DNp06$post_y, Posterior_LC6_DNp06$post_z, col = 'cyan', size= 10)

# LPLC1 > DNp06
points3d(Dorsal_LPLC1_DNp06$post_x, Dorsal_LPLC1_DNp06$post_y, Dorsal_LPLC1_DNp06$post_z, col = 'red', size= 10)
points3d(Ventral_LPLC1_DNp06$post_x, Ventral_LPLC1_DNp06$post_y, Ventral_LPLC1_DNp06$post_z, col = 'blue', size= 10)
points3d(Anterior_LPLC1_DNp06$post_x, Anterior_LPLC1_DNp06$post_y, Anterior_LPLC1_DNp06$post_z, col = 'purple', size= 10)
points3d(Posterior_LPLC1_DNp06$post_x, Posterior_LPLC1_DNp06$post_y, Posterior_LPLC1_DNp06$post_z, col = 'cyan', size= 10)

# LPLC2 > DNp06
points3d(Dorsal_LPLC2_DNp06$post_x, Dorsal_LPLC2_DNp06$post_y, Dorsal_LPLC2_DNp06$post_z, col = 'red', size= 10)
points3d(Ventral_LPLC2_DNp06$post_x, Ventral_LPLC2_DNp06$post_y, Ventral_LPLC2_DNp06$post_z, col = 'blue', size= 10)
points3d(Anterior_LPLC2_DNp06$post_x, Anterior_LPLC2_DNp06$post_y, Anterior_LPLC2_DNp06$post_z, col = 'purple', size= 10)
points3d(Posterior_LPLC2_DNp06$post_x, Posterior_LPLC2_DNp06$post_y, Posterior_LPLC2_DNp06$post_z, col = 'cyan', size= 10)



#(Figure 2B-F Left)

#VPNs to DNs
colnames(LC4_receptive_field)[1] <- "catmaid_id"
colnames(LC6_receptive_field)[1] <- "catmaid_id"
colnames(LC22_receptive_field)[1] <- "catmaid_id"
colnames(LPLC1_receptive_field)[1] <- "catmaid_id"
colnames(LPLC2_receptive_field)[1] <- "catmaid_id"
colnames(LPLC4_receptive_field)[1] <- "catmaid_id"

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
points3d(DNp01_syn_fly_points$post_x, DNp01_syn_fly_points$post_y, DNp01_syn_fly_points$post_z, col = DNp01_syn_fly_points$VPN_color, size= 10)

plot3d(DNp02_mesh/1000, col='black')
points3d(DNp02_syn_fly_points$post_x, DNp02_syn_fly_points$post_y, DNp02_syn_fly_points$post_z, col = DNp02_syn_fly_points$VPN_color, size= 10)

DNp03_skel <-read.neuron("datafiles/morphologyData/DNp03_morphData/DNp03_um_model.swc")

plot3d(DNp03_mesh/1000, col = 'black')
plot3d(DNp03_skel, col = 'black', WithNodes = FALSE, lwd = 3)
points3d(DNp03_syn_fly_points$post_x, DNp03_syn_fly_points$post_y, DNp03_syn_fly_points$post_z, col = DNp03_syn_fly_points$VPN_color, size= 7)

plot3d(DNp04_mesh/1000, col = 'black')
points3d(DNp04_syn_fly_points$post_x, DNp04_syn_fly_points$post_y, DNp04_syn_fly_points$post_z, col = DNp04_syn_fly_points$VPN_color, size= 10)

plot3d(DNp06_mesh/1000, col = 'black')
points3d(DNp06_syn_fly_points$post_x, DNp06_syn_fly_points$post_y, DNp06_syn_fly_points$post_z, col = DNp06_syn_fly_points$VPN_color, size= 10)



# Place scalebar near bottom-left-front corner with adjusted position
scalebar_length <- 50
margin <- 0.4  # x and y margin
depth_fraction <- 1.2  # how far forward along Z (0 = back, 1 = front)

xrange <- diff(bbox[1:2])
yrange <- diff(bbox[3:4])
zrange <- diff(bbox[5:6])

# Starting point of the scalebar
bar_start <- c(
  bbox[1] + margin * xrange,
  bbox[3] + margin * yrange,
  bbox[5] + depth_fraction * zrange
)

# End point along x-axis
bar_end <- bar_start + c(scalebar_length, 0, 0)

# Draw the scalebar line
segments3d(rbind(bar_start, bar_end), col = "black", lwd = 4)

# Add text label centered above the line
text3d(
  x = mean(c(bar_start[1], bar_end[1])),
  y = bar_start[2] + 0.02 * yrange,  # slight vertical offset
  z = bar_start[3],
  texts = paste0(scalebar_length, "um"),
  col = "black",
  cex = 1
)



#Saving fig image. 
#Adjust view point and then save as png
rgl.snapshot(filename = "LPLC2_quadratic_plane_mesh.png",fmt = "png")
