#Note to get the script to run will also need to load in the LC_startup code, which I have to also edit and update to contain only whats necessary

#Script is specific to figure 3.

#Dependencies 
{
library(tidyverse)
library(natverse)
library(alphashape3d) #alpha shape
library(alphahull) #2D alpha hull
library(RColorBrewer)
library(sf) #compute line-polygon intersections
library(neuprintr)
library(leaflet)
library(fields)
library(writexl)
library(stringi)
library(sp)
library(fafbseg)
}

#Select depending on the population being used
#LC4
a = -0.6; b = 1; c = 1.3; d = -265000 
#LC6
a = -0.6; b = 1; c = 1.32; d = -265000

#LC22
a = -0.62; b = 1.05; c = 1.3; d = -265000 

#LPLC1
a = -0.6; b = 1; c = 1.3; d = -265000 

#LPLC2
a = -0.6; b = 1; c = 1.3; d = -242500 

#LPLC4
a = -0.6; b = 1; c = 1.3; d = -265000 

#If working with an VPN without dendrites in the lobula plate ONLY (Lc4, LC6 etc)

{neu <- LC4

for (j in 1:length(neu)) {
  tar <- neu[[j]]
  xyz_ep <-  tar$d[tar$EndPoints, ] %>% xyzmatrix()
  xyz_bp = tar$d[tar$BranchPoints, ] %>% xyzmatrix()
  xyz_LO <- rbind(xyz_ep, xyz_bp) %>% 
    as_tibble() %>% 
    mutate(LO = a*X + b*Y + c*Z +d) %>%
    filter(LO > 0) %>%
    select(X,Y,Z)
  if (j == 1) {
    xyz_node <- xyz_LO
  } else {
    xyz_node <- bind_rows(xyz_node, xyz_LO)
  }
}
xyz_node <- data.matrix(xyz_node)

#Make alpha mesh
msh.a <- ashape3d(xyz_node, alpha = 20000, pert = T) # 20000 look ok
#plot3d(msh.a)
msh <- as.mesh3d(msh.a)
#plot3d(msh.a)
neu_lo <- nlapply(neu, subset, function(x) pointsinside(x, msh))
}


#If you are using a VPN population with dendrites in the Lobula plate
# need to use LO_msh for LPLCs and can use msh for just LCs, this below part stays commented out if only working with VPNs with no LP dendrites.
#If working with LPLC4, use the quadratic plane for LPLC1. 

#neu <- LPLC2
#neu_lo <- LPLC2

#For populations with lp dendrites will use a LO_msh, to select only lobula dendrites
#neu_lo <- nlapply(neu_lo, subset, function(x) pointsinside(x, LO_msh))
#plot3d(neu_lo)

#This below will select the points of the LPLC1 population with the above set transformations to remove parts of the tethers/axons and limit it to the lobula.
#neu_lo <- LPLC1

#From here you will select the points to calculate the quadratic plane. 
#for (j in 1:length(neu_lo)) {
#  tar <- neu_lo[[j]]
#  xyz_ep <-  tar$d[tar$EndPoints, ] %>% xyzmatrix()
#  xyz_bp = tar$d[tar$BranchPoints, ] %>% xyzmatrix()
#  xyz_LO <- rbind(xyz_ep, xyz_bp) %>% 
#    as_tibble() %>% 
#    mutate(LO = a*X + b*Y + c*Z + d) %>%
#    filter(LO > 0) %>%
#    select(X,Y,Z)
#  if (j == 1) {
#    xyz_node <- xyz_LO
#  } else {
#    xyz_node <- bind_rows(xyz_node, xyz_LO)
#  }
#}

#xyz_node <- data.matrix(xyz_node)

##Visualization of individual points
#{plot3d(xyz_node)
#plot3d(neu_lo)}


# - fit a 2nd order polynomial surface
# this takes some time!
polyfitorder <- 2 # 2nd order surface
gridMargin <- 30000 #add margin on the edge
xyz_node <- data.matrix(xyz_node)
X <- xyz_node[, 1];  Y <- xyz_node[, 2]; Z <- xyz_node[, 3]

# make an underlying grid to interpret the fit as a point set
dx2 <- 500 # grid spacing
dy2 <- 500
xx <- seq(range(X)[1] - gridMargin, range(X)[2] + gridMargin, by = dx2)
yy <- seq(range(Y)[1] - gridMargin, range(Y)[2] + gridMargin, by = dy2)
fitlm <- lm(Z ~ poly(X, Y, degree = polyfitorder, raw = TRUE)) #linear model fit
xygrid <- expand.grid(xx, yy)
xygrid <- setNames(data.frame(xygrid), c('X', 'Y'));
valfit <- predict(fitlm, xygrid) #generate values from the fit
xyz_lm <- cbind(xygrid, valfit)
#dist_min <- apply(xyz_lm, 1, function(pt) {min(rowSums(sweep(xyz_node, 2, pt) ^ 2))}) 
# keep points inside the lobula mesh
ii <- pointsinside(xyz_lm, LO_msh,rval='distance') > -3000
xyz_layer <- xyz_lm[ii,] # pts
#plot3d(neu_lo)


#If working with LPLC4 then run below after calculating the LPLC1 quadratic plane, otherwise skip this.
#neu_lo <- LPLC4

#For populations with lp dendrites will use a LO_msh, to select only lobula dendrites
#neu_lo <- nlapply(neu_lo, subset, function(x) pointsinside(x, LO_msh))


#(Figure 3 left panels)
### Visualizing the plane through the Lobula dendrites with skeletons
nopen3d()
par3d('windowRect' = c(10,10,1710,1710))
points3d(xyz_layer, color = "orange", alpha = 0.9, size = 2)
# points3d(dend_v, size = 2, col = 'gray90')
plot3d(neu[], col='gray80', soma = T, lwd=1, WithNodes = F)
plot3d(neu_lo[], col='gray80', soma = T, lwd=1, WithNodes = F)
plot3d(neu[[2]],  col= "red", lwd = 5, soma=T, WithNodes = F)
plot3d(neu[[8]],  col= "blue", lwd = 5, soma=T, WithNodes = F)
plot3d(nlapply(TM5[1], subset, function(x) pointsinside(x, msh,rval='distance')>-0.6e4), col = 'magenta', lwd = 4) #TM5
plot3d(nlapply(TM5[2], subset, function(x) pointsinside(x, msh,rval='distance')>-1.2e4), col = 'magenta', lwd = 4) 
arrow3d(axis_ori, axis_ori + axis_lat, theta = pi/6,n = 4, col="green", type = "rotation")
arrow3d(axis_ori, axis_ori + axis_dor, theta = pi/6,n = 4, col="magenta", type = "rotation")
arrow3d(axis_ori, axis_ori + axis_post, theta = pi/6,n = 4, col="blue", type = "rotation")
rgl.viewpoint(fov=0,zoom=0.8,userMatrix=rotationMatrix(170/180*pi,1,0,0) %*% rotationMatrix(30/180*pi,0,0,1) %*% rotationMatrix(-65/180*pi,0,1,0))
# plot3d(neu_target, lwd=5)

#Plotting of mesh and quadratic plane, same as above.
#change the ids below loaded from the previous file to read and plot meshes
VPN_ids = flywire_updateids(LC4_ids_mesh)
VPN_mesh=read_cloudvolume_meshes(VPN_ids)
# LC4_red = 23, blue = 4
# LC6_red = 61, blue = 40
# LC22_red = 26, blue = 18
# LPLC1_red = 23, blue = 47
# LPLC2_red = 1, blue = 8
# LPLC4_red = 46, blue = 24

plot3d(VPN_mesh[],  col= "red", lwd = 5, soma=T, WithNodes = F) 
plot3d(VPN_mesh[],  col= "blue", lwd = 5, soma=T, WithNodes = F)
plot3d(VPN_mesh[], col='gray80', soma = T, lwd=1, WithNodes = F)
points3d(xyz_layer, color = "orange", alpha = 0.9, size = 2)
plot3d(nlapply(TM5[1], subset, function(x) pointsinside(x, msh,rval='distance')>-0.6e4), col = 'magenta', lwd = 4) #TM5
plot3d(nlapply(TM5[2], subset, function(x) pointsinside(x, msh,rval='distance')>-1.2e4), col = 'magenta', lwd = 4) 
rgl.viewpoint(fov=0,zoom=0.8,userMatrix=rotationMatrix(170/180*pi,1,0,0) %*% rotationMatrix(30/180*pi,0,0,1) %*% rotationMatrix(-65/180*pi,0,1,0))

#Saving fig image. 
#rgl.snapshot(filename = "LPLC2_quadratic_plane_mesh.png",fmt = "png")



# this takes some time
row.names(xyz_layer) <-  seq(1, dim(xyz_layer)[1])
ind_pj <- list() #index of projected grid points
xyz_com <- list() # center-of-mass
xyz_pj_com <- list() # xyz of com projecting on grid
for (j in 1:length(neu_lo)){
  tar <- neu_lo[[j]]
  xyz_ep <-  tar$d[tar$EndPoints, ] %>% 
    xyzmatrix() %>%
    as_tibble() %>% 
    mutate(LO = a*X + b*Y + c*Z + d) %>%
    filter(LO > 0) %>%
    select(X,Y,Z) %>%
    data.matrix()
  xyz_bp <-  tar$d[tar$BranchPoints, ] %>% 
    xyzmatrix() %>%
    as_tibble() %>% 
    mutate(LO = a*X + b*Y + c*Z + d) %>%
    filter(LO > 0) %>%
    select(X,Y,Z) %>%
    data.matrix()
  
  # center-of-mass
  xyz_eb <- rbind(xyz_ep,xyz_bp)
  xyz_com[[j]] <- colSums(xyz_eb)/dim(xyz_eb)[1]
  
  # project dendrite end points and com to the fitted grid by shortest distance
  xyz_dend_pj <- rbind(xyz_com[[j]], xyz_ep) #append com at the beginning
  Nsigma <- 5 #exclude > 5 sigma
  com_dist <- as.matrix(dist(xyz_dend_pj))
  thrhd_dist <- sd(com_dist[,1])*Nsigma
  xyz_dend_pj <- cbind(xyz_dend_pj, com_dist[,1]) %>%
    as_tibble() %>%
    filter(V4 < thrhd_dist)%>%
    select(X,Y,Z) %>%
    data.matrix()
  
  ind_min <- c()
  ind_min <-apply(xyz_dend_pj, 1, function(pt) {which.min(rowSums(sweep(xyz_layer, 2, pt) ^ 2))}) # index in xyz_layer with min distance
  ind_com <- ind_min[1] #index of grid point that's closest to com
  ind_min <- unique(ind_min[-1]) #index of grid points
  xyz_pj_com[[j]] <- xyz_layer[ind_com,]
  ind_pj[[j]] <- row.names(xyz_layer[ind_min,])
}

# project onto a plane def by pc3 of xyz_layer
layer_pca <- prcomp(xyz_layer)
xyz_layer_rot <- t(layer_pca$rotation) %*% t(sweep(xyz_layer, 2, layer_pca$center,))
xyz_layer_rot <- t(xyz_layer_rot)
xy_layer_rot <- xyz_layer_rot[,c(1,2)]


# TM5 positions
TM5_c_xyz <- xyzmatrix(TM5[[1]]$d[match(TM5[[1]]$tags$"TM5 LO col", TM5[[1]]$d$PointNo),])
TM5_u_xyz <- xyzmatrix(TM5[[2]]$d[match(TM5[[2]]$tags$"TM5 LO col", TM5[[2]]$d$PointNo),])

# grid points near 2 TM5s
grid_c <- which.min(rowSums((sweep(xyz_layer, 2, c(TM5_c_xyz), "-"))^2))
grid_u <- which.min(rowSums((sweep(xyz_layer, 2, c(TM5_u_xyz), "-"))^2))


# - map equator and coord system 
center_new <- xyz_layer_rot[grid_c,1:2]
x_med_new <- as.numeric(center_new[1])
y_eq_new <- as.numeric(center_new[2])
angR_new <- acos((x_med_new - xyz_layer_rot[grid_u,1])/sqrt(sum((xyz_layer_rot[grid_c,1:2] - xyz_layer_rot[grid_u,1:2])^2)))

# 2d projection with dendrites, centered by 2 TM5
ang_2 <- pi/2 - angR_new
rot_2 <- matrix(c(cos(ang_2), sin(ang_2), -sin(ang_2), cos(ang_2)), ncol = 2)
xy_layer_align <- sweep(xy_layer_rot, MARGIN = 2, STATS = c(x_med_new, y_eq_new))
xy_layer_align <- t(rot_2 %*% t(xy_layer_align))



# reconstruct (x,y) projection - just keep x,y, then align
xy_pj <- list()
xy_pj_com <- list()
for (j in 1:length(ind_pj)){
  xy_tmp <- list()
  xy_ashape <- list()
  ii <- sort(as.integer(ind_pj[[j]]))
  xy_tmp <- xy_layer_align[ii, ]
  xy_ashape <- ashape(xy_tmp, alpha = 6000)
  xy_edge <- xy_ashape$edges[,1:6]
  xy_pj[[j]] <- list(xy=xy_tmp, ashape=xy_ashape, edge=xy_edge)
  
  xy_com_tmp <- unlist(xyz_pj_com[[j]]-layer_pca$center) %*% layer_pca$rotation
  xy_com_tmp <- xy_com_tmp[c(1,2)] - c(x_med_new, y_eq_new)
  xy_com_tmp <- t(rot_2 %*% xy_com_tmp)
  xy_pj_com[[j]] <- xy_com_tmp
}


### Plotting of the receptive field (Figure 3 middle) ###
windows(record = F, width = 8, height = 8)
# pdf(file = "LC4_2d.pdf", width = 8, height = 8,pointsize=12,family="Helvetica", useDingbats = F)
plot(xy_layer_align, col="orange", cex = 1, pch = ".",
     ylim = rev(range(xy_layer_align[,2])),
     xlim = (range(xy_layer_align[,1])),
     asp = 1,
     main = title("LC4 Receptive Field", cex=1.5, line=1),
     xlab = "Anterior-Posterior Coordinates in µm ",
     ylab = "Dorsal-Ventral Coordinates in µm")
for (j in 1:length(neu_lo)) {
twig <- neu_lo[[j]]
pp <- as.matrix(sweep(twig$d[,c("X","Y","Z")], 2, layer_pca$center)) %*% layer_pca$rotation
pp <- sweep(pp[,1:2], 2, STATS = c(x_med_new, y_eq_new))
twig$d[,c("X","Y")] <- sweep(t(rot_2 %*% t(pp)), 2, c(x_med_new, y_eq_new))
#plot(twig/1000,  col= col_com_2$colors[j], lwd = 2, soma=T, WithNodes = F, add = T)
plot(twig,  col= "gray80", lwd = 2, soma=F, WithNodes = F, add = T)
}
twig <- neu_lo[[30]]
pp <- as.matrix(sweep(twig$d[,c("X","Y","Z")], 2, layer_pca$center)) %*% layer_pca$rotation
pp <- sweep(pp[,1:2], 2, STATS = c(x_med_new, y_eq_new))
twig$d[,c("X","Y")] <- sweep(t(rot_2 %*% t(pp)), 2, c(x_med_new, y_eq_new))
plot(twig,  col= "red", lwd = 2, soma=T, WithNodes = F, add = T)
twig <- neu_lo[[4]]
pp <- as.matrix(sweep(twig$d[,c("X","Y","Z")], 2, layer_pca$center)) %*% layer_pca$rotation
pp <- sweep(pp[,1:2], 2, STATS = c(x_med_new, y_eq_new))
twig$d[,c("X","Y")] <- sweep(t(rot_2 %*% t(pp)), 2, c(x_med_new, y_eq_new))
plot(twig,  col= "blue", lwd = 2, soma=F, WithNodes = F, add = T)
lines(rbind(-xy_layer_align[grid_u,]*2.18, xy_layer_align[grid_u,]*1.8), lwd = 3, col = 'cyan')
points(matrix(c(x_med_new,y_eq_new), ncol =2), pch = 18, col = 'magenta', cex = 2)
points(matrix(xy_layer_align[grid_u,], ncol=2), pch = 18, col = 'magenta', cex = 2)
lines(c(50000,60000), c(-70000, -70000), col = 'black', lwd = 3)
text(x = 55000, -65000, labels = "10 µm")
points(matrix(unlist(xy_pj_com), ncol = 2, byrow = T), pch = 20, col = "black", cex = 1.5) 

## Setting up for transformation from lobula space to visual field
xy_poly <- list()
xy_bdpt <- matrix(ncol = 2, nrow = 0) #collect all LC4 polygon edge points to make a LO boundary
for (j in 1:length(ind_pj)) {
  ls_poly <- mkpoly(xy_pj[[j]]$ashape$edges)
  xy_poly[[j]] <- ls_poly[[1]][,3:4]
  xy_bdpt <- rbind(xy_bdpt, xy_pj[[j]]$edge[,3:4])
}


# use alpha-hull for the grid projection, 
xy_ashape_grid <- ashape(unique(xy_bdpt), alpha = 7500)
xy_grid_ahull <- mkpoly(xy_ashape_grid$edges)[[1]][,3:4]
xy_edge_grid <- xy_grid_ahull # hull edge points

grid_bdpt <- xy_edge_grid
grid_bdpt <- rbind(grid_bdpt, grid_bdpt[1,])

ymax <- max(xy_edge_grid[,2])
ymin <- min(xy_edge_grid[,2])

# VPN polygons and centers, in eye coordinates
poly_st <- st_polygon(list(data.matrix(grid_bdpt)))
xy_ori <- c(0,0)
R <- (ymax-ymin)*2
xy_com <- list() # com of the projected LC4 dendrite in eye coord

for (j in 1:length(xy_poly)) {
  pj_com <- xy_pj_com[[j]] %>% as.matrix()
  colnames(pj_com) <- c("x1", "y1")
  xy_poly[[j]] <- rbind(pj_com, xy_poly[[j]])
  xy_poly[[j]] %<>% 
    as_tibble() %>%
    mutate(phiC = acos((x1-x_med_new)/sqrt((x1-x_med_new)^2+(y1-y_eq_new)^2))*(-1)^(y1<y_eq_new) + 2*pi*(y1<y_eq_new)) %>%  #angle
    mutate(thetaC = NA) %>% #radius
    transmute(x1, y1, phiC, thetaC) %>%
    data.matrix()
  for (k in 1:dim(xy_poly[[j]])[1]) {
    alpha <- xy_poly[[j]][k,'phiC']
    line_st <- st_linestring(rbind(xy_ori, c(xy_ori[1] + R*cos(alpha), xy_ori[2] + R*sin(alpha))))
    int <- data.matrix(st_intersection(line_st, poly_st))[2,]
    xy_poly[[j]][k,'thetaC'] <- pi/2 * dist(rbind(xy_poly[[j]][k,1:2], xy_ori)) / dist(rbind(int,xy_ori))
  }
  # now turn inside out and a 90-rotation about x-axis to have the front edge on the x-z plane, 
  # angle wrt looking from behind the eye
  xy_poly[[j]] %<>%
    as_tibble() %>%
    # mutate(thetaC = pi - thetaC) %>% # turn inside out
    mutate(phiC = phiC ) %>% # align front edge to x-z plane
    mutate(x = sin(thetaC)*cos(phiC), y = sin(thetaC)*sin(phiC), z = cos(thetaC)) %>% #(x,y,z) coord
    mutate(xr = x, yr = z, zr = -y) %>% # +90 rotation about x-axis
    mutate(xrr = xr, yrr = yr, zrr = zr) %>% 
    mutate(theta = acos(zrr/sqrt(xrr^2+yrr^2+zrr^2)), phi = acos(xrr/sqrt(xrr^2+yrr^2))) %>%
    # mutate(theta_deg = 180 - theta/pi*180, phi_deg = phi/pi*180) %>% #the convention for x-y plot below
    mutate(theta_deg = theta/pi*180/180*diff(buchner_theta)+buchner_theta[1], 
           phi_deg = phi/pi*180/180*diff(buchner_phi)+buchner_phi[1]) %>% # longitude use buchner_phi
    # select(x1, y1, theta_deg, phi_deg) %>%
    data.matrix()
  
  xyMollweide <- Mollweide(xy_poly[[j]][,c("theta_deg", "phi_deg")])
  colnames(xyMollweide) <- c('xM','yM')
  xy_poly[[j]] <- cbind(xy_poly[[j]], xyMollweide)
  
  xy_com[[j]] <- xy_poly[[j]][1,]
  xy_poly[[j]] <- xy_poly[[j]][-1,]
}



# boundary in eye coord
grid_bdpt_tp <- grid_bdpt %>%
  as_tibble() %>%
  mutate(phiC = acos((x1-x_med_new)/sqrt((x1-x_med_new)^2+(y1-y_eq_new)^2))*(-1)^(y1<y_eq_new) + 2*pi*(y1<y_eq_new)) %>%  #angle
  mutate(thetaC = pi/2) %>% #radius
  mutate(x = sin(thetaC)*cos(phiC), y = sin(thetaC)*sin(phiC), z = cos(thetaC)) %>% #(x,y,z) coord
  mutate(xr = x, yr = z, zr = -y) %>% # +90 rotation about x-axis
  mutate(xrr = xr, yrr = yr, zrr = zr) %>% 
  mutate(theta = acos(zrr/sqrt(xrr^2+yrr^2+zrr^2)), phi = acos(xrr/sqrt(xrr^2+yrr^2))) %>%
  mutate(theta_deg = theta/pi*180/180*diff(buchner_theta)+buchner_theta[1], 
         phi_deg = phi/pi*180/180*diff(buchner_phi)+buchner_phi[1]) %>% # longitude use buchner_phi
  select(phi_deg, theta_deg) %>%
  data.matrix()


# sphere grid for computing areas
S2grid <- matrix(c(0,0), ncol = 2) # [theta phi] in rad
darc <- pi/180*1
for (ti in 1:179) {
  theta <- ti * darc / 1
  ddeg <- darc / sin(theta) / pi *180
  phiN <- floor(180 / ddeg)
  tpmat <- cbind(rep(theta/pi*180, 2*phiN+1), seq(-phiN, phiN, by = 1)*ddeg )
  S2grid <- rbind(S2grid, tpmat)
}
S2grid <- rbind(S2grid, c(180,0)) #south pole
darea <- 4*pi/dim(S2grid)[1]
S2grid <- cbind(S2grid, Mollweide(S2grid))
colnames(S2grid) <- c('t','p','xM','yM')

# background grid
bkgd_grid <- S2grid[, c('xM','yM')] # use S2 grid 


# set up background guidelines, here you can translate the backgrounds,to match the binocular overlap ~15 degress in overlap 
bkgd_eq <- Mollweide(cbind(rep(90,37), seq(-180, 180, by = 10)))
bkgd_eq_p45 <- Mollweide(cbind(rep(45,37), seq(-180, 180, by = 10)))
bkgd_eq_m45 <- Mollweide(cbind(rep(135,37), seq(-180, 180, by = 10)))

bkgd_eq_half <- Mollweide(cbind(rep(90,38), seq(0, 180, by = 10)))
bkgd_eq_p45_half <- Mollweide(cbind(rep(45,38), seq(0, 180, by = 10)))
bkgd_eq_m45_half <- Mollweide(cbind(rep(135,38), seq(0, 180, by = 10)))

bkgd_mer <- Mollweide(cbind(seq(0, 180, by = 10), rep(0,19)))
bkgd_mer_e <- Mollweide(cbind(seq(0, 180, by = 10), rep(90,19)))
bkgd_mer_ee <- Mollweide(cbind(seq(0, 180, by = 10), rep(180,19)))
bkgd_mer_w <- Mollweide(cbind(seq(0, 180, by = 10), rep(-90,19)))
bkgd_mer_ww <- Mollweide(cbind(seq(0, 180, by = 10), rep(-180,19)))


xy_bd_M <- matrix(ncol = 2)
xy_bd <- matrix(ncol = 2)
#bd_grid <- expand.grid(bd_phi, bd_theta) 

windows(record = F, width = 8, height = 8)

#Plotting of COMs of molliewide projection for looming stimulus estimation (Figure 3 right)
{plot(bkgd_grid, cex = 0.6, pch='', xlim = c(-0.5, pi), ylim = c(-1.5,1.5))
for (j in 1:length(xy_poly)) {
  polygon(xy_poly[[j]][,c("xM","yM")], angle = j*2, lwd = 0.5)
}
polygon(xy_poly[[4]][,c("xM","yM")],col = 'blue', density =20, angle = j*2, lwd =1)
polygon(xy_poly[[23]][,c("xM","yM")],col = 'red', density =20, angle = j*2, lwd =1)

lines(bkgd_mer, lwd = 2); lines(bkgd_mer_e, lwd =2, col = 'cyan'); lines(bkgd_mer_ee, lwd = 2);
lines(bkgd_eq_half, lwd =2); lines(bkgd_eq_m45_half, lwd =2 ); lines(bkgd_eq_p45_half, lwd =2)
points(matrix(c(1.12,0.85), ncol =2), pch = 18, col = 'magenta', cex = 2)
points(matrix(c(1.435,0.172), ncol =2), pch = 18, col = 'magenta', cex = 2)

for (j in 1:length(xy_poly)) {
  points(xy_com[[j]]['xM'], xy_com[[j]]['yM'],  col= 'black', cex = 1.5, pch = 20) #pch=1 circle, 32+j ASCII
}
}

######################################
#Calculating the area of each polygon 
plot(bkgd_grid, cex = 0.6, pch='', xlim = c(-0.5, pi), ylim = c(-1.5,1.5))
polygon(xy_poly[[1]][,c("xM","yM")], angle = j*2, lwd = 0.5)

VPN = c()
for (j in 1:length(xy_poly)) {
  polygon_coordinates <- xy_poly[[j]][, c("xM", "yM")]
  polygon <- Polygon(polygon_coordinates)
  
  # Create a Spatial Polygons object with a single polygon
  polygon <- Polygons(list(polygon), ID = "polygon1")
  area = polygon@area
  VPN <- c(VPN, area)
}

#Create a box and whisker plot
boxplot(VPN, main = "LC4 polygon normalized area", ylab = "Normalized Area", outline = TRUE, ylim = c(0, 1))

points(jitter(rep(1, length(VPN)), amount = 0.15), VPN, col = "red", pch = 16)

data <- data.frame(Area = VPN)

# Save the data frame to a CSV file
write.csv(data, "boxplot_data_LC_area.csv", row.names = FALSE)

######################################



