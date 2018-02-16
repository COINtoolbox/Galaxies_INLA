# Loads necessary R libraries and scripts: 'inla_fct.R':
library(INLA)
require(lattice)
require(classInt)
library(FITSio)
source("inla_fct.R")


# Read a STARLIGHT map from a FITS file and save x, y and mass variables:

massfile <- readFITS(file='data/NGC0309_mass.fits')
mass <- massfile$imDat[,]
dim <- dim(mass)
x <- array(0,dim=dim)
y <- array(0,dim=dim)
for (i in 1:dim[1]){
  for (j in 1:dim[2]){
    x[i,j] <- i
    y[i,j] <- j
  }
}


# Remove non-physical data:
valid <- which(!mass == 0)


# Run INLA pipeline on the mass map using a stationary model, elliptical distance, 
# and a voronoi mesh cutoff  of 1: 

massinla <- stationary_inla(x[valid],y[valid],mass[valid],xsize=dim[2],cutoff=1,
             ysize=dim[1],weight=mass[valid],shape='ellipse')


# Starry Night color palette

Van_Gogh <- c("#263C8B","#4E74A6", "#BDBF78", "#BFA524", "#2E231F")
colpal<-colorRampPalette(Van_Gogh)(100)
cutColor <- 100

# Natural colour scheme
fj5_mass <- classIntervals(massinla$out[valid], n = cutColor, style= "fisher")
fj5_ermass <- classIntervals(massinla$outsd[valid], n = cutColor, style= "fisher")



# Display original data
levelplot(massinla$image,col.regions=colpal,
          cut=cutColor,at=fj5_mass$brks)



# Display mean INLA output 
   levelplot(massinla$out, 
          col.regions=colpal,
          cut=cutColor,at=fj5_mass$brks
         )


# Display standard error around the mean 
  levelplot(massinla$outsd,
           col.regions=colpal,
            cut=cutColor,
            at=fj5_ermass$brks)


