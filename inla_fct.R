

### METRIC MEASUREMENT
residualCalculate <- function(mod,obs,ermod){
 mean <- mean(NaRV.omit((as.vector(mod)-as.vector(obs))^2/(as.vector(ermod))^2))
 median <- median(NaRV.omit((as.vector(mod)-as.vector(obs))^2/(as.vector(ermod))^2))
 std <- sd(NaRV.omit((as.vector(mod)-as.vector(obs))^2/(as.vector(ermod))^2))
 mad <- mad(NaRV.omit((as.vector(mod)-as.vector(obs))^2/(as.vector(ermod))^2))
 emean <- mean(NaRV.omit((as.vector(mod)-as.vector(obs))^2))
 emedian <- median(NaRV.omit((as.vector(mod)-as.vector(obs))^2))
 estd <- sd(NaRV.omit((as.vector(mod)-as.vector(obs))^2))
 emad <- mad(NaRV.omit((as.vector(mod)-as.vector(obs))^2))
 return(list(mean=mean,median=median,std=std,mad=mad,
             emean=emean,emedian=emedian,estd=estd,emad=emad))
}

## RESHAPE IMAGE WHEN ZOOM SET
zoom_fix <- function(img,zoom){
    mim <- melt(img)
    colnames(mim) <- c("x","y","value")
    xx <- mim$x/zoom 
    yy <- mim$y/zoom
    zz <- mim$value
    return(zz~xx*yy)
}



## -------- 1) FUNCTION FOR STATIONARY FITS
stationary_inla <- function(tx, ty, tpar, tepar, weight=1, name='age',  zoom = 1,
                            xsize=77,ysize=72,shape='ellipse',
                            p_range=c(2,0.2),p_sigma=c(2,0.2),cutoff=5){

    
x <- tx
y <- ty
par <- tpar
if (hasArg(tepar)) { epar <- tepar^2 } else { epar <- NULL}

# Create a mesh (tesselation) 
mesh <- inla.mesh.2d(cbind(x,y), max.n = 10, cutoff = cutoff)

#bookkeeeping
A <- inla.spde.make.A(mesh, loc=cbind(x,y))

#priors gaussian process (prior range: 20% less than 2, sigma 20 higher than 2)
spde <- inla.spde2.pcmatern(mesh=mesh,alpha=2, prior.range=p_range,prior.sigma=p_sigma) 

#calculate projection from model
projection <- inla.mesh.projector(mesh,xlim=c(0,ysize),ylim=c(0,xsize),
                                   dim=zoom*c(ysize+1,xsize+1))

#center
xcenter <- sum(x*weight)/sum(weight)
ycenter <- sum(y*weight)/sum(weight)
    
#radius fct
if(shape == 'radius'){
    radius <- sqrt((x-xcenter)^2 + (y-ycenter)^2)
    radius_2 <- (x-xcenter)^2 + (y-ycenter)^2

    #use parametric function of  radius+radius^2
    stk_rad <- inla.stack(data=list(par=par), A=list(A,1,1,1),
                          effects=list(i=1:spde$n.spde, m=rep(1,length(x)),
                                       radius=radius, radius_2=radius_2),tag='est')
 
     #caculate result model
    res_rad <- inla(par ~ 0 + m +radius +radius_2 +f(i, model=spde),
                    data=inla.stack.data(stk_rad),
                    control.predictor=list(A=inla.stack.A(stk_rad)),scale=epar)

     #porjection for radius
    projected_radius <- sqrt((rep(projection$x,each=length(projection$y))-xcenter)^2 +
                             (rep(projection$y,length(projection$x)) -ycenter)^2)
    projected_radius_2  <- (rep(projection$x,each=length(projection$y))-xcenter)^2 +
        (rep(projection$y,length(projection$x)) -ycenter)^2

    #output with matrix to include ellipse function
    output <- inla.mesh.project(inla.mesh.projector(mesh,
                                                    xlim=c(0,ysize),ylim=c(0,xsize),
                                                    dim=zoom*c(ysize+1,xsize+1)),
                                res_rad$summary.random$i$mean)+
        t(matrix(as.numeric(res_rad$summary.fixed$mean[1]+
                            res_rad$summary.fixed$mean[2]*projected_radius+
                            res_rad$summary.fixed$mean[3]*projected_radius_2),
                 nrow=zoom*(xsize+1),ncol=zoom*(ysize+1)))

    #output std with simple (no function)
    outputsd <- inla.mesh.project(inla.mesh.projector(mesh,xlim=c(0,ysize),
                                                      ylim=c(0,xsize),
                                                      dim=zoom*c(ysize+1,xsize+1)),
                                  res_rad$summary.random$i$sd)
 
#ellipse fct
}
else if(shape=='ellipse') {

    covar <- cov.wt(cbind(x,y), w=weight)
    eigens <- eigen(covar$cov)
    ellipse <- (cbind(x-xcenter,y-ycenter)%*%(eigens$vectors[,1]))^2/eigens$values[1] +
        (cbind(x-xcenter,y-ycenter)%*%(eigens$vectors[,2]))^2/eigens$values[2]
    ellipse_2 =  ellipse^2

    #use parametric function of ellipse & ellipse^2 
    stk_ell <- inla.stack(data=list(par=par), A=list(A,1,1,1),
                          effects=list(i=1:spde$n.spde,
                                       m=rep(1,length(x)),ellipse=ellipse,
                                       ellipse_2=ellipse_2),tag='est')
    #caculate result model
    res_ell <- inla(par ~ 0 + m +ellipse +ellipse_2 +f(i, model=spde),
                    data=inla.stack.data(stk_ell),
                    control.predictor=list(A=inla.stack.A(stk_ell)),scale=epar)

    #print restuls
    #print(res_rad$summary.fix)

    #prjection for ellipse
    px = rep(projection$x,each=length(projection$y))
    py = rep(projection$y,length(projection$x))
    projected_ellipse <- (cbind(px-xcenter,py-ycenter)%*%(eigens$vectors[,1]))^2/eigens$values[1] +
		  (cbind(px-xcenter,py-ycenter)%*%(eigens$vectors[,2]))^2/eigens$values[2]
    projected_ellipse_2 <- projected_ellipse^2

    #output with matrix to include ellipse function
    output <- inla.mesh.project(inla.mesh.projector(mesh,
                                                    xlim=c(0,ysize),ylim=c(0,xsize), #ch
                                                    dim=zoom*c(ysize+1,xsize+1)),#ch
                                res_ell$summary.random$i$mean)+
        t(matrix(as.numeric(res_ell$summary.fixed$mean[1]+
                            res_ell$summary.fixed$mean[2]*projected_ellipse+
                            res_ell$summary.fixed$mean[3]*projected_ellipse_2),
                 nrow=zoom*(xsize+1),ncol=zoom*(ysize+1)))#chan

    #output std with simple (no function)
    outputsd <- inla.mesh.project(inla.mesh.projector(mesh,xlim=c(0,ysize),
                                                      ylim=c(0,xsize),
                                                      dim=zoom*c(ysize+1,xsize+1)),
                                  res_ell$summary.random$i$sd)
}
    

#zoom    
if (zoom != 1){
    output <- zoom_fix(output,zoom)
    outputsd <- zoom_fix(outputsd,zoom)
}
    
    
#original data to compare
timage <- matrix(NA,nrow=ysize+1,ncol=xsize+1)
for (i in 1:length(x)) {timage[x[i],y[i]] <- par[i]}
terrimage = NULL
if (hasArg(tepar)) {
    terrimage <- matrix(NA,nrow=ysize+1,ncol=xsize+1)
    for (i in 1:length(x)) {terrimage[x[i],y[i]] <- epar[i]}
}
    

#sanity check (projected age compared to observed age as a funciton of age_error (should increase))
#plot(age_error,as.numeric(A%*%res$summary.random$i$mean)-age)
    return(list(out=output, image=timage, erimage=terrimage,outsd=outputsd))
}


### ----------- 2) FUNCTION FOR NON-STATIONARY FITS
nonstationary_inla <- function(tx, ty, tpar, tepar, weight=1, name='age',zoom=1,
                        xsize=77,ysize=72,cutoff=5,nbasis=2,degree=2,shape='ellipse'){

x <- tx
y <- ty
par <- tpar
if (hasArg(tepar)) { epar <- tepar^2 } else { epar <- NULL}

# Create a mesh (tesselation) 
mesh <- inla.mesh.2d(cbind(x,y), max.n = 10, cutoff = cutoff)

#bookkeeeping
A <- inla.spde.make.A(mesh, loc=cbind(x,y))

#inverse scale: degree=10, degree=10, n=2, n=2
basis.T <- inla.mesh.basis(mesh,type="b.spline", n=nbasis, degree=degree)
#inverse range
basis.K <- inla.mesh.basis(mesh,type="b.spline", n=nbasis, degree=degree)

spde <- inla.spde2.matern(mesh=mesh, alpha=2,
B.tau=cbind(0,basis.T,basis.K*0),B.kappa=cbind(0,basis.T*0,basis.K/2))

#calculate projection from model
projection <- inla.mesh.projector(mesh,xlim=c(0,xsize),ylim=c(0,ysize),dim=c(xsize+1,ysize+1))
    
#center
xcenter <- sum(x*weight)/sum(weight)
ycenter <- sum(y*weight)/sum(weight)    


#radius fct
if(shape == 'radius'){
    radius <- sqrt((x-xcenter)^2 + (y-ycenter)^2)
    radius_2 <- (x-xcenter)^2 + (y-ycenter)^2

    #use non-par fct of radius+radius_2
    stk_rad <- inla.stack(data=list(par=par), A=list(A,1,1,1),effects=list(i=1:spde$n.spde,radius=radius, radius_2=radius_2,m=rep(1,length(x))),tag='est')

    #calculate results
    res_rad <- inla(par ~ 0 + m +radius +radius_2 +f(i, model=spde),
         data=inla.stack.data(stk_rad), control.predictor=list(A=inla.stack.A(stk_rad)),
	scale=epar)

    #for radius
    projected_radius <- sqrt((rep(projection$x,each=length(projection$y))-xcenter)^2 +
		 (rep(projection$y,length(projection$x)) -ycenter)^2)
    projected_radius_2 <- (rep(projection$x,each=length(projection$y))-xcenter)^2 +
        (rep(projection$y,length(projection$x)) -ycenter)^2

    #output
    output <- inla.mesh.project(inla.mesh.projector(mesh,
	   xlim=c(0,xsize),ylim=c(0,ysize),dim=c(xsize+1,ysize+1)),res_rad$summary.random$i$mean)+
	   t(matrix(as.numeric(res_rad$summary.fixed$mean[1]+
	   res_rad$summary.fixed$mean[2]*projected_radius+
	   res_rad$summary.fixed$mean[3]*projected_radius_2),nrow=ysize+1,ncol=xsize+1))

    #output std with simple (no function)
    outputsd <- inla.mesh.project(inla.mesh.projector(mesh,xlim=c(0,xsize),
                  ylim=c(0,ysize),dim=c(xsize+1,ysize+1)),res_rad$summary.random$i$sd)

#ellipse fct
}
else if(shape=='ellipse') {  

    covar <- cov.wt(cbind(x,y), w=weight)
    eigens <- eigen(covar$cov)
    ellipse <- (cbind(x-xcenter,y-ycenter)%*%(eigens$vectors[,1]))^2/eigens$values[1] +(cbind(x-xcenter,y-ycenter)%*%(eigens$vectors[,2]))^2/eigens$values[2]
    ellipse_2 =  ellipse^2

    #nonpar fct of ellipse & ellipse^2
    stk_ell <- inla.stack(data=list(par=par), A=list(A,1,1,1),effects=list(i=1:spde$n.spde,ellipse=ellipse, ellipse_2=ellipse_2,m=rep(1,length(x))),tag='est')

    #result
    res_ell <- inla(par ~ 0 + m +ellipse +ellipse_2 +f(i, model=spde),
	data=inla.stack.data(stk_ell), control.predictor=list(A=inla.stack.A(stk_ell)),
	scale=epar)

    #print restuls
    #res_rad$summary.fix

    #for ellipse
    px = rep(projection$x,each=length(projection$y))
    py = rep(projection$y,length(projection$x))
    projected_ellipse <- (cbind(px-xcenter,py-ycenter)%*%(eigens$vectors[,1]))^2/eigens$values[1] +
		  (cbind(px-xcenter,py-ycenter)%*%(eigens$vectors[,2]))^2/eigens$values[2]
    projected_ellipse_2 <- projected_ellipse^2

    #output
    output <- inla.mesh.project(inla.mesh.projector(mesh,
	   xlim=c(0,xsize),ylim=c(0,ysize),dim=c(xsize+1,ysize+1)),res_ell$summary.random$i$mean)+
	   t(matrix(as.numeric(res_ell$summary.fixed$mean[1]+
	   res_ell$summary.fixed$mean[2]*projected_ellipse+
	   res_ell$summary.fixed$mean[3]*projected_ellipse_2),nrow=ysize+1,ncol=xsize+1))

    #output std with simple (no function)
    outputsd <- inla.mesh.project(inla.mesh.projector(mesh,xlim=c(0,xsize),ylim=c(0,ysize),dim=c(xsize+1,ysize+1)),res_ell$summary.random$i$sd)
}

#zoom    
if (zoom != 1){
    output <- zoom_fix(output,zoom)
    outputsd <- zoom_fix(outputsd,zoom)
}
    
    
#original data to compare
timage <- matrix(NA,nrow=ysize+1,ncol=xsize+1)
for (i in 1:length(x)) {timage[x[i],y[i]] <- par[i]}
if (hasArg(tepar)) {
terrimage <- matrix(NA,nrow=ysize+1,ncol=xsize+1)
for (i in 1:length(x)) {terrimage[x[i],y[i]] <- epar[i]}
}



return(list(out=output, image=timage, erimage=terrimage, outsd=outputsd))
}

#### --------------- 3) NON-PARAMETRIC FUNCTION
##Instead of predicting ourselves, we can pass a larger array with NA values that will be automatically predicted

nonparametric_inla <- function(tx, ty, tpar, tepar, weight=1, name='age',xsize=77,ysize=70){

tweight <- weight

#valid data
valid <- which(!(tpar <= 0) & !(tepar <= 0) & !(tepar == "NaN") & !(tweight <=0))
x <- c(tx[valid],rep(1:(xsize+1),each=(ysize+1)))
y <- c(ty[valid],rep(1:(ysize+1),xsize+1))
par <- c(tpar[valid],rep(NA,length=((xsize+1)*(ysize+1))))
epar <- c(tepar[valid],rep(NA,length=((xsize+1)*(ysize+1))))
weight <- c(tweight[valid],rep(NA,length=((xsize+1)*(ysize+1))))

#radius fct
xcenter <- sum(x[1:length(valid)]*weight[1:length(valid)])/sum(weight[1:length(valid)])
ycenter <- sum(y[1:length(valid)]*weight[1:length(valid)])/sum(weight[1:length(valid)])
radius <- sqrt((x-xcenter)^2 + (y-ycenter)^2)
radius_2 <- (x-xcenter)^2 + (y-ycenter)^2

#ellipse fct
covar <- cov.wt(cbind(x[1:length(valid)],y[1:length(valid)]), w=weight[1:length(valid)])
eigens <- eigen(covar$cov)
ellipse <- (cbind(x-xcenter,y-ycenter)%*%(eigens$vectors[,1]))^2/eigens$values[1] +(cbind(x-xcenter,y-ycenter)%*%(eigens$vectors[,2]))^2/eigens$values[2]
ellipse_2 =  ellipse^2

mesh <- inla.mesh.2d(cbind(x,y), max.n = 10, cutoff = 5)

A <- inla.spde.make.A(mesh, loc=cbind(x,y))

#inverse scale
basis.T <- inla.mesh.basis(mesh,type="b.spline", n=2, degree=2)
#inverse range
basis.K <- inla.mesh.basis(mesh,type="b.spline", n=2, degree=2)

spde <- inla.spde2.matern(mesh=mesh, alpha=2,
B.tau=cbind(0,basis.T,basis.K*0),B.kappa=cbind(0,basis.T*0,basis.K/2))

stk <- inla.stack(data=list(par=par), A=list(A,1,1,1),effects=list(i=1:spde$n.spde,
    ellipse=ellipse,ellipse_2=ellipse_2,m=rep(1,length(x))),tag='est')

##rw2=random_walk, ou=Ornstein-Uhlenbeck
## putting "compute=TRUE" will predict the rest of NA values

res <- inla(par ~ 0 + m + f(ellipse,model="ou") +f(i, model=spde), data=inla.stack.data(stk),
control.predictor=list(A=inla.stack.A(stk),compute=TRUE),scale=epar^2,quantiles=0.5)

#output <- inla.mesh.project(inla.mesh.projector(mesh,xlim=c(0,xsize),ylim=c(0,ysize),
#       dim=c(xsize+1,ysize+1)),t(matrix(as.numeric(res$summary.fitted.values$mean[(length(#valid)+1):length(x)]),nrow=ysize+1,ncol=xsize+1)))

output <- t(matrix(as.numeric(res$summary.fitted.values$mean[(length(valid)+1):length(x)]),nrow=ysize+1,ncol=xsize+1))



##chec radial nonpar fct
plot(res$summary.random$i$ID,res$summary.random$i$mean)


#original data to compare
timage <- matrix(NA,nrow=ysize+1,ncol=xsize+1)
for (i in 1:length(valid)) {timage[x[i],y[i]] <- par[i]}
terrimage <- matrix(NA,nrow=ysize+1,ncol=xsize+1)
for (i in 1:length(valid)) {terrimage[x[i],y[i]] <- epar[i]}


name='NGC2906-age-0.05'
ra=range(par[1:length(valid)])
png(filename=paste(name,'INLA_nonparametric.png',sep='_'))
par(mfrow=c(2,1))
image(timage,col=hsv(100:1/100*0.6),zlim=ra)
title(name)
image(output,col=hsv(100:1/100*0.6),zlim=ra)
title(paste(name,'INLA:ELLIPSE',sep='_'))
dev.off()
}



### QUICK AND DIRTY INTERPOLATION FROM ALBERTO

quickAndDirtyKrigforInterpolation2 <- function(x, y, z, zlab="Log(Age)") {
 require(fields)
 fit <- Krig(data.frame(x, y), z, theta=20)
 summary(fit)
 # Some diagnostic plots
 #set.panel(2,2)
 #plot(fit)
 quartz()
 surface( fit, type="I", extrap=FALSE, nx=200, ny=200, main=paste(zlab, "\n Simple GRSP Kriging"), zlab=zlab, col=terrain.colors(256), asp=1)
 points( fit$x, pch=19, cex=0.2 )
 return(fit)
}
quickAndDirtyKrigforInterpolationErrorPlot <- function(fit, zlab="Log(Age)") {
 quartz()
 sdSurface <- predictSurfaceSE(fit, extrap=FALSE)
 surface( sdSurface, type="I", extrap=FALSE, nx=200, ny=200, main=paste(zlab, "\n Simple GRSP Kriging Error Estimate"), col=terrain.colors(256), zlab=zlab, asp=1)
 points( fit$x, pch=19, cex=0.2 )
}

quickAndDirtyKrigforInterpolationWithErrorEstimate <- function(x, y, z, zlab="Log(Age)") {
 res <- quickAndDirtyKrigforInterpolation2(x, y, z, zlab=zlab)
 quickAndDirtyKrigforInterpolationErrorPlot(res, zlab=zlab)
}
