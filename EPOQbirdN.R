
library(data.table)
library(sp)
library(readxl)
library(glmmTMB)
library(mapview)
library(sf)
library(visreg)
library(geoR)
library(MASS)
library(mgcv)
library(scales)
library(raster)
library(velox)
library(mgcViz)
library(DHARMa)
library(FRutils)
library(viridisLite)
library(rasterVis)

windows()

##############################
### read data

d<-fread("C:/Users/rouf1703/Documents/Contrats/Antoinette/data/BDAllLU2.csv")

### I remove these habitats in case we want to change something in their measurement and it creates problems further down
d<-d[,!names(d)%in%c("SNE", "URB", "PAR", "PRA", "AGR", "FMX", "WAT", "WET","FCO", "FFE","FFM"),with=FALSE]


#d<-fread("C:/Users/rouf1703/Documents/Contrats/Antoinette/data/Extraction_LULC_EPOQ.txt")
d$lat<-as.numeric(gsub(",",".",d$Latitude))
d$lon<-as.numeric(gsub(",",".",d$Longitude))
d$date<-as.Date(paste(substr(d$Date,1,4),substr(d$Date,5,6),substr(d$Date,7,8),sep="-"))
d$jul<-as.integer(format(d$date,"%j"))
d$year<-as.integer(substr(d$date,1,4))
d$cell<-as.integer(factor(paste(d$lon,d$lat)))

g<-grep("Esp",names(d)) # which columns are species, to be used later

#################################
### use species name
codes<-read.csv("https://raw.githubusercontent.com/frousseu/EPOQbirdN/master/CodeEPOQ.csv")
names(d)[match(gsub("_","",codes[,"Nom_Var"]),names(d))]<-gsub(" ","_",codes[,"CodeFR"]) # replace species codes with names

species<-names(d)[g] # species names

###############################################
### change NAs to 0s (the data.table way)

for(sp in species){ 
  set(d,i=which(is.na(d[[sp]])),j=sp,value=0)
}

###################################
### checks nb records per species

n<-sapply(g,function(i){
  c(length(which(d[,..i]==0)),length(which(d[,..i]>0)))
})
o<-order(n[1,])
par(mar=c(10,6,2,2))
barplot(n[,o],names.arg=names(d)[g][o],border=NA,col=c("darkred","darkgreen"),las=2,ylab="Nombre de mentions")
legend("topleft",legend=c("  0",">0"),fill=c("darkred","darkgreen"),cex=3,bty="n",border=NA)

############################
### checks NAs per species

### now, supposed to be none

n<-sapply(g,function(i){
  length(which(is.na(d[,..i])))
})
o<-order(n)
par(mar=c(10,6,2,2))
barplot(n[o],names.arg=names(d)[g][o],border=NA,col="darkgreen",las=2,ylab="Nombre de NAs")


#######################################
### see locations

locs<-unique(d[,c("cell","lon","lat")])
coordinates(locs)<-~lon+lat
proj4string(locs)<-"+init=epsg:4326"
plot(locs,pch=16)
#mapview(locs)


#########################################
### look at data distribution
sp<-"BRCH"
hist(log(d[[sp]]))
plot(d$date,d[[sp]],col=gray(0.5,0.25),pch=16)


#####################################################
### Nombre de feuillets
x<-d[,.N,by=date]
plot(x$date,x$N,xaxt="n",xlab="Date",ylab="Nombre de feuillets",pch=16,col=gray(0.5,0.5))
r<-seq.Date(min(d$date),max(d$date),by="3 month")
axis.Date(1,at=r,las=2,format="%Y-%b")


#####################################################
### Nombre de feuillets par cell
x<-d[,.N,by=.(cell)]
locs$N<-x$N[match(locs$cell,x$cell)]
plot(locs,cex=sqrt(locs$N)/5,pch=16,col=gray(0.1,0.2),axes=TRUE)

### using mapview
mapview(locs,cex=sqrt(locs$N)/2,lwd=0.5)


#####################################################
### Nombre de feuillets par cell par date
d[,nbf:=.N,by=.(cell,date)]
hist(d$nbf,breaks=seq(0,100,by=1),xlim=range(d$nbf))


#####################################################
### Moyenne par espèces par cellule
sp<-d[,lapply(.SD,mean),by=cell,.SDcols = species]
locs@data<-cbind(locs@data,sp[match(locs$cell,sp$cell),])

par(mfrow=c(4,7),mar=c(0,0,0,0),oma=c(1,1,1,1))
for(i in seq_along(species)){
  plot(locs,pch=16,cex=rescale(sqrt(locs@data[,species[i]]),to=c(0.05,5)),col=gray(0.5,0.5))
  mtext(species[i],side=3,adj=0.1,line=-0.5,font=2,cex=1,col=gray(0,0.5))
}


#################################
### variog from previous models

### don't run, can be slow!

#dxs<-dx
#coordinates(dxs)<-~lon+lat
#proj4string(dxs)<-"+init=epsg:4326"
#dxs<-spTransform(dxs,CRS("+init=epsg:26917"))
#coords <- coordinates(dxs)
# v<-variog(coords=coords,data=resid(m),breaks=seq(0,50000,by=2000),max.dist=50000,bin.cloud=TRUE)
# par(mfrow=c(1,1))
# plot(v, main = "Variogram for spatial autocorrelation (LFY, fire occurrence)",type="b")



#######################################
### open canada raster data

### ACI code conversion
aci<-read.csv("https://raw.githubusercontent.com/frousseu/EPOQbirdN/master/CodeACI.csv",stringsAsFactors=FALSE)

### ACI data from Open Canada
# http://www.agr.gc.ca/atlas/data_donnees/agr/annualCropInventory/tif/2015/
r<-raster("C:/Users/rouf1703/Downloads/aci_2015_qc/aci_2015_qc.tif")
r<-crop(r,extend(extent(spTransform(locs,CRS(proj4string(r)))),c(5000,5000,5000,5000)))


### visualize LU types
# silenced because slow
#rr<-ratify(r)
#rr<-subs(r,aci[,c("ACI","Code")]) ### slow and memory hungry!
#cols<-c("grey90","cadetblue1","grey20","grey70","lightgoldenrod","chocolate","brown","goldenrod","chartreuse3","darkgreen","chartreuse1")
#levels(rr)[[1]]$cols<-cols
#plot(rr)
#levelplot(rr,maxpixels=1e6,col.regions=cols) # pixels are aggregated for speed, memory intensive...



##################################
### build grid with approcximate half size of cells in lat/lon

tablon<-table(diff(coordinates(locs)[,1]))
tablon<-tablon[names(tablon)!="0"]
tablon<-tablon[order(names(tablon))]

tablat<-table(diff(coordinates(locs)[,2]))
tablat<-tablat[names(tablat)!="0"]
tablat<-tablat[order(names(tablat))]

### w will be used to represent the half width in latlon coordinates of a cell
k<-1:4
w<-weighted.mean(abs(as.numeric(c(names(tablon)[k],names(tablon)[k]))),w=c(tablon[k],tablat[k]))/2 # some centroids are separated by 0.016, others by 0.017, impossible to make a regular grid. A weighted mean is used to make some kind of compromise

cells<-SpatialPolygons(apply(apply(coordinates(locs),1,function(r){cbind(r+c(w,w),r+c(w,-w),r+c(-w,-w),r+c(-w,w),r+c(w,w))}),2,function(v){Polygons(list(Polygon(matrix(v,ncol=2,byrow=TRUE))),ID=runif(1))}))
proj4string(cells)<-proj4string(locs)
#mapview(cells)
cells<-spTransform(cells,CRS(proj4string(r)))
#cells<-st_make_grid(st_as_sf(as(extent(r),"SpatialPolygons")),cell=1000)
#plot(cells)


####################################
### extract pixels
v<-velox(r)
e<-v$extract(cells)
le<-lapply(e,function(i){
  tab<-table(i[,1])	
  as.data.frame(rbind(tab/sum(tab)))
})


##########################################################
### bind cells data and fill unrepresented classes

lu<-rbindlist(le,fill=TRUE)
for(i in names(lu)){ 
  set(lu,i=which(is.na(lu[[i]])),j=i,value=0)
}
names(lu)<-aci$Code[match(names(lu),aci$ACI)]


############################################################
### combine columns of the same classes with an ugly hack (should use something with data.table
clu<-unique(names(lu))
temp<-list()
for(i in seq_along(clu)){
  m<-which(names(lu)%in%clu[i])
  if(length(m)>1){
    temp[[i]]<-rowSums(lu[,..m])
  }else{
    temp[[i]]<-lu[[m]]	
  }
}
lu<-do.call("data.table",temp)
names(lu)<-clu
lu$cell<-locs$cell
table(rowSums(lu[,1:(ncol(lu)-1)]),useNA="ifany")


####################################
### add data to the locs data and to the data

locs@data<-cbind(locs@data,lu)
plot(locs,pch=16,cex=locs$AGR*2)

d<-merge(d,lu)
d$FFM<-d$FFE+d$FMX


#####################################
### create prediction grid
#####################################

ext<-1 # 1 cell extension over the observation grid
inc<-1 # increase resolution of prediction grid by this factor compared to original grid
# need to fix grid ? uneven squares from EPOQ

gridlocs<-as.matrix(expand.grid(seq(min(coordinates(locs)[,1])-ext*2*w,max(coordinates(locs)[,1])+ext*2*w,by=2*w/inc),seq(min(coordinates(locs)[,2])-ext*2*w,max(coordinates(locs)[,2])+ext*2*w,by=2*w/inc)))
gridl<-SpatialPolygons(apply(apply(gridlocs,1,function(r){cbind(r+c(w,w),r+c(w,-w),r+c(-w,-w),r+c(-w,w),r+c(w,w))}),2,function(v){Polygons(list(Polygon(matrix(v,ncol=2,byrow=TRUE))),ID=runif(1))}))
proj4string(gridl)<-proj4string(locs)
grid<-spTransform(gridl,CRS(proj4string(r)))
grid<-SpatialPolygonsDataFrame(grid,data=data.frame(lon=coordinates(gridl)[,1],lat=coordinates(gridl)[,2]))
plot(cells,axes=TRUE,xlim=c(1800000,1820000),ylim=c(830000,850000),border="white",col="blue")
plot(grid,add=TRUE,border="red")
points(coordinates(grid),col="red")
points(spTransform(locs,CRS(proj4string(r))))

v<-velox(r)
e<-v$extract(grid)
le<-lapply(e,function(i){
  tab<-table(i[,1])	
  as.data.frame(rbind(tab/sum(tab)))
})

lu<-rbindlist(le,fill=TRUE)
for(i in names(lu)){ 
  set(lu,i=which(is.na(lu[[i]])),j=i,value=0)
}
names(lu)<-aci$Code[match(names(lu),aci$ACI)]

clu<-unique(names(lu))
temp<-list()
for(i in seq_along(clu)){
  ww<-which(names(lu)%in%clu[i])
  if(length(ww)>1){
    temp[[i]]<-rowSums(lu[,..ww])
  }else{
    temp[[i]]<-lu[[ww]]	
  }
}
lu<-do.call("data.table",temp)
names(lu)<-clu
table(rowSums(lu),useNA="ifany")
lu$lon<-coordinates(gridl)[,1]
lu$lat<-coordinates(gridl)[,2]
lu$FFM<-lu$FFE+lu$FMX
lu$cell<-1


#####################################################
### Models
#####################################################

### choose species and subset size
sp<-"PACO" 
dx<-d[sample(1:nrow(d),5000),]

f<-as.formula(paste(sp,"~","s(jul,bs=\"cc\")+s(lon,lat,bs=c(\"ds\"))+te(jul,lat,bs=c(\"cc\",\"ds\"))+s(FFM)+s(WAT)+s(URB)+s(AGR)+s(PAR)+s(cell,bs=\"re\")"))

dxs<-dx
coordinates(dxs)<-~lon+lat
proj4string(dxs)<-proj4string(locs)
dxs<-spTransform(dxs,CRS(proj4string(r)))

### model
m<-gam(f,data=dx,family=nb,select=FALSE,method="REML")
summary(m)
gam.check(m)

### response variable
barplot(table(model.frame(m)[,1]),col="darkgreen",border=NA)

### effects of single variables
va<-all.vars(f[[3]])
va<-va[va!="cell"]
par(mfrow=rep(ceiling(sqrt(length(va))),2),mar=c(4,4,3,3))
for(i in va){
  visreg(m,i,scale="response",xlab=i,rug=FALSE)
}
par(mfrow=c(1,1))


### effects of variables in interactions
visreg(m,"jul","lat",scale="response",overlay=TRUE)
visreg2d(m,"lon","lat",scale="response")


### compare simulated data to observed data using mgcvViz
b <- getViz(m, nsim = 100, post = TRUE, unconditional = TRUE)
check0D(b)+l_hist(binwidth=1) + scale_x_continuous(limits = c(-3, 5))


### simulate data using mgcViz
sims<-t(simulate(m,100))


### check residuals with DHARMa
dsims<-createDHARMa(simulatedResponse = sims, observedResponse = model.frame(m)[,1],fittedPredictedResponse = predict(m),integerResponse=TRUE)
plot(dsims,quantreg=FALSE)
testDispersion(dsims)


### compare proportions of zeros using simulated and observed data
sims0<-apply(sims,2,function(i){sum(i==0)/length(i)})
obs0<-sum(model.frame(m)[,1]==0)/length(model.frame(m)[,1])
hist(sims0,xlab="Proportions of 0s",border="white",col="grey70",main="",xlim=range(c(sims0,obs0)))
abline(v=obs0,lwd=2,col="red")


### compare simulated to observed data
plot(density(model.frame(m)[,1],adjust=2),xlim=c(0,min(50,max(model.frame(m)[,1]))),col="red",lwd=2)
invisible(
  lapply(1:ncol(sims),function(i){
    lines(density(sims[,i],adjust=2),col=gray(0,0.1))
  })
)



#######################################
### mapping predictions
#######################################

### choose julian dates and value to map
juls<-c(130,140,150,160)
vals<-c("pfit","pse")[1] # choose "pfit" [1] or "pse" [2] to visualize predictions or SE 

### generate values to predict
newdat<-as.data.frame(cbind(lu[rep(1:nrow(lu),length(juls)),],jul=rep(juls,each=nrow(lu))))
p<-predict(m,newdata=newdat,type="response",exclude = "s(cell)",se.fit=TRUE)
pfit<-as.vector(p$fit)
pse<-as.vector(p$se)
newdat<-cbind(newdat,data.frame(pfit,pse))

### choose color scale and split dates
newdat$p<-newdat[,vals]
newdat$col<-colo.scale(rescale(newdat$p,to=c(1,0)),inferno(100))
lnewdat<-split(newdat,newdat$jul)

### map predictions
par(mfrow=c(floor(sqrt(length(juls))),ceiling(sqrt(length(juls)))),mar=c(3,3,0,0),oma=c(0,0,0,5))
for(i in seq_along(lnewdat)){
  x<-lnewdat[[i]]
  plot(spTransform(grid,CRS(proj4string(locs))),col=x[,"col"],border=NA,axes=FALSE)
  mtext(paste(sp,format(as.Date(juls[i],origin="1970-01-01"),"%B-%d")),side=3,line=-3,adj=c(0.5),padj=0,font=2)
  #points(dxs$lon,dxs$lat,cex=10*(dxs[[sp]]+0.5)/max(dxs[[sp]]))
  #plot(locs,pch=1,cex=rescale(sqrt(locs@data[,sp]),to=c(0.05,5)),col=gray(0,0.75),add=TRUE)
  #grid$zcol<-x[,"pfit"]
  #mapview(spTransform(grid,CRS(proj4string(locs))),zcol="zcol")
}

### add legend
legend_image <- as.raster(matrix(inferno(100), ncol=1))
rec<-c(-72.6,45.01,-72.5,46.4)
par(new=TRUE,mfrow=c(1,1),mar=c(3,3,0,5),oma=c(0,0,0,0))
plot(0,0,type="n",xaxt="n",yaxt="n",bty="n")
rec<-c(par("usr")[2]+0.01,-1,1.13,1)
rasterImage(legend_image,rec[1],rec[2],rec[3],rec[4],xpd=TRUE,lwd=1)
brks<-20
at<-seq(min(newdat$p),max(newdat$p),length.out=brks)
labs<-rescale(c(at,range(newdat$p)),to=c(1,0))
labs<-at[1:brks]
text(x=rec[3],y= seq(rec[2],rec[4],l=brks), labels = round(labs,1),adj=c(-0.2,0.5),font=1,xpd=TRUE)


##############################################################
### map spatial effect assuming constant vars across cells
##############################################################

vals<-c("pfit","pse")[1] # choose "pfit" [1] or "pse" [2] to visualize predictions or SE 

### generate values to predict
newdat<-as.data.frame(cbind(lu[rep(1:nrow(lu),length(juls)),],jul=rep(juls,each=nrow(lu))))

newdat<-as.data.frame(as.list(colMeans(lu[,!names(lu)%in%c("lon","lat","cell","jul","0"),with=FALSE])))
newdat<-cbind(newdat[rep(1,nrow(lu)),],lu[,c("lat","lon")],jul=150,cell=1)

p<-predict(m,newdata=newdat,type="response",exclude = "s(cell)",se.fit=TRUE)
pfit<-as.vector(p$fit)
pse<-as.vector(p$se)
newdat<-cbind(newdat,data.frame(pfit,pse))

### choose color scale and split dates
newdat$p<-newdat[,vals]
newdat$col<-colo.scale(rescale(newdat$p,to=c(1,0)),inferno(100))

### map predictions
plot(spTransform(grid,CRS(proj4string(locs))),col=newdat[,"col"],border=NA,axes=FALSE)
mtext(paste(sp),side=3,line=-3,adj=c(0.5),padj=0,font=2)

legend_image <- as.raster(matrix(inferno(100), ncol=1))
rec<-c(-72.6,45.01,-72.5,46.4)
par(new=TRUE,mfrow=c(1,1),mar=c(3,3,0,5),oma=c(0,0,0,0))
plot(0,0,type="n",xaxt="n",yaxt="n",bty="n")
rec<-c(par("usr")[2]+0.01,-1,1.13,1)
rasterImage(legend_image,rec[1],rec[2],rec[3],rec[4],xpd=TRUE,lwd=1)
brks<-20
at<-seq(min(newdat$p),max(newdat$p),length.out=brks)
labs<-rescale(c(at,range(newdat$p)),to=c(1,0))
labs<-at[1:brks]
text(x=rec[3],y= seq(rec[2],rec[4],l=brks), labels = round(labs,1),adj=c(-0.2,0.5),font=1,xpd=TRUE)


#######################################
### interactive mapping
#######################################

map<-spTransform(grid,CRS(proj4string(locs)))
map$p<-pfit
mapviewOptions(raster.palette = grey.colors,vector.palette = colorRampPalette(rev(inferno(20))),na.color = "grey")
mapview(map,zcol="p",legend=TRUE,lwd=0,at = pretty(map$p,20))+mapview(locs,cex=sp,lwd=0,alpha.regions = 0.7, aplha = 0.2)


#####################################################
### test bs="re" effect on smooth across jul
#####################################################

jul<-1:365
ng<-1:100
newdat<-cbind(model.frame(m)[rep(ng,each=length(jul)),-(1:2)],jul=rep(jul,length(ng)))
lnd<-split(newdat,newdat$cell)
plot(0,0,xlim=c(0,365),ylim=c(0,15),type="n")
lvis<-lapply(lnd,function(i){
  p<-predict(m,newdata=i,type="response",se.fit=FALSE)
  lines(i$jul,as.vector(p))
})





