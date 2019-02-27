
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



##############################
### read data

d<-fread("C:/Users/User/Documents/Antoinette/Extraction_LULC_EPOQ.txt")
d$lat<-as.numeric(gsub(",",".",d$Latitude))
d$lon<-as.numeric(gsub(",",".",d$Longitude))
d$date<-as.Date(paste(substr(d$Date,1,4),substr(d$Date,5,6),substr(d$Date,7,8),sep="-"))
d$jul<-as.integer(format(d$date,"%j"))
d$year<-as.integer(substr(d$date,1,4))
d$cell<-as.integer(factor(paste(d$lon,d$lat)))

g<-grep("Esp_",names(d)) # which columns are species, to be used later

#################################
### use species name
codes<-read.csv("https://raw.githubusercontent.com/frousseu/EPOQbirdN/master/CodeEPOQ.csv")
names(d)[match(codes[,"Nom_Var"],names(d))]<-gsub(" ","_",codes[,"CodeFR"]) # replace species codes with names

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

#names(d)[grep("%",names(d))]<-gsub("%","P_",names(d)[grep("%",names(d))])
#locs<-unique(d[,c("cell","lon","lat",names(d)[grep("P_",names(d))])])
locs<-unique(d[,c("cell","lon","lat","%Agricole", "%BatiCommercial","%BatiResidentiel", "%Eau", "%MilieuNaturel", "%NonClassifie")])
coordinates(locs)<-~lon+lat
proj4string(locs)<-"+init=epsg:4326"
plot(locs,cex=as.numeric(gsub(",",".",locs$"%BatiResidentiel"))/100*2,pch=16)
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


#####################################################
### Simple models

sp<-"BRCH"

dx<-d[d$year%in%c(2002:2015),]
dx<-dx[sample(1:nrow(dx),10000),]
vars<-c("jul")

f<-as.formula(paste(sp,"~","s(jul,bs=\"cc\")+te(lon,lat)"))

m<-gam(f,data=dx,family=nb)
summary(m)

par(mfrow=c(1,1))
#visreg(m,scale="response")

visreg2d(m,"lon","lat",scale="response")
points(dx$lon,dx$lat,cex=10*(dx[[sp]]+0.5)/max(dx[[sp]]))


#####################################################
### multi species models 

dx<-d[d$year%in%c(2002:2015),]
dx<-dx[sample(1:nrow(dx),5000),]
dx<-melt(dx,measure.vars=species,variable.name = "species", value.name = "count")
dx<-dx[dx$Mangeoire==0,]
dx<-dx[dx$species%in%c("BRCH","BRFM","PAJA","MOCH","PACO","VIYR","VIME"),]

m<-gam(count~s(jul,by=species,bs="cc"),data=dx,family=nb)
summary(m)

par(mfrow=c(1,1))
visreg(m,"jul",by="species",scale="response",overlay=TRUE)




#################################
### variog from previous models

### don't run, can be slow!

dxs<-dx
coordinates(dxs)<-~lon+lat
proj4string(dxs)<-"+init=epsg:4326"
dxs<-spTransform(dxs,CRS("+init=epsg:26917"))
coords <- coordinates(dxs)
# v<-variog(coords=coords,data=resid(m),breaks=seq(0,50000,by=2000),max.dist=50000,bin.cloud=TRUE)
# par(mfrow=c(1,1))
# plot(v, main = "Variogram for spatial autocorrelation (LFY, fire occurrence)",type="b")



#######################################
### open canada raster data


### ACI code conversion
aci<-read.csv("https://raw.githubusercontent.com/frousseu/EPOQbirdN/master/CodeACI.csv",stringsAsFactors=FALSE)


### ACI data from Open Canada
# http://www.agr.gc.ca/atlas/data_donnees/agr/annualCropInventory/tif/2015/
r<-raster("C:/Users/User/Documents/Antoinette/aci_2015_qc/aci_2015_qc.tif")
r<-crop(r,extent(spTransform(locs,CRS(proj4string(r)))))
#plot(r)
#plot(spTransform(locs,CRS(proj4string(r))),add=TRUE,pch=1)


##################################
### build grid with approcximate half size of cells in lat/lon
w<-0.0165/2
cells<-SpatialPolygons(apply(apply(coordinates(locs),1,function(r){cbind(r+c(w,w),r+c(w,-w),r+c(-w,-w),r+c(-w,w),r+c(w,w))}),2,function(v){Polygons(list(Polygon(matrix(v,ncol=2,byrow=TRUE))),ID=runif(1))}))
proj4string(cells)<-proj4string(locs)
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
table(rowSums(lu),useNA="ifany")


####################################
### add data to the locs data and to the data

locs@data<-cbind(locs@data,lu)
plot(locs,pch=16,cex=locs$AGR*2)

d<-merge(d,lu)










