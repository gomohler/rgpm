# run this script to generate NIJ submission for RGPM.  Change TARGET_CRIME to target crime type.
# portland events data can be downloaded from https://nij.gov/funding/Pages/fy16-crime-forecasting-challenge-document.aspx#data
TARGET_CRIME="burglary"
data=read.csv("SADM/input_data/portland-events.csv")
library(lubridate)
library(reshape2)
library(rgl)
library(sp)
library(maptools)
library(raster)
library(rgeos)
library(ggmap)
library(rgdal)
library(glmnet)
library(randomForest)
set.seed(1)


pds=readShapeSpatial("SADM/input_data/portland-police-districts/Portland_Police_Districts.shp")
coordinates(data)<-~x+y
data$beat=over(data,pds)
nac=is.na(data$beat)
data=read.csv("SADM/input_data/portland-events.csv")
data=data[!nac[,1],]


minA=62500
maxA=360000
minTA=0.25*5280*5280
maxTA=.75*5280*5280
minH=125
maxH=500
PA=147.71*5280*5280

setupgrid<-function(data,h,w,shift,theta,fx,fy){
  
minx=min(data$x)
miny=min(data$y)
maxx=max(data$x)
maxy=max(data$y)
centerx=median(data$x)
centery=median(data$y)
Mx=ceiling((maxx-minx)/h)
My=ceiling((maxy-miny)/w)
Dx=ceiling(Mx*fx)*h
Dy=ceiling(My*fy)*w

grd <- GridTopology(cellcentre.offset=c(minx+shift-Dx,miny-Dy), cellsize=c(h,w), cells.dim=c(Mx+2*ceiling(Mx*fx),My+2*ceiling(My*fy)))
grdp <- as.SpatialPolygons.GridTopology(grd)
grid_rotate=elide(grdp, rotate = theta,center=c(centerx,centery))
return(grid_rotate)
}


setupfeatures<-function(data,YEAR,TARGET_CRIME){
data$datetime=ymd(as.character(data$date))
data$month=month(data$datetime)
data$split="valid"

UPPER=6
if(is.element(TARGET_CRIME,c("burglary","mvt"))){
UPPER=13
}

if(TARGET_CRIME=="all"){
data$split[data$year==YEAR&data$month>2&data$month<UPPER]="test"  
}else{
data$split[data$year==YEAR&data$month>2&data$month<UPPER&data$type==TARGET_CRIME]="test"
}

data$count=1
if(YEAR==2017){
  if(TARGET_CRIME=="all"){
  agg=aggregate(count~grid,data=data[data$year==2016&data$month>2&data$month<6,],FUN=sum)  
  }else{
  agg=aggregate(count~grid,data=data[data$year==2016&data$month>2&data$month<6&data$type==TARGET_CRIME,],FUN=sum)  
  }
}else{
agg=aggregate(count~grid,data=data[data$split=="test",],FUN=sum)
}
names(agg)=c("grid","count_test")

data$split[data$year<YEAR]="train"
data$split[data$year==YEAR&data$month<=2]="train"

ct=unique(data$type)
for(l in 1:4){
  grid_agg1=aggregate(count~grid,data=data[data$year==YEAR&data$month<=2&data$type==ct[l],],FUN=sum)
  names(grid_agg1)=c("grid",paste0("count1_",ct[l]))
  agg=merge(agg,grid_agg1,by="grid",all=TRUE)
  grid_agg2=aggregate(count~grid,data=data[data$year==(YEAR-1)&data$month>9&data$type==ct[l],],FUN=sum)
  names(grid_agg2)=c("grid",paste0("count2_",ct[l]))
  agg=merge(agg,grid_agg2,by="grid",all=TRUE)
  grid_agg3=aggregate(count~grid,data=data[data$year==(YEAR-1)&data$month<=9&data$type==ct[l],],FUN=sum)
  names(grid_agg3)=c("grid",paste0("count3_",ct[l]))
  agg=merge(agg,grid_agg3,by="grid",all=TRUE)
  grid_agg4=aggregate(count~grid,data=data[data$year<(YEAR-1)&data$type==ct[l],],FUN=sum)
  grid_agg4$count=grid_agg4$count/max(grid_agg4$count)
  names(grid_agg4)=c("grid",paste0("count4_",ct[l]))
  agg=merge(agg,grid_agg4,by="grid",all=TRUE)
  grid_agg5=aggregate(count~grid,data=data[data$split=="train"&data$type==ct[l]&data$month>2&data$month<6,],FUN=sum)
  grid_agg5$count=grid_agg5$count/max(grid_agg5$count)
  names(grid_agg5)=c("grid",paste0("count5_",ct[l]))
  agg=merge(agg,grid_agg5,by="grid",all=TRUE)
}

agg[is.na(agg)]<-0
return(agg)
}


PAIfun<-function(coef,data,TARGET_CRIME,fx,fy){
  h=coef[1]
  shift=coef[2]
  da=coef[3]
  theta=360*da
  
  minA=62500
  minTA=0.25*5280*5280
  PA=147.71*5280*5280

  w=minA/h
  Ncell=ceiling(minTA/(h*w))
  print("set up grid")
  
  grid_rotate<-setupgrid(data,h,w,shift,theta,fx,fy)
  
  print("assign grid cells to events")
  
  coordinates(data)<-~x+y
  data$grid=over(data,grid_rotate)
  
  print("calculate features")
  
  agg<-setupfeatures(data,2015,TARGET_CRIME)
  
  X=agg[,2:22]
  
  print("train random forest")
  model=randomForest(x=X[,2:21],y=log(X[,1]+1))
  agg$score=predict(model,as.matrix(X[,2:21]))
  agg=agg[order(-agg$score),] 
  
  PAI=sum(agg$count_test[1:Ncell])/sum(agg$count_test)*PA/(Ncell*h*w)
  print("current PAI:")
  print(PAI)
  if(h>500){PAI=0}
  if(h<125){PAI=0}
  
  return(PAI)
}

coefficients<-array(data=.1,dim=3)
coefficients[1]<-minH
coefficients[2]<-0
coefficients[3]<-0
est_par=optim(coefficients,PAIfun,data=data,TARGET_CRIME=TARGET_CRIME,fx=0,fy=0,control = list(maxit = 50,trace=1,fnscale=-1))

#validate

write.csv(est_par$par,paste0("~/Dropbox/Mohler/SADM/coefficients/",TARGET_CRIME,".csv"),row.names=F)
coef<-read.csv(paste0("~/Dropbox/Mohler/SADM/coefficients/",TARGET_CRIME,".csv"))

h=coef$x[1]
shift=coef$x[2]
da=coef$x[3]

theta=360*da


w=minA/h
Ncell=ceiling(minTA/(h*w))

grid_rotate<-setupgrid(data,h,w,shift,theta,.25,.25)
coordinates(data)<-~x+y
data$grid=over(data,grid_rotate)

agg<-setupfeatures(data,2015,TARGET_CRIME)

X=agg[,2:22]
model=randomForest(x=X[,2:21],y=log(X[,1]+1),ntree=2000)

agg<-setupfeatures(data,2016,TARGET_CRIME)

sort_crime=sort(agg$count_test,decreasing=TRUE)
X=agg[,2:22]
model2=randomForest(x=X[,2:21],y=log(X[,1]+1),ntree=2000)

agg$score2=predict(model,(X[,2:21]))
agg=agg[order(-agg$score2),] 
PAI=sum(agg$count_test[1:Ncell])/sum(agg$count_test)*PA/(Ncell*h*w)
print(PAI)
maxPAI=sum(sort_crime[1:Ncell])/sum(agg$count_test)*PA/(Ncell*h*w)
print(maxPAI)

spdf <- SpatialPolygonsDataFrame(grid_rotate, data=data.frame(id=row.names(grid_rotate), 
                                                       row.names=row.names(grid_rotate))) 
writeSpatialShape(x=spdf,paste0("SADM/validation/grid_",TARGET_CRIME))
write.csv(agg$grid[1:Ncell],paste0("SADM/validation/hotspot_id_",TARGET_CRIME,".csv"),row.names=F)

# create submission
agg<-setupfeatures(data,2017,TARGET_CRIME)

X=agg[,2:22]

agg$score3=predict(model2,(X[,2:21]))
agg=agg[order(-agg$score3),] 
PAI=sum(agg$count_test[1:Ncell])/sum(agg$count_test)*PA/(Ncell*h*w)
print(PAI)
sort_crime=sort(agg$count_test,decreasing=TRUE)
maxPAI=sum(sort_crime[1:Ncell])/sum(agg$count_test)*PA/(Ncell*h*w)
print(maxPAI)

spdf <- SpatialPolygonsDataFrame(grid_rotate, data=data.frame(id=row.names(grid_rotate), 
                                                              row.names=row.names(grid_rotate))) 

if(is.element(TARGET_CRIME,c("burglary","mvt"))){
  agg$grid=paste0("g",agg$grid)
}

writeSpatialShape(x=spdf,paste0("SADM/submission/grid_",TARGET_CRIME))
saveRDS(spdf, file = paste0("SADM/submission/spdf_",TARGET_CRIME,".rds"))
write.csv(agg[1:Ncell,c("grid","score3")],paste0("SADM/submission/hotspot_id_",TARGET_CRIME,".csv"),row.names=F)
write.csv(agg,paste0("SADM/submission/hotspot_score_",TARGET_CRIME,".csv"),row.names=F)




