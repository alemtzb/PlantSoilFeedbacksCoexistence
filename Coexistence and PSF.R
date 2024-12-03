rm(list=ls())
library(ggplot2)
#library(gridExtra)
library(gtools)
#library(plyr)
library(mgcv)
#load("~/betas supervivencia.Rdat")
setwd("~/Desktop/GitHub/PlantSoilFeedbacksCoexistence")
load("betas supervivencia.Rdat")
load("parametros.Rdata")
#this will load the required data from Martorell et al. 2021

# prepare data
b0p=rowMeans(b0)
for(i in 1:17) b0[,i]=b0p

#if short term effects are wanted, run:
matw=matrix(1,ncol=17,nrow=17)

#if long term effects are wanted, run:
matw=matrix(1,ncol=17,nrow=17)+wersp+wersp^2+wersp^3+wersp^4+wersp^5

# funtion to estimate a matrix of OSR for a given depth "prof" with a number of previous occupants "n" and long-term effects coefficient "ww"
calcmat=function(bersp,binter,prof,ww=matw,n=1) (bersp*n+binter*prof*n)*ww

###################################
# functions for pairwise analysis #
###################################

# estimate  strength of stabilizing and equalizing effects for a set of species "spid" at a given soil depth

stabiliz=function(prof,b1=bersp,b2=binter,ww=matw,spid=1:17){
	b1=b1[spid,spid]
	b2=b2[spid,spid]
	ww=ww[spid,spid]
	nsp=dim(b1)[1]
	mat=b1
	ersp=calcmat(b1,b2,prof, ww)
	for(i in 1:nsp){
		for(j in 1:nsp) mat[i,j]=-(ersp[i,i]-ersp[i,j]-ersp[j,i]+ersp[j,j])/2
	}
	mat
}

equaliz=function(prof,b1=bersp,b2=binter,ww=matw,spid=1:17){
	b1=b1[spid,spid]
	b2=b2[spid,spid]
	ww=ww[spid,spid]
	nsp=dim(b1)[1]
	mat=b1
	ersp=calcmat(b1,b2,prof,ww)
	for(i in 1:nsp){
		for(j in 1:nsp) mat[i,j]=(ersp[i,i]+ersp[j,i]-ersp[j,j]-ersp[i,j])/2
	}
	mat
}




# estimate stabilizing and equalizing effects for all pairs or species in spid and for
# all soil depths in vector "profs". Returns list of two arrays with stabilizing "ss" and
# equalizing "ee" effects.

pairprof=function(profs,b1=bersp,b2=binter,ww=matw,spid=1:17){
	b1=b1[spid,spid]
	b2=b2[spid,spid]
	ww=ww[spid,spid]
	nsp=dim(b1)[1]
	npf=length(profs)
	ss=array(dim=c(nsp,nsp,npf))
	ee=ss
	for(i in 1:npf){
		ss[,,i]=as.matrix(stabiliz(profs[i],b1,b2,ww))
		ee[,,i]=as.matrix(equaliz(profs[i],b1,b2,ww))
	}
	list("ss"=ss,"ee"=ee)
}


###########################################################
### How OCR relates to stability, fitness differences   ###
###  and their combined effect (Figure 1 in the ms)?    ###
###########################################################
ss8=stabiliz(8)
ee8=equaliz(8)
ss16=stabiliz(16)
ee16=equaliz(16)
ss24=stabiliz(24)
ee24=equaliz(24)

diag(ss8)=NA
diag(ss16)=NA
diag(ss24)=NA
diag(ee8)=NA
diag(ee16)=NA
diag(ee24)=NA


psf8=calcmat(bersp,binter,8)
psf16=calcmat(bersp,binter,16)
psf24=calcmat(bersp,binter,24)

diag(psf8)=NA
diag(psf16)=NA
diag(psf24)=NA

sqrttrans=function(x) sign(x)*sqrt(abs(x))

colu=colorRampPalette(list("orange","black"))(29)


par(mfrow=c(3,1))
#plot stabilization
xx8=sqrttrans(c(as.matrix(psf8[-18,-18])))
yy8=sqrttrans(c(as.matrix(ss8)))
xx16=sqrttrans(c(as.matrix(psf16[-18,-18])))
yy16=sqrttrans(c(as.matrix(ss16)))
xx24=sqrttrans(c(as.matrix(psf24[-18,-18])))
yy24=sqrttrans(c(as.matrix(ss24)))
mod8=gam(yy8~s(xx8))
xxx8=seq(min(xx8,na.rm=T),max(xx8,na.rm=T),length.out=25)
mod16=gam(yy16~s(xx16))
xxx16=seq(min(xx16,na.rm=T),max(xx16,na.rm=16),length.out=25)
mod24=gam(yy24~s(xx24))
xxx24=seq(min(xx24,na.rm=T),max(xx24,na.rm=T),length.out=25)

plot(xx8,yy8,xlab="occupancy-survival relationship",ylab="Stabilization",xaxt="n",xlim=c(min(c(xx8,xx16,xx24),na.rm=T),max(c(xx8,xx16,xx24),na.rm=T)),ylim=c(min(c(yy8,yy16,yy24),na.rm=T),max(c(yy8,yy16,yy24),na.rm=T)),col=colu[8])
axis(1,at=c(-sqrt(10),-sqrt(5),-1,0,1,sqrt(5),sqrt(10),sqrt(20),sqrt(30)),label=c(-10,-5,-1,0,1,5,10,20,30))
points(xx16,yy16,col=colu[16])
points(xx24,yy24,col=colu[24])

lines(xxx8,predict(mod8,data.frame(xx8=xxx8)),col="white",lwd=4)
lines(xxx16,predict(mod16,data.frame(xx16=xxx16)),col="white",lwd=4)
lines(xxx24,predict(mod24,data.frame(xx24=xxx24)),col="white",lwd=4)
lines(xxx8,predict(mod8,data.frame(xx8=xxx8)),col=colu[8])
lines(xxx16,predict(mod16,data.frame(xx16=xxx16)),col=colu[16])
lines(xxx24,predict(mod24,data.frame(xx24=xxx24)),col=colu[24])


#plot fitness differences
xx8=sqrttrans(c(as.matrix(psf8[-18,-18])))
yy8=abs(sqrttrans(c(as.matrix(ee8))))
xx16=sqrttrans(c(as.matrix(psf16[-18,-18])))
yy16=abs(sqrttrans(c(as.matrix(ee16))))
xx24=sqrttrans(c(as.matrix(psf24[-18,-18])))
yy24=abs(sqrttrans(c(as.matrix(ee24))))
mod8=gam(yy8~s(xx8))
xxx8=seq(min(xx8,na.rm=T),max(xx8,na.rm=T),length.out=25)
mod16=gam(yy16~s(xx16))
xxx16=seq(min(xx16,na.rm=T),max(xx16,na.rm=16),length.out=25)
mod24=gam(yy24~s(xx24))
xxx24=seq(min(xx24,na.rm=T),max(xx24,na.rm=T),length.out=25)

plot(xx8,yy8,xlab="occupancy-survival relationship",ylab="Fitness differences",xaxt="n",xlim=c(min(c(xx8,xx16,xx24),na.rm=T),max(c(xx8,xx16,xx24),na.rm=T)),ylim=c(min(c(yy8,yy16,yy24),na.rm=T),max(c(yy8,yy16,yy24),na.rm=T)),col=colu[8])
axis(1,at=c(-sqrt(10),-sqrt(5),-1,0,1,sqrt(5),sqrt(10),sqrt(20),sqrt(30)),label=c(-10,-5,-1,0,1,5,10,20,30))
points(xx16,yy16,col=colu[16])
points(xx24,yy24,col=colu[24])

lines(xxx8,predict(mod8,data.frame(xx8=xxx8)),col="white",lwd=4)
lines(xxx16,predict(mod16,data.frame(xx16=xxx16)),col="white",lwd=4)
lines(xxx24,predict(mod24,data.frame(xx24=xxx24)),col="white",lwd=4)
lines(xxx8,predict(mod8,data.frame(xx8=xxx8)),col=colu[8])
lines(xxx16,predict(mod16,data.frame(xx16=xxx16)),col=colu[16])
lines(xxx24,predict(mod24,data.frame(xx24=xxx24)),col=colu[24])


#plot net effects
xx8=sqrttrans(c(as.matrix(psf8[-18,-18])))
yy8=sqrttrans(c(as.matrix(ss8-abs(ee8))))
xx16=sqrttrans(c(as.matrix(psf16[-18,-18])))
yy16=sqrttrans(c(as.matrix(ss16-abs(ee16))))
xx24=sqrttrans(c(as.matrix(psf24[-18,-18])))
yy24=sqrttrans(c(as.matrix(ss24-abs(ee24))))
mod8=gam(yy8~s(xx8))
xxx8=seq(min(xx8,na.rm=T),max(xx8,na.rm=T),length.out=25)
mod16=gam(yy16~s(xx16))
xxx16=seq(min(xx16,na.rm=T),max(xx16,na.rm=16),length.out=25)
mod24=gam(yy24~s(xx24))
xxx24=seq(min(xx24,na.rm=T),max(xx24,na.rm=T),length.out=25)

plot(xx8,yy8,xlab="occupancy-survival relationship",ylab="Net effect",xaxt="n",xlim=c(min(c(xx8,xx16,xx24),na.rm=T),max(c(xx8,xx16,xx24),na.rm=T)),ylim=c(min(c(yy8,yy16,yy24),na.rm=T),max(c(yy8,yy16,yy24),na.rm=T)),col=colu[8])
axis(1,at=c(-sqrt(10),-sqrt(5),-1,0,1,sqrt(5),sqrt(10),sqrt(20),sqrt(30)),label=c(-10,-5,-1,0,1,5,10,20,30))
points(xx16,yy16,col=colu[16])
points(xx24,yy24,col=colu[24])

lines(xxx8,predict(mod8,data.frame(xx8=xxx8)),col="white",lwd=4)
lines(xxx16,predict(mod16,data.frame(xx16=xxx16)),col="white",lwd=4)
lines(xxx24,predict(mod24,data.frame(xx24=xxx24)),col="white",lwd=4)
lines(xxx8,predict(mod8,data.frame(xx8=xxx8)),col=colu[8])
lines(xxx16,predict(mod16,data.frame(xx16=xxx16)),col=colu[16])
lines(xxx24,predict(mod24,data.frame(xx24=xxx24)),col=colu[24])


###########################################################
### How depth relates to stability, fitness differences ###
###  and their combined effect (Figure 2 in the ms)?    ###
###########################################################


	profs=3:28
	todasprofs=pairprof(profs)
	tri=upper.tri(matrix(nrow=17,ncol=17))
	salss=matrix(nrow=5,ncol=26)
	salee=salss
	saldd=salss
	for(i in 1:26){
		ss=todasprofs$ss[,,i][tri]
		ee=abs(todasprofs$ee[,,i][tri])
		dd=ss-ee
		salss[,i]=quantile(ss,c(0.1,0.25,0.5,0.75,0.9))
		salee[,i]=quantile(ee,c(0.1,0.25,0.5,0.75,0.9))
		saldd[,i]=quantile(dd,c(0.1,0.25,0.5,0.75,0.9))
	}
	par(mfrow=c(3,1))
	plot(-1000,-1000,xlim=c(3,23),ylim=c(min(salss),max(salss)),ylab="stabilization")
	lines(c(-100,100),c(0,0),col="gray")
	lines(profs,salss[1,],lty=2)
	lines(profs,salss[2,],lty=1)
	lines(profs,salss[3,],lty=1,lwd=2)
	lines(profs,salss[4,],lty=1)
	lines(profs,salss[5,],lty=2)

	plot(-1000,-1000,xlim=c(3,23),ylim=c(0,max(salee)),ylab="fitness difference")
	lines(profs,salee[1,],lty=2)
	lines(profs,salee[2,],lty=1)
	lines(profs,salee[3,],lty=1,lwd=2)
	lines(profs,salee[4,],lty=1)
	lines(profs,salee[5,],lty=2)

	plot(-1000,-1000,xlim=c(3,23),ylim=c(min(saldd),max(saldd)),ylab="net effect")
	lines(c(-100,100),c(0,0),col="gray")
	lines(profs,saldd[1,],lty=2)
	lines(profs,saldd[2,],lty=1)
	lines(profs,saldd[3,],lty=1,lwd=2)
	lines(profs,saldd[4,],lty=1)
	lines(profs,saldd[5,],lty=2)




#######################################
# functions for multispecies analysis #
#######################################

# Calculate stability Ic and frequency of individuals per species
# (whose numbers are specified as a vector spid) in a given soil depth
Ic=function(prof,b1=bersp,b2=binter,ww=matw,spid=1:17){
	b1=b1[spid,spid]
	b2=b2[spid,spid]
	ww=ww[spid,spid]
	nsp=dim(b1)[1]
	sal=1:nsp
	names(sal)=names(b1)
	ersp=as.matrix(calcmat(b1,b2,prof,ww))
	dets=1:nsp
	for(i in 1:nsp){
		AA=ersp
		AA[,i]=1
		dets[i]=det(AA)
	}
	sal=dets/sum(dets)
	list("Ic"=(-1)^nsp*sum(dets),"feas"=sal)
}


#as Ic, but for several soil depths.
comprof=function(profs,b1=bersp,b2=binter,ww=matw,spid){
	b1=b1[spid,spid]
	b2=b2[spid,spid]
	ww=ww[spid,spid]
	nsp=dim(b1)[1]
	npf=length(profs)
	vIc=1:npf
	feas=matrix(nrow=nsp,ncol=npf)
	#row.names(feas)=names(b1)
	for(i in 1:npf){
		aa=Ic(profs[i],b1,b2,ww)
		vIc[i]=aa[[1]]
		feas[,i]=aa[[2]]
	}
	list("vIc"=vIc,"feas"=feas)
}

#Estimates Ic for all possible communities from those having just two species
#to that havinag all species. There are 131054 species so this takes a while.
#For avery community this returns a vector "species" containing the species "especies"
#in the community, the "Ic" value, the frequency of all species in the stable community "feas"
#and the pairwise stabilization and equalization metrics.

todascom=function(prof,b1=bersp,b2=binter,ww=matw){
	sal=list()
	for(com in 2:17){
		elems=combinations(17,com)
		for(i in 1:dim(elems)[1]){
			Ictemp=Ic(prof,spid=elems[i,])
			sal[[(length(sal) +1)]] = list("especies"=elems[i,],"Ic"=Ictemp$Ic,"feas"=Ictemp$feas)
		}
	}
	sal
}


###obtain statistics for all possible communities for all soil depths between 4
# and 28 cm.

################IMPORTANT################
#In what follows both cummulative (objects ending in l) or short term
# (ending in s). This convention is used throughout the code.

#We have #ed the following instructions because they consume a lot of time and
#runing them accidentally may freeze your computer. However, they need to be
#executed to produce the basic objects that will be needed below.
#matw=matrix(1,ncol=17,nrow=17)+wersp+wersp^2+wersp^3+wersp^4+wersp^5
#com4l=todascom(4)
#com8l=todascom(5)
#com4l=todascom(6)
#com8l=todascom(7)
#com4l=todascom(8)
#com8l=todascom(9)
#com10l=todascom(10)
#com11l=todascom(11)
#com12l=todascom(12)
#com13l=todascom(13)
#com14l=todascom(14)
#com15l=todascom(15)
#com16l=todascom(16)
#com17l=todascom(17)
#com18l=todascom(18)
#com19l=todascom(19)
#com20l=todascom(20)
#com21l=todascom(21)
#com22l=todascom(22)
#com23l=todascom(23)
#com24l=todascom(24)
#com25l=todascom(25)
#com26l=todascom(26)
#com27l=todascom(27)
#com28l=todascom(28)

#matw=matrix(1,ncol=17,nrow=17)
#com4s=todascom(4)
#com8s=todascom(5)
#com4s=todascom(6)
#com8s=todascom(7)
#com4s=todascom(8)
#com8s=todascom(9)
#com10s=todascom(10)
#com11s=todascom(11)
#com12s=todascom(12)
#com13s=todascom(13)
#com14s=todascom(14)
#com15s=todascom(15)
#com16s=todascom(16)
#com17s=todascom(17)
#com18s=todascom(18)
#com19s=todascom(19)
#com20s=todascom(20)
#com21s=todascom(21)
#com22s=todascom(22)
#com23s=todascom(23)
#com24s=todascom(24)
#com25s=todascom(25)
#com26s=todascom(26)
#com27s=todascom(27)
#com28s=todascom(28)



#find feasible subcomunities in an object "obj" produced by todascom
findfeas=function(obj) which(sapply(sapply(obj,"[[","feas"),min)>0)
#find stable subcomunities in an object produced by todascom
findestab=function(obj) which(sapply(obj,"[[","Ic")<0)
#find feasible and setable subcomunities in an object produced by todascom
findfe=function(obj) intersect(findfeas(obj),findestab(obj))



#remove all communites that are not SFC from the "com" objects generated before
com4l=com4l[findfe(com4l)]
com5l=com5l[findfe(com5l)]
com6l=com6l[findfe(com6l)]
com7l=com7l[findfe(com7l)]
com8l=com8l[findfe(com8l)]
com9l=com9l[findfe(com9l)]
com10l=com10l[findfe(com10l)]
com11l=com11l[findfe(com11l)]
com12l=com12l[findfe(com12l)]
com13l=com13l[findfe(com13l)]
com14l=com14l[findfe(com14l)]
com15l=com15l[findfe(com15l)]
com16l=com16l[findfe(com16l)]
com17l=com17l[findfe(com17l)]
com18l=com18l[findfe(com18l)]
com19l=com19l[findfe(com19l)]
com20l=com20l[findfe(com20l)]
com21l=com21l[findfe(com21l)]
com22l=com22l[findfe(com22l)]
com23l=com23l[findfe(com23l)]
com24l=com24l[findfe(com24l)]
com25l=com25l[findfe(com25l)]
com26l=com26l[findfe(com26l)]
com27l=com27l[findfe(com27l)]
com28l=com28l[findfe(com28l)]

com4s=com4s[findfe(com4s)]
com5s=com5s[findfe(com5s)]
com6s=com6s[findfe(com6s)]
com7s=com7s[findfe(com7s)]
com8s=com8s[findfe(com8s)]
com9s=com9s[findfe(com9s)]
com10s=com10s[findfe(com10s)]
com11s=com11s[findfe(com11s)]
com12s=com12s[findfe(com12s)]
com13s=com13s[findfe(com13s)]
com14s=com14s[findfe(com14s)]
com15s=com15s[findfe(com15s)]
com16s=com16s[findfe(com16s)]
com17s=com17s[findfe(com17s)]
com18s=com18s[findfe(com18s)]
com19s=com19s[findfe(com19s)]
com20s=com20s[findfe(com20s)]
com21s=com21s[findfe(com21s)]
com22s=com22s[findfe(com22s)]
com23s=com23s[findfe(com23s)]
com24s=com24s[findfe(com24s)]
com25s=com25s[findfe(com25s)]
com26s=com26s[findfe(com26s)]
com27s=com27s[findfe(com27s)]
com28s=com28s[findfe(com28s)]


####################################################
# fraction of all possible comunities with a given #
#     richness that are stable and feasible        #
#          and code for figure 3                   #
####################################################


# use an object derived from todascom from which communities that are not SFC
#have been removed to calculate the fraction of SFC
fracrich=function(obj){
	univers=2:17
	for(i in 2:17) univers[i-1]=dim(combinations(17,i))[1]
	aa=sapply(sapply(obj,"[[","especies"),length)
	hist(aa,breaks=1:17,plot=F)$counts/univers
}


#functions for plotting
plotpyram=function(fracs,prof,lims,med=T){
	dat=fracs[,prof/4]
	plot(-100,-100,ylim=c(1,11.5),xlim=c(-lims/2,lims/2),xaxt="n")
	for(i in 2:11){
		rect(-dat[i-1]/2,i-0.5,dat[i-1]/2,i+0.5,col=1)
	}
	if(med==T) points(0,medgrup(dat)+1,col="white",cex=2,pch=19)
}

todpyra=function(fracs,med){
	par(mfrow=c(1,8))
	lims=max(fracs)
	for(i in seq(4,28,4)) plotpyram(fracs,i,lims,med)
	plot(-100,-100,ylim=c(0,12),xlim=c(-lims/2,lims/2))
	rect(-0.05,2,0.05,2.2,col=1)
}

#Create plot for long term effects; change accordingly for short term.
#Note that the final panel contains a bar that indicates scale.

fracciones=matrix(nrow=16,ncol=7)
fracciones[,1]=fracrich(com4l)
fracciones[,2]=fracrich(com8l)
fracciones[,3]=fracrich(com12l)
fracciones[,4]=fracrich(com16l)
fracciones[,5]=fracrich(com20l)
fracciones[,6]=fracrich(com24l)
fracciones[,7]=fracrich(com28l)

todpyra(fracciones,F)

############################################################
#### Relationship between OSR and the richness in which ####
####         species tend to occur (Figure 4)           ####
############################################################



#eliminate from an object generated by todascom the subcomunities with n species
#that are subsets of others with N species.
#The use of this funtion only makes sense for the todascom object AFTER removing
#all communities that are not stable and feasible
killsubsets=function(n,N,obj){
	SFC=findfe(obj)
	obj=obj[SFC]
	rich=sapply(sapply(obj,"[[","especies"),length)
	whichn=obj[which(rich==n)]
	numnini=length(whichn)

	while(N>n&&length(whichn)>0){
		whichN=obj[which(rich==N)]
		numN=length(which(rich==N))
		numn=length(whichn)
		dataN=matrix(nrow=numN,ncol=N)
		for(i in 1:numN){
			dataN[i,]=whichN[[i]][[1]]
		}
		for(j in numn:1){
			compo=whichn[[j]][[1]]
			indic=FALSE
			i=1
			while(indic==F&&i<=numN){
				if(all(compo %in% dataN[i,])) {
					indic=TRUE
					whichn=whichn[-j]
					i=1
					} else {
						i=i+1
					}
			}
		}
		N=N-1
	}
	for(i in 1:numnini) obj=obj[-1]
	if(length(whichn)==0) obj
	else append(obj,whichn)
}

#As killsubsets, but will do for all possible n. maxim is the number of
#species of the richest stable and feasible community in the set.
killall=function(obj,maxim){
	for(i in 2:16){
		obj=killsubsets(i,maxim,obj)
	}
	obj
}


#creates a presence-absence matrix for all SFC that are not subsets of others
creamatcompa=function(objsub){
	comN=length(objsub)
	rich=sapply(sapply(objsub,"[[","especies"),length)
	sal=matrix(0,nrow=comN,ncol=17)
	for(i in 1:comN){
		sal[i,objsub[[i]][[1]]]=1
	}
	sal
}



#Use previous functions to obtain the presence absence matrices that we will require
subcom4l=killall(com4l,10)
subcom5l=killall(com5l,10)
subcom6l=killall(com6l,10)
subcom7l=killall(com7l,10)
subcom8l=killall(com8l,10)
subcom9l=killall(com9l,10)
subcom10l=killall(com10l,10)
subcom11l=killall(com11l,11)
subcom12l=killall(com12l,11)
subcom13l=killall(com13l,11)
subcom14l=killall(com14l,11)
subcom15l=killall(com15l,11)
subcom16l=killall(com16l,10)
subcom17l=killall(com17l,11)
subcom18l=killall(com18l,10)
subcom19l=killall(com19l,12)
subcom20l=killall(com20l,11)
subcom21l=killall(com21l,12)
subcom22l=killall(com22l,12)
subcom23l=killall(com23l,10)
subcom24l=killall(com24l,11)
subcom25l=killall(com25l,9)
subcom26l=killall(com26l,9)
subcom27l=killall(com27l,9)
subcom28l=killall(com28l,9)

subcom4s=killall(com4s,10)
subcom5s=killall(com5s,10)
subcom6s=killall(com6s,11)
subcom7s=killall(com7s,11)
subcom8s=killall(com8s,10)
subcom9s=killall(com9s,10)
subcom10s=killall(com10s,11)
subcom11s=killall(com11s,11)
subcom12s=killall(com12s,11)
subcom13s=killall(com13s,11)
subcom14s=killall(com14s,11)
subcom15s=killall(com15s,10)
subcom16s=killall(com16s,10)
subcom17s=killall(com17s,11)
subcom18s=killall(com18s,10)
subcom19s=killall(com19s,11)
subcom20s=killall(com20s,11)
subcom21s=killall(com21s,11)
subcom22s=killall(com22s,11)
subcom23s=killall(com23s,10)
subcom24s=killall(com24s,11)
subcom25s=killall(com25s,11)
subcom26s=killall(com26s,10)
subcom27s=killall(com27s,10)
subcom28s=killall(com28s,10)

matpa4l=creamatcompa(subcom4l)
matpa5l=creamatcompa(subcom5l)
matpa6l=creamatcompa(subcom6l)
matpa7l=creamatcompa(subcom7l)
matpa8l=creamatcompa(subcom8l)
matpa9l=creamatcompa(subcom9l)
matpa10l=creamatcompa(subcom10l)
matpa11l=creamatcompa(subcom11l)
matpa12l=creamatcompa(subcom12l)
matpa13l=creamatcompa(subcom13l)
matpa14l=creamatcompa(subcom14l)
matpa15l=creamatcompa(subcom15l)
matpa16l=creamatcompa(subcom16l)
matpa17l=creamatcompa(subcom17l)
matpa18l=creamatcompa(subcom18l)
matpa19l=creamatcompa(subcom19l)
matpa20l=creamatcompa(subcom20l)
matpa21l=creamatcompa(subcom21l)
matpa22l=creamatcompa(subcom22l)
matpa23l=creamatcompa(subcom23l)
matpa24l=creamatcompa(subcom24l)
matpa25l=creamatcompa(subcom25l)
matpa26l=creamatcompa(subcom26l)
matpa27l=creamatcompa(subcom27l)
matpa28l=creamatcompa(subcom28l)

matpa4s=creamatcompa(subcom4s)
matpa5s=creamatcompa(subcom5s)
matpa6s=creamatcompa(subcom6s)
matpa7s=creamatcompa(subcom7s)
matpa8s=creamatcompa(subcom8s)
matpa9s=creamatcompa(subcom9s)
matpa10s=creamatcompa(subcom10s)
matpa11s=creamatcompa(subcom11s)
matpa12s=creamatcompa(subcom12s)
matpa13s=creamatcompa(subcom13s)
matpa14s=creamatcompa(subcom14s)
matpa15s=creamatcompa(subcom15s)
matpa16s=creamatcompa(subcom16s)
matpa17s=creamatcompa(subcom17s)
matpa18s=creamatcompa(subcom18s)
matpa19s=creamatcompa(subcom19s)
matpa20s=creamatcompa(subcom20s)
matpa21s=creamatcompa(subcom21s)
matpa22s=creamatcompa(subcom22s)
matpa23s=creamatcompa(subcom23s)
matpa24s=creamatcompa(subcom24s)
matpa25s=creamatcompa(subcom25s)
matpa26s=creamatcompa(subcom26s)
matpa27s=creamatcompa(subcom27s)
matpa28s=creamatcompa(subcom28s)

#use creamatcompa object "objmat" to plot changes in the "selectivity" of species for
#communities with different richness, and return the correlations bewteen richness
# and selectivity. Using the slopes of the respective regression results in similar
# results as those presented in the main text.
plotselrich=function(objmat){
	comN=dim(objmat)[1]
	rich=rowSums(objmat)
	cuentarich=table(rich)
	riches=as.numeric(names(cuentarich))
	selecti=matrix(0,nrow=length(cuentarich),ncol=17)
	pends=1:length(cuentarich)
	esper=cuentarich*riches/17
	for(i in 1:17){
		obser=tapply(objmat[,i],rowSums(objmat),sum)
		selecti[,i]=obser/esper
		pends[i]=cor(selecti[,i],riches)
	}
	par(mfrow=c(6,3),mai = c(0.1, 0.1, 0.1, 0.1))
	for(i in 1:17){
		plot(c(1,18),c(1,1), col="gray",type="l",ylim=c(min(selecti)-.1,max(selecti)+.1),xlim=(c(min(rich),max(rich))),xaxt="n",xlab="",yaxt="n",ylab="")
		lines(riches,selecti[,i])
		text(x=min(riches)+1,y=max(selecti),names(bersp)[i])
	}
	pends
}

#check whether species that experience more positive or negative PSF occur in communities with
#different richness
#The following is for cummulative OSR, make the appropriate changes for short term effects.
#don't forget to set matw to the respective values

sel4=plotselrich(matpa4l)
sel8=plotselrich(matpa8l)
sel12=plotselrich(matpa12l)
sel16=plotselrich(matpa16l)
sel20=plotselrich(matpa20l)
sel24=plotselrich(matpa24l)
sel28=plotselrich(matpa28l)
sel5=plotselrich(matpa5l)
sel9=plotselrich(matpa9l)
sel13=plotselrich(matpa13l)
sel17=plotselrich(matpa17l)
sel21=plotselrich(matpa21l)
sel25=plotselrich(matpa25l)
sel6=plotselrich(matpa6l)
sel10=plotselrich(matpa10l)
sel14=plotselrich(matpa14l)
sel18=plotselrich(matpa18l)
sel22=plotselrich(matpa22l)
sel26=plotselrich(matpa26l)
sel7=plotselrich(matpa7l)
sel11=plotselrich(matpa11l)
sel15=plotselrich(matpa15l)
sel19=plotselrich(matpa19l)
sel23=plotselrich(matpa23l)
sel27=plotselrich(matpa27l)


psf4=calcmat(bersp,binter,4)
psf8=calcmat(bersp,binter,8)
psf12=calcmat(bersp,binter,12)
psf16=calcmat(bersp,binter,16)
psf20=calcmat(bersp,binter,20)
psf24=calcmat(bersp,binter,24)
psf5=calcmat(bersp,binter,5)
psf9=calcmat(bersp,binter,9)
psf13=calcmat(bersp,binter,13)
psf17=calcmat(bersp,binter,17)
psf21=calcmat(bersp,binter,21)
psf25=calcmat(bersp,binter,25)
psf6=calcmat(bersp,binter,6)
psf10=calcmat(bersp,binter,10)
psf14=calcmat(bersp,binter,14)
psf18=calcmat(bersp,binter,18)
psf22=calcmat(bersp,binter,22)
psf26=calcmat(bersp,binter,26)
psf7=calcmat(bersp,binter,7)
psf11=calcmat(bersp,binter,11)
psf15=calcmat(bersp,binter,15)
psf19=calcmat(bersp,binter,19)
psf23=calcmat(bersp,binter,23)
psf27=calcmat(bersp,binter,27)
psf28=calcmat(bersp,binter,28)


diag(psf4)=NA
diag(psf8)=NA
diag(psf12)=NA
diag(psf16)=NA
diag(psf20)=NA
diag(psf24)=NA
diag(psf28)=NA
diag(psf5)=NA
diag(psf9)=NA
diag(psf13)=NA
diag(psf17)=NA
diag(psf21)=NA
diag(psf25)=NA
diag(psf6)=NA
diag(psf10)=NA
diag(psf14)=NA
diag(psf18)=NA
diag(psf22)=NA
diag(psf26)=NA
diag(psf7)=NA
diag(psf11)=NA
diag(psf15)=NA
diag(psf19)=NA
diag(psf23)=NA
diag(psf27)=NA

# plot using mean OSR of the effects that each species receives from others
xx=c(-100,100)
colu=colorRampPalette(list("orange","black"))(29)
plot(-100,-100,ylim=c(-1,1),xlim=c(-.5,.2),xlab="mean OSR",ylab="richness preference")
for(i in 4:28){
	mod=lm(get(paste("sel",i,sep=""))~rowMeans(get(paste("psf",i,sep="")),na.rm=T))$coefficients
	lines(xx,mod[1]+mod[2]*xx,col=colu[i])
}
ampl=0.01
for(i in 4:28) lines(c(-1+ampl*i,-1+ampl*i),c(-.75,-.65),lwd=4,col=colu[i])
text(-.96,-.82,4)
text(-.72,-.82,28)

# plot using mean OSR of the effects that each species exerts on others
plot(-100,-100,ylim=c(-1,1),xlim=c(-.5,.2),xlab="mean OSR",ylab="richness preference")
for(i in 4:28){
	mod=lm(get(paste("sel",i,sep=""))~colMeans(get(paste("psf",i,sep="")),na.rm=T))$coefficients
	lines(xx,mod[1]+mod[2]*xx,col=colu[i])
}



##############################################################
### Changes in SFC over the depth gradient   (Figure 5) ######
##############################################################


#as creamatcompa, but the output has relative abundances instead of presences
creamatcom=function(objsub){
	comN=length(objsub)
	rich=sapply(sapply(objsub,"[[","especies"),length)
	sal=matrix(0,nrow=comN,ncol=17)
	for(i in 1:comN){
		sal[i,objsub[[i]][[1]]]=objsub[[i]][[3]]
	}
	sal
}

mat4l=creamatcom(subcom4l)
mat5l=creamatcom(subcom5l)
mat6l=creamatcom(subcom6l)
mat7l=creamatcom(subcom7l)
mat8l=creamatcom(subcom8l)
mat9l=creamatcom(subcom9l)
mat10l=creamatcom(subcom10l)
mat11l=creamatcom(subcom11l)
mat12l=creamatcom(subcom12l)
mat13l=creamatcom(subcom13l)
mat14l=creamatcom(subcom14l)
mat15l=creamatcom(subcom15l)
mat16l=creamatcom(subcom16l)
mat17l=creamatcom(subcom17l)
mat18l=creamatcom(subcom18l)
mat19l=creamatcom(subcom19l)
mat20l=creamatcom(subcom20l)
mat21l=creamatcom(subcom21l)
mat22l=creamatcom(subcom22l)
mat23l=creamatcom(subcom23l)
mat24l=creamatcom(subcom24l)
mat25l=creamatcom(subcom25l)
mat26l=creamatcom(subcom26l)
mat27l=creamatcom(subcom27l)
mat28l=creamatcom(subcom28l)

mat4s=creamatcom(subcom4s)
mat5s=creamatcom(subcom5s)
mat6s=creamatcom(subcom6s)
mat7s=creamatcom(subcom7s)
mat8s=creamatcom(subcom8s)
mat9s=creamatcom(subcom9s)
mat10s=creamatcom(subcom10s)
mat11s=creamatcom(subcom11s)
mat12s=creamatcom(subcom12s)
mat13s=creamatcom(subcom13s)
mat14s=creamatcom(subcom14s)
mat15s=creamatcom(subcom15s)
mat16s=creamatcom(subcom16s)
mat17s=creamatcom(subcom17s)
mat18s=creamatcom(subcom18s)
mat19s=creamatcom(subcom19s)
mat20s=creamatcom(subcom20s)
mat21s=creamatcom(subcom21s)
mat22s=creamatcom(subcom22s)
mat23s=creamatcom(subcom23s)
mat24s=creamatcom(subcom24s)
mat25s=creamatcom(subcom25s)
mat26s=creamatcom(subcom26s)
mat27s=creamatcom(subcom27s)
mat28s=creamatcom(subcom28s)

##add information to creamatcom objects "matcom" for subsequent analysis.
#Note that matcom has species abundances in columns 1 to 17
#in column 18 we add a unique identifier that is the same  for a community
#that is repeated in different soil depths, and in column 19 the depth
#in which the SFC occurs
matid=function(matcom,prof){
	nr=dim(matcom)[1]
	sal=matrix(nrow=nr,ncol=19)
	sal[,1:17]=matcom
	bb=2^(0:16)
	for(i in 1:nr) sal[i,18]=sum(matcom[i,]*bb)
	sal[,19]=prof
	sal
}


#paste all abundances matrices with their depths and identifiers into a single matrix.


supercoml=rbind(matid(mat4l,4),matid(mat5l,5),matid(mat6l,6),matid(mat7l,7),matid(mat8l,8),matid(mat9l,9),matid(mat10l,10),matid(mat11l,11),matid(mat12l,12),matid(mat13l,13),matid(mat14l,14),matid(mat15l,15),matid(mat16l,16),matid(mat17l,17),matid(mat18l,18),matid(mat19l,19),matid(mat20l,20),matid(mat21l,21),matid(mat22l,22),matid(mat23l,23),matid(mat24l,24),matid(mat25l,25),matid(mat26l,26),matid(mat27l,27),matid(mat28l,28))
supercoms=rbind(matid(mat4s,4),matid(mat5s,5),matid(mat6s,6),matid(mat7s,7),matid(mat8s,8),matid(mat9s,9),matid(mat10s,10),matid(mat11s,11),matid(mat12s,12),matid(mat13s,13),matid(mat14s,14),matid(mat15s,15),matid(mat16s,16),matid(mat17s,17),matid(mat18s,18),matid(mat19s,19),matid(mat20s,20),matid(mat21s,21),matid(mat22s,22),matid(mat23s,23),matid(mat24s,24),matid(mat25s,25),matid(mat26s,26),matid(mat27s,27),matid(mat28s,28))

#remove communities that are repeated in different soil depths
depurar=function(multid){
	multid=multid[order(multid[,18],decreasing=F),]
	for(i in (dim(multid)[1]-1):1){
		if(multid[i,18]==multid[i+1,18]) multid=multid[-(i+1),]
	}
	multid
}

#Conduct a PCA using the matrix from which repeated communities have been removed,
#and scale scores so that they are between 0 and 1 (This is required to color the plot).
#then assing the scores to the communities in supercom based on their identifiers.
ordenar=function(multid,objdep){
	pc=princomp(objdep[,1:17])
	objdep[,19]=(pc$scores[,1]-min(pc$scores[,1]))/max(pc$scores[,1]-min(pc$scores[,1]))
	ncom=dim(multid)[1]
	multid=cbind(multid,1:ncom)
	for(i in 1:dim(multid)[1]) multid[i,20]=objdep[which(objdep[,18]==multid[i,18]),19]
	cbind(multid,rowSums(multid[,1:17]))
}



#produce a stacked area plot of an "ordenar"object
graford=function(objgr){
	dat=data.frame("id"=as.numeric(as.factor(objgr[,20])),"prof"=objgr[,19])
	comsporprof=hist(dat$prof,breaks=3:28,plot=F)$counts
	frec=1:dim(dat)[1]
	for(i in 1:dim(dat)[1]) frec[i]=1/comsporprof[dat$prof[i]-3]
	dat=cbind(dat,frec)
	maxid=max(dat$id)
	for(i in unique(dat$prof)){
		tem=data.frame("id"=1:maxid,"prof"=rep(i,maxid),"frec"=rep(0,maxid))
		cual=dat[which(dat$prof==i),1]
		tem=tem[-cual,]
		dat=rbind(dat,tem)
	}
	scor=sort(unique(objgr[,20]))
	gamac=colorRamp(c("blue","white","red"))
	precol=gamac(scor)
	colores=rgb(precol[,1],precol[,2],precol[,3],maxColorValue=255)
	plox = ggplot(dat, aes(x=prof, y=frec, fill=as.factor(id)))
	plox + geom_area(show.legend = FALSE)+scale_fill_manual(values=colores)
}

#prepare data for a multi-panel stacked area plot
grafordb=function(objgr){

	dat=data.frame("id"=as.numeric(as.factor(objgr[,20])),"prof"=objgr[,19])
	comsporprof=hist(dat$prof,breaks=3:28,plot=F)$counts
	frec=1:dim(dat)[1]
	for(i in 1:dim(dat)[1]) frec[i]=1/comsporprof[dat$prof[i]-3]
	dat=cbind(dat,frec)
	maxid=max(dat$id)
	for(i in unique(dat$prof)){
		tem=data.frame("id"=1:maxid,"prof"=rep(i,maxid),"frec"=rep(0,maxid))
		cual=dat[which(dat$prof==i),1]
		tem=tem[-cual,]
		dat=rbind(dat,tem)
	}
	scor=sort(unique(objgr[,20]))
	gamac=colorRamp(c("blue","white","red"))
	precol=gamac(scor)
	colores=rgb(precol[,1],precol[,2],precol[,3],maxColorValue=255)
	list(ggplot(dat, aes(x=prof, y=frec, fill=as.factor(id))),colores)
}

#as graford, but produces graphs for communities having a richness of 5,6,
#7, 8 or 9 or more species (from top to bottom)
graford2=function(objgr){
	par(mfrow=c(5,1),mai=rep(0.05,4))
	plox5=grafordb(objgr[which(paragraf[,21]==5),])
	plox6=grafordb(objgr[which(paragraf[,21]==6),])
	plox7=grafordb(objgr[which(paragraf[,21]==7),])
	plox8=grafordb(objgr[which(paragraf[,21]==8),])
	plox9=grafordb(objgr[which(paragraf[,21]>8),])
	grid.arrange(plox5[[1]]+ geom_area(show.legend = FALSE)+scale_fill_manual(values=plox5[[2]]),
		plox6[[1]]+ geom_area(show.legend = FALSE)+scale_fill_manual(values=plox6[[2]]),
		plox7[[1]]+ geom_area(show.legend = FALSE)+scale_fill_manual(values=plox7[[2]]),
		plox8[[1]]+ geom_area(show.legend = FALSE)+scale_fill_manual(values=plox8[[2]]),
		plox9[[1]]+ geom_area(show.legend = FALSE)+scale_fill_manual(values=plox9[[2]]),
		ncol=1, nrow = 5)
}

depcoml=depurar(supercoml)
depcoms=depurar(supercoms)
paragraf=ordenar(supercoml,depcoml)

graford(paragraf)



