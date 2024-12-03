# prepare data
b0p=rowMeans(b0)
for(i in 1:17) b0[,i]=b0p

#if short term effects are wanted, run:
matw=matrix(1,ncol=17,nrow=17)

#if long term effects are wanted, run:
#matw=matrix(1,ncol=17,nrow=17)+wersp+wersp^2+wersp^3+wersp^4+wersp^5

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


# Set layout
layout(matrix(1:9, nrow = 3, byrow = FALSE), widths = c(1, 1, 1), heights = c(1.5, 1.5, 1.5))

# Adjust margins
par(mar = c(4, 4, 0, 0), oma = c(0.25, 0.25, 0.25, 0.25))

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

plot(xx8,yy8,
     ylab="Stabilization",
     xaxt="n",
     xlab=NA,
     xlim=c(min(c(xx8,xx16,xx24),na.rm=T),max(c(xx8,xx16,xx24),na.rm=T)),
     ylim=c(min(c(yy8,yy16,yy24),na.rm=T),max(c(yy8,yy16,yy24),na.rm=T)),
     col=colu[8])
axis(1,at=c(-sqrt(10),-sqrt(5),-1,0,1,sqrt(5),sqrt(10),sqrt(20),sqrt(30)),label=FALSE)
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


plot(xx8,yy8,
     xlab=NA,
     ylab="Fitness differences",
     xaxt="n",
     xlim=c(min(c(xx8,xx16,xx24),na.rm=T),max(c(xx8,xx16,xx24),na.rm=T)),
     ylim=c(min(c(yy8,yy16,yy24),na.rm=T),max(c(yy8,yy16,yy24),na.rm=T)),
     col=colu[8])
axis(1,at=c(-sqrt(10),-sqrt(5),-1,0,1,sqrt(5),sqrt(10),sqrt(20),sqrt(30)),label=FALSE)
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

plot(xx8,yy8,
     xlab="occupancy-survival relationship",
     ylab="Net effect",
     xaxt="n",
     xlim=c(min(c(xx8,xx16,xx24),na.rm=T),max(c(xx8,xx16,xx24),na.rm=T)),
     ylim=c(min(c(yy8,yy16,yy24),na.rm=T),max(c(yy8,yy16,yy24),na.rm=T)),
     col=colu[8])
axis(1,at=c(-sqrt(10),-sqrt(5),-1,0,1,sqrt(5),sqrt(10),sqrt(20),sqrt(30)),label=c(-10,-5,-1,0,1,5,10,20,30))
points(xx16,yy16,col=colu[16])
points(xx24,yy24,col=colu[24])

lines(xxx8,predict(mod8,data.frame(xx8=xxx8)),col="white",lwd=4)
lines(xxx16,predict(mod16,data.frame(xx16=xxx16)),col="white",lwd=4)
lines(xxx24,predict(mod24,data.frame(xx24=xxx24)),col="white",lwd=4)
lines(xxx8,predict(mod8,data.frame(xx8=xxx8)),col=colu[8])
lines(xxx16,predict(mod16,data.frame(xx16=xxx16)),col=colu[16])
lines(xxx24,predict(mod24,data.frame(xx24=xxx24)),col=colu[24])

####################


#if long term effects are wanted, run:
matw=matrix(1,ncol=17,nrow=17)+wersp+wersp^2+wersp^3+wersp^4+wersp^5

# funtion to estimate a matrix of OSR for a given depth "prof" with a number of previous occupants "n" and long-term effects coefficient "ww"
calcmat=function(bersp,binter,prof,ww=matw,n=1) (bersp*n+binter*prof*n)*ww

###################################
# functions for pairwise analysis #
###################################


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


#par(mfrow=c(3,1))
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

plot(xx8,yy8,
     xlab=NA,
     ylab="",
     xaxt="n",
     xlim=c(min(c(xx8,xx16,xx24),na.rm=T),max(c(xx8,xx16,xx24),na.rm=T)),
     ylim=c(min(c(yy8,yy16,yy24),na.rm=T),max(c(yy8,yy16,yy24),na.rm=T)),
     col=colu[8])
axis(1,at=c(-sqrt(10),-sqrt(5),-1,0,1,sqrt(5),sqrt(10),sqrt(20),sqrt(30)),label=FALSE)
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

plot(xx8,yy8,
     xlab=NA,
     ylab=NA,
     xaxt="n",
     xlim=c(min(c(xx8,xx16,xx24),na.rm=T),max(c(xx8,xx16,xx24),na.rm=T)),
     ylim=c(min(c(yy8,yy16,yy24),na.rm=T),max(c(yy8,yy16,yy24),na.rm=T)),
     col=colu[8])
axis(1,at=c(-sqrt(10),-sqrt(5),-1,0,1,sqrt(5),sqrt(10),sqrt(20),sqrt(30)),label=FALSE)
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

plot(xx8,yy8,
     xlab="occupancy-survival relationship",
     ylab=NA,
     xaxt="n",
     xlim=c(min(c(xx8,xx16,xx24),na.rm=T),max(c(xx8,xx16,xx24),na.rm=T)),
     ylim=c(min(c(yy8,yy16,yy24),na.rm=T),max(c(yy8,yy16,yy24),na.rm=T)),
     col=colu[8])
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
#par(mfrow=c(3,1))
plot(-1000,-1000,xlim=c(3,23),ylim=c(min(salss),max(salss)),xlab=NA,ylab=NA,xaxt="n")
axis(1, at = axTicks(1), labels = FALSE)
lines(c(-100,100),c(0,0),col="gray")
lines(profs,salss[1,],lty=2)
lines(profs,salss[2,],lty=1)
lines(profs,salss[3,],lty=1,lwd=2)
lines(profs,salss[4,],lty=1)
lines(profs,salss[5,],lty=2)

plot(-1000,-1000,xlim=c(3,23),ylim=c(0,max(salee)),xlab=NA,ylab=NA,xaxt="n")
lines(profs,salee[1,],lty=2)
lines(profs,salee[2,],lty=1)
lines(profs,salee[3,],lty=1,lwd=2)
lines(profs,salee[4,],lty=1)
lines(profs,salee[5,],lty=2)

plot(-1000,-1000,xlim=c(3,23),ylim=c(min(saldd),max(saldd)),xlab="Soil Depth (cm)",ylab=NA)
lines(c(-100,100),c(0,0),col="gray")
lines(profs,saldd[1,],lty=2)
lines(profs,saldd[2,],lty=1)
lines(profs,saldd[3,],lty=1,lwd=2)
lines(profs,saldd[4,],lty=1)
lines(profs,saldd[5,],lty=2)




