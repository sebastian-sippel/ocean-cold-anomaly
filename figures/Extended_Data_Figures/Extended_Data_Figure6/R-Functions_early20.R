###read a table containing a time series, first column = years
read.ts<-function(filename,sep=";",header=T){
  ind<-as.matrix(read.table(filename,sep=sep,header=header))
  data<-ts(ind[,-1],start=ind[1,1])
  return(data)
}

###plot multiple ts in one plot
#x is cbind(ts1,ts2,...)
pl.mts<-function(x,cols=1:dim(x)[2],xl=tsp(x)[1:2],yl=c(min(x,na.rm=T),max(x,na.rm=T)),lwd=1){
  if(7 %in% cols) cols[which(cols==7)]<-"orange"
  if(is.ts(x)==T){
    plot(x,plot.type="s",col=cols,xlim=xl,ylim=yl,lwd=lwd)
  }else{
    plot(ts(x,start=1),plot.type="s",col=cols,xlim=xl,ylim=yl,lwd=lwd)
  }
  
}

#### zeitreihe skalieren sodass innerhalb einer bestimmten periode mean=0 und sd=1
scaletoperiod<-function(data,start,end){
  if(length(dim(data))>0){
    sds<-apply(window(data,start=start,end=end),2,sd,na.rm=T)
  }else{
    sds<-sd(window(data,start=start,end=end),na.rm=T)
  }
  scaled<-ts(scale(data,scale=sds,center=F),start=start(data)[1])
  if(length(dim(scaled))>0){
    ms<-apply(window(scaled,start=start,end=end),2,mean,na.rm=T)
  }else{
    ms<-mean(window(scaled,start=start,end=end),na.rm=T)
  }
  scaled<-scale(scaled,scale=F,center=ms)
  return(scaled)
}

#add density curve as line to plot.   x is density object 
densityline<-function(x,col=1){
  par(new=T)
  plot(x$x,x$y,col=col,t="l",tick=F,labels=F,bty="n",xlab="",ylab="")
}