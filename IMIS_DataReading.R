##_______________________________________________________________
##|                                                             |  
##| Data reading                                                |
##|_____________________________________________________________|
input.data <- read.table(file=file.in, 
	col.names=list(
	"Sample","Depth", "Depth.err","Be10", "Be10.err", "Ne21", "Ne21.err","Al26", "Al26.err"
  ), row.names = NULL, skip=1)
if(sum(input.data$Be10)!=0 && is.numeric(sum(input.data$Be10))){Be10.flag='y'
           input.data$Be10[which(input.data$Be10==0)] <- NA; input.data$Be10.err[which(input.data$Be10==0)] <- NA
           } else {Be10.flag='n'}
if(sum(input.data$Ne21)!=0 && is.numeric(sum(input.data$Ne21))){Ne21.flag='y' 
           input.data$Ne21[which(input.data$Ne21==0)] <- NA; input.data$Ne21.err[which(input.data$Ne21==0)] <- NA
           lim.x.Ne[1]<-min(lim.x.Ne[1],input.data$Ne21[which(input.data$Ne21!=0)]); lim.x.Ne[2]<- max(lim.x.Ne[2],input.data$Ne21[which(!is.na(input.data$Ne21))])
           } else {Ne21.flag='n'}
if(sum(input.data$Al26)!=0 && is.numeric(sum(input.data$Al26))){Al26.flag='y'
           input.data$Al26[which(input.data$Al26==0)] <- NA; input.data$Al26.err[which(input.data$Al26==0)] <- NA
           lim.x.Al[1]<-min(lim.x.Al[1],input.data$Al26[which(input.data$Al26!=0)]); lim.x.Al[2]<-max(lim.x.Al[2],input.data$Al26)
           } else {Al26.flag='n'}
