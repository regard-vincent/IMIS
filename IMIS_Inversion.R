# # IMIS
# # INVERSION OF MULTI ISOTOPES IN SEDIMENTARY BASIN
# # --------------------
# # Vincent REGARD 
# # --------------------
# # 
# # Inversion module
# #
 
##_______________________________________________________________
##|                                                             |  
##| Initial data preparation                                    |
##|_____________________________________________________________|

n.step <- as.integer(max(c(t.b,na.omit(max.t.b)))/t.step)
if(n.step != (max(c(t.b,na.omit(max.t.b)))/t.step)){print(c("CAUTION: number of time steps is decimal"))}

	

##_______________________________________________________________
##|                                                             |  
##| CREATING A PRODUCTION vs DEPTH canevas                      |
##| calculation of production for various depths                |
##| in the program production is then interpolated              | 
##| (accelerates the process)                                   |
##|_____________________________________________________________|
if(progress.bar == 'txt'){print("Production vs. depth calculation");pb <- txtProgressBar(min = 0,  max = 4, char = "=",style=3)
}else if(progress.bar == 'win'){pb <- tkProgressBar(title = "Production vs. depth calculation", min = 0,  max = 4, width = 500)
}
z.prod <- seq(from=0,to=max(z.b),by=z.step.P )
if(Be10.flag=="y"){
	if(progress.bar == 'txt'){setTxtProgressBar(pb, 1)}else if(progress.bar == 'win'){setTkProgressBar(pb, 1, label="Be-10")}
	if(Prod.par=="lu.h"){PROD.Be.vs.z <- Production(z.prod*rho.b/10,abs(Latitude),Elevation,"Be",substr(Prod.par,4,4))
	                 # conversion z.lu (m) in g/cm2 : 1m <-> 100 x rho.b/1000 
		depth.prodBe <- rowSums(PROD.Be.vs.z)
	}else{
		suppressWarnings(PROD.Be.vs.z <- Production.z(as.numeric(Po.Be.b)[1],rho.b,z.prod,Att[1])) #neutrons
		for(i in 2:length(Po.Be.b)){
			suppressWarnings(PROD.Be.vs.z <- PROD.Be.vs.z + Production.z(as.numeric(Po.Be.b)[i],rho.b,z.prod,Att[i])) #muons
		}
		depth.prodBe <- PROD.Be.vs.z
	}
}
if(Al26.flag=="y"){
	if(progress.bar == 'txt'){setTxtProgressBar(pb, 2)}else if(progress.bar == 'win'){setTkProgressBar(pb, 2, label="Al-26")}
	
	if(Prod.par=="lu.h"){
		PROD.Al.vs.z <- Production(z.prod*rho.b/10,abs(Latitude),Elevation,"Al",substr(Prod.par,4,4))
		depth.prodAl <- rowSums(PROD.Al.vs.z)
	}else{
		suppressWarnings(PROD.Al.vs.z <- Production.z(as.numeric(Po.Al.b)[1],rho.b,z.prod,Att[1])) #neutrons
		for(i in 2:length(Po.Al.b)){
			suppressWarnings(PROD.Al.vs.z <- PROD.Al.vs.z + Production.z(as.numeric(Po.Al.b)[i],rho.b,z.prod,Att[i])) #muons
		}
		depth.prodAl <- PROD.Al.vs.z
	}
}
if(Ne21.flag=="y"){
	if(progress.bar == 'txt'){setTxtProgressBar(pb, 3)}else if(progress.bar == 'win'){setTkProgressBar(pb, 3, label="Ne-21")}
	if(Prod.par!="lu.e"){
		suppressWarnings(PROD.Ne.vs.z <- Production.z(as.numeric(Po.Ne.b)[1],rho.b,z.prod,Att[1])) #neutrons
		for(i in 2:length(Po.Ne.b)){
			suppressWarnings(PROD.Ne.vs.z <- PROD.Ne.vs.z + Production.z(as.numeric(Po.Ne.b)[i],rho.b,z.prod,Att[i])) #muons
		}
		depth.prodNe <- PROD.Ne.vs.z
	}else{
		PROD.Ne.vs.z <- Production(z.prod*rho.b/10,abs(Latitude),Elevation,"Ne",'h')
	                 # conversion z.lu (m) in g/cm2 : 1m <-> 100 x rho.b/1000 
		depth.prodNe <- rowSums(PROD.Ne.vs.z)
	}
}
if(progress.bar %in% c('txt','win')){close(pb)}


##_______________________________________________________________
##|                                                             |  
##| MATRIX setup. Size  n.step.var x n.iterations               |
##|    lines (1st index): Ages                                  |
##|    columns (2nd index): number of calcul iteration          |
##|_____________________________________________________________|


## generation scenarii output:                                [size]
## scenarii$breaks: bounds of different stages                [nbounds   x n.iterations]
## scenarii$t.step: duration of a time step during each stage [nbounds   x n.iterations]
## scenarii$n.step: number of steps during each stage         [nbounds-1]
scenarii <- Generate.tb()
n.step.cum <- c(0,cumsum(scenarii$n.step))   # cumulative sum, beginning at 0

## Creation erosion matrix a
MATRIX.epsilon.a <- Generate.eros.matrix()

## Creation matrix sr.b, of sedimentation rates during part b.
## verification it is possible
if(Repeat.step > min(scenarii$n.step)){
	Repeat.step <- min(scenarii$n.step)
	sprintf("nb of steps exceed minimum stage number of steps\n Repeat.step set to %i\n",Repeat.step) 
}

##########################
## Filling scenario (b)  #
##########################
z.b.mat <- time.matrix <- MATRIX.sr.b <- matrix(0,n.step,n.iterations)
if(Repeat.step==1){
	for(i.stage in 1:(length(t.b)-1)){
		MATRIX.sr.b[(n.step.cum[i.stage]+1):n.step.cum[i.stage+1],] <- t(matrix((z.b[i.stage]-z.b[i.stage+1]) /  # thickness
		          (scenarii$breaks[i.stage,]-scenarii$breaks[i.stage+1,]),n.iterations,scenarii$n.step[i.stage]))   # duration
		time.matrix[(n.step.cum[i.stage]+1):n.step.cum[i.stage+1],] <- t(matrix(scenarii$breaks[i.stage,],n.iterations,scenarii$n.step[i.stage]))+
		          (t(matrix(-scenarii$t.step[i.stage,],n.iterations,scenarii$n.step[i.stage])) * 
		          (matrix(1:(scenarii$n.step[i.stage]),scenarii$n.step[i.stage],n.iterations)))	           
	}
} else{
	n.step.breaks.cum <- matrix(0,Repeat.step*(length(n.step.cum)-1)+1, n.iterations) ; z.b.breaks <- n.step.breaks.cum # cf preceding
	n.step.breaks.cum[Repeat.step*(0:(length(n.step.cum)-1))+1,]   <- n.step.cum[1:length(n.step.cum)]    # indexes multiples of Repeat.step are those from n.stpep.cum
	                                                                                     # in between are inner breaks (calculated in the following)
	z.b.breaks[Repeat.step*(0:(length(n.step.cum)-1))+1,]   <- z.b[1:length(n.step.cum)]                  # same processing for z                                                             
	for(i.stage in 1:(length(t.b)-1)){
		n.step.breaks.cum[(Repeat.step*(i.stage-1))+(2:Repeat.step),] <- 
		               mapply(matrix(scenarii$n.step[i.stage],1,n.iterations),FUN=Sample.ordered)+n.step.cum[i.stage] 
		z.b.breaks[(Repeat.step*(i.stage-1))+(2:Repeat.step),] <- 
		               (mapply(matrix(999,1,n.iterations),FUN=Sample.ordered)/1000*(z.b[i.stage+1]-z.b[i.stage]))+z.b[i.stage] #idem for z
	}
	for(i.loc in 1:(dim(z.b.breaks)[1]-1)){
		for (j.it in 1:n.iterations){
			i.stage <- floor((i.loc+Repeat.step-1)/Repeat.step)
			MATRIX.sr.b[(n.step.breaks.cum[i.loc,j.it]+1):n.step.breaks.cum[i.loc+1,j.it],j.it] <- 
			          (z.b.breaks[i.loc+1,j.it]-z.b.breaks[i.loc,j.it]) /  # thickness
			          ((n.step.breaks.cum[i.loc,j.it]-n.step.breaks.cum[i.loc+1,j.it])*scenarii$t.step[i.stage,j.it])   # duration
			time.matrix[(n.step.breaks.cum[i.loc,j.it]+1):n.step.breaks.cum[i.loc+1,j.it],j.it] <- max(0,time.matrix[n.step.breaks.cum[i.loc,j.it],j.it]) +
			          (scenarii$t.step[i.stage,j.it]) * 
			          (1:(n.step.breaks.cum[i.loc+1,j.it]-n.step.breaks.cum[i.loc,j.it]))	           
	    }           
	}
	time.matrix <- max(t.b)-time.matrix
}
d.z.b <- (time.matrix[1:(n.step-1),]-time.matrix[2:n.step,]) * MATRIX.sr.b[1:(n.step-1),]
z.b.mat[2:n.step,] <- max(z.b)- apply(d.z.b,FUN=cumsum, MARGIN=2)
z.b.mat[1,] <- max(z.b) ; z.b.mat[which(z.b.mat<0)] <- 0 




##_______________________________________________________________
##|                                                             |  
##| Concentration calculation (direct)                          |
##|_____________________________________________________________|


# Hillslope erosion
#-------------------
if(Be10.flag=="y"){	CONCENTRATION.a.Be <- as.matrix(Concentration.Erosion(MATRIX.epsilon.a, sum(Po.Be.a), rho.a, L.Be))
	dim(CONCENTRATION.a.Be)<-c(n.step,n.iterations); CONCENTRATION.fin.Be <- CONCENTRATION.a.Be}
if(Ne21.flag=="y"){	CONCENTRATION.a.Ne <- as.matrix(Concentration.Erosion(MATRIX.epsilon.a, sum(Po.Ne.a), rho.a, 0))
	dim(CONCENTRATION.a.Ne)<-c(n.step,n.iterations); CONCENTRATION.fin.Ne <- CONCENTRATION.a.Ne}
if(Al26.flag=="y"){	CONCENTRATION.a.Al <- as.matrix(Concentration.Erosion(MATRIX.epsilon.a, sum(Po.Al.a), rho.a, L.Al))
	dim(CONCENTRATION.a.Al)<-c(n.step,n.iterations); CONCENTRATION.fin.Al <- CONCENTRATION.a.Al}

	
# Burial
#-------------------
if(Be10.flag=="y"){dC.b.Be <- CONCENTRATION.b.Be <- matrix(data = 0, nrow=n.step,ncol=n.iterations)}
if(Ne21.flag=="y"){dC.b.Ne <- CONCENTRATION.b.Ne <- matrix(data = 0, nrow=n.step,ncol=n.iterations)}
if(Al26.flag=="y"){dC.b.Al <- CONCENTRATION.b.Al <- matrix(data = 0, nrow=n.step,ncol=n.iterations)}

## i.step=1, initiation: hillslope concentration decreased by radioactive decay ----------------------------
CONCENTRATION.a.Be[1,] <- CONCENTRATION.a.Be[1,]*exp(-L.Be * n.step*scenarii$t.step[1,])
if(Al26.flag=="y"){CONCENTRATION.a.Al[1,] <- CONCENTRATION.a.Al[1,]*exp(-L.Al * n.step*scenarii$t.step[1,])}

# Loop -----------------------------------------------------------------------------------------------------
# progress bar
if(progress.bar == 'txt'){print("Concentration calculation");pb <- txtProgressBar(min = 0,  max = n.step, char = "=",style=3)
}else if(progress.bar == 'win'){pb <- tkProgressBar(title = "Concentration calculation", min = 0,  max = n.step, width = 500)
}

for(i.step in 2:n.step){
	index.stage <- max(which(n.step.cum < i.step)) # note this function needs the definition of i.step
	                                                          # the arg is istep+1, because we look at the former time step.
	# in the following, the former sample depth is given by : 
	depth <- z.b.mat[1:i.step,]- t(matrix(z.b.mat[i.step,],n.iterations,i.step))
	index.prod <-  round(depth/z.step.P)+1
	if(Ne21.flag=="y"){	
		dC.b.Ne[1:i.step,]<- depth.prodNe[index.prod] * t(matrix(scenarii$t.step[index.stage,],nrow=n.iterations,ncol=i.step))
		CONCENTRATION.b.Ne[1:i.step,] <- CONCENTRATION.b.Ne[1:i.step,] + dC.b.Ne[1:i.step,]
	}
	if(Be10.flag=="y"){	
		dC.b.Be[1:i.step,]<- depth.prodBe[index.prod] * t(matrix(scenarii$t.step[index.stage,],nrow=n.iterations,ncol=i.step))
		CONCENTRATION.b.Be[1:i.step,] <- CONCENTRATION.b.Be[1:i.step,]*exp(-L.Be*(scenarii$t.step[index.stage,])) + dC.b.Be[1:i.step,]
		CONCENTRATION.a.Be[i.step,] <- CONCENTRATION.a.Be[1,]*exp(-L.Be * (n.step+1-i.step)*scenarii$t.step[index.stage,])
	}
	if(Al26.flag=="y"){	
		dC.b.Al[1:i.step,]<- depth.prodAl[index.prod] * t(matrix(scenarii$t.step[index.stage,],nrow=n.iterations,ncol=i.step))
		CONCENTRATION.b.Al[1:i.step,] <- CONCENTRATION.b.Al[1:i.step,]*exp(-L.Al*(scenarii$t.step[index.stage,])) + dC.b.Al[1:i.step,]
		CONCENTRATION.a.Al[i.step,] <- CONCENTRATION.a.Al[1,]*exp(-L.Al * (n.step+1-i.step)*scenarii$t.step[index.stage,])
	}
	if(progress.bar == 'txt'){setTxtProgressBar(pb,  i.step)
	}else if(progress.bar == 'win'){setTkProgressBar(pb, i.step, label=paste( "Processing time step number", i.step, "of", n.step))
	}
}

if(progress.bar %in% c('txt','win')){close(pb)}

#remove(dC.b.Be, dC.b.Ne, dC.b.Al)
if(Be10.flag=="y"){CONCENTRATION.fin.Be <- CONCENTRATION.a.Be + CONCENTRATION.b.Be}
if(Ne21.flag=="y"){CONCENTRATION.fin.Ne <- CONCENTRATION.a.Ne + CONCENTRATION.b.Ne}
if(Al26.flag=="y"){CONCENTRATION.fin.Al <- CONCENTRATION.a.Al + CONCENTRATION.b.Al}



##_______________________________________________________________
##|                                                             |  
##| Chi2 calculation                                            |
##|_____________________________________________________________|

Chi2.Tout <- matrix(data = 0, nrow=1,ncol=n.iterations)
if(Be10.flag=="y"){
	Chi2.Be <- matrix(data = NA, nrow=1,ncol=n.iterations)
	for(i.iteration in 1:n.iterations){
		Chi2.Be[i.iteration] <- Chi2.data(input.data$Depth, input.data$Be10, input.data$Be10.err, 
		         CONCENTRATION.fin.Be[,i.iteration], z.b.mat[,i.iteration])
	}
	Chi2.Tout <- Chi2.Tout + Chi2.Be
}
if(Ne21.flag=="y"){
	Chi2.Ne <- matrix(data = NA, nrow=1,ncol=n.iterations)
	for(i.iteration in 1:n.iterations){
		Chi2.Ne[i.iteration] <- Chi2.data(input.data$Depth, input.data$Ne21, input.data$Ne21.err, CONCENTRATION.fin.Ne[,i.iteration], z.b.mat[,i.iteration])
	}
	Chi2.Tout <- Chi2.Tout + Chi2.Ne
}
if(Al26.flag=="y"){
	Chi2.Al <- matrix(data = NA, nrow=1,ncol=n.iterations)
	for(i.iteration in 1:n.iterations){
		Chi2.Al[i.iteration] <- Chi2.data(input.data$Depth, input.data$Al26, input.data$Al26.err, 
		          CONCENTRATION.fin.Al[,i.iteration], z.b.mat[,i.iteration])
	}
	Chi2.Tout <- Chi2.Tout + Chi2.Al
}

Best.Be <- which(Chi2.Be==min(Chi2.Be,na.rm=T)); Chi2.tout.inv <- Chi2.Be.inv <- Chi2.Be[Best.Be]/Chi2.Be
if(Al26.flag=="y"){
	if(progress.bar != 'no') {print(c("Chi2 10Be min", min(Chi2.Be,na.rm=T), "Chi2 26Al min",min(Chi2.Al,na.rm=T)))}
	Best.Al <- which(Chi2.Al == min(Chi2.Al,na.rm=T))
}
if(Ne21.flag=="y"){
	if(progress.bar != 'no') {print(c("Chi2 10Be min", min(Chi2.Be,na.rm=T), "Chi2 21Ne min",min(Chi2.Ne,na.rm=T)))}
	Best.Ne <- which(Chi2.Ne==min(Chi2.Ne,na.rm=T))
}
Best.tout <- which(Chi2.Tout==min(Chi2.Tout,na.rm=T));Worse.tout <- which(Chi2.Tout==max(Chi2.Tout,na.rm=T));
Chi2.tout.inv <- 1/Chi2.Tout; Chi2.tout.inv[is.na(Chi2.tout.inv)]<-0
Order.chi <- order(t(Chi2.tout.inv))



