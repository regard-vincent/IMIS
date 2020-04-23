##_______________________________________________________________
##|                                                             |  
##|   IMIS INVERSION PROGRAM: GRAPHIC OUTPUTS                   |
##|_____________________________________________________________|
##
## Graphical tuning parameters are described into IMIS_graphical_options.R
##
## Graphs drawn:
##    1- Depth vs age (name: 'DA')
##    2- Denudation rate vs time (name: 'DT')
##    3- Nuclide Concentration vs depth; Nuclides are 10Be (Blue), 21Ne (Red), 26Al(green)  (names: 'NB', 'NN', 'NA')
##    4- Sedimentation rate vs Age (name: 'SA')
##    5- Sedimentation rate vs depth (name: 'SD')
##    ALL (name: 'ALL')


##########################################
###                                    ###
###    Preparation                     ###
###                                    ###
##########################################

Space.top.legend <- 2      ### defines the size of the space for legend if nc.vs.z.legend.position == 't'
cm <- Chi2.tout.inv[Best.tout]
if(cm<0.9){Chi2.tout.inv <- Chi2.tout.inv*0.9/cm}
z.b.lim.im <- min(z.lim.im , max(z.b))
### Finding Delta_Chi depending on the degrees de freedom (obs) : qchisq(0.954, df=##) 95% ou qchisq(0.683, df=###)
Chi.Be95 <- min(Chi2.Be,na.rm=T) + qchisq(0.954, df=-1+length(na.omit(input.data$Be10)))
Chi.Be68 <- min(Chi2.Be,na.rm=T) + qchisq(0.683, df=-1+length(na.omit(input.data$Be10)))
Acceptable.Be95 <- which(Chi2.Be <= Chi.Be95)
Acceptable.Be68 <- which(Chi2.Be <= Chi.Be68)
N.Scenarii.Be.68 <- length(Acceptable.Be68); 
N.Scenarii.Be.95 <- length(Acceptable.Be95); 

if(Ne21.flag=="y"){
	Chi.Ne95   <- min(Chi2.Ne,na.rm=T)  + qchisq(0.954, df=-1+length(na.omit(input.data$Ne21)))
	Chi.Ne68   <- min(Chi2.Ne,na.rm=T)  + qchisq(0.683, df=-1+length(na.omit(input.data$Ne21)))
	Chi.Tout95 <- min(Chi2.Tout,na.rm=T)+ qchisq(0.954, df=-1+length(na.omit(input.data$Ne21))+length(na.omit(input.data$Be10)))
	Chi.Tout68 <- min(Chi2.Tout,na.rm=T)+ qchisq(0.683, df=-1+length(na.omit(input.data$Ne21))+length(na.omit(input.data$Be10)))
	Acceptable.Ne95   <- which(Chi2.Ne   <= Chi.Ne95);	Acceptable.Ne68 <- which(Chi2.Ne <= Chi.Ne68)
	Acceptable.Tout95 <- which(Chi2.Tout <= Chi.Tout95);Acceptable.Tout68 <- which(Chi2.Tout <= Chi.Tout68)
	N.Scenarii.Ne.68 <- length(Acceptable.Ne68); N.Scenarii.Ne.95 <- length(Acceptable.Ne95)
}   
if(Al26.flag=="y"){
	Chi.Al95   <- min(Chi2.Al,na.rm=T)  + qchisq(0.954, df=-1+length(na.omit(input.data$Al26)))
	Chi.Al68   <- min(Chi2.Al,na.rm=T)  + qchisq(0.683, df=-1+length(na.omit(input.data$Al26)))
	Chi.Tout95 <- min(Chi2.Tout,na.rm=T)+ qchisq(0.954, df=-1+length(na.omit(input.data$Al26))+length(na.omit(input.data$Be10)))
	Chi.Tout68 <- min(Chi2.Tout,na.rm=T)+ qchisq(0.683, df=-1+length(na.omit(input.data$Al26))+length(na.omit(input.data$Be10)))
	Acceptable.Al95   <- which(Chi2.Al   <= Chi.Al95);	Acceptable.Al68 <- which(Chi2.Al <= Chi.Al68)
	Acceptable.Tout95 <- which(Chi2.Tout <= Chi.Tout95);Acceptable.Tout68 <- which(Chi2.Tout <= Chi.Tout68)
	N.Scenarii.Al.68 <- length(Acceptable.Al68); N.Scenarii.Al.95 <- length(Acceptable.Al95)
}
if(Al26.flag=="y" | Ne21.flag=="y"){
	Chi.Tout95 <- min(Chi2.Tout,na.rm=T)+ qchisq(0.954, df=-1+length(na.omit(input.data$Ne21))+length(na.omit(input.data$Al26))+length(na.omit(input.data$Be10)))
	Chi.Tout68 <- min(Chi2.Tout,na.rm=T)+ qchisq(0.683, df=-1+length(na.omit(input.data$Ne21))+length(na.omit(input.data$Al26))+length(na.omit(input.data$Be10)))
	Acceptable.Tout95 <- which(Chi2.Tout <= Chi.Tout95);Acceptable.Tout68 <- which(Chi2.Tout <= Chi.Tout68)
	N.Scenarii.Tout.68 <- length(Acceptable.Tout68); N.Scenarii.Tout.95 <- length(Acceptable.Tout95)
}


if(length(which(test.lib=="raster"))==0)      {install.packages("raster")}       ; 	library(raster)
if(length(which(test.lib=="RColorBrewer"))==0){install.packages("RColorBrewer")} ; 	library(RColorBrewer)

suppressWarnings(if(Graphs.drawn %in% 'ALL'){Graphs.drawn=c('DA','DT','NB','NN','NA','SA','ST')})

### ===========================
###    PRINT RESULTS IN FILE
### ===========================
if(Ne21.flag=="y"){
	if(Al26.flag=="y"){
		write.table(t(c(min(Chi2.Be,na.rm=T),min(Chi2.Ne,na.rm=T),min(Chi2.Al,na.rm=T),min(Chi2.Tout,na.rm=T),N.Scenarii.Be.68,N.Scenarii.Be.95,N.Scenarii.Ne.68,N.Scenarii.Ne.95,N.Scenarii.Al.68,N.Scenarii.Al.95,N.Scenarii.Tout.68,N.Scenarii.Tout.95,info)), 
	                                file=file.out, col.names=FALSE, row.names=FALSE, append=TRUE, sep="\t")
    } else {
		write.table(t(c(min(Chi2.Be,na.rm=T),min(Chi2.Ne,na.rm=T),min(Chi2.Tout,na.rm=T),N.Scenarii.Be.68,N.Scenarii.Be.95,N.Scenarii.Ne.68,N.Scenarii.Ne.95,N.Scenarii.Tout.68,N.Scenarii.Tout.95,info)), 
	                                file=file.out, col.names=FALSE, row.names=FALSE, append=TRUE, sep="\t")
	}
} else {	             
	write.table(t(c(min(Chi2.Be,na.rm=T),min(Chi2.Al,na.rm=T),min(Chi2.Tout,na.rm=T),N.Scenarii.Be.68,N.Scenarii.Be.95,N.Scenarii.Al.68,N.Scenarii.Al.95,N.Scenarii.Tout.68,N.Scenarii.Tout.95,info)), 
	                                file=file.out, col.names=FALSE, row.names=FALSE, append=TRUE, sep="\t")
}

  
#########################################
##                                    ###
##    Depth vs age                    ###
##                                    ###
#########################################

if('DA' %in% Graphs.drawn){
	minmax <- matrix(0,2, n.step); minmax[1,] <- apply(z.b.mat,1,min); minmax[2,] <- apply(z.b.mat,1,max)
	timmax <- matrix(0,2, n.step); timmax[1,] <- apply(time.matrix,1,max); timmax[2,] <- apply(time.matrix,1,min)
	z.b.lim.im <- min(z.lim.im , max(z.b.mat[1,]))
	
	
	### in case option 'map', rasterisation is only computed once
	### ---------------------------------------------------------
	if(d.vs.t.type == 'm'){
		### raster definition
		e <- extent(c(0,n.step*t.step/1e6,-max(z.b),0))
		r <- raster(e, nrows=ceiling(max(z.b)/z.step.im), ncols=ceiling(n.step*t.step/t.step.im))
		
		### Interpolation from  points along each path at regularly spaced depths
		n.z.regular <- ceiling(max(z.b)/z.step.im)
		time.regular <- matrix(NA, nrow=n.z.regular, ncol=n.iterations)
		z.b.regular <- seq(res(r)[2]/2, -res(r)[2]/2 + res(r)[2]*n.z.regular,res(r)[2])
		for(i.z in 1:n.iterations){
			time.regular[,i.z] <- approx(z.b.mat[,i.z], time.matrix[,i.z],z.b.regular)$y
			}
		dim(time.regular) <- n.z.regular*n.iterations
		z.b.regularM <- matrix(z.b.regular,ncol=n.iterations,nrow=n.z.regular,byrow =F)
		dim(z.b.regularM) <- n.z.regular*n.iterations
		if (chi.scenario != 'i'){  ### default is 'n'
			### preparation raster map of Chi2 values
			time.z.chi.Be <- matrix(c(time.regular/1e6,-z.b.regularM,matrix(log10(Chi2.Be),ncol=n.iterations,nrow=n.z.regular,byrow =T)), ncol=3)
			if(Ne21.flag=="y"){time.z.chi.Ne <- matrix(c(time.regular/1e6,-z.b.regularM,matrix(log10(Chi2.Ne),ncol=n.iterations,nrow=n.z.regular,byrow =T)), ncol=3)}
			if(Al26.flag=="y"){time.z.chi.Al <- matrix(c(time.regular/1e6,-z.b.regularM,matrix(log10(Chi2.Al),ncol=n.iterations,nrow=n.z.regular,byrow =T)), ncol=3)}
			if(Al26.flag=="y" | Ne21.flag=="y"){
				time.z.chi <- matrix(c(time.regular/1e6,-z.b.regularM,matrix(log10(Chi2.Tout),ncol=n.iterations,nrow=n.z.regular,byrow =T)), ncol=3)
				rast.chi <- rasterize(na.omit(time.z.chi)[,1:2],r,(na.omit(time.z.chi))[,3], fun=function(x,...)min(x)) ### min
			} else {
				rast.chi <- rasterize(na.omit(time.z.chi.Be)[,1:2],r,(na.omit(time.z.chi)[,3]), fun=function(x,...)min(x)) ### min
			}
			rast.min <- summary(rast.chi)[1]; rast.q75 <- summary(rast.chi)[4]
			### raster plot preparation
			raster.coulour <- colorRampPalette(brewer.pal(name="YlOrBr",n=4))  ### https://moderndata.plot.ly/create-colorful-graphs-in-r-with-rcolorbrewer-and-plotly/
		}else{
			time.z.chi.inv <- matrix(c(time.regular/1e6,-z.b.regularM,matrix(Chi2.tout.inv,ncol=n.iterations,nrow=n.z.regular,byrow =T)), ncol=3)
			rast.chi.inv <- rasterize(na.omit(time.z.chi.inv)[,1:2],r,na.omit(time.z.chi.inv)[,3], fun=function(x,...)max(x)) ### max
			
			#### equalisation by z
			rast.chi.inv.eq <- r
			for(i.z in 1:dim(rast.chi.inv)[1]){ ### for eah z strip
				max.chi.z <- max(na.omit(rast.chi.inv[i.z,]))
				rast.chi.inv.eq[i.z,] <- rast.chi.inv[i.z,]/max.chi.z
			}
# 			raster.coulour <- colorRampPalette(brewer.pal(name="Greys",n=5))  ### https://moderndata.plot.ly/create-colorful-graphs-in-r-with-rcolorbrewer-and-plotly/
# 			raster.coulour <- colorRampPalette(brewer.pal(name="Oranges",n=5))  ### https://moderndata.plot.ly/create-colorful-graphs-in-r-with-rcolorbrewer-and-plotly/
			raster.coulour <- colorRampPalette(brewer.pal(name="YlOrBr",n=4))  ### https://moderndata.plot.ly/create-colorful-graphs-in-r-with-rcolorbrewer-and-plotly/
		}
	}
	
	
	### GRAPHICS
	for(i.graphe in 1:length(graphs)){
		### OPENING graph files
		if(graphs[i.graphe] =="pdf")      {pdf(file=paste(outpath,"/Graphe_scenarios.pdf",sep=""), width=w, height=h, pointsize=15, paper="a4r")
		}else if(graphs[i.graphe] =="png"){png(file=paste(outpath,"/Graphe_scenarios.png",sep=""), width=w, height=h, pointsize=15)
		}
	
		
		### PLOT
		if(d.vs.t.type == 'l'){ 
			par(mar=c(4,4,2,7))
			### draw scenario lines
			### the following lines draw a graph in which each scenario is drawn with a line whose color corresponds to scenario probability
			### ----------------------------------------------
			plot(time.matrix[,1]/1e6, -z.b.mat[,1],type='l',xlab='Time [Myrs]', ylab='Depth [m]',ylim=c(-z.b.lim.im ,0), col="grey", xaxp=c(0,100,50))
			polygon(c(timmax[1,1:n.step],timmax[2,n.step:1])/1e6, -c(minmax[1,1:n.step],minmax[2,n.step:1]), 
			           col="moccasin", border=NA)
			
			### draw the copper bearing layer 
			if(length(copper.layer.depth)>1){ polygon(max(timmax)*c(0,1,1,0)/1e6, -c(min(copper.layer.depth),min(copper.layer.depth),max(copper.layer.depth),max(copper.layer.depth)), col="paleturquoise2", border=NA)}
						for(i.iteration in 1:n.iterations){lines(time.matrix[,Order.chi[i.iteration]]/1e6, -z.b.mat[,Order.chi[i.iteration]], 
	         		col=grey(1-min(na.omit(Chi2.tout.inv[Order.chi[i.iteration]]), !is.na(Chi2.tout.inv[Order.chi[i.iteration]])) ))
	         	}
	         	### draw the scenarios
			spaced.col <- seq(1-Chi2.tout.inv[Worse.tout[1]], 1-Chi2.tout.inv[Best.tout[1]],by=(Chi2.tout.inv[Worse.tout[1]]-Chi2.tout.inv[Best.tout[1]])/10)
			color.legend(max(t.b)*1.05/1e6,-z.b.lim.im ,max(t.b)*1.07/1e6,0,c("low probability","","","","","high probability"),align="rb",grey(spaced.col),gradient="y")
			### best chi2
			if(Ne21.flag=="y"){lines(time.matrix[,Best.Ne]/1e6, -z.b.mat[,Best.Ne],col='red', lwd=2)}
			if(Al26.flag=="y"){lines(time.matrix[,Best.Al]/1e6, -z.b.mat[,Best.Al],col='chartreuse3', lwd=2)}
			lines(time.matrix[,Best.Be]/1e6, -z.b.mat[,Best.Be],col="blue", lwd=2)
			lines(time.matrix[,Best.tout]/1e6, -z.b.mat[,Best.tout],col="black", lwd=2)
			for(i.iteration in Best.tout){lines(time.matrix[,i.iteration]/1e6, -z.b.mat[,i.iteration],col='black', lwd=2)}
			polygon(c(timmax[1,1:n.step],timmax[2,n.step:1])/1e6, -c(minmax[1,1:n.step],minmax[2,n.step:1]), border="moccasin", lwd=4)
			### fixed points
			points(t.b/1e6,-z.b, pch=20, cex=2, col="gray60")
			for (i in (which(t.b.err!=0))){
			   Vect_err.y <- Vect_err.x <- 0; Vect_err.y[2] <- Vect_err.y[1] <- z.b[i]; Vect_err.x[1] <- t.b[i]+t.b.err[i]; Vect_err.x[2] <- t.b[i]-t.b.err[i]
			  lines(Vect_err.x/1e6,-Vect_err.y, lwd=10, col="gray60")
			}
			if(Ne21.flag=="y"){
				if(Al26.flag=="y"){
					legend(horiz=F, "topright", c("tested scenario", "best fit scenario", "best Be", "best Ne", "best Al", "constrains"), 
			              col=c("grey", "black", "dark blue","red", "chartreuse3", "gray40"), bg="white", lwd=c(1,2,2,2,2,10))  ##, bty='n'
			    } else {
					legend(horiz=F, "topright", c("tested scenario", "best fit scenario", "best Be", "best Ne", "constrains"), 
			              col=c("grey", "black", "dark blue","red", "gray40"), bg="white", lwd=c(1,2,2,2,10))  ##, bty='n'
				}
			} else {	             
					legend(horiz=F, "topright", c("tested scenario", "best fit scenario", "best Be", "best Al", "constrains"), 
			              col=c("grey", "black", "dark blue","chartreuse3", "gray40"), bg="white", lwd=c(1,2,2,2,10))  ##, bty='n'
			}
		
		}else{   #### if(d.vs.t.type == 'l'  --> d.vs.t.type='m'
		
			### draw scenario probabilities
			### the following lines draw a graph in which chi2 is mapped
			### ----------------------------------------------
			par(mar=c(4,4,2,5))
			if (chi.scenario != 'i'){
				### raster plot
				image(rast.chi, ylim=c(-z.b.lim.im ,0),xlab='Time [Myrs]', ylab='Depth [m]', xaxp=c(0,100,50))
				if(m.type==1){plot(rast.chi, col='white', colNA='light grey', add=T, legend=FALSE)}else{plot(rast.chi, col='grey90', colNA='white',add=T, legend=FALSE)}
				### draw the copper bearing layer 
				if(length(copper.layer.depth)>1){ polygon(max(timmax)*c(0,1,1,0)/1e6, -c(min(copper.layer.depth),min(copper.layer.depth),
					max(copper.layer.depth),max(copper.layer.depth)), col="paleturquoise2", border=NA)}

				if(m.type!=1){
					plot(rast.chi, zlim=c(0,Chi.Tout95),col='grey80', add=T, legend=FALSE)			
					plot(rast.chi, zlim=c(0,Chi.Tout68),col='grey70', add=T, legend=FALSE)			
				}					
				plot(rast.chi, col=rev(raster.coulour(255)), zlim=c(rast.min,rast.q75), add=T, alpha=0.75, 
					legend.args=list(text='log(Chi2)', side=4, line=2.5, cex=0.8))

	
				### Draw simple probabilty limits
				if(d.vs.t.levels=='y'){
					time.z.chi68 <- time.z.chi.Be[which(time.z.chi.Be[,3]<=log10(Chi.Be68)),]
				 	time.z.chi95 <- time.z.chi.Be[which(time.z.chi.Be[,3]<=log10(Chi.Be95)),]
				 	max.chi68.z.Be <- min.chi68.z.Be <- max.chi95.z.Be <- min.chi95.z.Be <- NA
					for(i.z in 1:n.z.regular){ ### for eah z strip
						max.chi68.z.Be[i.z] <- max( na.omit(time.z.chi68[ which(time.z.chi68[,2] == -z.b.regular[i.z]),1]))
						min.chi68.z.Be[i.z] <- min( na.omit(time.z.chi68[ which(time.z.chi68[,2] == -z.b.regular[i.z]),1]))
						max.chi95.z.Be[i.z] <- max( na.omit(time.z.chi95[ which(time.z.chi95[,2] == -z.b.regular[i.z]),1]))
						min.chi95.z.Be[i.z] <- min( na.omit(time.z.chi95[ which(time.z.chi95[,2] == -z.b.regular[i.z]),1]))
						}
			
					### draw the curves for Be-10
					labcurve( list( One=list( max.chi68.z.Be,-z.b.regular), Two=list( min.chi68.z.Be,-z.b.regular), Three=list( max.chi95.z.Be,-z.b.regular), Four= list( min.chi95.z.Be,-z.b.regular)),
							keys=c('68%conf','','95%conf',''),  keyloc="none", lty=c("dashed","dashed","twodash","twodash"), lwd=1,
							col="blue", ###xlim=c(10,15), ### ylim=c(-100-20),
							pl=TRUE, add=T) 
					### idem for other nuclides
					if(Ne21.flag=="y"){
						time.z.chi68 <- time.z.chi.Ne[which(time.z.chi.Ne[,3]<=log10(Chi.Ne68)),]
					 	time.z.chi95 <- time.z.chi.Ne[which(time.z.chi.Ne[,3]<=log10(Chi.Ne95)),]
					 	max.chi68.z.Ne <- min.chi68.z.Ne <- max.chi95.z.Ne <- min.chi95.z.Ne <- NA
						for(i.z in 1:n.z.regular){ ### for eah z strip
							max.chi68.z.Ne[i.z] <- max( na.omit(time.z.chi68[ which(time.z.chi68[,2] == -z.b.regular[i.z]),1]))
							min.chi68.z.Ne[i.z] <- min( na.omit(time.z.chi68[ which(time.z.chi68[,2] == -z.b.regular[i.z]),1]))
							max.chi95.z.Ne[i.z] <- max( na.omit(time.z.chi95[ which(time.z.chi95[,2] == -z.b.regular[i.z]),1]))
							min.chi95.z.Ne[i.z] <- min( na.omit(time.z.chi95[ which(time.z.chi95[,2] == -z.b.regular[i.z]),1]))
							}
						labcurve( list( One=list(max.chi68.z.Ne,-z.b.regular), Two=list( min.chi68.z.Ne,-z.b.regular), Three=list( max.chi95.z.Ne,-z.b.regular), Four= list( min.chi95.z.Ne,-z.b.regular)),
						keys=c('68%conf','','95%conf',''),  keyloc="none", lty=c("dashed","dashed","twodash","twodash"), lwd=1,
						col="red", ###xlim=c(10,15), ### ylim=c(-100-20),
						pl=TRUE, add=T) 
						}
					if(Al26.flag=="y"){
						time.z.chi68 <- time.z.chi.Al[which(time.z.chi.Al[,3]<=log10(Chi.Al68)),]
					 	time.z.chi95 <- time.z.chi.Al[which(time.z.chi.Al[,3]<=log10(Chi.Al95)),]
					 	max.chi68.z.Al <- min.chi68.z.Al <- max.chi95.z.Al <- min.chi95.z.Al <- NA
						for(i.z in 1:n.z.regular){ ### for eah z strip
							max.chi68.z.Al[i.z] <- max( na.omit(time.z.chi68[ which(time.z.chi68[,2] == -z.b.regular[i.z]),1]))
							min.chi68.z.Al[i.z] <- min( na.omit(time.z.chi68[ which(time.z.chi68[,2] == -z.b.regular[i.z]),1]))
							max.chi95.z.Al[i.z] <- max( na.omit(time.z.chi95[ which(time.z.chi95[,2] == -z.b.regular[i.z]),1]))
							min.chi95.z.Al[i.z] <- min( na.omit(time.z.chi95[ which(time.z.chi95[,2] == -z.b.regular[i.z]),1]))
							}
						labcurve( list( One=list( max.chi68.z.Al,-z.b.regular), Two=list( min.chi68.z.Al,-z.b.regular), Three=list( max.chi95.z.Al,-z.b.regular), Four= list( min.chi95.z.Al,-z.b.regular)),
						keys=c('68%conf','','95%conf',''),  keyloc="none", lty=c("dashed","dashed","twodash","twodash"), lwd=1,
						col="chartreuse3", ###xlim=c(10,15), ### ylim=c(-100-20),
						pl=TRUE, add=T) 
						}
					### idem for multinuclides chi2.
					if(Al26.flag=="y" | Ne21.flag=="y"){
						time.z.chi68 <- time.z.chi[which(time.z.chi[,3]<=log10(Chi.Tout68)),]
					 	time.z.chi95 <- time.z.chi[which(time.z.chi[,3]<=log10(Chi.Tout95)),]
					 	max.chi68.z <- min.chi68.z <- max.chi95.z <- min.chi95.z <- NA
						for(i.z in 1:n.z.regular){ ### for eah z strip
							max.chi68.z[i.z] <- max( na.omit(time.z.chi68[ which(time.z.chi68[,2] == -z.b.regular[i.z]),1]))
							min.chi68.z[i.z] <- min( na.omit(time.z.chi68[ which(time.z.chi68[,2] == -z.b.regular[i.z]),1]))
							max.chi95.z[i.z] <- max( na.omit(time.z.chi95[ which(time.z.chi95[,2] == -z.b.regular[i.z]),1]))
							min.chi95.z[i.z] <- min( na.omit(time.z.chi95[ which(time.z.chi95[,2] == -z.b.regular[i.z]),1]))
							}
						labcurve( list( One=list( max.chi68.z,-z.b.regular), Two=list( min.chi68.z,-z.b.regular), Three=list( max.chi95.z,-z.b.regular), Four= list( min.chi95.z,-z.b.regular)),
						keys=c('68%conf','','95%conf',''),  keyloc="none", lty=c("dashed","dashed","twodash","twodash"), lwd=2,
						col="gray20", ###xlim=c(10,15), ### ylim=c(-100-20),
						pl=TRUE, add=T) 
						}
				}
		
			}else{   #### if (chi.scenario != 'i')
				image(rast.chi.inv.eq, ylim=c(-z.b.lim.im ,0),xlab='Time [Myrs]', ylab='Depth [m]', xaxp=c(0,100,50))
				if(m.type==1){plot(rast.chi.inv.eq, col='white', colNA='light grey', add=T, legend=FALSE)}else{plot(rast.chi.inv.eq, col='grey90', colNA='white', add=T, legend=FALSE)}
				### draw the copper bearing layer 
				if(length(copper.layer.depth)>1){ polygon(max(timmax)*c(0,1,1,0)/1e6, -c(min(copper.layer.depth),min(copper.layer.depth),max(copper.layer.depth),max(copper.layer.depth)), col="paleturquoise2", border=NA)}
				if(m.type!=1){
						plot(rast.chi.inv.eq, zlim=c(1/Chi.Tout95,1),col='grey80', add=T, legend=FALSE)			
						plot(rast.chi.inv.eq, zlim=c(1/Chi.Tout68,1),col='grey70', add=T, legend=FALSE)						
				}					
				plot(rast.chi.inv.eq, zlim=c(rast.min/rast.q75,1),col=raster.coulour(50), add=T, 
					legend.args=list(text='Chi2_min / Chi2', side=4, line=2.5, cex=0.8))
				if(d.vs.t.levels=='y'){contour(rast.chi.inv.eq, levels=0.8, add=T, drawlabels=F, col="white")}
			} #### end if (chi.scenario != 'i')
			
			### best chi2
			if(Ne21.flag=="y"){lines(time.matrix[,Best.Ne]/1e6, -z.b.mat[,Best.Ne],col='red', lwd=2)}
			if(Al26.flag=="y"){lines(time.matrix[,Best.Al]/1e6, -z.b.mat[,Best.Al],col='chartreuse3', lwd=2)}
			lines(time.matrix[,Best.Be]/1e6, -z.b.mat[,Best.Be],col="blue", lwd=2)
			lines(time.matrix[,Best.tout]/1e6, -z.b.mat[,Best.tout],col="black", lwd=2)
			for(i.iteration in Best.tout){lines(time.matrix[,i.iteration]/1e6, -z.b.mat[,i.iteration],col='black', lwd=2)}
			### fixed points
			points(t.b/1e6,-z.b, pch=20, cex=2, col="gray40")
			for (i in (which(t.b.err!=0))){
				Vect_err.y <- Vect_err.x <- 0; Vect_err.y[2] <- Vect_err.y[1] <- z.b[i]; Vect_err.x[1] <- t.b[i]+t.b.err[i]; Vect_err.x[2] <- t.b[i]-t.b.err[i]
				lines(Vect_err.x/1e6,-Vect_err.y, lwd=10, col="gray40")
				if(Ne21.flag=="y"){
					if(Al26.flag=="y"){
						legend(horiz=F, "topright", c("best fit scenario", "best Be", "best Ne", "best Al", "constrains"), 
				              col=c("black", "dark blue","red", "chartreuse3", "gray40"), bg="white", lwd=c(2,2,2,2,10))  ##, bty='n'
				    	} else {
						legend(horiz=F, "topright", c("best fit scenario", "best Be", "best Ne", "constrains"), 
				              col=c("black", "dark blue","red", "gray40"), bg="white", lwd=c(2,2,2,10))  ##, bty='n'
					}
				} else {	             
						legend(horiz=F, "topright", c("best fit scenario", "best Be", "best Al", "constrains"), 
				              col=c("black", "dark blue","chartreuse3", "gray40"), bg="white", lwd=c(2,2,2,10))  ##, bty='n'
				}
			}  #### end for (i in (which(t.b.err!=0)))
		} #### end if(d.vs.t.type == 'l')
		### draw sample depths
		rug(-input.data$Depth,side=2,lwd=2)
		
		### CLOSING GRAPH FILES
		if(graphs[i.graphe] =="png" |  graphs[i.graphe] =="pdf"){dev.off()}
	} ## END for(i.graphe in 1:length(graphs))
}



#########################################
##                                    ###
##    Denudation rate vs time         ###
##                                    ###
#########################################

if('DT' %in% Graphs.drawn){
	for(i.graphe in 1:length(graphs)){
		if(graphs[i.graphe] =="pdf"){pdf(file=paste(outpath,"/Graphe_erosion.pdf",sep=""), width=w, height=h, pointsize=15, paper="a4r")
		}else if(graphs[i.graphe] =="png"){png(file=paste(outpath,"/Graphe_erosion.png",sep=""), width=w, height=h, pointsize=15)
		}

 		if(nc.vs.z.legend.position != 'ot'){par(mar=c(4,4,2,7)) } else {Prop.bas <- 2; par(mar=c(4,4,2+Space.top.legend,7))}

 		plot(time.matrix[,1]/1e6, MATRIX.epsilon.a[,1]*1e6,type='l',xlab='Time [Myrs]', ylab='Denudation rate [m/Ma]', col="white")
		for(i.iteration in 1:n.iterations){lines(time.matrix[,na.omit(Order.chi[i.iteration])]/1e6, MATRIX.epsilon.a[,na.omit(Order.chi[i.iteration])]*1e6, 
		            col=grey(1-min(na.omit(Chi2.tout.inv[Order.chi[i.iteration]]), !is.na(Chi2.tout.inv[Order.chi[i.iteration]])) )
		            )}
		# best chi2
		if(Ne21.flag=="y"){lines(time.matrix[,Best.Ne]/1e6, MATRIX.epsilon.a[,Best.Ne]*1e6,col="red",lwd=2)}
		if(Al26.flag=="y"){lines(time.matrix[,Best.Al]/1e6, MATRIX.epsilon.a[,Best.Al]*1e6,col='chartreuse3',lwd=2)}
		lines(time.matrix[,Best.Be]/1e6, MATRIX.epsilon.a[,Best.Be]*1e6,col="dark blue",lwd=2)
		lines(time.matrix[,Best.tout]/1e6, MATRIX.epsilon.a[,Best.tout]*1e6,col="black", lwd=2)
		for(i.iteration in Best.tout){points(time.matrix[,i.iteration]/1e6, MATRIX.epsilon.a[,i.iteration]*1e6,col='red',pch=19,cex=0.4)}
		
		# Legend
		if(Ne21.flag=="y"){
			if(Al26.flag=="y"){ ### Ne21 and Al26
				if(nc.vs.z.legend.position != 'ot'){
					legend(horiz=F, "topleft", c("tested scenario", "best fit scenario", "best Be", "best Ne", "best Al"), 
		               			col=c("grey", "black", "dark blue","red", "chartreuse3"), 
		               			bg="white", lwd=c(1,NULL, 2,2,2), pch=c(NULL,19, NULL,NULL,NULL), pt.cex=c(NULL,0.4, NULL,NULL,NULL)) 
	               		} else {
		               		legend(ncol=2, legend=c("tested scenario", "best fit scenario", "best Be", "best Ne", "best Al"), 
		               			col=c("grey", "black", "dark blue","red", "chartreuse3"), 
		               			bg="white", lwd=c(1,NULL, 2,2,2), pch=c(NULL,19, NULL,NULL,NULL), pt.cex=c(NULL,0.4, NULL,NULL,NULL),
		               			y=par()$usr[4]*1.01, x=10^(par()$usr[1]), yjust = 0, xpd = TRUE) 
	               		}
			} else { ### Ne21 and no Al26
				if(nc.vs.z.legend.position != 'ot'){
					legend(horiz=F, "topleft", c("tested scenario", "best fit scenario", "best Be", "best Ne"), 
		               			col=c("grey", "black", "dark blue","red"), 
		               			bg="white", lwd=c(1,NULL, 2,2), pch=c(NULL,19, NULL,NULL), pt.cex=c(NULL,0.4, NULL,NULL)) 
	               		}else{
		               		legend(ncol=2, legend=c("tested scenario", "best fit scenario", "best Be", "best Ne"), 
		               			col=c("grey", "black", "dark blue","red"), 
		               			bg="white", lwd=c(1,NULL, 2,2), pch=c(NULL,19, NULL,NULL), pt.cex=c(NULL,0.4, NULL,NULL),
		               			y=par()$usr[4]*1.01, x=10^(par()$usr[1]), yjust = 0, xpd = TRUE) 
	               		}
			}
		} else { ### no Ne21 and Al26	             
			if(nc.vs.z.legend.position != 'ot'){			
				legend(horiz=F, "topleft", c("tested scenario", "best fit scenario", "best Be", "best Al"), 
		               		col=c("grey", "black", "dark blue", "chartreuse3"), 
		               		bg="white", lwd=c(1,NULL, 2,2), pch=c(NULL,19,NULL,NULL), pt.cex=c(NULL,0.4, NULL,NULL,NULL))
	               	} else {
		               legend(ncol=2, legend=c("tested scenario", "best fit scenario", "best Be", "best Al"), 
		               		col=c("grey", "black", "dark blue", "chartreuse3"), 
		               		bg="white", lwd=c(1,NULL, 2,2), pch=c(NULL,19,NULL,NULL), pt.cex=c(NULL,0.4, NULL,NULL),
		               		y=par()$usr[4]*1.01, x=10^(par()$usr[1]), yjust = 0, xpd = TRUE) 
	               	}
		}
		
		spaced.col <- seq(1-Chi2.tout.inv[Worse.tout[1]], 1-Chi2.tout.inv[Best.tout[1]],by=(Chi2.tout.inv[Worse.tout[1]]-Chi2.tout.inv[Best.tout[1]])/10)
		color.legend(max(t.b)*1.05/1e6,0,max(t.b)*1.07/1e6,max(epsilon.a)*1e6,c("low probability","","","","","high probability"),align="rb",grey(spaced.col),gradient="y")
		
		if(graphs[i.graphe] =="png" |  graphs[i.graphe] =="pdf"){dev.off()}
	}
}

#########################################
##                                    ###
##    10Be Concentration vs depth     ###
##                                    ###
#########################################
if('NB' %in% Graphs.drawn){
	for(i.graphe in 1:length(graphs)){
		if(graphs[i.graphe] =="pdf"){pdf(file=paste(outpath,"/Graphe_10Be.pdf",sep=""), width=w, height=h, pointsize=15, paper="a4r")
		}else if(graphs[i.graphe] =="png"){png(file=paste(outpath,"/Graphe_10Be.png",sep=""), width=w, height=h, pointsize=15)
		}

 		if(nc.vs.z.legend.position != 'ot'){par(mar=c(4,4,2,2)) } else {par(mar=c(4,4,2+Space.top.legend,2))}
		
		min.x <- min(CONCENTRATION.fin.Be,input.data$Be10[which(input.data$Be10>0)]); max.x <- max(CONCENTRATION.fin.Be,input.data$Be10[which(input.data$Be10>0)])
		plot(CONCENTRATION.fin.Be[,1],  -z.b.mat[,1],  type='p', xlab='Be-10 Concentration [at/g]', ylab='depth [m]', pch=20, col="light blue", cex=0.2, log='x', xlim=c(min.x,max.x), ylim=c(-z.b.lim.im,0))
		for(i.iteration in 1:n.iterations){lines(CONCENTRATION.fin.Be[,i.iteration], -z.b.mat[,i.iteration], col="light blue")}
		# best chi2
		if(Ne21.flag=="y"){lines(CONCENTRATION.fin.Be[,Best.Ne], -z.b.mat[,Best.Ne], col="red", lwd=2, lty=3)}
		if(Al26.flag=="y"){lines(CONCENTRATION.fin.Be[,Best.Al], -z.b.mat[,Best.Al], col="chartreuse3", lwd=2, lty=3)}
		for(i.iteration in Best.Be){lines(CONCENTRATION.fin.Be[,i.iteration], -z.b.mat[,i.iteration],col='blue', lwd=2)}
		if(length(Best.Be)==0){lines(CONCENTRATION.fin.Be[,Best.Be], -z.b.mat[,Best.Be],col='blue', lwd=2)}
		lines(CONCENTRATION.fin.Be[,Best.tout], -z.b.mat[,Best.tout], col="black", lwd=2)
		
		points(input.data$Be10, -input.data$Depth, pch=21, bg='blue', col='blue')
		for (i in (1:length(input.data$Be10))){
		  Vect_err.y <- Vect_err.x <- 0
		  Vect_err.y[2] <- Vect_err.y[1] <- input.data$Depth[i]; Vect_err.x[1] <- input.data$Be10[i]+input.data$Be10.err[i]; Vect_err.x[2] <- input.data$Be10[i]-input.data$Be10.err[i]
		  lines(Vect_err.x,-Vect_err.y, lwd=2, col='blue')
		}
		
		# Legend
		if(Ne21.flag=="y" && Al26.flag!="y" ){
			if(nc.vs.z.legend.position != 'ot'){
				legend(horiz=F, "bottomright", c("tested scenario", "best fit scenario", "best Be", "best Ne"), 
		           			col=c("light blue", "black", "blue", "red"), bg="white", lwd=c(1,2,2,2), lty=c(1,1,1,3)) 
			} else {
		           	legend(ncol=2, legend=c("tested scenario", "best fit scenario", "best Be", "best Ne"), 
		           			col=c("light blue", "black", "blue", "red"), bg="white", lwd=c(1,2,2,2), lty=c(1,1,1,3),
		           			y=par()$usr[4]*1.5, x=10^(par()$usr[1]), yjust = 0, xpd = TRUE) 
			}
		}
		if(Ne21.flag!="y" && Al26.flag=="y" ){
			if(nc.vs.z.legend.position != 'ot'){
				legend(horiz=F, "bottomright", c("tested scenario", "best fit scenario", "best Be", "best Al"), 
		           			col=c("light blue", "black", "blue","chartreuse3"), bg="white", lwd=c(1,2,2,2), lty=c(1,1,1,3)) 
			} else {
		           	legend(ncol=2, legend=c("tested scenario", "best fit scenario", "best Be",  "best Al"), 
		           			col=c("light blue", "black", "blue", "chartreuse3"), bg="white", lwd=c(1,2,2,2), lty=c(1,1,1,3),
		           			y=par()$usr[4]*1.5, x=10^(par()$usr[1]), yjust = 0, xpd = TRUE) 
			}
		}
		if(Ne21.flag=="y" && Al26.flag=="y" ){
			if(nc.vs.z.legend.position != 'ot'){
				legend(horiz=F, "bottomright", c("tested scenario", "best fit scenario", "best Be", "best Ne", "best Al"), 
		           			col=c("light blue", "black", "blue", "red","chartreuse3"), bg="white", lwd=c(1,2,2,2,2), lty=c(1,1,1,3,3)) 
			} else {
		           	legend(ncol=3, legend=c("tested scenario", "best fit scenario", "best Be", "best Ne", "best Al"), 
		           			col=c("light blue", "black", "blue", "red","chartreuse3"), bg="white", lwd=c(1,2,2,2,2), lty=c(1,1,1,3,3),
		           			y=par()$usr[4]*1.5, x=10^(par()$usr[1]), yjust = 0, xpd = TRUE) 
			}
		}
		
		if(graphs[i.graphe] =="png" |  graphs[i.graphe] =="pdf"){dev.off()}
	}
}


#########################################
##                                    ###
##    26Al Concentration vs depth     ###
##                                    ###
#########################################
if('NA' %in% Graphs.drawn){
	for(i.graphe in 1:length(graphs)){
		if(Al26.flag=="y"){
		if(graphs[i.graphe] =="pdf"){pdf(file=paste(outpath,"Graphe_26Al.pdf",sep=""), width=w, height=h, pointsize=15, paper="a4r")
		}else if (graphs[i.graphe] =="png"){png(file=paste(outpath,"Graphe_26Al.png",sep=""), width=w, height=h, pointsize=15)
		}

 		if(nc.vs.z.legend.position != 'ot'){par(mar=c(4,4,2,2)) } else {Prop.bas <- 2; par(mar=c(4,4,2+Space.top.legend,2))}
		  
		min.x <- min(CONCENTRATION.fin.Al,input.data$Al26[which(input.data$Al26>0)]); max.x <- max(CONCENTRATION.fin.Al,input.data$Al26[which(input.data$Al26>0)])
		  
		plot(CONCENTRATION.fin.Al[,1],  -z.b.mat[,1],  type='p', xlab='Al-26 Concentration [at/g]', ylab='depth [m]', pch=20, col="darkolivegreen1", cex=0.2, 
		              log='x', xlim=c(min.x,max.x), ylim=c(-z.b.lim.im,0))
		for(i.iteration in 1:n.iterations){lines(CONCENTRATION.fin.Al[,i.iteration], -z.b.mat[,i.iteration], col="darkolivegreen1")}
		
		# best chi2
		lines(CONCENTRATION.fin.Al[,Best.Be], -z.b.mat[,Best.Be], col="blue", lwd=2, lty=3)
		if(Ne21.flag=="y"){lines(CONCENTRATION.fin.Al[,Best.Ne], -z.b.mat[,Best.Ne], col="red", lwd=2, lty=3)}
		for(i.iteration in Best.Al){lines(CONCENTRATION.fin.Al[,i.iteration], -z.b.mat[,i.iteration],col='chartreuse3', lwd=2)}
		if(length(Best.Al)==0){lines(CONCENTRATION.fin.Al[,Best.Al], -z.b.mat[,Best.Al],col='chartreuse3', lwd=2)}
		lines(CONCENTRATION.fin.Al[,Best.tout], -z.b.mat[,Best.tout], col="black", lwd=2)
		
		points(input.data$Al26, -input.data$Depth, pch=21, bg='chartreuse3', col='chartreuse3')
		for (i in (1:length(input.data$Be10))){
		  Vect_err.y <- Vect_err.x <- 0
		  Vect_err.y[2] <- Vect_err.y[1] <- input.data$Depth[i]; Vect_err.x[1] <- input.data$Al26[i]+input.data$Al26.err[i]; Vect_err.x[2] <- input.data$Al26[i]-input.data$Al26.err[i]
		  lines(Vect_err.x,-Vect_err.y, lwd=2, col='chartreuse3')
		}
		
		# Legend
		if(Ne21.flag=="y"){
			if(nc.vs.z.legend.position != 'ot'){
				legend(horiz=F, "bottomright", c("tested scenario", "best fit scenario","best Al", "best Be", "best Ne"), 
		           			col=c("darkolivegreen1","black","chartreuse3", "blue", "red"), bg="white", lwd=c(1,2,2,2,2), lty=c(1,1,1,3,3)) 
			} else {
				legend(ncol=2,  c("tested scenario", "best fit scenario","best Al", "best Be", "best Ne"), 
		           			col=c("darkolivegreen1","black","chartreuse3", "blue", "red"), bg="white", lwd=c(1,2,2,2,2), lty=c(1,1,1,3,3),
		           			y=par()$usr[4]*1.5, x=10^(par()$usr[1]), yjust = 0, xpd = TRUE) 
			}
		} else {
			if(nc.vs.z.legend.position != 'ot'){
				legend(horiz=F, "bottomright", c("tested scenario", "best fit scenario","best Al", "best Be"), 
		           			col=c("darkolivegreen1","black","chartreuse3", "blue", "red"), bg="white", lwd=c(1,2,2,2), lty=c(1,1,1,3)) 
			} else {
		           	legend(ncol=3, legend=c("tested scenario", "best fit scenario", "best Be", "best Ne"), 
		           			col=c("light blue", "black", "blue", "red","chartreuse3"), bg="white", lwd=c(1,2,2,2), lty=c(1,1,1,3),
		           			y=par()$usr[4]*1.5, x=10^(par()$usr[1]), yjust = 0, xpd = TRUE) 
			}
		}
		
		if(graphs[i.graphe] =="png" |  graphs[i.graphe] =="pdf"){dev.off()}
		}
	}
}
	
#########################################
##                                    ###
##    21Ne Concentration vs depth     ###
##                                    ###
#########################################
if('NN' %in% Graphs.drawn){
	for(i.graphe in 1:length(graphs)){
		if(Ne21.flag=="y"){
		if(graphs[i.graphe] =="pdf"){pdf(file=paste(outpath,"/Graphe_21Ne.pdf",sep=""), width=w, height=h, pointsize=15, paper="a4r")
		}else if (graphs[i.graphe] =="png"){png(file=paste(outpath,"/Graphe_21Ne.png",sep=""), width=w, height=h, pointsize=15)
		}

 		if(nc.vs.z.legend.position != 'ot'){par(mar=c(4,4,2,2)) } else {Prop.bas <- 2; par(mar=c(4,4,2+Space.top.legend,2))}
		
		min.x <- min(CONCENTRATION.fin.Ne,input.data$Ne21[which(input.data$Ne21>0)]); max.x <- max(CONCENTRATION.fin.Ne,input.data$Ne21[which(input.data$Ne21>0)])
		
		plot(CONCENTRATION.fin.Ne[,1],  -z.b.mat[,1],  type='p', xlab='Ne-21 Concentration [at/g]', ylab='depth [m]', pch=20, col="pink", cex=0.2, log='x', xlim=c(min.x,max.x), ylim=c(-z.b.lim.im,0))
		for(i.iteration in 1:n.iterations){lines(CONCENTRATION.fin.Ne[,i.iteration], -z.b.mat[,i.iteration], col="pink")}
		
		# best chi2
		lines(CONCENTRATION.fin.Ne[,Best.Be], -z.b.mat[,Best.Be], col="blue", lty=3)
		if(Al26.flag=="y"){lines(CONCENTRATION.fin.Ne[,Best.Al], -z.b.mat[,Best.Al], col="chartreuse3", lty=3)}
		for(i.iteration in Best.Ne){lines(CONCENTRATION.fin.Ne[,i.iteration], -z.b.mat[,i.iteration],col='red')}
		if(length(Best.Ne)==0){lines(CONCENTRATION.fin.Ne[,Best.Ne], -z.b.mat[,Best.Ne],col='red')}
		lines(CONCENTRATION.fin.Ne[,Best.tout], -z.b.mat[,Best.tout], col="black", lwd=2)
		
		points(input.data$Ne21, -input.data$Depth, pch=21, bg='brown', col='brown')
		for (i in (1:length(input.data$Be10))){
		  Vect_err.y <- Vect_err.x <- 0
		  Vect_err.y[2] <- Vect_err.y[1] <- input.data$Depth[i]; Vect_err.x[1] <- input.data$Ne21[i]+input.data$Ne21.err[i]; Vect_err.x[2] <- max(input.data$Ne21[i]-input.data$Ne21.err[i],1)
		  lines(Vect_err.x,-Vect_err.y, lwd=2, col='brown')
		}
		
		# Legend
		if(Al26.flag=="y"){
			if(nc.vs.z.legend.position != 'ot'){
				legend(horiz=F, "bottomright", c("tested scenario", "best fit scenario", "best Ne", "best Be", "best Al"), 
		           		col=c("pink", "black", "red", "blue","chartreuse3"), bg="white", lwd=c(1,2,2,2,2), lty=c(1,1,1,3,3))  
			} else {
				legend(ncol=2,  c("tested scenario", "best fit scenario", "best Ne", "best Be", "best Al"), 
		           		col=c("pink", "black", "red", "blue","chartreuse3"), bg="white", lwd=c(1,2,2,2,2), lty=c(1,1,1,3,3),
		           		y=par()$usr[4]*1.5, x=10^(par()$usr[1]), yjust = 0, xpd = TRUE) 
			}
		} else {
			if(nc.vs.z.legend.position != 'ot'){
				legend(horiz=F, "bottomright", c("tested scenario", "best fit scenario", "best Ne", "best Be"), 
		           		col=c("pink", "black", "red", "blue"), bg="white", lwd=c(1,2,2,2), lty=c(1,1,1,3)) 
			} else {
		           	legend(ncol=3, legend=c("tested scenario", "best fit scenario", "best Ne", "best Be"), 
		           			col=c("pink", "black", "red", "blue"), bg="white", lwd=c(1,2,2,2), lty=c(1,1,1,3),
		           			y=par()$usr[4]*1.5, x=10^(par()$usr[1]), yjust = 0, xpd = TRUE) 
			}
		}
			
		if(graphs[i.graphe] =="png" |  graphs[i.graphe] =="pdf"){dev.off()}
		}
	}
}

#########################################
##                                    ###
##    Sedimentation rate computation  ###
##                                    ###
#########################################
if('SA' %in% Graphs.drawn || 'SD' %in% Graphs.drawn){
	# SedRate <- (z.b.mat[1:(length(z.b.mat[,1])-1),]-z.b.mat[-1,]) / (time.matrix[1:(length(time.matrix[,1])-1),]-time.matrix[-1,])
	z.b.mat.lim <- z.b.mat; z.b.mat.lim[which(z.b.mat>z.lim.im)] <- NA
	SedRate <- (z.b.mat.lim[1:(length(z.b.mat.lim[,1])-1),]-z.b.mat.lim[-1,]) / (time.matrix[1:(length(time.matrix[,1])-1),]-time.matrix[-1,])
	SedRate[which(is.na(z.b.mat.lim[1:(length(z.b.mat.lim[,1])-1),]))] <- NA
	
	time.mean <- (time.matrix[-1,]+time.matrix[-length(time.matrix[,1]),])/2  ### sediment rates are calculated from 2 neighboring time data
	z.mean <- (z.b.mat[-1,]+z.b.mat[-length(time.matrix[,1]),])/2             ### sediment rates are calculated from 2 neighboring time data
	SedRate.summary <- summary(matrix(na.omit(SedRate), ncol=1),method="table")
	SedRate.min <- as.numeric(unlist(strsplit(as.array(SedRate.summary)[1],"Min.   :"))[2])
	SedRate.max <- as.numeric(unlist(strsplit(as.array(SedRate.summary)[6],"Max.   :"))[2])
	if(SedRate.min==0){SedRate.min <- 1/10* as.numeric(unlist(strsplit(as.array(SedRate.summary)[2],"1st Qu.:"))[2])}
	if(SedRate.max > max.sr && s.vs.t.log =='n'){SedRate.max <- max.sr}
}

#########################################
##                                    ###
##    Sedimentation rate vs Age       ###
##                                    ###
#########################################
if('SA' %in% Graphs.drawn){
	### Rasterisation is only computed once
	### raster definition
	SedRateLog <- log10(SedRate)
	if(s.vs.t.log == 'y'){
		min.sedRateLog <- floor(log10(SedRate.min)); max.sedRateLog <- ceiling(log10(SedRate.max))
		e <- extent(c(0,n.step*t.step/1e6,min.sedRateLog,max.sedRateLog))
		r <- raster(e, nrows=sr.n.im, ncols=ceiling(n.step*t.step/t.step.im))
	}else{     ####  Linear
		min.sedRate <- SedRate.min; max.sedRate <- SedRate.max
		e <- extent(c(0,n.step*t.step/1e6,min.sedRate,max.sedRate))
		r <- raster(e, nrows=sr.n.im, ncols=ceiling(n.step*t.step/t.step.im)) 
	}
		
	### Interpolation from  points along each path at regularly spaced depths
	n.time.regular <- ceiling(n.step*t.step/t.step.im)
	sr.regular <- matrix(NA, nrow=n.time.regular, ncol=n.iterations)
	time.regular <- seq(res(r)[1]/2, -res(r)[1]/2 + res(r)[1]*n.time.regular,res(r)[1])
	if(s.vs.t.log == 'y'){
		for(i.z in 1:n.iterations){sr.regular[,i.z] <- approx(time.mean[,i.z]/1e6,SedRateLog[,i.z],time.regular)$y}   ####  Log
	}else{  for(i.z in 1:n.iterations){sr.regular[,i.z] <- approx(time.mean[,i.z]/1e6,SedRate[,i.z]   ,time.regular)$y}   ####  Linear
	}
	dim(sr.regular) <- n.time.regular*n.iterations
	time.regularM <- matrix(time.regular,ncol=n.iterations,nrow=n.time.regular,byrow =F)
	dim(time.regularM) <- n.time.regular*n.iterations
		
	### preparation raster of Chi2 values
	sr.time.chi.Be <- matrix(c(sr.regular,time.regularM,matrix(log10(Chi2.Be),ncol=n.iterations,nrow=n.time.regular,byrow =T)), ncol=3)
	if(Ne21.flag=="y"){sr.time.chi.Ne <- matrix(c(sr.regular,time.regularM,matrix(log10(Chi2.Ne),ncol=n.iterations,nrow=n.time.regular,byrow =T)), ncol=3)}
	if(Al26.flag=="y"){sr.time.chi.Al <- matrix(c(sr.regular,time.regularM,matrix(log10(Chi2.Al),ncol=n.iterations,nrow=n.time.regular,byrow =T)), ncol=3)}
	if(Al26.flag=="y" | Ne21.flag=="y"){
		sr.time.chi <- matrix(c(sr.regular,time.regularM,matrix(log10(Chi2.Tout),ncol=n.iterations,nrow=n.time.regular,byrow =T)), ncol=3)
		rast.chi <- rasterize(na.omit(sr.time.chi)[,2:1],r,(na.omit(sr.time.chi)[,3]), fun=function(x,...)min(x)) ### min
	} else {
		rast.chi <- rasterize(na.omit(sr.time.chi.Be)[,2:1],r,(na.omit(sr.time.chi)[,3]), fun=function(x,...)min(x)) ### min
	}
	rast.min <- summary(rast.chi)[1]; rast.q75 <- summary(rast.chi)[4]
				
	
	for(i.graphe in 1:length(graphs)){
		if(graphs[i.graphe] =="pdf"){pdf(file=paste(outpath,"/Graphe_sedimentation.pdf",sep=""), width=w, height=h, pointsize=15, paper="a4r")
		}else if(graphs[i.graphe] =="png"){png(file=paste(outpath,"/Graphe_sedimentation.png",sep=""), width=w, height=h, pointsize=15)
		}
		par(mar=c(4,4,2,5))
		
		### draw scenario probabilities
		### the following lines draw a graph in which chi2 is mapped
		### ----------------------------------------------
			
		### raster plot
		raster.coulour <- colorRampPalette(brewer.pal(name="YlOrBr",n=4))  ### https://moderndata.plot.ly/create-colorful-graphs-in-r-with-rcolorbrewer-and-plotly/
	 	if(s.vs.t.log == 'y'){####  Log
			image(rast.chi, ylim=c(min.sedRateLog ,max.sedRateLog), xlab='Time [Myrs]', ylab='log(Sedimentation rate [m/yr])', xaxp=c(0,100,50))
	 	}else{   ####  Linear
	 		image(rast.chi, ylim=c(min.sedRate ,   max.sedRate   ), xlab='Time [Myrs]', ylab='Sedimentation rate [m/yr]', xaxp=c(0,100,50))
	 	}
		if(m.type==1){plot(rast.chi, col='white', colNA='light grey', add=T, legend=FALSE)}else{plot(rast.chi, col='grey90', colNA='white',add=T, legend=FALSE)}
		
		if(m.type!=1){
			plot(rast.chi, zlim=c(0,Chi.Tout95),col='grey80', add=T, legend=FALSE)			
			plot(rast.chi, zlim=c(0,Chi.Tout68),col='grey70', add=T, legend=FALSE)			
		}					
		plot(rast.chi, col=rev(raster.coulour(255)), zlim=c(rast.min,rast.q75), add=T, alpha=0.75, 
					legend.args=list(text='log(Chi2)', side=4, line=2.5, cex=0.8))
	
		if(s.vs.t.levels=='y'){
			### extract simple probabilty limits
			sr.time.chi68 <- sr.time.chi.Be[which(sr.time.chi.Be[,3]<=log10(Chi.Be68)),]
			sr.time.chi95 <- sr.time.chi.Be[which(sr.time.chi.Be[,3]<=log10(Chi.Be95)),]
			max.chi68.t.Be <- min.chi68.t.Be <- max.chi95.t.Be <- min.chi95.t.Be <- NA
			for(i.z in 1:n.time.regular){ ### for eah z strip
				max.chi68.t.Be[i.z] <- max( na.omit(sr.time.chi68[ which(sr.time.chi68[,2] == time.regular[i.z]),1]))
				min.chi68.t.Be[i.z] <- min( na.omit(sr.time.chi68[ which(sr.time.chi68[,2] == time.regular[i.z]),1]))
				max.chi95.t.Be[i.z] <- max( na.omit(sr.time.chi95[ which(sr.time.chi95[,2] == time.regular[i.z]),1]))
				min.chi95.t.Be[i.z] <- min( na.omit(sr.time.chi95[ which(sr.time.chi95[,2] == time.regular[i.z]),1]))
				}
			
			### draw the curves for Be-10
			labcurve( list( One=list(time.regular,max.chi68.t.Be), Two=list(time.regular,min.chi68.t.Be), Three=list(time.regular,max.chi95.t.Be), Four= list(time.regular,min.chi95.t.Be)),
					keys=c('68%conf','','95%conf',''),  keyloc="none", lty=c("dashed","dashed","twodash","twodash"), lwd=1,
					col="blue", ###xlim=c(10,15), ### ylim=c(-100-20),
					pl=TRUE, add=T) 
			### idem for other nuclides
			if(Ne21.flag=="y"){
				sr.time.chi68 <- sr.time.chi.Ne[which(sr.time.chi.Ne[,3]<=log10(Chi.Ne68)),]
			 	sr.time.chi95 <- sr.time.chi.Ne[which(sr.time.chi.Ne[,3]<=log10(Chi.Ne95)),]
			 	max.chi68.t.Ne <- min.chi68.t.Ne <- max.chi95.t.Ne <- min.chi95.t.Ne <- NA
				for(i.z in 1:n.time.regular){ ### for eah z strip
					max.chi68.t.Ne[i.z] <- max( na.omit(sr.time.chi68[ which(sr.time.chi68[,2] == time.regular[i.z]),1]))
					min.chi68.t.Ne[i.z] <- min( na.omit(sr.time.chi68[ which(sr.time.chi68[,2] == time.regular[i.z]),1]))
					max.chi95.t.Ne[i.z] <- max( na.omit(sr.time.chi95[ which(sr.time.chi95[,2] == time.regular[i.z]),1]))
					min.chi95.t.Ne[i.z] <- min( na.omit(sr.time.chi95[ which(sr.time.chi95[,2] == time.regular[i.z]),1]))
					}
				labcurve( list( One=list(time.regular,max.chi68.t.Ne), Two=list(time.regular,min.chi68.t.Ne), Three=list(time.regular,max.chi95.t.Ne), Four= list(time.regular,min.chi95.t.Ne)),
				keys=c('68%conf','','95%conf',''),  keyloc="none", lty=c("dashed","dashed","twodash","twodash"), lwd=1,
				col="red", ###xlim=c(10,15), ### ylim=c(-100-20),
				pl=TRUE, add=T) 
				}
			if(Al26.flag=="y"){
				sr.time.chi68 <- sr.time.chi.Al[which(sr.time.chi.Al[,3]<=log10(Chi.Al68)),]
			 	sr.time.chi95 <- sr.time.chi.Al[which(sr.time.chi.Al[,3]<=log10(Chi.Al95)),]
			 	max.chi68.t.Al <- min.chi68.t.Al <- max.chi95.t.Al <- min.chi95.t.Al <- NA
				for(i.z in 1:n.time.regular){ ### for eah z strip
					max.chi68.t.Al[i.z] <- max( na.omit(sr.time.chi68[ which(sr.time.chi68[,2] == time.regular[i.z]),1]))
					min.chi68.t.Al[i.z] <- min( na.omit(sr.time.chi68[ which(sr.time.chi68[,2] == time.regular[i.z]),1]))
					max.chi95.t.Al[i.z] <- max( na.omit(sr.time.chi95[ which(sr.time.chi95[,2] == time.regular[i.z]),1]))
					min.chi95.t.Al[i.z] <- min( na.omit(sr.time.chi95[ which(sr.time.chi95[,2] == time.regular[i.z]),1]))
					}
				labcurve( list( One=list(time.regular,max.chi68.t.Al), Two=list(time.regular,min.chi68.t.Al), Three=list(time.regular,max.chi95.t.Al), Four= list(time.regular,min.chi95.t.Al)),
				keys=c('68%conf','','95%conf',''),  keyloc="none", lty=c("dashed","dashed","twodash","twodash"), lwd=1,
				col="chartreuse3", ###xlim=c(10,15), ### ylim=c(-100-20),
				pl=TRUE, add=T) 
				}
			### idem for multinuclides chi2.
			if(Al26.flag=="y" | Ne21.flag=="y"){
				sr.time.chi68 <- sr.time.chi[which(sr.time.chi[,3]<=log10(Chi.Tout68)),]
			 	sr.time.chi95 <- sr.time.chi[which(sr.time.chi[,3]<=log10(Chi.Tout95)),]
			 	max.chi68.t <- min.chi68.t <- max.chi95.t <- min.chi95.t <- NA
				for(i.z in 1:n.time.regular){ ### for eah z strip
					max.chi68.t[i.z] <- max( na.omit(sr.time.chi68[ which(sr.time.chi68[,2] == time.regular[i.z]),1]))
					min.chi68.t[i.z] <- min( na.omit(sr.time.chi68[ which(sr.time.chi68[,2] == time.regular[i.z]),1]))
					max.chi95.t[i.z] <- max( na.omit(sr.time.chi95[ which(sr.time.chi95[,2] == time.regular[i.z]),1]))
					min.chi95.t[i.z] <- min( na.omit(sr.time.chi95[ which(sr.time.chi95[,2] == time.regular[i.z]),1]))
					}
				labcurve( list( One=list(time.regular,max.chi68.t), Two=list(time.regular,min.chi68.t), Three=list(time.regular,max.chi95.t), Four= list(time.regular,min.chi95.t)),
				keys=c('68%conf','','95%conf',''),  keyloc="none", lty=c("dashed","dashed","twodash","twodash"), lwd=2,
				col="gray20", ###xlim=c(10,15), ### ylim=c(-100-20),
				pl=TRUE, add=T) 
				}
		}
			
		### best chi2
	 	if(s.vs.t.log == 'y'){####  Log
			if(Ne21.flag=="y"){lines(time.mean[,Best.Ne]/1e6, SedRateLog[,Best.Ne],col='red', lwd=2)}
			if(Al26.flag=="y"){lines(time.mean[,Best.Al]/1e6, SedRateLog[,Best.Al],col='chartreuse3', lwd=2)}
			lines(time.mean[,Best.Be]/1e6  , SedRateLog[,Best.Be]  ,col="blue", lwd=2)
			lines(time.mean[,Best.tout]/1e6, SedRateLog[,Best.tout],col="black", lwd=2)
	 	}else{   ####  Linear
			if(Ne21.flag=="y"){lines(time.mean[,Best.Ne]/1e6, SedRate[,Best.Ne],col='red', lwd=2)}
			if(Al26.flag=="y"){lines(time.mean[,Best.Al]/1e6, SedRate[,Best.Al],col='chartreuse3', lwd=2)}
			lines(time.mean[,Best.Be]/1e6  , SedRate[,Best.Be]  ,col="blue", lwd=2)
			lines(time.mean[,Best.tout]/1e6, SedRate[,Best.tout],col="black", lwd=2)
	 	}
		### caption
		if(Ne21.flag=="y"){
			if(Al26.flag=="y"){
				legend(horiz=F, "topleft", c("best fit scenario", "best Be", "best Ne", "best Al"), 
		              col=c("black", "dark blue","red", "chartreuse3"), bg="white", lwd=c(2,2,2,2))  ##, bty='n'
		    } else {
				legend(horiz=F, "topleft", c("best fit scenario", "best Be", "best Ne"), 
		              col=c("black", "dark blue","red"), bg="white", lwd=c(2,2,2))  ##, bty='n'
			}
		} else {	             
				legend(horiz=F, "topleft", c("best fit scenario", "best Be", "best Al"), 
		              col=c("black", "dark blue","chartreuse3"), bg="white", lwd=c(2,2,2))  ##, bty='n'
		}
	
		if(graphs[i.graphe] =="png" |  graphs[i.graphe] =="pdf"){dev.off()}
	}
}

#########################################
##                                    ###
##    Sedimentation rate vs depth     ###
##                                    ###
#########################################
if('SD' %in% Graphs.drawn){
	### Rasterisation is only computed once
	### raster definition
	depth.mean <- (z.b.mat[-1,]+z.b.mat[-length(z.b.mat[,1]),])/2
	if(s.vs.z.log == 'y'){
		min.sedRateLog <- floor(log10(SedRate.min)); max.sedRateLog <- ceiling(log10(SedRate.max))
		e <- extent(c(0,max(z.b),min.sedRateLog,max.sedRateLog))
		r <- raster(e, nrows=sr.n.im, ncols=ceiling(max(z.b)/z.step.im))
	}else{     ####  Linear
		min.sedRate <- SedRate.min; max.sedRate <- SedRate.max
		e <- extent(c(0,max(z.b),min.sedRate,max.sedRate))
		r <- raster(e, nrows=sr.n.im, ncols=ceiling(max(z.b)/z.step.im))  
	}
	### Interpolation from  points along each path at regularly spaced depths
	n.z.regular <- ceiling(max(z.b)/z.step.im)
	sr.regular  <- matrix(NA, nrow=n.z.regular, ncol=n.iterations)
	z.b.regular <- seq(res(r)[1]/2, -res(r)[1]/2 + res(r)[1]*n.z.regular,res(r)[1])
	if(s.vs.z.log == 'y'){
		for(i.z in 1:n.iterations){sr.regular[,i.z] <- approx(z.mean[,i.z],SedRateLog[,i.z],z.b.regular)$y}   ####  Log
	}else{  for(i.z in 1:n.iterations){sr.regular[,i.z] <- approx(z.mean[,i.z],SedRate[,i.z]   ,z.b.regular)$y}   ####  Linear
	}
	dim(sr.regular) <- n.z.regular*n.iterations
	z.b.regularM <- matrix(z.b.regular,ncol=n.iterations,nrow=n.z.regular,byrow =F)
	dim(z.b.regularM) <- n.z.regular*n.iterations
	
	### preparation raster of Chi2 values
	sr.t.chi.Be <- matrix(c(sr.regular,z.b.regularM,matrix(log10(Chi2.Be),ncol=n.iterations,nrow=n.z.regular,byrow =T)), ncol=3)
	if(Ne21.flag=="y"){sr.t.chi.Ne <- matrix(c(sr.regular,z.b.regularM,matrix(log10(Chi2.Ne),ncol=n.iterations,nrow=n.z.regular,byrow =T)), ncol=3)}
	if(Al26.flag=="y"){sr.t.chi.Al <- matrix(c(sr.regular,z.b.regularM,matrix(log10(Chi2.Al),ncol=n.iterations,nrow=n.z.regular,byrow =T)), ncol=3)}
	if(Al26.flag=="y" | Ne21.flag=="y"){
		sr.t.chi <- matrix(c(sr.regular,z.b.regularM,matrix(log10(Chi2.Tout),ncol=n.iterations,nrow=n.z.regular,byrow =T)), ncol=3)
		rast.chi <- rasterize(na.omit(sr.t.chi)[,2:1],r,(na.omit(sr.t.chi)[,3]), fun=function(x,...)min(x)) ### min
	} else {
		rast.chi <- rasterize(na.omit(sr.t.chi.Be)[,2:1],r,(na.omit(sr.t.chi)[,3]), fun=function(x,...)min(x)) ### min
	}
	rast.min <- summary(rast.chi)[1]; rast.q75 <- summary(rast.chi)[4]
	
	for(i.graphe in 1:length(graphs)){
		if(graphs[i.graphe] =="pdf"){pdf(file=paste(outpath,"/Graphe_sedimentation_depth.pdf",sep=""), width=w, height=h, pointsize=15, paper="a4r")
		}else if(graphs[i.graphe] =="png"){png(file=paste(outpath,"/Graphe_sedimentation_depth.png",sep=""), width=w, height=h, pointsize=15)
		}
		
		par(mar=c(4,4,2,5))
		
		### draw scenario probabilities
		### the following lines draw a graph in which chi2 is mapped
		### ----------------------------------------------
		 		
		### raster plot
		raster.coulour <- colorRampPalette(brewer.pal(name="YlOrBr",n=4))  ### https://moderndata.plot.ly/create-colorful-graphs-in-r-with-rcolorbrewer-and-plotly/
	 	if(s.vs.z.log == 'y'){####  Log
			image(rast.chi, ylim=c(min.sedRateLog ,max.sedRateLog), xlim=c(0,z.b.lim.im ), xlab='Depth [m]', ylab='log(Sedimentation rate [m/yr])')
	 	}else{   ####  Linear
			image(rast.chi, ylim=c(min.sedRate    ,max.sedRate)   , xlim=c(0,z.b.lim.im ), xlab='Depth [m]', ylab='Sedimentation rate [m/yr]')
	 	}
		if(m.type==1){plot(rast.chi, col='white', colNA='light grey', add=T, legend=FALSE)}else{plot(rast.chi, col='grey90', colNA='white',add=T, legend=FALSE)}
		### draw the copper bearing layer 
	 	if(s.vs.z.log == 'y'){####  Log
			if(length(copper.layer.depth)>1){ polygon(c(min(copper.layer.depth),min(copper.layer.depth),max(copper.layer.depth),max(copper.layer.depth)),c(min.sedRateLog ,max.sedRateLog,max.sedRateLog,min.sedRateLog), col="paleturquoise2", border=NA)}
	 	}else{   ####  Linear
			if(length(copper.layer.depth)>1){ polygon(c(min(copper.layer.depth),min(copper.layer.depth),max(copper.layer.depth),max(copper.layer.depth)),c(min.sedRate ,max.sedRate,max.sedRate,min.sedRate), col="paleturquoise2", border=NA)}
	 	}
		if(m.type!=1){
			plot(rast.chi, zlim=c(0,Chi.Tout95),col='grey80', add=T, legend=FALSE)			
			plot(rast.chi, zlim=c(0,Chi.Tout68),col='grey70', add=T, legend=FALSE)			
		}					
		plot(rast.chi, col=rev(raster.coulour(255)), zlim=c(rast.min,rast.q75), add=T, alpha=0.75, legend.args=list(text='log(Chi2)', side=4, line=2.5, cex=0.8))
	
		if(s.vs.z.levels=='y'){
			### extract simple probabilty limits
			sr.t.chi68 <- sr.t.chi.Be[which(sr.t.chi.Be[,3]<=log10(Chi.Be68)),]
			sr.t.chi95 <- sr.t.chi.Be[which(sr.t.chi.Be[,3]<=log10(Chi.Be95)),]
			max.chi68.t.Be <- min.chi68.t.Be <- max.chi95.t.Be <- min.chi95.t.Be <- NA
			for(i.z in 1:n.z.regular){ ### for eah z strip
				max.chi68.t.Be[i.z] <- max( na.omit(sr.t.chi68[ which(sr.t.chi68[,2] == z.b.regular[i.z]),1]))
				min.chi68.t.Be[i.z] <- min( na.omit(sr.t.chi68[ which(sr.t.chi68[,2] == z.b.regular[i.z]),1]))
				max.chi95.t.Be[i.z] <- max( na.omit(sr.t.chi95[ which(sr.t.chi95[,2] == z.b.regular[i.z]),1]))
				min.chi95.t.Be[i.z] <- min( na.omit(sr.t.chi95[ which(sr.t.chi95[,2] == z.b.regular[i.z]),1]))
				}
			
			### draw the curves for Be-10
			labcurve( list( One=list(z.b.regular,max.chi68.t.Be), Two=list(z.b.regular,min.chi68.t.Be), Three=list(z.b.regular,max.chi95.t.Be), Four= list(z.b.regular,min.chi95.t.Be)),
					keys=c('68%conf','','95%conf',''),  keyloc="none", lty=c("dashed","dashed","twodash","twodash"), lwd=1,
					col="blue", ###xlim=c(10,15), ### ylim=c(-100-20),
					pl=TRUE, add=T) 
			### idem for other nuclides
			if(Ne21.flag=="y"){
				sr.t.chi68 <- sr.t.chi.Ne[which(sr.t.chi.Ne[,3]<=log10(Chi.Ne68)),]
			 	sr.t.chi95 <- sr.t.chi.Ne[which(sr.t.chi.Ne[,3]<=log10(Chi.Ne95)),]
			 	max.chi68.t.Ne <- min.chi68.t.Ne <- max.chi95.t.Ne <- min.chi95.t.Ne <- NA
				for(i.z in 1:n.z.regular){ ### for eah z strip
					max.chi68.t.Ne[i.z] <- max( na.omit(sr.t.chi68[ which(sr.t.chi68[,2] == z.b.regular[i.z]),1]))
					min.chi68.t.Ne[i.z] <- min( na.omit(sr.t.chi68[ which(sr.t.chi68[,2] == z.b.regular[i.z]),1]))
					max.chi95.t.Ne[i.z] <- max( na.omit(sr.t.chi95[ which(sr.t.chi95[,2] == z.b.regular[i.z]),1]))
					min.chi95.t.Ne[i.z] <- min( na.omit(sr.t.chi95[ which(sr.t.chi95[,2] == z.b.regular[i.z]),1]))
					}
				labcurve( list( One=list(z.b.regular,max.chi68.t.Ne), Two=list(z.b.regular,min.chi68.t.Ne), Three=list(z.b.regular,max.chi95.t.Ne), Four= list(z.b.regular,min.chi95.t.Ne)),
				keys=c('68%conf','','95%conf',''),  keyloc="none", lty=c("dashed","dashed","twodash","twodash"), lwd=1,
				col="red", ###xlim=c(10,15), ### ylim=c(-100-20),
				pl=TRUE, add=T) 
				}
			if(Al26.flag=="y"){
				sr.t.chi68 <- sr.t.chi.Al[which(sr.t.chi.Al[,3]<=log10(Chi.Al68)),]
			 	sr.t.chi95 <- sr.t.chi.Al[which(sr.t.chi.Al[,3]<=log10(Chi.Al95)),]
			 	max.chi68.t.Al <- min.chi68.t.Al <- max.chi95.t.Al <- min.chi95.t.Al <- NA
				for(i.z in 1:n.z.regular){ ### for eah z strip
					max.chi68.t.Al[i.z] <- max( na.omit(sr.t.chi68[ which(sr.t.chi68[,2] == z.b.regular[i.z]),1]))
					min.chi68.t.Al[i.z] <- min( na.omit(sr.t.chi68[ which(sr.t.chi68[,2] == z.b.regular[i.z]),1]))
					max.chi95.t.Al[i.z] <- max( na.omit(sr.t.chi95[ which(sr.t.chi95[,2] == z.b.regular[i.z]),1]))
					min.chi95.t.Al[i.z] <- min( na.omit(sr.t.chi95[ which(sr.t.chi95[,2] == z.b.regular[i.z]),1]))
					}
				labcurve( list( One=list(z.b.regular,max.chi68.t.Al), Two=list(z.b.regular,min.chi68.t.Al), Three=list(z.b.regular,max.chi95.t.Al), Four= list(z.b.regular,min.chi95.t.Al)),
				keys=c('68%conf','','95%conf',''),  keyloc="none", lty=c("dashed","dashed","twodash","twodash"), lwd=1,
				col="chartreuse3", ###xlim=c(10,15), ### ylim=c(-100-20),
				pl=TRUE, add=T) 
				}
			### idem for multinuclides chi2.
			if(Al26.flag=="y" | Ne21.flag=="y"){
				sr.t.chi68 <- sr.t.chi[which(sr.t.chi[,3]<=log10(Chi.Tout68)),]
			 	sr.t.chi95 <- sr.t.chi[which(sr.t.chi[,3]<=log10(Chi.Tout95)),]
			 	max.chi68.t <- min.chi68.t <- max.chi95.t <- min.chi95.t <- NA
				for(i.z in 1:n.z.regular){ ### for eah z strip
					max.chi68.t[i.z] <- max( na.omit(sr.t.chi68[ which(sr.t.chi68[,2] == z.b.regular[i.z]),1]))
					min.chi68.t[i.z] <- min( na.omit(sr.t.chi68[ which(sr.t.chi68[,2] == z.b.regular[i.z]),1]))
					max.chi95.t[i.z] <- max( na.omit(sr.t.chi95[ which(sr.t.chi95[,2] == z.b.regular[i.z]),1]))
					min.chi95.t[i.z] <- min( na.omit(sr.t.chi95[ which(sr.t.chi95[,2] == z.b.regular[i.z]),1]))
					}
				labcurve( list( One=list(z.b.regular,max.chi68.t), Two=list(z.b.regular,min.chi68.t), Three=list(z.b.regular,max.chi95.t), Four= list(z.b.regular,min.chi95.t)),
				keys=c('68%conf','','95%conf',''),  keyloc="none", lty=c("dashed","dashed","twodash","twodash"), lwd=2,
				col="gray20", ###xlim=c(10,15), ### ylim=c(-100-20),
				pl=TRUE, add=T) 
				}
		}
	
		### best chi2
	 	if(s.vs.t.log == 'y'){####  Log
			if(Ne21.flag=="y"){lines(z.mean[,Best.Ne], SedRateLog[,Best.Ne],col='red', lwd=2)}
			if(Al26.flag=="y"){lines(z.mean[,Best.Al], SedRateLog[,Best.Al],col='chartreuse3', lwd=2)}
			lines(z.mean[,Best.Be]  , SedRateLog[,Best.Be]  ,col="blue", lwd=2)
			lines(z.mean[,Best.tout], SedRateLog[,Best.tout],col="black", lwd=2)
	 	}else{   ####  Linear
			if(Ne21.flag=="y"){lines(z.mean[,Best.Ne], SedRate[,Best.Ne],col='red', lwd=2)}
			if(Al26.flag=="y"){lines(z.mean[,Best.Al], SedRate[,Best.Al],col='chartreuse3', lwd=2)}
			lines(z.mean[,Best.Be]  , SedRate[,Best.Be]  ,col="blue", lwd=2)
			lines(z.mean[,Best.tout], SedRate[,Best.tout],col="black", lwd=2)
		}
			
		### caption
		if(Ne21.flag=="y"){
			if(Al26.flag=="y"){
				legend(horiz=T, "bottom", c("best fit scenario", "best Be", "best Ne", "best Al"), 
		              col=c("black", "dark blue","red", "chartreuse3"), bg="white", lwd=c(2,2,2,2))  ##, bty='n'
		    } else {
				legend(horiz=T, "bottom", c("best fit scenario", "best Be", "best Ne"), 
		              col=c("black", "dark blue","red"), bg="white", lwd=c(2,2,2))  ##, bty='n'
			}
		} else {	             
				legend(horiz=T, "bottom", c("best fit scenario", "best Be", "best Al"), 
		              col=c("black", "dark blue","chartreuse3"), bg="white", lwd=c(2,2,2))  ##, bty='n'
		}
		### draw sample depths
		rug(input.data$Depth,side=3,lwd=2)
		
		if(graphs[i.graphe] =="png" |  graphs[i.graphe] =="pdf"){dev.off()}
	} ## end for(i.graphe)...
}

