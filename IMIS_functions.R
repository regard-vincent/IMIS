
####################
##                ##
##   IMIS Code    ##
##                ##
##   FUNCTIONS    ##
##                ##
####################

	
##_______________________________________________________________
##|                                                             |  
##|    Chi2                                                     |
##|_____________________________________________________________|
## 
## Calculates the Chi2 difference between model and data.
##
Chi2.data <- function(nature.z, nature.conc, nature.conc.err, model.conc, model.z){
	nb.data <- length(nature.conc)
	chi2 <- 0
	for(i.chi in 1:nb.data){
		if(nature.conc.err[i.chi]!=0){
			if(!is.numeric(nature.conc.err[i.chi]) || nature.conc.err[i.chi]==0){nature.conc.err[i.chi] <- 1}
			z1 <- min(model.z[which(model.z>nature.z[i.chi])]); z2 <- max(model.z[which(model.z<nature.z[i.chi])])
			if(z1==Inf){print(c("ERROR in the range of depth for Chi2 calculation (z=",nature.z[i.chi],")")); break  }
			if(z2==Inf){print(c("ERROR in the range of depth for Chi2 calculation (z=",nature.z[i.chi],")")); break  }
			dz1 <- (nature.z[i.chi]-z1)/(z2-z1) ; dz2 <- (nature.z[i.chi]-z2)/(z1-z2)
			chi2.i <- ((nature.conc[i.chi]-(model.conc[which(model.z==z1)]*dz2 + model.conc[which(model.z==z2)]*dz1))/nature.conc.err[i.chi])^2
			if(length(chi2.i)==0) {chi2 <- NA} else {chi2 <- chi2+chi2.i}
		}
	}
	chi2
	}

##_______________________________________________________________
##|                                                             |  
##|    Calculation of concentration during burial               |
##|_____________________________________________________________|
Concentration.Burial <- function(Time.f, PROD.f, rho.f, z.f, L.f){
  Conc.f <- matrix(0, nrow=length(z.f), ncol=length(Time.f))
  if(L.f !=0){
    for(i.f in (1:length(z.f))){
      for(j.f in (1:length(Time.f))){
        Conc.f[i.f,j.f] <- 1/L.f * sum(PROD.f*exp(-z.f[i.f]*rho.f/Att)*(1-exp(-L.f*Time.f[j.f]))) 
        }
      }
  } else {
    for(i.f in (1:length(z.f))){
      for(j.f in (1:length(Time.f))){
        Conc.f[i.f,j.f] <- sum(PROD.f*exp(-z.f[i.f]*rho.f/Att)*(Time.f[j.f])) 
      }
    }
  }
  Conc.f
  }


##_______________________________________________________________
##|                                                             |                
##|   Calculation of concentration for a steady state erosion   |
##|_____________________________________________________________|
Concentration.Erosion <- function(Eros.rate.f, PROD.f, rho.f, L.f){
	Conc.f <- seq(length=length(Eros.rate.f))*0
	for(i.f in (1:length(Eros.rate.f))){
		Conc.f[i.f] <- sum(PROD.f /(L.f + rho.f * Eros.rate.f[i.f]/Att)) 
		}
	return(Conc.f)
	}


##_______________________________________________________________
##|                                                             |  
##|    Calculation of concentration for a small time step       |
##|            delta.t at depth z.f                             |
##|_____________________________________________________________|
d.Concentration <- function(PROD.f, rho.f, z.f,delta.t){
	Conc.f <- seq(length=length(z.f))*0
	for(i.f in (1:length(z.f))){
		Conc.f[i.f] <- sum(PROD.f*exp(-rho.f/Att*z.f[i.f]))
		}
	return(Conc.f*delta.t)
	}

##_______________________________________________________________
##|                                                             |  
##|    Calculation of production (exponential mode) at depth z  |
##|_____________________________________________________________|
Production.z <- function(PROD_0.f, rho.f, z.f, attenuation){
	Prod.f <- PROD_0.f*exp(-rho.f/attenuation*z.f)
	return(Prod.f)
	}
	

##_______________________________________________________________
##|                                                             |  
##|    Generating erosion matrix                                      |
##|_____________________________________________________________|
Generate.erosion <- function(epsilon,timee,Law){
  erosmat <- matrix(data=NA,nrow=n.step)
  time.index <- ((n.step:1)-0.5)*t.step
  for(index in 1:n.step){
	  t1 <- min(timee[which(timee>time.index[index])])
	  if(Law=='step'){
	   	 erosmat[index]<-epsilon[which(timee==t1)]
  	  } else { 
	  	  if(Law!='cont'){
	  		print('Using default distribution: continuous')
  	  	  } 	
	  	  t2 <- max(timee[which(timee<time.index[index])])
	  	  e1 <- epsilon[which(timee==t1)] ; e2 <- epsilon[which(timee==t2)]	  
	  	  erosmat[index] <- ((e1*(time.index[index]-t2)) + (e2*(t1-time.index[index])))/(t1-t2)
  	  }
  }
  return(erosmat)
  }

Generate.eros.matrix <- function(){
  nb.elem.matrices <- n.step * n.iterations
  if(var.epsilon.a.law=='unif'){
  	MATRIX.epsilon.a <- as.matrix(1+ var.epsilon.a.percent/100*(Random.sr(var.epsilon.a.law,1,nb.elem.matrices)-0.5))
  } else if(var.epsilon.a.law=='chisq'){
  	median.obs <- median(Random.sr(var.epsilon.a.law,1,100000))
  	MATRIX.epsilon.a <- as.matrix(1+ var.epsilon.a.percent/100*((Random.sr(var.epsilon.a.law,1,nb.elem.matrices)/(2*median.obs))-0.5))
  }
  MATRIX.epsilon.a[which(MATRIX.epsilon.a<0)]<- -MATRIX.epsilon.a[which(MATRIX.epsilon.a<0)]    #negative values are impossible
  dim(MATRIX.epsilon.a)<-c(n.step,n.iterations)
  MATRIX.epsilon.a <- MATRIX.epsilon.a* as.vector(Generate.erosion(epsilon.a,t.a,var.a))
  return(MATRIX.epsilon.a)
}

    
##_______________________________________________________________
##|                                                             |  
##|    Generating t matrix                                      |
##|_____________________________________________________________|
Generate.tb <- function(){
	# ideal scenario
	t.step.ideal <- matrix(t.step,1,length(t.b)) 
	n.step.ideal <- matrix(1,1,length(t.b)-1)
	if(length(n.step.ideal)!=1){
		n.step.ideal[1:(length(n.step.ideal)-1)] <-  round((t.b[1:(length(n.step.ideal)-1)]-t.b[2:length(n.step.ideal)])/t.step)
		n.step.ideal[length(n.step.ideal)] <- round(max(t.b)/t.step)- sum(n.step.ideal[1:(length(n.step.ideal)-1)])
	}else{n.step.ideal[length(n.step.ideal)] <- round(max(t.b)/t.step)}

	# definition of the scenarii
	breaks.scenar <- t.b %*% matrix(1,1,n.iterations)+ runif(n.iterations*length(t.b.err), min=-1, max=1)*t.b.err 

   	# calculation of step duration
	t.step.scenar <- (breaks.scenar[1:(length(t.b)-1),] - breaks.scenar[2:length(t.b),]) /
    	   matrix(n.step.ideal,length(n.step.ideal),n.iterations)

	return(list(breaks=breaks.scenar,t.step=t.step.scenar,n.step=n.step.ideal))
	}


##_______________________________________________________________
##|                                                             |  
##|    Selecting distribution for sedimentation rates           |
##|_____________________________________________________________|
Random.sr <- function(Law,parameter,nb.values){
  if(Law=='unif'){
	   Sr.random <- runif(nb.values)
  } else { if(Law=='chisq'){
	  Sr.random <- rchisq(nb.values, df=parameter)
  } else {
	  print('Using default distribution: uniform')
	  Sr.random <- runif(nb.values)*2*median.value
  }}
  return(Sr.random)
  }

##_______________________________________________________________
##|                                                             |  
##|    sampling with ordering                                   |
##|_____________________________________________________________|
Sample.ordered <- function(number){
	a <- sample(number-1, size=Repeat.step-1) ## -1 because we want 3 stages between 2breaks
  	return(a[order(a)])
  }

