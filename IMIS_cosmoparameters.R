########################
##                    ##
##     COSMOGENICS    ## 
##     ENTRIES        ##
##                    ##
########################


### Desintegration and production constants ################
### ======================================= ################

L.Be <- 4.99e-7       		# 10Be desintegration constant [yr-1]; T1/2=1.39My
L.Al <- 9.8e-7         		# 26Al desintegration constant [yr-1]

Pslhl.Be <- 4.02   		# +/- 0.012  [at/g/yr] Borchers et al. 2016
R2610  <- 7.25			# Spallation production ratio 10Be/26Al. Borchers et al. 2016
Po.Ne.Be.Ra <- 4.23             # Ratio production 21Ne/10Be : 4.08+/-0.37[Balco et schuster, 2009]
		# 4.23+/-0.17 [Kober et al., 2011]


### Production parameters constants for various schemes  ################
### ==================================================== ################

### Attenuation length of Neutrons and various Muons [kg.m-2]
Att.gr <- c(1600 , 7380  , 26880 , 43600)     ## Granger and Muzikar, 2001
Att.br <- c(1600 , 15000 , 43200)             ## Braucher et al., 2011
### Productions ratios (total=1) ### Granger and Muzikar, 2001
PrRBe.gr <- c(NA,0.09/5 , 0.02/5 , 0.02/5 )      ## 10Be 
PrRAl.gr <- c(NA,0.72/30 ,0.16/30 ,0.19/30)      ## 26Al
### Productions Rates  ### Braucher et al., 2011
PrRBe.br <- c(NA,0.012  , 0.039  )               ## 10Be 
PrRAl.br <- c(NA,0.084  , 0.081  )               ## 26Al
### Neon production ratio Balco, Shuster, Blard, Zimmermann, Stone 2011
PrRNe.all<- c(NA, 0     , 0.036  )


### Calculation parameters   ################
### ======================== ################
z.step.P <- 0.01             # depth resolution of Production canevas [m]
### the caneva is a reference curve production vs depth. (Shortens calulation duration)
Po.Be.prop <- Po.Be.b <- Po.Be.a <- Att <- NULL; Po.Al.prop <- Po.Al.a <- Po.Al.b <- NULL #internal initiation







##_______________________________________________________________________________
##|                                                                              |  
##|    Function surface production followin LAL1991/STONE2000 SCALING factors    |
##|    Scaling coefficients as a function of latitude (Stone et al, 2000 - JGR)  |
##|______________________________________________________________________________|
Production.LalStone <- function(Elevation.f, Latitude.f){
	Lat.lalstone 	<- c(0,10,20,30,40,50,60);			
	A.lalstone		<- c(31.8518,34.3699,40.3153,42.0983,56.7733,69.0720,71.8733);
	B.lalstone		<- c(250.3193,258.4759,308.9894,512.6857,649.1343,832.4566,863.1927);
	C.lalstone		<- c(-0.083393,-0.089807,-0.106248,-0.120551,-0.160859,-0.199252,-0.207069);
	D.lalstone		<- c(7.4260e-5,7.9457e-5,9.4508e-5, 1.1752e-4,1.5463e-4,1.9391e-4,2.0127e-4);
	E.lalstone		<- c(-2.2397e-8,-2.3697e-8,-2.8234e-8,-3.8809e-8,-5.0330e-8,-6.3653e-8,-6.6043e-8);
	M.lalstone		<- c(0.587,0.600,0.678,0.833,0.933,1.000,1.000);
	PRES			<- 1013.25*exp(((-0.03417)/6.5e-3)*(log(288.15)-log(288.15-(6.5e-3*Elevation.f)))); #Mean atmospheric pressure
	a.lalstone		<- approx(Lat.lalstone,A.lalstone,abs(Latitude.f))$y;				#Interpolation of Stones empirical parameters according to input lat
	b.lalstone		<- approx(Lat.lalstone,B.lalstone,abs(Latitude.f))$y;
	c.lalstone		<- approx(Lat.lalstone,C.lalstone,abs(Latitude.f))$y;
	d.lalstone		<- approx(Lat.lalstone,D.lalstone,abs(Latitude.f))$y;
	e.lalstone		<- approx(Lat.lalstone,E.lalstone,abs(Latitude.f))$y;
	m.lalstone		<- approx(Lat.lalstone,M.lalstone,abs(Latitude.f))$y;
              		
	Fn_St <- (a.lalstone+b.lalstone*exp(-PRES/150)+c.lalstone*PRES+d.lalstone*PRES^2+e.lalstone*PRES^3)
	Fm_St <- m.lalstone*exp((1013.25-PRES)/242)

	return(c(Fn_St,Fm_St))
	}
	


### PRODUCTION PARAMETERS FOR BE or AL  [NEON IS PROCESSED HEREUNDER] ###
### ================================================================= ###

                    
### Computation of Lal / Stone coefficients in function of site latitude and elevation 
Production.LalStone(Elevation, Latitude)-> Prod.LalStone; Fn_St <- Prod.LalStone[1]; Fm_St <- Prod.LalStone[2]
Po.Be.a[1]  <- Pslhl.Be * Fn_St   		# Be production rate (neutrons) at earth's surface [at/g/yr], at source elevation (a)
Po.Al.a[1]  <- Pslhl.Be * Fn_St * R2610		# Al production rate (neutrons) at earth's surface [at/g/yr], at source elevation (a)

                 
### GRANGER
	if(Prod.par=="gr"){			# granger and Muzikar, 2001
	Att <- Att.gr
	# production SLHL (ratios)
	Po.Be.prop    <- PrRBe.gr
	Po.Be.prop[1] <- 1-sum(Po.Be.prop[2:4])        		
	Po.Al.prop    <- PrRAl.gr
 	Po.Al.prop[1] <- 1-sum(Po.Al.prop[2:4])        		

	# Lal /Stone Scaling
	Po.Be.a <- c(Fn_St, Fm_St , Fm_St, Fm_St )* Po.Be.prop * Pslhl.Be
	Po.Al.a <- c(Fn_St, Fm_St , Fm_St, Fm_St )* Po.Al.prop * Pslhl.Be
	
	Po.Be.prop <- Po.Be.a/sum(Po.Be.a)        		# calculated production ratio from different particles, 10Be
	Po.Al.prop <- Po.Al.a/sum(Po.Al.a)        		# calculated production ratio from different particles, 26Al
	}

	
## BRAUCHER
### BRAUCHER et al. 2011 (written after M. Lupker)
	if(Prod.par=="br11"|| Prod.par=="lu.h"){         
	Att <- Att.br
	Po.Be   <- PrRBe.br; Po.Al   <- PrRAl.br; 
	# production SLHL
	Po.Be.a <- c(Pslhl.Be         * Fn_St , Po.Be[2:3] * Fm_St )
	Po.Al.a <- c(Pslhl.Be * R2610 * Fn_St , Po.Al[2:3] * Fm_St )
	
	Po.Be.prop <- Po.Be.a/sum(Po.Be.a)        		# calculated production ratio from different particles, 10Be
	Po.Al.prop <- Po.Al.a/sum(Po.Al.a)        		# calculated production ratio from different particles, 26Al
	}

Po.Be.prop[1] <- 1 - sum(Po.Be.prop[-1])
	

### LUPKER
if(substr(Prod.par,1,2)=="lu"){
	source('DepthProd_Heisinger_Balco2017_vr.R') ### function written by M. Lupker to determine production at different depth
	if(Prod.par=="lu.e"){ Att <- Att.br }     			# this is simple case

	Po.Be.a <- as.data.frame(Production(0,abs(Latitude),Elevation,"Be",substr(Prod.par,4,4)))
	Po.Al.a <- as.data.frame(Production(0,abs(Latitude),Elevation,"Al",substr(Prod.par,4,4)))

	Po.Be.prop <- Po.Be.a/sum(Po.Be.a)        		# calculated production ratio from different particles, 10Be
	Po.Al.prop <- Po.Al.a/sum(Po.Al.a)        		# calculated production ratio from different particles, 26Al
}



##### PRODUCTION PARAMETERS FOR NEON ###
### ================================ ###

Po.Ne.prop <- Po.Ne.a <- Po.Ne.b <- NULL

### PROD PARAMETERS Balco, Shuster, Blard, Zimmermann, Stone 2011
if(Prod.par!="br11"){ 		
	Po.Ne.prop <- PrRNe.all
	Po.Ne.prop[1] <- 1 - sum(Po.Ne.prop[-1])	
}else{	Po.Ne.prop[1:2] <- c(1 - PrRNe.all[3], PrRNe.all[3])	}

Po.Ne.a  <- sum(Po.Be.a) * Po.Ne.Be.Ra * Po.Ne.prop




#### 
### FINALISATION AND HOMOGENEISATION OF DATA ################
### ======================================== ################
if(is.null(Po.Al.prop)){Po.Al.prop <- Po.Be.prop}
if(is.null(Po.Ne.prop)){Po.Ne.prop <- Po.Be.prop}

### to take into account the diversity of muons
### we verify length of attenuation parameters corresponds to that of prod parameters
### for example CRONUS gives only one muon prod rate whereas most production is using 2 muons with 2 attenuation lengths
if(Prod.par!="lu.h"){
	if(length(Att)!=length(Po.Be.a)){Po.Be.a[2:length(Att)] <- sum(Po.Be.a[-1]) * Po.Be.prop[-1]/sum(Po.Be.prop[-1])}
	if(length(Att)!=length(Po.Al.a)){Po.Al.a[2:length(Att)] <- sum(Po.Al.a[-1]) * Po.Al.prop[-1]/sum(Po.Al.prop[-1])}
	if(length(Att)!=length(Po.Ne.a)){suppressWarnings(Po.Ne.a[2:length(Att)] <- sum(Po.Ne.a[-1]) * Po.Ne.prop[-1]/sum(Po.Ne.prop[-1]))}
	if(progress.bar != 'no') {
		print(c("Prod (spallation, muons): Be=", Po.Be.a))
		if(Al26.flag=="y"){print(c("Prod (spallation, muons): Al=",Po.Al.a))}
		if(Ne21.flag=="y"){print(c("Prod (spallation, muons): Ne=",Po.Ne.a))}
		print(c("Attn.(spallation, muons)=", Att))
	}

}

Po.Be.b  <- Po.Be.a ; Po.Al.b  <- Po.Al.a ; Po.Ne.b  <- Po.Ne.a    ### Default is PROD(b) = PROD (a)
