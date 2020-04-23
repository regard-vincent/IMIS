##_______________________________________________________________
##|                                                             |  
##|   IMIS                                                      |  
##|   INVERSION OF MULTI ISOTOPES IN SEDIMENTARY BASIN          |
##|_____________________________________________________________|
##
# # --------------------
# # Vincent REGARD 
# # --------------------
# # PRINCIPLE
# # This code is analysing cosmogenic nuclide content of rocks collected deep into a sedimentary basin
# # 10Be, 21Ne and/or 26Al concentrations can be processed
# # A rock in the sedimentary pile had a story divided into two stages
# # (Stage a) Erosion and instantaneous transport to the basin; Erosion rate is epsilon.a, Production Po.a is constant
# # (Stage b) The sample has been sedimented and covered by younger sediments. 
# # The sediment rate is random (probability law sr.b1) and the final total sediment thickness is fixed (depth of sample collection)
# # 
# # 
# # R command to run the code
# # first, set the working directory by the command:
# # setwd("D:/Vincent/Beryl/10Be_26Al/PB_Inverse/")
# # then, run the code
# # source('IMIS_main.R')
# #
# # current version, april 2019.
# #
# # DEPENDENCIES
### This R Program needs the following R files ( --> indicate dependencies)
### IMIS_Inversion.R
### IMIS_input.R
###  --> IMIS_graphical_options.R
### IMIS_functions.R
### IMIS_DataReading.R
### IMIS_cosmoparameters.R
###  -->  DepthProd_Heisinger_Balco2017_vr.R
### IMIS_graphs.R

### Record of overall time duration
time.a <-  Sys.time()

### REQUIRED libraries
installed.packages() -> test.lib
# if(length(which(test.lib=="lattice"))==0){install.packages("lattice", quiet=T)}
# library(lattice) 
if(length(which(test.lib=="tcltk"))==0){install.packages("tcltk", quiet=T)}
library(tcltk) 
if(length(which(test.lib=="plotrix"))==0){install.packages("plotrix", quiet=T)}
library(plotrix) 
if(length(which(test.lib=="Hmisc"))==0){install.packages("Hmisc", quiet=T)}
library(Hmisc) 
if(length(which(test.lib=="raster"))==0){install.packages("raster", quiet=T)}
library(raster) 

outpath <- getwd()

##########################################
###                                    ###
###     RUN DIFFERENT PROGRAM FILES    ###
###                                    ###
##########################################

source('IMIS_input.R'); max.t.b <- max(t.b)
source('IMIS_functions.R') 
source('IMIS_DataReading.R')
source('IMIS_cosmoparameters.R')

#### compatibility checks 
Repeat.step -> Repeat.step.max
if(max(t.a)<max(t.b)){t.a[which(t.a==max(t.a))] <- max(t.b)}

##########################################
###                                    ###
### Preparation of output text file    ###
###                                    ###
##########################################
if(Ne21.flag=="y"){
	if(Al26.flag=="y"){
		write.table(t(c("Chi2 10Be min","Chi2 21Ne min","Chi2 26Al min","Chi2 min","N# Be 68%","N# Be 95%","N# Ne 68%","N# Ne 95%","N# Al 68%","N# Al 95%","N# All 68%","N# All 95%","remark")), file=file.out, col.names=FALSE, row.names=FALSE, append=FALSE, sep="\t")
    } else {
		write.table(t(c("Chi2 10Be min","Chi2 21Ne min","Chi2 min","N# Be 68%","N# Be 95%","N# Ne 68%","N# Ne 95%","N# All 68%","N# All 95%","remark")), file=file.out, col.names=FALSE, row.names=FALSE, append=FALSE, sep="\t")
	}
} else {	             
		write.table(t(c("Chi2 10Be min","Chi2 26Al min","Chi2 min","N# Be 68%","N# Be 95%","N# Al 68%","N# Al 95%","N# All 68%","N# All 95%","remark")), file=file.out, col.names=FALSE, row.names=FALSE, append=FALSE, sep="\t")
}

info <- ''


###############################################
###                                         ###
### Choice                                  ###
###    Single inversion                     ###
###    Inversion for a set of parameters    ###
###                                         ###
###############################################

###   1- Single inversion
# source('IMIS_Inversion.R')
# if(length(graphs>0) && graphs!="null" && graphs!=""){source('IMIS_graphs.R')}

##   2- Automatisation for investigating a wide range for the same parameter
## exemple is given for El Tesoro, varying the source erosion rate
auto.er.a <- c(seq(from=50, to=500, by=50), seq(from=600, to=1000, by=100))              
opt.graph <- graphs <- c("png")
for(auto in auto.er.a){
	if(auto == 500){graphs <- c("png","pdf")}else{graphs <- opt.graph}
	epsilon.a <- c(auto,auto,10,10)*1e-6 ### El Tesoro
	info <- paste("eros rate=  ",auto) 
	outpath <- paste(getwd(),"/erate_",auto,"/",sep=''); if (!dir.exists(outpath)){dir.create(outpath)}
	source('IMIS_Inversion.R')
	if(length(graphs>0) && graphs!="null" && graphs!=""){source('IMIS_graphs.R')}
}

if(Be10.flag=="y"){write.table(data.frame(z=z.prod, Be=depth.prodBe), paste("depth.prodBe_",Prod.par,".txt", sep=""), sep=',', row.names=F, col.names=F)}
if(Al26.flag=="y"){write.table(data.frame(z=z.prod, Ne=depth.prodAl), paste("depth.prodAl_",Prod.par,".txt", sep=""), sep=',', row.names=F, col.names=F)}
if(Ne21.flag=="y"){write.table(data.frame(z=z.prod, Ne=depth.prodNe), paste("depth.prodNe_",Prod.par,".txt", sep=""), sep=',', row.names=F, col.names=F)}

time.b <- Sys.time(); print(paste("End program. Duration=",format(time.b-time.a)))
