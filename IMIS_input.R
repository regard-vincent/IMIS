# # PROGRAM IMIS: main input
# # It is advised to prepare one such input file for each site
# # --------------------
# # Vincent REGARD 
# # --------------------
# #

############################
###    INPUT DATA FILE   ###
############################

file.in <- "ElTesoro.txt"; copper.layer.depth <- c(105.5,117.5)                # MAIN El Tesoro; copper.layer.depth only serves to colorize the depth


############################
###    OUTPUTS           ###
############################

## generic graphical parameters
graphs <- c("png") ### vector of graph types: "null", "R" (draw in R device), "pdf", and/or "png"  
            ### Note: may be long.
            ### CAREFUL: erase former files
Graphs.drawn <- c('DA', 'NB', 'NN', 'SA',  'SD', 'DT')   ### list of graphs drawn in 'DA' - Depth vs age; 'DT' - Denudation rate vs time; 
            ###    'NB', 'NN', 'NA' - 10Be, 21Ne, 26Al Concentration vs depth;
            ###    'SA' - Sedimentation rate vs Age; 'SD' - Sedimentation rate vs depth
            ###    'ALL' - all types of graphs drawn
progress.bar <- 'no'    ### 'txt' for text progress bar
                         ### 'win' for progress bar in external window
                         ### 'no' no progress bar [DEFAULT]
                         ### 'no' suppresses every information printed in device
            
### further graph options are tuned in IMIS_graphical_options.R
source('IMIS_graphical_options.R')

### File name for Chi2 outputs
file.out <- "Chi2_out.txt"

#################################################
### Cosmogenic nuclide parameters (general)   ###
#################################################

Prod.par <- "lu.br11"   # br11: Braucher et al. 2011
                     # gr : Granger and Muzikar, 2001
                     # lu.e: Same attenuation lengths and muons than Braucher et al. 2011 but with production for neutrons calculated from Lal/Stone
                     # lu.h: use of Heisinger Style written by M. Lupker
                     


############################################
### INVERSION first part: (a) [slopes]   ###
############################################
var.epsilon.a.percent <-0           # Variation in percent around central value
var.epsilon.a.law <- 'unif'         # Denudation rate distribution: either 'unif' (uniform) or 'chisq' (Chi2)  
epsilon.a <- c(500,500,10,10)*1e-6  # median denudation rate epsilon [m/yr], for each period
t.a  <- c(20,10,9.5,0)*1e6          # corresponding ages (yr)
var.a <- "cont"                     # Indicates how epsilon varies with time: 'cont' (continuous) or 'step' (step like)  
rho.a <- 2700                       # Source rock material density [kg/m3]


##############################################
### INVERSION second part: (b) [basin]     ###
##############################################
Latitude <- -23.5         	# Latitude [El Tesoro: 23.5° S ]
Elevation <- 2300       	# masl [El Tesoro: 2300]
rho.b <- 2300         		# deposit (b) material density [kg/m3]
t.step <- 5000 	                # time step for calculation [yr]
n.iterations <- 100           # number of tested scenarii  
var.sr.b.law <- 'unif'          # Sedimentation rate distribution: either 'unif' (uniform) or 'chisq' (Chi2)    
Repeat.step <- 7                # number of segment breaks between constraining points (see therefater)

# external constraints (constraining points)
z.b <- c(160,10,0)              # EL TESORO depths of known age t.b 
t.b <- c(20,9.5,0)*1e6	        # known ages[a]
t.b.err <- c(0,0.3,0)*1e6	# uncertainty on known ages[a]


