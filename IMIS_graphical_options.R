##_______________________________________________________________
##|                                                             |  
##|   IMIS INVERSION PROGRAM: GRAPHICAL OPTIONS                 |
##|_____________________________________________________________|
##
## Graphical options for graphic drawing (IMIS_graphs.R)
## Graphs drawn:
##    Depth vs age
##    Denudation rate vs time
##    Nuclide Concentration vs depth; Nuclides are 10Be (Blue), 21Ne (Red), 26Al(green)  
##    Sedimentation rate vs Age
##    Sedimentation rate vs depth


##########################################
###                                    ###
###    General parameters              ###
###                                    ###
##########################################

w <-1000; h <- w*2/3 ### Graph size and aspect

# specific options for graph depth vs time
# must be carefully chosen, because the drawing is time consuming.
d.vs.t.type <- 'm'       ### 'l' lines or 'm' map 
m.type <- 2             ### type of map. 2 possibilities: 1(grey background, white background explored paths) or 2 (no background, greys background explored paths, 3 levels of grey in function of outside 95% confidence, inside 95% and 68% confidence)	
t.step.im <- 250e3       ### mapping pixel size in yrs (Depth vs age and Sedimentation rate vs Age)
z.step.im <- 0.25        ### mapping pixel size in depth (Depth vs age and Sedimentation rate vs depth)
sr.n.im   <- 50         ### mapping pixel number in sed rate [advice: 50-200]
z.lim.im  <- 116          ### drawing limit, if different from max(z.b)
chi.scenario <- 'n' ### 'n' normal, chi, "i" chi inverse 
max.sr        <- NULL    # Default maximum sedimentation rate
max.sr        <- 1e-3    # manual(only for linear scale)


nc.vs.z.legend.position <- 'ibr' 		# Legend position 'ibr', inner [default] or 'ot' on top, outside the graph

#########################################
##                                    ###
##    Specific options                ###
##                                    ###
#########################################

##  Depth vs age  
##----------------
d.vs.t.levels <- 'n'     ### 'y' or 'n'. Draw levels (68% and 95% confidence for normal chi2, 0.8*max for inverse chi2)


##  Denudation rate vs age
##--------------------------
# currently no options


##  Nuclide Concentration vs depth
##----------------------------------
lim.x.Ne <- NULL # Default
lim.x.Al <- NULL # Default
#### manual limits can be set with the following lines
# lim.x.Ne <- c(2e6,5e7)                        #  extreme values for Ne (for graphs) 
# lim.x.Al <- c(1e4,10e7)                       # default extreme values for Al (for graphs) 



##  Sedimentation rate vs Age     
##------------------------------
s.vs.t.levels   <- 'n'     ### draw levels (68% and 95% confidence for normal chi2, 0.8*max for inverse chi2)
s.vs.t.log      <- 'n'   # logarithm scaling for sedimentation rate

##  Sedimentation rate vs depth     
##-------------------------------
s.vs.z.levels   <- 'n'     ### draw levels (68% and 95% confidence for normal chi2, 0.8*max for inverse chi2)
s.vs.z.log      <- 'n'   # logarithm scaling for sedimentation rate

