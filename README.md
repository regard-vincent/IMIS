<h1>    IMIS
	INVERSION OF MULTI ISOTOPES IN SEDIMENTARY BASIN
</h1>

</body>

<h2> 
        DESCRIPTION
</h2>

<b> PRINCIPLE </b> 

This a code running with R https://www.r-project.org/
This code is analysing cosmogenic nuclide content of rocks collected deep into a sedimentary basin
10Be, 21Ne and/or 26Al concentrations can be processed
A rock in the sedimentary pile had a story divided into two stages
(Stage a) Erosion and instantaneous transport to the basin; Erosion rate is epsilon.a, Production Po.a is constant
(Stage b) The sample has been sedimented and covered by younger sediments. 
The sediment rate is random (probability law sr.b1) and the final total sediment thickness is fixed (depth of sample collection)

<b> MAIN INPUT </b> 

Text file containing the following columns (separated by either space or tab): 
sample	z	+/-	10Be	+/-	21Ne	+/-	26Al	+/- 
sample: sample label
z: sample depth (m)
10Be, 21Ne, 26Al: cosmogenic nuclide concentration (at/g)
+/- refer to uncertainty for the previous column value

Code setup parameters are found in IMIS_input.R, IMIS_cosmoparameters.R, and IMIS_graphical_options.R

OUTPUTS
Various graphics: 1- Depth vs age; 2- Denudation rate vs time; 3- Nuclide Concentration vs depth;
      4- Sedimentation rate vs Age; 5- Sedimentation rate vs depth
Chi2 misfit values in output file.
Production vs. depth. For example file depth.prodBe_br11.txt contains two columns, the first is depth (m), the second the total production for 10Be (at/g)
      
      
_____________________________
|                            |
|   HOW TO RUN THE CODE      |
|   AND CODE WORKFLOW        |
|____________________________|

R command to run the code
first, set the working directory by the command setwd(YOUR WORKING DIRECTORY PATH HERE). Ex: setwd("D:/Vincent/Beryl/10Be_26Al/PB_Inverse/")
then, run the code source('IMIS_main.R').

The main input parameters can be tuned in IMIS_input.R
The main graphical options can be set in IMIS_graphical_options.R
Cosmogenic parameters can be changed in IMIS_cosmoparameters.R

IMIS WORKFLOW
The IMIS_main.R code calls the following codes (described in the following):
IMIS_input.R            | code containing the input parameters
IMIS_functions.R        | code containing various functions
IMIS_cosmoparameters.R  | code containing cosmogenic isotopes inputs and calculating cosmogenic production parameters. It can call DepthProd_Heisinger_Balco2017_vr.R
IMIS_DataReading.R      | code reading input data 
IMIS_Inversion.R        | main calculation code
IMIS_graphical_options.R| graphical options
IMIS_graphs.R           | code for graphical outputs 


____________________________________________
|                                           |
|   DESCRIPTION OF THE DIFFERENT CODES      |
|___________________________________________|

IMIS_main.R
-----------
Calls the different parts of the code.
It allows automatisation. In part "Choice Single inversion / Inversion for a set of parameters", release # signs in part 2
    This will create different folders to put output graphs.
    
IMIS_input.R. Contains the main input parameters
file.in                 | input file
copper.layer.depth      | indicates the interval colored in blue. Leave blank for no coloration
graphs                  | indicates what kind of file is written for graphs (possible: inside R, pdf or png)  
Graphs.drawn            | indicates the different graphs to be drawn (Depth vs age; Denudation rate vs time; Cosmogenic nuclide concentration vs depth; Sedimentation rate vs Age; Sedimentation rate vs depth)
file.out                | output file
Prod.par                | Cosmogenic scheme: br11 from Braucher et al. 2011, gr from Granger and Muzikar, 2001, lu.e: Lal/Stone surface production and Braucher et al. 2011 for attenuation lengths; lu.h: Lal/Stone surface production and Heisinger attenuation
var.epsilon.a.percent   | Slope (=stage a) variation in denudation rate
var.epsilon.a.law       | Random law for denudation rates in slopes (stage a): uniform or Chi squared laws  
t.a                     | Stage a/slope time step limits (yr)
epsilon.a               | Stage a/slope median denudation rates (m/yr)
var.a                   | Defines how denudation varies with time: continuous or step like  
rho.a                   | Source rock material density (kg/m3)
Latitude                | Site latitude
Elevation               | Site elevation (masl)
rho.b                   | Deposit material density (kg/m3)
t.step                  | Calculation time (yr)
n.iterations            | Number of tested scenarii  
var.sr.b.law            | Random law for sedimentation rate (stage b): uniform or Chi squared laws     
z.b                     | Basin (stage b) depths of known age t.b 
t.b                     | Stage b/basin time step limits (yr) 
t.b.err                 | uncertainty on known ages (yr)
Repeat.step             | Number of breaks in slope in between the t.b time step limits


IMIS_functions.R
----------------
code containing various functions


IMIS_cosmoparameters.R
----------------------
Contains Cosmo parameters and prepares a production vs depth dataset (matrix with a column for depth, the other for production)
In the following Production will not be recalculated but evaluated by linear interpolation of this dataset.
Calls DepthProd_Heisinger_Balco2017_vr.R for lu.e and lu.h
Main parameters
Pslhl.Be                | 10Be production at High Latitude Sea Level (not used for lu.e and lu.h)
R2610                   | Spallation production ratio 10Be/26Al
Po.Ne.Be.Ra             | Production Ratio 21Ne/10Be
Att.gr and Att.br       | Attenuation lengths (g/cm2)    
PrRBe, PrRAl, PrRNe     | production ratios or rates, depending on the scheme
z.step.P                | Resolution in depth of the production dataset (m)


IMIS_DataReading.R
------------------
Reads and setup the input data


IMIS_Inversion.R
----------------
This is the main calculation code. Each stage between two values of t.b is divided into Repeat.step time intervals with random sedimentation rate and duration.
For each time step, the cosmogenic nuclide concentration is computed following the production, depth and desintegration.
This is repeated n.iterations times (scenarios). Each scenario is characterized by the misfit between the record and the scenario given by the Chi-2 metric
The Chi-2 is calculated for each nuclide and with the different nuclide data.
The lowest Chi-2 is extracted (best fit)


IMIS_graphical_options.R
------------------------
Options for graphics
w, h                   | Graph size (external output)
d.vs.t.type            | Graph type. Lines: each scenari is draww, colored in function of Chi-2. Map: Chi-2 is rasterized; this raster is drawn 
m.type                 | 2 different kinds of raster maps
t.step.im              | For raster maps. Pixel size in yrs (Depth vs age and Sedimentation rate vs Age)
z.step.im              | For raster maps. Pixel size in depth (Depth vs age and Sedimentation rate vs depth)
sr.n.im                | For raster maps. Pixel number in sedimentation rate
z.lim.im               | For raster maps. Drawing limit in depth, if different from max(z.b)
chi.scenario           | Chi62 can be drawn or its inverse  
max.sr                 | Maximum sedimentation rate in graphs
d.vs.t.levels          | Depth vs age graph. Sets whether to draw levels (68% and 95% confidence for normal chi2, 0.8*max for inverse chi2)
nc.vs.z.legend.position| Nuclide Concentration vs depth. Legend position
lim.x.Ne and lim.x.Al  | Nuclide Concentration vs depth. Maximum concentration values in graphs
s.vs.t.levels          | Sedimentation rate vs Age. Sets whether to draw levels (68% and 95% confidence for normal chi2, 0.8*max for inverse chi2)
s.vs.t.log             | Sedimentation rate vs Age. Sets whether to use logarithm scaling for sedimentation rate
s.vs.z.levels          | Sedimentation rate vs depth. Sets whether to draw levels (68% and 95% confidence for normal chi2, 0.8*max for inverse chi2)
s.vs.z.log             | Sedimentation rate vs depth. Sets whether to use logarithm scaling for sedimentation rate                                  


IMIS_graphs.R          
-------------
Code for graphical outputs 
The graphs drawn are set by the option Graphs.drawn (IMIS_input.R)
Whether the graphs are written in files (pdf or png), or in R device is set by option graphs (IMIS_input.R)
Names of graphs:
Graphe_scenarios:           graph containing the depth vs time history and most probable (minimum Chi2 calculated on all nuclides) depth vs time paths
	In line mode, the darker the lines, the more probable the scenarios (Chi-squared value calculated on all nuclides)
	In map mode, a raster map indicating the most probable path crossing each pixel is drawn (Chi2 value). 
		The paths which are in the 25% less probable group are discarded.
		If m.type=2 (IMIS_graphical_options.R), the grey levels are drawn according to 68% and 95% confidence intevals.
	If chi.scenario is inverse, the values drawn are (minimum value of Chi_2)/Chi_2 ; thus values are in the range 0-1.
Graphe_erosion:             graph containing the denudation vs time history and most probable paths (see Graphe_scenarios for description)
Graphe_10Be, Graphe_26Al, Graphe_21Ne:  graphs containing the scenarios' concentration vs depth. The best fit scenarios are drawn. 
	In black line is the best fit considering all the cosmogenic nuclide data; in blue, the best 10Be fit; in green, the best 26Al fit; in red, the best 21Ne fit.
Graphe_sedimentation:       graph containing the Sedimentation rate vs time history and most probable paths (see Graphe_scenarios for description)
Graphe_sedimentation_depth: graph containing the Sedimentation rate vs depth history and most probable paths (see Graphe_scenarios for description)

