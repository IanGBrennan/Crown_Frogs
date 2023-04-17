#devtools::install_github(repo="nmatzke/BioGeoBEARS", force=T)
source("/Users/ianbrennan/Documents/GitHub/Crown_Frogs/Scripts/plot_BGB.R")

#######################################################
# SETUP -- libraries/BioGeoBEARS updates
#######################################################

# Load the package (after installation, see above).
library(optimx)   # optimx seems better than R's default optim()
library(GenSA)    # GenSA seems better than optimx (but slower) on 5+ parameters, 
# seems to sometimes fail on simple problems (2-3 parameters)
library(FD)       # for FD::maxent() (make sure this is up-to-date)
library(snow)     # (if you want to use multicore functionality; some systems/R versions prefer library(parallel), try either)
library(parallel)

#######################################################
# 2018-10-10 update: I have been putting the 
# updates on CRAN/GitHub
# You should use:
# rexpokit version 0.26.6 from CRAN
# cladoRcpp version 0.15 from CRAN
# BioGeoBEARS version 1.1 from GitHub, install with:
# library(devtools)
# devtools::install_github(repo="nmatzke/BioGeoBEARS")
#######################################################
library(rexpokit)
library(cladoRcpp)
library(BioGeoBEARS)


setwd("/Users/ianbrennan/Documents/GitHub/Crown_Frogs/Biogeography")
extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
list.files(extdata_dir)


##########################
# Time Stratified Analysis
##########################
# Note for interpreting 'x'
# If x is free, and inferred to be negative (which is usually the case), 
# then dispersal probability drops off with distance. If x=-1, then doubling 
# the distance halves the dispersal probability.  If x=-2, then dispersal 
# probability drops off with the square of distance.

# Tree File
trfn = "CrownFrogs_BGB_Fossils.tre"
tr <- read.tree(trfn)
# Geography File, specify which Hypothesis!
geogfn = "CrownFrogs_BGB_H2.txt"

# Set your tip ranges
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
# Maximum range size observed:
#max(rowSums(dfnums_to_numeric(tipranges@df)))
max_range_size = 3



#######################################################
# Run DEC
#######################################################

# Intitialize a default model (DEC model)
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
# Give BioGeoBEARS the location of the phylogeny Newick file
BioGeoBEARS_run_object$trfn = trfn
# Give BioGeoBEARS the location of the geography text file
BioGeoBEARS_run_object$geogfn = geogfn
# Input the maximum range size
BioGeoBEARS_run_object$max_range_size = max_range_size

# Min branch length to be considered a tip (below threshold is an ancestor)
BioGeoBEARS_run_object$min_branchlength = 0.000001 
# Allow null ranges
BioGeoBEARS_run_object$include_null_range = TRUE

# Provide your time slices
BioGeoBEARS_run_object$timesfn = "CrownFrogs_timeperiods.txt"
# Provide your distance matrices
BioGeoBEARS_run_object$distsfn = "CrownFrogs_distances_scaled.txt"

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 8    # use multiple cores in parallel

# Sparse matrix for huge state spaces (600+)
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE, cut_fossils=FALSE)
# The stratified tree is described in this table:
#BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Run this to check inputs. Read the error messages if you get them!
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

# For a slow analysis, run once, then set runslow=FALSE to just 
# load the saved result.
runslow = TRUE
resfn = "CrownFrogs_DEC_H2.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  save(res, file=resfn)
  resDEC = res
} else {
  # Loads to "res"
  load(resfn)
  resDEC = res
}

# Plot the time-stratified DEC results
analysis_titletxt ="DEC CrownFrogs H2"
# Setup
results_object = resDEC
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
# States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", 
                                label.offset=0.45, tipcex=0.3, statecex=0.3, splitcex=0.3, titlecex=0.5, 
                                plotsplits=FALSE, cornercoords_loc=scriptdir, include_null_range=TRUE, 
                                tr=tr, tipranges=tipranges)
# Pie chart
plot_BGB(results_object, analysis_titletxt, plotwhat="pie", 
         label.offset=0.45, tipcex=0.3, statecex=0.5, splitcex=0.001, titlecex=0.5, 
         plotsplits=FALSE, cornercoords_loc=scriptdir, include_null_range=TRUE, 
         tr=tr, tipranges=tipranges, 
         #skiplabels = T,
         #simplify_piecharts = F,
         tipboxes_TF = T,
         pie_tip_statecex = 0.3,
         plotlegend = T, legend_cex = 0.2)


#######################################################
# Run DEC+J
#######################################################
# for notes on each command, see the DEC model above
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
BioGeoBEARS_run_object$timesfn = "CrownFrogs_timeperiods.txt"
BioGeoBEARS_run_object$distsfn = "CrownFrogs_distances_scaled.txt"
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 8
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE, cut_fossils=FALSE)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
# Set up DEC+J model
# Get the ML parameter values from the 2-parameter nested model
# (this will ensure that the 3-parameter model always does at least as good)
dstart = resDEC$outputs@params_table["d","est"]
estart = resDEC$outputs@params_table["e","est"]
jstart = 0.0001
# Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart
# Add j as a free parameter
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart
# Check all is good
check_BioGeoBEARS_run(BioGeoBEARS_run_object)
# Run DEC+J
resfn = "CrownFrogs_DEC_j_H2.Rdata"
runslow = TRUE
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  
  resDECj = res
} else {
  # Loads to "res"
  load(resfn)
  resDECj = res
}

# Plot the time-stratified DEC results
analysis_titletxt ="DEC+J CrownFrogs H2"
# Setup
results_object = resDECj
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
# States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", 
                                label.offset=0.45, tipcex=0.3, statecex=0.3, splitcex=0.3, titlecex=0.5, 
                                plotsplits=FALSE, cornercoords_loc=scriptdir, include_null_range=TRUE, 
                                tr=tr, tipranges=tipranges)
# Pie chart
plot_BGB(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", 
         label.offset=0.45, tipcex=0.3, statecex=0.5, splitcex=0.001, titlecex=0.5, 
         plotsplits=F, cornercoords_loc=scriptdir, include_null_range=TRUE, 
         tr=tr, tipranges=tipranges, 
         #skiplabels = T,
         #simplify_piecharts = F,
         tipboxes_TF = T,
         #plotlegend = T, legend_cex = 0.2,
         pie_tip_statecex = 0.3)

#######################################################
# Run DEC+X
#######################################################

# Intitialize a default model (DEC model)
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
# Give BioGeoBEARS the location of the phylogeny Newick file
BioGeoBEARS_run_object$trfn = trfn
# Give BioGeoBEARS the location of the geography text file
BioGeoBEARS_run_object$geogfn = geogfn
# Input the maximum range size
BioGeoBEARS_run_object$max_range_size = max_range_size

# Min branch length to be considered a tip (below threshold is an ancestor)
BioGeoBEARS_run_object$min_branchlength = 0.000001 
# Allow null ranges
BioGeoBEARS_run_object$include_null_range = TRUE

# Provide your time slices
BioGeoBEARS_run_object$timesfn = "CrownFrogs_timeperiods.txt"
# Provide your distance matrices
BioGeoBEARS_run_object$distsfn = "CrownFrogs_distances_scaled.txt"

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 8    # use multiple cores in parallel

# Sparse matrix for huge state spaces (600+)
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE, cut_fossils=FALSE)
# The stratified tree is described in this table:
#BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Add x as a free parameter (function of distance)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","est"] = 0.0

# Run this to check inputs. Read the error messages if you get them!
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

# For a slow analysis, run once, then set runslow=FALSE to just 
# load the saved result.
runslow = TRUE
resfn = "CrownFrogs_DEC_x_H1.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  save(res, file=resfn)
  resDECx = res
} else {
  # Loads to "res"
  load(resfn)
  resDECx = res
}

# Plot the time-stratified DEC results
analysis_titletxt ="DEC Time-Stratified Dispersal (Scaled) CrownFrogs H1"
# Setup
results_object = resDECx
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
# States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", 
                                label.offset=0.45, tipcex=0.3, statecex=0.3, splitcex=0.3, titlecex=0.5, 
                                plotsplits=FALSE, cornercoords_loc=scriptdir, include_null_range=TRUE, 
                                tr=tr, tipranges=tipranges)
# Pie chart
plot_BGB(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", 
         label.offset=0.45, tipcex=0.3, statecex=0.5, splitcex=0.001, titlecex=0.5, 
         plotsplits=FALSE, cornercoords_loc=scriptdir, include_null_range=TRUE, 
         tr=tr, tipranges=tipranges, 
         #skiplabels = T,
         simplify_piecharts = T,
         tipboxes_TF = T,
         pie_tip_statecex = 0.3,
         plotlegend = T, legend_cex = 0.2)


#######################################################
# Run DEC+J+X
#######################################################
# for notes on each command, see the DEC model above
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
BioGeoBEARS_run_object$timesfn = "CrownFrogs_timeperiods.txt"
BioGeoBEARS_run_object$distsfn = "CrownFrogs_distances_scaled.txt"
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 8
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE, cut_fossils=FALSE)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
# Set up DEC+J model
# Get the ML parameter values from the 2-parameter nested model
# (this will ensure that the 3-parameter model always does at least as good)
dstart = resDEC$outputs@params_table["d","est"]
estart = resDEC$outputs@params_table["e","est"]
jstart = 0.0001
xstart = resDEC$outputs@params_table["x","est"]
# Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart
# Add j as a free parameter
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart
# Add x as a free parameter
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","init"] = xstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","est"] = xstart
# Check all is good
check_BioGeoBEARS_run(BioGeoBEARS_run_object)
# Run DEC+J
resfn = "CrownFrogs_DEC_j_x_H1.Rdata"
runslow = TRUE
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  
  resDECjx = res
} else {
  # Loads to "res"
  load(resfn)
  resDECjx = res
}

# Plot the time-stratified DEC results
analysis_titletxt ="DEC+J Time-Stratified Dispersal (Scaled) CrownFrogs H1"
# Setup
results_object = resDECjx
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
# States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", 
                                label.offset=0.45, tipcex=0.3, statecex=0.3, splitcex=0.3, titlecex=0.5, 
                                plotsplits=FALSE, cornercoords_loc=scriptdir, include_null_range=TRUE, 
                                tr=tr, tipranges=tipranges)
# Pie chart
plot_BGB(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", 
         label.offset=0.45, tipcex=0.3, statecex=0.5, splitcex=0.001, titlecex=0.5, 
         plotsplits=FALSE, cornercoords_loc=scriptdir, include_null_range=TRUE, 
         tr=tr, tipranges=tipranges, 
         #skiplabels = T,
         #simplify_piecharts = F,
         tipboxes_TF = T,
         pie_tip_statecex = 0.3,
         plotlegend = T, legend_cex = 0.2)


#####################################################
#####################################################
#####################################################

#######################################################
# Run DEC+X+W
#######################################################
# more info: https://groups.google.com/g/biogeobears/c/BmTiLh13PO4
# more info: https://groups.google.com/g/biogeobears/c/FVsXo4Z38Z8/m/Wo-u27rGBQAJ

# Intitialize a default model (DEC model)
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
# Give BioGeoBEARS the location of the phylogeny Newick file
BioGeoBEARS_run_object$trfn = trfn
# Give BioGeoBEARS the location of the geography text file
BioGeoBEARS_run_object$geogfn = geogfn
# Input the maximum range size
BioGeoBEARS_run_object$max_range_size = max_range_size

# Min branch length to be considered a tip (below threshold is an ancestor)
BioGeoBEARS_run_object$min_branchlength = 0.000001 
# Allow null ranges
BioGeoBEARS_run_object$include_null_range = TRUE

# Provide your time slices
BioGeoBEARS_run_object$timesfn = "CrownFrogs_timeperiods.txt"
# Provide your distance matrices
BioGeoBEARS_run_object$distsfn = "CrownFrogs_distances_scaled.txt"
# Provide your dispersal matrices
BioGeoBEARS_run_object$dispersal_multipliers_fn = "CrownFrogs_dispersal.txt"

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 8    # use multiple cores in parallel

# Sparse matrix for huge state spaces (600+)
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE, cut_fossils=FALSE)
# The stratified tree is described in this table:
#BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Add x as a free parameter (function of distance)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","est"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","min"] = -3
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","max"] = 0

# Add w as a free parameter (function of dispersal probability)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["w","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["w","init"] = 0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["w","est"] = 0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["w","min"] = 0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["w","max"] = 3

# Run this to check inputs. Read the error messages if you get them!
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

# For a slow analysis, run once, then set runslow=FALSE to just 
# load the saved result.
runslow = TRUE
resfn = "CrownFrogs_DEC_x_w_H2_Fossils.Rdata"
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  save(res, file=resfn)
  resDECxw = res
} else {
  # Loads to "res"
  load(resfn)
  resDECxw = res
}

# Plot the time-stratified DEC results
analysis_titletxt ="DEC Time-Stratified +x (Scaled) +w CrownFrogs w/Fossils H2"
# Setup
results_object = resDECxw
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
# States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", 
                                label.offset=0.45, tipcex=0.3, statecex=0.3, splitcex=0.3, titlecex=0.5, 
                                plotsplits=FALSE, cornercoords_loc=scriptdir, include_null_range=TRUE, 
                                tr=tr, tipranges=tipranges)
# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", 
                         label.offset=0.45, tipcex=0.3, statecex=0.7, splitcex=0.7, titlecex=0.5, 
                         plotsplits=FALSE, cornercoords_loc=scriptdir, include_null_range=TRUE, 
                         tr=tr, tipranges=tipranges)
plot_BGB(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", 
         label.offset=0.45, tipcex=0.3, statecex=0.5, splitcex=0.001, titlecex=0.5, 
         plotsplits=FALSE, cornercoords_loc=scriptdir, include_null_range=TRUE, 
         tr=tr, tipranges=tipranges, 
         #skiplabels = T,
         #simplify_piecharts = F,
         tipboxes_TF = T,
         pie_tip_statecex = 0.3,
         plotlegend = T, legend_cex = 0.2)


#######################################################
# Run DEC+J+X+W
#######################################################
# for notes on each command, see the DEC model above
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
BioGeoBEARS_run_object$timesfn = "CrownFrogs_timeperiods.txt"
BioGeoBEARS_run_object$distsfn = "CrownFrogs_distances_scaled.txt"
BioGeoBEARS_run_object$dispersal_multipliers_fn = "CrownFrogs_dispersal.txt"
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = 8
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)
BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE, cut_fossils=FALSE)
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run
# Set up DEC+J model
# Get the ML parameter values from the 2-parameter nested model
# (this will ensure that the 3-parameter model always does at least as good)
dstart = resDECxw$outputs@params_table["d","est"]
estart = resDECxw$outputs@params_table["e","est"]
jstart = 0.0001
xstart = resDECxw$outputs@params_table["x","est"]
wstart = resDECxw$outputs@params_table["w","est"]
# Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart
# Add j as a free parameter
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart
# Add x as a free parameter
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","init"] = xstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","est"] = xstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","min"] = -3
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","max"] = 0
# Add w as a free parameter (function of dispersal probability)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["w","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["w","init"] = 0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["w","est"] = 0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["w","min"] = 0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["w","max"] = 3
# Check all is good
check_BioGeoBEARS_run(BioGeoBEARS_run_object)
# Run DEC+J
resfn = "CrownFrogs_DEC_j_x_w_H2_Fossils.Rdata"
runslow = TRUE
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  
  resDECjxw = res
} else {
  # Loads to "res"
  load(resfn)
  resDECjxw = res
}

# Plot the time-stratified DEC results
analysis_titletxt ="DEC+J Time-Stratified +x (Scaled) +w CrownFrogs w/Fossils H1"
# Setup
results_object = resDECjxw
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
# States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", 
                                label.offset=0.45, tipcex=0.3, statecex=0.3, splitcex=0.3, titlecex=0.5, 
                                plotsplits=FALSE, cornercoords_loc=scriptdir, include_null_range=TRUE, 
                                tr=tr, tipranges=tipranges)
# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", 
                         label.offset=0.45, tipcex=0.3, statecex=0.7, splitcex=0.7, titlecex=0.5, 
                         plotsplits=FALSE, cornercoords_loc=scriptdir, include_null_range=TRUE, 
                         tr=tr, tipranges=tipranges)
plot_BGB(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", 
         label.offset=0.45, tipcex=0.3, statecex=0.5, splitcex=0.001, titlecex=0.5, 
         plotsplits=FALSE, cornercoords_loc=scriptdir, include_null_range=TRUE, 
         tr=tr, tipranges=tipranges, 
         #skiplabels = T,
         simplify_piecharts = T,
         tipboxes_TF = T,
         pie_tip_statecex = 0.3,
         plotlegend = T, legend_cex = 0.2)








#######################################################
# Make a table that is the abbreviations key
#######################################################
# Here, the file with the "key" is a tab-delimited text file -- 
# change to your file location if you don't want Hawaii
keyfn = "CrownFrogs_RangesKey.txt"

keydf = read.table(keyfn, header=TRUE, sep="\t", stringsAsFactors=FALSE, fill=TRUE)
keydf

keydf2 = keydf[,c("ab1", "desc")]
names(keydf2) = c("Abbr", "Description")
keydf2

library(gridExtra)	# for grid.table() function
grid.table(keydf2, show.rownames=FALSE)


#######################################################
# Plot legend for SOME states/ranges 
# (especially since the widespread ranges are just
#  combinations of the primary colors used for 
#  single-area ranges)
#######################################################

pdffn = "colors_legend_some_v1.pdf"
pdf(pdffn, width=6, height=6)


# Subset to just some ranges (since there are sooo many combinations)
states_to_put_in_legend = c(1,2,3,4,5,6,7,8)
colors_list_for_states_subset = colors_list_for_states[states_to_put_in_legend]
possible_ranges_list_txt_subset = possible_ranges_list_txt[states_to_put_in_legend]

legend_ncol=NULL
legend_cex=1.5
colors_legend(possible_ranges_list_txt_subset, colors_list_for_states_subset, legend_ncol=legend_ncol, legend_cex=legend_cex)

dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)

