#plot_BGB(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", 
#                         label.offset=0.45, tipcex=0.3, statecex=0.5, splitcex=0.001, titlecex=0.5, 
#                         plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, 
#                         tr=tr, tipranges=tipranges, 
#                         #skiplabels = T,
#                         #simplify_piecharts = F,
#                         tipboxes_TF = T,
#                         pie_tip_statecex = 0.3,
#                         plotlegend = T, legend_cex = 0.2)
#
plot_BGB <- function (results_object, analysis_titletxt = NULL, addl_params = list(), 
                      plotwhat = "text", label.offset = NULL, tipcex = 0.8, statecex = 0.7, 
                      splitcex = 0.6, titlecex = 0.8, plotsplits = TRUE, plotlegend = FALSE, 
                      legend_ncol = NULL, legend_cex = 1, cornercoords_loc = "auto", 
                      tr = NULL, tipranges = NULL, if_ties = "takefirst", pie_tip_statecex = 0.7, 
                      juststats = FALSE, xlab = "Millions of years ago", root.edge = TRUE, 
                      colors_list_for_states = NULL, skiptree = FALSE, show.tip.label = TRUE, 
                      tipcol = "black", dej_params_row = NULL, plot_max_age = NULL, 
                      skiplabels = FALSE, plot_stratum_lines = TRUE, include_null_range = NULL, 
                      plot_null_range = FALSE, simplify_piecharts = FALSE, tipboxes_TF = TRUE, 
                      tiplabel_adj = c(0.5), no.margin = FALSE, xlims = NULL, ylims = NULL) 
{
  junk = "\n\t# manual_ranges_txt=NULL, \n\t# @manual_ranges_txt If you dont want to use the default text for each range, produced\n\t# by areas_list_to_states_list_new(), specify the list here.\n\n\t\n\tscriptdir = \"/Dropbox/_njm/__packages/BioGeoBEARS_setup/inst/extdata/a_scripts/\"\n\tplot_BioGeoBEARS_results(results_object, analysis_titletxt=NULL, addl_params=list(), plotwhat=\"text\", label.offset=NULL, tipcex=0.8, statecex=0.8, splitcex=0.8, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=NULL, tipranges=NULL)\n\t\n\t# Defaults\n\taddl_params=list(\"j\"); plotwhat=\"text\"; label.offset=0.45; tipcex=0.7; statecex=0.7; splitcex=0.6; titlecex=0.8; plotsplits=TRUE; cornercoords_loc=scriptdir; include_null_range=TRUE; tr=tr; tipranges=tipranges; juststats = FALSE; plotlegend=FALSE; \txlab=\"Millions of years ago\"; if_ties=\"takefirst\"\n\t\n\t\n\t# Setup\nresults_object = resDEC\nanalysis_titletxt =\"BioGeoBEARS DEC on Mariana M1v4_unconstrained\"\naddl_params=list(\"j\"); plotwhat=\"text\"; label.offset=0.45; tipcex=0.7; statecex=0.7; splitcex=0.6; titlecex=0.8; plotsplits=TRUE; cornercoords_loc=scriptdir; include_null_range=TRUE; tr=tr; tipranges=tipranges\njuststats=FALSE; plotlegend=FALSE; \txlab=\"Millions of years ago\"; if_ties=\"takefirst\"\nshow.tip.label=TRUE\ntipcol=\"black\"; dej_params_row=NULL; plot_max_age=NULL; skiplabels=FALSE; \ncolors_list_for_states=NULL\nskiptree=FALSE\ninclude_null_range=NULL\nplot_stratum_lines=TRUE\n\tplot_null_range = FALSE\n\t"
  if (is.null(include_null_range) == TRUE) {
    include_null_range = results_object$inputs$include_null_range
  }
  results_object$inputs$include_null_range = include_null_range
  tmp_fg = par("fg")
  par(fg = "black")
  BioGeoBEARS_run_object = results_object$inputs
  if (is.null(tr)) {
    tr = check_trfn(trfn = BioGeoBEARS_run_object$trfn)
  }
  tr_pruningwise = reorder(tr, "pruningwise")
  tips = 1:length(tr_pruningwise$tip.label)
  nodes = (length(tr_pruningwise$tip.label) + 1):(length(tr_pruningwise$tip.label) + 
                                                    tr_pruningwise$Nnode)
  if (is.null(tipranges)) {
    if (BioGeoBEARS_run_object$use_detection_model == FALSE) {
      tipranges = getranges_from_LagrangePHYLIP(lgdata_fn = np(BioGeoBEARS_run_object$geogfn))
    }
    if (BioGeoBEARS_run_object$use_detection_model == TRUE) {
      if (BioGeoBEARS_run_object$use_detection_model == 
          TRUE) {
        tipranges = tipranges_from_detects_fn(detects_fn = BioGeoBEARS_run_object$detects_fn)
      }
    }
  }
  areas = getareas_from_tipranges_object(tipranges)
  areas
  numareas = length(areas)
  numareas
  if (!is.na(results_object$inputs$max_range_size)) {
    max_range_size = results_object$inputs$max_range_size
  }
  else {
    max_range_size = length(areas)
  }
  max_range_size
  if (is.null(results_object$inputs$states_list)) {
    numstates = numstates_from_numareas(numareas = length(areas), 
                                        maxareas = max_range_size, include_null_range = results_object$inputs$include_null_range)
    numstates
    states_list_areaLetters = areas_list_to_states_list_new(areas, 
                                                            maxareas = max_range_size, include_null_range = results_object$inputs$include_null_range)
    states_list_0based_index = rcpp_areas_list_to_states_list(areas, 
                                                              maxareas = max_range_size, include_null_range = results_object$inputs$include_null_range)
  }
  else {
    states_list_0based_index = results_object$inputs$states_list
  }
  param_ests = extract_params_from_BioGeoBEARS_results_object(results_object, 
                                                              returnwhat = "table", addl_params = addl_params, paramsstr_digits = 4)
  if (juststats == TRUE) {
    return(param_ests)
  }
  else {
    paramstr = extract_params_from_BioGeoBEARS_results_object(results_object, 
                                                              returnwhat = "string", addl_params = addl_params, 
                                                              paramsstr_digits = 4)
  }
  param_names = extract_params_from_BioGeoBEARS_results_object(results_object, 
                                                               returnwhat = "param_names", addl_params = addl_params, 
                                                               paramsstr_digits = 4)
  if (is.null(analysis_titletxt)) {
    tmptxt = results_object$inputs$description
    if (any(is.null(tmptxt), tmptxt == "", tmptxt == "defaults", 
            tmptxt == "default")) {
      analysis_titletxt = ""
    }
    else {
      analysis_titletxt = results_object$inputs$description
    }
  }
  if (is.null(dej_params_row)) {
    analysis_titletxt = paste(analysis_titletxt, "\n", "ancstates: global optim, ", 
                              max_range_size, " areas max. ", paramstr, sep = "")
    analysis_titletxt
  }
  else {
    dej_params_row
    brate_col_TF = names(dej_params_row) == "brate"
    brate_col = (1:length(dej_params_row))[brate_col_TF]
    biogeog_params = dej_params_row[1:(brate_col - 1)]
    biogeog_param_names = names(dej_params_row)[1:(brate_col - 
                                                     1)]
    equals_col = "="
    tmpcols = cbind(biogeog_param_names, equals_col, unlist(biogeog_params))
    tmpcols
    txtrows = apply(X = tmpcols, MARGIN = 1, FUN = paste, 
                    sep = "", collapse = "")
    txtrows
    biogeog_params_txt = paste(txtrows, sep = "", collapse = "; ")
    biogeog_params_txt
    titletxt2 = bquote(paste(.(max_range_size), " areas max., ", 
                             .(biogeog_params_txt), "; ", lambda, "=", .(dej_params_row$brate), 
                             "; ", mu, "=", .(dej_params_row$drate), "; ", alpha, 
                             "=", .(dej_params_row$rangesize_b_exponent), "; ", 
                             omega, "=", .(dej_params_row$rangesize_d_exponent), 
                             "", sep = ""))
  }
  leftright_nodes_matrix = get_leftright_nodes_matrix_from_results(tr_pruningwise)
  marprobs = results_object$ML_marginal_prob_each_state_at_branch_bottom_below_node
  left_ML_marginals_by_node = marprobs[leftright_nodes_matrix[, 
                                                              2], ]
  right_ML_marginals_by_node = marprobs[leftright_nodes_matrix[, 
                                                               1], ]
  right_ML_marginals_by_node
  if (is.null(dim(left_ML_marginals_by_node))) {
    left_ML_marginals_by_node = matrix(data = left_ML_marginals_by_node, 
                                       nrow = 1)
  }
  if (is.null(dim(right_ML_marginals_by_node))) {
    right_ML_marginals_by_node = matrix(data = right_ML_marginals_by_node, 
                                        nrow = 1)
  }
  relprobs_matrix = results_object$ML_marginal_prob_each_state_at_branch_top_AT_node
  if (length(nodes) > 1) {
    relprobs_matrix_for_internal_states = relprobs_matrix[nodes, 
    ]
  }
  else {
    relprobs_matrix_for_internal_states = relprobs_matrix[nodes, 
    ]
    relprobs_matrix_for_internal_states = matrix(data = relprobs_matrix_for_internal_states, 
                                                 nrow = 1, ncol = ncol(relprobs_matrix))
  }
  relprobs_matrix
  if (is.null(states_list_0based_index)) {
    statenames = areas_list_to_states_list_new(areas, maxareas = max_range_size, 
                                               include_null_range = results_object$inputs$include_null_range, 
                                               split_ABC = FALSE)
    ranges_list = as.list(statenames)
    statenames
  }
  else {
    ranges_list = states_list_0based_to_ranges_txt_list(state_indices_0based = states_list_0based_index, 
                                                        areanames = areas)
    ranges_list
    statenames = unlist(ranges_list)
    statenames
  }
  MLprobs = get_ML_probs(relprobs_matrix)
  MLprobs
  MLstates = get_ML_states_from_relprobs(relprobs_matrix, statenames, 
                                         returnwhat = "states", if_ties = if_ties)
  if (is.null(colors_list_for_states)) {
    colors_matrix = get_colors_for_numareas(length(areas))
    colors_list_for_states = mix_colors_for_states(colors_matrix, 
                                                   states_list_0based_index, plot_null_range = results_object$inputs$include_null_range)
    colors_list_for_states
  }
  if (is.null(ranges_list)) {
    possible_ranges_list_txt = areas_list_to_states_list_new(areas, 
                                                             maxareas = max_range_size, split_ABC = FALSE, include_null_range = results_object$inputs$include_null_range)
  }
  else {
    possible_ranges_list_txt = ranges_list
  }
  cols_byNode = rangestxt_to_colors(possible_ranges_list_txt, 
                                    colors_list_for_states, MLstates)
  if (plotlegend == TRUE) {
    colors_legend(possible_ranges_list_txt, colors_list_for_states, 
                  legend_ncol = legend_ncol, legend_cex = legend_cex)
  }
  if (root.edge == FALSE) {
    tr$root.edge = 0
  }
  if (root.edge == TRUE) {
    if (is.null(tr$root.edge) == TRUE) {
      tr$root.edge = 0
    }
  }
  if (is.null(label.offset)) {
    label.offset = 0.007 * (get_max_height_tree(tr) + tr$root.edge)
  }
  if (show.tip.label == TRUE) {
    if (is.null(plot_max_age)) {
      max_x = 1.25 * (get_max_height_tree(tr) + tr$root.edge)
      min_x = 0
    }
    else {
      nontree_part_of_x = plot_max_age - (get_max_height_tree(tr) + 
                                            tr$root.edge)
      max_x = 1.25 * (get_max_height_tree(tr) + tr$root.edge)
      min_x = -1 * nontree_part_of_x
    }
  }
  else {
    if (is.null(plot_max_age)) {
      max_x = 1.05 * (get_max_height_tree(tr) + tr$root.edge)
      min_x = 0
    }
    else {
      nontree_part_of_x = plot_max_age - (get_max_height_tree(tr) + 
                                            tr$root.edge)
      max_x = 1.05 * (get_max_height_tree(tr) + tr$root.edge)
      min_x = -1 * nontree_part_of_x
    }
  }
  max_tree_x = 1 * (get_max_height_tree(tr) + tr$root.edge)
  if (is.null(xlims)) {
    xlims = c(min_x, max_x)
  }
  else {
    xlims = xlims
  }
  nodecoords = node_coords(tr, tmplocation = cornercoords_loc, 
                           root.edge = root.edge)
  max_tree_x = max(nodecoords$x)
  if (is.null(plot_max_age)) {
    xticks_desired_lims = c(0, max_tree_x)
  }
  else {
    xticks_desired_lims = c(0, plot_max_age)
  }
  xticks_desired = pretty(xticks_desired_lims)
  xaxis_ticks_locs = max_tree_x - xticks_desired
  if (skiptree != TRUE) {
    plot(tr_pruningwise, x.lim = xlims, y.lim = ylims, show.tip.label = FALSE, 
         label.offset = label.offset, cex = tipcex, no.margin = no.margin, 
         root.edge = root.edge, edge.width = 0.5)
    if (show.tip.label == TRUE) {
      tiplabels_to_plot = sapply(X = tr_pruningwise$tip.label, 
                                 FUN = substr, start = 1, stop = 30)
      if (skiplabels == FALSE) {
        tiplabels(text = tiplabels_to_plot, tip = tips, 
                  cex = tipcex, adj = 0, bg = "white", frame = "n", 
                  pos = 4, offset = label.offset, col = tipcol)
      }
    }
    axis(side = 1, at = xaxis_ticks_locs, labels = xticks_desired)
    mtext(text = xlab, side = 1, line = 2, cex = titlecex)
  }
  if (plotwhat == "text") {
    par(fg = tmp_fg)
    if (skiplabels == FALSE) {
      nodelabels(text = MLstates[nodes], node = nodes, 
                 bg = cols_byNode[nodes], cex = statecex, fr="c")
      tiplabels(text = MLstates[tips], tip = tips, bg = cols_byNode[tips], 
                cex = statecex, adj = tiplabel_adj)
    }
    par(fg = "black")
  }
  if (plotwhat == "pie") {
    par(fg = tmp_fg)
    if (skiplabels == FALSE) {
      if (simplify_piecharts == TRUE) {
        colnums_to_keep_in_probs = NULL
        probs = results_object$ML_marginal_prob_each_state_at_branch_top_AT_node
        probs2 = probs
        maxprob = rep(0, nrow(probs))
        other = rep(0, nrow(probs))
        num_to_keep = 1
        cat("\nSince simplify_piecharts==TRUE, reducing prob pie charts to (most probable, other)...\n")
        for (i in 1:nrow(probs)) {
          cat(i, " ", sep = "")
          tmprow = probs[i, ]
          positions_highest_prob_to_lowest = rev(order(tmprow))
          positions_to_keep = positions_highest_prob_to_lowest[1:num_to_keep]
          colnums_to_keep_in_probs = c(colnums_to_keep_in_probs, 
                                       positions_to_keep)
          keepTF = rep(FALSE, length(tmprow))
          keepTF[positions_to_keep] = TRUE
          otherTF = keepTF == FALSE
          other[i] = sum(tmprow[otherTF])
          tmprow[otherTF] = 0
          probs2[i, ] = tmprow
        }
        cat("\n")
        colnums_to_keep_in_probs_in_order = sort(unique(colnums_to_keep_in_probs))
        probs3 = cbind(probs2[, colnums_to_keep_in_probs_in_order], 
                       other)
        probs3 = probs3[nodes, ]
        newcols = c(colors_list_for_states[colnums_to_keep_in_probs_in_order], 
                    "white")
        nodelabels(pie = probs3, node = nodes, piecol = newcols, 
                   cex = statecex)
      }
      else {
        nodelabels(pie = relprobs_matrix_for_internal_states, 
                   node = nodes, piecol = colors_list_for_states, 
                   cex = statecex)
      }
      if (tipboxes_TF == TRUE) {
        tiplabels(text = MLstates[tips], tip = tips, 
                  bg = cols_byNode[tips], cex = pie_tip_statecex, 
                  adj = tiplabel_adj, fr="c")
      }
    }
    par(fg = "black")
  }
  if (skiptree != TRUE) {
    if (titlecex > 0) {
      par(cex.main = titlecex)
      title(analysis_titletxt)
      if (!is.null(dej_params_row)) {
        title(titletxt2, line = 1)
      }
    }
  }
  if (plotsplits == TRUE) {
    if (cornercoords_loc == "manual") {
      stoptxt = cat("\nNOTE: To plot splits, this function needs to access the function 'plot_phylo3_nodecoords'.\n", 
                    "The function is modified from an APE function, and cannot be directly included in the package,\n", 
                    "due to some C code that does not meet CRAN standards. To solve this, give plot_BioGeoBEARS_results\n", 
                    "a 'cornercoords_loc' string that gives the directory of plot_phylo3_nodecoords.R.  Typically this\n", 
                    "can be found via: ", "tmp=np(system.file(\"extdata/a_scripts\", package=\"BioGeoBEARS\"))\n", 
                    "then: list.files(tmp); print(tmp)\n", sep = "")
      plotsplits = FALSE
    }
  }
  if (plotsplits == TRUE) {
    coords = corner_coords(tr, tmplocation = cornercoords_loc, 
                           root.edge = root.edge)
    coords
    relprobs_matrix = left_ML_marginals_by_node
    if (plotwhat == "text") {
      MLprobs = get_ML_probs(relprobs_matrix)
      MLprobs
      MLstates = get_ML_states_from_relprobs(relprobs_matrix, 
                                             statenames, returnwhat = "states", if_ties = if_ties)
      MLstates
      length(MLstates)
      cols_byNode = rangestxt_to_colors(possible_ranges_list_txt, 
                                        colors_list_for_states, MLstates)
      par(fg = tmp_fg)
      if (skiplabels == FALSE) {
        cornerlabels(text = MLstates, coords = coords$leftcorns, 
                     bg = cols_byNode, cex = splitcex)
      }
      par(fg = "black")
    }
    if (plotwhat == "pie") {
      par(fg = tmp_fg)
      cornerpies(pievals = relprobs_matrix, coords$leftcorns, 
                 piecol = colors_list_for_states, cex = splitcex)
      par(fg = "black")
    }
    relprobs_matrix = right_ML_marginals_by_node
    if (plotwhat == "text") {
      MLprobs = get_ML_probs(relprobs_matrix)
      MLprobs
      MLstates = get_ML_states_from_relprobs(relprobs_matrix, 
                                             statenames, returnwhat = "states", if_ties = if_ties)
      MLstates
      length(MLstates)
      cols_byNode = rangestxt_to_colors(possible_ranges_list_txt, 
                                        colors_list_for_states, MLstates)
      par(fg = tmp_fg)
      if (skiplabels == FALSE) {
        cornerlabels(text = MLstates, coords = coords$rightcorns, 
                     bg = cols_byNode, cex = splitcex)
      }
      par(fg = "black")
    }
    if (plotwhat == "pie") {
      par(fg = tmp_fg)
      cornerpies(pievals = relprobs_matrix, coords$rightcorns, 
                 piecol = colors_list_for_states, cex = splitcex)
      par(fg = "black")
    }
  }
  if (((is.null(BioGeoBEARS_run_object$timeperiods) == FALSE)) && 
      (plot_stratum_lines == TRUE)) {
    timeperiods = BioGeoBEARS_run_object$timeperiods
    line_positions_on_plot = add_statum_boundaries_to_phylo_plot(tr, 
                                                                 timeperiods = timeperiods, lty = "dashed", col = "gray50", 
                                                                 plotlines = TRUE)
  }
  param_ests = matrix(data = param_ests, nrow = 1)
  param_ests = adf2(param_ests)
  param_ests = dfnums_to_numeric(param_ests)
  names(param_ests) = c("LnL", "nparams", param_names)
  return(param_ests)
}
