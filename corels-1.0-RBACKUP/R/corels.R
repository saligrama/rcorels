corels <- function(data_fname,label_fname,meta_fname,curiosity_policy=2,max_nodes=10000, regularization=0.01, verbosity=list("progress"), map_type=2, ablation=0, calculate_size=FALSE, latex_out=FALSE) {
    params_list <- list(curiosity_policy, max_nodes, regularization, verbosity, map_type, ablation, calculate_size, latex_out)
    .Call('_run_corels', PACKAGE='corels', 0, 0, params_list, data_fname, label_fname, meta_fname)
}
