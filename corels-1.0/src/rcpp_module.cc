#include <assert.h>
#include <string.h>

#include <Rcpp.h>

#include "params.h"
#include "run.hh"
#include "utils.hh"
#include "alloc.hh"

#define BUFSZ 512

int load_params(Rcpp::StringVector *params_list, run_params_t *params) {
    // since all the inputs are strings (otherwise Rcpp won't recognize anything), convert them back to ints later
    params->curiosity_policy = std::stoi(Rcpp::as<std::string> (params_list[0][0]));
    params->max_num_nodes = std::stoi(Rcpp::as<std::string> (params_list[0][1]));
    params->c = std::stod(Rcpp::as<std::string> (params_list[0][2]));
    const char* vstring = (Rcpp::as<std::string> (params_list[0][3])).c_str();
    char *vopt = NULL;
    char *vcopy = strdup(vstring);
    while ((vopt = strsep(&vcopy, ",")) != NULL) {
        if (strcmp(vopt, "progress") == 0)
            params->v_progress = 1;
        if (strcmp(vopt, "log") == 0)
            params->v_log = 1;
        if (strcmp(vopt, "rule") == 0)
            params->v_rule = 1;
        if (strcmp(vopt, "label") == 0)
            params->v_label = 1;
        if (strcmp(vopt, "samples") == 0)
            params->v_samples = 1;
        if (strcmp(vopt, "silent") == 0)
            params->v_silent = 1;
    }
    free(vcopy);
    params->map_type = std::stoi(Rcpp::as<std::string> (params_list[0][4]));
    params->ablation = std::stoi(Rcpp::as<std::string> (params_list[0][5]));
    params->calculate_size = std::stoi(Rcpp::as<std::string> (params_list[0][6]));
    params->latex_out = std::stoi(Rcpp::as<std::string> (params_list[0][7]));
    return 0;
}

int load_data(const char *data_file, const char *label_file, const char *minor_file,
	               int *ret_nsamples, int *ret_nrules, int *ret_nlabels, int *ret_nmeta, rule_t **rules, rule_t **labels, rule_t **meta) {
    int nrules, nlabels, nsamples_lchk, nsamples_mchk, nmeta, ret;

    /* Load data. */
    if ((ret = rules_init(data_file, ret_nrules, ret_nsamples, rules, 1)) != 0) {
        printf("%s\n", data_file);
        return (ret);
    }

    /* Load labels. */
    if ((ret = rules_init(label_file, &nlabels, &nsamples_lchk, labels, 0)) != 0) {
        rules_free(*rules, nrules, 1);
        return (ret);
    }

    if (minor_file && strlen(minor_file) > 0 && (ret = rules_init(minor_file, ret_nmeta, &nsamples_mchk, meta, 0)) != 0) {
        rules_free(*rules, nrules, 1);
        rules_free(*labels, nlabels, 0);
        return (ret);
    }

    assert(nlabels == 2);
    assert(nsamples_lchk == *ret_nsamples);
    return (0);
}

// workaround so we don't use goto
Rcpp::List error (char* error_txt) {
    Rcpp::Rcerr << error_txt << std::endl;
    return Rcpp::List::create();
}

// [[Rcpp::export]]
Rcpp::List _train (Rcpp::StringVector param_list, std::string data_fname, std::string label_fname, std::string minor_fname) {
    const char *df = data_fname.c_str();
    const char *lf = label_fname.c_str();
    const char *mf = minor_fname.c_str();

    run_params_t params;
    set_default_params(&params);
    if (load_params(&param_list, &params) != 0)
        return Rcpp::List::create();

    char error_txt[BUFSZ];

    if(!df || strlen(df) == 0) {
        snprintf(error_txt, BUFSZ, "training data file must be a valid file path");
        return error(error_txt);
    }
    if(!lf || strlen(lf) == 0) {
        snprintf(error_txt, BUFSZ, "label file must be a valid file path");
        return error(error_txt);
    }
    if(!params.opt_fname || !strlen(params.opt_fname)) {
        snprintf(error_txt, BUFSZ, "optimal rulelist file must be a valid file path");
        return error(error_txt);
    }
    if(!params.log_fname || !strlen(params.log_fname)) {
        snprintf(error_txt, BUFSZ, "log file must be a valid file path");
        return error(error_txt);
    }
    if (params.max_num_nodes < 0) {
        snprintf(error_txt, BUFSZ, "maximum number of nodes must be positive");
        return error(error_txt);
    }
    if (params.c < 0.0) {
        snprintf(error_txt, BUFSZ, "regularization constant must be postitive");
        return error(error_txt);
    }
    if (params.map_type > 2 || params.map_type < 0) {
        snprintf(error_txt, BUFSZ, "symmetry-aware map must be (0|1|2)");
        return error(error_txt);
    }
    if (params.curiosity_policy < 0 || params.curiosity_policy > 4) {
        snprintf(error_txt, BUFSZ, "you must specify a curiosity type (0|1|2|3|4)");
        return error(error_txt);
    }

    if (load_data(df, lf, mf, &params.nsamples, &params.nrules, &params.nlabels, &params.nmeta,
                    &params.rules, &params.labels, &params.meta) != 0) {
        snprintf(error_txt, BUFSZ, "error loading data");
        return error(error_txt);
    }

    tracking_vector<unsigned short, DataStruct::Tree> rulelist;
    tracking_vector<bool, DataStruct::Tree> preds;
    run_corels(&params, rulelist, preds);

    Rcpp::IntegerVector rule_v;
    for (int i = 0; i < rulelist.size(); i++)
        rule_v.push_back(rulelist[i]);

    Rcpp::LogicalVector pred_v;
    for (int i = 0; i < preds.size(); i++)
        pred_v.push_back(preds[i]);

    rules_free(params.rules, params.nrules, 1);
    rules_free(params.labels, params.nlabels, 0);

    if (params.meta)
        rules_free(params.meta, params.nmeta, 0);

    return(Rcpp::List::create(Rcpp::Named("rule")=rule_v, Rcpp::Named("pred")=pred_v));
}

// [[Rcpp::export]]
RcppExport SEXP corels_train(SEXP paramSEXP, SEXP dataSEXP, SEXP labelSEXP, SEXP minorSEXP) {
    BEGIN_RCPP
    Rcpp::StringVector params = Rcpp::as<Rcpp::StringVector>(paramSEXP);
    std::string data = Rcpp::as<std::string>(dataSEXP);
    std::string label = Rcpp::as<std::string>(labelSEXP);
    std::string minor = Rcpp::as<std::string>(minorSEXP);

    return Rcpp::wrap(_train(params, data, label, minor));
    END_RCPP
}
