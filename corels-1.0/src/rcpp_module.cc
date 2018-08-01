#include <Rcpp.h>
#include <assert.h>

#include "params.h"
#include "run.hh"
#include "utils.hh"
#include "alloc.hh"

#define BUFSZ 512

void load_params(Rcpp::List *params_list, run_params_t *params) {
    Rcpp::NumericVector nv;
    Rcpp::IntegerVector iv;
    Rcpp::LogicalVector lv;
    Rcpp::List verbosity;

    iv = params_list[0];
    params->curiosity_policy = iv[0];
    iv = params_list[1];
    params->max_num_nodes = iv[0];
    nv = params_list[2];
    params->c = nv[0];
    verbosity = params_list[3];
    for (int i = 0; i < verbosity.size(); i++) {
        strcat(params->vstring, verbosity[i]);
    }
    iv = params_list[4];
    params->map_type = iv[0];
    iv = params_list[5];
    params->ablation = iv[0];
    lv = params_list[6];
    params->calculate_size = lv[0];
    lv = params_list[7];
    params->latex_out = lv[0];
}

int load_data(const char *data_file, const char *label_file, const char *minor_file,
	               int *ret_nsamples, int *ret_nrules, int *ret_nlabels, int *ret_nmeta, rule_t **rules, rule_t **labels, rule_t **meta) {
    int nrules, nlabels, nsamples_lchk, nsamples_mchk, nmeta, ret;

    /* Load data. */
    if ((ret = rules_init(data_file, ret_nrules, ret_nsamples, rules, 1)) != 0)
        return (ret);

    /* Load labels. */
    if ((ret = rules_init(label_file, &nlabels, &nsamples_lchk, labels, 0)) != 0) {
        rules_free(*rules, nrules, 1);
        return (ret);
    }

    if (minor_file && (ret = rules_init(minor_file, ret_nmeta, &nsamples_mchk, meta, 0)) != 0) {
        rules_free(*rules, nrules, 1);
        rules_free(*labels, nlabels, 0);
        return (ret);
    }

    assert(nlabels == 2);
    assert(nsamples_lchk == *ret_nsamples);
    return (0);
}

// [[Rcpp::export]]
Rcpp::List _train (Rcpp::List param_list, Rcpp::CharacterVector data_fname, Rcpp::CharacterVector label_fname, Rcpp::CharacterVector minor_fname) {
    const char *df = Rcpp::as<std::string>(data_fname).c_str();
    const char *lf = Rcpp::as<std::string>(label_fname).c_str();
    const char *mf = (minor_fname) ? Rcpp::as<std::string>(minor_fname).c_str() : NULL;

    run_params_t params;
    set_default_params(&params);
    load_params(&param_list, &params);

    char error_txt[BUFSZ];

    if(!df || strlen(df) == 0) {
        snprintf(error_txt, BUFSZ, "training data file must be a valid file path");
        goto error;
    }
    if(!lf || strlen(lf) == 0) {
        snprintf(error_txt, BUFSZ, "label file must be a valid file path");
        goto error;
    }
    if(!params.opt_fname || !strlen(params.opt_fname)) {
        snprintf(error_txt, BUFSZ, "optimal rulelist file must be a valid file path");
        goto error;
    }
    if(!params.log_fname || !strlen(params.log_fname)) {
        snprintf(error_txt, BUFSZ, "log file must be a valid file path");
        goto error;
    }
    if (params.max_num_nodes < 0) {
        snprintf(error_txt, BUFSZ, "maximum number of nodes must be positive");
        goto error;
    }
    if (params.c < 0.0) {
        snprintf(error_txt, BUFSZ, "regularization constant must be postitive");
        goto error;
    }
    if (params.map_type > 2 || params.map_type < 0) {
        snprintf(error_txt, BUFSZ, "symmetry-aware map must be (0|1|2)");
        goto error;
    }
    if (params.curiosity_policy < 0 || params.curiosity_policy > 4) {
        snprintf(error_txt, BUFSZ, "you must specify a curiosity type (0|1|2|3|4)");
        goto error;
    }

    if (load_data(df, lf, mf, &params.nsamples, &params.nrules, &params.nlabels, &params.nmeta,
                    &params.rules, &params.labels, &params.meta) != 0) {
        snprintf(error_txt, BUFSZ, "error loading data");
        goto error;
    }

error:
    ::Rf_error(error_txt);
    return NULL;

    tracking_vector<unsigned short, DataStruct::Tree> rulelist;
    tracking_vector<bool, DataStruct::Tree> preds;

    run_corels(params, rulelist, preds);

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
    
    Rcpp::DataFrame rs = Rcpp::DataFrame::create(Rcpp::Named("V1")=rule_v, Rcpp::Named("V2")=pred_v);

    return(Rcpp::List::create(Rcpp::Named("rs")=rs));
}

RcppExport SEXP corels_train(SEXP paramSEXP, SEXP dataSEXP, SEXP labelSEXP, SEXP minorSEXP) {
    BEGIN_RCPP
    Rcpp::traits::input_parameter<Rcpp::List>::type params(paramSEXP);
    Rcpp::traits::input_parameter<Rcpp::CharacterVector>::type data(dataSEXP);
    Rcpp::traits::input_parameter<Rcpp::CharacterVector>::type label(labelSEXP);
    Rcpp::traits::input_parameter<Rcpp::CharacterVector>::type minor(minorSEXP);

    return Rcpp::wrap(_train(params, data, label, minor));
    END_RCPP
}
