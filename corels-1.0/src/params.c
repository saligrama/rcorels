#include "params.h"

void set_default_params(run_params_t* out)
{
    out->opt_fname = (char*)D_OPT_FNAME;
    out->log_fname = (char*)D_LOG_FNAME;
    out->max_num_nodes = D_MAX_NUM_NODES;
    out->c = D_C;
    out->curiosity_policy = D_CURIOSITY_POLICY;
    out->map_type = D_MAP_TYPE;
    out->freq = D_FREQ;
    out->ablation = D_ABLATION;
    out->calculate_size = D_CALCULATE_SIZE;
    out->latex_out = D_LATEX_OUT;

    out->v_progress = D_V_PROGRESS;
    out->v_log = D_V_LOG;
    out->v_rule = D_V_RULE;
    out->v_label = D_V_LABEL;
    out->v_samples = D_V_SAMPLES;
    out->v_silent = D_V_SILENT;

    out->nrules = 0;
    out->nlabels = 0;
    out->nsamples = 0;
    out->nmeta = 0;
    out->rules = NULL;
    out->labels = NULL;
    out->meta = NULL;
}
