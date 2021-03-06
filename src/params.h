#pragma once

#include "rule.h"

#define D_OPT_FNAME          "corels-opt.txt"
#define D_LOG_FNAME          "corels.txt"
#define D_MAX_NUM_NODES      10000
#define D_C                  0.01
#define D_VSTRING            "progress"
#define D_CURIOSITY_POLICY   2
#define D_MAP_TYPE           1
#define D_FREQ               1000
#define D_ABLATION           0
#define D_CALCULATE_SIZE     0
#define D_LATEX_OUT          0
#define D_V_PROGRESS         1
#define D_V_LOG              0
#define D_V_RULE             0
#define D_V_LABEL            0
#define D_V_SAMPLES          0
#define D_V_SILENT           0

#ifdef __cplusplus
extern "C" {
#endif

typedef struct run_params {
    char* opt_fname;
    char* log_fname;
    int max_num_nodes;
    double c;
    int curiosity_policy;
    int map_type;
    int freq;
    int ablation;
    int calculate_size;
    int latex_out;

    int v_progress;
    int v_log;
    int v_rule;
    int v_label;
    int v_samples;
    int v_silent;

    int nrules;
    int nlabels;
    int nsamples;
    int nmeta;

    rule_t* rules;
    rule_t* labels;
    rule_t* meta;
} run_params_t;

void set_default_params(run_params_t* out);

#ifdef __cplusplus
}
#endif
