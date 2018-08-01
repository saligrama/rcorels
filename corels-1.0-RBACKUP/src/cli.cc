#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include <map>
#include <set>

#include "run.hh"
#include "params.h"

#define BUFSZ 512

int main(int argc, char *argv[]) {
    const char usage[] = "USAGE: %s [-b] "
        "[-n max_num_nodes] [-r regularization] [-v (rule|label|samples|progress|log|silent)] "
        "-c (1|2|3|4) -p (0|1|2) [-f logging_frequency] "
        "-a (0|1|2) [-s] [-L latex_out] "
        "data.out data.label [data.minor]\n\n"
        "%s\n";

    extern char *optarg;

    run_params_t params;
    set_default_params(&params);

    char ch;
    bool error = false;
    char error_txt[BUFSZ];

    bool run_bfs = false;
    bool run_curiosity = false;
    bool vstring_malloc = false;

    /* only parsing happens here */
    while ((ch = getopt(argc, argv, "bsLc:p:v:n:r:f:a:")) != -1) {
        switch (ch) {
        case 'b':
            run_bfs = true;
            params.curiosity_policy = 0;
            break;
        case 's':
            params.calculate_size = true;
            break;
        case 'c':
            run_curiosity = true;
            params.curiosity_policy = atoi(optarg);
            break;
        case 'L':
            params.latex_out = true;
            break;
        case 'p':
            params.map_type = atoi(optarg);
            break;
        case 'v':
            vstring_malloc = true;
            params.vstring = (char*)malloc(sizeof(char) * strlen(optarg));
            strcpy(params.vstring, optarg);
            break;
        case 'n':
            params.max_num_nodes = atoi(optarg);
            break;
        case 'r':
            params.c = atof(optarg);
            break;
        case 'f':
            params.freq = atoi(optarg);
            break;
        case 'a':
            params.ablation = atoi(optarg);
            break;
        default:
            error = true;
            snprintf(error_txt, BUFSZ, "unknown option: %c", ch);
        }
    }
    if (params.max_num_nodes < 0) {
        error = true;
        snprintf(error_txt, BUFSZ, "number of nodes must be positive");
    }
    if (params.c < 0) {
        error = true;
        snprintf(error_txt, BUFSZ, "regularization constant must be postitive");
    }
    if (params.map_type > 2 || params.map_type < 0) {
        error = true;
        snprintf(error_txt, BUFSZ, "symmetry-aware map must be (0|1|2)");
    }
    if ((run_bfs + run_curiosity) != 1) {
        error = true;
        snprintf(error_txt, BUFSZ,
                "you must use exactly one of (-b | -c)");
    }
    if (argc < 2 + optind) {
        error = true;
        snprintf(error_txt, BUFSZ,
                "you must specify data files for rules and labels");
    }
    if (run_curiosity && !((params.curiosity_policy >= 1) && (params.curiosity_policy <= 4))) {
        error = true;
        snprintf(error_txt, BUFSZ,
                "you must specify a curiosity type (1|2|3|4)");
    }

    if (error) {
        fprintf(stderr, usage, argv[0], error_txt);
        return 1;
    }

    argc -= optind;
    argv += optind;

    int nsamples_chk;
    if(rules_init(argv[0], &params.nrules, &params.nsamples, &params.rules, 1) != 0) {
        fprintf(stderr, "could not load out file at path '%s'\n", argv[0]);
        return 1;
    }
    if(rules_init(argv[1], &params.nlabels, &nsamples_chk, &params.labels, 0) != 0) {
        fprintf(stderr, "could not load label file at path '%s'\n", argv[1]);
        rules_free(params.rules, params.nrules, 1);
        return 1;
    }

    int nmeta, nsamples_check;
    if (argc == 3) {
        if(rules_init(argv[2], &nmeta, &nsamples_check, &params.meta, 0) != 0) {
            fprintf(stderr, "could not load minority file at path '%s'\n", argv[2]);
            rules_free(params.rules, params.nrules, 1);
            rules_free(params.labels, params.nlabels, 0);
            return 1;
        }
    }

    if(params.nsamples != nsamples_chk) {
        fprintf(stderr, "number of samples in out file (%d) and label file (%d) must match\n", params.nsamples, nsamples_chk);
        rules_free(params.rules, params.nrules, 1);
        rules_free(params.labels, params.nlabels, 0);
        rules_free(params.meta, nmeta, 0);
        return 1;
    }

    if(params.meta && params.nsamples != nsamples_check) {
        fprintf(stderr, "number of samples in out file (%d) and minority file (%d) must match\n", params.nsamples, nsamples_check);
        rules_free(params.rules, params.nrules, 1);
        rules_free(params.labels, params.nlabels, 0);
        rules_free(params.meta, nmeta, 0);
        return 1;
    }

    std::map<int, std::string> curiosity_map;
    curiosity_map[1] = "curiosity";
    curiosity_map[2] = "curious_lb";
    curiosity_map[3] = "curious_obj";
    curiosity_map[4] = "dfs";

    char froot[BUFSZ];
    params.opt_fname = (char*)malloc(sizeof(char) * BUFSZ);
    params.log_fname = (char*)malloc(sizeof(char) * BUFSZ);
    const char* pch = strrchr(argv[0], '/');
    snprintf(froot, BUFSZ, "../logs/for-%s-%s%s-%s-%s-removed=%s-max_num_nodes=%d-c=%.7f-f=%d",
            pch ? pch + 1 : "",
            run_bfs ? "bfs" : "",
            run_curiosity ? curiosity_map[params.curiosity_policy].c_str() : "",
            (params.map_type == 1) ? "with_prefix_perm_map" :
                (params.map_type == 2 ? "with_captured_symmetry_map" : "no_pmap"),
            params.meta ? "minor" : "no_minor",
            params.ablation ? ((params.ablation == 1) ? "support" : "lookahead") : "none",
            params.max_num_nodes, params.c, params.freq);
    snprintf(params.log_fname, BUFSZ, "%s.txt", froot);
    snprintf(params.opt_fname, BUFSZ, "%s-opt.txt", froot);

    run_corels(params);

    if(vstring_malloc)
        free(params.vstring);

    free(params.opt_fname);
    free(params.log_fname);

    if (params.meta) {
        rules_free(params.meta, nmeta, 0);
    }

    rules_free(params.rules, params.nrules, 1);
    rules_free(params.labels, params.nlabels, 0);

    return 0;
}
