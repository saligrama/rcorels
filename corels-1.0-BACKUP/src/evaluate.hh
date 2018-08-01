#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "queue.hh"
#include "evaluate.hh"

#include "rule.h"



/**
    Stores an optimal rule list, need to be paired with a data_t that contains all the rule data
**/
typedef struct rulelist {

    int nrules;

    unsigned short * ids; // ids of the optimal rule list
    int * predictions; // predictions of the optimal rule list

    int default_prediction;

    double c;
} rulelist_t;



/**
    Verbosity usage is common, so explained here:

        value               behavior

        >0              Only fatal errors

        >1              Only most important data
                        (usually, the final results of the function)

        >2              Diverse progress data

**/



/**
    randomize_rule randomizes the truthtable of a rule for a given number of bits, and
    randomize_data randomizes all the rules in a data object and the labels, by first randomizing
    the 0 label and then setting the 1 label to the complement of that
**/

#ifdef GMP

void
randomize_rule(rule_t * rule, int nsamples, gmp_randstate_t state);

void
randomize_data(data_t * data, gmp_randstate_t state);

#else

void
randomize_rule(rule_t * rule, int nsamples);

void
randomize_data(data_t * data);

#endif

/**
    Runs a number of tests, where for each test the following happens:
    1) A random dataset is generated
    2) Its optimal rule list is determined with CORELS
    3) The objective of the optimal rule list determined by CORELS is checked to see if it indeed
       has the objective that is also outputted by CORELS as the optimal objective. This is done
       with the evaluate_data function.
    4) The optimal rule list is determined by the below brute force algorithm, that checks all
       possible rule lists and predictions (up to a certain rulelist length).
    5) The objective determined by brute force is checked to see if it is greater than or equal to that
       determined by CORELS (it could be greater than CORELS because you can set the the brute force algorithm
       to only check rule lists up to a certain length [otherwise the time needed gets ridiculous]). If
       max_num_nodes is set to a small value, CORELS could get a higher objective than brute force even if
       brute force checked only small length rulelists, and this would trigger an error and exit the loop
       (possibly needs improvement here).

    Parameters:
        num_iters - number of tests
        num_rules - number of rules for the random datasets
        num_samples - number of samples for the random datasets
        c - length constant
        b_max_list_len - maximum length of rulelists checked by the brute force algorithm. Set this to 0 to completely
                         skip (disable) any brute force calculations
        ablation - ablation used when running CORELS
        q_cmp - comparison function used by the queue when running CORELS
        node_type - the type of node used by the tree (if using curious_cmp for the comparison function,
                    this must be set to "curious" in order for the tree to actually use curiosity, otherwise usually it's "node")
        useCapturedPMap - whether the captured permutation map should be used as opposed to the prefix permutation map
        max_num_nodes - maximum number of nodes allowed by CORELS when running
        epsilon - tolerance for floating-point comparisons
        seed - random seed for generating the random datasets
        v - verbosity

    Returns:
        On success:
            0

        On failure:
            1
**/
int
run_random_tests(size_t num_iters, int num_rules, int num_samples, double c, int b_max_list_len,
                     int ablation, std::function<bool(Node*, Node*)> q_cmp, const char* node_type, bool useCapturedPMap,
                     size_t max_num_nodes, double epsilon, unsigned long seed, int v);


/**
    Used by run_random_tests, outputs a reader-friendly error dump when one of the tests in run_random_tests fails
    and the loop stops running. It prints the optimal rule list determined by CORELS and brute force, and the objectives
    determined by CORELS, the evaluate function when checking the rule list outputted by CORELS, and brute force

    Paramters:
        data - contains info about rules and labels for the dataset that caused the error
        corels_opt_list - rule ids of optimal rulelist determined by CORELS
        corels_opt_preds - predictions of optimal rulelist determined by CORELS
        brute_opt_list - rule ids of optimal rulelist determined by brute force
        brute_opt_preds - predictions of optimal rulelist determined by brute force
        output_brute - whether the brute force info should be outputted
                       (set to false by run_random_tests when it is given b_max_list_len = 0)
        corels_obj - optimal objective determined by CORELS
        eval_check_obj - objective of optimal rulelist determined by CORELS, calculated by evaluate_data
        brute_obj - optimal objective determined by brute force
        v - verbosity
**/
void
output_error(data_t data, tracking_vector<unsigned short, DataStruct::Tree> corels_opt_list,
                 tracking_vector<bool, DataStruct::Tree> corels_opt_preds,
                 tracking_vector<unsigned short, DataStruct::Tree> brute_opt_list,
                 tracking_vector<bool, DataStruct::Tree> brute_opt_preds, bool output_brute, double corels_obj,
                 double eval_check_obj, double brute_obj, int v);


/**
    Loads the data from a model, out, label, and minor file into a data and rulelist struct

    Parameters:
        out - pointer to struct in which to store data info about all the rules
        opt_out - pointer to struct in which to store the rulelist indicated by the model file
        model_file - file containing optimal rule list and predictions
                     if NULL, only the general data is loaded and not the optimal rule list
        out_file - .out file containing data
        label_file - .label file
        v - verbosity
**/
int
data_init(data_t * out, rulelist_t * opt_out, const char * model_file, const char * out_file, const char * label_file, int v);



/**
    Called in data_init, loads the ids, predictions, nrules, and default prediction information
    of an optimal rule list stored in model_file

    Parameters:
        out - pointer to struct in which to store rule list
        data - information about all the rules and labels
        model_file - file containing optimal rule list and predictions
        v - verbosity

    Returns:
        On success:
            0

        On failure:
            1
**/
int
data_init_model(rulelist_t * out, data_t data, const char * model_file, int v);



// Frees rules and labels arrays
void data_free(data_t model);

// Frees ids and predictions arrays
void rulelist_free(rulelist_t rulelist);




/**
    Calculates the optimal objective from given data by checking every possible rule list and prediction permutation
    Then, it stores the information of the optimal rule list in opt_list

    Parameters:
        data - contains info about the rule and label data
        opt_list - where the optimal rule list is stored
        max_list_len - Max length of rule lists checked
        c - length constant
        v - verbosity

    Returns:
        On success:
            (double) minimum objective for all the data

        On failure:
            -1.0
**/
double
obj_brute(data_t data, rulelist_t * opt_list, int max_list_len, double c, int v);



/**
    Recursive helper function for finding all the possible rule lists, simply takes a prefix and adds every
    possible rule to it (excluding the ones already in the prefix) and evaluates the rulelists created by adding each
    rule with every possible prediction and default prediction as well. If a rulelist with an objective smaller than
    the current minimum objective is found, the minimum objective and optimal rulelist are updated.

    Parameters:
        data - contains rule and label info for the dataset being searched
        min_obj - current minumum objective
        opt_list - current optimal rulelist
        prefix - rulelist prefix to which every remaning rule is added and the subsequently created rulelists evaluated
        max_list_len - maximum length of rulelists that are checked
        c - length constant
        v - verbosity
**/
void
_obj_brute_helper(data_t data, double * min_obj, rulelist_t * opt_list, rulelist_t prefix, int max_list_len, double c, int v);




/**
    Finds the objective of the rule list stored in the given model file when evaluated with the given data files

    Parameters:
        model_file - file containing optimal rule list and predictions
        out_file - .out file containing data
        label_file - .label file
        c - length constant
        v - verbosity

    Returns:
        On success:
            (double) objective of rule list evaluated with given data

        On failure:
            -1.0
**/
double
evaluate(const char * model_file, const char * out_file, const char * label_file, VECTOR *total_captured_correct, double c, int v);




/**
        Same as before, except it takes a data_t object preloaded with the rule and label info
        and a rulelist_t preloaded with the rulelist being evaluated
**/
double
evaluate_data(data_t data, rulelist_t list, VECTOR *total_captured_correct, double c, int v);
