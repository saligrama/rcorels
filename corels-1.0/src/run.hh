#pragma once

#include "params.h"
#include "alloc.hh"

#ifdef __cplusplus
extern "C" {
#endif

double run_corels (run_params_t params, tracking_vector<unsigned short, DataStruct::Tree> rulelist, tracking_vector<bool, DataStruct::Tree> preds);

#ifdef __cplusplus
}
#endif
