#ifndef AMG_EXT_H
#define AMG_EXT_H

#ifdef __cplusplus
extern "C" {
#endif

#include "../lsd_ns_supg_ext/lsh_ns_supg_ext_intf.h"

struct AMGAlgorithmData {
	int coarsening_algorithm;
	int interpolation_algorithm;
	double strength_threshold;
	int levels_number;
};

extern struct AMGAlgorithmData amg_algorithm_data;

void init_amg(int coarsening_algorithm, int interpolation_algorithm, double strength_threshold, int levels_number);
extern double amg_get_global_diff(double diff);
extern void amg_project_solution_to_level(int level_id);
extern void amg_project_solution_from_level(int level_id);

#ifdef __cplusplus
}
#endif


#endif
