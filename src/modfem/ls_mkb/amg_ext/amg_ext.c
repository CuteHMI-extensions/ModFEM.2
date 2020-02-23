#include "amg_ext.h"

struct AMGAlgorithmData amg_algorithm_data;

void init_amg(int coarsening_algorithm, int interpolation_algorithm, double strength_threshold, int levels_number)
{
	amg_algorithm_data.coarsening_algorithm = coarsening_algorithm;
	amg_algorithm_data.interpolation_algorithm = interpolation_algorithm;
	amg_algorithm_data.strength_threshold = strength_threshold;
	amg_algorithm_data.levels_number = levels_number;
}
