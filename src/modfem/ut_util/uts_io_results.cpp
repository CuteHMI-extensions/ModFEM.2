
#include <vector>
#include <sstream>
#include <iomanip>
#include <ctime>
#include <fstream>
#include <limits>
#include <libconfig.h>
#include <map>
#include <algorithm>


#include <modfem/uth_io_results.h>
#include <modfem/uth_log.h>
#include <modfem/uth_system.h>


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////// UTR_IO_RESULT (file with infomrations)

const char        *       utv_io_result_column_names[RESULT_LAST]
	= { "STEP",
					"NLIN",
					"  TIME",
					"DTIME",
					"CFL_MIN",
					"CFL_MAX",
					"CFL_AVG",
					"N_DOFS",
					"NORM_PLAS_FL",
					"NORM_NS_SUPG",
					"NORM_NS_SUPG_NONL",
					"NORM_HEAT",
					"NORM_HEAT_NONL",
					"FILES",
					"TIME_TO_SOL", // time-to-solution
					"TIME_DDM", // domain decomposition
					"TIME_PCM", // communication time
					"TIME_MMPM", // parallel mesh module time
					"TIME_APPM", // parallel aproximation module
					"TIME_SIM",
					"CONTROL_POINTS"
			};

const int utv_io_result_column_width[RESULT_LAST]
= {
	4,
	4,
	6,
	6,
	6,
	6,
	6,
	6,
	6,
	6,
	6,
	6,
	6,
	6,
	6,
	6,
	6,
	6,
	6,
};

struct utt_io_result_ctrl_pt
{
	static int n_solutions;
	double x, y, z;
	std::string name;
	/// map from Mesh_id -> elem_id in this mesh with this point.
	std::map<int, int>   mesh2elem;
	std::vector<double> solution_values;
};

int utt_io_result_ctrl_pt::n_solutions = 1;

std::vector<std::string> utv_io_result_record(RESULT_LAST);
std::vector<bool>   utv_io_result_columns(RESULT_LAST);
std::vector<double> utv_io_result_times(RESULT_LAST);
std::string         utv_io_result_filename("result_test.csv");


std::vector<utt_io_result_ctrl_pt> utv_io_result_ctrl_pts;
std::string                        utv_io_result_ctrl_pts_filename("control_points.csv");

#ifdef __cplusplus
extern "C"
{
#endif


// Forward declaration.
int utr_io_result_ctrl_pts_write_values_to_file_and_proceed();
int utr_io_result_ctrl_pts_write_columns(const int Keep_old_data); // 1+ = keep old data in file
// 0 = clear file

int utr_io_result_set_filename(const char * Workdir, const char * Filename)
{
	utv_io_result_filename = Workdir;
	utv_io_result_filename.append("/");
	utv_io_result_filename.append(Filename);

	return 0;
}

void utr_io_result_clear_rows()
{
	for (int i = 0; i < RESULT_LAST; ++i) {
		utv_io_result_record[i].clear();
	}
}

int utr_io_result_clear_columns()
{
	utr_io_result_clear_rows();
	for (int i = 0; i < RESULT_LAST; ++i) {
		utv_io_result_columns[i] = false;
		utv_io_result_times[i] = 0.0;
	}
	return 0;
}

int	utr_io_result_add_column(utt_io_result_name column)
{
	utv_io_result_columns[column] = true;
	return 0;
}

int utr_io_result_add_value(utt_io_result_name name, const char * value)
{
	if (! utv_io_result_record[name].empty() ) {
		utv_io_result_record[name].append(",");
	}
	utv_io_result_record[name].append(value);
	return 0;
}

int utr_io_result_add_value_int(utt_io_result_name name, int value)
{
	std::ostringstream strs;
	strs << std::setw(utv_io_result_column_width[name]) << value;
	return utr_io_result_add_value(name, strs.str().c_str());
}

int utr_io_result_add_value_double(utt_io_result_name name, double value)
{
	std::ostringstream strs;
	strs << std::setw(utv_io_result_column_width[name]) << value;
	return utr_io_result_add_value(name, strs.str().c_str());
}

int utr_io_result_cumulative_timer_start(utt_io_result_name name)
{
	utv_io_result_times[name] = time_clock();
	return 0;
}

int utr_io_result_cumulative_timer_stop(utt_io_result_name name)
{
	const double dt = time_clock() - utv_io_result_times[name];
	double prev_val = 0.0;

	mf_check(dt >= 0.0, "Measured negative time of execution!");

	if (! utv_io_result_record[name].empty() ) {
		std::stringstream strs;
		strs << utv_io_result_record[name];
		strs >> prev_val;
		utv_io_result_record[name].clear();
	}

	return utr_io_result_add_value_double(name, prev_val + dt);

}

int utr_io_result_write_columns_to_file(const int Keep_old_data, // 1+ = keep old data in file
		// 0 = clear file
		const std::string & Filename)
{
	std::ios_base::openmode mode = std::fstream::out;

	if (Keep_old_data != 0) {
		mode |= std::fstream::app;
	}

	std::ofstream   file( Filename.c_str(), mode );

	if (Keep_old_data == 0) { // write to empty file
		file << " TIMESTAMP ";
		for (int i = 0; i < RESULT_LAST; ++i) {
			if (utv_io_result_columns.at(i) == true) {
				file << ", " << utv_io_result_column_names[i];
			}
		}
	}

	file << std::endl << std::endl;
	file.close();


	utr_io_result_ctrl_pts_write_columns(Keep_old_data);


	return 0;
}

int	utr_io_result_write_columns(const int Keep_old_data // 1+ = keep old data in file
		// 0 = clear file
)
{
	return utr_io_result_write_columns_to_file(Keep_old_data, utv_io_result_filename);
}

int utr_io_result_write_values_to_file_and_proceed(const std::string & Filename)
{
	utr_io_result_ctrl_pts_write_values_to_file_and_proceed();

	char date[64];
	time_t t = time(0);
	strftime(date, sizeof(date), "%Y-%m-%d %H:%M:%S", gmtime(&t));
	std::ofstream   file(Filename.c_str(), std::fstream::out | std::fstream::app);
	file << date ;
	for (int i = 0; i < RESULT_LAST; ++i) {
		if (utv_io_result_columns[i] == true) {
			// if field is multivalue, then wrap it with ""
			if (utv_io_result_record[i].find(',') != std::string::npos) {
				file <<  ", " << "\"" << utv_io_result_record[i] << "\"";
			}
			else {
				file <<  ", " << utv_io_result_record[i] ;
			}
			utv_io_result_record[i].clear();
			utv_io_result_times[i] = 0.0;
		}
	}
	file << std::endl;
	file.close();

	return 0;
}

int utr_io_result_write_values_and_proceed()
{
	return utr_io_result_write_values_to_file_and_proceed(utv_io_result_filename);
}

int utr_io_result_read_row_from_file(const std::string & Filename, const int row_nr, double row[RESULT_LAST])
{
	int n_read = 0;
	std::ifstream   file(Filename.c_str());

	if (file.is_open()) {

		std::string     line;
		for (int l = 0; l < row_nr; ++l) {
			std::getline(file, line);
		}

		std::stringstream sline(line);

		// Ignore timestamp;
		std::string cell;
		std::getline(sline, cell, ',');
		//ss.ignore(20, ',');

		for (int v = 0; v < RESULT_LAST; ++v) {
			if (utv_io_result_columns[v]) {
				std::getline(sline, cell, ',');
				std::stringstream scell(cell);
				scell >> row[v];
				++n_read;
				mf_log_info("Read %lf", row[v]);
			}
			//ss.ignore(line.size(), ',');
		}

		file.close();
	}

	if (n_read == 0) {
		mf_log_err("Unable to gather data from file %s.", Filename.c_str());
	}

	return n_read;
}

int utr_io_result_gather_and_avg_all_files(const int N_procs, const int My_proc_ID,
		const int Row_nr, const int Keep_data)
{
	// Assumnig other process filenames will be like ours,
	// but with other proc id.
	std::stringstream ss1;
	ss1 << My_proc_ID;
	const std::string & my_proc_id_sting = ss1.str();

	// create file for avg values
	std::string avg_filename = utv_io_result_filename;
	// find last last symbol matching out proc number
	int proc_id_pos = avg_filename.find_last_of(my_proc_id_sting);

	avg_filename.replace(proc_id_pos,
			my_proc_id_sting.size(),
			"AVG");

	// write column names
	utr_io_result_write_columns_to_file(Keep_data, avg_filename);
	// from each process csv file, gather values and add to avg
	double min_val[RESULT_LAST];
	double avg_val[RESULT_LAST];
	double max_val[RESULT_LAST];
	for (int i = 0; i < RESULT_LAST; ++i) {
		min_val[i] = std::numeric_limits<double>::max();
		avg_val[i] = 0.0;
		max_val[i] = std::numeric_limits<double>::min();
	}


	for (int p = 1; p <= N_procs; ++p) {
		std::stringstream ss2;
		ss2 << p;
		const std::string & proc_id_sting = ss2.str();

		std::string p_filename = utv_io_result_filename;
		// find last last symbol matching out proc number
		int proc_id_pos = p_filename.find_last_of(my_proc_id_sting);

		p_filename.replace(proc_id_pos,
				my_proc_id_sting.size(),
				proc_id_sting);

		mf_log_info("Gathering avg. data from file %s.", p_filename.c_str() );

		double p_val[RESULT_LAST] = {0.0};
		utr_io_result_read_row_from_file(p_filename, Row_nr, p_val);

		for (int i = 0; i < RESULT_LAST; ++i) {
			//mf_log_info("Before: min %lf, avg %lf, max %lf",min_val[i],avg_val[i],max_val[i]);
			min_val[i] = std::min(min_val[i], p_val[i]);
			max_val[i] = std::max(max_val[i], p_val[i]);
			avg_val[i] += p_val[i] / double(N_procs);
			//mf_log_info("After: min %lf, avg %lf, max %lf",min_val[i],avg_val[i],max_val[i]);
		}
	}

	utr_io_result_add_value(RESULT_STEP, "Min.");
	utr_io_result_add_value_int(RESULT_STEP, Row_nr);
	for (int i = 0; i < RESULT_LAST; ++i) {
		if (utv_io_result_columns[i]) {
			utr_io_result_add_value_double(utt_io_result_name(i), min_val[i]);
		}
	}
	utr_io_result_write_values_to_file_and_proceed(avg_filename);

	utr_io_result_add_value(RESULT_STEP, "Avg.");
	utr_io_result_add_value_int(RESULT_STEP, Row_nr);
	for (int i = 0; i < RESULT_LAST; ++i) {
		if (utv_io_result_columns[i]) {
			utr_io_result_add_value_double(utt_io_result_name(i), avg_val[i]);
		}
	}
	utr_io_result_write_values_to_file_and_proceed(avg_filename);

	utr_io_result_add_value(RESULT_STEP, "Max.");
	utr_io_result_add_value_int(RESULT_STEP, Row_nr);
	for (int i = 0; i < RESULT_LAST; ++i) {
		if (utv_io_result_columns[i]) {
			utr_io_result_add_value_double(utt_io_result_name(i), max_val[i]);
		}
	}
	utr_io_result_write_values_to_file_and_proceed(avg_filename);

	return 0;
}

/////////////////////////////////////////////////////////////////////////////////
/// For getting values at given points in all time steps.

///
/// \brief utr_io_result_ctrl_pts_read_file - read control points data from input file
/// \return
///
int utr_io_result_ctrl_pts_read_file(const char * Work_dir, const char * Filename)
{
	config_t cfg;
	config_setting_t * root_setting;
	const char * str;
	char ctrl_pts_file[300];
	int ret_val = EXIT_SUCCESS;

	//printf("%s\n",Work_dir);
	//printf("%s\n",Filename);
	sprintf(ctrl_pts_file, "%s/%s", Work_dir, Filename);
	//printf("%s\n",material_file);

	config_init(&cfg);

	/* Read the file. If there is an error, report it and exit. */
	if (! config_read_file(&cfg, ctrl_pts_file)) {
		//fprintf(Interactive_output, "Material file error: %s:%d - %s\n", config_error_file(&cfg),
		//        config_error_line(&cfg), config_error_text(&cfg));
		//mf_log_err("Could not open file! (%s) ", mapping_file);
		mf_log_err("Control points file error: %s:%d - %s\n", config_error_file(&cfg),
				config_error_line(&cfg), config_error_text(&cfg));
		config_destroy(&cfg);
		return (-1);
	}

	root_setting = config_lookup(&cfg, "control_points");

	if (root_setting != NULL)
	{
		int count = config_setting_length(root_setting);
		int setting_length;

		mf_log_info("Number of control points assigments in %s: %d", Filename, count);

		for (int i = 0; i < count; ++i) {
			utt_io_result_ctrl_pt ctrl_pt;
			config_setting_t * ctrl_pt_setting = config_setting_get_elem(root_setting, i);
			config_setting_t * setting;

			int matnum = 0;

			// COORDS
			if (NULL != (setting = config_setting_get_member(ctrl_pt_setting, "coords"))) {
				if (config_setting_is_array(setting) == CONFIG_TRUE) {
					ctrl_pt.x = (double)config_setting_get_float_elem(setting, 0);
					ctrl_pt.y = (double)config_setting_get_float_elem(setting, 1);
					ctrl_pt.z = (double)config_setting_get_float_elem(setting, 2);
				}
				else {
					mf_fatal_err("In 'control_points' in file %s 'coords' is not an array (e.g.[0.0,2.0,4.0]) in assigment %d ",
							ctrl_pts_file, i);
				}
			}
			else {
				mf_fatal_err("In 'control_points' in file %s missing 'coords' array in assigment %d ",
						ctrl_pts_file, i);
			}

			// name
			if (NULL != (setting = config_setting_get_member(ctrl_pt_setting, "name"))) {
				ctrl_pt.name = config_setting_get_string(setting);

				if (ctrl_pt.name.empty()) {
					mf_fatal_err("In 'control_points' in file %s 'name' is empty in assigment %d ",
							ctrl_pts_file, i);
				}
			}
			else {
				mf_fatal_err("In 'control_points' in file %s missing 'name' in assigment %d ",
						ctrl_pts_file, i);
			}

			utv_io_result_ctrl_pts.push_back(ctrl_pt);

		} //!for
	}//!root_setting

	mf_log_info("'control_points' read! ");

	config_destroy(&cfg);

	return (ret_val);

}


int utr_io_result_ctrl_pts_nr()
{
	return utv_io_result_ctrl_pts.size();
}

int utr_io_result_ctrl_pt_get_coords(const int Point_nr, double * Coords)
{
	if (Coords != NULL) {
		Coords[0] = utv_io_result_ctrl_pts[Point_nr].x;
		Coords[1] = utv_io_result_ctrl_pts[Point_nr].y;
		Coords[2] = utv_io_result_ctrl_pts[Point_nr].z;
	}
	return 1;
}

int utr_io_result_ctrl_pt_get_elemID(const int Pt_nr, const int Mesh_id)
{
	if (utv_io_result_ctrl_pts[Pt_nr].mesh2elem.find(Mesh_id)
			== utv_io_result_ctrl_pts[Pt_nr].mesh2elem.end()) {
		mf_fatal_err("Control point element not set!");
	}
	return utv_io_result_ctrl_pts[Pt_nr].mesh2elem[Mesh_id];
}

int utr_io_result_ctrl_pt_set_elemID(const int Pt_nr, const int Mesh_id, const int ElemID)
{
	mf_check(ElemID > 0, "Wrong element ID (%d)!", ElemID );
	utv_io_result_ctrl_pts[Pt_nr].mesh2elem[Mesh_id] = ElemID;
	return 1;
}

int utr_io_result_add_ctrl_pt_values(const int Point_nr, const int N_val, double * Values)
{
	if (Values != NULL) {
		utv_io_result_ctrl_pts[Point_nr].solution_values.insert(
				utv_io_result_ctrl_pts[Point_nr].solution_values.end(),
				Values, Values + N_val);
	}
	return 1;
}

int utr_io_result_ctrl_pts_set_sol_len(const int N_solutions)
{
	utt_io_result_ctrl_pt::n_solutions = N_solutions;
	return 1;
}

int utr_io_result_ctrl_pts_write_columns(const int Keep_old_data) // 1+ = keep old data in file
// 0 = clear file)
{
	if (utv_io_result_columns[RESULT_CONTROL_POINTS]
			&& !utv_io_result_ctrl_pts.empty()) {
		std::ios_base::openmode mode = std::fstream::out;

		if (Keep_old_data != 0) {
			mode |= std::fstream::app;
		}

		std::ofstream   file( utv_io_result_ctrl_pts_filename.c_str(), mode );

		if (Keep_old_data == 0) { // write to empty file

			std::stringstream row1, row2;
			// First row.
			row1 << "STEP, TIME[s], ";
			row2 << ", , ";
			for (int i = 0; i < utv_io_result_ctrl_pts.size(); ++i) {
				for (int j = 0; j < utt_io_result_ctrl_pt::n_solutions; ++j) {
					row1 << utv_io_result_ctrl_pts[i].name << ", ";
					row2 << "Sol " << j << ", ";
				}
			}

			mf_log_info("Creating file %s with control points data.",
					utv_io_result_ctrl_pts_filename.c_str());
			file << row1.str() << "\n" << row2.str();
		}

		file << std::endl ;
		file.close();
	}
	return 0;
}

int utr_io_result_ctrl_pts_write_values_to_file_and_proceed()
{
	if (utv_io_result_columns[RESULT_CONTROL_POINTS]
			&& !utv_io_result_ctrl_pts.empty() ) {
		std::ofstream   file(utv_io_result_ctrl_pts_filename.c_str(), std::fstream::out | std::fstream::app);

		file << utv_io_result_record[RESULT_STEP] << ", " << utv_io_result_record[RESULT_CUR_TIME] << ", ";

		for (int i = 0; i < utv_io_result_ctrl_pts.size(); ++i) {
			const int len = std::min(utt_io_result_ctrl_pt::n_solutions,
							(int)utv_io_result_ctrl_pts[i].solution_values.size());
			for (int j = 0; j < len; ++j) {
				file << utv_io_result_ctrl_pts[i].solution_values[j] << ", ";
			}
			utv_io_result_ctrl_pts[i].solution_values.clear();
		}
		file << std::endl;
		file.close();
	}
	return 0;
}

int utr_io_result_ctrl_pts_set_filename(const char * Workdir, const char * Filename)
{
	utv_io_result_ctrl_pts_filename = Workdir;
	utv_io_result_ctrl_pts_filename.append("/");
	utv_io_result_ctrl_pts_filename.append(Filename);

	return 0;
}


/////////////////////////////////////////////////////////////////////////////////


#ifdef __cplusplus
}
#endif
