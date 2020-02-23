#ifndef UTH_IO_RESULTS_H_
#define UTH_IO_RESULTS_H_
#ifdef __cplusplus
extern "C"
{
#endif

/**
 \defgroup UTM_IO_RESULTS Input-Ouput Results Utilities
 \ingroup UTM
 @{
 *
 * Example control_points block (can be in any libconfig file):
 *
 * control_points: (
 * {
 * coords = [0.0, 0.0, 0.0];
 * name = "Point 01";
 * },
 * {
 * coords = [0.2, 0.4, 0.6];
 * name = "Point H_68";
 * }
 * );
 */

///  UTR_IO_RESULT (file with infomrations)

typedef enum {
    RESULT_STEP = 0,
    RESULT_NON_LIN_STEP,
    RESULT_CUR_TIME,
    RESULT_DT,
    RESULT_CFL_MIN,
    RESULT_CFL_MAX,
    RESULT_CFL_AVG,
    RESULT_N_DOFS,
    RESULT_NORM_PLAST_FLOW,
    RESULT_NORM_NS_SUPG,
    RESULT_NORM_NS_SUPG_NONL,
    RESULT_NORM_NS_MIXED,
    RESULT_NORM_NS_MIXED_NONL,
    RESULT_NORM_HEAT,
    RESULT_NORM_HEAT_NONL,
    RESULT_OUT_FILENAME,
    RESULT_TIME_TO_SOL, // time-to-solution
    RESULT_TIME_DDM, // domain decomposition
    RESULT_TIME_PCM, // parallel communication time
    RESULT_TIME_MMPM, // parallel mesh module time
    RESULT_TIME_APPM, // parallel aproximation module
    RESULT_TIME_SIM, // solver interface module
    RESULT_CONTROL_POINTS,
    RESULT_LAST
} utt_io_result_name;

int utr_io_result_set_filename(const char* Workdir,const char* Filename);

int utr_io_result_clear_columns();

void utr_io_result_clear_rows();

int	utr_io_result_add_column(utt_io_result_name column);


int	utr_io_result_write_columns(const int Keep_old_data // 1+ = keep old data in file
                                                        // 0 = clear file
                                );

int utr_io_result_add_value(utt_io_result_name name, const char* value);

int utr_io_result_add_value_int(utt_io_result_name name, int value);

int utr_io_result_add_value_double(utt_io_result_name name, double value);

int utr_io_result_cumulative_timer_start(utt_io_result_name name);
int utr_io_result_cumulative_timer_stop(utt_io_result_name name);

int utr_io_result_write_values_and_proceed();
int utr_io_result_ctrl_pts_write_values_to_file_and_proceed();

/////////////////////////////////////////////////////////////////////////////////
/// For statistics data in multi-processor computations
int utr_io_result_gather_and_avg_all_files(const int N_procs, const int My_proc_ID,
                                           const int Row_nr, const int Keep_data);

/////////////////////////////////////////////////////////////////////////////////
/// For getting values at given points in all time steps.

///
/// \brief utr_io_results_ctrl_pts_read_file - read control points data from input file
/// \return
///
int utr_io_result_ctrl_pts_read_file(const char* Work_dir, const char* Filename);
int utr_io_result_ctrl_pts_set_sol_len(const int N_solutions);
int utr_io_result_ctrl_pts_set_filename(const char* Workdir,const char* Filename);
int utr_io_result_add_ctrl_pt_values(const int Point_nr, const int N_val, double* Values);

int utr_io_result_ctrl_pts_nr();
int utr_io_result_ctrl_pt_get_coords(const int Point_nr, double* Coords);
int utr_io_result_ctrl_pt_set_elemID(const int Pt_nr, const int Mesh_id, const int ElemID);
int utr_io_result_ctrl_pt_get_elemID(const int Pt_nr, const int Mesh_id);


/////////////////////////////////////////////////////////////////////////////////

/** @}*/

#ifdef __cplusplus
}
#endif

/** @} */ // end of group


#endif //UTH_IO_RESULTS_H_
