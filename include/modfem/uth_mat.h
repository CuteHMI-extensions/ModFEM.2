#ifndef UTS_MAT_H
#define UTS_MAT_H

#ifndef ANALYTIC_EXPRESSIONS
#define ANALYTIC_EXPRESSIONS
#endif

#include <stdlib.h>

#ifdef __cplusplus
extern "C"
{
#endif


/**
 \defgroup UTM_MAT Materials Utilities
 \ingroup UTM
*
* The materials in ModFEM are defined as \link utr_material_query stucture of different material properties \endlink.
* The materials submodule of the utilities module is responsible for \link utr_mat_read reading the material file(s)  \endlink
* and \link utr_material_query providing access to the exisiting database \endlink.
* The default and preferred way of accessing the materials is to make use of \link utr_material_query material query function \endlink with \link mmr_el_groupID groupID \endlink
* as a parameter.
*
* The correct way of using the material database is as follows:
*
* At program initialisation:
* 1. \link utr_mat_read Read the material file(s) \endlink at least once.
*   - the \link utc_mat_database_filename main material base \endlink will be automaticly read
*   - provided material file should contain information about \link utr_mat_get_matID group-to-material \endlink and \link utr_mat_get_blockID group-to-block \endlink information.
*     example blocks with those informations could look like this:
*
*    group_to_materials_assignment: (
*    {
*    material_number = 10;
*    groups = (15,11);
*    },
*    {
*    material_number = 5;
*    groups = ( 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 );
*    });
*
*
* During computations:
* 2. Get the \link mmr_el_groupID groupID of the element \endlink you are analysing.
* 3. Create \link utt_material_query_params structure with params \endlink of the material query.
* - obligatory paramteres are:
*   + \link utt_material_query_params.group_idx groupID \endlink of the element you are querying
*   + \link utt_material_query_params.temperature current temperature \endlink if you want to use material properties depending on temperature
* - Additional parameters are:
*   + \link utt_query_type type of query \endlink
*   + \link utt_material_query_params.xg coordinates of point of interest \endlink
* 4. Create \link utt_material_query_result query result \endlink and initialize it with \link utr_mat_init_query_result \endlink.
* 5. Select material properties that you want to obtain:
*  - for each desired field of initialized \link utt_material_query_result query result structure \endlink set it to the value \link UTC_MAT_QUERY_RESULT_REQUIRED \endlink.
* 6. Finally call \link utr_material_query \endlink providing query parameters from step 3 and query result from step 4-5 as arguments.
* 7. As function successfully returns fields from the \link utt_material_query_result result \endlink previously assigned with
* \link UTC_MAT_QUERY_RESULT_REQUIRED \endlink will store information about selected material properties.
*
* Examples of the material query usage can be found in \link pdr_heat_material_query heat module \endlink, \link pdr_ns_supg_material_query ns_supg module \endlink
* and \link pdr_ns_supg_heat_vof_material_query ns_supg-heat-vof supermodule \endlink.
*
* NOTE: The materials submodule provides also \link utr_mat_get_material direct access \endlink to the \link utt_material_data raw material data \endlink,
* but it is intended to be an additional feature for special purposes,
* and using this other way is much more cumbersome during computations.
*
 @{
 */

/// Filename of the default built-in material database linked with executable.
/// The default location of this file is the directory of the executable.
/// If given problem require to provide additional materials, then the material file
/// given in problem file should be modyfied.
extern const char* utc_mat_database_filename;

/** internal structure filled from materials file - not for direct access!*/
#ifndef PHASE_TRANSFORMATION
typedef struct {

  int matnum; ///< REQUIRED: arbitrary number of material. Should be bigger then \link utc_mat_database_max_reserved_id \endlink if material defined outside file \link utc_mat_database_filename \endlink.
  char name[80];  ///< REQUIRED: name of material
  char is_fluid;    ///< (optional) default: "FLUID", possible: "FLUID","NOT_FLUID"

  double *Tfor_dynamic_viscosity;
  double *atT_dynamic_viscosity;
  void*	expression_of_dynamic_viscosity_from_T;
  int dynamic_viscosity_num;

  double *Tfor_density;
  double *atT_density;
  void*	expression_of_density_from_T;
  int density_num;

  double *Tfor_thermal_conductivity;
  double *atT_thermal_conductivity;
  int thermal_conductivity_num;

  double *Tfor_specific_heat;
  double *atT_specific_heat;
  int specific_heat_num;

  double *Tfor_enthalpy;
  double *atT_enthalpy;
  int enthalpy_num;

  double *Tfor_thermal_expansion_coefficient;
  double *atT_thermal_expansion_coefficient;
  int thermal_expansion_coefficient_num;

  double *Tfor_electrical_resistivity;
  double *atT_electrical_resistivity;
  int electrical_resistivity_num;

  double *Tfor_VOF;
  double *atT_VOF;
  int VOF_num;

  double *Tfor_dg_dT;
  double *atT_dg_dT;
  int dg_dT_num;

  double temp_solidus;
  double temp_liquidus;
  double temp_vaporization;

  double latent_heat_of_fusion;
  double latent_heat_of_vaporization;

  double composition;

  double *TSfor_flow_stress;	// plastic Strain for flow stress (Temperature = constant)
  double *atTS_flow_stress;		// flow stress at plastic Strain (Temperature = constant)
  void*	expression_of_flow_stress_from_S;	// flow_stress = f(strain_intensity)
  int flow_stress_num;

  double Young_mod;
  double Poisson_rat;
  double bulk_mod;
  double shear_mod;

} utt_material_data;
#else
typedef struct {

  int matnum;
  char name[20];
  char is_fluid;    ///< (optional) default: "FLUID", possible: "FLUID","NOT_FLUID"

  double *Tfor_dynamic_viscosity;
  double *atT_dynamic_viscosity;
  void*	expression_of_dynamic_viscosity_from_T;
  int dynamic_viscosity_num;

  double *Tfor_density;
  double *atT_density;
  void*	expression_of_density_from_T;
  int density_num;

  double *Tfor_thermal_conductivity;
  double *atT_thermal_conductivity;
  int thermal_conductivity_num;

  double *Tfor_specific_heat;
  double *atT_specific_heat;
  int specific_heat_num;

  double *Tfor_enthalpy;
  double *atT_enthalpy;
  int enthalpy_num;

  double *Tfor_thermal_expansion_coefficient;
  double *atT_thermal_expansion_coefficient;
  int thermal_expansion_coefficient_num;

  double *Tfor_electrical_resistivity;
  double *atT_electrical_resistivity;
  int electrical_resistivity_num;

  double *Tfor_VOF;
  double *atT_VOF;
  int VOF_num;

  double *Tfor_dg_dT;
  double *atT_dg_dT;
  int dg_dT_num;

  double temp_solidus;
  double temp_liquidus;
  double temp_vaporization;

  double latent_heat_of_fusion;
  double latent_heat_of_vaporization;

  double composition;

  double *TSfor_flow_stress;	// plastic Strain for flow stress (Temperature = constant)
  double *atTS_flow_stress;		// flow stress at plastic Strain (Temperature = constant)
  void*	expression_of_flow_stress_from_S;	// flow_stress = f(strain_intensity)
  int flow_stress_num;

  double Young_mod;
  double Poisson_rat;
  double bulk_mod;
  double shear_mod;

  double *tfor_T_Fs;	// time to estimate the temperature of phase transformation beginning
  double *att_T_Fs;		// estimated temperature of phase transformation beginning
  void*	expression_of_T_Fs_from_t;
  int T_Fs_num;

  double *tfor_T_Ff;
  double *att_T_Ff;
  void*	expression_of_T_Ff_from_t;
  int T_Ff_num;

  double *tfor_T_Ps;
  double *att_T_Ps;
  void*	expression_of_T_Ps_from_t;
  int T_Ps_num;

  double *tfor_T_Pf;
  double *att_T_Pf;
  void*	expression_of_T_Pf_from_t;
  int T_Pf_num;

  double *tfor_T_Bs;
  double *att_T_Bs;
  void*	expression_of_T_Bs_from_t;
  int T_Bs_num;

  double *tfor_T_Bf;
  double *att_T_Bf;
  void*	expression_of_T_Bf_from_t;
  int T_Bf_num;

  double *tfor_T_Ms;
  double *att_T_Ms;
  void*	expression_of_T_Ms_from_t;
  int T_Ms_num;

  double *tfor_T_Mf;
  double *att_T_Mf;
  void*	expression_of_T_Mf_from_t;
  int T_Mf_num;

  double *Vcfor_M_Vc;
  double *atVc_M_Vc;
  int M_Vc_num;

  double *Vcfor_B_Vc;
  double *atVc_B_Vc;
  int B_Vc_num;

  double *Vcfor_F_Vc;
  double *atVc_F_Vc;
  int F_Vc_num;

  double *Vcfor_P_Vc;
  double *atVc_P_Vc;
  int P_Vc_num;

  double *Tfor_H_A_F;
  double *atT_H_A_F;
  void*	expression_of_H_A_F_from_T;
  int H_A_F_num;

  double *Tfor_H_A_P;
  double *atT_H_A_P;
  void*	expression_of_H_A_P_from_T;
  int H_A_P_num;

  double *Tfor_H_A_B;
  double *atT_H_A_B;
  void*	expression_of_H_A_B_from_T;
  int H_A_B_num;

  double *Tfor_H_A_M;
  double *atT_H_A_M;
  void*	expression_of_H_A_M_from_T;
  int H_A_M_num;

  double tB_lim;
  double tM_lim;

} utt_material_data;
#endif

typedef enum {
    QUERY_NODE, ///< if query concerns exact mesh node (vertex) of element
    QUERY_POINT ///< if query concerns any arbitrary selected point within given element
} utt_query_type;

typedef enum {
    NOT_FLUID=0,
    FLUID=1
} utt_is_fluid;

/** input for pdr_material_query */
typedef struct{
  char * name; ///< if searching by idx (idx >= 0) set to anything
  int group_idx; ///< set -1 if you want to search by name
  double temperature; ///< needed for temperature-depended material properties
  double strain_intensity;
  double xg[3]; ///< coordinates of the point
  utt_query_type query_type; ///< \link utt_query_type mesh node or any point \endlink
  int cell_id;  ///< element ID
  int node_id;  ///< node (vertex) ID
  double time; /// needed to query for current simulation time
  double Vc;  // query for cooling rate
  //in future: elem no., coordinates etc.
} utt_material_query_params;

/// \brief Result of the \link utr_material_query \endlink .
/// NOTE: if parameter is required, it have to be set to \link UTC_MAT_QUERY_RESULT_REQUIRED \endlink before executing  \link utr_material_query \endlink !.
/// By default all parameters are set to  \link UTC_MAT_QUERY_RESULT_NOT_NEEDED \endlink

#ifndef PHASE_TRANSFORMATION
typedef struct{
    char name[80];
    char is_fluid;
    double dynamic_viscosity;
    double density;
    double thermal_conductivity;
    double specific_heat;
    double enthalpy;
    double thermal_expansion_coefficient;
    double electrical_resistivity;
    double VOF;
    double dg_dT;
    double temp_solidus;
    double temp_liquidus;
    double temp_vaporization;
    double latent_heat_of_fusion;
    double latent_heat_of_vaporization;

    double composition;
    double flow_stress;
    double Young_mod;
    double Poisson_rat;
    double bulk_mod;
    double shear_mod;
} utt_material_query_result;
#else
typedef struct{
    char name[20];
    char is_fluid;
    double dynamic_viscosity;
    double density;
    double thermal_conductivity;
    double specific_heat;
    double enthalpy;
    double thermal_expansion_coefficient;
    double electrical_resistivity;
    double VOF;
    double dg_dT;
    double temp_solidus;
    double temp_liquidus;
    double temp_vaporization;
    double latent_heat_of_fusion;
    double latent_heat_of_vaporization;

    double composition;
    double flow_stress;		// yield stress
    double Young_mod;
    double Poisson_rat;
    double bulk_mod;
    double shear_mod;
	
	double T_Fs;
	double T_Ff;
	double T_Ps;
	double T_Pf;
	double T_Bs;
	double T_Bf;
	double T_Ms;
	double T_Mf;
	
	double M_Vc;
	double B_Vc;
	double F_Vc;
	double P_Vc;
	
	double H_A_F;
	double H_A_P;
	double H_A_B;
	double H_A_M;
	
	double tB_lim;
	double tM_lim;
	
} utt_material_query_result;
#endif

/// Used to indicate not needed parameters in \link utt_material_query_result \endlink .
#define UTC_MAT_QUERY_RESULT_NOT_NEEDED -1.0
/// Used to indicate required parameters in \link utt_material_query_result \endlink .
#define UTC_MAT_QUERY_RESULT_REQUIRED -3.0

/// For the purpose of backward-compatibility and support of auto-assigning to groups, a few additional flags are provided.
typedef enum {
    UTE_GROUP_ID_NOT_ASSIGNED = -1, ///< to indicate, that given group ID is not connected with any material or block
    UTE_GROUP_ID_DEFAULT_BLOCK = 1  ///< to indicate auto-assignment for so-called "default block" which is created when no blocks have been defined
} ute_mat_groupID_flag;

typedef enum {
    UTE_SUCCESS = EXIT_SUCCESS,
    UTE_FAIL = -1,
    UTE_REF_TEMP_MUST_BE_GEQ_ZERO = -2,
    UTE_NOT_NEED_MATERIAL_DATABASE = 1
} ute_mat_read_result;

// Main functions

ute_mat_read_result utr_mat_read(const char *Work_dir,
                 const char *Filename,
                 FILE *Interactive_output);

void utr_mat_init_query_result(utt_material_query_result* result);

/** \brief Main function for queering for material data.
 *
 * \param Params
 * \param Result
 * \return success code
 *
 */
int utr_material_query(const utt_material_query_params *Params,
                         utt_material_query_result *Result);

/** \brief To get material number associated with given element group number.
 *
 * \param groupID int see \link mmr_el_groupID \endlink
 * \return int material number (ID)
 *
 */
int utr_mat_get_matID(int groupID);



/** \brief To get total number of materials in current database.
 *
 * \return int
 *
 */
int utr_mat_get_n_materials();

/**
 * @brief utr_mat_get_materials_IDs To get array of all known materialIDs
 * @param materialIDs array of size returned by \link utr_mat_get_n_materials \endlink
 * @return number of written material IDs into array \param materialIDs or UTE_FAIL.
 */
int utr_mat_get_materials_IDs(int *materialIDs);


/** \brief To permanently delete stored material and free resources.
 *
 * \return void
 *
 */
void utr_mat_clear_material(utt_material_data* mat);


/** \brief To permanently delete all stored materials.
 *
 * \return void
 *
 */
void utr_mat_clear_all();


// Additional functions

//int utr_read_mat_and_block_mapID(config_setting_t *root_setting);

/** \brief To gain access to raw material data.
 *  NOTE: this is NOT the same as calling \link utr_material_query \endlink !
 *
 * \param groupID const int see \link mmr_el_groupID \endlink
 * \return const utt_material_data*
 *
 */
const utt_material_data* utr_mat_get_material(const int groupID);

/** \brief  To gain access to raw material data.
 *  NOTE: this is NOT the same as calling \link utr_material_query \endlink !
 *
 * \param matID int see see \link utr_mat_get_matID \endlink
 * \return const utt_material_data*
 *
 */
const utt_material_data *utr_mat_get_material_by_matID(int matID);


/** @} */ // end of group

#ifdef __cplusplus
}
#endif

#endif //UTS_MAT_H
