#ifndef UTS_MAT_CPP
#define UTS_MAT_CPP

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libconfig.h>
#include <set>
#include <string>
#include <vector>
#include <algorithm>
#include <map>

#include <modfem/uth_log.h>
#include <modfem/uth_mat.h>
#include <modfem/uth_bc.h>

#include "fastmathparser/exprtk.hpp"

#ifdef __cplusplus
extern "C"
{
#endif


const char * utc_mat_database_filename = "materials_database.dat";


const int utc_mat_database_max_reserved_id = 100;

// LOCAL STRUCUTES
/* main structure containing materials data */
typedef struct {
	typedef std::map<int, utt_material_data> utt_material_map;
	std::vector<int> material_ids;
	utt_material_map material_data;
} utt_materials;



// LOCAL VARIABLES
utt_materials utv_materials;
std::set<std::string> utv_read_material_files;

//
std::vector<int>    utv_group2material;

void utr_mat_init_material(utt_material_data * mat)
{
	if (mat != NULL) {
		memset(mat, 0, sizeof(utt_material_data));
	}
}

void utr_mat_clear_material(utt_material_data * mat)
{
	if (mat != NULL) {
		mat->is_fluid = FLUID;
		UTM_SAFE_FREE_PTR(mat->Tfor_density);
		UTM_SAFE_FREE_PTR(mat->atT_density);

		UTM_SAFE_FREE_PTR(mat->Tfor_thermal_conductivity);
		UTM_SAFE_FREE_PTR(mat->atT_thermal_conductivity);

		UTM_SAFE_FREE_PTR(mat->Tfor_specific_heat);
		UTM_SAFE_FREE_PTR(mat->atT_specific_heat);

		UTM_SAFE_FREE_PTR(mat->Tfor_enthalpy);
		UTM_SAFE_FREE_PTR(mat->atT_enthalpy);

		UTM_SAFE_FREE_PTR(mat->Tfor_thermal_expansion_coefficient);
		UTM_SAFE_FREE_PTR(mat->atT_thermal_expansion_coefficient);


		UTM_SAFE_FREE_PTR(mat->Tfor_electrical_resistivity);
		UTM_SAFE_FREE_PTR(mat->atT_electrical_resistivity);

		UTM_SAFE_FREE_PTR(mat->Tfor_VOF);
		UTM_SAFE_FREE_PTR(mat->atT_VOF);

		UTM_SAFE_FREE_PTR(mat->Tfor_dg_dT);
		UTM_SAFE_FREE_PTR(mat->atT_dg_dT);

		UTM_SAFE_FREE_PTR(mat->TSfor_flow_stress);
		UTM_SAFE_FREE_PTR(mat->atTS_flow_stress);

#ifdef PHASE_TRANSFORMATION
		UTM_SAFE_FREE_PTR(mat->tfor_T_Fs);
		UTM_SAFE_FREE_PTR(mat->att_T_Fs);

		UTM_SAFE_FREE_PTR(mat->tfor_T_Ff);
		UTM_SAFE_FREE_PTR(mat->att_T_Ff);

		UTM_SAFE_FREE_PTR(mat->tfor_T_Ps);
		UTM_SAFE_FREE_PTR(mat->att_T_Ps);

		UTM_SAFE_FREE_PTR(mat->tfor_T_Pf);
		UTM_SAFE_FREE_PTR(mat->att_T_Pf);

		UTM_SAFE_FREE_PTR(mat->tfor_T_Bs);
		UTM_SAFE_FREE_PTR(mat->att_T_Bs);

		UTM_SAFE_FREE_PTR(mat->tfor_T_Bf);
		UTM_SAFE_FREE_PTR(mat->att_T_Bf);

		UTM_SAFE_FREE_PTR(mat->tfor_T_Ms);
		UTM_SAFE_FREE_PTR(mat->att_T_Ms);

		UTM_SAFE_FREE_PTR(mat->tfor_T_Mf);
		UTM_SAFE_FREE_PTR(mat->att_T_Mf);

		UTM_SAFE_FREE_PTR(mat->Vcfor_M_Vc);
		UTM_SAFE_FREE_PTR(mat->atVc_M_Vc);

		UTM_SAFE_FREE_PTR(mat->Vcfor_B_Vc);
		UTM_SAFE_FREE_PTR(mat->atVc_B_Vc);

		UTM_SAFE_FREE_PTR(mat->Vcfor_F_Vc);
		UTM_SAFE_FREE_PTR(mat->atVc_F_Vc);

		UTM_SAFE_FREE_PTR(mat->Vcfor_P_Vc);
		UTM_SAFE_FREE_PTR(mat->atVc_P_Vc);

		UTM_SAFE_FREE_PTR(mat->Tfor_H_A_F);
		UTM_SAFE_FREE_PTR(mat->atT_H_A_F);

		UTM_SAFE_FREE_PTR(mat->Tfor_H_A_P);
		UTM_SAFE_FREE_PTR(mat->atT_H_A_P);

		UTM_SAFE_FREE_PTR(mat->Tfor_H_A_B);
		UTM_SAFE_FREE_PTR(mat->atT_H_A_B);

		UTM_SAFE_FREE_PTR(mat->Tfor_H_A_M);
		UTM_SAFE_FREE_PTR(mat->atT_H_A_M);

#endif
	}
}

ute_mat_read_result utr_mat_read_assigments( const char * Work_dir,
		const char * Filename,
		FILE * Interactive_output)
{
	config_t cfg;
	config_setting_t * root_setting;
	const char * str;
	char mapping_file[300];


	ute_mat_read_result ret_val = UTE_SUCCESS;

	//printf("%s\n",Work_dir);
	//printf("%s\n",Filename);
	sprintf(mapping_file, "%s/%s", Work_dir, Filename);
	//printf("%s\n",material_file);

	config_init(&cfg);

	/* Read the file. If there is an error, report it and exit. */
	if (! config_read_file(&cfg, mapping_file)) {
		//fprintf(Interactive_output, "Material file error: %s:%d - %s\n", config_error_file(&cfg),
		//        config_error_line(&cfg), config_error_text(&cfg));
		//mf_log_err("Could not open file! (%s) ", mapping_file);
		mf_log_err("Material file error: %s:%d - %s\n", config_error_file(&cfg),
				config_error_line(&cfg), config_error_text(&cfg));
		config_destroy(&cfg);
		return (UTE_FAIL);
	}

	root_setting = config_lookup(&cfg, "group_to_materials_assignment");
	if (root_setting != NULL)
	{
		int count = config_setting_length(root_setting);
		int setting_length;

		mf_log_info("Number of material assigments in %s: %d", Filename, count);

		if (count > utv_group2material.size()) {
			utv_group2material.resize(count + 1, UTE_GROUP_ID_NOT_ASSIGNED);

			for (int i = 0; i < count; ++i) {
				config_setting_t * material = config_setting_get_elem(root_setting, i);
				config_setting_t * setting;

				int matnum = 0;

				// MATERIAL NUMBER
				if (NULL != (setting = config_setting_get_member(material, "material_number"))) {
					matnum  = (int)config_setting_get_int(setting);
				}
				else {
					mf_fatal_err("In 'group_to_materials_assignment' in file %s missing 'material_number' in assigment %d ",
							mapping_file, i);
				}

				// groups
				if (NULL != (setting = config_setting_get_member(material, "groups"))) {

					if (config_setting_is_list(setting) == CONFIG_TRUE) {
						setting_length = config_setting_length(setting);
						for (int j = 0; j < setting_length; ++j) {
							int  groupID = config_setting_get_int_elem(setting, j);

							if (groupID >= utv_group2material.size()) {
								utv_group2material.resize(groupID + 1, UTE_GROUP_ID_NOT_ASSIGNED);
							}

							if (utv_group2material[groupID] != UTE_GROUP_ID_NOT_ASSIGNED) {
								mf_log_warn("Duplicate 'group_to_materials_assignment' for groupID=%d ", groupID);
							}

							utv_group2material[groupID] = matnum;

							mf_log_info("Assigning group %d to material %d", groupID, matnum);
						}
					}
					else {
						int  groupID = config_setting_get_int(setting);

						if (groupID >= utv_group2material.size()) {
							utv_group2material.resize(groupID + 1, UTE_GROUP_ID_NOT_ASSIGNED);
						}

						if (utv_group2material[groupID] != UTE_GROUP_ID_NOT_ASSIGNED) {
							mf_log_warn("Duplicate 'group_to_materials_assignment' for groupID=%d ", groupID);
						}

						utv_group2material[groupID] = matnum;

						mf_log_info("Assigning group %d to material %d", groupID, matnum);

					}
				}
				else {
					mf_log_err("In 'group_to_materials_assignment' in file %s missing 'groups' in assigment %d ", mapping_file, i);
				}

			}
		}
		mf_log_info("'group_to_materials_assignment' read! ");
	}

	config_destroy(&cfg);

	return (ret_val);
}


ute_mat_read_result utr_mat_read_material_file( const char * Work_dir,
		const char * Filename,
		FILE * Interactive_output)
{


	utr_mat_read_assigments(Work_dir, Filename, Interactive_output);
	utr_bc_read_block_assigments(Work_dir, Filename, Interactive_output);
	//utr_mat_read_assigments(Work_dir,"materials_database.dat",Interactive_output);


	utv_read_material_files.insert(Filename);

	config_t cfg;
	config_setting_t * root_setting;
	const char * str;
	char material_file[300];

	int constant_thermal_conductivity = 1;
	int constant_density = 1;
	int constant_specific_heat = 1;
	int no_other_parameters = 1;
	int constant_viscosity = 1;
	int constant_flow_stress = 1;

#ifdef PHASE_TRANSFORMATION
	int constant_T_Fs = 1;
	int constant_T_Ff = 1;
	int constant_T_Ps = 1;
	int constant_T_Pf = 1;
	int constant_T_Bs = 1;
	int constant_T_Bf = 1;
	int constant_T_Ms = 1;
	int constant_T_Mf = 1;

	int constant_M_Vc = 1;
	int constant_B_Vc = 1;
	int constant_F_Vc = 1;
	int constant_P_Vc = 1;

	int constant_H_A_F = 1;
	int constant_H_A_P = 1;
	int constant_H_A_B = 1;
	int constant_H_A_M = 1;

#endif

	ute_mat_read_result ret_val = UTE_SUCCESS;

	//printf("%s\n",Work_dir);
	//printf("%s\n",Filename);
	sprintf(material_file, "%s/%s", Work_dir, Filename);
	//printf("%s\n",material_file);

	mf_log_info("Materials configuration file: %s", material_file);

	config_init(&cfg);

	/* Read the file. If there is an error, report it and exit. */
	if (! config_read_file(&cfg, material_file)) {
		//fprintf(Interactive_output, "Material file error: %s:%d - %s\n", config_error_file(&cfg),
		//        config_error_line(&cfg), config_error_text(&cfg));
		mf_log_err("Could not open material file! (%s) ", material_file);
		config_destroy(&cfg);
		return (UTE_FAIL);
	}

	root_setting = config_lookup(&cfg, "materials");
	if (root_setting != NULL)
	{
		utt_materials * Materials_db = &utv_materials;



		int count = config_setting_length(root_setting);
		int i, j;
		int setting_length;
		const char * name;

		mf_log_info("Number of materials: %d", count);

		for (i = 0; i < count; ++i) {
			config_setting_t * material = config_setting_get_elem(root_setting, i);
			config_setting_t * setting;

			// MATERIAL NUMBER
			if (NULL != (setting = config_setting_get_member(material, "material_number"))) {
				int matnum  = (int)config_setting_get_int(setting);

				if ( Materials_db->material_ids.end() !=
						std::find( Materials_db->material_ids.begin(),
								Materials_db->material_ids.end(),
								matnum ) ) {
					mf_fatal_err("Material with id (%d) already exist!", matnum );
				}

//////////////////////////////  //////////////////////////////  //////////////////////////////
/// Built-in material database uses id (0-utc_mat_database_max_reserved_id).
/// Users can define own materials with ids bigger than utc_mat_database_max_reserved_id.
//////////////////////////////  //////////////////////////////  //////////////////////////////
				if (strcmp(Filename, utc_mat_database_filename)) {
					if (matnum <= utc_mat_database_max_reserved_id) {
						mf_log_err("Problem material file (%s) uses id (%d) reserved for built-in material database (0-%d)",
								Filename, matnum, utc_mat_database_max_reserved_id);
					}
				}
				else { // this is main material database
					if (matnum > utc_mat_database_max_reserved_id) {
						mf_log_err("Material from database uses id (%d) reserved for user materials(bigger than %d)",
								matnum, utc_mat_database_max_reserved_id);
					}
				}
//////////////////////////////  //////////////////////////////  //////////////////////////////
				Materials_db->material_ids.push_back(matnum);
			}
			else {
				int matnum = Materials_db->material_data.size(); // because we wont materialid=0 fot only 1 material
				mf_log_warn("Material without 'material_number' asssigned to material %d", matnum);
				Materials_db->material_ids.push_back(matnum);
			}

			utt_material_data & mat_data = Materials_db->material_data[ Materials_db->material_ids.back() ];
			utr_mat_init_material(& mat_data);

			mat_data.matnum = Materials_db->material_ids.back();

			// MATERIAL NAME
			config_setting_lookup_string(material, "name", &name);

			if (strlen(name) >= sizeof(mat_data.name)) {
				mf_fatal_err("Material name too long (%d). Should be max %d",
						(int)strlen(name),
						(int)sizeof(mat_data.name));
			}

			strcpy(mat_data.name, name);

			// IS_FLUID
			if (NULL != (setting = config_setting_get_member(material, "is_fluid"))) {

				const char * buf_is_fluid = NULL;
				config_setting_lookup_string(material,  "is_fluid", & buf_is_fluid );

				if (0 == strcmp(buf_is_fluid, "NOT_FLUID")) {
					mat_data.is_fluid = NOT_FLUID;
				}
				else if (0 == strcmp(buf_is_fluid, "FLUID")) {
					mat_data.is_fluid = FLUID;
				}
				else {
					mf_fatal_err("Incorrect value of 'is_fluid=%s' in file %s (correct are: \"FLUID\" or \"NOT_FLUID\". ",
							buf_is_fluid,
							Filename);
				}
			}
			else { // default:
				mat_data.is_fluid = FLUID;
			}

			// DENSITY
			if (NULL != (setting = config_setting_get_member(material, "density"))) {

				if (config_setting_is_list(setting) == CONFIG_TRUE)
				{
					setting_length = config_setting_length(setting);
					mat_data.Tfor_density =
							(double *)malloc(sizeof(double) * setting_length);
					mat_data.atT_density =
							(double *)malloc(sizeof(double) * setting_length);
					mat_data.density_num = setting_length;
					for (j = 0; j < setting_length; ++j)
					{
						mat_data.Tfor_density[j] =
								(double)config_setting_get_float_elem(config_setting_get_elem(setting, j), 0);
						mat_data.atT_density[j] =
								(double)config_setting_get_float_elem(config_setting_get_elem(setting, j), 1);
					}

					//              // for temperature dependance reference temperature must be > 0
					//              //TODO: move to check problem module
					//              if(pdv_heat_problem.ctrl.ref_temperature <= 0.0){
					//                fprintf(Interactive_output, "ref_temperature <= 0.0 for temperature dependent visocisty. Exiting \n");
					//                config_destroy(&cfg);
					//                exit(-1);
					//                //return(EXIT_FAILURE);
					//              }
					ret_val = UTE_REF_TEMP_MUST_BE_GEQ_ZERO;

					constant_density = 0;

				}
				else if (NULL != config_setting_get_string(setting))
				{
#ifdef ANALYTIC_EXPRESSIONS
					mat_data.density_num = 1;
					mat_data.Tfor_density = (double *)malloc(sizeof(double));
					mat_data.expression_of_density_from_T = (void *) new exprtk::expression<double>();

					std::string expression_string(config_setting_get_string(setting));
					exprtk::symbol_table<double>	symbol_table;
					exprtk::expression<double>	&	expression = *((exprtk::expression<double> *)mat_data.expression_of_density_from_T);

					symbol_table.add_variable("T", mat_data.Tfor_density[0]);
					symbol_table.add_constants();
					expression.register_symbol_table(symbol_table);

					exprtk::parser<double>			parser;
					parser.compile(expression_string, expression);
#endif
				}
				else
				{
					mat_data.density_num = 1;
					mat_data.Tfor_density = (double *)malloc(sizeof(double));
					mat_data.atT_density = (double *)malloc(sizeof(double));
					mat_data.Tfor_density[0] = -1.0;
					mat_data.atT_density[0] =
							(double)config_setting_get_float(setting);


				}
			}

			// FLOW_STRESS
			if (NULL != (setting = config_setting_get_member(material, "flow_stress"))) {

				if (config_setting_is_list(setting) == CONFIG_TRUE)
				{
					setting_length = config_setting_length(setting);
					mat_data.TSfor_flow_stress =
							(double *)malloc(sizeof(double) * setting_length);
					mat_data.atTS_flow_stress =
							(double *)malloc(sizeof(double) * setting_length);
					mat_data.flow_stress_num = setting_length;
					for (j = 0; j < setting_length; ++j)
					{
						mat_data.TSfor_flow_stress[j] =
								(double)config_setting_get_float_elem(config_setting_get_elem(setting, j), 0);
						mat_data.atTS_flow_stress[j] =
								(double)config_setting_get_float_elem(config_setting_get_elem(setting, j), 1);
					}

					//              // for temperature dependance reference temperature must be > 0
					//              //TODO: move to check problem module
					//              if(pdv_heat_problem.ctrl.ref_temperature <= 0.0){
					//                fprintf(Interactive_output, "ref_temperature <= 0.0 for temperature dependent visocisty. Exiting \n");
					//                config_destroy(&cfg);
					//                exit(-1);
					//                //return(EXIT_FAILURE);
					//              }
					ret_val = UTE_REF_TEMP_MUST_BE_GEQ_ZERO;

					constant_flow_stress = 0;

				}
				else if (NULL != config_setting_get_string(setting))
				{
#ifdef ANALYTIC_EXPRESSIONS
					mat_data.flow_stress_num = 1;
					mat_data.TSfor_flow_stress = (double *)malloc(sizeof(double));
					mat_data.expression_of_flow_stress_from_S = (void *) new exprtk::expression<double>();

					std::string expression_string(config_setting_get_string(setting));
					exprtk::symbol_table<double>	symbol_table;
					exprtk::expression<double>	&	expression = *((exprtk::expression<double> *)mat_data.expression_of_flow_stress_from_S);

					symbol_table.add_variable("eps_i", mat_data.TSfor_flow_stress[0]);	// eps_i strain intensity
					symbol_table.add_constants();
					expression.register_symbol_table(symbol_table);

					exprtk::parser<double>			parser;
					parser.compile(expression_string, expression);
#endif
				}
				else
				{
					mat_data.flow_stress_num = 1;
					mat_data.TSfor_flow_stress = (double *)malloc(sizeof(double));
					mat_data.atTS_flow_stress = (double *)malloc(sizeof(double));
					mat_data.TSfor_flow_stress[0] = -1.0;
					mat_data.atTS_flow_stress[0] =
							(double)config_setting_get_float(setting);


				}
			}

			// THERMAL_COND.
			if (NULL != (setting = config_setting_get_member(material, "thermal_conductivity"))) {

				if (config_setting_is_list(setting) == CONFIG_TRUE)
				{
					setting_length = config_setting_length(setting);
					mat_data.Tfor_thermal_conductivity =
							(double *) malloc(sizeof(double) * setting_length);
					mat_data.atT_thermal_conductivity =
							(double *) malloc(sizeof(double) * setting_length);
					mat_data.thermal_conductivity_num = setting_length;
					for (j = 0; j < setting_length; ++j)
					{
						mat_data.Tfor_thermal_conductivity[j] =
								(double)config_setting_get_float_elem(config_setting_get_elem(setting, j), 0);
						mat_data.atT_thermal_conductivity[j]  =
								(double)config_setting_get_float_elem(config_setting_get_elem(setting, j), 1);
					}

//              // for temperature dependance reference temperature must be > 0
//              if(pdv_heat_problem.ctrl.ref_temperature <= 0.0){
//                fprintf(Interactive_output, "ref_temperature <= 0.0 for temperature dependent thermal_conductivity. Exiting \n");
//                config_destroy(&cfg);
//                exit(-1);
//                //return(EXIT_FAILURE);
//              }
					ret_val = UTE_REF_TEMP_MUST_BE_GEQ_ZERO;

					constant_thermal_conductivity = 0;

				}
				else
				{
					mat_data.thermal_conductivity_num = 1;
					mat_data.Tfor_thermal_conductivity =
							(double *) malloc(sizeof(double));
					mat_data.atT_thermal_conductivity =
							(double *) malloc(sizeof(double));
					mat_data.Tfor_thermal_conductivity[0] = -1.0;
					mat_data.atT_thermal_conductivity[0] =
							(double)config_setting_get_float(setting);
				}
			}

			// SPECIFIC HEAT
			if (NULL != (setting = config_setting_get_member(material, "specific_heat"))) {

				if (config_setting_is_list(setting) == CONFIG_TRUE)
				{
					setting_length = config_setting_length(setting);
					mat_data.Tfor_specific_heat =
							(double *) malloc(sizeof(double) * setting_length);
					mat_data.atT_specific_heat =
							(double *) malloc(sizeof(double) * setting_length);
					mat_data.specific_heat_num = setting_length;
					for (j = 0; j < setting_length; ++j)
					{
						mat_data.Tfor_specific_heat[j] =
								(double)config_setting_get_float_elem(config_setting_get_elem(setting, j), 0);
						mat_data.atT_specific_heat[j]  =
								(double)config_setting_get_float_elem(config_setting_get_elem(setting, j), 1);
					}

					// for temperature dependance reference temperature must be > 0
//              if(pdv_heat_problem.ctrl.ref_temperature <= 0.0){
//                fprintf(Interactive_output, "ref_temperature <= 0.0 for temperature dependent specific heat. Exiting \n");
//                config_destroy(&cfg);
//                exit(-1);
//                //return(EXIT_FAILURE);
//              }
					ret_val = UTE_REF_TEMP_MUST_BE_GEQ_ZERO;

					constant_specific_heat = 0;

				}
				else
				{
					mat_data.specific_heat_num = 1;
					mat_data.Tfor_specific_heat =
							(double *) malloc(sizeof(double));
					mat_data.atT_specific_heat =
							(double *) malloc(sizeof(double));
					mat_data.Tfor_specific_heat[0] = -1.0;
					mat_data.atT_specific_heat[0]  =
							(double)config_setting_get_float(setting);
				}
			}

#ifdef PHASE_TRANSFORMATION
			// MARTENSITE_START
			if (NULL != (setting = config_setting_get_member(material, "T_Ms"))) {

				if (config_setting_is_list(setting) == CONFIG_TRUE)
				{
					setting_length = config_setting_length(setting);
					mat_data.tfor_T_Ms =
							(double *)malloc(sizeof(double) * setting_length);
					mat_data.att_T_Ms =
							(double *)malloc(sizeof(double) * setting_length);
					mat_data.T_Ms_num = setting_length;
					for (j = 0; j < setting_length; ++j)
					{
						mat_data.tfor_T_Ms[j] =
								(double)config_setting_get_float_elem(config_setting_get_elem(setting, j), 0);
						mat_data.att_T_Ms[j] =
								(double)config_setting_get_float_elem(config_setting_get_elem(setting, j), 1);
					}

					//              // for temperature dependance reference temperature must be > 0
					//              //TODO: move to check problem module
					//              if(pdv_heat_problem.ctrl.ref_temperature <= 0.0){
					//                fprintf(Interactive_output, "ref_temperature <= 0.0 for temperature dependent visocisty. Exiting \n");
					//                config_destroy(&cfg);
					//                exit(-1);
					//                //return(EXIT_FAILURE);
					//              }
					ret_val = UTE_REF_TEMP_MUST_BE_GEQ_ZERO;

					constant_T_Ms = 0;

				}
				else if (NULL != config_setting_get_string(setting))
				{
#ifdef ANALYTIC_EXPRESSIONS
					mat_data.T_Ms_num = 1;
					mat_data.tfor_T_Ms = (double *)malloc(sizeof(double));
					mat_data.expression_of_T_Ms_from_t = (void *) new exprtk::expression<double>();

					std::string expression_string(config_setting_get_string(setting));
					exprtk::symbol_table<double>	symbol_table;
					exprtk::expression<double>	&	expression = *((exprtk::expression<double> *)mat_data.expression_of_T_Ms_from_t);

					symbol_table.add_variable("t", mat_data.tfor_T_Ms[0]);	// current simulation time
					symbol_table.add_constants();
					expression.register_symbol_table(symbol_table);

					exprtk::parser<double>			parser;
					parser.compile(expression_string, expression);
#endif
				}
				else
				{
					mat_data.T_Ms_num = 1;
					mat_data.tfor_T_Ms = (double *)malloc(sizeof(double));
					mat_data.att_T_Ms = (double *)malloc(sizeof(double));
					mat_data.tfor_T_Ms[0] = -1.0;
					mat_data.att_T_Ms[0] =
							(double)config_setting_get_float(setting);


				}
			}

			// MARTENSITE_FINISH
			if (NULL != (setting = config_setting_get_member(material, "T_Mf"))) {

				if (config_setting_is_list(setting) == CONFIG_TRUE)
				{
					setting_length = config_setting_length(setting);
					mat_data.tfor_T_Mf =
							(double *)malloc(sizeof(double) * setting_length);
					mat_data.att_T_Mf =
							(double *)malloc(sizeof(double) * setting_length);
					mat_data.T_Mf_num = setting_length;
					for (j = 0; j < setting_length; ++j)
					{
						mat_data.tfor_T_Mf[j] =
								(double)config_setting_get_float_elem(config_setting_get_elem(setting, j), 0);
						mat_data.att_T_Mf[j] =
								(double)config_setting_get_float_elem(config_setting_get_elem(setting, j), 1);
					}

					//              // for temperature dependance reference temperature must be > 0
					//              //TODO: move to check problem module
					//              if(pdv_heat_problem.ctrl.ref_temperature <= 0.0){
					//                fprintf(Interactive_output, "ref_temperature <= 0.0 for temperature dependent visocisty. Exiting \n");
					//                config_destroy(&cfg);
					//                exit(-1);
					//                //return(EXIT_FAILURE);
					//              }
					ret_val = UTE_REF_TEMP_MUST_BE_GEQ_ZERO;

					constant_T_Mf = 0;

				}
				else if (NULL != config_setting_get_string(setting))
				{
#ifdef ANALYTIC_EXPRESSIONS
					mat_data.T_Mf_num = 1;
					mat_data.tfor_T_Mf = (double *)malloc(sizeof(double));
					mat_data.expression_of_T_Mf_from_t = (void *) new exprtk::expression<double>();

					std::string expression_string(config_setting_get_string(setting));
					exprtk::symbol_table<double>	symbol_table;
					exprtk::expression<double>	&	expression = *((exprtk::expression<double> *)mat_data.expression_of_T_Mf_from_t);

					symbol_table.add_variable("t", mat_data.tfor_T_Mf[0]);	// current simulation time
					symbol_table.add_constants();
					expression.register_symbol_table(symbol_table);

					exprtk::parser<double>			parser;
					parser.compile(expression_string, expression);
#endif
				}
				else
				{
					mat_data.T_Mf_num = 1;
					mat_data.tfor_T_Mf = (double *)malloc(sizeof(double));
					mat_data.att_T_Mf = (double *)malloc(sizeof(double));
					mat_data.tfor_T_Mf[0] = -1.0;
					mat_data.att_T_Mf[0] =
							(double)config_setting_get_float(setting);


				}
			}

			// BAINITE_START
			if (NULL != (setting = config_setting_get_member(material, "T_Bs"))) {

				if (config_setting_is_list(setting) == CONFIG_TRUE)
				{
					setting_length = config_setting_length(setting);
					mat_data.tfor_T_Bs =
							(double *)malloc(sizeof(double) * setting_length);
					mat_data.att_T_Bs =
							(double *)malloc(sizeof(double) * setting_length);
					mat_data.T_Bs_num = setting_length;
					for (j = 0; j < setting_length; ++j)
					{
						mat_data.tfor_T_Bs[j] =
								(double)config_setting_get_float_elem(config_setting_get_elem(setting, j), 0);
						mat_data.att_T_Bs[j] =
								(double)config_setting_get_float_elem(config_setting_get_elem(setting, j), 1);
					}

					//              // for temperature dependance reference temperature must be > 0
					//              //TODO: move to check problem module
					//              if(pdv_heat_problem.ctrl.ref_temperature <= 0.0){
					//                fprintf(Interactive_output, "ref_temperature <= 0.0 for temperature dependent visocisty. Exiting \n");
					//                config_destroy(&cfg);
					//                exit(-1);
					//                //return(EXIT_FAILURE);
					//              }
					ret_val = UTE_REF_TEMP_MUST_BE_GEQ_ZERO;

					constant_T_Bs = 0;

				}
				else if (NULL != config_setting_get_string(setting))
				{
#ifdef ANALYTIC_EXPRESSIONS
					mat_data.T_Bs_num = 1;
					mat_data.tfor_T_Bs = (double *)malloc(sizeof(double));
					mat_data.expression_of_T_Bs_from_t = (void *) new exprtk::expression<double>();

					std::string expression_string(config_setting_get_string(setting));
					exprtk::symbol_table<double>	symbol_table;
					exprtk::expression<double>	&	expression = *((exprtk::expression<double> *)mat_data.expression_of_T_Bs_from_t);

					symbol_table.add_variable("t", mat_data.tfor_T_Bs[0]);	// current simulation time
					symbol_table.add_constants();
					expression.register_symbol_table(symbol_table);

					exprtk::parser<double>			parser;
					parser.compile(expression_string, expression);
#endif
				}
				else
				{
					mat_data.T_Bs_num = 1;
					mat_data.tfor_T_Bs = (double *)malloc(sizeof(double));
					mat_data.att_T_Bs = (double *)malloc(sizeof(double));
					mat_data.tfor_T_Bs[0] = -1.0;
					mat_data.att_T_Bs[0] =
							(double)config_setting_get_float(setting);


				}
			}

			// BAINITE_FINISH
			if (NULL != (setting = config_setting_get_member(material, "T_Bf"))) {

				if (config_setting_is_list(setting) == CONFIG_TRUE)
				{
					setting_length = config_setting_length(setting);
					mat_data.tfor_T_Bf =
							(double *)malloc(sizeof(double) * setting_length);
					mat_data.att_T_Bf =
							(double *)malloc(sizeof(double) * setting_length);
					mat_data.T_Bf_num = setting_length;
					for (j = 0; j < setting_length; ++j)
					{
						mat_data.tfor_T_Bf[j] =
								(double)config_setting_get_float_elem(config_setting_get_elem(setting, j), 0);
						mat_data.att_T_Bf[j] =
								(double)config_setting_get_float_elem(config_setting_get_elem(setting, j), 1);
					}

					//              // for temperature dependance reference temperature must be > 0
					//              //TODO: move to check problem module
					//              if(pdv_heat_problem.ctrl.ref_temperature <= 0.0){
					//                fprintf(Interactive_output, "ref_temperature <= 0.0 for temperature dependent visocisty. Exiting \n");
					//                config_destroy(&cfg);
					//                exit(-1);
					//                //return(EXIT_FAILURE);
					//              }
					ret_val = UTE_REF_TEMP_MUST_BE_GEQ_ZERO;

					constant_T_Bf = 0;

				}
				else if (NULL != config_setting_get_string(setting))
				{
#ifdef ANALYTIC_EXPRESSIONS
					mat_data.T_Bf_num = 1;
					mat_data.tfor_T_Bf = (double *)malloc(sizeof(double));
					mat_data.expression_of_T_Bf_from_t = (void *) new exprtk::expression<double>();

					std::string expression_string(config_setting_get_string(setting));
					exprtk::symbol_table<double>	symbol_table;
					exprtk::expression<double>	&	expression = *((exprtk::expression<double> *)mat_data.expression_of_T_Bf_from_t);

					symbol_table.add_variable("t", mat_data.tfor_T_Bf[0]);	// current simulation time
					symbol_table.add_constants();
					expression.register_symbol_table(symbol_table);

					exprtk::parser<double>			parser;
					parser.compile(expression_string, expression);
#endif
				}
				else
				{
					mat_data.T_Bf_num = 1;
					mat_data.tfor_T_Bf = (double *)malloc(sizeof(double));
					mat_data.att_T_Bf = (double *)malloc(sizeof(double));
					mat_data.tfor_T_Bf[0] = -1.0;
					mat_data.att_T_Bf[0] =
							(double)config_setting_get_float(setting);


				}
			}

			// Heat of the Austenite to Bainite transformation
			if (NULL != (setting = config_setting_get_member(material, "H_A_B"))) {

				if (config_setting_is_list(setting) == CONFIG_TRUE)
				{
					setting_length = config_setting_length(setting);
					mat_data.Tfor_H_A_B =
							(double *)malloc(sizeof(double) * setting_length);
					mat_data.atT_H_A_B =
							(double *)malloc(sizeof(double) * setting_length);
					mat_data.H_A_B_num = setting_length;
					for (j = 0; j < setting_length; ++j)
					{
						mat_data.Tfor_H_A_B[j] =
								(double)config_setting_get_float_elem(config_setting_get_elem(setting, j), 0);
						mat_data.atT_H_A_B[j] =
								(double)config_setting_get_float_elem(config_setting_get_elem(setting, j), 1);
					}

					//              // for temperature dependance reference temperature must be > 0
					//              //TODO: move to check problem module
					//              if(pdv_heat_problem.ctrl.ref_temperature <= 0.0){
					//                fprintf(Interactive_output, "ref_temperature <= 0.0 for temperature dependent visocisty. Exiting \n");
					//                config_destroy(&cfg);
					//                exit(-1);
					//                //return(EXIT_FAILURE);
					//              }
					ret_val = UTE_REF_TEMP_MUST_BE_GEQ_ZERO;

					constant_H_A_B = 0;

				}
				else if (NULL != config_setting_get_string(setting))
				{
#ifdef ANALYTIC_EXPRESSIONS
					mat_data.H_A_B_num = 1;
					mat_data.Tfor_H_A_B = (double *)malloc(sizeof(double));
					mat_data.expression_of_H_A_B_from_T = (void *) new exprtk::expression<double>();

					std::string expression_string(config_setting_get_string(setting));
					exprtk::symbol_table<double>	symbol_table;
					exprtk::expression<double>	&	expression = *((exprtk::expression<double> *)mat_data.expression_of_H_A_B_from_T);

					symbol_table.add_variable("T", mat_data.Tfor_H_A_B[0]);	// current simulation time
					symbol_table.add_constants();
					expression.register_symbol_table(symbol_table);

					exprtk::parser<double>			parser;
					parser.compile(expression_string, expression);
#endif
				}
				else
				{
					mat_data.H_A_B_num = 1;
					mat_data.Tfor_H_A_B = (double *)malloc(sizeof(double));
					mat_data.atT_H_A_B = (double *)malloc(sizeof(double));
					mat_data.Tfor_H_A_B[0] = -1.0;
					mat_data.atT_H_A_B[0] =
							(double)config_setting_get_float(setting);


				}
			}

			// Heat of the Austenite to Martensite transformation
			if (NULL != (setting = config_setting_get_member(material, "H_A_M"))) {

				if (config_setting_is_list(setting) == CONFIG_TRUE)
				{
					setting_length = config_setting_length(setting);
					mat_data.Tfor_H_A_M =
							(double *)malloc(sizeof(double) * setting_length);
					mat_data.atT_H_A_M =
							(double *)malloc(sizeof(double) * setting_length);
					mat_data.H_A_M_num = setting_length;
					for (j = 0; j < setting_length; ++j)
					{
						mat_data.Tfor_H_A_M[j] =
								(double)config_setting_get_float_elem(config_setting_get_elem(setting, j), 0);
						mat_data.atT_H_A_M[j] =
								(double)config_setting_get_float_elem(config_setting_get_elem(setting, j), 1);
					}

					//              // for temperature dependance reference temperature must be > 0
					//              //TODO: move to check problem module
					//              if(pdv_heat_problem.ctrl.ref_temperature <= 0.0){
					//                fprintf(Interactive_output, "ref_temperature <= 0.0 for temperature dependent visocisty. Exiting \n");
					//                config_destroy(&cfg);
					//                exit(-1);
					//                //return(EXIT_FAILURE);
					//              }
					ret_val = UTE_REF_TEMP_MUST_BE_GEQ_ZERO;

					constant_H_A_M = 0;

				}
				else if (NULL != config_setting_get_string(setting))
				{
#ifdef ANALYTIC_EXPRESSIONS
					mat_data.H_A_M_num = 1;
					mat_data.Tfor_H_A_M = (double *)malloc(sizeof(double));
					mat_data.expression_of_H_A_M_from_T = (void *) new exprtk::expression<double>();

					std::string expression_string(config_setting_get_string(setting));
					exprtk::symbol_table<double>	symbol_table;
					exprtk::expression<double>	&	expression = *((exprtk::expression<double> *)mat_data.expression_of_H_A_M_from_T);

					symbol_table.add_variable("T", mat_data.Tfor_H_A_M[0]);	// current simulation time
					symbol_table.add_constants();
					expression.register_symbol_table(symbol_table);

					exprtk::parser<double>			parser;
					parser.compile(expression_string, expression);
#endif
				}
				else
				{
					mat_data.H_A_M_num = 1;
					mat_data.Tfor_H_A_M = (double *)malloc(sizeof(double));
					mat_data.atT_H_A_M = (double *)malloc(sizeof(double));
					mat_data.Tfor_H_A_M[0] = -1.0;
					mat_data.atT_H_A_M[0] =
							(double)config_setting_get_float(setting);


				}
			}

			// Heat of the Austenite to Ferrite transformation
			if (NULL != (setting = config_setting_get_member(material, "H_A_F"))) {

				if (config_setting_is_list(setting) == CONFIG_TRUE)
				{
					setting_length = config_setting_length(setting);
					mat_data.Tfor_H_A_F =
							(double *)malloc(sizeof(double) * setting_length);
					mat_data.atT_H_A_F =
							(double *)malloc(sizeof(double) * setting_length);
					mat_data.H_A_F_num = setting_length;
					for (j = 0; j < setting_length; ++j)
					{
						mat_data.Tfor_H_A_F[j] =
								(double)config_setting_get_float_elem(config_setting_get_elem(setting, j), 0);
						mat_data.atT_H_A_F[j] =
								(double)config_setting_get_float_elem(config_setting_get_elem(setting, j), 1);
					}

					//              // for temperature dependance reference temperature must be > 0
					//              //TODO: move to check problem module
					//              if(pdv_heat_problem.ctrl.ref_temperature <= 0.0){
					//                fprintf(Interactive_output, "ref_temperature <= 0.0 for temperature dependent visocisty. Exiting \n");
					//                config_destroy(&cfg);
					//                exit(-1);
					//                //return(EXIT_FAILURE);
					//              }
					ret_val = UTE_REF_TEMP_MUST_BE_GEQ_ZERO;

					constant_H_A_F = 0;

				}
				else if (NULL != config_setting_get_string(setting))
				{
#ifdef ANALYTIC_EXPRESSIONS
					mat_data.H_A_F_num = 1;
					mat_data.Tfor_H_A_F = (double *)malloc(sizeof(double));
					mat_data.expression_of_H_A_F_from_T = (void *) new exprtk::expression<double>();

					std::string expression_string(config_setting_get_string(setting));
					exprtk::symbol_table<double>	symbol_table;
					exprtk::expression<double>	&	expression = *((exprtk::expression<double> *)mat_data.expression_of_H_A_F_from_T);

					symbol_table.add_variable("T", mat_data.Tfor_H_A_F[0]);	// current simulation time
					symbol_table.add_constants();
					expression.register_symbol_table(symbol_table);

					exprtk::parser<double>			parser;
					parser.compile(expression_string, expression);
#endif
				}
				else
				{
					mat_data.H_A_F_num = 1;
					mat_data.Tfor_H_A_F = (double *)malloc(sizeof(double));
					mat_data.atT_H_A_F = (double *)malloc(sizeof(double));
					mat_data.Tfor_H_A_F[0] = -1.0;
					mat_data.atT_H_A_F[0] =
							(double)config_setting_get_float(setting);


				}
			}

			// Heat of the Austenite to Pearlite transformation
			if (NULL != (setting = config_setting_get_member(material, "H_A_P"))) {

				if (config_setting_is_list(setting) == CONFIG_TRUE)
				{
					setting_length = config_setting_length(setting);
					mat_data.Tfor_H_A_P =
							(double *)malloc(sizeof(double) * setting_length);
					mat_data.atT_H_A_P =
							(double *)malloc(sizeof(double) * setting_length);
					mat_data.H_A_P_num = setting_length;
					for (j = 0; j < setting_length; ++j)
					{
						mat_data.Tfor_H_A_P[j] =
								(double)config_setting_get_float_elem(config_setting_get_elem(setting, j), 0);
						mat_data.atT_H_A_P[j] =
								(double)config_setting_get_float_elem(config_setting_get_elem(setting, j), 1);
					}

					//              // for temperature dependance reference temperature must be > 0
					//              //TODO: move to check problem module
					//              if(pdv_heat_problem.ctrl.ref_temperature <= 0.0){
					//                fprintf(Interactive_output, "ref_temperature <= 0.0 for temperature dependent visocisty. Exiting \n");
					//                config_destroy(&cfg);
					//                exit(-1);
					//                //return(EXIT_FAILURE);
					//              }
					ret_val = UTE_REF_TEMP_MUST_BE_GEQ_ZERO;

					constant_H_A_P = 0;

				}
				else if (NULL != config_setting_get_string(setting))
				{
#ifdef ANALYTIC_EXPRESSIONS
					mat_data.H_A_P_num = 1;
					mat_data.Tfor_H_A_P = (double *)malloc(sizeof(double));
					mat_data.expression_of_H_A_P_from_T = (void *) new exprtk::expression<double>();

					std::string expression_string(config_setting_get_string(setting));
					exprtk::symbol_table<double>	symbol_table;
					exprtk::expression<double>	&	expression = *((exprtk::expression<double> *)mat_data.expression_of_H_A_P_from_T);

					symbol_table.add_variable("T", mat_data.Tfor_H_A_P[0]);	// current simulation time
					symbol_table.add_constants();
					expression.register_symbol_table(symbol_table);

					exprtk::parser<double>			parser;
					parser.compile(expression_string, expression);
#endif
				}
				else
				{
					mat_data.H_A_P_num = 1;
					mat_data.Tfor_H_A_P = (double *)malloc(sizeof(double));
					mat_data.atT_H_A_P = (double *)malloc(sizeof(double));
					mat_data.Tfor_H_A_P[0] = -1.0;
					mat_data.atT_H_A_P[0] =
							(double)config_setting_get_float(setting);


				}
			}

			// final volume fraction of MARTENSITE after cooling
			if (NULL != (setting = config_setting_get_member(material, "M_Vc"))) {

				if (config_setting_is_list(setting) == CONFIG_TRUE) {
					setting_length = config_setting_length(setting);
					mat_data.Vcfor_M_Vc = (double *) malloc(sizeof(double) * setting_length);
					mat_data.atVc_M_Vc  = (double *) malloc(sizeof(double) * setting_length);
					mat_data.M_Vc_num = setting_length;
					for (j = 0; j < setting_length; ++j) {
						mat_data.Vcfor_M_Vc[j] = (double)config_setting_get_float_elem(config_setting_get_elem(setting, j), 0);
						mat_data.atVc_M_Vc[j]  = (double)config_setting_get_float_elem(config_setting_get_elem(setting, j), 1);
					}

					ret_val = UTE_REF_TEMP_MUST_BE_GEQ_ZERO;

					constant_M_Vc = 0;

				} else {
					mat_data.M_Vc_num = 1;
					mat_data.Vcfor_M_Vc = (double *) malloc(sizeof(double));
					mat_data.atVc_M_Vc = (double *) malloc(sizeof(double));
					mat_data.Vcfor_M_Vc[0] = -1.0;
					mat_data.atVc_M_Vc[0] = (double)config_setting_get_float(setting);
				}
			}

			// final volume fraction of BAINITE after cooling
			if (NULL != (setting = config_setting_get_member(material, "B_Vc"))) {

				if (config_setting_is_list(setting) == CONFIG_TRUE) {
					setting_length = config_setting_length(setting);
					mat_data.Vcfor_B_Vc = (double *) malloc(sizeof(double) * setting_length);
					mat_data.atVc_B_Vc  = (double *) malloc(sizeof(double) * setting_length);
					mat_data.B_Vc_num = setting_length;
					for (j = 0; j < setting_length; ++j) {
						mat_data.Vcfor_B_Vc[j] = (double)config_setting_get_float_elem(config_setting_get_elem(setting, j), 0);
						mat_data.atVc_B_Vc[j]  = (double)config_setting_get_float_elem(config_setting_get_elem(setting, j), 1);
					}

					ret_val = UTE_REF_TEMP_MUST_BE_GEQ_ZERO;

					constant_B_Vc = 0;

				} else {
					mat_data.B_Vc_num = 1;
					mat_data.Vcfor_B_Vc = (double *) malloc(sizeof(double));
					mat_data.atVc_B_Vc = (double *) malloc(sizeof(double));
					mat_data.Vcfor_B_Vc[0] = -1.0;
					mat_data.atVc_B_Vc[0] = (double)config_setting_get_float(setting);
				}
			}

			// final volume fraction of FERRITE after cooling
			if (NULL != (setting = config_setting_get_member(material, "F_Vc"))) {

				if (config_setting_is_list(setting) == CONFIG_TRUE) {
					setting_length = config_setting_length(setting);
					mat_data.Vcfor_F_Vc = (double *) malloc(sizeof(double) * setting_length);
					mat_data.atVc_F_Vc  = (double *) malloc(sizeof(double) * setting_length);
					mat_data.F_Vc_num = setting_length;
					for (j = 0; j < setting_length; ++j) {
						mat_data.Vcfor_F_Vc[j] = (double)config_setting_get_float_elem(config_setting_get_elem(setting, j), 0);
						mat_data.atVc_F_Vc[j]  = (double)config_setting_get_float_elem(config_setting_get_elem(setting, j), 1);
					}

					ret_val = UTE_REF_TEMP_MUST_BE_GEQ_ZERO;

					constant_F_Vc = 0;

				} else {
					mat_data.F_Vc_num = 1;
					mat_data.Vcfor_F_Vc = (double *) malloc(sizeof(double));
					mat_data.atVc_F_Vc = (double *) malloc(sizeof(double));
					mat_data.Vcfor_F_Vc[0] = -1.0;
					mat_data.atVc_F_Vc[0] = (double)config_setting_get_float(setting);
				}
			}

			// final volume fraction of PERLITE after cooling
			if (NULL != (setting = config_setting_get_member(material, "P_Vc"))) {

				if (config_setting_is_list(setting) == CONFIG_TRUE) {
					setting_length = config_setting_length(setting);
					mat_data.Vcfor_P_Vc = (double *) malloc(sizeof(double) * setting_length);
					mat_data.atVc_P_Vc  = (double *) malloc(sizeof(double) * setting_length);
					mat_data.P_Vc_num = setting_length;
					for (j = 0; j < setting_length; ++j) {
						mat_data.Vcfor_P_Vc[j] = (double)config_setting_get_float_elem(config_setting_get_elem(setting, j), 0);
						mat_data.atVc_P_Vc[j]  = (double)config_setting_get_float_elem(config_setting_get_elem(setting, j), 1);
					}

					ret_val = UTE_REF_TEMP_MUST_BE_GEQ_ZERO;

					constant_P_Vc = 0;

				} else {
					mat_data.P_Vc_num = 1;
					mat_data.Vcfor_P_Vc = (double *) malloc(sizeof(double));
					mat_data.atVc_P_Vc = (double *) malloc(sizeof(double));
					mat_data.Vcfor_P_Vc[0] = -1.0;
					mat_data.atVc_P_Vc[0] = (double)config_setting_get_float(setting);
				}
			}

#endif
			/////////////////////////// NS_SUPG/////////////////////////
			// DYNAMIC_VISCOSITTY

			if (NULL != (setting = config_setting_get_member(material, "dynamic_viscosity"))) {

				if (config_setting_is_list(setting) == CONFIG_TRUE)
				{
					setting_length = config_setting_length(setting);
					//printf("\n\tAS: dynamic viscosity - setting_length =%d\n", setting_length);
					mat_data.Tfor_dynamic_viscosity =
							(double *) malloc(sizeof(double) * setting_length);
					mat_data.atT_dynamic_viscosity =
							(double *) malloc(sizeof(double) * setting_length);
					mat_data.dynamic_viscosity_num = setting_length;
					for (j = 0; j < setting_length; ++j)
					{
						mat_data.Tfor_dynamic_viscosity[j] =
								(double)config_setting_get_float_elem(
										config_setting_get_elem(setting, j), 0);
						mat_data.atT_dynamic_viscosity[j]  =
								(double)config_setting_get_float_elem(
										config_setting_get_elem(setting, j), 1);
					}

					// for temperature dependance reference temperature must be > 0
					//if(pdv_ns_supg_problem.ctrl.ref_temperature <= 0.0){
//                fprintf(Interactive_output, "ref_temperature <= 0.0 for temperature dependent visocisty. Exiting \n");
//                config_destroy(&cfg);
//                exit(-1);
//                //return(EXIT_FAILURE);
					ret_val = UTE_REF_TEMP_MUST_BE_GEQ_ZERO;


					constant_viscosity = 0;

				}
				else
				{

					mat_data.dynamic_viscosity_num = 1;
					mat_data.Tfor_dynamic_viscosity =
							(double *) malloc(sizeof(double));
					mat_data.atT_dynamic_viscosity =
							(double *) malloc(sizeof(double));
					mat_data.Tfor_dynamic_viscosity[0] = -1.0;
					mat_data.atT_dynamic_viscosity[0]  =
							(double)config_setting_get_float(setting);

				}

			}


		}

		// REGISTERING deleting function
		atexit(utr_mat_clear_all);

		config_destroy(&cfg);

		if (constant_thermal_conductivity == 1 && constant_density == 1
				&& constant_specific_heat == 1 && no_other_parameters == 1 && count == 1) {

			if (ret_val == UTE_SUCCESS) {
				ret_val = UTE_NOT_NEED_MATERIAL_DATABASE;
				mf_log_warn("Materials database not needed!");
			}
		}

	}

	mf_log_info("Materials read! ");

	return (ret_val);
}

// INTERFACE FUNCTIONS
void utr_mat_init_query_result(utt_material_query_result * result) {
	if (result != NULL) {
		result->is_fluid = FLUID;
		result->dynamic_viscosity = UTC_MAT_QUERY_RESULT_NOT_NEEDED;
		result->density = UTC_MAT_QUERY_RESULT_NOT_NEEDED;
		result->thermal_conductivity = UTC_MAT_QUERY_RESULT_NOT_NEEDED;
		result->specific_heat = UTC_MAT_QUERY_RESULT_NOT_NEEDED;
		result->enthalpy = UTC_MAT_QUERY_RESULT_NOT_NEEDED;
		result->thermal_expansion_coefficient = UTC_MAT_QUERY_RESULT_NOT_NEEDED;
		result->electrical_resistivity = UTC_MAT_QUERY_RESULT_NOT_NEEDED;
		result->VOF = UTC_MAT_QUERY_RESULT_NOT_NEEDED;
		result->dg_dT = UTC_MAT_QUERY_RESULT_NOT_NEEDED;
		result->temp_solidus = UTC_MAT_QUERY_RESULT_NOT_NEEDED;
		result->temp_liquidus = UTC_MAT_QUERY_RESULT_NOT_NEEDED;
		result->temp_vaporization = UTC_MAT_QUERY_RESULT_NOT_NEEDED;
		result->latent_heat_of_fusion = UTC_MAT_QUERY_RESULT_NOT_NEEDED;
		result->latent_heat_of_vaporization = UTC_MAT_QUERY_RESULT_NOT_NEEDED;
		result->flow_stress = UTC_MAT_QUERY_RESULT_NOT_NEEDED;
#ifdef PHASE_TRANSFORMATION
		result->T_Fs = UTC_MAT_QUERY_RESULT_NOT_NEEDED;
		result->T_Ff = UTC_MAT_QUERY_RESULT_NOT_NEEDED;
		result->T_Ps = UTC_MAT_QUERY_RESULT_NOT_NEEDED;
		result->T_Pf = UTC_MAT_QUERY_RESULT_NOT_NEEDED;
		result->T_Bs = UTC_MAT_QUERY_RESULT_NOT_NEEDED;
		result->T_Bf = UTC_MAT_QUERY_RESULT_NOT_NEEDED;
		result->T_Ms = UTC_MAT_QUERY_RESULT_NOT_NEEDED;
		result->T_Mf = UTC_MAT_QUERY_RESULT_NOT_NEEDED;

		result->M_Vc = UTC_MAT_QUERY_RESULT_NOT_NEEDED;
		result->B_Vc = UTC_MAT_QUERY_RESULT_NOT_NEEDED;
		result->F_Vc = UTC_MAT_QUERY_RESULT_NOT_NEEDED;
		result->P_Vc = UTC_MAT_QUERY_RESULT_NOT_NEEDED;

		result->H_A_F = UTC_MAT_QUERY_RESULT_NOT_NEEDED;
		result->H_A_P = UTC_MAT_QUERY_RESULT_NOT_NEEDED;
		result->H_A_B = UTC_MAT_QUERY_RESULT_NOT_NEEDED;
		result->H_A_M = UTC_MAT_QUERY_RESULT_NOT_NEEDED;
#endif
	}
}


ute_mat_read_result utr_mat_read(const char * Work_dir,
		const char * Filename,
		FILE * Interactive_output)
{
	// If file already read ignore it.
	if ( utv_read_material_files.find(Filename) != utv_read_material_files.end()) {
		return UTE_SUCCESS;
	}
	else if (utv_read_material_files.empty()) {
		// read main database material file
		utr_mat_read_material_file(Work_dir, utc_mat_database_filename, Interactive_output);
	}

	return utr_mat_read_material_file(Work_dir, Filename, Interactive_output);

//    utr_mat_read_material_file(Work_dir,Filename,Interactive_output);
//    return utr_mat_read_material_file(Work_dir,utc_mat_database_filename,Interactive_output);


}




int utr_material_query(const utt_material_query_params * Params,
		utt_material_query_result * Result)
{
	int i;
	//int material_idx = -1;
	const utt_material_data * material = utr_mat_get_material(Params->group_idx);



	double v1, v2, t1, t2;
	double a, b;

	strcpy(Result->name, material->name);
	//result->queryparams = params;

	Result->is_fluid = material->is_fluid;

	if (Result->dynamic_viscosity == UTC_MAT_QUERY_RESULT_REQUIRED) {
		//viscosity
		if (material->dynamic_viscosity_num == 1)
			Result->dynamic_viscosity = material->atT_dynamic_viscosity[0];
		else {



			if (Params->temperature <= material->Tfor_dynamic_viscosity[0])
				Result->dynamic_viscosity = material->atT_dynamic_viscosity[0];
			else if (Params->temperature >=
					material->Tfor_dynamic_viscosity[material->dynamic_viscosity_num - 1])
				Result->dynamic_viscosity =
						material->atT_dynamic_viscosity[material->dynamic_viscosity_num - 1];
			else
				for (i = 0; i < material->dynamic_viscosity_num - 1; ++i) {
					if (Params->temperature <= material->Tfor_dynamic_viscosity[i + 1]) {
						t1 = material->Tfor_dynamic_viscosity[i];
						t2 = material->Tfor_dynamic_viscosity[i + 1];
						v1 = material->atT_dynamic_viscosity[i];
						v2 = material->atT_dynamic_viscosity[i + 1];
						a = (v1 - v2) / (t1 - t2);
						b = v1 - a * t1;
						Result->dynamic_viscosity = a * (Params->temperature) + b;
						break;
					}
				}
		}
	}

	//density
	if (Result->density == UTC_MAT_QUERY_RESULT_REQUIRED ) {
		if (material->density_num == 1) {
			if (material->expression_of_density_from_T != NULL) {
#ifdef ANALYTIC_EXPRESSIONS
				material->Tfor_density[0] = Params->temperature;
				// calculating value using given formula
				Result->density = ((exprtk::expression<double> *)material->expression_of_density_from_T)->value();
#endif
			}
			else {
				Result->density = material->atT_density[0];
			}
		}
		else {
			if (Params->temperature <= material->Tfor_density[0])
				Result->density = material->atT_density[0];
			else if (Params->temperature >=
					material->Tfor_density[material->density_num - 1])
				Result->density = material->atT_density[material->density_num - 1];
			else
				for (i = 0; i < material->density_num - 1; ++i) {
					if (Params->temperature <= material->Tfor_density[i + 1]) {
						t1 = material->Tfor_density[i];
						t2 = material->Tfor_density[i + 1];
						v1 = material->atT_density[i];
						v2 = material->atT_density[i + 1];
						a = (v1 - v2) / (t1 - t2);
						b = v1 - a * t1;
						Result->density = a * (Params->temperature) + b;
						break;
					}
				}
		}
	}

	//flow_stress
	if (Result->flow_stress == UTC_MAT_QUERY_RESULT_REQUIRED ) {
		if (material->flow_stress_num == 1) {
			if (material->expression_of_flow_stress_from_S != NULL) {
#ifdef ANALYTIC_EXPRESSIONS
				material->TSfor_flow_stress[0] = Params->strain_intensity;
				// calculating value using given formula
				Result->flow_stress = ((exprtk::expression<double> *)material->expression_of_flow_stress_from_S)->value();
#endif
			}
			else {
				Result->flow_stress = material->atTS_flow_stress[0];
			}
		}
		else {
			if (Params->strain_intensity <= material->TSfor_flow_stress[0])
				Result->flow_stress = material->atTS_flow_stress[0];
			else if (Params->strain_intensity >=
					material->TSfor_flow_stress[material->flow_stress_num - 1])
				Result->flow_stress = material->atTS_flow_stress[material->flow_stress_num - 1];
			else
				for (i = 0; i < material->flow_stress_num - 1; ++i) {
					if (Params->strain_intensity <= material->TSfor_flow_stress[i + 1]) {
						t1 = material->TSfor_flow_stress[i];
						t2 = material->TSfor_flow_stress[i + 1];
						v1 = material->atTS_flow_stress[i];
						v2 = material->atTS_flow_stress[i + 1];
						a = (v1 - v2) / (t1 - t2);
						b = v1 - a * t1;
						Result->flow_stress = a * (Params->strain_intensity) + b;
						break;
					}
				}
		}
	}

#ifdef PHASE_TRANSFORMATION

	//T_Ms
	if (Result->T_Ms == UTC_MAT_QUERY_RESULT_REQUIRED ) {
		if (material->T_Ms_num == 1) {
			if (material->expression_of_T_Ms_from_t != NULL) {
#ifdef ANALYTIC_EXPRESSIONS
				material->tfor_T_Ms[0] = Params->time;		// current simulation time
				// calculating value using given formula
				Result->T_Ms = ((exprtk::expression<double> *)material->expression_of_T_Ms_from_t)->value();
#endif
			}
			else {
				Result->T_Ms = material->att_T_Ms[0];
			}
		}
		else {
			if (Params->time <= material->tfor_T_Ms[0])
				Result->T_Ms = material->att_T_Ms[0];
			else if (Params->time >=
					material->tfor_T_Ms[material->T_Ms_num - 1])
				Result->T_Ms = material->att_T_Ms[material->T_Ms_num - 1];
			else
				for (i = 0; i < material->T_Ms_num - 1; ++i) {
					if (Params->time <= material->tfor_T_Ms[i + 1]) {
						t1 = material->tfor_T_Ms[i];
						t2 = material->tfor_T_Ms[i + 1];
						v1 = material->att_T_Ms[i];
						v2 = material->att_T_Ms[i + 1];
						a = (v1 - v2) / (t1 - t2);
						b = v1 - a * t1;
						Result->T_Ms = a * (Params->time) + b;
						break;
					}
				}
		}
	}

	//T_Mf
	if (Result->T_Mf == UTC_MAT_QUERY_RESULT_REQUIRED ) {
		if (material->T_Mf_num == 1) {
			if (material->expression_of_T_Mf_from_t != NULL) {
#ifdef ANALYTIC_EXPRESSIONS
				material->tfor_T_Mf[0] = Params->time;		// current simulation time
				// calculating value using given formula
				Result->T_Mf = ((exprtk::expression<double> *)material->expression_of_T_Mf_from_t)->value();
#endif
			}
			else {
				Result->T_Mf = material->att_T_Mf[0];
			}
		}
		else {
			if (Params->time <= material->tfor_T_Mf[0])
				Result->T_Mf = material->att_T_Mf[0];
			else if (Params->time >=
					material->tfor_T_Mf[material->T_Mf_num - 1])
				Result->T_Mf = material->att_T_Mf[material->T_Mf_num - 1];
			else
				for (i = 0; i < material->T_Mf_num - 1; ++i) {
					if (Params->time <= material->tfor_T_Mf[i + 1]) {
						t1 = material->tfor_T_Mf[i];
						t2 = material->tfor_T_Mf[i + 1];
						v1 = material->att_T_Mf[i];
						v2 = material->att_T_Mf[i + 1];
						a = (v1 - v2) / (t1 - t2);
						b = v1 - a * t1;
						Result->T_Mf = a * (Params->time) + b;
						break;
					}
				}
		}
	}

	//T_Bs
	if (Result->T_Bs == UTC_MAT_QUERY_RESULT_REQUIRED ) {
		if (material->T_Bs_num == 1) {
			if (material->expression_of_T_Bs_from_t != NULL) {
#ifdef ANALYTIC_EXPRESSIONS
				material->tfor_T_Bs[0] = Params->time;		// current simulation time
				// calculating value using given formula
				Result->T_Bs = ((exprtk::expression<double> *)material->expression_of_T_Bs_from_t)->value();
#endif
			}
			else {
				Result->T_Bs = material->att_T_Bs[0];
			}
		}
		else {
			if (Params->time <= material->tfor_T_Bs[0])
				Result->T_Bs = material->att_T_Bs[0];
			else if (Params->time >=
					material->tfor_T_Bs[material->T_Bs_num - 1])
				Result->T_Bs = material->att_T_Bs[material->T_Bs_num - 1];
			else
				for (i = 0; i < material->T_Bs_num - 1; ++i) {
					if (Params->time <= material->tfor_T_Bs[i + 1]) {
						t1 = material->tfor_T_Bs[i];
						t2 = material->tfor_T_Bs[i + 1];
						v1 = material->att_T_Bs[i];
						v2 = material->att_T_Bs[i + 1];
						a = (v1 - v2) / (t1 - t2);
						b = v1 - a * t1;
						Result->T_Bs = a * (Params->time) + b;
						break;
					}
				}
		}
	}

	//T_Bf
	if (Result->T_Bf == UTC_MAT_QUERY_RESULT_REQUIRED ) {
		if (material->T_Bf_num == 1) {
			if (material->expression_of_T_Bf_from_t != NULL) {
#ifdef ANALYTIC_EXPRESSIONS
				material->tfor_T_Bf[0] = Params->time;		// current simulation time
				// calculating value using given formula
				Result->T_Bf = ((exprtk::expression<double> *)material->expression_of_T_Bf_from_t)->value();
#endif
			}
			else {
				Result->T_Bf = material->att_T_Bf[0];
			}
		}
		else {
			if (Params->time <= material->tfor_T_Bf[0])
				Result->T_Bf = material->att_T_Bf[0];
			else if (Params->time >=
					material->tfor_T_Bf[material->T_Bf_num - 1])
				Result->T_Bf = material->att_T_Bf[material->T_Bf_num - 1];
			else
				for (i = 0; i < material->T_Bf_num - 1; ++i) {
					if (Params->time <= material->tfor_T_Bf[i + 1]) {
						t1 = material->tfor_T_Bf[i];
						t2 = material->tfor_T_Bf[i + 1];
						v1 = material->att_T_Bf[i];
						v2 = material->att_T_Bf[i + 1];
						a = (v1 - v2) / (t1 - t2);
						b = v1 - a * t1;
						Result->T_Bf = a * (Params->time) + b;
						break;
					}
				}
		}
	}

	//M_Vc
	if (Result->M_Vc == UTC_MAT_QUERY_RESULT_REQUIRED) {
		if (material->M_Vc_num == 1)
			Result->M_Vc = material->atVc_M_Vc[0];
		else {
			if (Params->Vc <= material->Vcfor_M_Vc[0])
				Result->M_Vc = material->atVc_M_Vc[0];
			else if (Params->Vc >=
					material->Vcfor_M_Vc[material->M_Vc_num - 1])
				Result->M_Vc =
						material->atVc_M_Vc[material->M_Vc_num - 1];
			else
				for (i = 0; i < material->M_Vc_num - 1; ++i) {
					if (Params->Vc <= material->Vcfor_M_Vc[i + 1]) {
						t1 = material->Vcfor_M_Vc[i];
						t2 = material->Vcfor_M_Vc[i + 1];
						v1 = material->atVc_M_Vc[i];
						v2 = material->atVc_M_Vc[i + 1];
						a = (v1 - v2) / (t1 - t2);
						b = v1 - a * t1;
						Result->M_Vc = a * (Params->Vc) + b;
						break;
					}
				}
		}
	}

	//B_Vc
	if (Result->B_Vc == UTC_MAT_QUERY_RESULT_REQUIRED) {
		if (material->B_Vc_num == 1)
			Result->B_Vc = material->atVc_B_Vc[0];
		else {
			if (Params->Vc <= material->Vcfor_B_Vc[0])
				Result->B_Vc = material->atVc_B_Vc[0];
			else if (Params->Vc >=
					material->Vcfor_B_Vc[material->B_Vc_num - 1])
				Result->B_Vc =
						material->atVc_B_Vc[material->B_Vc_num - 1];
			else
				for (i = 0; i < material->B_Vc_num - 1; ++i) {
					if (Params->Vc <= material->Vcfor_B_Vc[i + 1]) {
						t1 = material->Vcfor_B_Vc[i];
						t2 = material->Vcfor_B_Vc[i + 1];
						v1 = material->atVc_B_Vc[i];
						v2 = material->atVc_B_Vc[i + 1];
						a = (v1 - v2) / (t1 - t2);
						b = v1 - a * t1;
						Result->B_Vc = a * (Params->Vc) + b;
						break;
					}
				}
		}
	}

	//F_Vc
	if (Result->F_Vc == UTC_MAT_QUERY_RESULT_REQUIRED) {
		if (material->F_Vc_num == 1)
			Result->F_Vc = material->atVc_F_Vc[0];
		else {
			if (Params->Vc <= material->Vcfor_F_Vc[0])
				Result->F_Vc = material->atVc_F_Vc[0];
			else if (Params->Vc >=
					material->Vcfor_F_Vc[material->F_Vc_num - 1])
				Result->F_Vc =
						material->atVc_F_Vc[material->F_Vc_num - 1];
			else
				for (i = 0; i < material->F_Vc_num - 1; ++i) {
					if (Params->Vc <= material->Vcfor_F_Vc[i + 1]) {
						t1 = material->Vcfor_F_Vc[i];
						t2 = material->Vcfor_F_Vc[i + 1];
						v1 = material->atVc_F_Vc[i];
						v2 = material->atVc_F_Vc[i + 1];
						a = (v1 - v2) / (t1 - t2);
						b = v1 - a * t1;
						Result->F_Vc = a * (Params->Vc) + b;
						break;
					}
				}
		}
	}

	//P_Vc
	if (Result->P_Vc == UTC_MAT_QUERY_RESULT_REQUIRED) {
		if (material->P_Vc_num == 1)
			Result->P_Vc = material->atVc_P_Vc[0];
		else {
			if (Params->Vc <= material->Vcfor_P_Vc[0])
				Result->P_Vc = material->atVc_P_Vc[0];
			else if (Params->Vc >=
					material->Vcfor_P_Vc[material->P_Vc_num - 1])
				Result->P_Vc =
						material->atVc_P_Vc[material->P_Vc_num - 1];
			else
				for (i = 0; i < material->P_Vc_num - 1; ++i) {
					if (Params->Vc <= material->Vcfor_P_Vc[i + 1]) {
						t1 = material->Vcfor_P_Vc[i];
						t2 = material->Vcfor_P_Vc[i + 1];
						v1 = material->atVc_P_Vc[i];
						v2 = material->atVc_P_Vc[i + 1];
						a = (v1 - v2) / (t1 - t2);
						b = v1 - a * t1;
						Result->P_Vc = a * (Params->Vc) + b;
						break;
					}
				}
		}
	}

	//H_A_M
	if (Result->H_A_M == UTC_MAT_QUERY_RESULT_REQUIRED ) {
		if (material->H_A_M_num == 1) {
			if (material->expression_of_H_A_M_from_T != NULL) {
#ifdef ANALYTIC_EXPRESSIONS
				material->Tfor_H_A_M[0] = Params->temperature;		// current temperature
				// calculating value using given formula
				Result->H_A_M = ((exprtk::expression<double> *)material->expression_of_H_A_M_from_T)->value();
#endif
			} else {
				Result->H_A_M = material->atT_H_A_M[0];
			}
		} else {
			if (Params->temperature <= material->Tfor_H_A_M[0])
				Result->H_A_M = material->atT_H_A_M[0];
			else if (Params->temperature >= material->Tfor_H_A_M[material->H_A_M_num - 1])
				Result->H_A_M = material->atT_H_A_M[material->H_A_M_num - 1];
			else
				for (i = 0; i < material->H_A_M_num - 1; ++i) {
					if (Params->temperature <= material->Tfor_H_A_M[i + 1]) {
						t1 = material->Tfor_H_A_M[i];
						t2 = material->Tfor_H_A_M[i + 1];
						v1 = material->atT_H_A_M[i];
						v2 = material->atT_H_A_M[i + 1];
						a = (v1 - v2) / (t1 - t2);
						b = v1 - a * t1;
						Result->H_A_M = a * (Params->temperature) + b;
						break;
					}
				}
		}
	}

	//H_A_B
	if (Result->H_A_B == UTC_MAT_QUERY_RESULT_REQUIRED ) {
		if (material->H_A_B_num == 1) {
			if (material->expression_of_H_A_B_from_T != NULL) {
#ifdef ANALYTIC_EXPRESSIONS
				material->Tfor_H_A_B[0] = Params->temperature;		// current temperature
				// calculating value using given formula
				Result->H_A_B = ((exprtk::expression<double> *)material->expression_of_H_A_B_from_T)->value();
#endif
			} else {
				Result->H_A_B = material->atT_H_A_B[0];
			}
		} else {
			if (Params->temperature <= material->Tfor_H_A_B[0])
				Result->H_A_B = material->atT_H_A_B[0];
			else if (Params->temperature >= material->Tfor_H_A_B[material->H_A_B_num - 1])
				Result->H_A_B = material->atT_H_A_B[material->H_A_B_num - 1];
			else
				for (i = 0; i < material->H_A_B_num - 1; ++i) {
					if (Params->temperature <= material->Tfor_H_A_B[i + 1]) {
						t1 = material->Tfor_H_A_B[i];
						t2 = material->Tfor_H_A_B[i + 1];
						v1 = material->atT_H_A_B[i];
						v2 = material->atT_H_A_B[i + 1];
						a = (v1 - v2) / (t1 - t2);
						b = v1 - a * t1;
						Result->H_A_B = a * (Params->temperature) + b;
						break;
					}
				}
		}
	}

	//H_A_F
	if (Result->H_A_F == UTC_MAT_QUERY_RESULT_REQUIRED ) {
		if (material->H_A_F_num == 1) {
			if (material->expression_of_H_A_F_from_T != NULL) {
#ifdef ANALYTIC_EXPRESSIONS
				material->Tfor_H_A_F[0] = Params->temperature;		// current temperature
				// calculating value using given formula
				Result->H_A_F = ((exprtk::expression<double> *)material->expression_of_H_A_F_from_T)->value();
#endif
			} else {
				Result->H_A_F = material->atT_H_A_F[0];
			}
		} else {
			if (Params->temperature <= material->Tfor_H_A_F[0])
				Result->H_A_F = material->atT_H_A_F[0];
			else if (Params->temperature >= material->Tfor_H_A_F[material->H_A_F_num - 1])
				Result->H_A_F = material->atT_H_A_F[material->H_A_F_num - 1];
			else
				for (i = 0; i < material->H_A_F_num - 1; ++i) {
					if (Params->temperature <= material->Tfor_H_A_F[i + 1]) {
						t1 = material->Tfor_H_A_F[i];
						t2 = material->Tfor_H_A_F[i + 1];
						v1 = material->atT_H_A_F[i];
						v2 = material->atT_H_A_F[i + 1];
						a = (v1 - v2) / (t1 - t2);
						b = v1 - a * t1;
						Result->H_A_F = a * (Params->temperature) + b;
						break;
					}
				}
		}
	}

	//H_A_P
	if (Result->H_A_P == UTC_MAT_QUERY_RESULT_REQUIRED ) {
		if (material->H_A_P_num == 1) {
			if (material->expression_of_H_A_P_from_T != NULL) {
#ifdef ANALYTIC_EXPRESSIONS
				material->Tfor_H_A_P[0] = Params->temperature;		// current temperature
				// calculating value using given formula
				Result->H_A_P = ((exprtk::expression<double> *)material->expression_of_H_A_P_from_T)->value();
#endif
			} else {
				Result->H_A_P = material->atT_H_A_P[0];
			}
		} else {
			if (Params->temperature <= material->Tfor_H_A_P[0])
				Result->H_A_P = material->atT_H_A_P[0];
			else if (Params->temperature >= material->Tfor_H_A_P[material->H_A_P_num - 1])
				Result->H_A_P = material->atT_H_A_P[material->H_A_P_num - 1];
			else
				for (i = 0; i < material->H_A_P_num - 1; ++i) {
					if (Params->temperature <= material->Tfor_H_A_P[i + 1]) {
						t1 = material->Tfor_H_A_P[i];
						t2 = material->Tfor_H_A_P[i + 1];
						v1 = material->atT_H_A_P[i];
						v2 = material->atT_H_A_P[i + 1];
						a = (v1 - v2) / (t1 - t2);
						b = v1 - a * t1;
						Result->H_A_P = a * (Params->temperature) + b;
						break;
					}
				}
		}
	}

#endif

	//thermal_conductivity
	if (Result->thermal_conductivity == UTC_MAT_QUERY_RESULT_REQUIRED) {
		if (material->thermal_conductivity_num == 1)
			Result->thermal_conductivity = material->atT_thermal_conductivity[0];
		else {
			if (Params->temperature <= material->Tfor_thermal_conductivity[0])
				Result->thermal_conductivity = material->atT_thermal_conductivity[0];
			else if (Params->temperature >=
					material->Tfor_thermal_conductivity[material->thermal_conductivity_num - 1])
				Result->thermal_conductivity =
						material->atT_thermal_conductivity[material->thermal_conductivity_num - 1];
			else
				for (i = 0; i < material->thermal_conductivity_num - 1; ++i) {
					if (Params->temperature <= material->Tfor_thermal_conductivity[i + 1]) {
						t1 = material->Tfor_thermal_conductivity[i];
						t2 = material->Tfor_thermal_conductivity[i + 1];
						v1 = material->atT_thermal_conductivity[i];
						v2 = material->atT_thermal_conductivity[i + 1];
						a = (v1 - v2) / (t1 - t2);
						b = v1 - a * t1;
						Result->thermal_conductivity = a * (Params->temperature) + b;
						break;
					}
				}
		}
	}

	//specific_heat
	if (Result->specific_heat == UTC_MAT_QUERY_RESULT_REQUIRED) {
		if (material->specific_heat_num == 1)
			Result->specific_heat = material->atT_specific_heat[0];
		else {
			if (Params->temperature <= material->Tfor_specific_heat[0])
				Result->specific_heat = material->atT_specific_heat[0];
			else if (Params->temperature >=
					material->Tfor_specific_heat[material->specific_heat_num - 1])
				Result->specific_heat =
						material->atT_specific_heat[material->specific_heat_num - 1];
			else
				for (i = 0; i < material->specific_heat_num - 1; ++i) {
					if (Params->temperature <= material->Tfor_specific_heat[i + 1]) {
						t1 = material->Tfor_specific_heat[i];
						t2 = material->Tfor_specific_heat[i + 1];
						v1 = material->atT_specific_heat[i];
						v2 = material->atT_specific_heat[i + 1];
						a = (v1 - v2) / (t1 - t2);
						b = v1 - a * t1;
						Result->specific_heat = a * (Params->temperature) + b;
						break;
					}
				}
		}
	}

	//thermal_expansion_coefficient
	if (Result->thermal_expansion_coefficient == UTC_MAT_QUERY_RESULT_REQUIRED) {
		if (material->thermal_expansion_coefficient_num == 1)
			Result->thermal_expansion_coefficient =
					material->atT_thermal_expansion_coefficient[0];
		else {
			if (Params->temperature <= material->Tfor_thermal_expansion_coefficient[0])
				Result->thermal_expansion_coefficient =
						material->atT_thermal_expansion_coefficient[0];
			else if (Params->temperature >=
					material->Tfor_thermal_expansion_coefficient[material->thermal_expansion_coefficient_num - 1])
				Result->thermal_expansion_coefficient =
						material->atT_thermal_expansion_coefficient[material->thermal_expansion_coefficient_num - 1];
			else
				for (i = 0; i < material->thermal_expansion_coefficient_num - 1; ++i) {
					if (Params->temperature <=
							material->Tfor_thermal_expansion_coefficient[i + 1]) {
						t1 = material->Tfor_thermal_expansion_coefficient[i];
						t2 = material->Tfor_thermal_expansion_coefficient[i + 1];
						v1 = material->atT_thermal_expansion_coefficient[i];
						v2 = material->atT_thermal_expansion_coefficient[i + 1];
						a = (v1 - v2) / (t1 - t2);
						b = v1 - a * t1;
						Result->thermal_expansion_coefficient = a * (Params->temperature) + b;
						break;
					}
				}
		}
	}

	//electrical_resistivity
	if (Result->electrical_resistivity == UTC_MAT_QUERY_RESULT_REQUIRED) {
		if (material->electrical_resistivity_num == 1)
			Result->electrical_resistivity = material->atT_electrical_resistivity[0];
		else {
			if (Params->temperature <= material->Tfor_electrical_resistivity[0])
				Result->electrical_resistivity = material->atT_electrical_resistivity[0];
			else if (Params->temperature >=
					material->Tfor_electrical_resistivity[material->electrical_resistivity_num - 1])
				Result->electrical_resistivity =
						material->atT_electrical_resistivity[material->electrical_resistivity_num - 1];
			else
				for (i = 0; i < material->electrical_resistivity_num - 1; ++i) {
					if (Params->temperature <=
							material->Tfor_electrical_resistivity[i + 1]) {
						t1 = material->Tfor_electrical_resistivity[i];
						t2 = material->Tfor_electrical_resistivity[i + 1];
						v1 = material->atT_electrical_resistivity[i];
						v2 = material->atT_electrical_resistivity[i + 1];
						a = (v1 - v2) / (t1 - t2);
						b = v1 - a * t1;
						Result->electrical_resistivity = a * (Params->temperature) + b;
						break;
					}
				}
		}
	}

	//enthalpy
	if (Result->enthalpy == UTC_MAT_QUERY_RESULT_REQUIRED) {
		if (material->enthalpy_num == 1)
			Result->enthalpy = material->atT_enthalpy[0];
		else {
			if (Params->temperature <= material->Tfor_enthalpy[0])
				Result->enthalpy = material->atT_enthalpy[0];
			else if (Params->temperature >=
					material->Tfor_enthalpy[material->enthalpy_num - 1])
				Result->enthalpy = material->atT_enthalpy[material->enthalpy_num - 1];
			else
				for (i = 0; i < material->enthalpy_num - 1; ++i) {
					if (Params->temperature <= material->Tfor_enthalpy[i + 1]) {
						t1 = material->Tfor_enthalpy[i];
						t2 = material->Tfor_enthalpy[i + 1];
						v1 = material->atT_enthalpy[i];
						v2 = material->atT_enthalpy[i + 1];
						a = (v1 - v2) / (t1 - t2);
						b = v1 - a * t1;
						Result->enthalpy = a * (Params->temperature) + b;
						break;
					}
				}
		}
	}

	//VOF */
	if (Result->VOF == UTC_MAT_QUERY_RESULT_REQUIRED) {
		if (material->VOF_num == 1)
			Result->VOF = material->atT_VOF[0];
		else {
			if (Params->temperature <= material->Tfor_VOF[0])
				Result->VOF = material->atT_VOF[0];
			else if (Params->temperature >= material->Tfor_VOF[material->VOF_num - 1])
				Result->VOF = material->atT_VOF[material->VOF_num - 1];
			else
				for (i = 0; i < material->VOF_num - 1; ++i) {
					if (Params->temperature <= material->Tfor_VOF[i + 1]) {
						t1 = material->Tfor_VOF[i];
						t2 = material->Tfor_VOF[i + 1];
						v1 = material->atT_VOF[i];
						v2 = material->atT_VOF[i + 1];
						a = (v1 - v2) / (t1 - t2);
						b = v1 - a * t1;
						Result->VOF = a * (Params->temperature) + b;
						break;
					}
				}
		}
	}

	//dg_dT */
	if (Result->dg_dT == UTC_MAT_QUERY_RESULT_REQUIRED) {
		if (material->dg_dT_num == 1)
			Result->dg_dT = material->atT_dg_dT[0];
		else {
			if (Params->temperature <= material->Tfor_dg_dT[0])
				Result->dg_dT = material->atT_dg_dT[0];
			else if (Params->temperature >= material->Tfor_dg_dT[material->dg_dT_num - 1])
				Result->dg_dT = material->atT_dg_dT[material->dg_dT_num - 1];
			else
				for (i = 0; i < material->dg_dT_num - 1; ++i) {
					if (Params->temperature <= material->Tfor_dg_dT[i + 1]) {
						t1 = material->Tfor_dg_dT[i];
						t2 = material->Tfor_dg_dT[i + 1];
						v1 = material->atT_dg_dT[i];
						v2 = material->atT_dg_dT[i + 1];
						a = (v1 - v2) / (t1 - t2);
						b = v1 - a * t1;
						Result->dg_dT = a * (Params->temperature) + b;
						break;
					}
				}
		}
	}

	if (Result->temp_solidus == UTC_MAT_QUERY_RESULT_REQUIRED) {
		Result->temp_solidus = material->temp_solidus;
	}
	if (Result->temp_liquidus  == UTC_MAT_QUERY_RESULT_REQUIRED) {
		Result->temp_liquidus = material->temp_liquidus;
	}
	if (Result->temp_vaporization  == UTC_MAT_QUERY_RESULT_REQUIRED) {
		Result->temp_vaporization = material->temp_vaporization;
	}
	if (Result->latent_heat_of_fusion == UTC_MAT_QUERY_RESULT_REQUIRED) {
		Result->latent_heat_of_fusion = material->latent_heat_of_fusion;
	}
	if (Result->latent_heat_of_vaporization == UTC_MAT_QUERY_RESULT_REQUIRED) {
		Result->latent_heat_of_vaporization = material->latent_heat_of_vaporization;
	}
	return 0;
}



int utr_mat_get_matID(int groupID)
{
	int material_idx = UTE_GROUP_ID_NOT_ASSIGNED;


	if (groupID < utv_group2material.size()) {
		material_idx = utv_group2material[groupID];
	}

	if (material_idx ==  UTE_GROUP_ID_NOT_ASSIGNED) {
		mf_log_warn("Requested 'material_number' for groupID=%d which is not assigned!", groupID );



		for (int i = 0; i < utv_materials.material_data.size(); ++i) {
			if (utv_materials.material_ids[i] == groupID) {
				material_idx = utv_materials.material_ids[i];
				break;
			}
		}



		if (groupID >= utv_group2material.size()) {
			#pragma omp critical
			{
				utv_group2material.resize(groupID + 1, UTE_GROUP_ID_NOT_ASSIGNED);
			}
		}
		// Assigning permanently to avoid next searches.
		utv_group2material[groupID] = material_idx;
		mf_log_info("Group %d assigned to material %d", groupID, material_idx);

	}
	//mf_log_info("Group_id(%d), Material_id(%d)",groupID,material_idx);
	if (material_idx == UTE_GROUP_ID_NOT_ASSIGNED) {
		mf_log_err("Material not found for given Group_id (%d) ", groupID);
	}

	return material_idx;

}

int utr_mat_get_n_materials()
{
	return utv_materials.material_data.size();
}

/**
 * @brief utr_mat_get_materials_IDs To get array of all known materialIDs
 * @param materialIDs array of size returned by \link utr_mat_get_n_materials \endlink
 * @return number of written material IDs into array \param materialIDs
 */
int utr_mat_get_materials_IDs(int * materialIDs)
{
	int n_written_mats = UTE_FAIL;
	if (materialIDs != NULL) {
		utt_materials::utt_material_map::const_iterator it = utv_materials.material_data.begin();
		for (n_written_mats = 0; n_written_mats < utv_materials.material_data.size(); ++n_written_mats) {
			materialIDs[n_written_mats] = it->first;    //< this is material ID
			++it;
		}
	}
	return n_written_mats;
}
void utr_mat_clear_all()
{
	utt_materials::utt_material_map::iterator it = utv_materials.material_data.begin();

	while (it != utv_materials.material_data.end()) {
		utr_mat_clear_material(& (it->second) );
		++it;
	}
	utv_materials.material_data.clear();
	utv_materials.material_ids.clear();
	utv_group2material.clear();
	utv_read_material_files.clear();

}

const utt_material_data * utr_mat_get_material(const int groupID)
{
	mfp_check_debug(utr_mat_get_matID(groupID) >= 0, "Material id incorrect(%d)", utr_mat_get_matID(groupID) );

//    mfp_check_debug( ( utv_materials.material_data.find(utr_mat_get_matID(groupID)) !=  utv_materials.material_data.end() ),
//                    "Material id incorrect(%d)",
//                     utr_mat_get_matID(groupID) );

	return & ( utv_materials.material_data[utr_mat_get_matID(groupID)] );
}

const utt_material_data * utr_mat_get_material_by_matID(int matID)
{
	mfp_check_debug(matID >= 0, "Material id incorrect(%d)", matID);
	mfp_check_debug(matID <= *std::max_element(utv_materials.material_ids.begin(),
					utv_materials.material_ids.end()),
			"Material id incorrect(%d)", matID);

	return & ( utv_materials.material_data[matID] );
}

#ifdef __cplusplus
}
#endif

#endif //UTS_MAT_CPP
