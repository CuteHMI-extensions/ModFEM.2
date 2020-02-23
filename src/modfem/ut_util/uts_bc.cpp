
#include <vector>
#include <libconfig.h>
#include <algorithm>

#include <modfem/uth_log.h>    //USES
#include <modfem/uth_bc.h> // IMPELMENTS
#include <modfem/uth_mat.h> // ute_mat_groupID_flag


// HERE
//std::vector<utt_bc_assignment>  new_bc_nums;

std::vector<int>    utv_group2block;
std::vector<double>    utv_temp2block;
std::vector<utt_bc_assignment> utv_bc2insert;

bool operator==(const utt_bc_assignment & f, const utt_bc_assignment & s)
{
	// return TRUE if equal
	return 0 == memcmp(&f, &s, sizeof(utt_bc_assignment));
}

#ifdef __cplusplus
extern "C"
{
#endif

/**
 \defgroup UTM_BC Boundary Condition Utilities
 \ingroup UTM
*/

void utr_bc_clear_all()
{
	//new_bc_nums.clear();
	utv_group2block.clear();
	utv_temp2block.clear();
	utv_bc2insert.clear();
}

int utr_bc_read_to_insert_assigments(const char * Work_dir,
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
		mf_log_err("Could not open file! (%s) ", mapping_file);
		config_destroy(&cfg);
		return (UTE_FAIL);
	}

	root_setting = config_lookup(&cfg, "bc_to_insert");
	if (root_setting != NULL)
	{
		int count = config_setting_length(root_setting);
		int setting_length;

		mf_log_info("Number of bc_to_insert assigments in %s: %d", Filename, count);

		if (count > 0) {
			for (int i = 0; i < count; ++i) {
				config_setting_t * bc_to_insert = config_setting_get_elem(root_setting, i);
				config_setting_t * setting;

				utt_bc_assignment bc2insert_assigment;
				bc2insert_assigment.BC_num = UTE_MESH_BC_UNDEFINED;
				bc2insert_assigment.Ids[0] = UTE_MESH_BC_UNDEFINED;
				bc2insert_assigment.Ids[1] = UTE_MESH_BC_UNDEFINED;
				bc2insert_assigment.Id_type = UTE_MESH_BC_GROUP;


				// BC NUMBER
				if (NULL != (setting = config_setting_get_member(bc_to_insert, "bcnum"))) {
					bc2insert_assigment.BC_num  = (int)config_setting_get_int(setting);
				}
				else {
					mf_fatal_err("In 'bc_to_insert' assigment in file %s missing 'bcnum' in assigment %d ",
							mapping_file, i);
				}
				// What concerns ID1 and ID2.
				const char * concerns_string = NULL;
				setting = config_setting_get_member(bc_to_insert, "concerns");
				if (NULL != setting) {
					concerns_string = config_setting_get_string(setting);
					if (0 == strcmp("material", concerns_string )) {
						bc2insert_assigment.Id_type = UTE_MESH_BC_MATERIAL;
					}
					else if (0 == strcmp("block", concerns_string )) {
						bc2insert_assigment.Id_type = UTE_MESH_BC_BLOCK;
					}
					else if (0 != strcmp("group", concerns_string )) {
						mf_log_err("'bc_to_insert' wtih invalid 'concerns' parameter found (%s). Assuming IDs are referencing to group IDs.",
								concerns_string);
					}
				}
				else {
					mf_log_warn("'bc_to_insert' boundary condition wtih no 'concerns' parameter found. Assuming IDs are referencing to group IDs");
				}


				if (NULL != config_setting_get_member(bc_to_insert, "ID1")) {
					bc2insert_assigment.Ids[0] = config_setting_get_int(config_setting_get_member(bc_to_insert, "ID1"));
				}
				else {
					mf_fatal_err("'contact' boundary condition wtih no 'ID1' parameter found.");
				}


				if (NULL != config_setting_get_member(bc_to_insert, "ID2")) {
					bc2insert_assigment.Ids[1] = config_setting_get_int(config_setting_get_member(bc_to_insert, "ID2"));
				}
				else {
					mf_fatal_err("'contact' boundary condition wtih no 'ID2' parameter found.");
				}

				// if no such assigment exist already add this one.
				if (std::find(utv_bc2insert.begin(), utv_bc2insert.end(), bc2insert_assigment) == utv_bc2insert.end()) {
					utv_bc2insert.push_back(bc2insert_assigment);
					mf_log_info("'bc_to_insert' with {bcnum=%d;ID1=%d;ID2=%d;concerns=%s;} read from file %s.",
							bc2insert_assigment.BC_num, bc2insert_assigment.Ids[0], bc2insert_assigment.Ids[1], concerns_string, Filename);
				}
				else {
					mf_log_warn("Ignoring 'bc_to_insert' with {bcnum=%d;ID1=%d;ID2=%d;concerns=%s;} from file %s: it duplicates existing bc_to_insert assigment!",
							bc2insert_assigment.BC_num, bc2insert_assigment.Ids[0], bc2insert_assigment.Ids[1], concerns_string, Filename);
				}

			}
		}
		mf_log_info("'bc_to_insert' read with %d assigments! ", count);
	}


	config_destroy(&cfg);

	// REGISTERING deleting function
	atexit(utr_bc_clear_all);

	return (ret_val);

}


int utr_bc_to_insert_n_assigments()
{
	return utv_bc2insert.size();
}

int utr_bc_to_insert_completed()
{
	utv_bc2insert.clear();
	return 0;
}

const utt_bc_assignment * utr_bc_to_insert_get_assigments()
{
	return utv_bc2insert.data();
}

int utr_bc_get_blockID(int groupID)
{
	if (groupID <= 0) {
		mf_log_err("GroupID should be > 0 (is %d)", groupID);
	}
	int blockID = UTE_GROUP_ID_NOT_ASSIGNED;

	if (groupID < utv_group2block.size()) {
		blockID = utv_group2block[groupID];
	}

	if (blockID ==  UTE_GROUP_ID_NOT_ASSIGNED) {
		mf_log_warn("Requested blockID for groupID=%d which is not assigned!", groupID);
		blockID = UTE_GROUP_ID_DEFAULT_BLOCK;
		mf_log_info("Group %d assigned to default block %d", groupID, blockID);
	}

	return blockID;
}

double utr_temp2block_id(int id) {

	double ID = -1.0;

	if (id < utv_temp2block.size()) {
		ID = utv_temp2block[id];
	}

	if (ID ==  -1.0) {
		mf_log_warn("Requested utr_temp2block_id for id =%d which is not assigned!", id);
	}

	return ID;

}

int utr_temp2block_size() {

	return utv_temp2block.size();

}

int utr_bc_read_block_assigments(const char * Work_dir,
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
		mf_log_err("Could not open file! (%s) ", mapping_file);
		config_destroy(&cfg);
		return (UTE_FAIL);
	}

	root_setting = config_lookup(&cfg, "group_to_blocks_assignment");
	if (root_setting != NULL)
	{
		int count = config_setting_length(root_setting);
		int setting_length;

		mf_log_info("Number of block assigments in %s: %d", Filename, count);

		if (count > utv_group2block.size()) {
			utv_group2block.resize(count + 1, UTE_GROUP_ID_NOT_ASSIGNED);

			for (int i = 0; i < count; ++i) {
				config_setting_t * group_to_blocks = config_setting_get_elem(root_setting, i);
				config_setting_t * setting;

				int blocknum = 0;
				double startTemp = 0.0;

				// BLOCK NUMBER
				if (NULL != (setting = config_setting_get_member(group_to_blocks, "block_number"))) {
					blocknum  = (int)config_setting_get_int(setting);
				}
				else {
					mf_fatal_err("In 'group_to_blocks_assignment' in file %s missing 'block_number' in assigment %d ",
							mapping_file, i);
				}

				// Init Temp
				if (NULL != (setting = config_setting_get_member(group_to_blocks, "init_temp"))) {

					startTemp = 0.0;
					startTemp  = (double)config_setting_get_float(setting);
					utv_temp2block.push_back(double(blocknum));
					utv_temp2block.push_back(startTemp);

				}
				else {
					mf_log_warn("In 'group_to_blocks_assignment' in file %s missing 'init_temp' in assigment %d ",
							mapping_file, i);
				}

				// groups
				if (NULL != (setting = config_setting_get_member(group_to_blocks, "groups"))) {

					if (config_setting_is_list(setting) == CONFIG_TRUE) {
						setting_length = config_setting_length(setting);
						for (int j = 0; j < setting_length; ++j) {
							int  groupID = config_setting_get_int_elem(setting, j);

							if (groupID >= utv_group2block.size()) {
								utv_group2block.resize(groupID + 1, UTE_GROUP_ID_NOT_ASSIGNED);
							}

							if (utv_group2block[groupID] != UTE_GROUP_ID_NOT_ASSIGNED) {
								mf_log_warn("Duplicate 'group_to_blocks_assignment' for groupID=%d ", groupID);
							}

							utv_group2block[groupID] = blocknum;

							mf_log_info("Assigning group %d to block %d", groupID, blocknum);
						}
					}
					else {
						int  groupID = config_setting_get_int(setting);

						if (groupID >= utv_group2block.size()) {
							utv_group2block.resize(groupID + 1, UTE_GROUP_ID_NOT_ASSIGNED);
						}

						if (utv_group2block[groupID] != UTE_GROUP_ID_NOT_ASSIGNED) {
							mf_log_warn("Duplicate 'group_to_blocks_assignment' for groupID=%d ", groupID);
						}

						utv_group2block[groupID] = blocknum;

						mf_log_info("Assigning group %d to block %d", groupID, blocknum);

					}
				}
				else {
					mf_log_err("In 'group_to_blocks_assignment' in file %s missing 'groups' in assigment %d ", mapping_file, i);
				}

			}
		}
		mf_log_info("'group_to_blocks_assignment' read! ");
	}


	config_destroy(&cfg);

	// REGISTERING deleting function
	atexit(utr_bc_clear_all);

	return (ret_val);
}

int utr_bc_get_n_block_assignments()
{
	return utv_group2block.size();
}

//int utr_bc_get_assignment(const int number, utt_bc_assignment * bc_num_to_set)
//{
//    mf_check_mem(bc_num_to_set);

//    memcpy(bc_num_to_set, & new_bc_nums[number], sizeof(utt_bc_assignment) );

//    return 0;
//}


/** @} */ // end of group

#ifdef __cplusplus
}
#endif
