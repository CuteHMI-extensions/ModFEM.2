#ifndef UTH_BC_H_
#define UTH_BC_H_

#include "uth_mesh.h"

#ifdef __cplusplus
extern "C"
{
#endif

/**
 \defgroup UTM_BC Boundary Condition Utilities
 \ingroup UTM
*
* The boundary conditions submodule of the utilities module is responsible for:
* 1. \link utr_bc_read_block_assigments
* reading the block to group assigment from bc file(s)  \endlink
* and \link utr_bc_get_blockID providing access to the block-to-group mapping. \endlink.
*
* The correct way of using the material database is as follows:
*
*
* 1.1. At program initialisation: \link utr_bc_read_block_assigments Read the bc file(s) with mapping(s)  \endlink at least once.
* - example mapping block:
*    group_to_blocks_assignment: (
*    {
*    block_number = 1;
*    groups = (15,11);
*    },
*    {
*    block_number = 35;
*    groups = ( 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 );
*    });
*
* 1.2. During computations: call \link utr_bc_get_blockID \endlink to access provided mappings.
*
* 2. Inserting \link utt_bc_assignment boundary conditions into mesh at run-time \endlink.
* 'bo_to_insert' block is expected to be in any boundary conditions file.
* If there are given two ID's, then in mesh between them (with respect to their kind), a new boundary condition with number 'bcnum' will appear.
* This, however, generates a split in mesh, introducing a new 'contact' plane between selected ID's of given kind.
*
* bc_to_insert: (
*    {
*	bcnum = 2;
*	concerns="block";
*	ID1 = 1;
*	ID2 = 2;
*	}
*  );
*
*/


//////////////////// BC TO INSERT ////////////////////////////

typedef struct {
    int  BC_num;
    int  Ids[2];
    utt_mesh_bc_type  Id_type; /// see \link utt_mesh_bc_type \endlink;
} utt_bc_assignment;


int utr_bc_read_to_insert_assigments(const char *Work_dir,
                           const char *Filename,
                           FILE *Interactive_output);

int utr_bc_to_insert_n_assigments();

const utt_bc_assignment* utr_bc_to_insert_get_assigments( );

int utr_bc_to_insert_completed();


//////////////////// BLOCKS ////////////////////////////
///
/// \brief utr_bc_read_assigments
/// \param Work_dir
/// \param Filename
/// \param Interactive_output
/// \return
///
int utr_bc_read_block_assigments(const char *Work_dir,
                           const char *Filename,
                           FILE *Interactive_output);

/** \brief To get block number associated with given element group number.
 *
 * \param groupID int see \link mmr_el_groupID \endlink
 * \return int block number (ID)
 *
 */
int utr_bc_get_blockID(int groupID);

/** \brief
 *
 * \return int
 *
 */
int utr_bc_get_n_block_assignments();

///** \brief
// *
// * \param number const int
// * \param bc_num_to_set utt_bc_assignment*
// * \return int
// *
// */
//int utr_bc_get_assignment(const int number, utt_bc_assignment * bc_num_to_set);

double utr_temp2block_id(int id);
int utr_temp2block_size();

void utr_bc_clear_all();

/** @} */ // end of group

#ifdef __cplusplus
}
#endif

#endif // UTH_BC_H_
