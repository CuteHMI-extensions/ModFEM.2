/************************************************************************
File apph_intf.h - interface of parallel approximation modules of the code.

Procedures:
  appr_get_ent_owner - to return owning process(or) ID
  appr_get_ent_id_at_owner

  appr_init_exchange_tables - to initialize data structure related to exchange 
                             of dofs
  appr_free_exchange_tables - to free data structure related to exchange 
                             of dofs

  appr_get_nr_dofs_owned - utility to calculate the number of dofs owned at
                          the finest level

Global functions that can be (should be?) parallelized
  appr_create_exchange_tables - to create lists of dof structures
                              exchanged between pairs of processors.
  appr_exchange_dofs - to exchange dofs between processors
  appr_sol_vec_norm - to compute a norm of global vector in parallel
  appr_sol_sc_prod - to compute a scalar product of two global vectors 

------------------------------  			
History:        
	08.2008 - Krzysztof Banas, initial version		
*************************************************************************/

#ifndef _apph_intf_
#define _apph_intf_

#ifdef __cplusplus
extern "C" {
#endif

/** @defgroup APPR Parallel Approximation
 *
 *  @{
 */

/**--------------------------------------------------------
appr_get_ent_owner - to return owning process(or) ID
---------------------------------------------------------*/
extern int appr_get_ent_owner( 
                   /** returns: >=0 -success code - owner ID, <0 -error code */
  int Field_id,    /** in: approximation field ID  */
  int Ent_type,      /** in: type of mesh entity */
  int Ent_id         /** in: mesh entity ID */
  );

/**--------------------------------------------------------
appr_get_ent_id_at_owner
---------------------------------------------------------*/
extern int appr_get_ent_id_at_owner( 
        /** returns: >=0 -success code - entity local ID, <0 -error code */
  int Field_id,    /** in: approximation field ID  */
  int Ent_type,      /** in: type of mesh entity */
  int Ent_id         /** in: mesh entity ID */
  );


/**--------------------------------------------------------
  appr_init_exchange_tables - to initialize data structure related to exchange 
                             of dofs
---------------------------------------------------------*/
extern int appr_init_exchange_tables(/* returns: ID of created tables */
  int Nr_proc,      /* number of process(or)s */
  int My_proc_id,   /* executing process(or)s ID */
  int Field_id,
  int Nreq_offset,
  int Nreq,
  int Nr_levels    /* in: number of levels for each field */
  );

/**--------------------------------------------------------
  appr_create_exchange_tables - to create lists of dofs
                              exchanged between pairs of processors.
  REMARK. a simplified setting is used, where dof structures are modified
    (updated) only by their owners and data send to other processors.
    Thus both exchange1 and exchange2 groups of arrays are empty and
    only send and receive arrays are filled. In the routine
    appr_create_exchange_tables each processor constructs lists of dof
    structures needed for solver (dofrecv group) and sends them to owning 
    processors. It also receives from other processors requests for owned 
    dof structures needed by other processors (dofsend group).
---------------------------------------------------------*/
extern int appr_create_exchange_tables( 
                   /** returns: >=0 -success code, <0 -error code */
  int Table_index,    /* in: exchange table index  */
  int Level_id,    /** in: level ID */
  int Nr_dof_ent,  /** in: number of DOF entities in the level */
  /** all four subsequent arrays are indexed by block IDs with 1(!!!) offset */
  int* L_dof_ent_type,/** in: list of DOF entities associated with DOF blocks */
  int* L_dof_ent_id,  /** in: list of DOF entities associated with DOF blocks */
  int* L_bl_nrdof,    /** in: list of nrdofs for each dof block */
  int* L_bl_posg,     /** in: list of positions within the global */
                      /**     vector of dofs for each dof block */
  int* L_elem_to_bl,  /** in: list of DOF blocks associated with DOF entities */
  int* L_face_to_bl,  /** in: list of DOF blocks associated with DOF entities */
  int* L_edge_to_bl,  /** in: list of DOF blocks associated with DOF entities */
  int* L_vert_to_bl  /** in: list of DOF blocks associated with DOF entities */
  );

/**--------------------------------------------------------
  appr_exchange_dofs - to exchange dofs between processors
---------------------------------------------------------*/
extern int appr_exchange_dofs(
  int Table_index,    /* in: exchange table index  */
  int Level_id,    /** in: level ID */
  double* Vec_dofs  /** in: vector of dofs to be exchanged */
  );

/**--------------------------------------------------------
  appr_sol_vec_norm - to compute a norm of global vector in parallel
---------------------------------------------------------*/
extern double appr_sol_vec_norm( /** returns: L2 norm of global Vector */
  int Table_index,    /* in: exchange table index  */
  int Level_id,    /** in: level ID */
  int Nrdof,            /** in: number of vector components */
  double* Vector        /** in: local part of global Vector */
  );

/**--------------------------------------------------------
  appr_sol_sc_prod - to compute a scalar product of two global vectors 
---------------------------------------------------------*/
extern double appr_sol_sc_prod(
		   /** retruns: scalar product of Vector1 and Vector2 */
  int Table_index,    /* in: exchange table index  */
  int Level_id,    /** in: level ID */
  int Nrdof,           /** in: number of vector components */
  double* Vector1,     /** in: local part of global Vector */
  double* Vector2      /** in: local part of global Vector */
  );


/**--------------------------------------------------------
  appr_get_nr_dofs_owned - utility to calculate the number of dofs owned at
                          the finest level
---------------------------------------------------------*/
extern int appr_get_nr_dofs_owned( 
  int Table_index    /* in: exchange table index  */
  );

/**--------------------------------------------------------
  appr_free_exchange_tables - to free data structure related to exchange 
                             of dofs
---------------------------------------------------------*/
extern int appr_free_exchange_tables( 
		  /** returns: >=0 -success code, <0 -error code */
  int Table_index    /* in: exchange table index  */
  );

/** @} */ // end of group


#ifdef __cplusplus
}
#endif

#endif
