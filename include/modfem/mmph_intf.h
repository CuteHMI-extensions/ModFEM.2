/************************************************************************
File mmph_intf.h - interface of parallel mesh manipulation modules of the code.

  mmpr_init_mesh - to initialize the parallel mesh data structure

  mmpr_el_owner - to get an owning processor identifier of an element
  mmpr_el_set_owner - to set an owning processor identifier of an element
  mmpr_fa_owner - to get an owning processor identifier of a face
  mmpr_fa_set_owner - to set an owning processor identifier of a face
  mmpr_ed_owner - to get an owning processor identifier of an edge
  mmpr_ed_set_owner - to set an owning processor identifier of an edge
  mmpr_ve_owner - to get an owning processor identifier of a vertex
  mmpr_ve_set_owner - to set an owning processor identifier of a vertex
  mmpr_el_id_at_owner - to get an local identifier of an element
  mmpr_el_set_id_at_owner - to set an local identifier of an element
  mmpr_fa_id_at_owner - to get an local identifier of a face
  mmpr_fa_set_id_at_owner - to set an local identifier of a face
  mmpr_ed_id_at_owner - to get an local identifier of an edge
  mmpr_ed_set_id_at_owner - to set an local identifier of an edge
  mmr_ve_id_at_owner - to get an local identifier of a vertex
  mmpr_ve_set_id_at_owner - to set an local identifier of a vertex

  mmpr_create_subdomains - to decompose the mesh and create subdomains

  mmpr_mesh_ent_ovl - to create tables (for each mesh entity type) of local
                     numbers and storing processors for each mesh entity owned
                     by a given processor but stored also on other processors
  mmpr_update_ipid - to update interprocessor ID for all mesh entities
                    using interprocessor communication
  mmpr_free_mesh_ovl - to free data structure with overlap info
  mmpr_update_ref_list - to update list of refined elements due to irregularity
                      constraint using inter-processor communication
  mmpr_exchange_list_ref - to update the list of refined elements
                          using inter-processor communication
  mmpr_balance_load - to ... (what you think?)

 ------------------------------
History:
	08.2008 - Krzysztof Banas, initial version
*************************************************************************/

#ifndef _mmph_intf_
#define _mmph_intf_


#ifdef __cplusplus
extern "C"{
#endif


/** @defgroup MMP Parallel Mesh
 *
 *  @{
 */

/*** CONSTANTS ***/
#define MMPC_MAX_NUM_SUB 1024     /** maximal number of processors/processes */
//#define MMPC_MAX_NUM_SUB 128     /** maximal number of processors/processes */

/*** FUNCTION DECLARATIONS - headers for external functions ***/

/**-----------------------------------------------------------
  mmpr_init_mesh - to initialize the parallel mesh data structure
------------------------------------------------------------*/
extern int mmpr_init_mesh(  /** returns: >0 - Mesh ID, <0 - error code */
  int Control,
  int Mesh_id,
  int Nr_proc,
  int My_proc_id
  );

/**-----------------------------------------------------------
  mmpr_check_mesh - to check the parallel mesh data structure
------------------------------------------------------------*/
extern int mmpr_check_mesh(  /** returns: >0 - Mesh ID, <0 - error code */
  int Mesh_id
  );


/**-------------------------------------------------------------------------
  mmpr_el_owner - to get an owning processor identifier of an element
---------------------------------------------------------------------------*/
extern int mmpr_el_owner(
	       /** returns: >0 -an owning processor identifier of an element */
	         /**          <0 - error code */
  int Mesh_id,	 /** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int El         /** in: global (within the subdomain) element ID */
  );


/**-------------------------------------------------------------------------
  mmpr_el_set_owner - to set an owning processor identifier of an element
---------------------------------------------------------------------------*/
extern int mmpr_el_set_owner( /** returns: >=0 - success code, <0 - error code */
  int Mesh_id,	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int El,	/** in: global (within the subdomain) element ID */
  int Owner /** in: owning processor identifier of the element */
  );


/**-------------------------------------------------------------------------
  mmpr_fa_owner - to get an owning processor identifier of a face
---------------------------------------------------------------------------*/
extern int mmpr_fa_owner(
		 /** returns: >0 - an owning processor identifier of a face */
	          /**          <0 - error code */
  int Mesh_id,	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Fa        /** in: global (within the subdomain)  ID */
  );

/**-------------------------------------------------------------------------
  mmpr_fa_set_owner - to set an owning processor identifier of a face
---------------------------------------------------------------------------*/
extern int mmpr_fa_set_owner(  /** returns: >=0 - success code, <0 - error code */
  int Mesh_id,  /** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Fa,       /** in: global (within the subdomain) face ID */
  int Owner      /** in: owning processor identifier of the face */
  );

/**-------------------------------------------------------------------------
  mmpr_ed_owner - to get an owning processor identifier of an edge
---------------------------------------------------------------------------*/
extern int mmpr_ed_owner(
                 /** returns: >0 an owning processor identifier of an edge */
	          /**          <0 - error code */
  int Mesh_id,    /** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Ed          /** in: global (within the subdomain) edge ID */
  );

/**-------------------------------------------------------------------------
  mmpr_ed_set_owner - to set an owning processor identifier of an edge
---------------------------------------------------------------------------*/
extern int mmpr_ed_set_owner(  /** returns: >=0 - success code, <0 - error code */
  int Mesh_id,  /** in: mesh ID or 0 (MMPC_CUR_MESH_ID) for the current mesh */
  int Ed,       /** in: global (within the subdomain) edge ID */
  int Owner      /** in: owning processor identifier of the edge */
  );


/**-------------------------------------------------------------------------
  mmpr_ve_owner - to get an owning processor identifier of a vertex
---------------------------------------------------------------------------*/
extern int mmpr_ve_owner(
		 /** returns: >0 - an owning processor identifier of a vertex */
	         /**          <0 - error code */
  int Mesh_id,  /** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Ve        /** in: global (within the subdomain) vertex ID */
  );


/**-------------------------------------------------------------------------
  mmpr_ve_set_owner - to set an owning processor identifier of a vertex
---------------------------------------------------------------------------*/
extern int mmpr_ve_set_owner(  /** returns: >=0 - success code, <0 - error code */
  int Mesh_id,  /** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Ve,       /** in: global (within the subdomain) vertex ID */
  int Owner      /** in: owning processor identifier of the vertex */
  );


/**-------------------------------------------------------------------------
  mmpr_el_id_at_owner - to get an local identifier of an element
---------------------------------------------------------------------------*/
extern int mmpr_el_id_at_owner(/** returns: >0 - an local identifier of an element */
	         /**          <0 - error code */
  int Mesh_id,	 /** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int El         /** in: global (within the subdomain) element ID */
  );

/**-------------------------------------------------------------------------
  mmpr_el_set_id_at_owner - to set an local identifier of an element
---------------------------------------------------------------------------*/
extern int mmpr_el_set_id_at_owner(
                    /** returns: >=0 - success code, <0 - error code */
  int Mesh_id,	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int El,	/** in: global (within the subdomain) element ID */
  int Id_At_Owner /** in: local identifier of the element */
  );


/**-------------------------------------------------------------------------
  mmpr_fa_id_at_owner - to get an local identifier of a face
---------------------------------------------------------------------------*/
extern int mmpr_fa_id_at_owner( /** returns: >0 - an local identifier of a face */
	          /**          <0 - error code */
  int Mesh_id,	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Fa        /** in: global (within the subdomain)  ID */
  );


/**-------------------------------------------------------------------------
  mmpr_fa_set_id_at_owner - to set an local identifier of a face
---------------------------------------------------------------------------*/
extern int mmpr_fa_set_id_at_owner(/** returns: >=0 - success code, <0 - error code*/
  int Mesh_id,  /** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Fa,       /** in: global (within the subdomain) face ID */
  int Id_At_Owner      /** in: local identifier of the face */
  );


/**-------------------------------------------------------------------------
  mmpr_ed_id_at_owner - to get an local identifier of an edge
---------------------------------------------------------------------------*/
extern int mmpr_ed_id_at_owner( /** returns: >0 an local identifier of an edge */
	          /**          <0 - error code */
  int Mesh_id,    /** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Ed          /** in: global (within the subdomain) edge ID */
  );

/**-------------------------------------------------------------------------
  mmpr_ed_set_id_at_owner - to set an local identifier of an edge
---------------------------------------------------------------------------*/
extern int mmpr_ed_set_id_at_owner(/** returns: >=0 - success code, <0 - error code*/
  int Mesh_id,  /** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Ed,       /** in: global (within the subdomain) edge ID */
  int Id_At_Owner      /** in: local identifier of the edge */
  );

/**-------------------------------------------------------------------------
  mmr_ve_id_at_owner - to get an local identifier of a vertex
---------------------------------------------------------------------------*/
extern int mmpr_ve_id_at_owner(/** returns: >0 - a local identifier of a vertex */
	         /**          <0 - error code */
  int Mesh_id,  /** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Ve        /** in: global (within the subdomain) vertex ID */
  );

/**-------------------------------------------------------------------------
  mmpr_ve_set_id_at_owner - to set an local identifier of a vertex
---------------------------------------------------------------------------*/
extern int mmpr_ve_set_id_at_owner(
                   /** returns: >=0 - success code, <0 - error code */
  int Mesh_id,  /** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Ve,       /** in: global (within the subdomain) vertex ID */
  int Id_At_Owner      /** in: local identifier of the vertex */
  );


/**--------------------------------------------------------
  mmpr_create_subdomains - to decompose the mesh and create subdomains
---------------------------------------------------------*/
extern int mmpr_create_subdomains( /** returns: >=0 - success code, */
				  /**           <0 - error code */
  int Mesh_id,	     /** in: mesh ID */
  int Control        /** in: ??? */
  );

/**-----------------------------------------------------------
  mmpr_init_ref - to initialize the process of refinement
------------------------------------------------------------*/
extern int mmpr_init_ref(  /** returns: >=0 - success code, <0 - error code */
  int Mesh_id	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  );

/**-----------------------------------------------------------
  mmr_is_ready_for_proj_dof_ref - to check if mesh module is ready
                                  for dofs projection
------------------------------------------------------------*/
extern int mmpr_is_ready_for_proj_dof_ref(
  int Mesh_id
);

/**-----------------------------------------------------------
  mmpr_final_ref - to finalize the process of parallel refinement
------------------------------------------------------------*/
extern int mmpr_final_ref(  /** returns: >=0 - success code, <0 - error code */
  int Mesh_id 	/** in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  );

/**--------------------------------------------------------
  mmpr_mesh_ent_ovl - to create tables (for each mesh entity type) of local
                     numbers and storing processors for each mesh entity owned
                     by a given processor but stored also on other processors
---------------------------------------------------------*/
extern int mmpr_mesh_ent_ovl( /** returns: >=0 - success code, <0 - error code */
  int Mesh_id   /** in: mesh ID */
  );

/**--------------------------------------------------------
  mmpr_update_ipid - to update interprocessor ID for all mesh entities
                    using interprocessor communication
---------------------------------------------------------*/
extern int mmpr_update_ipid( /** returns: >=0 - processor ID, <0 - error code */
  int Mesh_id 	/** in: mesh ID */
  );

/**--------------------------------------------------------
  mmpr_free_mesh_ovl - to free data structure with overlap info
----------------------------------------------------------*/
extern int mmpr_free_mesh_ovl(
  int Mesh_id     /** in: mesh ID */
  );


/**-----------------------------------------------------------
  mmpr_balance_load - to ... (what you think?)
------------------------------------------------------------*/
extern int mmpr_balance_load(
  int Mesh_id  /** in: mesh ID */
  );


  // NEW PROCEDURES NOT SUPPORTED BY MMPD_PRISM

/**-----------------------------------------------------------
mmpr_update_ref_list - to update list of refined elements due to irregularity
                      constraints using inter-processor communication
------------------------------------------------------------*/
extern int mmpr_update_ref_list(
                   /** returns: >=0 - success code, <0 - error code */
  const int Mesh_id,     /** in: mesh ID */
  int *Nr_ref,     /** in/out: number of refined elements */
  int **List_ref   /** in/out: list of refined elements */
  );

/**-----------------------------------------------------------
  mmpr_ipc_n_connected_procs - to get number of neighburing subdomains assigned to other procs.
------------------------------------------------------------*/
extern int mmpr_n_neighbouring_procs( // returns number of procs connected with caller proc
        const int Mesh_id //in: mesh id
        );

/**-----------------------------------------------------------
  mmpr_interproc_connectivity - to get detailed info about neighbouring procs
------------------------------------------------------------*/
extern int mmpr_neighbouring_procs( // returns number of procs connected with caller proc
        const int Mesh_id,  //in: mesh id
        int * Neigh_proc_ids //out:[mmpr_ipc_n_connected_procs]
                                   // connected negibouring procs ids
        );

/**-----------------------------------------------------------
  mmpr_ipc_faces_connected_to_proc - to get array of face (ids) connected with given proc id
------------------------------------------------------------*/
extern int mmpr_ipc_faces_connected_to_proc( // returns number of faces
        const int Mesh_id, //in: mesh id
        const int Proc_id, //in: proc id to check
        int * Faces //out: [upto Connected_proc_n_links for Proc_id]
        );

/**-----------------------------------------------------------
  mmpr_ipc_fa_crossboundary_neig - to get info about element on the other side of boundary
------------------------------------------------------------*/
extern int mmpr_ipc_other_subdomain_neigh(const int Mesh_id,
        const int Fa,
        int *Neig_el_id_at_owner,
        int *Neig_el_owner
        );

/** @} */ // end of group


#ifdef __cplusplus
}
#endif



#endif
