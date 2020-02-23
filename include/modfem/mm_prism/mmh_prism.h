/************************************************************************
File mmh_prism.h - internal information for mesh manipulation module
				   for prismatic elements

Contains:
  - constants
  - data types
  - global variables (for the whole module)
  - function headers

------------------------------
History:
		02.2002 - Krzysztof Banas, initial version
*************************************************************************/

#ifndef _mmh_prism_
#define _mmh_prism_


/*** CONSTANTS ***/

#define MMC_MAX_NUM_MESH   10   /* maximal number of meshes */

/* Other */
extern const int MMC_SAME_ORIENT;   /* indicator for the same orientation */
extern const int MMC_OPP_ORIENT;    /* indicator for the opposite orientation */


/*** DATA TYPES ***/

/*** Mesh data structure (for prismatic elements but with provisions
						  for hexahedral and tetrahedral elements)   ***/

/* structure with mesh parameters */
typedef struct {
  int    nrno;        /* number of active nodes */
  int    nmno;        /* maximal index of initiated node */
  int    mxno;        /* maximal available index of node */
  int    pfno;        /* pointer to first free node space*/
  int    nred;        /* number of active edges */
  int    nmed;        /* maximal index of initiated edge */
  int    mxed;        /* maximal available index of edge */
  int    pfed;        /* pointer to first free edge space*/
  int    nrfa;        /* number of active faces */
  int    nmfa;        /* maximal index of initiated face */
  int    mxfa;        /* maximal available index of face */
  int    pffa;        /* pointer to first free face space */
  int    nrel;        /* number of active elements */
  int    nmel;        /* maximal index of initiated element */
  int    mxel;        /* maximal available index of element */
  int    pfel;        /* pointer to first free element space */
  int    maxgen;      /* maximum generation level for elements */
  int    maxgendiff;  /* maximum difference of generation levels */
					  /* for neighboring elements */
} mmt_meshp;

/* node structure */
typedef struct {
  double   x;        /* x-coor of the node (or -1e11 for free spaces ) */
  double   y;        /* y-coor of the node (or next free space ID) */
  double   z;        /* z-coor of the node */
} mmt_nodes;

/* edge structure */
typedef struct {
  int type;     /*  1 (MMC_EDGE) - active edge */
				/* <0 - inactive edge, number of attempted divisions */
				/*  0 (MMC_FREE) - free space */
  int node[2];  /* for active edges: node's IDs */
				/* for free spaces: node[0] - the next free space */
				/* for inactive edges: son's IDs */
  int* elems;   /* temporary table (not saved in dump files) storing */
				/* the list of elements containing the edge */
} mmt_edges;

/* face structure */
typedef struct {
  int type;       /* type and status of the face: */
				  /*   1 (MMC_EDGE) - 1D face for 2D elements, */
				  /*   2 (MMC_QUAD) - quadrilateral */
				  /*   3 (MMC_TRIA) - triangle, */
				  /*   0 (MMC_FREE) - free space */
				  /*  >0 - active, not refined */
				  /*  <0 - inactive, refined */
  int bc;         /* for boundary faces: bc flag for the face (obtained */
				  /*   from mesh generator and interpreted by the problem */
				  /*   dependent module) */
				  /* for interelement faces: the shift in node numbering */
				  /*   between the face and its second neighbor (if it has */
				  /*   different orientation than the face,  */
				  /*   0 if only possible and for the same orientation) */
	/* WARNING: for quadrilateral faces only shift=1 is currently implemented */
				  /* for free spaces: next free space's ID */
  int edge[4];    /* edges' IDs - in proper order */
				  /*  >0 - edge with the same orientation as face */
				  /*  <0 - edge with opposite orientation */
  int neig[2];    /* neighboring elements' global IDs */
				  /* (first neighbor has the same orientation */
				  /*  as the face and ordering of its nodes */
				  /*  defines ordering of nodes on the face) */
				  /*  >0 - ID of equal size neighbor */
				  /*  -1 (MMC_BIG_NGB) - indicates big neighbor */
				  /*   0 (MMC_BOUNDARY) - boundary (always second neighbor) */
				  /*  -2 (MMC_SUB_BND) - inter-subdomain boundary */
  int *sons;      /* sons */
} mmt_faces;

/* element structure */
typedef struct {
  int    type;    /* 7 (MMC_TETRA) - tetrahedron */
				  /* 5 (MMC_PRISM) - prism */
				  /* 6 (MMC_BRICK) - brick */
				  /* 0 (MMC_FREE) - free space (face and sons not allocated) */
				  /* >0 - active, not refined */
				  /* <0 - inactive, refined */
  int    mate;    /* material ID */
				  /* for free space: ID of the next free structure */
  int    fath;    /* parent's global ID */
				  /* 0 (MMC_NO_FATH) - if the element has no parent */
  int    refi;    /* type of the last refinement */
				  /* 0 (MMC_NOT_REF) - sons not allocated */
				  /* 1 (MMC_REF_TYP_ISO) - isotropic refinement */
  int    *face;   /* faces' global IDs (sign=orientation)*/
  int    *sons;   /* children */
} mmt_elems;

/* structure with mesh data */
typedef struct {
  int     space_dim;        /* number of space dimensions */
  int       monitor;        /* printing flag */
  mmt_meshp    parm;        /* structure with mesh parameters */
  mmt_nodes   *node;        /* pointer to array of nodes */
  mmt_edges   *edge;        /* pointer to array of edges */
  mmt_faces   *face;        /* pointer to array of faces */
  mmt_elems   *elem;        /* pointer to array of elements */
} mmt_mesh;



/*** GLOBAL VARIABLES for the whole module ***/

extern int       mmv_nr_meshes;   /* the number of meshes in the problem */
extern int       mmv_cur_mesh_id;              /* ID of the current mesh */
extern mmt_mesh  mmv_meshes[MMC_MAX_NUM_MESH];        /* array of meshes */


/*** FUNCTIONS DECLARATIONS - headers for internal functions ***/

/* in file mms_prism_io.c:
  mmr_read_mesh - to dump-in mesh data stored by previous runs
  mmr_write_mesh - to dump-out mesh data in the standard HP_FEM format
*/

/**--------------------------------------------------------
  mmr_read_mesh - to dump-in mesh data stored by previous runs
---------------------------------------------------------*/
extern int mmr_read_mesh( /* returns: >=0 - success code, <0 - error code */
  int Mesh_id,	  /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  char *Filename  /* in: name of the file to read mesh data */
  );

/**--------------------------------------------------------
  mmr_write_mesh - to dump-out mesh data in the standard HP_FEM format
---------------------------------------------------------*/
extern int mmr_write_mesh( /* returns: >=0 - success code, <0 - error code */
  int Mesh_id,    /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  char *Filename  /* in: name of the file to write mesh data */
  );

/* in file mms_prism_io.c:
  mmr_import_mesh_grad - to read mesh data from input file created by
	"gradmesh" mesh generator (jkucwaj@zms.pk.edu.pl)
		 and generate structured 3D mesh of prismatic elements
*/

/**--------------------------------------------------------
mmr_import_mesh_grad - to read mesh data from input file created by
	"gradmesh" and generate structured 3D mesh of prismatic elements
---------------------------------------------------------*/
int mmr_import_mesh_grad( /* returns: >=0 - success code, <0 - error code */
  int Mesh_id,	  /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  char *Filename  /* in: name of the file to read mesh data */
  );

/* in file mms_prism_datstr.c:
  mmr_select_mesh - to select the proper mesh
  mmr_get_mesh_i_params - to return mesh parameters
*/


/**--------------------------------------------------------
  mmr_select_mesh - to select the proper mesh
---------------------------------------------------------*/
extern mmt_mesh* mmr_select_mesh( /* returns pointer to the chosen mesh */
  int Mesh_id    /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  );

/**--------------------------------------------------------
  mmr_get_mesh_i_params - to return mesh parameters
---------------------------------------------------------*/
int mmr_get_mesh_i_params(  /* returns: >=0 - integer mesh parameter,
									 <0 - error code  */
  mmt_mesh* Mesh,  /* in: pointer to the mesh data structure */
  int Num          /* in: parameter number in control structure */
  );

/* in file mms_prism_util.c:
  mmr_create_edge_elems - to create for each edge a list of elements
				  to which it belongs
  mmr_delete_edge_elems - to delete for each edge a list of elements
				  to which it belongs

*/

/**-----------------------------------------------------------
  mmr_create_edge_elems - to create for each edge a list of elements
				  to which it belongs
------------------------------------------------------------*/
int mmr_create_edge_elems( /* returns: 1-success, <=0-failure */
  int  Mesh_id,	 /* in: ID of the mesh to be used or 0 for the current mesh */
  int Max_edge_id /* in: the range of edge IDs to consider or 0 for default */
  );

/**-----------------------------------------------------------
  mmr_delete_edge_elems - to delete for each edge a list of elements
				  to which it belongs
------------------------------------------------------------*/
int mmr_delete_edge_elems( /* returns: 1-success, <=0-failure */
  int  Mesh_id,	 /* in: ID of the mesh to be used or 0 for the current mesh */
  int Max_edge_id /* in: the range of edge IDs to consider or 0 for default */
  );


/* in file mms_prism_ref.c:
  mmr_divide_el8_p - to break a prismatic element into 8 sons
  mmr_clust_el8_p - to cluster back a family of eight prisms
*/

extern int mmr_divide_el8_p( /* returns: >=0 - success code, <0 - error code */
  int Mesh_id,    /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int El        /* in: element ID */
);

extern int mmr_clust_el8_p( /* returns: >=0 - success code, <0 - error code */
  int Mesh_id,    /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int El        /* in: element ID */
);

/* in file mms_prism_util.c:
  mmr_clust_fa4_t - to cluster back a family of 4 triangular faces
  mmr_clust_fa4_q - to cluster back a family of 4 quadrilateral faces
  mmr_divide_face4_t - to break a triangular face into 4 sons
  mmr_divide_face4_q - to divide a quadrilateral face into 4 sons
*/

extern int mmr_divide_face4_t( /* returns: >=0 - success code, <0 - error code */
  int  Mesh_id,	   /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int  Fa,         /* in: face to be divided */
  int* Face_sons,  /* out: face's sons */
  int* Sons_edges, /* out: created new edges */
  int* New_nodes   /* out: created new nodes */
);

extern int mmr_divide_face4_q(/* returns: >=0 - success code, <0 - error code */
  int  Mesh_id,	   /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int  Fa,         /* in: face to be divided */
  int* Face_sons,  /* out: sons */
  int* Sons_edges, /* out: created new edges */
  int* New_nodes   /* out: created new nodes */
);

extern int mmr_clust_fa4_t( /* returns: >=0 - success code, <0 - error code */
  int  Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Face,      /* in: face ID */
  int* Face_sons /* in (optional): face sons' IDs */
);

extern int mmr_clust_fa4_q( /* returns: >=0 - success code, <0 - error code */
  int  Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Face,      /* in: face ID */
  int* Face_sons /* in (optional): face sons' IDs */
);


/* in file mms_util.c:
  mmr_chk_list - list manipulation
  mmr_ivector - to allocate space for vectors
  mmr_vec3_prod - to compute vector product of 3D vectors
  mmr_vec3_mxpr - to compute mixed vector product of 3D vectors
  mmr_vec3_length - to compute length of a 3D vector
*/

/**--------------------------------------------------------
mmr_chk_list - to check whether Num is on the list List
	with length Ll
---------------------------------------------------------*/
extern int mmr_chk_list(	/* returns: */
			/* >0 - position on the list */
					/* 0 - not found on the list */
	int Num, 	/* number to be checked */
	int* List, 	/* list of numbers */
	int Ll		/* length of the list */
	);

/**--------------------------------------------------------
 mmr_ivector - to allocate an integer vector: name[0..ncom-1]:
				  name=mmr_ivector(ncom,error_text)
---------------------------------------------------------*/
int *mmr_ivector(/* returns: pointer to array of integers */
	int ncom, 		/* in: number of components */
	char error_text[]	/* in: text to print in case of error */
	);

/**--------------------------------------------------------
mmr_vec3_prod - to compute vector product of 3D vectors
---------------------------------------------------------*/
extern void mmr_vec3_prod(
	double* vec_a, 	/* in: vector a */
	double* vec_b, 	/* in: vector b */
	double* vec_c	/* out: vector product axb */
	);

/**--------------------------------------------------------
mmr_vec3_mxpr - to compute mixed vector product of 3D vectors
---------------------------------------------------------*/
extern double mmr_vec3_mxpr( /* returns: mixed product [a,b,c] */
	double* vec_a, 	/* in: vector a */
	double* vec_b, 	/* in: vector b */
	double* vec_c	/* in: vector c */
	);

/**--------------------------------------------------------
mmr_vec3_length - to compute length of a 3D vector
---------------------------------------------------------*/
extern double mmr_vec3_length(	/* returns: vector length */
	double* vec	/* in: vector */
	);

#endif
