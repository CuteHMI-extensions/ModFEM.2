/************************************************************************
File mms_prism_datstr.c - definitions of mesh data structure routines

Contains routines:
  mmr_select_mesh - to select the proper mesh
  mmr_get_mesh_i_params - to return mesh parameters
  mmr_el_status - to return element status (active, inactive, free space)
  mmr_el_type - to return element type
  mmr_el_groupID - to return material number for element
  mmr_el_type_ref - to return element's type of refinement
  mmr_el_faces - to get faces and big neighbors of an element
  mmr_el_eq_neig - to get equal size (or larger) neighbors of an element
  mmr_el_node_coor - to get the coordinates of element's vertices
  mmr_el_edges - to get the list of element's edges
  mmr_el_fam - to return family information for an element
  mmr_el_gen - to return generation level for an element
  mmr_el_hsize - to compute a characteristic linear size for an element
  mmr_fa_status - to return face status (active, inactive, free space)
  mmr_fa_type - to return face type (triangle, quad, free space)
  mmr_fa_bc - to get the boundary condition type of the face
  mmr_fa_sub_bnd - to indicate face is on the intersubdomain boundary
  mmr_fa_set_sub_bnd - to set that face is on the intersubdomain boundary
  mmr_fa_edges - to return a list of face's edges
  mmr_fa_eq_neig - to return a list of face's neighbors and
		corresponding neighbors' sides (only equal size
		neighbors considered)
  mmr_fa_neig - to return a list of face's neighbors and
		corresponding neighbors' sides (for active
		faces all neighbors are active)
  mmr_fa_node_coor - to get the list and coordinates of faces's nodes
  mmr_fa_elem_coor - to find coordinates within neighboring
		elements for a point on face
  mmr_fa_area - to compute the area of face and vector normal
  mmr_fa_fam - to return face's family information
  mmr_edge_nodes - to return edge node's numbers
  mmr_edge_status - to return edge status (active, inactive, free space)
  mmr_node_status - to return node status (active, free space)
  mmr_node_coor - to return node coordinates
  mmr_move_node - moves node by vector
  mmr_set_node_coor - set node coordinates


Terminology convention:
  local - refers to an individual entity (e.g. local numbering of element's
		  nodes)
  global - refers to IDs within the domain (sequential) or subdomain (parallel)
  inter-processor - refers to IDs within the domain, for entities distributed
					among different processors, inter-processor ID includes
					information on entity's ownership

------------------------------
History:
	  02.2002 - Krzysztof Banas, initial version
*************************************************************************/

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<memory.h> // or string.h

/* interface for the mesh manipulation module */
#include <modfem/mmh_intf.h>

/* mesh manipulation data structure and headers for internal routines */
#include <modfem/mm_prism/mmh_prism.h>

/*---------------------------------------------------------
  mmr_select_mesh - to select the proper mesh
---------------------------------------------------------*/
mmt_mesh * mmr_select_mesh( /* returns pointer to the chosen mesh */
		/* to avoid errors if input is not valid */
		/* it returns the pointer to the current mesh */
		int Mesh_id    /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
)
{

	/* select the proper mesh from the array of meshes */
	if ( Mesh_id == MMC_CUR_MESH_ID ) {
		return (&mmv_meshes[mmv_cur_mesh_id - 1]);
	}
	else if ( Mesh_id > 0 && Mesh_id <= mmv_nr_meshes ) {
		return (&mmv_meshes[Mesh_id - 1]);
	}
	else {
		return (&mmv_meshes[mmv_cur_mesh_id - 1]);
		/* alternative:  return(NULL);   */
	}
}

/*---------------------------------------------------------
  mmr_get_mesh_i_params - to return mesh parameters
---------------------------------------------------------*/
int mmr_get_mesh_i_params(  /* returns: >=0 - integer mesh parameter,
									 <0 - error code  */
		mmt_mesh * Mesh, /* in: pointer to the mesh data structure */
		int Num          /* in: parameter number in control structure */
)
{

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	if (Num == 1) return (Mesh->parm.nrno);
	else if (Num == 2) return (Mesh->parm.nmno);
	else if (Num == 5) return (Mesh->parm.nred);
	else if (Num == 6) return (Mesh->parm.nmed);
	else if (Num == 9) return (Mesh->parm.nrfa);
	else if (Num == 10) return (Mesh->parm.nmfa);
	else if (Num == 13) return (Mesh->parm.nrel);
	else if (Num == 14) return (Mesh->parm.nmel);
	else if (Num == 21) return (Mesh->parm.maxgen);
	else if (Num == 22) return (Mesh->parm.maxgendiff);
	else {
		printf("Wrong parameter number in mesh_i_params!");
		exit(1);
	}

	/* error condition - that point should not be reached */
//	return -1;
}

/*---------------------------------------------------------
mmr_elem_structure - to return elem structure (e.g. for sending)
---------------------------------------------------------*/
int mmr_elem_structure( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id,	   /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int El,  	   /* in: elem ID */
		int * Elem_struct /* out: elem structure in the form of integer array */
)
{

	/* auxiliary variables */
	mmt_mesh * mesh;
	int i;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	Elem_struct[0] = mesh->elem[El].type;
	//Elem_struct[1]=mesh->elem[El].ipid;
	Elem_struct[2] = mesh->elem[El].mate;
	Elem_struct[3] = mesh->elem[El].fath;
	Elem_struct[4] = mesh->elem[El].refi;
	if (abs(mesh->elem[El].type) == MMC_PRISM) {
		for (i = 0; i < 5; i++) Elem_struct[5 + i] = mesh->elem[El].face[i];
	}
	else if (abs(mesh->elem[El].type) == MMC_BRICK) {
		for (i = 0; i < 6; i++) Elem_struct[5 + i] = mesh->elem[El].face[i];
	}
	else if (abs(mesh->elem[El].type) == MMC_TETRA) {
		for (i = 0; i < 4; i++) Elem_struct[5 + i] = mesh->elem[El].face[i];
	}
	if (mesh->elem[El].type < 0) {
		for (i = 0; i < 8; i++) {
			Elem_struct[11 + i] = mesh->elem[El].sons[i];
		}
	}

	return 1;

}

/*---------------------------------------------------------
mmr_el_status - to return element status (active, inactive, free space)
---------------------------------------------------------*/
int mmr_el_status( /* returns element status: */
		/* +1 (MMC_ACTIVE)   - active element */
		/*  0 (MMC_FREE)     - free space */
		/* -1 (MMC_INACTIVE) - inactive (refined) element */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int El  	/* in: element ID */
)
{

	/* auxiliary variables */
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	if (mesh->elem[El].type > 0) return (MMC_ACTIVE);
	else if (mesh->elem[El].type < 0) return (MMC_INACTIVE);
	else return (MMC_FREE);

}

/*---------------------------------------------------------
mmr_el_type - to return element type
---------------------------------------------------------*/
int mmr_el_type( /* returns element type or <0 - error code */
		/*	 7 (MMC_TETRA) - tetrahedron */
		/*	 5 (MMC_PRISM) - prism */
		/*	 6 (MMC_BRICK) - brick */
		/*	 0 (MMC_FREE)  - free space */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int El  	/* in: element ID */
)
{

	/* auxiliary variables */
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	return (abs(mesh->elem[El].type));
}

/*---------------------------------------------------------
mmr_el_groupID - to return material number for element
---------------------------------------------------------*/
int mmr_el_groupID( /* returns material flag for element, <0 - error code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int El  	/* in: element ID */
)
{

	/* auxiliary variables */
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	return (mesh->elem[El].mate);
}

/*---------------------------------------------------------
mmr_el_set_mate - to set material number for element
---------------------------------------------------------*/
int mmr_el_set_groupID( /* sets material flag for element, <0 - error code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int El,  	/* in: element ID */
		int Mat_id    /* in: material ID */
)
{

	/* auxiliary variables */
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	mesh->elem[El].mate = Mat_id;

	return 1;
}

/*---------------------------------------------------------
mmr_el_type_ref - to return element's type of refinement
---------------------------------------------------------*/
int mmr_el_type_ref( /* returns element's type of refinement, <0 - error code */
		/*         MMC_NOT_REF - not refined */
		/*         MMC_REF_ISO - isotropic refinement */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int El  	/* in: element ID */
)
{

	/* auxiliary variables */
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	return (mesh->elem[El].refi);
}

/*---------------------------------------------------------
mmr_el_faces - to get faces of an element
---------------------------------------------------------*/
int mmr_el_faces( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int El,     	/* in: element ID */
		int * Faces,  	/* out: list of faces */
		/*	(Faces[0] - number of faces) */
		int * Orient	/* out: orientation for each face */
		/*	+1 (MMC_SAME_ORIENT) - the same as element */
		/*	-1 (MMC_OPP_ORIENT) - opposite */
)
{
	/* auxiliary variables */
	mmt_mesh * mesh;
	int i, el_type, num_faces = 0;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	el_type = mmr_el_type(Mesh_id, El);
	if (el_type == MMC_BRICK) num_faces = 6;
	else if (el_type == MMC_PRISM) num_faces = 5;
	else if (el_type == MMC_TETRA) num_faces = 4;

	Faces[0] = num_faces;
	for (i = 1; i <= num_faces; i++) {
		Faces[i] = abs(mesh->elem[El].face[i - 1]);
		if (Orient != NULL) {
			if (mesh->elem[El].face[i - 1] > 0) Orient[i] = MMC_SAME_ORIENT;
			else Orient[i] = MMC_OPP_ORIENT;
		}
	}

	return 1;
}

/*---------------------------------------------------------
mmr_el_eq_neig - to get equal size neighbors of an element
---------------------------------------------------------*/
int mmr_el_eq_neig( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int El,     	/* in: element ID */
		int * Neig,   	/* out: list of equal size neighbors */
		/*  >0 - equal size neighbor ID */
		/*  -1 (MMC_BIG_NGB) - big neighbor */
		/*   0 (MMC_BOUNDARY) - boundary (always second neighbor)*/
		int * Neig_sides /* out: list of sides of neighbors */
)
{

	/* auxiliary variables */
	int i, el_faces[MMC_MAXELFAC + 1] ;
	int face_neig[2]; 	/* face's neighbors and their */
	int fa_neig_sides[2];	/* side numbers */

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* get element faces */
	mmr_el_faces(Mesh_id, El, el_faces, NULL);

	if (Neig != NULL) {
		/* for each face */
		Neig[0] = el_faces[0];
		for (i = 1; i <= Neig[0]; i++) {
			/* get its equal size neighbors */
			mmr_fa_eq_neig(Mesh_id, el_faces[i], face_neig, fa_neig_sides, NULL);
			if (face_neig[0] == El) {
				Neig[i] = face_neig[1];
				if (Neig_sides != NULL) Neig_sides[i] = fa_neig_sides[1];
			}
			else {
				Neig[i] = face_neig[0];
				if (Neig_sides != NULL) Neig_sides[i] = fa_neig_sides[0];
			}
		}
	}

	return 1;
}


/*---------------------------------------------------------
mmr_el_node_coor - to get the coordinates of element's nodes
---------------------------------------------------------*/
int mmr_el_node_coor( 	/* returns number of nodes */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int El,  	  /* in: element ID */
		int * Nodes,  	 /* out (optional Nodes==NULL switches off the option): */
		/*    list of vertex node IDs */
		/*	(Nodes[0] - number of nodes) */
		double * Node_coor /* out (optional Node_ccor==NULL switches off the option): */
		/*    coordinates of element vertices */
)
{

	/* auxiliary variables */
	int i, num_node = 0, el_type, node_shift;
	double fa_node_coor[3 * MMC_MAXFAVNO];
	int fa_node[MMC_MAXFAVNO + 1];
	int el_faces[MMC_MAXELFAC + 1], face_orient[MMC_MAXELFAC + 1];
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	el_type = mmr_el_type(Mesh_id, El);
	if (el_type == MMC_BRICK) {
		num_node = 8;
	}
	else if (el_type == MMC_PRISM) {
		num_node = 6;
	}
	else if (el_type == MMC_TETRA) {
		num_node = 4;
	}

	if (Nodes == NULL && Node_coor == NULL) return (num_node);

	// the only considered in fact
	if (el_type == MMC_PRISM) {

		/* take just two bases */
		mmr_el_faces(Mesh_id, El, el_faces, face_orient);

		/*kbw
		printf("element %d, \nfaces:",El);
		for(i=0;i<el_faces[0];i++) printf("  %d",el_faces[i+1]);
		printf("\n");
		kew*/

		/* for the first base */
		/* get nodes, their coordinates and possible shift */
		if (Node_coor != NULL) {
			mmr_fa_node_coor(Mesh_id, el_faces[1], fa_node, fa_node_coor);
		}
		else {
			mmr_fa_node_coor(Mesh_id, el_faces[1], fa_node, NULL);
		}

		/*kbw
		printf("element %d, face 1 - global %d, orient %d\nnodes:",El,el_faces[1]
		,face_orient[1]);
		for(i=0;i<fa_node[0];i++) printf("  %d",fa_node[i+1]);
		printf("\n");
		kew*/


		/* if orientation the same - shift==0 */
		if (face_orient[1] == MMC_SAME_ORIENT) {
			if (Nodes != NULL) {
				Nodes[1] = fa_node[1];
				Nodes[2] = fa_node[3];
				Nodes[3] = fa_node[2];
			}
			if (Node_coor != NULL) {
				for (i = 0; i < 3; i++) {
					Node_coor[i] = fa_node_coor[i];
					Node_coor[3 + i] = fa_node_coor[6 + i];
					Node_coor[6 + i] = fa_node_coor[3 + i];
				}
			}
		}
		else {
			/* different orientation - get shift and substitute */
			node_shift = -(mesh->face[el_faces[1]].bc);
			if (Nodes != NULL) {
				Nodes[1] = fa_node[node_shift + 1];
				Nodes[2] = fa_node[(node_shift + 1) % 3 + 1];
				Nodes[3] = fa_node[(node_shift + 2) % 3 + 1];
			}
			if (Node_coor != NULL) {
				for (i = 0; i < 3; i++) {
					Node_coor[i] = fa_node_coor[3 * node_shift + i];
					Node_coor[3 + i] = fa_node_coor[3 * ((node_shift + 1) % 3) + i];
					Node_coor[6 + i] = fa_node_coor[3 * ((node_shift + 2) % 3) + i];
				}
			}
		}

		/* for the second base */
		/* get nodes, their coordinates and possible shift */
		if (Node_coor != NULL) {
			mmr_fa_node_coor(Mesh_id, el_faces[2], fa_node, fa_node_coor);
		}
		else {
			mmr_fa_node_coor(Mesh_id, el_faces[2], fa_node, NULL);
		}

		/*kbw
		printf("element %d, face 2 - global %d, orient %d\nnodes:",El,el_faces[2]
		,face_orient[2]);
		for(i=0;i<fa_node[0];i++) printf("  %d",fa_node[i+1]);
		printf("\n");
		kew*/


		/* if orientation the same - shift==0 */
		if (face_orient[2] == MMC_SAME_ORIENT) {
			if (Nodes != NULL) {
				Nodes[4] = fa_node[1];
				Nodes[5] = fa_node[2];
				Nodes[6] = fa_node[3];
			}
			if (Node_coor != NULL) {
				for (i = 0; i < 3; i++) {
					Node_coor[9 + i] = fa_node_coor[i];
					Node_coor[12 + i] = fa_node_coor[3 + i];
					Node_coor[15 + i] = fa_node_coor[6 + i];
				}
			}
		}
		else {
			/* different orientation - get shift and substitute */
			node_shift = -mesh->face[el_faces[2]].bc;
			if (Nodes != NULL) {
				Nodes[4] = fa_node[node_shift + 1];
				Nodes[6] = fa_node[(node_shift + 1) % 3 + 1];
				Nodes[5] = fa_node[(node_shift + 2) % 3 + 1];
			}
			if (Node_coor != NULL) {
				for (i = 0; i < 3; i++) {
					Node_coor[9 + i] = fa_node_coor[3 * node_shift + i];
					Node_coor[15 + i] = fa_node_coor[3 * ((node_shift + 1) % 3) + i];
					Node_coor[12 + i] = fa_node_coor[3 * ((node_shift + 2) % 3) + i];
				}
			}
		}

		/* checking */
#ifdef DEBUG_MMM
		if (Nodes != NULL) {
			iaux = 0;
			/* for three side faces */
			for (ifa = 0; ifa < 3; ifa++) {

				/* get nodes, their coordinates and possible shift */
				mmr_fa_node_coor(Mesh_id, el_faces[3 + ifa], fa_node, fa_node_coor);

				/*kbw
				printf("element %d, face %d - global %d, orient %d\nnodes:",El,3+ifa,
				el_faces[3+ifa],face_orient[3+ifa]);
				for(i=0;i<fa_node[0];i++) printf("  %d",fa_node[i+1]);
				printf("\n");
				kew*/

				if (ifa == 0) {
					if (face_orient[3 + ifa] == MMC_SAME_ORIENT) {
						if (fa_node[1] != Nodes[1]) iaux = 1;
						if (fa_node[2] != Nodes[2]) iaux = 2;
						if (fa_node[3] != Nodes[5]) iaux = 5;
						if (fa_node[4] != Nodes[4]) iaux = 4;
					}
					else {
						node_shift = -mesh->face[el_faces[3 + ifa]].bc;
						if (node_shift != 1) {
							printf("Face_shift!=1 for vertical, quadrilateral face %d !\n", el_faces[3 + ifa]);
							printf("Use half structured meshes with shift==1 or modify\n");
							printf("data structure routines.\n");
							continue;
						}
						if (fa_node[1] != Nodes[2]) iaux = 2;
						if (fa_node[2] != Nodes[1]) iaux = 1;
						if (fa_node[3] != Nodes[4]) iaux = 4;
						if (fa_node[4] != Nodes[5]) iaux = 5;
					}
				}
				else if (ifa == 1) {
					if (face_orient[3 + ifa] == MMC_SAME_ORIENT) {
						if (fa_node[1] != Nodes[2]) iaux = 2;
						if (fa_node[2] != Nodes[3]) iaux = 3;
						if (fa_node[3] != Nodes[6]) iaux = 6;
						if (fa_node[4] != Nodes[5]) iaux = 5;
					}
					else {
						node_shift = -mesh->face[el_faces[3 + ifa]].bc;
						if (node_shift != 1) {
							printf("Face_shift!=1 for vertical, quadrilateral face %d !\n", el_faces[3 + ifa]);
							printf("Use half structured meshes with shift==1 or modify\n");
							printf("data structure routines.\n");
							continue;
						}
						if (fa_node[1] != Nodes[3]) iaux = 3;
						if (fa_node[2] != Nodes[2]) iaux = 2;
						if (fa_node[3] != Nodes[5]) iaux = 5;
						if (fa_node[4] != Nodes[6]) iaux = 6;
					}
				}
				else if (ifa == 2) {
					if (face_orient[3 + ifa] == MMC_SAME_ORIENT) {
						if (fa_node[1] != Nodes[3]) iaux = 3;
						if (fa_node[2] != Nodes[1]) iaux = 1;
						if (fa_node[3] != Nodes[4]) iaux = 4;
						if (fa_node[4] != Nodes[6]) iaux = 6;
					}
					else {
						node_shift = -mesh->face[el_faces[3 + ifa]].bc;
						if (node_shift != 1) {
							printf("Face_shift!=1 for vertical, quadrilateral face %d !\n", el_faces[3 + ifa]);
							printf("Use half structured meshes with shift==1 or modify\n");
							printf("data structure routines.\n");
							continue;
						}
						if (fa_node[1] != Nodes[1]) iaux = 1;
						if (fa_node[2] != Nodes[3]) iaux = 3;
						if (fa_node[3] != Nodes[6]) iaux = 6;
						if (fa_node[4] != Nodes[4]) iaux = 4;
					}
				}
			}

			if (iaux) {
				printf("At least node %d wrong for element %d\n", iaux, El);
			}

		}
#endif

	}
	else if (el_type == MMC_TETRA) {
		num_node = 4;
	}


	if (Nodes != NULL) Nodes[0] = num_node;

	return (num_node);
}

/*---------------------------------------------------------
mmr_el_edges - to get the list of element's edges
---------------------------------------------------------*/
int mmr_el_edges( 	/* returns the number of edges or error code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int El,  	  /* in: element ID */
		int * Edges   	 /* out : list of edge IDs */
		/*	(Edges[0] - number of edges) */
)
{

	/* auxiliary variables */
	int ifa, nr_fa_edges, num_edges = 0, el_type;
	int fa_edges[5];
	int el_faces[7], face_orient[7];
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	el_type = mmr_el_type(Mesh_id, El);
	if (el_type == MMC_BRICK) {
		num_edges = 12;
	}
	else if (el_type == MMC_PRISM) {
		num_edges = 9;

		mmr_el_faces(Mesh_id, El, el_faces, face_orient);

		/*kbw
		printf("element %d, \nfaces:",El);
		for(i=0;i<el_faces[0];i++) printf("  %d",el_faces[i+1]);
		printf("\n");
		kew*/

		/* we consider only vertical quadrilateral faces */
		for (ifa = 3; ifa <= el_faces[0]; ifa++) {

			nr_fa_edges = mmr_fa_edges(Mesh_id, el_faces[ifa], fa_edges, NULL);

			/*kbw
			printf("element %d, face %d - global %d, orient %d, nr_edges %d\nedges:",
				   El,ifa,el_faces[ifa],face_orient[ifa], nr_fa_edges);
			for(i=0;i<nr_fa_edges;i++) printf("  %d",fa_edges[i]);
			printf("\n");
			kew*/

			/* we take always the first (bottom) edge and the third (top) edge */
			/* as horizontal edges */
			Edges[ifa - 2] = fa_edges[0];
			Edges[ifa + 1] = fa_edges[2];

			/* depending on the orientation we take the second or the fourth edge */
			/* as vertical edges */
			if (face_orient[ifa] == MMC_SAME_ORIENT) {
				Edges[ifa + 4] = fa_edges[3];
			}
			else {
				Edges[ifa + 4] = fa_edges[1];
			}

			/*kbw
				  printf("element edges: %d, %d, %d\n",
					 Edges[ifa-2],Edges[ifa+1],Edges[ifa+4]);
			kew*/


		}

	}
	else if (el_type == MMC_TETRA) {
		num_edges = 6;
	}

	Edges[0] = num_edges;

	/*kbw
	printf("element %d, \nedges:",El);
	for(i=0;i<Edges[0];i++) printf("  %d",Edges[i+1]);
	printf("\n");
	kew*/


	return (num_edges);
}

/*---------------------------------------------------------
mmr_el_fam - to return family information for an element
---------------------------------------------------------*/
int mmr_el_fam( /* returns: element father ID or 0 (MMC_NO_FATH) */
		/*         for initial mesh element or <0 - error code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int El,       /* in: element ID (1...) */
		int * Elsons, /* out: list of element's sons */
		/* 	Elsons[0] - number of sons */
		int * Type	/* out: type of refinement */
		/*         MMC_NOT_REF - not refined */
		/*         MMC_REF_ISO - isotropic refinement */
)
{

	/* auxiliary variables */
	int i, nrelsons, el_type, ref_type;
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	if (Elsons != NULL) {

		if (mesh->elem[El].sons == NULL) {
			Elsons[0] = 0;
			ref_type = MMC_NOT_REF;
		}
		else {

			/* initialize */
			nrelsons = 0;

			/* ref_type = refinement type for BRICK and PRISM */
			ref_type = mmr_el_type_ref(Mesh_id, El);

			el_type = mmr_el_type(Mesh_id, El);
			if (el_type == MMC_BRICK) {
			}
			else if (el_type == MMC_PRISM) {

				if (ref_type == MMC_REF_ISO) nrelsons = 8; /* isotropic refinement */

			}
			else if (el_type == MMC_TETRA) {
			}

			Elsons[0] = nrelsons;


			/* substitute list of sons */
			for (i = 1; i <= nrelsons; i++)
				Elsons[i] = mesh->elem[El].sons[i - 1];

		} /* end if sons specified */
	} /* end if sons required */
	else if (Type != NULL) {
		/* if type of refinement required only! */
		if (mesh->elem[El].sons == NULL) ref_type = MMC_NOT_REF;
		else ref_type = mmr_el_type_ref(Mesh_id, El);
	}

	if (Type != NULL) *Type = ref_type;

	/* return father's number */
	return (mesh->elem[El].fath);
}

//---------------------------------------------------------
// mmr_el_fam_all - to return all recursive family information for an element
//---------------------------------------------------------
extern int mmr_el_fam_all( /* returns: element father ID or 0 (MMC_NO_FATH) for */
		/*          initial mesh element or <0 - error code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int El,       /* in: element ID (1...) */
		int * Elsons  /* out: list of all recursive element's sons */
		/* 	Elsons[0] - number of sons */
)
{

	mmt_mesh * mesh = mmr_select_mesh(Mesh_id);

	if (Elsons != NULL) {

		int sons[MMC_MAXELSONS + 1] = {0}, ison = 1, i = 1;
		const int n_grandsons_size = MMC_MAXELSONS * (mmr_get_max_gen(Mesh_id) - 1) + 1;
		int * grandsons = (int *) calloc( n_grandsons_size, sizeof(int));

		if (grandsons == NULL) {
			printf("Error: bad_calloc");
			exit(-1);
		}

		mmr_el_fam( Mesh_id, El, sons, NULL);
		Elsons[0] = sons[0];
		memcpy(Elsons + 1, sons + 1, sons[0]*sizeof(int));
		i += sons[0];

		for (; ison <= sons[0]; ++ison) {
			mmr_el_fam_all( Mesh_id, sons[ison], grandsons );
			Elsons[0] += grandsons[0];
			memcpy(Elsons + i, grandsons + 1, grandsons[0]*sizeof(int));
		}

		if (grandsons[0] > n_grandsons_size) {
			printf("Error: Array exceedeed!");
			exit(-1);
		}

		if (grandsons != NULL) {
			free(grandsons);
			grandsons = NULL;
		}
	}
	return (mesh->elem[El].fath);
}


/*---------------------------------------------------------
mmr_el_gen - to return generation level for an element
---------------------------------------------------------*/
int mmr_el_gen(   /* returns: El's generation ID */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int El        /* in: element ID (1...) */
)
{

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	if (mmr_el_fam(Mesh_id, El, NULL, NULL) == MMC_NO_FATH) return (0);
	else return (1 + mmr_el_gen(Mesh_id, mmr_el_fam(Mesh_id, El, NULL, NULL)));

}

/*---------------------------------------------------------
  mmr_el_ancestor - to find the ancestor of an element
					with generation level Ilev
---------------------------------------------------------*/
int mmr_el_ancestor( /* returns: >0 - ancestor ID, <0 - error code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int El,       /* in: element ID (1...) */
		int Ilev      /* in: level ID */
)
{

	int ancel;

	if (mmr_el_gen(Mesh_id, El) <= Ilev) return (El);
	else {
		ancel = mmr_el_ancestor(Mesh_id, mmr_el_fam(Mesh_id, El, NULL, NULL), Ilev);
	}

	return (ancel);
}

/*---------------------------------------------------------
mmr_el_hsize - to compute a characteristic linear size for an element
	(for linear and multi-linear 3D elements)
---------------------------------------------------------*/
double mmr_el_hsize(   /* returns: element size */
		int Mesh_id,	  /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int El,         /* in: element ID (1...) */
		double * Size_x, /* out: for anizotropic elements */
		double * Size_y,
		double * Size_z
)
{

	/* auxiliary variables */
	double hsize, node_coor[3 * MMC_MAXELVNO];
	int el_type, nodes[MMC_MAXELVNO + 1];
	double v1[3], v2[3], v3[3];
	int i;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* find element type */
	el_type = mmr_el_type(Mesh_id, El);

	/* get the coordinates of element nodes in the right order */
	mmr_el_node_coor(Mesh_id, El, nodes, node_coor);

	if (el_type == MMC_BRICK) {
	}
	else if (el_type == MMC_PRISM) {

		/* form three vectors from nodes coordinates */
		for (i = 0; i < 3; i++) {
			v1[i] = node_coor[3 + i] - node_coor[i];
			v2[i] = node_coor[6 + i] - node_coor[i];
			v3[i] = node_coor[9 + i] - node_coor[i];
		}

		/* compute hsize as third root of volume (computed as mixed product) */
		hsize = pow(fabs(mmr_vec3_mxpr(v1, v2, v3)), 1.0 / 3.0);

	}
	else if (el_type == MMC_TETRA) {
	}

	return (hsize);
}

/*---------------------------------------------------------
mmr_face_structure - to return face structure (e.g. for sending)
---------------------------------------------------------*/
int mmr_face_structure( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id,	   /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Fa,  	   /* in: face ID */
		int * Face_struct /* out: face structure in the form of an integer array */
)
{

	/* auxiliary variables */
	mmt_mesh * mesh;
	int i, nr_edge;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	Face_struct[0] = mesh->face[Fa].type;
	//Face_struct[1]=mesh->face[Fa].ipid;
	Face_struct[2] = mesh->face[Fa].bc;
	if (abs(Face_struct[0]) == MMC_QUAD) nr_edge = 4;
	else if (abs(Face_struct[0]) == MMC_TRIA) nr_edge = 3;
	for (i = 0; i < nr_edge; i++) Face_struct[3 + i] = mesh->face[Fa].edge[i];
	Face_struct[7] = mesh->face[Fa].neig[0];
	Face_struct[8] = mesh->face[Fa].neig[1];
	if (mesh->face[Fa].type < 0) {
		for (i = 0; i < 4; i++) {
			Face_struct[9 + i] = mesh->face[Fa].sons[i];
		}
	}

	return 1;

}

/*---------------------------------------------------------
mmr_fa_status - to return face status (active, inactive, free space)
---------------------------------------------------------*/
int mmr_fa_status( /* returns face status: */
		/*	 1 (MMC_ACTIVE)   - active face */
		/*	 0 (MMC_FREE)     - free space */
		/*	-1 (MMC_INACTIVE) - inactive face */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Fa  	/* in: face ID */
)
{

	/* auxiliary variables */
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	if (mesh->face[Fa].type > 0) return (MMC_ACTIVE);
	else if (mesh->face[Fa].type < 0) return (MMC_INACTIVE);
	else return (MMC_FREE);
}

/*---------------------------------------------------------
mmr_fa_type - to return face type (triangle, quad, free space)
---------------------------------------------------------*/
int mmr_fa_type( /* returns face type: */
		/*	 3 (MMC_TRIA) - triangle */
		/*	 4 (MMC_QUAD) - quadrilateral */
		/*	 0 (MMC_FREE) - free space */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Fa  	/* in: face ID */
)
{

	/* auxiliary variables */
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	return (abs(mesh->face[Fa].type));
}


/*--------------------------------------------------------------------------
 mmr_fa_bc - to get the boundary condition flag for a face
---------------------------------------------------------------------------*/
int mmr_fa_bc( /* returns: bc flag for a face */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Fa	      /* in: global face ID */
)
{
	/* auxiliary variables */
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	return (mesh->face[Fa].bc);

}


/*--------------------------------------------------------------------------
 mmr_fa_set_bc - to set the boundary condition flag for a face
---------------------------------------------------------------------------*/
int mmr_fa_set_bc( /* returns: bc flag for a face */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Fa,	      /* in: global face ID */
		int BC
)
{
	/* auxiliary variables */
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	mesh->face[Fa].bc = BC;

	return (mesh->face[Fa].bc);

}


/*---------------------------------------------------------
  mmr_fa_sub_bnd - to indicate face is on the intersubdomain boundary
---------------------------------------------------------*/
int mmr_fa_sub_bnd( /* returns: 1 - true, 0 - false */
		int Mesh_id,    /* in: mesh ID */
		int Face_id     /* in: face ID */
)
{
	/* auxiliary variables */
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	if (    mesh->face[Face_id].neig[0] == MMC_SUB_BND
			|| mesh->face[Face_id].neig[1] == MMC_SUB_BND ) return 1;
	else  return (0);
}

/*---------------------------------------------------------
  mmr_fa_set_sub_bnd - to set that the face is on the intersubdomain boundary
---------------------------------------------------------*/
int mmr_fa_set_sub_bnd( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id,    /* in: mesh ID */
		int Face_id,    /* in: face ID */
		int Side_id     /* in: side ID ??? */
)
{
	/* auxiliary variables */
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	mesh->face[Face_id].neig[Side_id] = MMC_SUB_BND;

	return 1;
}

/*---------------------------------------------------------
mmr_fa_edges - to return a list of face's edges
---------------------------------------------------------*/
int mmr_fa_edges( /* returns: >=0 - number of edges, <0 - error code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Fa,         /* in:  global face ID */
		int * Fa_edges, /* out: face edges IDs */
		int * Ed_orient /* out: edges orientation: */
		/*	+1 (MMC_SAME_ORIENT) - the same as face */
		/*	-1 (MMC_OPP_ORIENT)  - opposite */
)
{

	/* auxiliary variables */
	int i, fa_type, num_edges;
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	fa_type = mmr_fa_type(Mesh_id, Fa);

#ifdef DEBUG_MMM
	if (fa_type == MMC_FREE) {
		printf("no edges for free space\n");
		exit(1);
	}
#endif

	if (fa_type == MMC_TRIA) num_edges = 3;
	else num_edges = 4;

	/* rewrite edges from data structure */
	if (Fa_edges != NULL) {
		for (i = 0; i < num_edges; i++) {
			Fa_edges[i] = abs(mesh->face[Fa].edge[i]);
		}
		if (Ed_orient != NULL) {
			for (i = 0; i < num_edges; i++) {
				if (mesh->face[Fa].edge[i] > 0) Ed_orient[i] = MMC_SAME_ORIENT;
				else Ed_orient[i] = MMC_OPP_ORIENT;
			}
		}
	}

	/*kbw
	printf("In mmr_fa_edges: face %d\n", Fa);
	printf("edges:");
	for(i=0;i<num_edges;i++) printf("  %d", Ed_orient[i]*Fa_edges[i]);
	printf("\n");
	kew*/

	return (num_edges);
}

/*---------------------------------------------------------
mmr_fa_eq_neig - to return a list of face's neighbors and
		corresponding neighbors' sides (only equal size
		neighbors considered)
---------------------------------------------------------*/
void mmr_fa_eq_neig(
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Fa,       /* in:  global face ID */
		int * Fa_neig, /* out: face neighbors */
		/* (first neighbor has the same orientation */
		/*  as the face and ordering of its nodes */
		/*  defines ordering of nodes on the face) */
		/*  >0 - equal size neighbor ID */
		/*  -1 (MMC_BIG_NGB) - big neighbor */
		/*   0 (MMC_BOUNDARY) - boundary (always second neighbor)*/
		/*  -2 (MMC_SUB_BND) - inter-subdomain boundary */
		int * Neig_sides,     /* out: side local IDs for equal size neighbors */
		int * Node_shift	/* out: the difference in positions between */
		/*	first face's node for first and second */
		/* 	neighbor (usage in mmr_fa_node_coor) */
)
{

	/* auxiliary variables */
	int el_faces[MMC_MAXELFAC + 1];
	int iaux;
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

#ifdef DEBUG_MMM
	if (mesh->face[Fa].type == MMC_FREE) {
		printf("no neighbors for free space %d\n", Fa);
		exit(1);
	}
#endif

	/* rewrite equal size neighbors from data structure */
	if (Fa_neig != NULL) {
		Fa_neig[0] = mesh->face[Fa].neig[0];
		Fa_neig[1] = mesh->face[Fa].neig[1];
	}

	/* find which neighbor's side the face forms */
	if (Fa_neig != NULL && Neig_sides != NULL) {
		if (Fa_neig[0] > 0) {
			mmr_el_faces(Mesh_id, Fa_neig[0], el_faces, NULL);
			for (iaux = 0; iaux < el_faces[0]; iaux++) {
				if (el_faces[iaux + 1] == Fa) {
					Neig_sides[0] = iaux;
					break;
				}
			}
		}
		else Neig_sides[0] = 0;
		if (Fa_neig[1] > 0) {
			mmr_el_faces(Mesh_id, Fa_neig[1], el_faces, NULL);
			for (iaux = 0; iaux < el_faces[0]; iaux++) {
				if (el_faces[iaux + 1] == Fa) {
					Neig_sides[1] = iaux;
					break;
				}
			}
		}
		else Neig_sides[1] = 0;
	}

	if (Node_shift != NULL) {
		*Node_shift = -mesh->face[Fa].bc;
		if (*Node_shift < 0) *Node_shift = 0;
	}

	return;
}

/*---------------------------------------------------------
mmr_fa_neig - to return a list of face's neighbors and
		corresponding neighbors' sides (for active
		faces all neighbors are active)
---------------------------------------------------------*/
void mmr_fa_neig(
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Fa,               /* in:  global face ID */
		int * Fa_neig,        /* out: face neighbors */
		/* (first neighbor has the same orientation */
		/*  as the face and ordering of its nodes */
		/*  defines ordering of nodes on the face) */
		/*  >0 - equal size neighbor ID */
		/*  <0 - big neighbor ID */
		/*   0 (MMC_BOUNDARY) - boundary (always 2nd neighbor)*/
		int * Neig_sides,     /* out: side IDs for neighbors */
		int * Node_shift,     /* out: the difference in positions between */
		/*	first face's node for first and second */
		/* 	neighbor (usage in mmr_fa_node_coor) */
		int * Diff_gen,	/* out: generation difference between neighbors */
		double * Acoeff,      /* out: coefficients of linear transformation... */
		double * Bcoeff	/* 	between face coordinates and */
		/* 	coordinates on an ancestor face  */
		/* 	being a side of big neighbor  */
)
{

	/* auxiliary variables */
	int ineig, side, neig, iaux;
	int igen, ison, diff_gen, fa_type;
	int el_faces[MMC_MAXELFAC + 1];
	int side_fath, side_old, side_sons[5], node_mid;
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

#ifdef DEBUG_MMM
	if (mesh->face[Fa].type == MMC_FREE) {
		printf("no neighbors for free space %d in fa_neig\n", Fa);
		exit(-1);
	}
	if (Fa_neig == NULL) {
		printf("no space for array Fa_neig in fa_neig\n");
		exit(-1);
	}
	if (Fa == 0) {
		printf("asking for neighbors of %d in fa_neig\n", Fa);
		exit(-1);
	}
#endif

	/* find equal size neighbors and their sides */
	mmr_fa_eq_neig(Mesh_id, Fa, Fa_neig, Neig_sides, Node_shift);

	/* initialize generation difference for big neighbor case */
	diff_gen = 0;

	/* for each face side */
	for (ineig = 0; ineig < 2; ineig++) {

#ifdef DEBUG_MMM
		if (Fa_neig[ineig] == MMC_SUB_BND) {
			printf("looking for neighbor across intersubdomain boundary!!!\n");
			exit(-1);
		}
#endif

		/* for big neighbor case - find its number */
		side = Fa;
		neig = Fa_neig[ineig];
		while (neig == MMC_BIG_NGB) {
			diff_gen++;
			side = mmr_fa_fam(Mesh_id, side, NULL, NULL);
			neig = mesh->face[side].neig[ineig];
		}
		if (Fa_neig[ineig] == MMC_BIG_NGB) {
			Fa_neig[ineig] = -neig;
			if (Neig_sides != NULL) {
				mmr_el_faces(Mesh_id, neig, el_faces, NULL);
				for (iaux = 0; iaux < el_faces[0]; iaux++) {
					if (el_faces[iaux + 1] == side) {
						Neig_sides[ineig] = iaux;
						break;
					}
				}
			}
		}

	}/* loop over ineig */

	if (Diff_gen != NULL) *Diff_gen = diff_gen;

	/* in case of big neighbor find position in ancestor older by diff_gen */
	if (Acoeff != NULL && diff_gen > 0) {

		/*kbw
		if(Neig_sides!=NULL) printf("Face %d, neighbors %d (side %d), %d  (side %d)- generation difference %d\n",
		Fa, Fa_neig[0], Neig_sides[0], Fa_neig[1], Neig_sides[1], diff_gen);
		else printf("Face %d, neighbors %d, %d - generation difference %d\n",
		Fa, Fa_neig[0], Fa_neig[1], diff_gen);
		kew*/

		fa_type = mmr_fa_type(Mesh_id, Fa);
		side = Fa;
		side_fath = mmr_fa_fam(Mesh_id, side, NULL, NULL);

		/* to speed up computations for 1-irregular meshes */
		if (diff_gen == 1) {

			mmr_fa_fam(Mesh_id, side_fath, side_sons, &node_mid);

			ison = mmr_chk_list(side, &side_sons[1], side_sons[0]);
			if (fa_type == MMC_QUAD && node_mid > 0) {
				Acoeff[0] = 0.5;
				Acoeff[1] = 0.0;
				Acoeff[2] = 0.0;
				Acoeff[3] = 0.5;
				if (ison == 1) {
					Bcoeff[0] = -0.5;
					Bcoeff[1] = -0.5;
				}
				else if (ison == 2) {
					Bcoeff[0] = 0.5;
					Bcoeff[1] = -0.5;
				}
				else if (ison == 3) {
					Bcoeff[0] = 0.5;
					Bcoeff[1] = 0.5;
				}
				else if (ison == 4) {
					Bcoeff[0] = -0.5;
					Bcoeff[1] = 0.5;
				}

			}
			else if (fa_type == MMC_TRIA) {

				if (ison < 4) {

					Acoeff[0] = 0.5;
					Acoeff[1] = 0.0;
					Acoeff[2] = 0.0;
					Acoeff[3] = 0.5;

					if (ison == 1) {
						Bcoeff[0] = 0.0;
						Bcoeff[1] = 0.0;
					}
					else if (ison == 2) {
						Bcoeff[0] = 0.5;
						Bcoeff[1] = 0.0;
					}
					else if (ison == 3) {
						Bcoeff[0] = 0.0;
						Bcoeff[1] = 0.5;
					}
				}
				else {

					Acoeff[0] = -0.5;
					Acoeff[1] =  0.0;
					Acoeff[2] =  0.0;
					Acoeff[3] = -0.5;
					Bcoeff[0] =  0.5;
					Bcoeff[1] =  0.5;

				}
			}
		} /* end if one generation difference between faces */
		else {

			/* in a loop over generations */
			Acoeff[0] = 1.0;
			Acoeff[1] = 0.0;
			Acoeff[2] = 0.0;
			Acoeff[3] = 1.0;
			Bcoeff[0] = 0.0;
			Bcoeff[1] = 0.0;

			for (igen = 0; igen < diff_gen; igen++) {

				side_old = mmr_fa_fam(Mesh_id, side_fath, side_sons, &node_mid);

				/* check which son face is */
				ison = mmr_chk_list(side, &side_sons[1], side_sons[0]);

				/*kbw
				printf("Father %d, son %d - number %d,node_mid %d\n",
				side_fath,side,ison,node_mid);
				printf("Sons (%d): ",side_sons[0]);
				for(i=1;i<=side_sons[0];i++) printf(" %d",side_sons[i]);
				printf("\n");
				kew*/

				/* update coefficients A and B */

				/* find transformation coefficients to father element coordinates */
				if (fa_type == MMC_QUAD && node_mid > 0) {

					/* find coordinates within father face */
					Acoeff[0] = 0.5 * Acoeff[0];
					Acoeff[1] = 0.5 * Acoeff[1];
					Acoeff[2] = 0.5 * Acoeff[2];
					Acoeff[3] = 0.5 * Acoeff[3];
					if (ison == 1) {
						Bcoeff[0] = 0.5 * Bcoeff[0] - 0.5;
						Bcoeff[1] = 0.5 * Bcoeff[1] - 0.5;
					}
					else if (ison == 2) {
						Bcoeff[0] = 0.5 * Bcoeff[0] + 0.5;
						Bcoeff[1] = 0.5 * Bcoeff[1] - 0.5;
					}
					else if (ison == 3) {
						Bcoeff[0] = 0.5 * Bcoeff[0] + 0.5;
						Bcoeff[1] = 0.5 * Bcoeff[1] + 0.5;
					}
					else if (ison == 4) {
						Bcoeff[0] = 0.5 * Bcoeff[0] - 0.5;
						Bcoeff[1] = 0.5 * Bcoeff[1] + 0.5;
					}

				}
				else if (fa_type == MMC_TRIA) {

					if (ison < 4) {

						Acoeff[0] = 0.5 * Acoeff[0];
						Acoeff[1] = 0.5 * Acoeff[1];
						Acoeff[2] = 0.5 * Acoeff[2];
						Acoeff[3] = 0.5 * Acoeff[3];

						if (ison == 1) {
							Bcoeff[0] = 0.5 * Bcoeff[0];
							Bcoeff[1] = 0.5 * Bcoeff[1];
						}
						else if (ison == 2) {
							Bcoeff[0] = 0.5 * Bcoeff[0] + 0.5;
							Bcoeff[1] = 0.5 * Bcoeff[1];
						}
						else if (ison == 3) {
							Bcoeff[0] = 0.5 * Bcoeff[0];
							Bcoeff[1] = 0.5 * Bcoeff[1] + 0.5;
						}
					}
					else {

						Acoeff[0] = -0.5 * Acoeff[0];
						Acoeff[1] = -0.5 * Acoeff[1];
						Acoeff[2] = -0.5 * Acoeff[2];
						Acoeff[3] = -0.5 * Acoeff[3];
						Bcoeff[0] = -0.5 * Bcoeff[0] + 0.5;
						Bcoeff[1] = -0.5 * Bcoeff[1] + 0.5;

					}
				}
				/*kbw
				printf("Acoeff %lf %lf %lf %lf\n",
				Acoeff[0],Acoeff[1],Acoeff[2],Acoeff[3]);
				printf("Bcoeff %lf %lf\n",
				Bcoeff[0],Bcoeff[1]);
				kew*/

				side = side_fath;
				side_fath = side_old;

			} /* end loop over generations */

		} /* end if more than one generation difference between faces */

	} /* end if big neighbor */

}

/*---------------------------------------------------------
mmr_fa_node_coor - to get the list and coordinates of faces's nodes
---------------------------------------------------------*/
int mmr_fa_node_coor( /* returns: number of nodes for a face */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Fa,  	     /* in: face ID */
		int * Nodes,  	    /* out: list of vertex node IDs */
		double * Node_coor /* out: coordinates of face vertices */
)
{

	/* auxiliary variables */
	int i, ino, num_node, fa_type, ied, nodes[2];
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	/* find the type of the face */
	fa_type = mmr_fa_type(Mesh_id, Fa);
	if (fa_type == MMC_QUAD) num_node = 4;
	else if (fa_type == MMC_TRIA) num_node = 3;
	if (Nodes != NULL) Nodes[0] = num_node;

	/* in a loop over edges */
	for (i = 0; i < num_node; i++) {
		ied = mesh->face[Fa].edge[i];
		mmr_edge_nodes(Mesh_id, abs(ied), nodes);
		if (ied > 0) ino = nodes[0];
		else ino = nodes[1];
		if (Nodes != NULL) Nodes[i + 1] = ino;
		if (Node_coor != NULL) {
			Node_coor[3 * i] = mesh->node[ino].x;
			Node_coor[3 * i + 1] = mesh->node[ino].y;
			Node_coor[3 * i + 2] = mesh->node[ino].z;
		}

		/*checking*/
#ifdef DEBUG_MMM
		/*kbw
		printf("face %d, edge %d - global %d, chosen node - %d\n",
		Fa,i,ied,ino);
		kew*/
#endif

	}

	return (num_node);
}

/*---------------------------------------------------------
mmr_fa_elem_coor - to find coordinates within neighboring
		elements for a point on face
---------------------------------------------------------*/
void mmr_fa_elem_coor(
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		double * Xloc,	       /* in: local coordinates on a face */
		/* !!!!! for triangular faces coordinates are standard [0,1][0,1] */
		/* !!!!! for quadrilateral faces coordinates are [-1,1][-1,1] - */
		/* !!!!! which means that they do not conform to element coordinates */
		int * Fa_neig,	       /* in: face neighbors (<0 - big) */
		/*	first - same orientation */
		/*	second - opposite orientation */
		int * Neig_sides,     /* in: which side face is for neighbors */
		int Node_shift,       /* in: the difference in positions between */
		/*	first face's node for first and second */
		/* 	neighbor (usage in mmr_fa_node_coor) */
		double * Acoeff,      /* in: coefficients of transformation between ... */
		double * Bcoeff,	/* in: face coord and big neighb face coord */
		double * Xneig        /* out: local coordinates for neighbors */
)
{

	/* auxiliary variables */
	int neig, ineig, el_type;
	double coor[2];

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* loop over neighbors */
	for (ineig = 0; ineig < 2; ineig++) {

		if (Fa_neig[ineig] != 0) {

			neig = abs(Fa_neig[ineig]);
			el_type = mmr_el_type(Mesh_id, neig);

			/* equal size neighbor - coordinate is just pased */
			if (Fa_neig[ineig] > 0) {
				coor[0] = Xloc[0];
				coor[1] = Xloc[1];
			}
			else {
				/* transform coeficients */
				coor[0] = Acoeff[0] * Xloc[0] + Acoeff[1] * Xloc[1] + Bcoeff[0];
				coor[1] = Acoeff[2] * Xloc[0] + Acoeff[3] * Xloc[1] + Bcoeff[1];
			}

			if (el_type == MMC_BRICK) {
			}
			else if (el_type == MMC_PRISM) {

				if (Neig_sides[ineig] == 0) {

					Xneig[3 * ineig + 2] = -1.0;

					if (ineig == 0) {
						/* first neighbor with the same orientation */

						Xneig[3 * ineig]   = coor[1];
						Xneig[3 * ineig + 1] = coor[0];

					}
					else {
						/* second neighbor - possibly shifted and with the opposite orientation */

						if (Node_shift == 0) {
							Xneig[3 * ineig]   = coor[0];
							Xneig[3 * ineig + 1] = coor[1];
						}
						if (Node_shift == 1) {
							Xneig[3 * ineig]   = 1 - coor[0] - coor[1];
							Xneig[3 * ineig + 1] = coor[1];
						}
						if (Node_shift == 2) {
							Xneig[3 * ineig]   = coor[0];
							Xneig[3 * ineig + 1] = 1 - coor[0] - coor[1];
						}

					}

				}
				else if (Neig_sides[ineig] == 1) {

					Xneig[3 * ineig + 2] = 1.0;

					if (ineig == 0) {
						/* first neighbor with the same orientation */

						Xneig[3 * ineig]   = coor[0];
						Xneig[3 * ineig + 1] = coor[1];

					}
					else {
						/* second neighbor - possibly shifted and with the opposite orientation */

						if (Node_shift == 0) {
							Xneig[3 * ineig]   = coor[1];
							Xneig[3 * ineig + 1] = coor[0];
						}
						if (Node_shift == 1) {
							Xneig[3 * ineig]   = 1 - coor[0] - coor[1];
							Xneig[3 * ineig + 1] = coor[0];
						}
						if (Node_shift == 2) {
							Xneig[3 * ineig]   = coor[1];
							Xneig[3 * ineig + 1] = 1 - coor[0] - coor[1];
						}

					}
				}

				else if (Neig_sides[ineig] == 2) {

					Xneig[3 * ineig + 1] = 0.0;
					Xneig[3 * ineig + 2] = coor[1];

					if (ineig == 0) {
						/* first neighbor with the same orientation */

						Xneig[3 * ineig]   = 0.5 * (coor[0] + 1.0);

					}
					else {
						/* second neighbor - with the opposite orientation */

						Xneig[3 * ineig]   = 0.5 * (-coor[0] + 1.0);

					}
				}

				else if (Neig_sides[ineig] == 3) {

					Xneig[3 * ineig + 2] = coor[1];

					if (ineig == 0) {
						/* first neighbor with the same orientation */

						Xneig[3 * ineig]   = 0.5 * (-coor[0] + 1.0);;

					}
					else {
						/* second neighbor - with the opposite orientation */

						Xneig[3 * ineig]   = 0.5 * (coor[0] + 1.0);

					}

					Xneig[3 * ineig + 1] = 1.0 - Xneig[3 * ineig];

				}

				else if (Neig_sides[ineig] == 4) {

					Xneig[3 * ineig] = 0.0;
					Xneig[3 * ineig + 2] = coor[1];

					if (ineig == 0) {
						/* first neighbor with the same orientation */

						Xneig[3 * ineig + 1]   = 0.5 * (-coor[0] + 1.0);

					}
					else {
						/* second neighbor - with the opposite orientation */

						Xneig[3 * ineig + 1]   = 0.5 * (coor[0] + 1.0);

					}
				}

			}
			else if (el_type == MMC_TETRA) {
			}
		}
	}

	/*kbw
	#ifdef DEBUG_MMM
	printf("In elem coord for a point %lf, %lf on a face\n",
	coor[0],coor[1]);
	printf("Neigh1 %d, side %d\n", Fa_neig[0], Neig_sides[0]);
	printf("Coor: %lf, %lf, %lf\n", Xneig[0], Xneig[1], Xneig[2]);
	if(Fa_neig[1]!=0){
	printf("Neigh2 %d, side %d, node_shift %d\n",
	Fa_neig[1], Neig_sides[1], Node_shift);
	printf("Coor: %lf, %lf, %lf\n", Xneig[3], Xneig[4], Xneig[5]);
	}
	#endif
	kew*/


	return;
}

/*---------------------------------------------------------
mmr_fa_area - to compute the area of face and vector normal
---------------------------------------------------------*/
void mmr_fa_area(
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Fa,		/* in: global face ID */
		double * Area,	       /* out: face area */
		double * Vec_norm     /* out: normal vector */
)
{

	/* auxiliary variables */
	double node_coor[3 * MMC_MAXFAVNO], daux;
	double vec_a[3], vec_b[3], vec_c[3], vec_d[3], vec_e[3], vec_f[3];
	int i, fa_type, nodes[MMC_MAXFAVNO + 1];

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* find face type */
	fa_type = mmr_fa_type(Mesh_id, Fa);

	/* get the coordinates of face nodes in the right order */
	mmr_fa_node_coor(Mesh_id, Fa, nodes, node_coor);


	if (fa_type == MMC_QUAD) {

		for (i = 0; i < 3; i++) {
			vec_a[i] = node_coor[3 + i] - node_coor[i];
			vec_b[i] = node_coor[6 + i] - node_coor[3 + i];
			vec_c[i] = node_coor[9 + i] - node_coor[6 + i];
			vec_d[i] = node_coor[i] - node_coor[9 + i];
		}
		mmr_vec3_prod(vec_a, vec_b, vec_e);
		mmr_vec3_prod(vec_c, vec_d, vec_f);
		daux = mmr_vec3_length(vec_e);
		if (Area != NULL) *Area = 0.5 * ( daux + mmr_vec3_length(vec_f) );
		/* Vec_norm valid only if faces are plane by assumption */
		if (Vec_norm != NULL) {
			for (i = 0; i < 3; i++)
				Vec_norm[i] = vec_e[i] / daux;

#ifdef DEBUG_MMM
			for (i = 0; i < 3; i++) {
				if (fabs(Vec_norm[i] - vec_f[i] / mmr_vec3_length(vec_f)) > 1e-9) {
					printf("Error 3275 in mmr_fa_area %.12lf != %.12lf\n",
							Vec_norm[i], vec_f[i] / mmr_vec3_length(vec_f));
				}
			}
#endif

		}
	}
	else if (fa_type == MMC_TRIA) {

		for (i = 0; i < 3; i++) {
			vec_a[i] = node_coor[3 + i] - node_coor[i];
			vec_b[i] = node_coor[6 + i] - node_coor[3 + i];
		}
		mmr_vec3_prod(vec_a, vec_b, vec_e);
		daux = mmr_vec3_length(vec_e);
		if (Area != NULL) *Area = 0.5 * daux;
		if (Vec_norm != NULL)
			for (i = 0; i < 3; i++)
				Vec_norm[i] = vec_e[i] / daux;

	}

	return;
}


/*--------------------------------------------------------------------------
mmr_fa_fam - to return face's family information
---------------------------------------------------------------------------*/
int mmr_fa_fam( /* returns: face's father ID or 0 (MMC_NO_FATH) if there is */
		/*          no father or <0 for error code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Fa,      	/* in: global face ID */
		int * Fasons,	/* out: sons */
		/* 	Fasons[0] - number of sons */
		int * Node_mid	/* out: node in the middle (if any) */
)
{

	/* auxiliary variables */
	int i, iaux, iside, ineig, neig, ref_kind;
	int el_fath, fa_fath, neig_type, neig_sons[MMC_MAXELSONS + 1];
	int fath_sons[MMC_MAXELSONS + 1], nodes[2];
	int node_shift, face_neig[2], neig_sides[2];
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

#ifdef DEBUG_MMM
	/* check whether face is active */
	if (mmr_fa_status(Mesh_id, Fa) == MMC_FREE) {
		printf("No father for free space %d!\n", Fa);
		getchar();
		getchar();
		getchar();
		return -1;
	}
	/* check consistency */
	if (Fasons == NULL && Node_mid != NULL) {
		printf("Middle node requested in fa_fam with no space for sons!\n");
		return -1;
	}
#endif

	/* new data structure */
	if (Fasons != NULL && mmr_fa_status(Mesh_id, Fa) == MMC_INACTIVE) {
		Fasons[0] = 4;
		for (i = 1; i <= 4; i++) {
			Fasons[i] = mesh->face[Fa].sons[i - 1];
		}
	}
	if (Node_mid != NULL && mmr_fa_status(Mesh_id, Fa) == MMC_INACTIVE) {
		if (mmr_fa_type(Mesh_id, Fa) == MMC_TRIA) {
			*Node_mid = 0;
		}
		else {
			i = mesh->face[Fasons[1]].edge[1];
			mmr_edge_nodes(Mesh_id, abs(i), nodes);
			if (i > 0) *Node_mid = nodes[1];
			else *Node_mid = nodes[0];
		}
	}



	/* find equal size face neighbors */
	mmr_fa_eq_neig(Mesh_id, Fa, face_neig, neig_sides, &node_shift);

	/*kbw
	  if(Fa==1980||Fa==1968){
		printf("face %d, neighbors %d %d, sides %d %d, shift %d\n",
		   Fa, face_neig[0], face_neig[1], neig_sides[0],
		   neig_sides[1], node_shift);
	  }
	kew*/

	/* choose neighbor - better refined one, since we may need its sons as well */
	if (face_neig[0] <= 0) {

		if (face_neig[1] <= 0) {

#ifdef DEBUG_MMM
			printf("empty face on intersubdomain boundary - unknown father\n");
			if (face_neig[1] == face_neig[0]) {
				printf("both neighbors indicated %d - error exiting\n", face_neig[1]);
				exit(-1);
			}
#endif

			return (0);
		}
		else {
			ineig = 1;
		}
	}
	else {
		ineig = 0;
		if (face_neig[1] > 0 && mmr_el_status(Mesh_id, face_neig[0]) != MMC_INACTIVE) {
			ineig = 1;
		}
	}
	neig = face_neig[ineig];

	/* check the element type of a neighbor */
	neig_type = mmr_el_type(Mesh_id, neig);

	/* find side */
	iside = neig_sides[ineig];

	/* find its father and sons */
	el_fath = mmr_el_fam(Mesh_id, neig, NULL, NULL);

	/* substitute proper face */
	if (el_fath == MMC_NO_FATH) fa_fath = MMC_NO_FATH;
	else {

		/* if face has a big neighbor or is on the boundary - situation is clear */
		if (face_neig[0] <= 0 || face_neig[1] <= 0) {
			if (el_fath > 0) fa_fath = abs(mesh->elem[el_fath].face[iside]);
			else fa_fath = MMC_NO_FATH;
		}
		else {
			/* otherwise we have to check whether it has father at all */

			/* check which son neig is */
			mmr_el_fam(Mesh_id, el_fath, fath_sons, &i);

			if (neig_type == MMC_PRISM && i == MMC_REF_ISO) {

				iaux = mmr_chk_list(neig, &fath_sons[1], fath_sons[0]);

				if ((iaux == 4 || iaux == 8) && iside >= 2) fa_fath = MMC_NO_FATH;
				else if ((iaux <= 4 && iside == 1) || (iaux > 4 && iside == 0)) fa_fath = MMC_NO_FATH;
				else if ((iaux <= 4 && iside == iaux % 3 + 2) ||
						(iaux > 4 && iside == (iaux - 4) % 3 + 2)) fa_fath = MMC_NO_FATH;
				else fa_fath = abs(mesh->elem[el_fath].face[iside]);

				/*kbw
					if(Fa==1980||Fa==1968){
					  printf("el_fath sons:");
					  for(i=0;i<fath_sons[0];i++) printf("  %d",fath_sons[i+1]);
					  printf("\n");
					  printf("neig %d is %dth son and side is %d - el_fath %d\n",
						 neig,iaux,iside,fa_fath);
					}
				kew*/

			}
		}
	}

	/*kbw
	  if(Fa==1980||Fa==1968){
		printf("neighbor %d - element %d, side %d, father %d, face_father %d\n",
		   ineig, neig, iside, el_fath, fa_fath);
		if(Fasons!=NULL){
		  printf("neig sons:");
		  for(i=0;i<neig_sons[0];i++) printf("  %d",neig_sons[i+1]);
		  printf("\n");
		}
	  }
	kew*/

#ifdef DEBUG_MMM
	if (face_neig[0] > 0 && face_neig[1] > 0) {
		/*simple check*/
		if ((mmr_el_fam(Mesh_id, face_neig[0], NULL, NULL) ==
						mmr_el_fam(Mesh_id, face_neig[1], NULL, NULL)) && fa_fath != MMC_NO_FATH) {
			printf("error in fa_fam, internal face %d has father!\n", Fa);
			exit(-1);
		}
	}
#endif

	return (fa_fath);



	/************** OLD ALGORITHM *********************************/

#ifdef DEBUG_MMM

	/* if sons required */
	if (Fasons != NULL) {

		el_fath = mmr_el_fam(Mesh_id, neig, neig_sons, &ref_kind);
		/* check whether any of neighboring elements has been refined */
		if (mmr_el_status(Mesh_id, neig) != MMC_INACTIVE) {
			printf("no sons for face %d in fa_faces\n", Fa);
			exit(-1);
		}
		else {

			/* for isotropic refinements of prisms */
			if (neig_type == MMC_PRISM && ref_kind == MMC_REF_ISO) {

				/* substitute proper faces - orientation!*/
				if (ineig == 0) {
					/* first neighbor - same orientation */
					if (iside == 0) {
						Fasons[0] = 4;
						Fasons[1] = mesh->elem[neig_sons[1]].face[iside];
						Fasons[2] = mesh->elem[neig_sons[3]].face[iside];
						Fasons[3] = mesh->elem[neig_sons[2]].face[iside];
						Fasons[4] = mesh->elem[neig_sons[4]].face[iside];
					}
					else if (iside == 1) {
						Fasons[0] = 4;
						Fasons[1] = mesh->elem[neig_sons[5]].face[iside];
						Fasons[2] = mesh->elem[neig_sons[6]].face[iside];
						Fasons[3] = mesh->elem[neig_sons[7]].face[iside];
						Fasons[4] = mesh->elem[neig_sons[8]].face[iside];
					}
					else if (iside == 2) {
						Fasons[0] = 4;
						Fasons[1] = mesh->elem[neig_sons[1]].face[iside];
						Fasons[2] = mesh->elem[neig_sons[2]].face[iside];
						Fasons[3] = mesh->elem[neig_sons[6]].face[iside];
						Fasons[4] = mesh->elem[neig_sons[5]].face[iside];
					}
					else if (iside == 3) {
						Fasons[0] = 4;
						Fasons[1] = mesh->elem[neig_sons[2]].face[iside];
						Fasons[2] = mesh->elem[neig_sons[3]].face[iside];
						Fasons[3] = mesh->elem[neig_sons[7]].face[iside];
						Fasons[4] = mesh->elem[neig_sons[6]].face[iside];
					}
					else if (iside == 4) {
						Fasons[0] = 4;
						Fasons[1] = mesh->elem[neig_sons[3]].face[iside];
						Fasons[2] = mesh->elem[neig_sons[1]].face[iside];
						Fasons[3] = mesh->elem[neig_sons[5]].face[iside];
						Fasons[4] = mesh->elem[neig_sons[7]].face[iside];
					}
				} /* end if first neighbor */
				else {
					/* second neighbor - opposite orientation */
					if (iside == 0) {
						Fasons[0] = 4;
						Fasons[node_shift + 1] = -mesh->elem[neig_sons[1]].face[iside];
						Fasons[(node_shift + 1) % 3 + 1] = -mesh->elem[neig_sons[2]].face[iside];
						Fasons[(node_shift + 2) % 3 + 1] = -mesh->elem[neig_sons[3]].face[iside];
						Fasons[4] = -mesh->elem[neig_sons[4]].face[iside];
					}
					else if (iside == 1) {
						Fasons[0] = 4;
						Fasons[node_shift + 1] = -mesh->elem[neig_sons[5]].face[iside];
						Fasons[(node_shift + 1) % 3 + 1] = -mesh->elem[neig_sons[7]].face[iside];
						Fasons[(node_shift + 2) % 3 + 1] = -mesh->elem[neig_sons[6]].face[iside];
						Fasons[4] = -mesh->elem[neig_sons[8]].face[iside];
					}
					else if (iside == 2) {
						if (node_shift == 1) {
							Fasons[0] = 4;
							Fasons[1] = -mesh->elem[neig_sons[2]].face[iside];
							Fasons[2] = -mesh->elem[neig_sons[1]].face[iside];
							Fasons[3] = -mesh->elem[neig_sons[5]].face[iside];
							Fasons[4] = -mesh->elem[neig_sons[6]].face[iside];
						}
						else {
							printf("Face_shift!=1 for vertical, quadrilateral face %d !\n", neig);
							printf("Use half structured meshes with shift==1 or modify\n");
							printf("data strucures routines.\n");
							return -1;
						}
					}
					else if (iside == 3) {
						if (node_shift == 1) {
							Fasons[0] = 4;
							Fasons[1] = -mesh->elem[neig_sons[3]].face[iside];
							Fasons[2] = -mesh->elem[neig_sons[2]].face[iside];
							Fasons[3] = -mesh->elem[neig_sons[6]].face[iside];
							Fasons[4] = -mesh->elem[neig_sons[7]].face[iside];
						}
						else {
							printf("Face_shift!=1 for vertical, quadrilateral face %d !\n", neig);
							printf("Use half structured meshes with shift==1 or modify\n");
							printf("data strucures routines.\n");
							return -1;
						}
					}
					else if (iside == 4) {
						if (node_shift == 1) {
							Fasons[0] = 4;
							Fasons[1] = -mesh->elem[neig_sons[1]].face[iside];
							Fasons[2] = -mesh->elem[neig_sons[3]].face[iside];
							Fasons[3] = -mesh->elem[neig_sons[7]].face[iside];
							Fasons[4] = -mesh->elem[neig_sons[5]].face[iside];
						}
						else {
							printf("Face_shift!=1 for vertical, quadrilateral face %d !\n", neig);
							printf("Use half structured meshes with shift==1 or modify\n");
							printf("data strucures routines.\n");
							return -1;
						}
					}

				}

				/* substitute middle node (if exists and is requested) */
				if (Node_mid != NULL) {
					if (iside < 2) {
						*Node_mid = 0;
					}
					else {
						if (ineig == 0) {
							/* first neighbor - same orientation */
							i = mesh->face[Fasons[1]].edge[1];
							mmr_edge_nodes(Mesh_id, abs(i), nodes);
							if (i > 0) *Node_mid = nodes[1];
							else *Node_mid = nodes[0];
						}
						else {
							i = mesh->face[Fasons[1]].edge[1];
							mmr_edge_nodes(Mesh_id, abs(i), nodes);
							if (i > 0) *Node_mid = nodes[1];
							else *Node_mid = nodes[0];
						}
					}
				}


			} /* end if PRISM and ISO refinement */
			else {
				printf("mmr_fa_fam not implemented for this element and refinement kind !\n");
				return -1;
			}

		}
	}

	return (fa_fath);

#endif

	/************** END OF OLD ALGORITHM *********************************/

}


/*---------------------------------------------------------
mmr_edge_nodes - to return edge node's numbers
---------------------------------------------------------*/
int mmr_edge_nodes( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Ed,		/* in: edge ID */
		int * Edge_nodes	/* out: IDs of edge nodes */
)
{

	/* auxiliary variables */
	int iaux, nodes[2];
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* if(Ed==15443){ */
	/*   printf("Edge %d\n",Ed); */
	/*   getchar(); */
	/* } */


	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	/* if edge is active just substitute its nodes */
	if (mesh->edge[Ed].type > 0) {
		Edge_nodes[0] = mesh->edge[Ed].node[0];
		Edge_nodes[1] = mesh->edge[Ed].node[1];
	}
	else if (mesh->edge[Ed].type < 0) {

		/* collect nodes from two edge sons */
		iaux = mmr_edge_nodes(Mesh_id, mesh->edge[Ed].node[0], nodes);
		if (iaux < 0) return -1;
		Edge_nodes[0] = nodes[0];
		iaux = mmr_edge_nodes(Mesh_id, mesh->edge[Ed].node[1], nodes);
		if (iaux < 0) return -1;
		Edge_nodes[1] = nodes[1];

	}
	else {

		printf("No nodes for free space %d\n", Ed);
		exit(-1);

	}

	return 1;
}

/*---------------------------------------------------------
mmr_edge_sons - to return edge son's numbers
---------------------------------------------------------*/
int mmr_edge_sons( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Ed,		/* in: edge ID */
		int * Edge_sons,	/* out: IDs of edge sons */
		int * Edge_mid_node	/** out: IDs of edge mid-node son */
)
{

	/* auxiliary variables */
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	/* if edge is active just substitute its sons */
	if (mesh->edge[Ed].type < 0) {
		Edge_sons[0] = mesh->edge[Ed].node[0];
		Edge_sons[1] = mesh->edge[Ed].node[1];
	}
	else {

		printf("No sons for active edge or free space %d\n", Ed);
		exit(-1);

	}

	if (Edge_mid_node != NULL) {

		printf("returning Edge_mid_node not implemented in mmr_edge_sons for prism\n");
		exit(-1);

	}

	return 1;
}

/*---------------------------------------------------------
mmr_edge_status - to return edge status (active, inactive, free space)
---------------------------------------------------------*/
int mmr_edge_status( /* returns edge status: */
		/* +1 (MMC_ACTIVE)   - active edge */
		/*  0 (MMC_FREE)     - free space */
		/* -1 (MMC_INACTIVE) - inactive (refined) edge */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Ed  	/* in: edge ID */
)
{

	/* auxiliary variables */
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	if (mesh->edge[Ed].type > 0) return (MMC_ACTIVE);
	else if (mesh->edge[Ed].type < 0) return (MMC_INACTIVE);
	else return (MMC_FREE);

}

/*---------------------------------------------------------
mmr_edge_structure - to return edge structure (e.g. for sending)
---------------------------------------------------------*/
int mmr_edge_structure( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id,	   /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Ed,  	   /* in: edge ID */
		int * Edge_struct /* out: edge structure in the form of integer array */
)
{

	/* auxiliary variables */
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	Edge_struct[0] = mesh->edge[Ed].type;
	Edge_struct[1] = mesh->edge[Ed].node[0];
	Edge_struct[2] = mesh->edge[Ed].node[1];
	//Edge_struct[3]=mesh->edge[Ed].ipid;

	return 1;

}

/*---------------------------------------------------------
mmr_set_edge_type - to set type for an edge
---------------------------------------------------------*/
int mmr_set_edge_type( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id,	   /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Edge_id,     /* in: edge ID */
		int Type         /* in: edge type (the number of attempted subdivisions !!!)*/
)
{

	/* auxiliary variables */
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	mesh->edge[Edge_id].type = Type;

	return 1;

}

/*---------------------------------------------------------
mmr_set_edge_fam - to set family data for an edge
---------------------------------------------------------*/
int mmr_set_edge_fam( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id,	   /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Edge_id,     /* in: edge ID */
		int Son1,        /* in: first son ID */
		int Son2         /* in: second son ID */
)
{

	/* auxiliary variables */
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

#ifdef DEBUG_MMM
	if (mesh->edge[Edge_id].type > 0) {
		printf("Cannot set sons' IDs for active edge %d\n", Edge_id);
		exit(-1);
	}
#endif

	if (Son1 > 0) mesh->edge[Edge_id].node[0] = Son1;
	if (Son2 > 0) mesh->edge[Edge_id].node[1] = Son2;

	return 1;

}

/*---------------------------------------------------------
mmr_set_face_fam - to set family data for an face
---------------------------------------------------------*/
int mmr_set_face_fam( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id,	   /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Face_id,     /* in: face ID */
		int * Sons      /* in: sons IDs */
)
{

	/* auxiliary variables */
	mmt_mesh * mesh;
	const int num_sons = 4;
	int ison;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

#ifdef DEBUG_MMM
	if (mesh->face[Face_id].type > 0) {
		printf("Cannot set sons' IDs for active face %d\n", Face_id);
		exit(-1);
	}
#endif

	for (ison = 0; ison < num_sons; ison++) {
		if (Sons[ison] != 0) mesh->face[Face_id].sons[ison] = Sons[ison];
	}

	return 1;

}

/*---------------------------------------------------------
mmr_set_face_neig - to set neighbors data for a face
---------------------------------------------------------*/
int mmr_set_face_neig( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id,	   /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Face_id,     /* in: face ID */
		int Neig1,       /* in: first neig ID */
		int Neig2,        /* in: second neig ID */
		int Neig1Type,    /** in: first neig type */
		int Neig2Type     /** in: second neig type */
)
{

	/* auxiliary variables */
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	if (Neig1 != 0) mesh->face[Face_id].neig[0] = Neig1;
	if (Neig2 != 0) mesh->face[Face_id].neig[1] = Neig2;

	return 1;

}

/*---------------------------------------------------------
mmr_set_elem_fam - to set family data for an elem
---------------------------------------------------------*/
int mmr_set_elem_fam( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id,	   /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Elem_id,     /* in: elem ID */
		int Fath,        /* in: father ID */
		int * Sons       /* in: sons IDs */
)
{

	/* auxiliary variables */
	mmt_mesh * mesh;
	const int num_sons = 8;
	int ison;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

#ifdef DEBUG_MMM
	if (mesh->elem[Elem_id].type > 0) {
		printf("Cannot set sons' IDs for active elem %d\n", Elem_id);
		exit(-1);
	}
#endif

	if (Fath > 0) mesh->elem[Elem_id].fath = Fath;

	for (ison = 0; ison < num_sons; ison++) {
		if (Sons[ison] > 0) mesh->elem[Elem_id].sons[ison] = Sons[ison];
	}

	return 1;

}

/*---------------------------------------------------------
mmr_set_elem_fath - to set family data for an elem
---------------------------------------------------------*/
int mmr_set_elem_fath( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id,	   /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Elem_id,     /* in: elem ID */
		int Fath        /* in: father ID */
)
{

	/* auxiliary variables */
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	if (Fath > 0) mesh->elem[Elem_id].fath = Fath;

	return 1;

}

/*---------------------------------------------------------
mmr_node_status - to return node status (active, inactive, free space)
---------------------------------------------------------*/
int mmr_node_status( /* returns node status: */
		/* +1 (MMC_ACTIVE)   - active node */
		/*  0 (MMC_FREE)     - free space */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Node  	/* in: node ID */
)
{

	/* auxiliary variables */
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	if (mesh->node[Node].x < -1e10) return (MMC_FREE);
	else return (MMC_ACTIVE);

}

/*---------------------------------------------------------
mmr_node_coor - to return node coordinates
---------------------------------------------------------*/
int mmr_node_coor( /* returns success (>=0) or error (<0) code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Node,  	/* in: node ID */
		double * Coor /* out: coordinates */
)
{

	/* auxiliary variables */
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	if (mesh->node[Node].x < -1e10) return -1;
	else {

		Coor[0] = mesh->node[Node].x;
		Coor[1] = mesh->node[Node].y;
		Coor[2] = mesh->node[Node].z;

	}
	return 1;
}

/*---------------------------------------------------------
mmr_move_node - moves node by vector
---------------------------------------------------------*/
int mmr_move_node( /* returns success (>=0) or error (<0) code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Node,  	/* in: node ID */
		double * Dx  	/* in: displacement (vector) */
)
{

	/* auxiliary variables */
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	if (mesh->node[Node].x < -1e10) return -1;
	else {
		mesh->node[Node].x += Dx[0];
		mesh->node[Node].y += Dx[1];
		mesh->node[Node].z += Dx[2];
	}

	return 1;

}

/*---------------------------------------------------------
mmr_set_node_coor - set node coordinates
---------------------------------------------------------*/
int mmr_set_node_coor( /* returns success (>=0) or error (<0) code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Node,  	/* in: node ID */
		double * Xcoor  	/* in: node coordinates */
)
{

	/* auxiliary variables */
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	mesh->node[Node].x = Xcoor[0];
	mesh->node[Node].y = Xcoor[1];
	mesh->node[Node].z = Xcoor[2];

	return 1;

}

/*---------------------------------------------------------
  mmr_el_fa_nodes - to get list local face nodes indexes in elem
---------------------------------------------------------*/
extern int mmr_el_fa_nodes( // returns face type flag
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int El,		// in: global elem ID
		int Fa,		/* in: local face number in elem El */
		int * fa_nodes		/* out: list of local indexes of face nodes */
)
{
	int type = MMC_QUAD;
	switch (Fa) {
		case 0: {
			fa_nodes[0] = 0;
			fa_nodes[1] = 1;
			fa_nodes[2] = 2;
			type = MMC_TRIA;
		}
		break;
		case 1: {
			fa_nodes[0] = 3;
			fa_nodes[1] = 4;
			fa_nodes[2] = 5;
			type = MMC_TRIA;
		}
		break;
		case 2: {
			fa_nodes[0] = 0;
			fa_nodes[1] = 1;
			fa_nodes[2] = 4;
			fa_nodes[3] = 3;
		}
		break;
		case 3: {
			fa_nodes[0] = 1;
			fa_nodes[1] = 2;
			fa_nodes[2] = 5;
			fa_nodes[3] = 4;
		}
		break;
		case 4: {
			fa_nodes[0] = 2;
			fa_nodes[1] = 0;
			fa_nodes[2] = 3;
			fa_nodes[3] = 5;
		}
		break;
	}//!switch(Fa)
	return type;
}

//END

