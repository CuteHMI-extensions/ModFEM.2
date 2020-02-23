/************************************************************************
File mms_prism_util.c - utility routines for mesh modification

Contains routines:
  mmr_loc_loc - to compute local coordinates within an element,
	given local coordinates within an element of the same family
  mmr_divide_face4_t - to break a triangular face into 4 sons
  mmr_divide_face4_q - to divide a quadrilateral face into 4 sons
  mmr_create_face - to create a new face
  mmr_divide_edge - to divide an edge into two sons
  mmr_create_edge - to create a new edge structure
  mmr_create_node - to create a node
  mmr_create_edge_elems - to create (or update) a list for each edge of elements
				  to which it belongs
  mmr_edge_elems - to return IDs of elements containing the edge

  mmr_delete_edge_elems - to delete for each edge a list of elements
				  to which it belongs

  mmr_clust_fa4_t - to cluster back family of 4 triangular faces
  mmr_clust_fa4_q - to cluster back family of 4 quadrilateral faces
  mmr_del_face    - to free a face structure
  mmr_clust_edge  - to cluster two edges back into their father
  mmr_del_edge    - to free an edge structure
  mmr_del_node    - to free a node structure

------------------------------
History:
	02.2002 - Krzysztof Banas, initial version
*************************************************************************/

#include<stdlib.h>
#include<stdio.h>
#include<math.h>

/* interface for the mesh manipulation module */
#include <modfem/mmh_intf.h>

/* mesh manipulation data structure and headers for internal routines */
#include <modfem/mm_prism/mmh_prism.h>

/* internal utility procedures */
void mut_mat3vec(
		double * m1, 	/* in: matrix (stored by rows as a vector!) */
		double * v1, 	/* in: vector */
		double * v2	/* out: resulting vector */
);
void mut_mat3mat(
		double * m1,	/* in: matrix */
		double * m2,	/* in: matrix */
		double * m3	/* out: matrix m1*m2 */
);
int mut_chk_list(	/* returns: */
		/* >0 - position on the list */
		/* 0 - not found on the list */
		int Num, 	/* number to be checked */
		int * List, 	/* list of numbers */
		int Ll		/* length of the list */
);

/*---------------------------------------------------------
mmr_loc_loc - to compute local coordinates within an element,
	given local coordinates within an element of the same family
---------------------------------------------------------*/
int mmr_loc_loc(/* returns: 1 - success, 0 - failure */
		int Mesh_id,   /* in: mesh ID */
		int El_from, 	/* in: element number */
		double * X_from, /* in: local element coordinates */
		int El_to, 	/* in: another element number */
		double * X_to	/* out: local another element coordinates */
)
{

	/* local variables */
	int igen;	/* generation counter */
	int gen_from, gen_to; /* element's generations */
	int gen_diff;	/* generation difference between element and ancestor */
	int father, elsons[MMC_MAXELSONS + 1]; /* family information */
	int ison;	/* which son are you? */
	int type_ref;	/* refinement type: h2, h3, h4 */

	/* auxiliary variables */
	int i, nel, nel_old;
	double amat[9], atrans[9], acoef[9];
	double bvec[3], btrans[3], bcoef[3];

	/*++++++++++++++++ executable statements ++++++++++++++++*/


	/* find generation levels for both elements */
	gen_from = mmr_el_gen(Mesh_id, El_from);
	gen_to = mmr_el_gen(Mesh_id, El_to);
	gen_diff = abs(gen_from - gen_to);

	/*kbw
		printf("In mmr_loc_loc: El_from %d (gen %d), El_to %d (gen %d)\n",
		   El_from, gen_from, El_to, gen_to);

	kew*/

	/* starting with the younger element */
	if (gen_to < gen_from) {
		nel = El_from;
		nel_old = El_to;
	}
	else {
		nel = El_to;
		nel_old = El_from;
	}

	if (gen_diff > 1) {
		for (i = 0; i < 9; i++) amat[i] = 0.0;
		for (i = 0; i < 3; i++) bvec[i] = 0.0;
		amat[0] = 1.0;
		amat[4] = 1.0;
		amat[8] = 1.0;
	}

	/* in a loop over generations */
	for (igen = 0; igen < gen_diff; igen++) {


		/* initialize for further generations */
		if (igen > 0) {
			nel = father;
		}

		/* find father */
		father = mmr_el_fam(Mesh_id, nel, NULL, NULL);

		/* check family informarion */
		mmr_el_fam(Mesh_id, father, elsons, &type_ref);

		/* check which son are you */
		ison = mut_chk_list(nel, &elsons[1], elsons[0]);

		/*kbw
			printf("Father %d, son %d - number %d\n",father,nel,ison);
			printf("Sons (%d): ",elsons[0]);
			for(i=1;i<=elsons[0];i++) printf(" %d",elsons[i]);
			printf("\n");
		kew*/

		/* find transformation coefficients to father element coordinates */
		if (gen_to < gen_from) {

			atrans[0] = 0.5;
			atrans[1] = 0;
			atrans[2] = 0;
			atrans[3] = 0;
			atrans[4] = 0.5;
			atrans[5] = 0;
			atrans[6] = 0;
			atrans[7] = 0;
			atrans[8] = 0.5;

			if (ison == 1) {

				btrans[0] =  0;
				btrans[1] =  0;
				btrans[2] = -0.5;

			}
			else if (ison == 2) {

				btrans[0] =  0.5;
				btrans[1] =  0;
				btrans[2] = -0.5;

			}
			else if (ison == 3) {

				btrans[0] =  0;
				btrans[1] =  0.5;
				btrans[2] = -0.5;

			}
			else if (ison == 4) {

				atrans[0] = -0.5;
				atrans[4] = -0.5;

				btrans[0] =  0.5;
				btrans[1] =  0.5;
				btrans[2] = -0.5;

			}
			else if (ison == 5) {

				btrans[0] =  0;
				btrans[1] =  0;
				btrans[2] =  0.5;

			}
			else if (ison == 6) {

				btrans[0] =  0.5;
				btrans[1] =  0;
				btrans[2] =  0.5;

			}
			else if (ison == 7) {

				btrans[0] =  0;
				btrans[1] =  0.5;
				btrans[2] =  0.5;

			}
			else if (ison == 8) {

				atrans[0] = -0.5;
				atrans[4] = -0.5;

				btrans[0] =  0.5;
				btrans[1] =  0.5;
				btrans[2] =  0.5;

			}

		}
		else {

			atrans[0] = 2;
			atrans[1] = 0;
			atrans[2] = 0;
			atrans[3] = 0;
			atrans[4] = 2;
			atrans[5] = 0;
			atrans[6] = 0;
			atrans[7] = 0;
			atrans[8] = 2;

			if (ison == 1) {

				btrans[0] =  0;
				btrans[1] =  0;
				btrans[2] =  1;

			}
			else if (ison == 2) {

				btrans[0] = -1;
				btrans[1] =  0;
				btrans[2] =  1;

			}
			else if (ison == 3) {

				btrans[0] =  0;
				btrans[1] = -1;
				btrans[2] =  1;

			}
			else if (ison == 4) {

				atrans[0] = -2;
				atrans[4] = -2;

				btrans[0] =  1;
				btrans[1] =  1;
				btrans[2] =  1;

			}
			else if (ison == 5) {

				btrans[0] =  0;
				btrans[1] =  0;
				btrans[2] = -1;

			}
			else if (ison == 6) {

				btrans[0] = -1;
				btrans[1] =  0;
				btrans[2] = -1;

			}
			else if (ison == 7) {

				btrans[0] =  0;
				btrans[1] = -1;
				btrans[2] = -1;

			}
			else if (ison == 8) {

				atrans[0] = -2;
				atrans[4] = -2;

				btrans[0] =  1;
				btrans[1] =  1;
				btrans[2] = -1;

			}

		}

		if (gen_diff > 1) {
			mut_mat3mat(atrans, amat, acoef);
			for (i = 0; i < 9; i++) amat[i] = acoef[i];
			mut_mat3vec(atrans, bvec, bcoef);
			for (i = 0; i < 3; i++) bvec[i] = btrans[i] + bcoef[i];
		}


	} /* end loop over generations: igen */

	/* if final element is identical with one sent */
	if (father == nel_old) {

		if (gen_diff > 1) {
			mut_mat3vec(amat, X_from, X_to);
			for (i = 0; i < 3; i++) X_to[i] += bvec[i];
		}
		else {
			mut_mat3vec(atrans, X_from, X_to);
			for (i = 0; i < 3; i++) X_to[i] += btrans[i];
		}


		if ( X_to[0] < 0.0 || X_to[1] < 0.0 ||
				X_to[0] + X_to[1] > 1.0 ||
				X_to[2] < -1.0 || X_to[2] > 1.0 ) return (0);


	}
	else {
		printf("Elements %d and %d do not overlap in loc_loc\n", El_from, El_to );
		exit(-1);
	}

	/*kbw
	printf("In loc-loc: from: %d, to: %d\n",El_from,El_to);
	printf("In loc-loc: older: %d, younger: %d\n",nel_old,nel);
	printf("Coordinates from: %lf %lf %lf\n",
	X_from[0],X_from[1],X_from[2]);
	printf("Coefficients : \n");
	for(i=0;i<9;i++) printf("%20.15lf",amat[i]);
	printf("\n");
	for(i=0;i<3;i++) printf("%20.15lf",bvec[i]);
	printf("\n");
	printf("Coordinates to: %lf %lf %lf\n",
	X_to[0],X_to[1],X_to[2]);
	kew*/

	return (1);
}

/*------------------------------------------------------------
  mmr_create_element - to create an element and fill its data structure
------------------------------------------------------------*/
int mmr_create_element( /* returns: ID of the created element (<=0 - failure) */
		int  Mesh_id,	 /* in: ID of the mesh to be used or 0 for the current mesh */
		int  Type,     /* in: type for the face */
		int  Mate,     /* in: material indicator */
		int  Fath,     /* in: father element ID */
		int  Refi,     /* in: refinement type indicator */
		int * Faces,   /* in: list of faces' IDs */
		int * Sons    /* in: list of sons (only for inactive elements, Type<0) */
		/*||begin||*/
		//int Ipid       /* in: inter-processor ID for new element */
		/*||end||*/
)
{

	/* local variables */
	mmt_mesh * mesh;
	const int nr_sons = 8;
	int i, iaux, iel, pel, ifa, nr_faces;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh( Mesh_id);

	/* find free data structure */
	/*take as a pointer value stored at pfel*/
	pel = mesh->parm.pfel;
	if (pel > mesh->parm.mxel) {
		printf("Not enough space for next element, Nrel= %d\n",
				mesh->parm.nrel);
		exit(-1);
	}

	/* if space already visited */
	if (pel < mesh->parm.nmel + 1) {

		/* small check - whether we have space for new element */
#ifdef DEBUG_MMM
		if (mesh->elem[pel].type != MMC_FREE) {
			printf("create_element: old pointer %d points not into free space\n",
					pel);
			exit(-1);
		}
#endif

		/* subsitute for pfel value stored in mate */
		mesh->parm.pfel = mesh->elem[pel].mate;

	} /* end if: found free space */
	else {

		/*small check - whether we are at the end of available space */
#ifdef DEBUG_MMM
		if (pel != mesh->parm.nmel + 1) {
			printf("create_element - error in pointers\n");
			exit(-1);
		}
#endif

		/* update control variables */
		mesh->parm.pfel++;
		mesh->parm.nmel++;

	}

	mesh->parm.nrel++;

	/* substitute parameters */
	mesh->elem[pel].type = Type;
	mesh->elem[pel].mate = Mate;
	mesh->elem[pel].fath = Fath;
	mesh->elem[pel].refi = Refi;
	iaux = abs(Type);
	if (iaux == MMC_BRICK) nr_faces = 6;
	else if (iaux == MMC_PRISM) nr_faces = 5;
	else if (iaux == MMC_TETRA) nr_faces = 4;
	mesh->elem[pel].face = mmr_ivector(nr_faces, "element faces in create_elem");
	for (i = 0; i < nr_faces; i++) mesh->elem[pel].face[i] = Faces[i];
	if (Type > 0) {
		mesh->elem[pel].sons = NULL;
	}
	else if (Type < 0) {
		mesh->elem[pel].sons = mmr_ivector(nr_sons, "element sons in create_elem");
		if (Sons != NULL) {
			for (i = 0; i < nr_sons; i++) mesh->elem[pel].sons[i] = Sons[i];
		}
	}
	else {
		printf("create_element - wrong type\n");
		exit(-1);
	}

	/*||begin||*/
	//mesh->elem[pel].ipid = Ipid;
	/*||end||*/

	/* check partial consistency */
#ifdef DEBUG_MMM
	{
		int el_node[MMC_MAXELVNO + 1];
		mmr_el_node_coor(Mesh_id, pel, el_node, NULL);
	}
#endif

	return (pel);
}


/*------------------------------------------------------------
  mmr_divide_face4_t - to break a triangular face into 4 sons
	Warning: face neighbours are not specified (except the case Ipids!=NULL,
	where it is assumed that one neighbor is big and the second is
	inter-subdomain boundary)
------------------------------------------------------------*/
int mmr_divide_face4_t( /* returns: 1 - success, <=0 - failure */
		int  Mesh_id,	   /* in: ID of the mesh to be used or 0 for the current mesh */
		int  Fa,         /* in: face to be divided */
		int * Face_sons, /* out: face's sons */
		int * Sons_edges, /* out: created new edges */
		int * New_nodes /* out: created new nodes */
		/*||begin||*/
		//int* Ipids       /* in: list of IPIDs for new entities */
		/*||end||*/
)
{

	/* local variables */
	mmt_mesh * mesh;
	int edge_sons[2], int_edge[3], out_edge[6];
	int iaux, ied, ison, edge, nodes[2], bc_flag;
	//, new_ipid, ipid_count;

	const int num_vert = 3;
	const int num_sons = 4;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh( Mesh_id);

	/* first check whether face is active */
	if (mmr_fa_status( Mesh_id, Fa) != MMC_ACTIVE) {
#ifdef DEBUG_MMM
		printf("not active face or free space %d in divide_face4_t!\n", Fa);
#endif
		return (-1);
	}

	/*kbwSS
	  if(Fa==669||Fa==739){
		int i;
		printf("!!!!!!!!!!!!!!!!Dividing face %d , type %d, bc %d, neig1 %d, neig2 %d, owner %d\nedges:",
		   Fa,mesh->face[Fa].type,
		   mesh->face[Fa].bc,
		   mesh->face[Fa].neig[0],
		   mesh->face[Fa].neig[1], mmr_get_owner(Mesh_id,MMC_FACE,Fa ));
		for(i=0;i<num_vert;i++){
		  printf("  %d",mesh->face[Fa].edge[i]);
		}
		printf("\n");
	  }
	kew*/

	/*||begin||*/
	/*   if(Ipids==NULL){ */
	/*     /\* prepare ipid - inter-processor ID for newly created entities *\/ */
	/*     if(mesh->face[Fa].ipid==0){ */
	/*       new_ipid = 0; */
	/*       mmv_ref_loc.face[0]++; */
	/*       mmv_ref_loc.face[mmv_ref_loc.face[0]]=Fa; */
	/*       mmv_ref_loc.face[0]++; */
	/*       mmv_ref_loc.face[mmv_ref_loc.face[0]]=mmr_fa_type(Mesh_id, Fa); */
	/*     } */
	/*     else if(mesh->face[Fa].ipid>0){ */
	/*       new_ipid = -1; */
	/*       mmv_ref_ali.face[0]++; */
	/*       mmv_ref_ali.face[mmv_ref_ali.face[0]]=Fa; */
	/*       mmv_ref_ali.face[0]++; */
	/*       mmv_ref_ali.face[mmv_ref_ali.face[0]]=mmr_fa_type(Mesh_id, Fa); */
	/*     } */
	/*     else{ */
	/*       printf("external mesh entity %d in divide_fa\n",Fa); */
	/*       exit(-1); */
	/*     } */
	/*   } */
	/*   else{ */
	/*     ipid_count=0; */
	/*   } */
	/*||end||*/

	/* divide edges, prepare list of edges for sons */
	for (ied = 0; ied < num_vert; ied++) {

		edge = mesh->face[Fa].edge[ied];

		if (mesh->edge[abs(edge)].type > 0) {

#ifdef DEBUG_MMM
			/*       if(Ipids!=NULL){ */
			/* 	printf("dividing face %d on inter-processor boundary: not divided edges!\n", */
			/* 	       Fa); */
			/* 	exit(-1); */
			/*       } */
#endif

			/* active edge, we have to divide it */
			iaux = mmr_divide_edge(Mesh_id, abs(edge), edge_sons, &New_nodes[ied]);
			//,NULL);
			if (iaux <= 0) return (-1);

			/*kbw
			printf("For face %d divided edge %d (node in the middle %d)\n",
			Fa,edge,New_nodes[ied]);
			printf("two sons: %d (nodes %d, %d) and %d (nodes %d, %d)\n",
			edge_sons[0],mesh->edge[edge_sons[0]].node[0],
			mesh->edge[edge_sons[0]].node[1],
			edge_sons[1],mesh->edge[edge_sons[1]].node[0],
			mesh->edge[edge_sons[1]].node[1]);
			kew*/

		}
		else if (mesh->edge[abs(edge)].type < 0) {

			/* increase the counter for attempted edge divisions */
			mesh->edge[abs(edge)].type--;

			/* inactive edge - we just collect its sons and middle node */
			edge_sons[0] = mesh->edge[abs(edge)].node[0];
			edge_sons[1] = mesh->edge[abs(edge)].node[1];

			/*kbw
			printf("0For face %d found inactive edge %d, two sons: %d and %d \n",
			Fa,edge,edge_sons[0],edge_sons[1]);
			kew*/

			mmr_edge_nodes( Mesh_id, edge_sons[0], nodes);
			New_nodes[ied] = nodes[1];

			/*kbw
			printf("For face %d found inactive edge %d, (node in the middle %d)\n",
			Fa,edge,New_nodes[ied]);
			printf("two sons: %d (nodes %d, %d) and %d \n",
			edge_sons[0],nodes[0],nodes[1],edge_sons[1]);
			kew*/


		}
		else {
			printf("divide_face4_t - error in face edges\n");
			exit(-1);
		}

		/* put edge sons on a list of outer edges */
		if (edge > 0) {
			out_edge[2 * ied]   =  edge_sons[0];
			out_edge[2 * ied + 1] =  edge_sons[1];
		}
		else {
			out_edge[2 * ied]   = -edge_sons[1];
			out_edge[2 * ied + 1] = -edge_sons[0];
		}

	}

	/* create three edges joining three middle nodes */
	for (ied = 0; ied < num_vert; ied++) {

		int_edge[ied] = mmr_create_edge( Mesh_id, MMC_EDGE,
						New_nodes[(ied + 2) % num_vert], New_nodes[ied]);
		//, new_ipid);
		if (int_edge[ied] <= 0) return (-1);

		/*||begin||*/
		/*     if(Ipids==NULL){ */
		/*       if(new_ipid==0){ */
		/* 	mmv_ref_loc.face[0]++; */
		/* 	mmv_ref_loc.face[mmv_ref_loc.face[0]]=int_edge[ied]; */
		/*       } */
		/*       else if(new_ipid==-1){ */
		/* 	mmv_ref_ali.face[0]++; */
		/* 	mmv_ref_ali.face[mmv_ref_ali.face[0]]=int_edge[ied]; */
		/*       } */
		/*     } */
		/*     else{ */
		/*       mmr_ed_set_ipid(Mesh_id, int_edge[ied], Ipids[ipid_count]); */
		/*       ipid_count++; */
		/*     } */
		/*||end||*/

	}

	/* create new faces */
	bc_flag = mesh->face[Fa].bc;
	for (ison = 0; ison < num_sons; ison++) {

		/* edge numbers */

		if (ison < num_sons - 1) {

			/* two edges are inherited from father, one is newly created */
			Sons_edges[3 * ison + ison] = out_edge[2 * ison];
			Sons_edges[3 * ison + (ison + 1) % num_vert] = -int_edge[ison];
			Sons_edges[3 * ison + (ison + 2) % num_vert] =
					out_edge[2 * ((ison + 2) % num_vert) + 1];

		}
		else {

			/* fourth son */
			Sons_edges[3 * ison + 0] = int_edge[2];
			Sons_edges[3 * ison + 1] = int_edge[0];
			Sons_edges[3 * ison + 2] = int_edge[1];

		}

		Face_sons[ison] = mmr_create_face(Mesh_id, MMC_TRIA, bc_flag,
						&Sons_edges[3 * ison], NULL, NULL);
		//, new_ipid);

		/*||begin||*/
		/*     if(Ipids==NULL){ */
		/*       if(new_ipid==0){ */
		/* 	mmv_ref_loc.face[0]++; */
		/* 	mmv_ref_loc.face[mmv_ref_loc.face[0]]=Face_sons[ison]; */
		/*       } */
		/*       else if(new_ipid==-1){ */
		/* 	mmv_ref_ali.face[0]++; */
		/* 	mmv_ref_ali.face[mmv_ref_ali.face[0]]=Face_sons[ison]; */
		/*       } */
		/*     } */
		/*     else{ */
		/*       mmr_fa_set_ipid(Mesh_id, Face_sons[ison], Ipids[ipid_count]); */
		/*       ipid_count++; */

		/*       if(mesh->face[Fa].neig[0]==MMC_SUB_BND){ */
		/* 	mesh->face[Face_sons[ison]].neig[0]=MMC_SUB_BND; */
		/* 	mesh->face[Face_sons[ison]].neig[1]=-1; */
		/*       } */
		/*       else if(mesh->face[Fa].neig[1]==MMC_SUB_BND){ */
		/* 	mesh->face[Face_sons[ison]].neig[1]=MMC_SUB_BND; */
		/* 	mesh->face[Face_sons[ison]].neig[0]=-1; */
		/*       } */
		/*       else{ */
		/* 	printf("dividing face with Ipids, not on inter-subdomain boundary!\n"); */
		/* 	exit(-1); */
		/*       } */
		/*     } */
		/*||end||*/

		/*kbw
			if(Fa==68){
			  int i;
			  printf("Son %d: face %d , type %d, bc %d, neig1 %d, neig2 %d\nedges:",
				 ison,Face_sons[ison],mesh->face[Face_sons[ison]].type,
				 mesh->face[Face_sons[ison]].bc,
				 mesh->face[Face_sons[ison]].neig[0],
				 mesh->face[Face_sons[ison]].neig[1] );
			  for(i=0;i<num_vert;i++){
			printf("  %d",mesh->face[Face_sons[ison]].edge[i]);
			  }
			  printf("\n");
			  printf("ipid: %d\n", mmr_fa_ipid(Mesh_id,Face_sons[ison]));
			}
		kew*/

	}

	/* update parameters for the father */
	mesh->parm.nrfa--;
	mesh->face[Fa].type *= -1;
	mesh->face[Fa].sons = (int *) malloc(num_sons * sizeof(int));
	for (ison = 0; ison < num_sons; ison++) mesh->face[Fa].sons[ison] = Face_sons[ison];

	return (1);
}

/*------------------------------------------------------------
mmr_divide_face4_q - to divide a quadrilateral face into 4 sons
------------------------------------------------------------*/
int mmr_divide_face4_q( /* returns: 1 - success, <=0 - failure */
		int  Mesh_id,	   /* in: ID of the mesh to be used or 0 for the current mesh */
		int  Fa,         /* in: face to be divided */
		int * Face_sons, /* out: sons */
		int * Sons_edges, /* out: created new edges */
		int * New_nodes /* out: created new nodes */
		/*||begin||*/
		//int* Ipids       /* in: list of IPIDs for new entities */
		/*||end||*/
)
{
	/* local variables */
	mmt_mesh * mesh;
	int edge_sons[2], int_edge[4], out_edge[8];
	int i, iaux, ied, ison, edge, nodes[2], bc_flag;
	//, new_ipid, ipid_count;
	double daux, faux, eaux;

	const int num_vert = 4;
	const int num_sons = 4;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh( Mesh_id);

	/* first check whether face is active */
	if (mmr_fa_status( Mesh_id, Fa) != MMC_ACTIVE) {
#ifdef DEBUG_MMM
		printf("not active face %d in divide_face4_t!\n", Fa);
#endif
		return (-1);
	}

	/*kbw
	  if(Fa==669||Fa==739){
		printf("!!!!!!!!!!!!!!!!Dividing face %d , type %d, bc %d, neig1 %d, neig2 %d, owner %d\nedges:",
		   Fa,mesh->face[Fa].type,
		   mesh->face[Fa].bc,
		   mesh->face[Fa].neig[0],
		   mesh->face[Fa].neig[1], mmr_get_owner(Mesh_id,MMC_FACE,Fa) );
		for(i=0;i<num_vert;i++){
		  printf("  %d",mesh->face[Fa].edge[i]);
		}
		printf("\n");
	  }
	kew*/

	/*||begin||*/
	/*   if(Ipids==NULL){ */
	/*   /\* prepare ipid - inter-processor ID for newly created entities *\/ */
	/*     if(mesh->face[Fa].ipid==0){ */
	/*       new_ipid = 0; */
	/*       mmv_ref_loc.face[0]++; */
	/*       mmv_ref_loc.face[mmv_ref_loc.face[0]]=Fa; */
	/*       mmv_ref_loc.face[0]++; */
	/*       mmv_ref_loc.face[mmv_ref_loc.face[0]]=mmr_fa_type(Mesh_id, Fa); */
	/*     } */
	/*     else if(mesh->face[Fa].ipid>0){ */
	/*       new_ipid = -1; */
	/*       mmv_ref_ali.face[0]++; */
	/*       mmv_ref_ali.face[mmv_ref_ali.face[0]]=Fa; */
	/*       mmv_ref_ali.face[0]++; */
	/*       mmv_ref_ali.face[mmv_ref_ali.face[0]]=mmr_fa_type(Mesh_id, Fa); */
	/*     } */
	/*     else{ */
	/*       printf("external mesh entity %d in divide_fa\n",Fa); */
	/*       exit(-1); */
	/*     } */
	/*   } */
	/*   else{ */
	/*     ipid_count=0; */
	/*   } */
	/*||end||*/

	/* divide edges, prepare list of edges for sons */
	for (ied = 0; ied < num_vert; ied++) {

		edge = mesh->face[Fa].edge[ied];

		if (mesh->edge[abs(edge)].type > 0) {

#ifdef DEBUG_MMM
			/*       if(Ipids!=NULL){ */
			/* 	printf("dividing face %d on inter-processor boundary: not divided edges!\n", */
			/* 	       Fa); */
			/* 	exit(-1); */
			/*       } */
#endif

			/* active edge, we have to divide it */
			iaux = mmr_divide_edge(Mesh_id, abs(edge), edge_sons, &New_nodes[ied]);
			//,NULL);
			if (iaux <= 0) return (-1);

			/*kbw
				  if(Fa==68){
				printf("For face %d divided edge %d (node in the middle %d)\n",
					   Fa,edge,New_nodes[ied]);
				printf("two sons: %d (nodes %d, %d) and %d (nodes %d, %d)\n",
					   edge_sons[0],mesh->edge[edge_sons[0]].node[0],
					   mesh->edge[edge_sons[0]].node[1],
					   edge_sons[1],mesh->edge[edge_sons[1]].node[0],
					   mesh->edge[edge_sons[1]].node[1]);
				  }
			kew*/

		}
		else if (mesh->edge[abs(edge)].type < 0) {

			/* increase the counter for attempted edge divisions */
			mesh->edge[abs(edge)].type--;

			/* inactive edge - we just collect its sons and middle node */
			edge_sons[0] = mesh->edge[abs(edge)].node[0];
			edge_sons[1] = mesh->edge[abs(edge)].node[1];

			mmr_edge_nodes( Mesh_id, edge_sons[0], nodes);
			New_nodes[ied] = nodes[1];

			/*kbw
				  if(Fa==68){
				printf("For face %d found inactive edge %d, (node in the middle %d)\n",
					   Fa,edge,New_nodes[ied]);
				printf("two sons: %d (nodes %d, %d) and %d \n",
					   edge_sons[0],nodes[0],nodes[1],edge_sons[1]);
				  }
			kew*/


		}
		else {
			printf("divide_face4_q - error in face edges\n");
			exit(-1);
		}

		/* put edge sons on a list of outer edges */
		if (edge > 0) {
			out_edge[2 * ied]   =  edge_sons[0];
			out_edge[2 * ied + 1] =  edge_sons[1];
		}
		else {
			out_edge[2 * ied]   = -edge_sons[1];
			out_edge[2 * ied + 1] = -edge_sons[0];
		}

	}

	/* create a new node in the middle of a face */
	daux = 0.0;
	faux = 0.0;
	eaux = 0.0;
	for (ied = 0; ied < num_vert; ied++) {
		daux += mesh->node[New_nodes[ied]].x;
		faux += mesh->node[New_nodes[ied]].y;
		eaux += mesh->node[New_nodes[ied]].z;
	}
	New_nodes[4] = mmr_create_node(Mesh_id, daux / 4.0, faux / 4.0, eaux / 4.0);
	//, new_ipid);

#ifdef DEBUG_MMM
	if (New_nodes[4] <= 0) {
		printf("cannot create node in divide face\n");
		exit(-1);
	}
#endif

	/*||begin||*/
	/*   if(Ipids==NULL){ */
	/*     if(new_ipid==0){ */
	/*       mmv_ref_loc.face[0]++; */
	/*       mmv_ref_loc.face[mmv_ref_loc.face[0]]=New_nodes[4]; */
	/*     } */
	/*     else if(new_ipid==-1){ */
	/*       mmv_ref_ali.face[0]++; */
	/*       mmv_ref_ali.face[mmv_ref_ali.face[0]]=New_nodes[4]; */
	/*     } */
	/*   } */
	/*   else{ */
	/*     mmr_ve_set_ipid(Mesh_id, New_nodes[4], Ipids[ipid_count]); */
	/*     ipid_count++; */
	/*   } */
	/*||end||*/

	/* for each old edge create a new one joining it with node in the middle */
	for (ied = 0; ied < num_vert; ied++) {

		int_edge[ied] = mmr_create_edge(Mesh_id, MMC_EDGE,
						New_nodes[ied], New_nodes[4]);
		//, new_ipid);
		if (int_edge[ied] <= 0) return (-1);

		/*||begin||*/
		/*     if(Ipids==NULL){ */
		/*       if(new_ipid==0){ */
		/* 	mmv_ref_loc.face[0]++; */
		/* 	mmv_ref_loc.face[mmv_ref_loc.face[0]]=int_edge[ied]; */
		/*       } */
		/*       else if(new_ipid==-1){ */
		/* 	mmv_ref_ali.face[0]++; */
		/* 	mmv_ref_ali.face[mmv_ref_ali.face[0]]=int_edge[ied]; */
		/*       } */
		/*     } */
		/*     else{ */
		/*       mmr_ed_set_ipid(Mesh_id, int_edge[ied], Ipids[ipid_count]); */
		/*       ipid_count++; */
		/*     } */
		/*||end||*/

	}

	/* create new faces */
	bc_flag = mesh->face[Fa].bc;
	for (ison = 0; ison < num_sons; ison++) {
		/* edge numbers */
		/* two edges are inherited from father, two are newly created */
		Sons_edges[4 * ison + ison] = out_edge[2 * ison];
		Sons_edges[4 * ison + (ison + 1) % num_vert] = int_edge[ison];
		Sons_edges[4 * ison + (ison + 2) % num_vert] = -int_edge[(ison + 3) % num_vert];
		Sons_edges[4 * ison + (ison + 3) % num_vert] = out_edge[2 * ((ison + 3) % num_vert) + 1];

		Face_sons[ison] = mmr_create_face(Mesh_id, MMC_QUAD, bc_flag,
						&Sons_edges[4 * ison], NULL, NULL);
		//, new_ipid);

		/*||begin||*/
		/*     if(Ipids==NULL){ */
		/*       if(new_ipid==0){ */
		/* 	mmv_ref_loc.face[0]++; */
		/* 	mmv_ref_loc.face[mmv_ref_loc.face[0]]=Face_sons[ison]; */
		/*       } */
		/*       else if(new_ipid==-1){ */
		/* 	mmv_ref_ali.face[0]++; */
		/* 	mmv_ref_ali.face[mmv_ref_ali.face[0]]=Face_sons[ison]; */
		/*       } */
		/*     } */
		/*     else{ */
		/*       mmr_fa_set_ipid(Mesh_id, Face_sons[ison], Ipids[ipid_count]); */
		/*       ipid_count++; */

		/*       if(mesh->face[Fa].neig[0]==MMC_SUB_BND){ */
		/* 	mesh->face[Face_sons[ison]].neig[0]=MMC_SUB_BND; */
		/* 	mesh->face[Face_sons[ison]].neig[1]=-1; */
		/*       } */
		/*       else if(mesh->face[Fa].neig[1]==MMC_SUB_BND){ */
		/* 	mesh->face[Face_sons[ison]].neig[1]=MMC_SUB_BND; */
		/* 	mesh->face[Face_sons[ison]].neig[0]=-1; */
		/*       } */
		/*       else{ */
		/* 	printf("dividing face with Ipids, not on inter-subdomain boundary!\n"); */
		/* 	exit(-1); */
		/*       } */
		/*     } */
		/*||end||*/

		/*kbw
			if(Fa==68){
			  int i;
			  printf("Son %d: face %d , type %d, bc %d, neig1 %d, neig2 %d\nedges:",
				 ison,Face_sons[ison],mesh->face[Face_sons[ison]].type,
				 mesh->face[Face_sons[ison]].bc,
				 mesh->face[Face_sons[ison]].neig[0],
				 mesh->face[Face_sons[ison]].neig[1] );
			  for(i=0;i<num_vert;i++){
			printf("  %d",mesh->face[Face_sons[ison]].edge[i]);
			  }
			  printf("\n");
			  printf("ipid: %d\n", mmr_fa_ipid(Mesh_id,Face_sons[ison]));
			}
		kew*/

	}

	/* update parameters for the father */
	mesh->parm.nrfa--;
	mesh->face[Fa].type *= -1;
	mesh->face[Fa].sons = (int *) malloc(num_sons * sizeof(int));
	for (ison = 0; ison < num_sons; ison++) mesh->face[Fa].sons[ison] = Face_sons[ison];

	return (1);
}

/*------------------------------------------------------------
mmr_create_face - to create a new face
------------------------------------------------------------*/
int mmr_create_face( /* returns: number of the created face (<=0 - failure) */
		int  Mesh_id,	 /* in: ID of the mesh to be used or 0 for the current mesh */
		int  Type,     /* in: type for the face */
		int  Flag_bc,  /* in: bc flag for the face */
		int * Edges,   /* in: edges for the new face */
		int * Neig,    /* in (optional): neighbors (or NULL) */
		int * Sons     /* in (optional): sons (if Type<0) or NULL
/*||begin||*/
		//int Ipid       /* in: inter-processor ID the new face */
		/*||end||*/
)
{

	/* local variables */
	mmt_mesh * mesh;
	int i, pfa, ied, nr_edges;
	const int num_sons = 4;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh( Mesh_id);

	/* create a new face*/
	/*take as a pointer to its space value stored at pffa*/
	pfa = mesh->parm.pffa;
	if (pfa > mesh->parm.mxfa) {
		printf("Not enough space for next face, Nrfa= %d\n",
				mesh->parm.nrfa);
		exit(-1);
	}

	/* if space already visited */
	if (pfa < mesh->parm.nmfa + 1) {

		/* small check - whether we have space for new face */
#ifdef DEBUG_MMM
		if (mesh->face[pfa].type != MMC_FREE) {
			printf("Create_face_t: old pointer %d points not into free space\n",
					pfa);
			return (-1);
		}
#endif

		/* subsitute for pffa value stored in bc */
		mesh->parm.pffa = mesh->face[pfa].bc;

	} /* end if: found free space */
	else {

		/*small check - whether we are at the end of available space */
#ifdef DEBUG_MMM
		if (pfa != mesh->parm.nmfa + 1) {
			printf("break4t - error in pointers\n");
			return (-1);
		}
#endif

		/* update control variables */
		mesh->parm.pffa++;
		mesh->parm.nmfa++;

	}

	mesh->parm.nrfa++;

	/* substitute type, bc and edges */
	mesh->face[pfa].type = Type;
	mesh->face[pfa].bc = Flag_bc;
	if (abs(Type) == MMC_TRIA) nr_edges = 3;
	else nr_edges = 4;
	for (ied = 0; ied < nr_edges; ied++) mesh->face[pfa].edge[ied] = Edges[ied];

	if (Neig != NULL) {
		mesh->face[pfa].neig[0] = Neig[0];
		mesh->face[pfa].neig[1] = Neig[1];
	}

	if (Type > 0) {
		mesh->face[pfa].sons = NULL;
	}
	else {
		mesh->face[pfa].sons = (int *) malloc(num_sons * sizeof(int));
		if (Sons != NULL) {
			for (i = 0; i < num_sons; i++) mesh->face[pfa].sons[i] = Sons[i];
		}
	}

	/*||begin||*/
	//mesh->face[pfa].ipid=Ipid;
	/*||end||*/

	/*kbw
	  if(pfa==669||pfa==739){
	printf("Created new face %d, type %d, bc %d, owner %d\n edges:",
	pfa, mesh->face[pfa].type, mesh->face[pfa].bc, mmr_get_owner(Mesh_id,MMC_FACE,pfa));
	for(ied=0;ied<nr_edges;ied++) printf("  %d", mesh->face[pfa].edge[ied]);
	printf("\n");
	  }
	kew*/

	return (pfa);
}

/*------------------------------------------------------------
mmr_divide_edge - to divide an edge into two sons
------------------------------------------------------------*/
int mmr_divide_edge( /* returns: 1 - success, <=0 - failure */
		int  Mesh_id,	  /* in: ID of the mesh to be used or 0 for the current mesh */
		int  Edge,      /* in: edge to be divided */
		int * Edge_sons, /* out: two sons */
		int * Node_mid /* out: node in the middle */
		/*||begin||*/
		//int* Ipids      /* in: list of IPIDs for new entities */
		/*||end||*/
)
{

	/* local variables */
	mmt_mesh * mesh;
	int iaux, ino1, ino2, ison, node;
	//, new_ipid, ipid_count;
	double daux, faux, eaux;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh( Mesh_id);

	/*kbw
	if(Edge==47){
	printf("Dividing edge %d, type %d, nodes %d, %d\n",
	Edge, mesh->edge[Edge].type,
	mesh->edge[Edge].node[0],mesh->edge[Edge].node[1]);
	}
	kew*/

	/*||begin||*/
	/*   if(Ipids==NULL){ */
	/*     /\* prepare ipid - inter-processor ID for newly created entities *\/ */
	/*     if(mesh->edge[Edge].ipid==0){ */
	/*       new_ipid = 0; */
	/*       mmv_ref_loc.edge[0]++; */
	/*       mmv_ref_loc.edge[mmv_ref_loc.edge[0]]=Edge; */
	/*     } */
	/*     else if(mesh->edge[Edge].ipid>0){ */
	/*       new_ipid = -1; */
	/*       mmv_ref_ali.edge[0]++; */
	/*       mmv_ref_ali.edge[mmv_ref_ali.edge[0]]=Edge; */
	/*     } */
	/*     else{ */
	/*       printf("external mesh entity %d in divide_edge\n",Edge); */
	/*       exit(-1); */
	/*     } */
	/*   } */
	/*   else{ */
	/*     ipid_count=0; */
	/*   } */
	/*||end||*/

	/* create a new node in the middle */
	iaux = mesh->edge[Edge].node[0];
	daux = 0.5 * mesh->node[iaux].x;
	faux = 0.5 * mesh->node[iaux].y;
	eaux = 0.5 * mesh->node[iaux].z;
	iaux = mesh->edge[Edge].node[1];
	daux += 0.5 * mesh->node[iaux].x;
	faux += 0.5 * mesh->node[iaux].y;
	eaux += 0.5 * mesh->node[iaux].z;
	node = mmr_create_node(Mesh_id, daux, faux, eaux);
	//, new_ipid);
	if (node < 0) return (-1);

	/*||begin||*/
	/*   if(Ipids==NULL){ */
	/*     if(new_ipid==0){ */
	/*       mmv_ref_loc.edge[0]++; */
	/*       mmv_ref_loc.edge[mmv_ref_loc.edge[0]]=node; */
	/*     } */
	/*     else if(new_ipid==-1){ */
	/*       mmv_ref_ali.edge[0]++; */
	/*       mmv_ref_ali.edge[mmv_ref_ali.edge[0]]=node; */
	/*     } */
	/*   } */
	/*   else{ */
	/*     mmr_ve_set_ipid(Mesh_id, node, Ipids[ipid_count]); */
	/*     ipid_count++; */
	/*   } */
	/*||end||*/

	/* for each created edge */
	for (ison = 0; ison < 2; ison++) {

		/* substitute node numbers */
		if (ison == 0) {
			ino1 = mesh->edge[Edge].node[0];
			ino2 = node;
		}
		else {
			ino1 = node;
			ino2 = mesh->edge[Edge].node[1];
		}

		Edge_sons[ison] = mmr_create_edge(Mesh_id, MMC_EDGE, ino1, ino2);
		//, new_ipid);
		if (Edge_sons[ison] < 0) return (-1);

		/*||begin||*/
		/*     if(Ipids==NULL){ */
		/*       if(new_ipid==0){ */
		/* 	mmv_ref_loc.edge[0]++; */
		/* 	mmv_ref_loc.edge[mmv_ref_loc.edge[0]]=Edge_sons[ison]; */
		/*       } */
		/*       else if(new_ipid==-1){ */
		/* 	mmv_ref_ali.edge[0]++; */
		/* 	mmv_ref_ali.edge[mmv_ref_ali.edge[0]]=Edge_sons[ison]; */
		/*       } */
		/*     } */
		/*     else{ */
		/*       mmr_ed_set_ipid(Mesh_id, Edge_sons[ison], Ipids[ipid_count]); */
		/*       ipid_count++; */
		/*     } */
		/*||end||*/

	}

	/* update information for divided face: indicate it is inactive */
	mesh->parm.nred--;
	/* initialize the counter for attempted divisions */
	mesh->edge[Edge].type = -1;
	mesh->edge[Edge].node[0] = Edge_sons[0];
	mesh->edge[Edge].node[1] = Edge_sons[1];

	/*kbw
	if(Edge==47){
	printf("Divided edge %d, type %d, sons %d, %d\n",
	Edge, mesh->edge[Edge].type,
	mesh->edge[Edge].node[0],mesh->edge[Edge].node[1]);
	}
	kew*/

	*Node_mid = node;

	return (1);

}

/*------------------------------------------------------------
mmr_create_edge - to create a new edge structure
------------------------------------------------------------*/
int mmr_create_edge( /* returns: ID of a new edge (<=0 - failure)*/
		int  Mesh_id,	 /* in: ID of the mesh to be used or 0 for the current mesh */
		int Type,      /* in: type indicator (MMC_EDGE or */
		/*     number of attempted divisions for inactive edges0 */
		int  Node1,    /* in: nodes for the new edge */
		int  Node2
		/*||begin||*/
		//int Ipid       /* in: inter-processor ID the new edge */
		/*||end||*/
)
{

	/* local variables */
	mmt_mesh * mesh;
	int edge;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh( Mesh_id);

	/* create a new edge*/
	/*take as a pointer to its space value stored at pfed*/
	edge = mesh->parm.pfed;
	if (edge > mesh->parm.mxed) {
		printf("Not enough space for next edge, Nred = %d\n",
				mesh->parm.nred);
		exit(-1);
	}

	/* if space already visited */
	if (edge < mesh->parm.nmed + 1) {

		/* small check - whether we have space for new edge */
#ifdef DEBUG_MMM
		if (mesh->edge[edge].type != MMC_FREE) {
			printf("Create_edge: old pointer %d points not into free space\n",
					edge);
			exit(-1);
		}
#endif

		/* subsitute for pfed value stored in node[0] */
		mesh->parm.pfed = mesh->edge[edge].node[0];

	} /* end if: found free space */
	else {

		/*small check - whether we are at the end of available space */
#ifdef DEBUG_MMM
		if (edge != mesh->parm.nmed + 1) {
			printf("create_new_edge - error in pointers\n");
			exit(-1);
		}
#endif

		/* update control variables */
		mesh->parm.pfed++;
		mesh->parm.nmed++;

	}

	mesh->parm.nred++;

	/* substitute type and nodes */
	mesh->edge[edge].type = Type;
	mesh->edge[edge].node[0] = Node1;
	mesh->edge[edge].node[1] = Node2;

	/* initialize edge[].elems structure */
	mesh->edge[edge].elems = NULL;

	/*||begin||*/
	//mesh->edge[edge].ipid=Ipid;
	/*||end||*/

	/*kbw
	printf("Created edge %d: type %d, nodes %d, %d\n",
	edge, mesh->edge[edge].type,
	mesh->edge[edge].node[0],mesh->edge[edge].node[1]);
	kew*/

	return (edge);
}


/*------------------------------------------------------------
mmr_create_node - to create a node at point (Xcoor,Ycoor,Zcoor)
------------------------------------------------------------*/
int mmr_create_node( /* returns: node number of created node (<=0 - failure) */
		int  Mesh_id,	 /* in: ID of the mesh to be used or 0 for the current mesh */
		double Xcoor,  /* in: coordinates of new node */
		double Ycoor,
		double Zcoor
		/*||begin||*/
		//int Ipid       /* in: inter-processor ID the new node */
		/*||end||*/
)
{

	mmt_mesh * mesh;
	int pno;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh( Mesh_id);

	/* create a new node */
	/*take as a pointer to its space value stored at pfno*/
	pno = mesh->parm.pfno;
	if (pno > mesh->parm.mxno) {
		printf("Not enough space for next node, Nrno= %d\n",
				mesh->parm.nrno);
		exit(-1);
	}

	/* if space already visited */
	if (pno < mesh->parm.nmno + 1) {

#ifdef DEBUG_MMM
		/* small check - whether we have space for new node */
		if (mesh->node[pno].x > -1e10) {
			printf("Create node: old pointer %d points not into free space\n",
					pno);
			return (-1);
		}
#endif

		/* substitute for pfno value stored in y coordinate */
		mesh->parm.pfno = mesh->node[pno].y;

	} /* end if: found free space */
	else {

#ifdef DEBUG_MMM
		/*small check - whether we are at the end of available space */
		if (pno != mesh->parm.nmno + 1) {
			printf("create node - error in pointers\n");
			return (-1);
		}
#endif

		/* update control variables */
		mesh->parm.pfno++;
		mesh->parm.nmno++;

	}

	mesh->parm.nrno++;

	/* substitute values */
	mesh->node[pno].x = Xcoor;
	mesh->node[pno].y = Ycoor;
	mesh->node[pno].z = Zcoor;

	/*||begin||*/
	//mesh->node[pno].ipid=Ipid;
	/*||end||*/

	/*kbw
	printf("Created node %d: coord: %lf, %lf %lf\n",
	pno,mesh->node[pno].x,
	mesh->node[pno].y,mesh->node[pno].z);
	kew*/

	return (pno);
}


/*------------------------------------------------------------
  mmr_create_edge_elems - to create (or recreate if Max_edge_id > 0) for each
				  edge a list of elements to which it belongs
------------------------------------------------------------*/
int mmr_create_edge_elems( /* returns: 1-success, <=0-failure */
		int  Mesh_id,	 /* in: ID of the mesh to be used or 0 for the current mesh */
		int Max_edge_id /* in: the range of edge IDs to consider or 0 for default */
		/* if Max_edge_id==0 it is assumed that there are no structures to free !!! */
)
{
	/* local variables */
	mmt_mesh * mesh;
	int ned, nel, ied, el_edges[13];

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* deallocate the space for elems structure for active and inactive edges */
	mmr_delete_edge_elems(mmv_cur_mesh_id, Max_edge_id);

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh( Mesh_id);

	/* in a loop over all edges create elems tables */
	ned = 0;
	while ((ned = mmr_get_next_edge_all(Mesh_id, ned)) != 0) {

		mesh->edge[ned].elems = (int *) malloc((MMC_MAX_EDGE_ELEMS + 1) * sizeof(int));
		mesh->edge[ned].elems[0] = 0;

	}

	/* in a loop over all elements */
	nel = 0;
	while ((nel = mmr_get_next_elem_all(Mesh_id, nel)) != 0) {

		/* get element's edges */
		mmr_el_edges(Mesh_id, nel, el_edges);

		/* put the element on elems lists */
		for (ied = 1; ied <= el_edges[0]; ied++) {

			ned = el_edges[ied];
			mesh->edge[ned].elems[0]++;

#ifdef DEBUG_MMM
			if (mesh->edge[ned].elems[0] > MMC_MAX_EDGE_ELEMS) {
				printf("Increase MMC_MAX_EDGE_ELEMS. At least to %d\n",
						mesh->edge[ned].elems[0]);
				exit(-1);
			}
#endif

			mesh->edge[ned].elems[mesh->edge[ned].elems[0]] = nel;

		}

	}

#ifdef DEBUG_MMM
	ned = 0;
	while ((ned = mmr_get_next_edge_all(Mesh_id, ned)) != 0) {
		int i, iel, edge_elems[MMC_MAX_EDGE_ELEMS + 1];
		mmr_edge_elems(Mesh_id, ned, edge_elems);
		for (iel = 1; iel <= edge_elems[0]; iel++) {
			nel = edge_elems[iel];
			mmr_el_edges(Mesh_id, nel, el_edges);
			if (mut_chk_list(ned, &el_edges[1], el_edges[0]) == 0) {
				printf("Error 2345398 in mmr_create_edge_elems. Exiting!\n");
				printf("Edge %d, nr_elems %d, el %d, nr_edges %d\nedge_elems: ",
						ned, edge_elems[0], nel, el_edges[0]);
				for (i = 1; i <= edge_elems[0]; i++) {
					printf("%d  ", edge_elems[i]);
				}
				printf("\nel_edges: ");
				for (ied = 1; ied < el_edges[0]; ied++) {
					printf("%d  ", el_edges[ied]);
				}
				printf("\n");
				exit(-1);
			}
		}
	}
#endif


	return (1);
}

/*---------------------------------------------------------
mmr_edge_elems - to return IDs of elements containing the edge
---------------------------------------------------------*/
extern int mmr_edge_elems( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int Ed,		/* in: edge ID */
		int * Edge_elems	/* out: IDs of elements containing the edge */
		/*      Edge_elems[0] - the number of elements */
)
{

	/* local variables */
	mmt_mesh * mesh;
	int i;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh( Mesh_id);

	for (i = 0; i <= mesh->edge[Ed].elems[0]; i++) {

		Edge_elems[i] = mesh->edge[Ed].elems[i];

	}

	return (1);
}

/*------------------------------------------------------------
  mmr_delete_edge_elems - to delete for each edge a list of elements
				  to which it belongs
------------------------------------------------------------*/
int mmr_delete_edge_elems( /* returns: 1-success, <=0-failure */
		int  Mesh_id,	 /* in: ID of the mesh to be used or 0 for the current mesh */
		int Max_edge_id /* in: the range of edge IDs to consider or 0 for default */
		/* if Max_edge_id==0 it is assumed that there are no structures to free !!! */
)
{
	/* local variables */
	mmt_mesh * mesh;
	int ned, iaux;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh( Mesh_id);

	if (Max_edge_id == 0) {
		iaux = mmr_get_max_edge_id(Mesh_id);
		for (ned = 1; ned <= iaux; ned++) mesh->edge[ned].elems = NULL;
	}
	else {
		for (ned = 1; ned <= Max_edge_id; ned++) {
			if (mesh->edge[ned].elems != NULL) free(mesh->edge[ned].elems);
			mesh->edge[ned].elems = NULL;
		}
	}

	return (1);
}

/*------------------------------------------------------------
mmr_clust_fa4_t - to cluster back family of 4 triangular faces
------------------------------------------------------------*/
int mmr_clust_fa4_t( /* returns: 1-success, <=0-failure */
		int  Mesh_id,	 /* in: ID of the mesh to be used or 0 for the current mesh */
		int Face,      /* in: father face ID */
		int * Face_sons /* in (optional): face sons' IDs */
)
{

	/* local variables */
	mmt_mesh * mesh;
	int iaux, ied, edge, face_sons[5], ison;

	const int nr_edges = 3;
	const int nr_sons = 4;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh( Mesh_id);

	/* check whether face is inactive */
#ifdef DEBUG_MMM
	if (mmr_fa_status( Mesh_id, Face) != MMC_INACTIVE) {
		printf("not an inactive face %d in cluster face\n", Face);
		exit(-1);
	}
#endif

	/* find sons */
	if (Face_sons == NULL) mmr_fa_fam( Mesh_id, Face, face_sons, NULL);
	else {

		face_sons[0] = nr_sons;
		for (ison = 1; ison <= nr_sons; ison++) {
			face_sons[ison] = Face_sons[ison];
		}

	}

	/* check whether all sons are active */
#ifdef DEBUG_MMM
	for (ison = 1; ison <= nr_sons; ison++) {
		if (mmr_fa_status( Mesh_id, face_sons[ison]) != MMC_ACTIVE) {
			printf("not an active face son %d in cluster face %d\n",
					face_sons[ison], Face);
			exit(-1);
		}
	}
#endif

	/*kbw
	  if(Face==669||Face==739){
		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Clustering back face %d, type %d, bc %d, owner %d\n",Face,mesh->face[Face].type, mesh->face[Face].bc, mmr_get_owner(Mesh_id,MMC_FACE,Face ));
	}
	kew*/
	/*kbw
	  if(Face==669||Face==739){
	  for(ison=1;ison<=nr_sons;ison++){
	printf("Clustering back face %d, type %d, bc %d\n",
	face_sons[ison],mesh->face[face_sons[ison]].type,
	mesh->face[face_sons[ison]].bc);
	printf("edges: ");
	for(ied=0;ied<nr_edges;ied++) printf(" %d ",
	mesh->face[face_sons[ison]].edge[ied]);
	printf("\nneighbors %d %d\n",mesh->face[face_sons[ison]].neig[0],
	mesh->face[face_sons[ison]].neig[1]);
	  }
	}
	kew*/

	/* cluster back six boundary edges of a family (if possible) */
	for (ied = 0; ied < nr_edges; ied++) {
		edge = abs(mesh->face[Face].edge[ied]);


		/*||begin||*/
		/* if edge is owned */
		//if(mesh->edge[edge].ipid==0){
		/*||end||*/
		if (mesh->edge[edge].type < -1) {
			/* increase the counter of attempted divisions */
			mesh->edge[edge].type++;

			/*kbw
			if(Face>0){
			printf("Attempted clustering of edge %d, type %d\n",
			edge,mesh->edge[edge].type);
			}
			kew*/

		}
		else {

			/* cluster back edge */
			iaux = mmr_clust_edge(Mesh_id, edge);
			if (iaux <= 0) return (-1);

			/*||begin||*/
			/* put on the list of clustered edges */
			//mmv_del_loc.edge[0]++;
			//mmv_del_loc.edge[mmv_del_loc.edge[0]]=edge;
			/*||end||*/

		}
		/*||begin||*/
		//} /* end if owned edge */
		/*||end||*/
	}

	/* delete three edges of the fourth son */
	for (ied = 0; ied < nr_edges; ied++) {
		edge = abs(mesh->face[face_sons[4]].edge[ied]);
		iaux = mmr_del_edge(Mesh_id, edge);
		if (iaux <= 0) return (-1);
	}

	/* delete sons */
	for (ison = 1; ison <= nr_sons; ison++) {
		mmr_del_face(Mesh_id, face_sons[ison]);
	}

	/* change parameters for the father */
	mesh->parm.nrfa++;
	mesh->face[Face].type *= -1;
	free(mesh->face[Face].sons);
	mesh->face[Face].sons = NULL;

	/*kbw
	if(Face>0){
	printf("Clustered back TO face %d, type %d, bc %d\n",
	Face,mesh->face[Face].type,
	mesh->face[Face].bc);
	printf("edges: ");
	for(ied=0;ied<nr_edges;ied++) printf(" %d ",
	mesh->face[Face].edge[ied]);
	printf("\nneighbors %d %d\n",mesh->face[Face].neig[0],
	mesh->face[Face].neig[1]);
	}
	kew*/

	return (1);
}

/*------------------------------------------------------------
mmr_clust_fa4_q - to cluster back family of 4 quadrilateral faces
------------------------------------------------------------*/
int mmr_clust_fa4_q( /* returns: 1-success, <=0-failure */
		int  Mesh_id,	 /* in: ID of the mesh to be used or 0 for the current mesh */
		int Face,      /* in: face ID */
		int * Face_sons /* in (optional): face sons' IDs */
)
{


	/* local variables */
	mmt_mesh * mesh;
	int iaux, ied, edge, face_sons[5], ison;

	const int nr_edges = 4;
	const int nr_sons = 4;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh( Mesh_id);

	/* check whether face is inactive */
#ifdef DEBUG_MMM
	if (mmr_fa_status( Mesh_id, Face) != MMC_INACTIVE) {
		printf("not an inactive face %d in cluster face\n", Face);
		exit(-1);
	}
#endif

	/* find sons */
	if (Face_sons == NULL) mmr_fa_fam( Mesh_id, Face, face_sons, NULL);
	else {

		face_sons[0] = nr_sons;
		for (ison = 1; ison <= nr_sons; ison++) {
			face_sons[ison] = Face_sons[ison];
		}

	}

	/* check whether all sons are active */
#ifdef DEBUG_MMM
	for (ison = 1; ison <= nr_sons; ison++) {
		if (mmr_fa_status( Mesh_id, face_sons[ison]) != MMC_ACTIVE) {
			printf("not an active face son %d in cluster face %d\n",
					face_sons[ison], Face);
			exit(-1);
		}
	}
#endif

	/*kbw
	  if(Face==669||Face==739){
		printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Clustering back face %d, type %d, bc %d, owner %d\n",Face,mesh->face[Face].type, mesh->face[Face].bc, mmr_get_owner(Mesh_id,MMC_FACE,Face ));
	}
	/*kbw
	  if(Face==669||Face==739){
		for(ison=1;ison<=nr_sons;ison++){
		  printf("Clustering back face %d, type %d, bc %d\n",
			 face_sons[ison],mesh->face[face_sons[ison]].type,
			 mesh->face[face_sons[ison]].bc);
		  printf("edges: ");
		  for(ied=0;ied<nr_edges;ied++) printf(" %d ",
				 mesh->face[face_sons[ison]].edge[ied]);
		  printf("\nneighbors %d %d\n",mesh->face[face_sons[ison]].neig[0],
			 mesh->face[face_sons[ison]].neig[1]);
	  }
	}
	kew*/

	/* cluster back eight boundary edges of a family (if possible) */
	for (ied = 0; ied < nr_edges; ied++) {
		edge = abs(mesh->face[Face].edge[ied]);

		/*||begin||*/
		/* if edge is not owned */
		//if(mesh->edge[edge].ipid==0){
		/*||end||*/
		if (mesh->edge[edge].type < -1) {
			/* decrease the counter of attempted divisions */
			mesh->edge[edge].type++;

			/*kbw
			if(Face>0){
			printf("Attempted clustering of edge %d, type %d\n",
			edge,mesh->edge[edge].type);
			}
			kew*/

		}
		else {

			/* cluster back edge */
			iaux = mmr_clust_edge(Mesh_id, edge);
			if (iaux <= 0) return (-1);

			/*||begin||*/
			/* 	/\* put on the list of clustered edges *\/ */
			/* 	mmv_del_loc.edge[0]++; */
			/* 	mmv_del_loc.edge[mmv_del_loc.edge[0]]=edge; */
			/*||end||*/

		}
		/*||begin||*/
		//}
		/*||end||*/
	}

	/* delete node in the middle */
	edge = abs(mesh->face[face_sons[1]].edge[1]);
	iaux = mmr_del_node(Mesh_id, mesh->edge[edge].node[1]);
	if (iaux <= 0) return (-1);

	/* delete four internal edges */
	for (ied = 0; ied < nr_edges; ied++) {
		edge = abs(mesh->face[face_sons[ied + 1]].edge[(ied + 1) % nr_edges]);
		iaux = mmr_del_edge(Mesh_id, edge);
		if (iaux <= 0) return (-1);
	}

	/* delete sons */
	for (ison = 1; ison <= nr_sons; ison++) {
		mmr_del_face(Mesh_id, face_sons[ison]);
	}

	/* change parameters for the father */
	mesh->parm.nrfa++;
	mesh->face[Face].type *= -1;
	free(mesh->face[Face].sons);
	mesh->face[Face].sons = NULL;

	/*kbw
	if(Face>0){
	printf("Clustered back TO face %d, type %d, bc %d\n",
	Face,mesh->face[Face].type,
	zmesh->face[Face].bc);
	printf("edges: ");
	for(ied=0;ied<nr_edges;ied++) printf(" %d ",
	mesh->face[Face].edge[ied]);
	printf("\nneighbors %d %d\n",mesh->face[Face].neig[0],
	mesh->face[Face].neig[1]);
	}
	kew*/

	return (1);
}

/*------------------------------------------------------------
mmr_del_face    - to free a face structure
------------------------------------------------------------*/
int mmr_del_face( /* returns: 1-success, <=0-failure */
		int  Mesh_id,	/* in: ID of the mesh to be used or 0 for the current mesh */
		int Face      /* in: face ID */
)
{

	/* local variables */
	mmt_mesh * mesh;
	int face_type, nr_edges, ied;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh( Mesh_id);

	face_type = mmr_fa_type( Mesh_id, Face);
	if (face_type == MMC_TRIA) nr_edges = 3;
	else nr_edges = 4;

	/*kbw
	  if(Face==669||Face==739){
	printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Deleting face %d, type %d, bc %d, owner %d\n",
	Face,mesh->face[Face].type,
	mesh->face[Face].bc, mmr_get_owner(Mesh_id,MMC_FACE,Face ));
	printf("edges: ");
	for(ied=0;ied<nr_edges;ied++) printf(" %d ",
	mesh->face[Face].edge[ied]);
	printf("\nneighbors %d %d\n",mesh->face[Face].neig[0],
	mesh->face[Face].neig[1]);
	}
	kew*/

	if (mmr_fa_status(Mesh_id, Face) == MMC_ACTIVE) {
		/* decrease the number of active faces */
		mesh->parm.nrfa--;
	}

	mesh->face[Face].type = MMC_FREE;
	mesh->face[Face].bc = mesh->parm.pffa;
	for (ied = 0; ied < nr_edges; ied++) mesh->face[Face].edge[ied] = 0;
	mesh->face[Face].neig[0] = 0;
	mesh->face[Face].neig[1] = 0;

	/* put Face as the first available free space */
	mesh->parm.pffa = Face;

	/*kbw
	printf("Deleted face %d, nrfa %d, nmfa %d, pffa %d\n",
	Face, mesh->parm.nrfa, mesh->parm.nmfa, mesh->parm.pffa);
	kew*/

	return (1);
}

/*------------------------------------------------------------
mmr_clust_edge  - to cluster two edges back into their father
------------------------------------------------------------*/
int mmr_clust_edge( /* returns: 1-success, <=0-failure */
		int  Mesh_id,	/* in: ID of the mesh to be used or 0 for the current mesh */
		int Edge      /* in: father edge ID */
)
{

	/* local variables */
	mmt_mesh * mesh;
	int iaux, son1, son2;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh( Mesh_id);

	/* check whether Edge is inactive */
#ifdef DEBUG_MMM
	if (mesh->edge[Edge].type >= 0) {
		printf("not an inactive edge or free space %d in cluster edge\n", Edge);
		getchar();
		getchar();
		getchar();
		exit(-1);
	}
#endif

	/* check whether edge can be clustered */
	/*||begin||*/
	//if(mesh->edge[Edge].ipid==0&&mesh->edge[Edge].type<-1){
	/*||end||*/
#ifdef DEBUG_MMM
	if (mesh->edge[Edge].type < -1) {
		printf("edge %d cannot be clustered in cluster edge\n", Edge);
		return (0);
	}
#endif

	/* find sons, check whether they are active */
	son1 = mesh->edge[Edge].node[0];
	son2 = mesh->edge[Edge].node[1];

#ifdef DEBUG_MMM
	if (mesh->edge[son1].type <= 0) {
		printf("not an active edge son %d in cluster edge %d\n", son1, Edge);
		return (0);
	}

	if (mesh->edge[son2].type <= 0) {
		printf("not an active edge son %d in cluster edge %d\n", son2, Edge);
		return (0);
	}
#endif

	/*kbw
	printf("Clustering edge %d, type %d, sons %d %d\n",
		Edge,mesh->edge[Edge].type,
		mesh->edge[Edge].node[0],mesh->edge[Edge].node[1]);
	printf("Son1 %d, type %d, nodes %d %d\n",
		son1,mesh->edge[son1].type,
		mesh->edge[son1].node[0],mesh->edge[son1].node[1]);
	printf("Son2 %d, type %d, nodes %d %d\n",
		son2,mesh->edge[son2].type,
		mesh->edge[son2].node[0],mesh->edge[son2].node[1]);
	kew*/

	/* substitute nodes for father Edge */
	mesh->edge[Edge].node[0] = mesh->edge[son1].node[0];
	mesh->edge[Edge].node[1] = mesh->edge[son2].node[1];

	/* delete node between sons */
	iaux = mmr_del_node(Mesh_id, mesh->edge[son1].node[1]);
	if (iaux <= 0) return (-1);

	iaux = mmr_del_edge(Mesh_id, son1);
	if (iaux <= 0) return (-1);
	iaux = mmr_del_edge(Mesh_id, son2);
	if (iaux <= 0) return (-1);

	/* change parameters for father Edge */
	mesh->parm.nred++;
	mesh->edge[Edge].type = MMC_EDGE;

	/*kbw
	printf("Clustered TO edge %d, type %d, nodes %d %d\n",
		Edge,mesh->edge[Edge].type,
		mesh->edge[Edge].node[0],mesh->edge[Edge].node[1]);
	kew*/

	return (1);
}

/*------------------------------------------------------------
mmr_del_edge    - to free an edge structure
------------------------------------------------------------*/
int mmr_del_edge( /* returns: 1-success, <=0-failure */
		int  Mesh_id,	/* in: ID of the mesh to be used or 0 for the current mesh */
		int Edge      /* in: edge ID */
)
{

	/* local variables */
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh( Mesh_id);

	/*kbw
	printf("Deleting edge %d, type %d, nodes %d %d\n",
		Edge,mesh->edge[Edge].type,
		mesh->edge[Edge].node[0],mesh->edge[Edge].node[1]);
	kew*/

	if (mmr_edge_status(Mesh_id, Edge) == MMC_ACTIVE) {
		/* decrease the number of active edges */
		mesh->parm.nred--;
	}

	/* delete edge */
	mesh->edge[Edge].type = MMC_FREE;
	/* store pointer to next free space in node[0] */
	mesh->edge[Edge].node[0] = mesh->parm.pfed;
	mesh->edge[Edge].node[1] = 0;

	/* free edge[].elems structure */
#ifdef DEBUG_MM
	if (mesh->edge[Edge].elems == NULL) {
		printf("Error 24248 in mmr_del_edge (edge %d, elems==NULL). Exiting!\n",
				ned);
	}
#endif
	free(mesh->edge[Edge].elems);
	mesh->edge[Edge].elems = NULL;

	/* put Edge as the first available free space */
	mesh->parm.pfed = Edge;

	/*kbw
	printf("Deleted edge %d, nred %d, nmed %d, pfed %d\n",
	Edge, mesh->parm.nred, mesh->parm.nmed, mesh->parm.pfed);
	kew*/

	return (1);
}


/*------------------------------------------------------------
mmr_del_node    - to free a node structure
------------------------------------------------------------*/
int mmr_del_node( /* returns: 1-success, <=0-failure */
		int  Mesh_id,	/* in: ID of the mesh to be used or 0 for the current mesh */
		int Node      /* in: node ID */
)
{

	/* local variables */
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh( Mesh_id);

	/*kbw
	printf("Deleting node %d, coord: x=%lf, y=%lf, z=%lf\n",
	Node,mesh->node[Node].x,
	mesh->node[Node].y,mesh->node[Node].z);
	kew*/

	/* delete node */
	mesh->node[Node].x = -1e11;
	/* store pointer to next free space in y-coordinate */
	mesh->node[Node].y = mesh->parm.pfno;
	mesh->node[Node].z = 0.0;
	/* put Node as the first available free space */
	mesh->parm.pfno = Node;
	/* decrease the number of active nodes */
	mesh->parm.nrno--;

	/*kbw
	printf("Deleted node %d, nrno %d, nmno %d, pfno %d\n",
	Node, mesh->parm.nrno, mesh->parm.nmno, mesh->parm.pfno);
	kew*/

	return (1);
}


/*------------------------------------------------------------
mmr_del_elem    - to free an element structure
------------------------------------------------------------*/
int mmr_del_elem( /* returns: 1-success, <=0-failure */
		int  Mesh_id,	/* in: ID of the mesh to be used or 0 for the current mesh */
		int Elem      /* in: element ID */
)
{

	/* local variables */
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh( Mesh_id);

	/*kbw
	#ifdef DEBUG_MMM
	  if(mmr_el_status(Mesh_id,Elem)==MMC_ACTIVE){
		printf("deleting active elem %d, nr_elem %d\n",
		   Elem, mmr_get_nr_elem(Mesh_id));
	  }
	  else if(mmr_el_status(Mesh_id,Elem)==MMC_INACTIVE){
		printf("deleting inactive elem %d, nr_elem %d\n",
		   Elem, mmr_get_nr_elem(Mesh_id));
	  }
	  else{
		printf("deleting free space %d, error!!!!\n",
		   Elem);
		exit(-1);
	  }
	#endif
	kew*/

	if (mmr_el_status(Mesh_id, Elem) == MMC_ACTIVE) {
		/* decrease the number of active elements */
		mesh->parm.nrel--;
	}

	mesh->elem[Elem].type = MMC_FREE;
	mesh->elem[Elem].mate = mesh->parm.pfel;
	/* put Elem as the first available free space */
	mesh->parm.pfel = Elem;

	/* free not necessary storage space */
	free(mesh->elem[Elem].face);
	if (mesh->elem[Elem].sons != NULL) {
		free(mesh->elem[Elem].sons);
		mesh->elem[Elem].sons = NULL;
	}

	/* auxiliary substitutions */
	mesh->elem[Elem].refi = MMC_NOT_REF;
	mesh->elem[Elem].fath = MMC_NO_FATH;

	return (1);
}

/* internal utility procedures */

/*---------------------------------------------------------
mut_mat3vec - to compute matrix vector product in 3D space
---------------------------------------------------------*/
void mut_mat3vec(
		double * m1, 	/* in: matrix (stored by rows as a vector!) */
		double * v1, 	/* in: vector */
		double * v2	/* out: resulting vector */
)
{

	v2[0] = m1[0] * v1[0] + m1[1] * v1[1] + m1[2] * v1[2] ;
	v2[1] = m1[3] * v1[0] + m1[4] * v1[1] + m1[5] * v1[2] ;
	v2[2] = m1[6] * v1[0] + m1[7] * v1[1] + m1[8] * v1[2] ;
	return;
}

/*---------------------------------------------------------
mut_mat3mat - to compute matrix matrix product in 3D space
	(all matrices are stored by rows as vectors!)
---------------------------------------------------------*/
void mut_mat3mat(
		double * m1,	/* in: matrix */
		double * m2,	/* in: matrix */
		double * m3	/* out: matrix m1*m2 */
)
{

	m3[0] = m1[0] * m2[0] + m1[1] * m2[3] + m1[2] * m2[6] ;
	m3[1] = m1[0] * m2[1] + m1[1] * m2[4] + m1[2] * m2[7] ;
	m3[2] = m1[0] * m2[2] + m1[1] * m2[5] + m1[2] * m2[8] ;

	m3[3] = m1[3] * m2[0] + m1[4] * m2[3] + m1[5] * m2[6] ;
	m3[4] = m1[3] * m2[1] + m1[4] * m2[4] + m1[5] * m2[7] ;
	m3[5] = m1[3] * m2[2] + m1[4] * m2[5] + m1[5] * m2[8] ;

	m3[6] = m1[6] * m2[0] + m1[7] * m2[3] + m1[8] * m2[6] ;
	m3[7] = m1[6] * m2[1] + m1[7] * m2[4] + m1[8] * m2[7] ;
	m3[8] = m1[6] * m2[2] + m1[7] * m2[5] + m1[8] * m2[8] ;

}

/*---------------------------------------------------------
mut_chk_list - to check whether Num is on the list List
	with length Ll
---------------------------------------------------------*/
int mut_chk_list(	/* returns: */
		/* >0 - position on the list */
		/* 0 - not found on the list */
		int Num, 	/* number to be checked */
		int * List, 	/* list of numbers */
		int Ll		/* length of the list */
)
{

	int i, il;

	for (i = 0; i < Ll; i++) {
		if ((il = List[i]) == 0) break;
		/* found on the list on (i+1) position */
		if (Num == il) return (i + 1);
	}
	/* not found on the list */
	return (0);
}

