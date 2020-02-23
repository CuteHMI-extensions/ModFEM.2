/************************************************************************
File mms_prism_ref.c - routines for breaking and clustering prismatic
					   elements

Contains routines:
  mmr_divide_el8_p - to break a prismatic element into 8 sons
  mmr_clust_el8_p - to cluster back a family of eight prisms

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

/*------------------------------------------------------------
  mmr_divide_el8_p - to break a prismatic element into 8 sons
------------------------------------------------------------*/
int mmr_divide_el8_p( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id,    /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int El        /* in: element ID */
)
{

	/* auxiliary variables */
	mmt_mesh * mesh;
//	int el_node[MMC_MAXELVNO + 1];      /* list of nodes of El */
	int sons_edges[16], sons_new_nodes[16], face_sons[5];
	int new_nodes[5], new_edges[3], orient[20];
	int hor_edges[6], vert_edges[6], base_edges[6], edges[4];
	int hor_faces[4], vert_faces[6], out_hface[8], out_vface[12];
	int i, iaux, ifa, ison, face, ifaneig, shift, pel = 5; //, iprint = 5;
	/*||begin||*/
	//int new_ipid;
	/*||end||*/

	const int num_face = 5;
	const int num_sons = 8;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	/*kbw
	  if(El==8){
		int ino,num_dof_scal,init=0;
		double dofs_fath[10], dofs_fath_2[10];
		printf("Dividing element %d, type %d, mate %d, fath %d, refi %d\n",
		   El,mesh->elem[El].type,mesh->elem[El].mate,mesh->elem[El].fath,
		   mesh->elem[El].refi);
		printf("faces: ");
		for(ino=0;ino<5;ino++) printf(" %d ",
					  mesh->elem[El].face[ino]);
		printf("\nnodes: ");
		mmr_el_node_coor(Mesh_id,El,el_node,NULL);
		for(ino=0;ino<el_node[0];ino++) printf(" %d ",el_node[ino+1]);
		printf("\n");
	//printf("dofs_1: \n");
	//for(ino=0;ino<num_dof_scal;ino++) printf("%20.15lf",dofs_fath[ino]);
	//printf("\n");
	//if(init>1){
	//  printf("dofs_2: ");
	//  for(ino=0;ino<num_dof_scal;ino++) printf("%20.15lf",dofs_fath_2[ino]);
	//  printf("\n");
	//}
		getchar();
	  }
	kew*/


	/*||begin||*/
	/* prepare ipid - inter-processor ID for newly created entities */
	/*   if(mesh->elem[El].ipid==0){ */
	/*     new_ipid = 0; */
	/*     mmv_ref_loc.elem[0]++; */
	/*     mmv_ref_loc.elem[mmv_ref_loc.elem[0]]=El; */
	/*     mmv_ref_loc.elem[0]++; */
	/*     mmv_ref_loc.elem[mmv_ref_loc.elem[0]]=mmr_el_type(Mesh_id, El);; */
	/*   } */
	/*   else if(mesh->elem[El].ipid>0){ */
	/*     new_ipid = -1; */
	/*     mmv_ref_ali.elem[0]++; */
	/*     mmv_ref_ali.elem[mmv_ref_ali.elem[0]]=El; */
	/*     mmv_ref_ali.elem[0]++; */
	/*     mmv_ref_ali.elem[mmv_ref_ali.elem[0]]=mmr_el_type(Mesh_id, El);; */
	/*   } */
	/*   else{ */
	/*     printf("external mesh entity %d in divide_el\n",El); */
	/*     exit(-1); */
	/*   } */
	/*||end||*/


	/* divide faces, prepare list of faces for new elements */
	for (ifa = 0; ifa < num_face; ifa++) {

		face = mesh->elem[El].face[ifa];

		if (mmr_fa_status(Mesh_id, abs(face)) == MMC_ACTIVE) {

			/*kbw
			  if(abs(face)==669||abs(face)==739){
				printf("dividing elem %d, ACTIVE face %d\n", El, abs(face));
			  }
			kew*/

			/* active face: we have to divide it and create new faces */
			if (ifa < 2) {

				iaux = mmr_divide_face4_t(Mesh_id, abs(face), face_sons,
								sons_edges, sons_new_nodes);
				//, NULL);

				/* the same orientation for the element and the face */
				if (face > 0) {

					if (ifa == 0) {
						/* put face sons on a list of future elements' faces */
						out_hface[0]  =  face_sons[0];
						out_hface[1]  =  face_sons[2];
						out_hface[2]  =  face_sons[1];
						out_hface[3]  =  face_sons[3];
						/* put face edges on a list of future faces' edges */
						base_edges[0] = -sons_edges[10];
						base_edges[1] = -sons_edges[9];
						base_edges[2] = -sons_edges[11];
					}
					else {
						/* put face sons on a list of future elements' faces */
						out_hface[4]  =  face_sons[0];
						out_hface[5]  =  face_sons[1];
						out_hface[6]  =  face_sons[2];
						out_hface[7]  =  face_sons[3];
						/* put face edges on a list of future faces' edges */
						base_edges[3] = sons_edges[10];
						base_edges[4] = sons_edges[11];
						base_edges[5] = sons_edges[9];
					}
				}
				/* the opposite orientation for the element and the face */
				else if (face < 0) {
					shift = -mesh->face[-face].bc;
					if (ifa == 0) {
						/* put face sons on a list of future elements' faces */
						out_hface[0]  = -face_sons[shift];
						out_hface[1]  = -face_sons[(shift + 1) % 3];
						out_hface[2]  = -face_sons[(shift + 2) % 3];
						out_hface[3]  = -face_sons[3];
						/* put face edges on a list of future faces' edges */
						base_edges[0] = sons_edges[9 + (shift + 1) % 3];
						base_edges[1] = sons_edges[9 + (shift + 2) % 3];
						base_edges[2] = sons_edges[9 + shift];
					}
					else {
						/* put face sons on a list of future elements' faces */
						out_hface[4]  = -face_sons[shift];
						out_hface[5]  = -face_sons[(shift + 2) % 3];
						out_hface[6]  = -face_sons[(shift + 1) % 3];
						out_hface[7]  = -face_sons[3];
						/* put face edges on a list of future faces' edges */
						base_edges[3] = -sons_edges[9 + (shift + 1) % 3];
						base_edges[4] = -sons_edges[9 + shift];
						base_edges[5] = -sons_edges[9 + (shift + 2) % 3];
					}
				}
			}
			else {

				iaux = mmr_divide_face4_q(Mesh_id, abs(face), face_sons,
								sons_edges, sons_new_nodes);
				//, NULL);

				new_nodes[ifa] = sons_new_nodes[4];
				/* collect necessary horizontal and vertical faces */
				/* the same orientation for the element and the face */
				if (face > 0) {
					/* put face sons on a list of future elements' faces */
					out_vface[4 * (ifa - 2)]    =  face_sons[0];
					out_vface[4 * (ifa - 2) + 1]  =  face_sons[1];
					out_vface[4 * (ifa - 2) + 2]  =  face_sons[2];
					out_vface[4 * (ifa - 2) + 3]  =  face_sons[3];
					/* put face edges on a list of future faces' edges */
					/* a bit of magic */
					hor_edges[2 * (ifa - 2)]    = -sons_edges[2];
					hor_edges[2 * (ifa - 2) + 1]  =  sons_edges[8];
					vert_edges[2 * (ifa - 2)]   =  sons_edges[1];
					vert_edges[2 * (ifa - 2) + 1] = -sons_edges[11];
				}
				/* the opposite orientation for the element and the face */
				else if (face < 0) {
					shift = -mesh->face[-face].bc;
					if (shift == 1) {
						/* put face sons on a list of future elements' faces */
						out_vface[4 * (ifa - 2)]    = -face_sons[1];
						out_vface[4 * (ifa - 2) + 1]  = -face_sons[0];
						out_vface[4 * (ifa - 2) + 2]  = -face_sons[3];
						out_vface[4 * (ifa - 2) + 3]  = -face_sons[2];
						/* put face edges on a list of future faces' edges */
						/* a bit of magic */
						hor_edges[2 * (ifa - 2)]    = -sons_edges[8];
						hor_edges[2 * (ifa - 2) + 1]  =  sons_edges[2];
						vert_edges[2 * (ifa - 2)]   =  sons_edges[1];
						vert_edges[2 * (ifa - 2) + 1] = -sons_edges[11];
					}
					else {
#ifdef DEBUG_MMM
						printf("Face_shift!=1 for vertical, quadrilateral face %d !\n", -face);
						printf("Use half structured meshes with shift==1 or modify\n");
						printf("mesh modification routines.\n");
#endif
					}
				}
			}

			/*kbw
				  if(El==8){
				if(ifa<2){
				  printf("For element %d divided face %d\n",
					 El, abs(face));
				  printf("sons:");
				  for(i=0;i<4;i++) printf("  %d", face_sons[i]);
				  printf("\n");
				  printf("sons' edges:");
				  for(i=0;i<12;i++) printf("  %d", sons_edges[i]);
				  printf("\n");
				  printf("new nodes:");
				  for(i=0;i<4;i++) printf("  %d", sons_new_nodes[i]);
				  printf("\n");
				}
				else{
				  printf("For element %d divided face %d, (node in the middle %d)\n",
					 El, abs(face), new_nodes[ifa]);
				  printf("sons:");
				  for(i=0;i<4;i++) printf("  %d", face_sons[i]);
				  printf("\n");
				  printf("sons' edges:");
				  for(i=0;i<16;i++) printf("  %d", sons_edges[i]);
				  printf("\n");
				  printf("new nodes:");
				  for(i=0;i<5;i++) printf("  %d", sons_new_nodes[i]);
				  printf("\n");
				}
				  }
			kew*/

#ifdef DEBUG_MMM
			if (iaux <= 0) {
				printf("error in dividing faces\n");
				exit(-1);
			}
#endif

			/* update the other neighbors for new faces */
			if (face > 0) {

				/* if orientation the same - El is the first neighbor */

				if (mesh->face[abs(face)].neig[1] == MMC_BOUNDARY) {
					/* if the second neighbor is the boundary */
					for (i = 0; i < 4; i++) mesh->face[face_sons[i]].neig[1] = MMC_BOUNDARY;
				}
				else if (mesh->face[abs(face)].neig[1] == MMC_SUB_BND) {
					/* if the second neighbor is the intersubdomain boundary */
					for (i = 0; i < 4; i++) mesh->face[face_sons[i]].neig[1] = MMC_SUB_BND;
				}
				else {
					for (i = 0; i < 4; i++) mesh->face[face_sons[i]].neig[1] = MMC_BIG_NGB;
				}

			}
			else {

				/* for opposite orientation - El is the second neighbor */

				if (mesh->face[abs(face)].neig[0] == MMC_SUB_BND) {
					/* if the second neighbor is the intersubdomain boundary */
					for (i = 0; i < 4; i++) mesh->face[face_sons[i]].neig[0] = MMC_SUB_BND;
				}
				else {
					for (i = 0; i < 4; i++) mesh->face[face_sons[i]].neig[0] = MMC_BIG_NGB;
				}

			}

		} /* end if active face */
		else if (mmr_fa_status(Mesh_id, abs(face)) == MMC_INACTIVE) {

			/*kbw
			  if(abs(face)==669||abs(face)==739){
				printf("dividing elem %d, INACTIVE face %d\n", El, abs(face));
			  }
			kew*/

			/* inactive face - we already have its sons */
			mmr_fa_fam(Mesh_id, abs(face), face_sons, &new_nodes[ifa]);

			if (ifa < 2) {
				/* get necessary edges of sons */

				mmr_fa_edges(Mesh_id, face_sons[4], sons_edges, orient);

				if (face > 0) {
					if (ifa == 0) {
						/* put face sons on a list of future elements' faces */
						out_hface[0]  =  face_sons[1];
						out_hface[1]  =  face_sons[3];
						out_hface[2]  =  face_sons[2];
						out_hface[3]  =  face_sons[4];
						/* put face edges on a list of future faces' edges */
						base_edges[0] = -orient[1] * sons_edges[1];
						base_edges[1] = -orient[0] * sons_edges[0];
						base_edges[2] = -orient[2] * sons_edges[2];
					}
					else {
						/* put face sons on a list of future elements' faces */
						out_hface[4]  =  face_sons[1];
						out_hface[5]  =  face_sons[2];
						out_hface[6]  =  face_sons[3];
						out_hface[7]  =  face_sons[4];
						/* put face edges on a list of future faces' edges */
						base_edges[3] = orient[1] * sons_edges[1];
						base_edges[4] = orient[2] * sons_edges[2];
						base_edges[5] = orient[0] * sons_edges[0];
					}
				}
				else if (face < 0) {
					shift = -mesh->face[-face].bc;
					if (ifa == 0) {
						/* put face sons on a list of future elements' faces */
						out_hface[0]  = -face_sons[shift + 1];
						out_hface[1]  = -face_sons[(shift + 1) % 3 + 1];
						out_hface[2]  = -face_sons[(shift + 2) % 3 + 1];
						out_hface[3]  = -face_sons[4];
						/* put face edges on a list of future faces' edges */
						base_edges[0] = orient[(shift + 1) % 3] * sons_edges[(shift + 1) % 3];
						base_edges[1] = orient[(shift + 2) % 3] * sons_edges[(shift + 2) % 3];
						base_edges[2] = orient[shift] * sons_edges[shift];
					}
					else {
						/* put face sons on a list of future elements' faces */
						out_hface[4]  = -face_sons[shift + 1];
						out_hface[5]  = -face_sons[(shift + 2) % 3 + 1];
						out_hface[6]  = -face_sons[(shift + 1) % 3 + 1];
						out_hface[7]  = -face_sons[4];
						/* put face edges on a list of future faces' edges */
						base_edges[3] = -orient[(shift + 1) % 3] * sons_edges[(shift + 1) % 3];
						base_edges[4] = -orient[shift] * sons_edges[shift];
						base_edges[5] = -orient[(shift + 2) % 3] * sons_edges[(shift + 2) % 3];
					}
				}
			}
			else if (ifa >= 2) {
				if (face > 0) {
					/* put face sons on a list of future elements' faces */
					out_vface[4 * (ifa - 2)]    =  face_sons[0 + 1];
					out_vface[4 * (ifa - 2) + 1]  =  face_sons[1 + 1];
					out_vface[4 * (ifa - 2) + 2]  =  face_sons[2 + 1];
					out_vface[4 * (ifa - 2) + 3]  =  face_sons[3 + 1];
				}
				else if (face < 0) {
					shift = -mesh->face[-face].bc;
					if (shift == 1) {
						/* put face sons on a list of future elements' faces */
						out_vface[4 * (ifa - 2)]    = -face_sons[1 + 1];
						out_vface[4 * (ifa - 2) + 1]  = -face_sons[0 + 1];
						out_vface[4 * (ifa - 2) + 2]  = -face_sons[3 + 1];
						out_vface[4 * (ifa - 2) + 3]  = -face_sons[2 + 1];
					}
					else {
#ifdef DEBUG_MMM
						printf("Face_shift!=1 for vertical, quadrilateral face %d !\n", -face);
						printf("Use half structured meshes with shift==1 or modify\n");
						printf("mesh modification routines.\n");
#endif
					}
				}

				mmr_fa_edges(Mesh_id, face_sons[1], sons_edges, orient);

				/*kbw
				if(El==8){
				  printf(" face %d, son0 %d\n", abs(face), face_sons[1]);
				  printf("sons' edges:");
				  for(i=0;i<4;i++) printf("  %d", orient[i]*sons_edges[i]);
				  printf("\n");
				}
				kew*/

				if (face > 0) {
					/* a bit of magic */
					hor_edges[2 * (ifa - 2)]    = -orient[2] * sons_edges[2];
					vert_edges[2 * (ifa - 2)]   =  orient[1] * sons_edges[1];
				}
				else if (face < 0) {
					shift = -mesh->face[-face].bc;
					/* a bit of magic */
					if (shift == 1) {
						hor_edges[2 * (ifa - 2) + 1]  =  orient[2] * sons_edges[2];
						vert_edges[2 * (ifa - 2)]   =  orient[1] * sons_edges[1];
					}
					else {
#ifdef DEBUG_MMM
						printf("Face_shift!=1 for vertical, quadrilateral face %d !\n", -face);
						printf("Use half structured meshes with shift==1 or modify\n");
						printf("mesh modification routines.\n");
#endif
					}
				}

				mmr_fa_edges(Mesh_id, face_sons[3], sons_edges, orient);

				/*kbw
					if(El==8){
					  printf(" face %d, son2 %d\n", abs(face), face_sons[3]);
					  printf("sons' edges:");
					  for(i=0;i<4;i++) printf("  %d", orient[i]*sons_edges[i]);
					  printf("\n");
					}
				kew*/

				if (face > 0) {
					/* a bit of magic */
					hor_edges[2 * (ifa - 2) + 1]  =  orient[0] * sons_edges[0];
					vert_edges[2 * (ifa - 2) + 1] = -orient[3] * sons_edges[3];
				}
				else if (face < 0) {
					shift = -mesh->face[-face].bc;
					/* a bit of magic */
					if (shift == 1) {
						hor_edges[2 * (ifa - 2)]    = -orient[0] * sons_edges[0];
						vert_edges[2 * (ifa - 2) + 1] = -orient[3] * sons_edges[3];
					}
					else {
#ifdef DEBUG_MMM
						printf("Face_shift!=1 for vertical, quadrilateral face %d !\n", -face);
						printf("Use half structured meshes with shift==1 or modify\n");
						printf("mesh modification routines.\n");
#endif
					}
				}
			}


			/*kbw
				  if(El==8){
				printf("For element %d found inactive face %d, (node in the middle %d)\n",
					   El, abs(face), new_nodes[ifa]);
				printf("sons with nodes:\n");
				for(i=0;i<4;i++) {
				  mmr_fa_node_coor(Mesh_id,face_sons[i+1],el_node,NULL);
				  printf("  %d:  %d  %d  %d  %d", face_sons[i+1], el_node[1]
					 , el_node[2], el_node[3], el_node[4]);
				}
				printf("\n");
			}
			kew*/

		}
		else {
			printf("break8_p - error in element faces\n");
			exit(-1);
		}


	} /* end loop over faces: ifa */

	/*kbw
		if(El==8){
		  printf("base edges:");
		  for(i=0;i<6;i++) printf("  %d", base_edges[i]);
		  printf("\n");
		  printf("horizontal edges:");
		  for(i=0;i<6;i++) printf("  %d", hor_edges[i]);
		  printf("\n");
		  printf("vertical edges:");
		  for(i=0;i<6;i++) printf("  %d", vert_edges[i]);
		  printf("\n");
		  printf("outer horizontal faces:");
		  for(i=0;i<8;i++) printf("  %d", out_hface[i]);
		  printf("\n");
		  printf("outer vertical faces:");
		  for(i=0;i<12;i++) printf("  %d", out_vface[i]);
		  printf("\n");
		}
	kew*/

	/* create new edges inside old element */
	for (i = 0; i < 3; i++) {

		new_edges[i] = mmr_create_edge(Mesh_id, MMC_EDGE, new_nodes[(i + 2) % 3 + 2],
						new_nodes[i + 2]);
		//, new_ipid);

		/*||begin||*/
		/*     if(new_ipid==0){ */
		/*       mmv_ref_loc.elem[0]++; */
		/*       mmv_ref_loc.elem[mmv_ref_loc.elem[0]]=new_edges[i]; */
		/*     } */
		/*     else if(new_ipid== -1){ */
		/*       mmv_ref_ali.elem[0]++; */
		/*       mmv_ref_ali.elem[mmv_ref_ali.elem[0]]=new_edges[i]; */
		/*     } */
		/*||end||*/

	}

	/* create new faces inside old element */
	/* horizontal faces */
	for (ison = 0; ison < 3; ison++) {

		edges[ison] =  hor_edges[2 * ison];
		edges[(ison + 1) % 3] = -new_edges[ison];
		edges[(ison + 2) % 3] =  hor_edges[2 * ((ison + 2) % 3) + 1];

		iaux = 0; /* no shift between upper and lower element */
		hor_faces[ison] = mmr_create_face(Mesh_id, MMC_TRIA, iaux, edges,
						NULL, NULL);
		//, new_ipid);

		/*||begin||*/
		/*     if(new_ipid==0){ */
		/*       mmv_ref_loc.elem[0]++; */
		/*       mmv_ref_loc.elem[mmv_ref_loc.elem[0]]=hor_faces[ison]; */
		/*     } */
		/*     else if(new_ipid== -1){ */
		/*       mmv_ref_ali.elem[0]++; */
		/*       mmv_ref_ali.elem[mmv_ref_ali.elem[0]]=hor_faces[ison]; */
		/*     } */
		/*||end||*/

	}

	/* fourth face */
	iaux = 0; /* no shift between upper and lower element */
	edges[0] = new_edges[2];
	edges[1] = new_edges[0];
	edges[2] = new_edges[1];
	hor_faces[3] = mmr_create_face(Mesh_id, MMC_TRIA, iaux, edges,
					NULL, NULL);
	//, new_ipid);

	/*||begin||*/
	/*   if(new_ipid==0){ */
	/*     mmv_ref_loc.elem[0]++; */
	/*     mmv_ref_loc.elem[mmv_ref_loc.elem[0]]=hor_faces[3]; */
	/*   } */
	/*   else if(new_ipid== -1){ */
	/*     mmv_ref_ali.elem[0]++; */
	/*     mmv_ref_ali.elem[mmv_ref_ali.elem[0]]=hor_faces[3]; */
	/*   } */
	/*||end||*/

	for (ifa = 0; ifa < 3; ifa++) {
		/* vertical faces */
		iaux = -1; /* standard half structured mesh - shift between faces -1 */
		edges[0] =  base_edges[ifa];
		edges[1] =  vert_edges[2 * ifa];
		edges[2] = -new_edges[ifa];
		edges[3] = -vert_edges[2 * ((ifa + 2) % 3)];
		vert_faces[2 * ifa] = mmr_create_face(Mesh_id, MMC_QUAD, iaux, edges,
						NULL, NULL);
		//, new_ipid);

		iaux = -1; /* standard half structured mesh - shift between faces -1 */
		edges[0] =  new_edges[ifa];
		edges[1] =  vert_edges[2 * ifa + 1];
		edges[2] = -base_edges[ifa + 3];
		edges[3] = -vert_edges[2 * ((ifa + 2) % 3) + 1];
		vert_faces[2 * ifa + 1] = mmr_create_face(Mesh_id, MMC_QUAD, iaux, edges,
						NULL, NULL);
		//, new_ipid);

		/*||begin||*/
		/*     if(new_ipid==0){ */
		/*       mmv_ref_loc.elem[0]++; */
		/*       mmv_ref_loc.elem[mmv_ref_loc.elem[0]]=vert_faces[2*ifa]; */
		/*       mmv_ref_loc.elem[0]++; */
		/*       mmv_ref_loc.elem[mmv_ref_loc.elem[0]]=vert_faces[2*ifa+1]; */
		/*     } */
		/*     else if(new_ipid== -1){ */
		/*       mmv_ref_ali.elem[0]++; */
		/*       mmv_ref_ali.elem[mmv_ref_ali.elem[0]]=vert_faces[2*ifa]; */
		/*       mmv_ref_ali.elem[0]++; */
		/*       mmv_ref_ali.elem[mmv_ref_ali.elem[0]]=vert_faces[2*ifa+1]; */
		/*     } */
		/*||end||*/

	}

	/* create data structures for sons */
	mesh->elem[El].refi = MMC_REF_ISO; /* indicate isotropic refinement */
	mesh->elem[El].sons = mmr_ivector(num_sons, "sons in divide_el8_p");

	/* for each new element */
	for (ison = 0; ison < num_sons; ison++) {

		/*take as a pointer to its space value stored at pfel*/
		pel = mesh->parm.pfel;
		if (pel > mesh->parm.mxel) {
			printf("Not enough space for next element, Nrel= %d\n",
					mesh->parm.nrel);
			exit(-1);
		}

		/* if space already visited */
		if (pel < mesh->parm.nmel + 1) {

			/* small check - whether we have space for new element */
			if (mesh->elem[pel].type != MMC_FREE) {
				printf("BREAK8: old pointer %d points not into free space\n",
						pel);
				exit(-1);
			}

			/* subsitute for pfel value stored in mate */
			mesh->parm.pfel = mesh->elem[pel].mate;

		} /* end if: found free space */
		else {

			/*small check - whether we are at the end of available space */
			if (pel != mesh->parm.nmel + 1) {
				printf("break4 - error in pointers\n");
				exit(-1);
			}

			/* update control variables */
			mesh->parm.pfel++;
			mesh->parm.nmel++;

		}

		mesh->parm.nrel++;

		/* substitute parameters */
		mesh->elem[El].sons[ison] = pel;

		mesh->elem[pel].type = mesh->elem[El].type;
		mesh->elem[pel].mate = mesh->elem[El].mate;
		mesh->elem[pel].fath = El;
		mesh->elem[pel].refi = MMC_NOT_REF;
		mesh->elem[pel].sons = NULL;
		mesh->elem[pel].face = mmr_ivector(num_face, "element faces in read data");

		/*||begin||*/
		/*     mesh->elem[pel].ipid = new_ipid; */
		/*     if(new_ipid==0){ */
		/*       mmv_ref_loc.elem[0]++; */
		/*       mmv_ref_loc.elem[mmv_ref_loc.elem[0]]=pel; */
		/*     } */
		/*     else if(new_ipid== -1){ */
		/*       mmv_ref_ali.elem[0]++; */
		/*       mmv_ref_ali.elem[mmv_ref_ali.elem[0]]=pel; */
		/*     } */
		/*||end||*/


		/* face numbers */
		if (ison < 3) {
			/* lower base */
			mesh->elem[pel].face[0] = out_hface[ison];
			/* upper base */
			mesh->elem[pel].face[1] = hor_faces[ison];
			/* sides */
			mesh->elem[pel].face[2 + ison] = out_vface[4 * ison];
			mesh->elem[pel].face[2 + (ison + 1) % 3] = -vert_faces[2 * ison];
			mesh->elem[pel].face[2 + (ison + 2) % 3] = out_vface[4 * ((ison + 2) % 3) + 1];
		}
		else if (ison == 3) {
			/* lower base */
			mesh->elem[pel].face[0] = out_hface[ison];
			/* upper base */
			mesh->elem[pel].face[1] = hor_faces[3];
			/* sides */
			mesh->elem[pel].face[2] = vert_faces[4];
			mesh->elem[pel].face[3] = vert_faces[0];
			mesh->elem[pel].face[4] = vert_faces[2];
		}
		else if (ison < 7) {
			ifa = ison - 4;
			/* lower base */
			mesh->elem[pel].face[0] = -hor_faces[ifa];
			/* upper base */
			mesh->elem[pel].face[1] = out_hface[4 + ifa];
			/* sides */
			mesh->elem[pel].face[2 + ifa] = out_vface[4 * ifa + 3];
			mesh->elem[pel].face[2 + (ifa + 1) % 3] = -vert_faces[2 * ifa + 1];
			mesh->elem[pel].face[2 + (ifa + 2) % 3] = out_vface[4 * ((ifa + 2) % 3) + 2];
		}
		else if (ison == 7) {
			/* lower base */
			mesh->elem[pel].face[0] = -hor_faces[3];
			/* upper base */
			mesh->elem[pel].face[1] = out_hface[7];
			/* sides */
			mesh->elem[pel].face[2] = vert_faces[5];
			mesh->elem[pel].face[3] = vert_faces[1];
			mesh->elem[pel].face[4] = vert_faces[3];
		}

		/* update information on face neighbors */
		for (ifa = 0; ifa < num_face; ifa++) {

			face = mesh->elem[pel].face[ifa];
			/* if face number is positive pel is first neighbor */
			if (face > 0) ifaneig = 0;
			else ifaneig = 1;
			mesh->face[abs(face)].neig[ifaneig] = pel;

		} /* for all faces of new elements */

		/* check consistency */
#ifdef DEBUG_MMM
		mmr_el_node_coor(Mesh_id, pel, el_node, NULL);
#endif

		/*kbw
			//if(El==8){
			  int ino;
			  printf("son %d - element %d, type %d, mate %d, fath %d, refi %d\n",
				 ison,pel,mesh->elem[pel].type,mesh->elem[pel].mate,
				 mesh->elem[pel].fath, mesh->elem[pel].refi);
			  printf("faces: ");
			  for(ino=0;ino<5;ino++) printf(" %d ",
							mesh->elem[pel].face[ino]);
			  printf("\nnodes: ");
			  for(ino=0;ino<el_node[0];ino++) printf(" %d ",el_node[ino+1]);
			  printf("\n");
			  getchar();
			  //}
		kew*/

	} /* end loop over sons */

	/* update parameters for the father */
	mesh->parm.nrel--;
	mesh->elem[El].type *= -1;

	return (1);
}

/*------------------------------------------------------------
mmr_clust_el8_p - to cluster back a family of eight prisms
------------------------------------------------------------*/
int mmr_clust_el8_p( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id,    /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		int El          /* in: element ID */
)
{

	/* auxiliary variables */
	mmt_mesh * mesh;
	int el_node[MMC_MAXELVNO + 1];      /* list of nodes of El */
	int fath_sons[MMC_MAXELSONS];
	int face_sons[5];
	int i, j, iaux, jaux, ifa, ison, son, max_gen, max_gen_diff;
	int face, face_neig, fath, el_type, mate, ifason;

	const int iprint = MMC_PRINT_INFO;

	const int num_face = 5;
	const int num_sons = 8;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* first check whether element is active */
	if (mmr_el_status(Mesh_id, El) != MMC_ACTIVE) {
#ifdef DEBUG_MMM
		if (iprint > MMC_PRINT_ERRORS) {
			printf("not active element or free space %d in clust_el8_p!\n", El);
		}
#endif
		return (-2);
	}

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	//i=21; max_gen=mmr_get_mesh_i_params(mesh,i);
	max_gen = mmr_get_max_gen(Mesh_id);
	//i=22; max_gen_diff=mmr_get_mesh_i_params(mesh,i);
	max_gen_diff = mmr_get_max_gen_diff(Mesh_id);

	/* find father */
	fath = mesh->elem[El].fath;
	el_type = mmr_el_type(Mesh_id, fath);

	/* small check */
	if (fath == MMC_NO_FATH) {
#ifdef DEBUG_MMM
		printf("Attempt to clust initial mesh element; in clust_el8_p\n");
#endif
		return (-3);
	}

	/* check whether all sons are active and create list of sons */
	for (ison = 0; ison < num_sons; ison++) {
		son = mesh->elem[fath].sons[ison];
		fath_sons[ison] = son;

#ifdef DEBUG_MMM
		if (mmr_el_status(Mesh_id, son) == MMC_INACTIVE) {
			printf("not active brother %d of %d in clust_el8_p!\n", son, El);
			return (-1);
		}
		if (mesh->elem[son].sons != NULL) {
			printf("son %d of %d with sons in clust_el8_p!\n", son, El);
			return (-1);
		}
		/* find material number */
		if (ison == 0) mate = mesh->elem[son].mate;
		else {
			if (mate != mesh->elem[son].mate) {
				printf("Clustering elements with different material types\n");
				return (-1);
			}
		}
		if (mate != mesh->elem[fath].mate) {
			printf("Clustering elements with different material type than father\n");
			return (-1);
		}
#endif


		/* check whether derefinement will not cause too big generation difference */

		/* simple test for max_gen_diff==1 */
		if (max_gen_diff == 1) {
			for (ifa = 0; ifa < num_face; ifa++) {
				face = abs(mesh->elem[son].face[ifa]);

				/*kbw
				  if(abs(face)==669||abs(face)==739){
						  printf("-father %d, son %d, face %d, status %d, type %d, bc %d\n",
							 fath, son, face, mmr_fa_status(Mesh_id,face),
							 mmr_fa_type(Mesh_id,face), mmr_fa_bc(Mesh_id,face));
				  }
				kew*/

				if (mmr_fa_status(Mesh_id, face) == MMC_INACTIVE) {

					if (mmr_fa_sub_bnd(Mesh_id, face)) {

						int ifason;

#ifdef DEBUG_MMM
						printf("inactive face %d on int.sub.bnd. - for clustering\n",
								face);
						getchar();
						getchar();
						getchar();
						getchar();
						getchar();
						getchar();
#endif

						mmr_fa_fam( Mesh_id, face, face_sons, NULL);

#ifdef DEBUG_MMM
						for (ifason = 1; ifason <= face_sons[0]; ifason++) {
							printf("son %d, status %d\n",
									face_sons[ifason], mmr_fa_status( Mesh_id, face_sons[ifason]));
						}
#endif
						if (mmr_fa_type(Mesh_id, face) == MMC_TRIA) {
							mmr_clust_fa4_t(Mesh_id, face, NULL);
						}
						else mmr_clust_fa4_q(Mesh_id, face, NULL);

					}
					else {

#ifdef DEBUG_MMM
						if (iprint > MMC_PRINT_ERRORS)
							printf("In deref: too big generation difference for element %d across the face %d!\n",
									son, abs(face));
#endif
						return (-1);

					}
				}

			}

		}
		else {

#ifdef DEBUG_MMM
			printf("clust_el_8p: max_gen_diff barrier not implemented for specified case!\n");
#endif

		}

	}

	/*kbw
	if(El==200){
	  int ino;
	  for(ison=0;ison<num_sons;ison++){
		son=mesh->elem[fath].sons[ison];
	printf("Clustering back element %d, type %d, mate %d fath %d, refi %d\n",
	son,mesh->elem[son].type,mesh->elem[son].mate,mesh->elem[son].fath,
	mesh->elem[son].refi);
	printf("faces: ");
	for(ino=0;ino<5;ino++) printf(" %d ",
	mesh->elem[son].face[ino]);
	printf("\nnodes: ");
	mmr_el_node_coor(Mesh_id,son,el_node,NULL);
	for(ino=0;ino<el_node[0];ino++) printf(" %d ",el_node[ino+1]);
	printf("\n");
	  }
	}
	kew*/

	/* checking */
#ifdef DEBUG_MMM
	if (mesh->elem[fath].refi != MMC_REF_ISO) {
		printf("attempt to cluster back isotropically element %d refined anisotropically\n", fath);
		exit(-1);
	}
#endif

	/*kbw
	if(El==200){
	printf("Clustering back TO element %d, type %d, mate %d, fath %d, refi %d\n",
	fath,mesh->elem[fath].type,mesh->elem[fath].mate,mesh->elem[fath].fath,
	mesh->elem[fath].refi);
	printf("faces: ");
	for(ino=0;ino<5;ino++) printf(" %d ",
	mesh->elem[fath].face[ino]);
	printf("\nnodes: ");
	mmr_el_node_coor(Mesh_id,fath,el_node,NULL);
	for(ino=0;ino<el_node[0];ino++) printf(" %d ",el_node[ino+1]);
	printf("\n");
	printf("dofs_1: \n");
	for(ino=0;ino<num_dof_scal;ino++) printf("%20.15lf",dofs_fath[ino]);
	printf("\n");
	if(init>1){
	  printf("dofs_2: ");
	  for(ino=0;ino<num_dof_scal;ino++) printf("%20.15lf",dofs_fath_2[ino]);
	  printf("\n");
	}
	}
	kew*/

	/* cluster back horizontal and vertical boundary faces of a family */

	/* for each face of big element */
	for (ifa = 0; ifa < num_face; ifa++) {
		face = abs(mesh->elem[fath].face[ifa]);

		if (mesh->face[face].neig[0] == fath)
			face_neig = mesh->face[face].neig[1];
		else face_neig = mesh->face[face].neig[0];


		/* if face lies on the intersubdoamin boundary */
		if (face_neig == MMC_SUB_BND) {

			/*kbw
				  if(face==669||face==739){
				printf("clustering back elem %d, face %d on inters. bound., owner %d \n",
					   fath, face, mmr_get_owner(Mesh_id,MMC_FACE,face));
				  }
			kew*/

			/* just in case - create list of sons
			mesh->face[face].sons = (int *) malloc(4*sizeof(int));
			mmr_fa_fam(Mesh_id, face, face_sons, NULL);
			for(i=1;i<=face_sons[0];i++){
			mesh->face[face].sons[i-1]=face_sons[i];
				 }
				 */

			/* checking */
#ifdef DEBUG_MMM
			/*       if(mmr_fa_ipid(Mesh_id, face)==0){ */
			/* 	printf("owned face %d on intersubdomain boundary\n",fath);  */
			/* 	exit(-1); */
			/*       } */
#endif

			/* it will be checked later whether this face should be clustered */

		}
		else {
			/* if neighbor across face is active or big or boundary */
			if (face_neig <= 0 ||
					mmr_el_status(Mesh_id, face_neig) == MMC_ACTIVE) {

				/*kbw
					if(abs(face)==669||abs(face)==739){
					  printf("clustering back elem %d, face %d, ACTIVE neighbor %d\n",
						 fath, abs(face), face_neig );
					}
				kew*/

				if (ifa < 2) mmr_clust_fa4_t(Mesh_id, face, NULL);
				else mmr_clust_fa4_q(Mesh_id, face, NULL);

				/*||begin||*/
				/* if face is owned */
				/* 	if(mesh->face[face].ipid==0){ */
				/* 	  /\* put on the list of clustered faces *\/ */
				/* 	  mmv_del_loc.face[0]++; */
				/* 	  mmv_del_loc.face[mmv_del_loc.face[0]]=face; */
				/* 	} */
				/*||end||*/

			}
			else if (face_neig > 0 && mmr_el_status(Mesh_id, face_neig) == MMC_INACTIVE) {
				/*kbw
					if(abs(face)==669||abs(face)==739){
					  printf("clustering back elem %d, face %d, INACTIVE neighbor %d\n",
						 fath, abs(face), face_neig );
					}
				kew*/
				/* the neighbor is divided - inform its sons on change of neighbors */
				mmr_fa_fam(Mesh_id, face, face_sons, NULL);
				for (i = 1; i <= face_sons[0]; i++) {
					iaux = face_sons[i];
					for (j = 0; j < num_sons; j++) {
						if (fath_sons[j] == mesh->face[iaux].neig[0]) {
							mesh->face[iaux].neig[0] = MMC_BIG_NGB;
							break;
						}
						else if (fath_sons[j] == mesh->face[iaux].neig[1]) {
							mesh->face[iaux].neig[1] = MMC_BIG_NGB;
							break;
						}
					}
				}
			}
			else {
				printf("error in neighbors for element %d in clust_el8_p!\n", El);
				exit(-1);
			}

		}

		/*kbw
			if(El==200){
			  printf("\nafter ifa %d, sons:",face);
			  for(ifason=0;ifason<num_sons;ifason++){
			printf("%d (%d), ", fath_sons[ifason], mesh->elem[fath].sons[ifason]);
			  }
			}
		kew*/

	}

	for (ison = 0; ison < num_sons; ison++) {
		son = fath_sons[ison];

		if (ison < 3) {

			/* delete horizontal internal faces */
			mmr_del_face(Mesh_id, abs(mesh->elem[son].face[1]));

		}
		else if (ison == 3) {

			/* delete horizontal internal face */
			face = abs(mesh->elem[son].face[1]);
			/* delete internal edges */
			for (i = 0; i < 3; i++) mmr_del_edge(Mesh_id,
						abs(mesh->face[face].edge[i]));
			mmr_del_face(Mesh_id, face);

			/* delete vertical internal faces */
			mmr_del_face(Mesh_id, abs(mesh->elem[son].face[2]));
			mmr_del_face(Mesh_id, abs(mesh->elem[son].face[3]));
			mmr_del_face(Mesh_id, abs(mesh->elem[son].face[4]));

		}
		else if (ison == 7) {

			/* delete vertical internal faces */
			mmr_del_face(Mesh_id, abs(mesh->elem[son].face[2]));
			mmr_del_face(Mesh_id, abs(mesh->elem[son].face[3]));
			mmr_del_face(Mesh_id, abs(mesh->elem[son].face[4]));

		}

		/* inactivate sons */
		mesh->elem[son].type = MMC_FREE;
		mesh->elem[son].mate = mesh->parm.pfel;
		mesh->elem[son].refi = MMC_NOT_REF;
		free(mesh->elem[son].face);
		/*||begin||*/
		/* no meaning */
		/*     mesh->elem[son].ipid = -2; */
		/*||end||*/

		/* put son as the first available free space */
		mesh->parm.pfel = son;
		/* decrease the number of active elements */
		mesh->parm.nrel--;


		/*kbw
			if(El==200){
			  printf("\nafter son %d, sons:",son);
			  for(ifason=0;ifason<num_sons;ifason++){
			printf("%d (%d), ", fath_sons[ifason], mesh->elem[fath].sons[ifason]);
			  }
			}
		kew*/

	}

	/* change parameters for the father */
	mesh->parm.nrel++;
	mesh->elem[fath].type *= -1;

	/*kbw
	if(fath>0){
	printf("Clustered back TO element %d, type %d, mate %d, fath %d, refi %d\n",
	fath,mesh->elem[fath].type,mesh->elem[fath].mate,mesh->elem[fath].fath,
	mesh->elem[fath].refi);
	printf("faces: ");
	for(ino=0;ino<5;ino++) printf(" %d ",
	mesh->elem[fath].face[ino]);
	printf("\nnodes: ");
	mmr_el_node_coor(Mesh_id,fath,el_node,NULL);
	for(ino=0;ino<el_node[0];ino++) printf(" %d ",el_node[ino+1]);
	printf("\n");
	}
	kew*/

	return (0);
}

