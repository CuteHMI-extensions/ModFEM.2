/************************************************************************
File mms_prism_input_grad.c - interface between prismatic element mesh
						manipulation module and the "gradmesh" mesh generator
						(contact: jkucwaj@zms.pk.edu.pl)

Contains routines:
  mmr_import_mesh_grad - to read mesh data from input file created by
	"gradmesh" mesh generator (jkucwaj@zms.pk.edu.pl)
		 and generate structured 3D mesh of prismatic elements

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

/*---------------------------------------------------------
mmr_import_mesh_grad - to read mesh data from input file created by
	"gradmesh" and generate structured 3D mesh of prismatic elements
---------------------------------------------------------*/
int mmr_import_mesh_grad( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id,	  /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		char * Filename /* in: name of the file to read mesh data */
)
{
	/* auxiliary variables */
	FILE * fp;
	int i, iaux, istr, num_el_lay, top_bc, bottom_bc;
	int nface, nr_nodel, ino, nodel[4], neig, ifa, ied, i_lay;
	int num_nod_jk, num_el_jk, nr_hor_face, nr_ver_edge, mat_num;
	double z_bottom, z_top, dz, x_temp, y_temp;
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* open the intput file */
	fp = fopen(Filename, "r");
	if (fp == NULL) {
		printf("Not found file '%s' with \"gradmesh\" mesh data!!! Exiting\n", Filename);
		exit(-1);
	}

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	/* read dimensions for mesh arrays */
	fscanf(fp, "%d %d %d %d\n", &mesh->parm.mxno, &mesh->parm.mxed,	&mesh->parm.mxfa, &mesh->parm.mxel );
	printf("TOUT: DIMENSIONS FOR MESH ARRAYS: %d %d %d %d\n", mesh->parm.mxno, mesh->parm.mxed,	mesh->parm.mxfa, mesh->parm.mxel);


	/* read the number of stored node structures, the bottom nad top */
	/* z coordinate,  the number of elements across the thickness */
	/* the flags for top and bottom boundary conditions */
	fscanf(fp, "%d %lg %lg %d %d %d\n", &num_nod_jk, &z_bottom, &z_top, &num_el_lay, &bottom_bc, &top_bc);
	mesh->parm.nmno = num_nod_jk * (num_el_lay + 1);
	mesh->parm.nrno = mesh->parm.nmno;
	mesh->parm.pfno = mesh->parm.nmno + 1;
	dz = (z_top - z_bottom) / num_el_lay;
#ifdef DEBUG_MMM
	if (dz < 0) {
		printf("wrong bottom and top coordinates for layers\n");
		return -1;
	}
#endif
	/* number of vertical edges */
	nr_ver_edge = num_nod_jk * num_el_lay;

	/* adjust the number of allocated structures */
	if (mesh->parm.mxno < mesh->parm.nmno)
		mesh->parm.mxno = mesh->parm.nmno;

	/* allocate space for nodes' structures */
	mesh->node = (mmt_nodes *)malloc((mesh->parm.mxno + 1) * sizeof(mmt_nodes));

	/* adjust the number of allocated edge structures */
	/* - at least 4*(num_el_lay+1) times as much as nodes in one layer */
	if (mesh->parm.mxed < 4 * num_nod_jk * (num_el_lay + 1))
		mesh->parm.mxed = 4 * num_nod_jk * (num_el_lay + 1);

	/* allocate space for edges' structures */
	mesh->edge =
			(mmt_edges *)malloc((mesh->parm.mxed + 1) * sizeof(mmt_edges));

	printf("TOUT: num_nod_jk: %d \n", num_nod_jk);

	/* read nodes structures */
	for (istr = 1; istr <= num_nod_jk; istr++) {
		fscanf(fp, "%lg %lg", &x_temp, &y_temp);

		for (i_lay = 0; i_lay <= num_el_lay; i_lay++) {

			/* set interprocessor ID to zero */
			//mesh->node[istr+i_lay*num_nod_jk].ipid = 0;

			mesh->node[istr + i_lay * num_nod_jk].x = x_temp;
			mesh->node[istr + i_lay * num_nod_jk].y = y_temp;
			mesh->node[istr + i_lay * num_nod_jk].z = z_bottom + i_lay * dz;

			printf("TOUT: x: %lg, y: %lg, z: %lg\n", x_temp, y_temp, z_bottom + i_lay * dz);

			/*kbw
			printf("node %d, x - %lf, y - %lf, z - %lf\n",
			istr+i_lay*num_nod_jk,mesh->node[istr+i_lay*num_nod_jk].x,
			mesh->node[istr+i_lay*num_nod_jk].y,mesh->node[istr+i_lay*num_nod_jk].z);
			kew*/

			/* create vertical edges */
			if (i_lay < num_el_lay) {

				/* set interprocessor ID to zero */
				//mesh->edge[istr+i_lay*num_nod_jk].ipid = 0;

				mesh->edge[istr + i_lay * num_nod_jk].type = MMC_EDGE;
				mesh->edge[istr + i_lay * num_nod_jk].node[0] = istr + i_lay * num_nod_jk;
				mesh->edge[istr + i_lay * num_nod_jk].node[1] = istr + (i_lay + 1) * num_nod_jk;

				printf("TOUT: Edge: 0: %d, 1: %d \n", mesh->edge[istr + i_lay * num_nod_jk].node[0], mesh->edge[istr + i_lay * num_nod_jk].node[1]);
			}
		}
	}

	/* report */
#ifdef DEBUG_MMM
	printf("\nMesh module run in DEBUG mode!\n\n");
	printf("Nodes   : allocated %d structures, read %d structures for %d nodes\n",
			mesh->parm.mxno, mesh->parm.nmno, mesh->parm.nrno);
#endif

	/* read the number of stored element structures */
	fscanf(fp, "%d\n", &num_el_jk);
	mesh->parm.nmel = num_el_jk * num_el_lay;

	/* adjust the number of allocated element structures */
	if (mesh->parm.mxel < mesh->parm.nmel)
		mesh->parm.mxel = mesh->parm.nmel;

	/* allocate space for elements' structures */
	mesh->elem = (mmt_elems *)malloc((mesh->parm.mxel + 1) * sizeof(mmt_elems));

	/* adjust the number of allocated face structures */
	/* - at least three times as much as elements */
	if (mesh->parm.mxfa < 3 * num_el_jk * (num_el_lay + 1))
		mesh->parm.mxfa = 3 * num_el_jk * (num_el_lay + 1);

	/* allocate space for faces' structures */
	mesh->face = (mmt_faces *)malloc((mesh->parm.mxfa + 1) * sizeof(mmt_faces));

	/* set parameters for elements */
	mesh->parm.nrel = mesh->parm.nmel;
	mesh->parm.pfel = mesh->parm.nmel + 1;
	/* number of horizontal faces */
	nr_hor_face = num_el_jk * (num_el_lay + 1);

	/* read element structures and create face structures */
	nface = 0;
	for (istr = 1; istr <= num_el_jk; istr++) {

		/* read material type indicator */
		fscanf(fp, "%d", &iaux);
		mat_num = abs(iaux);

		if (iaux > 0) {
			nr_nodel = 3;
		}
		else {
			nr_nodel = 4;
		}

		/* read nodes to temporary array */
		for (ino = 0; ino < nr_nodel; ino++) {
			fscanf(fp, "%d", &nodel[ino]);
		}

		/* kbw
		printf("JK tria %d, vertices %d %d %d\n",
		   istr, nodel[0], nodel[1], nodel[2]);
		printf("x1 %lf, y1 %lf, x2 %lf, y2 %lf, x3 %lf, y3 %lf\n",
		   mesh->node[nodel[0]].x, mesh->node[nodel[0]].y,
		   mesh->node[nodel[1]].x, mesh->node[nodel[1]].y,
		   mesh->node[nodel[2]].x, mesh->node[nodel[2]].y);
		kew*/

#ifdef DEBUG_MMM
		x_temp = mesh->node[nodel[1]].x - mesh->node[nodel[0]].x;
		y_temp = mesh->node[nodel[1]].y - mesh->node[nodel[0]].y;
		daux = mesh->node[nodel[2]].x - mesh->node[nodel[0]].x;
		faux = mesh->node[nodel[2]].y - mesh->node[nodel[0]].y;
		if (x_temp * faux - y_temp * daux < 0.0) {
			printf("Clockwise orientation in JK file (%lf) !!!\n",
					x_temp * faux - y_temp * daux);
		}
#endif

		for (i_lay = 0; i_lay < num_el_lay; i_lay++) {

			/* set interprocessor ID to zero */
			//mesh->elem[istr+i_lay*num_el_jk].ipid = 0;

			if (nr_nodel == 3) {
				mesh->elem[istr + i_lay * num_el_jk].type = MMC_PRISM;
				i = 5;
				mesh->elem[istr + i_lay * num_el_jk].face =
						mmr_ivector(i, "element faces in read data");
			}
			else {
				mesh->elem[istr + i_lay * num_el_jk].type = MMC_BRICK;
				printf("El_mate (A.K.A. mat_num) ==0 in mesh_jk.dat. Brick not implemeted. Exiting.");
				exit(-1);
			}

			mesh->elem[istr + i_lay * num_el_jk].mate = mat_num;
			mesh->elem[istr + i_lay * num_el_jk].fath = MMC_NO_FATH;
			mesh->elem[istr + i_lay * num_el_jk].refi = MMC_NOT_REF;
			mesh->elem[istr + i_lay * num_el_jk].sons = NULL;

			/* create horizontal faces */
			if (i_lay == 0)
				mesh->face[istr + i_lay * num_el_jk].bc = bottom_bc;
			/* we assume no shift between neighboring elements */
			else mesh->face[istr + i_lay * num_el_jk].bc = 0;
			if (i_lay == num_el_lay - 1)
				mesh->face[istr + (i_lay + 1)*num_el_jk].bc = top_bc;
			/* we assume no shift between neighboring elements */
			else mesh->face[istr + (i_lay + 1)*num_el_jk].bc = 0;

			if (nr_nodel == 3) {
				/* for PRISM elements */

				/* set interprocessor ID to zero */
				//mesh->face[istr+i_lay*num_el_jk].ipid = 0;
				//if(i_lay==num_el_lay-1)
				//mesh->face[istr+(i_lay+1)*num_el_jk].ipid = 0;

				/* create horizontal triangular faces */
				mesh->face[istr + i_lay * num_el_jk].type = MMC_TRIA;
				mesh->face[istr + i_lay * num_el_jk].sons = NULL;
				if (i_lay == num_el_lay - 1) {
					mesh->face[istr + (i_lay + 1)*num_el_jk].type = MMC_TRIA;
					mesh->face[istr + (i_lay + 1)*num_el_jk].sons = NULL;
				}

				if (i_lay == 0)
					mesh->elem[istr + i_lay * num_el_jk].face[0] =
							istr + i_lay * num_el_jk;
				else
					mesh->elem[istr + i_lay * num_el_jk].face[0] =
							-(istr + i_lay * num_el_jk);
				mesh->elem[istr + i_lay * num_el_jk].face[1] =
						istr + (i_lay + 1) * num_el_jk;

				if (i_lay == 0) {
					mesh->face[istr + i_lay * num_el_jk].neig[0] =
							istr + i_lay * num_el_jk;
					mesh->face[istr + i_lay * num_el_jk].neig[1] = MMC_BOUNDARY;
				}
				else {
					mesh->face[istr + i_lay * num_el_jk].neig[0] =
							istr + (i_lay - 1) * num_el_jk;
					mesh->face[istr + i_lay * num_el_jk].neig[1] =
							istr + i_lay * num_el_jk;
				}
				/* for the last layer of faces */
				if (i_lay == num_el_lay - 1) {
					mesh->face[istr + (i_lay + 1)*num_el_jk].neig[0] =
							istr + i_lay * num_el_jk;
					mesh->face[istr + (i_lay + 1)*num_el_jk].neig[1] = MMC_BOUNDARY;
				}
			}
			if (nr_nodel == 4) {
				/* for BRICK elements */
			}

		}

		/* create or find edges for layer zero */
		for (ino = 0; ino < nr_nodel; ino++) {
			fscanf(fp, "%d", &neig);

			/* for internal faces */
			if (neig > 0) {

				/* check whether a face joining two neighbors already exists */
				for (iaux = 0; iaux < nface; iaux++) {
					ifa = nr_hor_face + iaux * num_el_lay + 1;

					/* if face exists */
					if ((mesh->face[ifa].neig[0] == istr &&
									mesh->face[ifa].neig[1] == neig) ||
							(mesh->face[ifa].neig[1] == istr &&
									mesh->face[ifa].neig[0] == neig)) {

						/* first element creates face with proper orientation */
						for (i_lay = 0; i_lay < num_el_lay; i_lay++) {
							mesh->elem[istr + i_lay * num_el_jk].face[2 + ino] =
									-(ifa + i_lay);
						}
						neig = 0;
						break;
					}

				}

				/* if face not found, create a new one */
				if (neig) {

					nface++;
					if (num_el_jk * (num_el_lay + 1) + nface * num_el_lay > mesh->parm.mxfa) {
						printf("Too much faces, increase mxfa!\n");
						exit(-1);
					}
					for (i_lay = 0; i_lay < num_el_lay; i_lay++) {

						ifa = nr_hor_face + (nface - 1) * num_el_lay + i_lay + 1;
						mesh->elem[istr + i_lay * num_el_jk].face[2 + ino] =
								ifa;


						/* set interprocessor ID to zero */
						//mesh->face[ifa].ipid = 0;

						/* type - quadrilateral */
						mesh->face[ifa].type = MMC_QUAD;
						mesh->face[ifa].sons = NULL;
						/* shift between neighboring elements - 1 node */
						mesh->face[ifa].bc = -1;
						mesh->face[ifa].neig[0] = istr + i_lay * num_el_jk;
						mesh->face[ifa].neig[1] = neig + i_lay * num_el_jk;


						/* create horizontal edges */
						ied = nr_ver_edge + (nface - 1) * (num_el_lay + 1) + i_lay + 1;

						/* set interprocessor ID to zero */
						//mesh->edge[ied].ipid = 0;

						mesh->edge[ied].type = MMC_EDGE;
						mesh->edge[ied].node[0] =
								nodel[ino] + i_lay * num_nod_jk;
						mesh->edge[ied].node[1] =
								nodel[(ino + 1) % nr_nodel] + i_lay * num_nod_jk;

						if (i_lay == num_el_lay - 1) {

							/* set interprocessor ID to zero */
							//mesh->edge[ied+1].ipid = 0;

							mesh->edge[ied + 1].type = MMC_EDGE;
							mesh->edge[ied + 1].node[0] =
									nodel[ino] + (i_lay + 1) * num_nod_jk;
							mesh->edge[ied + 1].node[1] =
									nodel[(ino + 1) % nr_nodel] + (i_lay + 1) * num_nod_jk;
						}

						/* assign edges to faces */
						mesh->face[ifa].edge[0] = ied;
						mesh->face[ifa].edge[1] =
								nodel[(ino + 1) % nr_nodel] + i_lay * num_nod_jk;
						mesh->face[ifa].edge[2] = -(ied + 1);
						mesh->face[ifa].edge[3] =
								-(nodel[ino] + i_lay * num_nod_jk);

						/* kbw
								  printf("JK tria %d, neig %d\n",
									 istr, neig);
						kew*/

						/*kbw
								  printf("created face %d , type %d, bc %d, neig1 %d, neig2 %d\nedges:",
									 ifa,
									 mesh->face[ifa].type,
									 mesh->face[ifa].bc,
									 mesh->face[ifa].neig[0],
									 mesh->face[ifa].neig[1] );
								  for(i=0;i<4;i++){
								printf("  %d",mesh->face[ifa].edge[i]);
								  }
								  printf("\n");
						kew*/

					} /* end loop over layers i_lay */

				}

			} /* end if internal face */

			/* for boundary faces */
			else if (neig < 0) {

				/* create new boundary face */
				nface++;
				if (num_el_jk * (num_el_lay + 1) + nface * num_el_lay > mesh->parm.mxfa) {
					printf("Too much faces, increase mxfa!\n");
					exit(-1);
				}

				for (i_lay = 0; i_lay < num_el_lay; i_lay++) {
					ifa = nr_hor_face + (nface - 1) * num_el_lay + i_lay + 1;
					mesh->elem[istr + i_lay * num_el_jk].face[2 + ino] =
							ifa;


					/* set interprocessor ID to zero */
					//mesh->face[ifa].ipid = 0;

					/* type - quadrilateral */
					mesh->face[ifa].type = MMC_QUAD;
					mesh->face[ifa].sons = NULL;
					/* boundary condition flag */
					mesh->face[ifa].bc = -neig;
					mesh->face[ifa].neig[0] = istr + i_lay * num_el_jk;
					mesh->face[ifa].neig[1] = MMC_BOUNDARY;

					/* create horizontal edges */
					ied = nr_ver_edge + (nface - 1) * (num_el_lay + 1) + i_lay + 1;

					/* set interprocessor ID to zero */
					//mesh->edge[ied].ipid = 0;

					mesh->edge[ied].type = MMC_EDGE;
					mesh->edge[ied].node[0] =
							nodel[ino] + i_lay * num_nod_jk;
					mesh->edge[ied].node[1] =
							nodel[(ino + 1) % nr_nodel] + i_lay * num_nod_jk;

					if (i_lay == num_el_lay - 1) {

						/* set interprocessor ID to zero */
						//mesh->edge[ied+1].ipid = 0;

						mesh->edge[ied + 1].type = MMC_EDGE;
						mesh->edge[ied + 1].node[0] =
								nodel[ino] + (i_lay + 1) * num_nod_jk;
						mesh->edge[ied + 1].node[1] =
								nodel[(ino + 1) % nr_nodel] + (i_lay + 1) * num_nod_jk;

					}

					/* assign edges to faces */
					mesh->face[ifa].edge[0] = ied;
					mesh->face[ifa].edge[1] =
							nodel[(ino + 1) % nr_nodel] + i_lay * num_nod_jk;
					mesh->face[ifa].edge[2] = -(ied + 1);
					mesh->face[ifa].edge[3] =
							-(nodel[ino] + i_lay * num_nod_jk);

					/*kbw
					printf("created boundary face %d , type %d, bc %d, neig1 %d, neig2 %d\nedges:",
					ifa,
					mesh->face[ifa].type,
					mesh->face[ifa].bc,
					mesh->face[ifa].neig[0],
					mesh->face[ifa].neig[1] );
					for(i=0;i<4;i++){
					  printf("  %d",mesh->face[ifa].edge[i]);
					}
					printf("\n");
					kew*/

				} /* end of loop over layers i_lay */

			} /* end if boundary face */
		} /* end loop over faces - ino */


		for (ino = 0; ino < nr_nodel; ino++) {
			/* loop again over layers and update edges for horizontal faces */
			for (i_lay = 0; i_lay < num_el_lay; i_lay++) {

				if (nr_nodel == 3) {
					/* for PRISMS */

					ifa = istr + i_lay * num_el_jk;
					iaux = mesh->elem[istr + i_lay * num_el_jk].face[2 + ino];
					if (iaux > 0) ied = mesh->face[iaux].edge[0];
					else ied = -mesh->face[-iaux].edge[0];
					if (i_lay == 0) {
						mesh->face[ifa].edge[2 - ino] = -ied;
					}
					else {
						mesh->face[ifa].edge[ino] = ied;
					}

					if (i_lay == num_el_lay - 1) {
						ifa = istr + (i_lay + 1) * num_el_jk;
						iaux = mesh->elem[istr + i_lay * num_el_jk].face[2 + ino];
						if (iaux > 0) ied = mesh->face[iaux].edge[2];
						else ied = -mesh->face[-iaux].edge[2];
						mesh->face[ifa].edge[ino] = -ied;
					}

					/*checking*/
#ifdef DEBUG_MMM
					if (i_lay > 0) {
						ifa = istr + i_lay * num_el_jk;
						iaux = mesh->elem[istr + (i_lay - 1) * num_el_jk].face[2 + ino];
						if (iaux > 0) {
							ied = mesh->face[iaux].edge[2];
							if (mesh->face[ifa].edge[ino] != -ied) {
								printf("Wrong edge %d for face %d, %d != %d\n",
										ino, ifa, mesh->face[ifa].edge[ino], ied);
							}
						}
					}
#endif

					/*kbw
					if(ino==nr_nodel-1){
					ifa=istr+i_lay*num_el_jk;
					printf("created face %d , type %d, bc %d, neig1 %d, neig2 %d\nedges:",ifa,
					mesh->face[ifa].type,
					mesh->face[ifa].bc,
					mesh->face[ifa].neig[0],
					mesh->face[ifa].neig[1] );
					for(i=0;i<nr_nodel;i++){
					  printf("  %d",mesh->face[ifa].edge[i]);
					}
					printf("\n");
					if(i_lay==num_el_lay-1){
					printf("created face %d , type %d, bc %d, neig1 %d, neig2 %d\nedges:",
					ifa+num_el_jk,
					mesh->face[ifa+num_el_jk].type,
					mesh->face[ifa+num_el_jk].bc,
					mesh->face[ifa+num_el_jk].neig[0],
					mesh->face[ifa+num_el_jk].neig[1] );
					for(i=0;i<nr_nodel;i++){
					  printf("  %d",mesh->face[ifa+num_el_jk].edge[i]);
					}
					printf("\n");
					}
					}
					kew*/

				}
				else if (nr_nodel == 4) {
					/* for BRICKS */

				}

			} /* end of loop over layers: i_lay */
		} /* end of loop over faces: ino */


		/*kbw
		for(i_lay=0;i_lay<num_el_lay;i_lay++){
		int el_node[MMC_MAXELVNO+1];
		printf("created element %d, type %d, fath %d,  mate %d\n",
		istr+i_lay*num_el_jk,
		mesh->elem[istr+i_lay*num_el_jk].type,
		mesh->elem[istr+i_lay*num_el_jk].fath,
		mesh->elem[istr+i_lay*num_el_jk].mate);
		printf("faces: ");
		for(ino=0;ino<5;ino++) printf(" %d ",
		mesh->elem[istr+i_lay*num_el_jk].face[ino]);
		printf("\nnodes: ");
		mmr_el_node_coor(Mesh_id,istr+i_lay*num_el_jk,el_node,NULL);
		for(ino=0;ino<el_node[0];ino++) printf(" %d ",el_node[ino+1]);
		printf("\n");
		}
		getchar();
		kew*/

	} /* end loop over elements */

	/* set parameters for edges */
	mesh->parm.nred = num_nod_jk * num_el_lay + nface * (num_el_lay + 1);
	mesh->parm.nmed = mesh->parm.nred;
	mesh->parm.pfed = mesh->parm.nred + 1;

	/* set parameters for faces */
	mesh->parm.nrfa = num_el_jk * (num_el_lay + 1) + nface * num_el_lay;
	mesh->parm.nmfa = mesh->parm.nrfa;
	mesh->parm.pffa = mesh->parm.nrfa + 1;

	/* report */
#ifdef DEBUG_MMM
	printf("Edges   : allocated %d structures, read %d structures for %d edges\n",
			mesh->parm.mxed, mesh->parm.nmed, mesh->parm.nred);

	printf("Faces   : allocated %d structures, read %d structures for %d faces\n",
			mesh->parm.mxfa, mesh->parm.nmfa, mesh->parm.nrfa);

	printf("Elements: allocated %d structures, read %d structures for %d elements\n",
			mesh->parm.mxel, mesh->parm.nmel, mesh->parm.nrel);
#endif

	/* close file with data for a given mesh*/
	fclose(fp);

	return (1);
}
