/************************************************************************
File mms_prism_io_dump.c - input/output routines for prismatic elements' mesh

Contains routines:
  mmr_read_mesh - to dump-in mesh data stored by previous runs
  mmr_write_mesh - to dump-out mesh data in the standard HP_FEM format

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
  mmr_read_mesh - to dump-in mesh data stored by previous runs
---------------------------------------------------------*/
int mmr_read_mesh( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id,	  /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		char * Filename /* in: name of the file to read mesh data */
)
{
	/* auxiliary variables */
	FILE * fp;
	int i, iaux, istr, ison;
	mmt_mesh * mesh;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* open the input file */
	fp = fopen(Filename, "r");
	if (fp == NULL) {
		printf("Not found file '%s' with mesh data!!! Exiting.\n", Filename);
		exit(-1);
	}

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);

	/* read dimensions for mesh arrays */
	fscanf(fp, "%d %d %d %d\n",
			&mesh->parm.mxno, &mesh->parm.mxed,
			&mesh->parm.mxfa, &mesh->parm.mxel );

	/* read the number of stored node structures and the pointer to */
	/* the first free space */
	fscanf(fp, "%d %d\n", &mesh->parm.nmno, &mesh->parm.pfno);

	/* adjust the number of allocated structures */
	if (mesh->parm.mxno < mesh->parm.nmno)
		mesh->parm.mxno = mesh->parm.nmno;

	/* allocate space for nodes' structures */
	mesh->node =
			(mmt_nodes *)malloc((mesh->parm.mxno + 1) * sizeof(mmt_nodes));

	/* read nodes structures */
	iaux = 0;
	for (istr = 1; istr <= mesh->parm.nmno; istr++) {

		/* set interprocessor ID to zero */
		//mesh->node[istr].ipid = 0;

		fscanf(fp, "%lg", &mesh->node[istr].x);

		if (mesh->node[istr].x < -1e10) {
			fscanf(fp, "%lg\n", &mesh->node[istr].y);

		}
		else {
			fscanf(fp, "%lg %lg\n",
					&mesh->node[istr].y, &mesh->node[istr].z);
			iaux++;
		}

	}

	/* set the number of active nodes */
	mesh->parm.nrno = iaux;

	/* report */
#ifdef DEBUG_MMM
	printf("Mesh %d:\n", Mesh_id);
	printf("Nodes   : allocated %d structures, read %d structures for %d nodes\n",
			mesh->parm.mxno, mesh->parm.nmno, mesh->parm.nrno);
#endif
	/* read the number of stored edge structures and the pointer to */
	/* first free space */
	fscanf(fp, "%d %d\n", &mesh->parm.nmed, &mesh->parm.pfed);

	/* adjust the number of allocated structures */
	if (mesh->parm.mxed < mesh->parm.nmed)
		mesh->parm.mxed = mesh->parm.nmed;

	/* allocate space for nodes' structures */
	mesh->edge =
			(mmt_edges *)malloc((mesh->parm.mxed + 1) * sizeof(mmt_edges));

	/* read edges structures */
	iaux = 0;
	for (istr = 1; istr <= mesh->parm.nmed; istr++) {

		/* set interprocessor ID to zero */
		//mesh->edge[istr].ipid = 0;

		fscanf(fp, "%d %d %d\n",
				&mesh->edge[istr].type,
				&mesh->edge[istr].node[0], &mesh->edge[istr].node[1]);
		if (mesh->edge[istr].type > 0) iaux++;
	}

	/* set the number of active edges */
	mesh->parm.nred = iaux;

	/* report */
#ifdef DEBUG_MMM
	printf("Edges   : allocated %d structures, read %d structures for %d edges\n",
			mesh->parm.mxed, mesh->parm.nmed, mesh->parm.nred);
#endif
	/* read the number of stored face structures and the pointer to */
	/* first free space */
	fscanf(fp, "%d %d\n", &mesh->parm.nmfa, &mesh->parm.pffa);

	/* adjust the number of allocated structures */
	if (mesh->parm.mxfa < mesh->parm.nmfa)
		mesh->parm.mxfa = mesh->parm.nmfa;

	/* allocate space for faces' structures */
	mesh->face =
			(mmt_faces *)malloc((mesh->parm.mxfa + 1) * sizeof(mmt_faces));

	/* read face structures */
	iaux = 0;
	for (istr = 1; istr <= mesh->parm.nmfa; istr++) {

		/* set interprocessor ID to zero */
		//mesh->face[istr].ipid = 0;

		fscanf(fp, "%d %d %d %d\n",
				&mesh->face[istr].type, &mesh->face[istr].bc,
				&mesh->face[istr].neig[0], &mesh->face[istr].neig[1] );
		if (abs(mesh->face[istr].type) == MMC_TRIA) {
			fscanf(fp, "%d %d %d\n",
					&mesh->face[istr].edge[0], &mesh->face[istr].edge[1],
					&mesh->face[istr].edge[2]);
		}
		else if (abs(mesh->face[istr].type) == MMC_QUAD) {
			fscanf(fp, "%d %d %d %d\n",
					&mesh->face[istr].edge[0], &mesh->face[istr].edge[1],
					&mesh->face[istr].edge[2], &mesh->face[istr].edge[3]);
		}

		/* initialize pointer to sons */
		mesh->face[istr].sons = NULL;

		if (mesh->face[istr].type > 0) {
			iaux++;
		}
		/* read sons for inactive faces */
		else if (mesh->face[istr].type < 0) {
			ison = 4;
			mesh->face[istr].sons =
					mmr_ivector(ison, "face sons in read data");
			for (i = 0; i < ison; i++) {
				fscanf(fp, "%d", &mesh->face[istr].sons[i]);
			}
			fscanf(fp, "\n");
		}

	}

	/* set the number of active faces */
	mesh->parm.nrfa = iaux;

	/* report */
#ifdef DEBUG_MMM
	printf("Faces   : allocated %d structures, read %d structures for %d faces\n",
			mesh->parm.mxfa, mesh->parm.nmfa, mesh->parm.nrfa);
#endif
	/* read the number of stored element structures and the pointer to */
	/* first free space */
	fscanf(fp, "%d %d\n", &mesh->parm.nmel, &mesh->parm.pfel);

	/* adjust the number of allocated structures */
	if (mesh->parm.mxel < mesh->parm.nmel)
		mesh->parm.mxel = mesh->parm.nmel;

	/* allocate space for elements' structures */
	mesh->elem = (mmt_elems *)malloc((mesh->parm.mxel + 1) * sizeof(mmt_elems));

	/* read element structures */
	iaux = 0;
	for (istr = 1; istr <= mesh->parm.nmel; istr++) {

		/* set interprocessor ID to zero */
		//mesh->elem[istr].ipid = 0;

		fscanf(fp, "%d %d %d %d\n", &mesh->elem[istr].type,
				&mesh->elem[istr].mate,
				&mesh->elem[istr].fath, &mesh->elem[istr].refi);

		/* initialize pointer to sons */
		mesh->elem[istr].sons = NULL;

		/* for active and inactive elements */
		if (abs(mesh->elem[istr].type) == MMC_BRICK) {
		}
		else if (abs(mesh->elem[istr].type) == MMC_PRISM) {

			i = 5;
			mesh->elem[istr].face = mmr_ivector(i, "element faces in read data");
			fscanf(fp, "%d %d %d %d %d\n", &mesh->elem[istr].face[0],
					&mesh->elem[istr].face[1], &mesh->elem[istr].face[2],
					&mesh->elem[istr].face[3], &mesh->elem[istr].face[4] );

		}
		else if (abs(mesh->elem[istr].type) == MMC_TETRA) {
		}

		/* if element not a free space */
		if (mesh->elem[istr].type != MMC_FREE) {
			/* for inactive elements */
			if (mesh->elem[istr].type < 0) {
				/* refinement kind */
				if (mesh->elem[istr].refi == MMC_REF_ISO) {
					ison = 8;
				}
				else {
					printf("Unknown refinement type for inactive element %d\n",
							istr);
					exit(-1);
				}
				mesh->elem[istr].sons =
						mmr_ivector(ison, "sons in read data");
				for (i = 0; i < ison; i++) {
					fscanf(fp, "%d", &mesh->elem[istr].sons[i]);
				}
				fscanf(fp, "\n");
			}
			else {
				/* for active elements */
				iaux++;
#ifdef DEBUG_MMM
				if (mesh->elem[istr].refi != MMC_NOT_REF) {
					printf("Refinement type indicated for active element %d\n",
							istr);
					exit(-1);
				}
#endif
			}
		}

	}

	/* set the number of active elements */
	mesh->parm.nrel = iaux;

	/* report */
#ifdef DEBUG_MMM
	printf("Elements: allocated %d structures, read %d structures for %d elements\n",
			mesh->parm.mxel, mesh->parm.nmel, mesh->parm.nrel);
#endif

	/* close file with data for a given mesh*/
	fclose(fp);

	return (1);
}

/*---------------------------------------------------------
  mmr_write_mesh - to dump-out mesh data in the standard HP_FEM format
---------------------------------------------------------*/
int mmr_write_mesh( /* returns: >=0 - success code, <0 - error code */
		int Mesh_id,    /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
		char * Filename /* in: name of the file to write mesh data */
)
{

	/* auxiliary variables */
	mmt_mesh * mesh;
	FILE * fp;
	int i, istr, ison;

	/*++++++++++++++++ executable statements ++++++++++++++++*/

	/* open the output file */
	fp = fopen(Filename, "w");
	if (fp == NULL) {
		printf("Cannot open file '%s' for mesh data\n", Filename);
		return (-1);
	}

	/* select the proper mesh data structure */
	mesh = mmr_select_mesh(Mesh_id);


	/* write dimensions for mesh arrays */
	fprintf(fp, "%d %d %d %d\n",
			mesh->parm.mxno, mesh->parm.mxed,
			mesh->parm.mxfa, mesh->parm.mxel );

	fprintf(fp, "%d %d\n", mesh->parm.nmno, mesh->parm.pfno);
	for (istr = 1; istr <= mesh->parm.nmno; istr++) {
		if (mesh->node[istr].x < -1e10) {
			fprintf(fp, "%.12lg %.12lg\n",
					mesh->node[istr].x, mesh->node[istr].y);
		}
		else {
			fprintf(fp, "%.12lg %.12lg %.12lg\n",
					mesh->node[istr].x, mesh->node[istr].y,
					mesh->node[istr].z);
		}
	}

	fprintf(fp, "%d %d\n", mesh->parm.nmed, mesh->parm.pfed);
	for (istr = 1; istr <= mesh->parm.nmed; istr++) {
		fprintf(fp, "%d %d %d\n",
				mesh->edge[istr].type,
				mesh->edge[istr].node[0], mesh->edge[istr].node[1]);
	}

	fprintf(fp, "%d %d\n", mesh->parm.nmfa, mesh->parm.pffa);
	for (istr = 1; istr <= mesh->parm.nmfa; istr++) {
		fprintf(fp, "%d %d %d %d\n",
				mesh->face[istr].type, mesh->face[istr].bc,
				mesh->face[istr].neig[0], mesh->face[istr].neig[1] );
		if (abs(mesh->face[istr].type) == MMC_TRIA) {
			fprintf(fp, "%d %d %d\n",
					mesh->face[istr].edge[0], mesh->face[istr].edge[1],
					mesh->face[istr].edge[2]);
		}
		else if (abs(mesh->face[istr].type) == MMC_QUAD) {
			fprintf(fp, "%d %d %d %d\n",
					mesh->face[istr].edge[0], mesh->face[istr].edge[1],
					mesh->face[istr].edge[2], mesh->face[istr].edge[3]);
		}
		if (mesh->face[istr].type < 0) {
			fprintf(fp, "%d %d %d %d\n",
					mesh->face[istr].sons[0], mesh->face[istr].sons[1],
					mesh->face[istr].sons[2], mesh->face[istr].sons[3]);
		}

	}

	fprintf(fp, "%d %d\n", mesh->parm.nmel, mesh->parm.pfel);
	for (istr = 1; istr <= mesh->parm.nmel; istr++) {

		fprintf(fp, "%d %d %d %d\n", mesh->elem[istr].type,
				mesh->elem[istr].mate, mesh->elem[istr].fath, mesh->elem[istr].refi );

		/* for active and inactive elements */
		if (abs(mesh->elem[istr].type) == MMC_BRICK) {
		}
		else if (abs(mesh->elem[istr].type) == MMC_PRISM) {

			fprintf(fp, "%d %d %d %d %d\n", mesh->elem[istr].face[0],
					mesh->elem[istr].face[1], mesh->elem[istr].face[2],
					mesh->elem[istr].face[3], mesh->elem[istr].face[4] );

		}
		else if (abs(mesh->elem[istr].type) == MMC_TETRA) {
		}

		/* for inactive elements */
		if (mesh->elem[istr].type < 0) {
			/* refinement kind */
			if (mesh->elem[istr].refi == MMC_REF_ISO) {
				ison = 8;
			}
			else {
				printf("Unknown refinement type for inactive element %d\n",
						istr);
				exit(-1);
			}
			for (i = 0; i < ison; i++) {
				fprintf(fp, "%d ", mesh->elem[istr].sons[i]);
			}
			fprintf(fp, "\n");
		}
	}


	/* close file with mesh data */
	fclose(fp);

	return 1;
}

