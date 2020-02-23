#ifndef LAH_AMG_INTERFACE_H
#define LAH_AMG_INTERFACE_H

#ifdef __cplusplus
extern "C" {
#endif

/*---------------------------------------------------------
  lar_allocate_SM_and_LV - to allocate space for stiffness matrix and load vector
---------------------------------------------------------*/
extern int lar_allocate_SM_and_LV_petsc( // returns: matrix index in itv_matrices array
  int Max_SM_size, /* maximal size of element stiffness matrix */
  int Nrblocks,    /* in: number of DOF blocks */
  int Nrdof_glob,  /* in: total number of DOFs */
  int* Nrdofbl,	   /* in: list of numbers of dofs in a block */
  int* Posglob,	   /* in: list of global numbers of first dof */
  int* Nroffbl,	   /* in: list of numbers of off diagonal blocks */
  // if Nroffbl[iblock]<=0 - iblock is a ghost block with no rows associated
  int** L_offbl	   /* in: list of lists of off diagonal blocks */
	);

extern int lar_initialize_SM_and_LV_petsc(int Matrix_id, int Comp_type);

extern double lar_get_storage_petsc(int Matrix_id);

/*------------------------------------------------------------
  lar_assemble_SM_and_LV - to assemble entries to the global stiffness matrix
                           and the global load vector using the provided local
                           stiffness matrix and load vector
------------------------------------------------------------*/
extern int lar_assemble_SM_and_LV_petsc(
                         /* returns: >=0 - success code, <0 - error code */
  int Matrix_id,   /* in: matrix ID */
  int Comp_type,         /* in: indicator for the scope of computations: */
                         /*   ITC_SOLVE - solve the system */
                         /*   ITC_RESOLVE - resolve for the new rhs vector */
  int Nr_dof_bl,         /* in: number of global dof blocks */
                         /*     associated with the local stiffness matrix */
  int* L_bl_id,          /* in: list of dof blocks' IDs */
  int* L_bl_nrdof,       /* in: list of blocks' numbers of dof */
  double* Stiff_mat,     /* in: stiffness matrix stored columnwise */
  double* Rhs_vect,      /* in: rhs vector */
  char* Rewr_dofs         /* in: flag to rewrite or sum up entries */
                         /*   'T' - true, rewrite entries when assembling */
                         /*   'F' - false, sum up entries when assembling */
  );

extern int lar_allocate_preconditioner_petsc(int Matrix_id);

extern int lar_fill_preconditioner_petsc(int Matrix_id);

extern int lar_free_preconditioner_petsc(int Matrix_id);

extern int lar_free_SM_and_LV_petsc(int Matrix_id);

extern void lar_compute_residual_petsc(int Matrix_id, int Use_rhs, int Ini_zero, int Ndof, double* X, double* B, double* V);

extern void lar_perform_BJ_or_GS_iterations_petsc(
  int Matrix_id,   /* in: matrix ID */
  int Use_rhs,	/* in: 0 - no rhs, 1 - with rhs */
  int Ini_zero,	/* in: flag for zero initial guess */
  int Nr_prec,  /* in: number of preconditioner iterations */
  int Ndof,	/* in: number of unknowns (components of v*) */
  double* V,	/* in,out: vector of unknowns updated */
		/* during the loop over subdomains */
  double* B	/* in:  the rhs vector, if NULL take rhs */
		/*      from block data structure */
	);

extern void lar_perform_rhsub_petsc(int Matrix_id, int Ndof,	double* V, double* B);

extern int lar_block_print_matrix_petsc(int Matrix_id);

#ifdef __cplusplus
}
#endif

#endif
