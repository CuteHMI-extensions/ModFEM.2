/************************************************************************
File pdh_control_intf.h - interface for control parameters and some general 
                          purpose functions from problem dependent module 

Contains declarations of constants and interface routines:

  pdr_module_introduce - to return the problem module's name
  pdr_get_problem_structure - to return pointer to problem structure
  pdr_problem_name - to return the problem's name
  pdr_ctrl_i_params - to return one of control parameters
  pdr_ctrl_d_params - to return one of control parameters
  pdr_adapt_i_params - to return parameters of adaptation
  pdr_adapt_d_params - to return parameters of adaptation
  pdr_time_i_params - to return parameters of time integration
  pdr_time_d_params - to return parameters of time integration
  pdr_set_time_i_params - to set parameters of time integration
  pdr_set_time_d_params - to set parameters of time integration
  pdr_lins_i_params - to return parameters of linear system solution
  pdr_lins_d_params - to return parameters of linear system solution
  pdr_change_data - to change some of control data 

------------------------------  			
History:
	02.2002 - Krzysztof Banas, initial version		
*************************************************************************/

#ifndef _pdh_control_intf_
#define _pdh_control_intf_

#include <stdio.h>

#ifdef __cplusplus
extern "C" 
{
#endif

/** @defgroup PRD_CTRL Problem Control
 *  @ingroup PRD
 *  @{
 */

/**-----------------------------------------------------------
  pdr_get_problem_structure - to return pointer to problem structure
------------------------------------------------------------*/
extern void* pdr_get_problem_structure(int Problem_id);

/*------------------------------------------------------------
  pdr_problem_name - to return the problem's name
------------------------------------------------------------*/
int pdr_problem_name(
                  /* returns: >=0 - success code, <0 - error code */
  int Problem_id,  /* in: problem ID or 
		         PDC_USE_CURRENT_PROBLEM_ID for the current problem */
  char* Problem_name /* out: the name of the problem solved */
  );

/**--------------------------------------------------------
  pdr_ctrl_i_params - to return one of control parameters
  //if(Num==1) return(problem->ctrl.name); - changed for string
  else if(Num==2) return(problem->ctrl.mesh_id);
  else if(Num==3) return(problem->ctrl.field_id);
  else if(Num==4) return(problem->ctrl.init_bvp);
  else if(Num==5) return(problem->ctrl.nreq);
  else if(Num==9) return(problem->ctrl.nr_sol);
  else if(Num==6) return(problem->ctrl.nr_mat);
  else if(Num==7) return(problem->ctrl.slope);
  else if(Num==8) return(problem->ctrl.base);
---------------------------------------------------------*/

extern int pdr_ctrl_i_params( /** returns: integer problem parameter */
  int Problem_id,	/** in: problem ID  */
  int Num         /** in: parameter number in control structure */
  );

/**--------------------------------------------------------
  pdr_ctrl_d_params - to return one of control parameters
  if(Num==8) return(problem->ctrl.coeff1);
  else if(Num==9) return(problem->ctrl.coeff2);
  else if(Num==10) return(problem->ctrl.coeff3);
---------------------------------------------------------*/
extern double pdr_ctrl_d_params( /** returns: integer problem parameter */
	int Problem_id,	/** in: problem ID  */
	int Num         /** in: parameter number in control structure */
	);

/**--------------------------------------------------------
  pdr_adapt_i_params - to return parameters of adaptation
  if(Num==1) return(problem->adpt.type);
  else if(Num==2) return(problem->adpt.interval);
  else if(Num==3) return(problem->adpt.maxgen);
  else if(Num==4) return(problem->adpt.maxgendiff);
  else if(Num==7) return(problem->adpt.monitor);
---------------------------------------------------------*/
extern int pdr_adapt_i_params( /** returns: integer adaptation parameter */
	int Problem_id,	/** in: data structure to be used  */
	int Num         /** in: parameter number in control structure */
	);

/**--------------------------------------------------------
  pdr_adapt_d_params - to return parameters of adaptation
  if(Num==5) return(problem->adpt.eps);
  else if(Num==6) return(problem->adpt.ratio);
---------------------------------------------------------*/
extern double pdr_adapt_d_params( /** returns: real adaptation parameter */
	int Problem_id,	/** in: data structure to be used  */
	int Num         /** in: parameter number in control structure */
	);

/**-----------------------------------------------------------
pdr_time_i_params - to return parameters of timeation
------------------------------------------------------------*/
extern int pdr_time_i_params(int Problem_id, int Num);

/**-----------------------------------------------------------
pdr_time_d_params - to return parameters of timeation
------------------------------------------------------------*/
extern double pdr_time_d_params(int Problem_id, int Num);

/**--------------------------------------------------------
pdr_set_time_i_params - to change parameters of time discretization
---------------------------------------------------------*/
extern void pdr_set_time_i_params( 
        int Problem_id,	     /** in: data structure to be used  */
	int Num,             /** in: parameter number in control structure */
	int Value            /** in: parameter value */
				   );

/**--------------------------------------------------------
pdr_set_time_d_params - to change parameters of time discretization
---------------------------------------------------------*/
extern void pdr_set_time_d_params( 
        int Problem_id,	     /** in: data structure to be used  */
	int Num,             /** in: parameter number in control structure */
	double Value         /** in: parameter value */
			    );

/**--------------------------------------------------------
pdr_lins_i_params - to return parameters of linear equations solver
---------------------------------------------------------*/
extern int pdr_lins_i_params( /** returns: integer linear solver parameter */
	int Problem_id,	/** in: data structure to be used  */
	int Num         /** in: parameter number in control structure */
		       );

/**--------------------------------------------------------
pdr_lins_d_params - to return parameters of linear equations solver
---------------------------------------------------------*/
extern double pdr_lins_d_params( /** returns: real linear solver parameter */
	int Problem_id,	/** in: data structure to be used  */
	int Num         /** in: parameter number in control structure */
			  );

/**--------------------------------------------------------
pdr_change_data - to change some of control data 
---------------------------------------------------------*/
extern void pdr_change_data(
	int Problem_id	/** in: data structure to be used  */
		     );

/** @} */ // end of group


#ifdef __cplusplus
}
#endif

#endif
