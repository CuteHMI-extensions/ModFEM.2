/***********************************************************************
 * File: mod_fem_viewer.h
 * 2011 Paweł Macioł, Cracow university of Technology
 *
 * This header file is required to connect to the mod_fem_viwer
 * visualization module
 ***********************************************************************/

#ifndef _MOD_FEM_VIEWER_H_
#define _MOD_FEM_VIEWER_H_

/** @defgroup VIEWER Viewer
 *
 *  @{
 */

#ifdef MOD_FEM_VIEWER
# ifdef __cplusplus
#   define femviewer_export extern "C"
# else
#   define femviewer_export extern
# endif
#else
# define femviewer_export
#endif

#ifdef MOD_FEM_VIEWER
# define _USE_GRAPHICS_MODULE	1
#endif

#include<stdio.h>

#ifdef _USE_GRAPHICS_MODULE
femviewer_export
int init_mod_fem_viewer(int argc, char **argv,FILE* out);
#else
//inline
int init_mod_fem_viewer(int argc,char** argv,FILE* out) {
	return fprintf(out,"Graphic module is not linked!\nPlease compile with \"-DMOD_FEM_VIEWER\" option and run again.\n");
}
#endif

/** @} */ // end of group


#endif /** _MOD_FEM_VIEWER_H_ */
