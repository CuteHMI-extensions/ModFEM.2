#include <modfem/uth_log.h>


#ifdef PARALLEL
#include <stdio.h>
#include "pch_intf.h"
#endif


#ifdef __cplusplus
extern "C"
{
#endif

/// For logging support.
FILE * utv_log_out = stderr;

#ifdef PARALLEL
const char * utr_log_pre()
{
	static char * pre_txt = NULL;

	if (pre_txt == NULL) {
		if ( pcr_my_proc_id() != -1 ) {
			pre_txt = new char[4];
			sprintf(pre_txt, "%d", pcr_my_proc_id());
		}
	}

	return pre_txt == NULL ? "" : pre_txt;
}
#else
const char * utr_log_pre()
{
	return "";
}
#endif

#ifdef __cplusplus
}
#endif
