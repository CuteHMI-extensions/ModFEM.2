#ifndef __uth_err_h__
#define __uth_err_h__

/// \file Contains informations about COMPUTATIONAL errors during solution
/// and from solved problems domain.
/// It is not suitable for handling programming,system or other kind of errors.
/// K.Michalik - 04.2016 initial version


#ifdef __cplusplus
extern "C"
{
#endif

typedef enum {
    NO_ERROR = 0,
    ERR_OUTFLOW_BACKWATER,
    ERR_FINAL_VELOCITY_MAG_EXCEEDED
} utt_err_code;


#ifdef __cplusplus
}
#endif

#endif //! __uth_err_h__
