#ifndef _PROFTIMER_HPP_
#define _PROFTIMER_HPP_

/** \addtogroup MMM_HYBRID Hybrid Mesh Module
 *  \ingroup MM
 *  @{
 */

#ifdef WIN32
#define WIN
#endif

#ifdef WIN
	#include <windows.h>
#else // LINUX
	#include <time.h>
#endif
// very simple, high-precision (at least in theory) timer for timing API calls


struct ProfTimer {
	void Start(void)
	{
#ifdef WIN
		QueryPerformanceCounter(&mTimeStart);
#else
		clock_gettime(CLOCK_REALTIME,&mTimeStart);
#endif
	}

	void Stop(void)
	{
#ifdef WIN
		QueryPerformanceCounter(&mTimeStop);
#else
		clock_gettime(CLOCK_REALTIME, &mTimeStop);
#endif
	}

	double GetDurationInSecs(void) const
	{
#ifdef WIN
		LARGE_INTEGER freq;
		QueryPerformanceFrequency(&freq);
		return ((double)(mTimeStop.QuadPart-mTimeStart.QuadPart)/(double)freq.QuadPart);
#else
		return double(mTimeStop.tv_sec - mTimeStart.tv_sec) + (double(mTimeStop.tv_nsec - mTimeStart.tv_nsec)/double(1000000000));
#endif
	}

#ifdef WIN
	LARGE_INTEGER mTimeStart;
	LARGE_INTEGER mTimeStop;
#else
	timespec mTimeStart;
	timespec mTimeStop;
#endif
};


/** @}*/

#endif // _PROFTIMER_HPP_
