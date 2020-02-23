/******************************************************************************
Plik uts_time.cpp - time measurments for Windows



------------------------------
Historia:
	  02.2015 - Kazimierz Michalik, kamich@agh.edu.pl
******************************************************************************/
/* interfejs */
#include <modfem/uth_system.h>


// system
#define WINDOWS_LEAN_AND_MEAN
#include<windows.h>
#include<stdio.h>
#include<time.h>

#ifdef _OPENMP
#include<omp.h>
#endif

#ifdef __cplusplus
extern "C"
{
#endif

// pomocnicze zmienne
static LARGE_INTEGER iFrequency;
static LARGE_INTEGER iLastQuery;
static clock_t ct;

/*---------------------------------------------------------
  czas_C - zwraca czas CPU w sekundach od rozpoczecia dzialania lub -1 jesli
		   funkcja clock() nie jest prawidlowo zaimplementowana w systemie
---------------------------------------------------------*/
double time_C()
{
	clock_t timer;

	timer = clock();

	if(timer==-1) return( (double) timer);
	else return ( (double) timer / (double) CLOCKS_PER_SEC );
}
double czas_C() { return time_C(); }
/*---------------------------------------------------------
inicjuj_czas - inicjuje pomiar czasu nadajac wartosci
			   poczatkowe odpowiednim zmiennym
---------------------------------------------------------*/
void time_init()
{
	QueryPerformanceCounter  ( &iLastQuery);
	QueryPerformanceFrequency( &iFrequency);

	ct = clock();
}
void inicjuj_czas() { time_init(); }

/*---------------------------------------------------------
  czas_zegara - procedura zwraca czas zegarowy w sekundach od momentu
				inicjacji (jako liczbe podwojnej precyzji)
---------------------------------------------------------*/
double time_clock()
{

#ifdef _OPENMP

 return(omp_get_wtime());

#else

	time_t daytime;

	time( &daytime);
	return( (double)daytime); // resolution only in whole seconds

#endif
}

double czas_zegara() { return time_clock(); }
/*---------------------------------------------------------
  czas_CPU - procedura zwraca czas CPU w sekundach od momentu
			 inicjacji (jako liczbe podwojnej precyzji)
 ---------------------------------------------------------*/
double time_CPU()
{
	LARGE_INTEGER iCurrentQuery;
	double		  cputime;

	QueryPerformanceCounter( &iCurrentQuery);
	cputime = (double)(iCurrentQuery.QuadPart - iLastQuery.QuadPart) / (double) iFrequency.QuadPart;
	cputime /= 1e3;

	return( cputime);
}
double czas_CPU() { return time_CPU(); }
/*---------------------------------------------------------
 drukuj_czas - pomiar i wydruk czasu CPU i zegarowego
			   od momentu inicjacji pomiaru czasu
 ---------------------------------------------------------*/
void time_print()
{
	LARGE_INTEGER iCurrentQuery;
	clock_t	timer;
	time_t daytime;
	double stdtime, cputime;

	timer = clock();
	QueryPerformanceCounter( &iCurrentQuery);
	time( &daytime);

	stdtime = (double)(timer - ct) / (double) CLOCKS_PER_SEC;
	cputime = (double)(iCurrentQuery.QuadPart - iLastQuery.QuadPart) / (double)(iFrequency.QuadPart);
	cputime /= 1e3;

	printf("czas standardowy = %lf [s]\n",stdtime);
	printf("czas CPU         = %lf [s]\n",cputime);
	printf("czas zegarowy    = %lf [s]\n",daytime);

}
void drukuj_czas() { time_print(); }

#ifdef __cplusplus
}
#endif
