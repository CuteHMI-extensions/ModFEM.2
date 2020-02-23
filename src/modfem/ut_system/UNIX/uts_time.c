/******************************************************************************
Plik pomiar_czasu.c - procedury pomiaru czasu dla roznych systemow operacyjnych

Procedury wykorzystujace funkcje standardowe C:
  time_C - zwraca czas CPU w sekundach od rozpoczecia dzialania lub -1 jesli
		   funkcja clock() nie jest prawidlowo zaimplementowana w systemie

Procedury UNIX-owe:
  time_init   - inicjacja pomiaru czasu
  time_print    - pomiar i wydruk czasu CPU i zegarowego od momentu inicjacji

// ONLY time_CPU and time_clock are thread safe !!!
  time_CPU       - zwrocenie czasu CPU od momentu inicjacji
  time_clock    - zwrocenie czasu zegarowego od momentu inicjacji

Sposob uzycia:
  zgodnie z interfejsem w pliku pomiar_czasu.h

------------------------------
Historia:
	  02.2002 - Krzysztof Banas, kbanas@pk.edu.pl
******************************************************************************/

#include<stdlib.h>
#include<stdio.h>
#include<time.h>
#include<sys/time.h>
#include<sys/resource.h>

/* interfejs */
#include <modfem/uth_system.h>

#ifdef _OPENMP
#include<omp.h>
#endif

/* zmienne statyczne do wielokrotnego wykorzystania procedur pomiaru */
static clock_t ct;
static struct rusage rp;
static struct timeval tp;
static struct timezone tzp;

/*---------------------------------------------------------
  czas_C - zwraca czas CPU w sekundach od rozpoczecia dzialania lub -1 jesli
		   funkcja clock() nie jest prawidlowo zaimplementowana w systemie
---------------------------------------------------------*/
double time_C()
{
  clock_t time;

  time = clock();

  if(time==-1) return( (double) time);
  else return ( (double) time / (double) CLOCKS_PER_SEC );
}

/*---------------------------------------------------------
inicjuj_czas - inicjuje pomiar czasu nadajac wartosci
			   poczatkowe strukturom tp, tzp i rp
---------------------------------------------------------*/
void time_init()
{
  ct = clock();
  getrusage(RUSAGE_SELF, &rp);
  gettimeofday(&tp,&tzp);
}


/*---------------------------------------------------------
  czas_zegara - procedura zwraca czas zegarowy w sekundach od momentu
				inicjacji (jako liczbe podwojnej precyzji)
---------------------------------------------------------*/
double time_clock()
{

#ifdef _OPENMP

  return(omp_get_wtime());

#else

  struct timeval tk;
  struct timezone tzp;
  double daytime;

  gettimeofday(&tk, &tzp);

  daytime=(tk.tv_usec)/1e6+tk.tv_sec;

  return(daytime);

#endif

}

/*---------------------------------------------------------
  czas_CPU - procedura zwraca czas CPU w sekundach od momentu
			 inicjacji (jako liczbe podwojnej precyzji)
 ---------------------------------------------------------*/
double time_CPU()
{

  struct rusage rk;
  double cputime;

  getrusage(RUSAGE_SELF, &rk);

  cputime = (rk.ru_utime.tv_usec)/1e6;
  cputime += rk.ru_utime.tv_sec;

  return(cputime);
}

/*---------------------------------------------------------
 drukuj_czas - pomiar i wydruk czasu CPU i zegarowego
			   w sekundach od momentu inicjacji pomiaru czasu
 ---------------------------------------------------------*/
void time_print()
{

  clock_t time;
  struct rusage rk;
  struct timeval tk;
  double stdtime, cputime,daytime;

  time = clock();
  getrusage(RUSAGE_SELF, &rk);
  gettimeofday(&tk, &tzp);

  stdtime = (double) (time - ct) / (double)  CLOCKS_PER_SEC ;

  cputime = ( rk.ru_utime.tv_usec - rp.ru_utime.tv_usec ) / 1.0e6 ;
  cputime += rk.ru_utime.tv_sec - rp.ru_utime.tv_sec;

  daytime = ( tk.tv_usec - tp.tv_usec ) / 1.0e6 + tk.tv_sec - tp.tv_sec ;

  printf("time C library   = %lf\n",stdtime);
  printf("time CPU         = %lf\n",cputime);
  printf("time wallclock   = %lf\n",daytime);

}
