/*
James Clark
	hseries.c
	
	Harmonic series sum threshold
	
	Helps approximate Sum(1/x, 1, N) > M, for an arbitrary M
	
	Relies on very large number formats provided by the C library. Implementation
	may vary from platform to platform.
	
	Programmed in C (C99)
	Currently compiles and works with MinGW gcc for Windows (64-bit).
	
	To-do: 
		*Testing for UNIX
		*Custom number storage for extremely large values
			Built-in types are good up to about M=22 or slightly more, results in
			a value of N greater than 2 billion
		*Possibly multithreading harmonicSeriesEx() to improve performance
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

/* Verbose mode. Remove if you don't like to see lots of numbers */
#define VERBOSE

/* Arbitarily small margin of error */
#define DELTA 0.000000001L

/* ratio of e^M/N ~ 0.5604 to 0.5741 through some analysis.
   As M gets larger, N/e^M approaches a value approximately 5.61454.
   
   The ratios used for the initial guess are slightly higher so that
   our guess approaches from above.
*/
#define EM_N_RATIO 0.564L
#define EM_N_RATIO_MEDIUM 0.5618L
#define EM_N_RATIO_LARGE 0.56147L
#define EM_N_RATIO_VERY_LARGE 0.56146L
#define EM_N_RATIO_HUGE 0.5614596L
#define EM_N_RATIO_COLOSSAL 0.5614595L


/* Initial harmonic series sum, for reporting later */
long double G_sum = 0.0L;


static inline long double harmonicSeriesEx (long nStart, long nEnd)
/* Returns the sum of the harmonic series beginning from an arbitrary begin and end point */
{
	long double sum = 0.0L;
	
	if ((nStart == 0 || nEnd == 0) || (nStart > nEnd))
		return 0.0L;
	
	#ifdef VERBOSE
		printf("Processing harmonic series...\n");
	#endif
	
	for (long i = nStart; i <= nEnd; i++)
		sum += 1.0L/i;
		
	return sum;
}


static long double harmonicSeries (long n)
/* Returns the sum of the harmonic series from 1 to an arbitrary end point */
{
	return harmonicSeriesEx(1, n);
}


static unsigned long hseriesThreshold (long double M)
/* Solves Sum(1/x, 1, N) > M for an arbitrary M */
{
	unsigned long n = 0L;
	long double sum = 0.0L;
	
	#ifdef VERBOSE
		long double diff = 0.0L;
		long double r;
		long num_tests = 0L;
	#endif
	
	/*
  	   Determine an initial guess, will be high. Cast from long double to 
	   unsigned long throws away the decimal portion.
	   N gets large very quickly, so adjust the ratio used to improve our guess
	*/
	if (M >= 20.0L)
		n = (unsigned long)(expl(M) * EM_N_RATIO_COLOSSAL);
	else if (M >= 18.0L)
		n = (unsigned long)(expl(M) * EM_N_RATIO_HUGE);
	else if (M >= 16.0L)
		n = (unsigned long)(expl(M) * EM_N_RATIO_VERY_LARGE);
	else if (M >= 12.0L)
		n = (unsigned long)(expl(M) * EM_N_RATIO_LARGE);
	else if (M >= 9.0L)
		n = (unsigned long)(expl(M) * EM_N_RATIO_MEDIUM);
	else n = (unsigned long)(expl(M) * EM_N_RATIO); 
	
	sum = harmonicSeries(n);
	G_sum = sum; /* Save for reporting later */
	
	while ((sum - M) > -DELTA && n > 0L)
	{	
		num_tests++;
		
		#ifdef VERBOSE
			diff = fabsl(M - sum);
			r = n/expl(M);
			printf(" Guess %lu, sum = %.8Lf, diff = %.8Lg, n/e^M = %.8Lg\n", n, sum, diff, r);
		#endif
		
		G_sum = sum;  /* Save for reporting later */
		sum -= 1.0L/n;   /* Faster than summing 1..n-1 again */
		n--;
	}
	
	#ifdef VERBOSE
		diff = M - sum;
		r = n/expl(M);
		printf(" Guess %lu, sum = %.8Lf, diff = %.8Lg, n/e^M = %.8Lg\n", n, sum, diff, r);
		
		printf("    Total number of guesses: %d\n\n", num_tests + 1);
	#endif
	
	return n + 1UL;
}


int main (int argc, char* argv[])
{
	if (argc != 2)
	{
		printf("Usage: %s <NUMBER>\n\n", argv[0]);
		return -1;
	}
	
	char* end; /* To capture where conversion to double stopped */
	long double M = strtold(argv[1], &end);
	
	/* Check for bad inputs */
	if (*end != '\0') /* Entire string was not successfully converted */
	{
		printf("Invalid input.\nUsage: %s <NUMBER>\n\n", argv[0]);
		return -1;
	}
	else if (M < 0.0L) /* Harmonic series can only be positive, though N = 0 is technically correct */
	{
		printf("Number must be greater than or equal to zero.\n\n");
		return -1;
	}
	
	/* Input is successful */
	unsigned long N = hseriesThreshold(M);

	printf("Sum(1/n, 1, N) > %.8Lf, when N >= %d\n\n", M, N);
	printf("Sum(1/n, 1, %lu) ~ %.8Lf\n", N, G_sum);
	
	return 0;
}