#include <stdio.h> //Used for standard I/O.
#include <stdlib.h> //Used for exit().
#include <math.h> //Used for sqrt and log functions.
#include <sys/time.h>//Used for time calculations.
#include <unistd.h>//Used for getpid(),sleep(),fork() etc.
#include <sys/wait.h>//Used for waitpid().
#include <string.h>//Used for argument handling (strcmp()).
#include <omp.h>
double get_wtime(void){ //Function for timestamps.
    struct timeval t;
    gettimeofday(&t, NULL);
    return (double)t.tv_sec + (double)t.tv_usec*1.0e-6;
}

double f(double x){ //F function used for calculation of log*sqrt.
       return log(x)*sqrt(x);
}

double proc_calc(int ifk, int procn,double a, double b){ //Child process function.
    unsigned long const n = 1e9;//Every process calculates a part of the result for quicker solving times
    const double dx = (b-a)/n;
    double S = 0;

    for (unsigned long i = ifk; i < (n); i+=procn) { //i+=process_number for division of the total addition.
         double xi = a + (i + 0.5)*dx;
            S += f(xi);
    }
            S *= dx;
            return S;
}

int main(int argc, char *argv[]){
	int threads_numb;
	double a,b,t0,t1;
	double result;
	int thread_id;
	system("clear");
	printf("Integral log(x)*sqrt(x) with openmp threads\nGive number of lower limit\n");
	scanf("%lf", &a);
	printf("\nGive upper limit\n");
	scanf("%lf", &b);
	printf("\nGive number of threads\n");
	scanf("%d", &threads_numb);
	t0 = get_wtime();
	#pragma omp parallel reduction (+:result)num_threads(threads_numb)
	{
		thread_id=omp_get_num_threads();
		result=0.0;
		result = proc_calc(thread_id,threads_numb,a,b);
	}
	t1=get_wtime();

	printf("\nResult of integral log(x)*sqrt(x) with limits %.4f-%.4f is %.6f", a , b , result);
	printf("\nTotal time taken %.4f", t1-t0);
	return 0;
}

