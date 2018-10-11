#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include<conio.h>
#include <malloc.h>
#include<omp.h>
#include <cstdlib>
#include <ctime>
#include <Windows.h>

#define  Max(a,b) ((a)>(b)?(a):(b))
int  N; /*20*/
double   maxeps = 0.1e-7;
int itmax = 10;
int i, j, k;
double eps;

double* A;
double* B;

int do_it = 1;
int processor_amount;

void relax();
void resid();
void init();
void verify();

#ifdef __GLIBC__
#if !defined(_POSIX_C_SOURCE) && !defined(_POSIX_SOURCE)
typedef struct _SYSTEMTIME {
	WORD wYear;
	WORD wMonth;
	WORD wDayOfWeek;
	WORD wDay;
	WORD wHour;
	WORD wMinute;
	WORD wSecond;
	WORD wMilliseconds;
} SYSTEMTIME, *PSYSTEMTIME;
#endif
#endif


int get_value() {

	char ch = 'a';
	int a = 0;
	while (ch != '\n')
	{
		ch = getchar();
		if (((ch >= '0') && (ch <= '9')) || (ch == '\n'))
		{
			if ((ch >= '0') && (ch <= '9'))
				a = a * 10 + (ch - '0');
		}
		else
		{
			printf("error");
			goto towards_end;
		}
	}
	return a;
towards_end:
	return -1000;
}

void wtime(double *t) {
	static int sec = -1;
	//struct timeval /*TIMEVAL*/ tv;
	SYSTEMTIME tv;
	/*gettimeofday(&tv, (void *)0);*/
	GetLocalTime(&tv);
	if (sec < 0) sec = tv.wSecond;
	*t = (tv.wSecond - sec) + 1.0e-6*tv.wMilliseconds;
}

SYSTEMTIME wtimeMiliSeconds() {
	SYSTEMTIME tv;
	/*gettimeofday(&tv, (void *)0);*/
	GetLocalTime(&tv);
	return tv;
}

int main() {
	int it, thread_number;
	double time0, time1; 

	printf("Set itmax\n");
	itmax = get_value();

	printf("Set processor amount\n");
	processor_amount = get_value();

	printf("Set N - array size\n");
	N = get_value();

	if ((itmax >= 0) && (processor_amount >= 0) && (N >= 0)) {
		init();
		if (do_it == 1) {
		//	time0 = clock();
			
			time0 = omp_get_wtime (); 
			wtime(&time0);
			SYSTEMTIME curTimeStart = wtimeMiliSeconds();

#pragma omp parallel num_threads(processor_amount)
			{
#pragma omp for schedule(dynamic)
				for (it = 1; it <= itmax; it++) {
					eps = 0.;
					relax();
					resid();
					printf("it=%4i eps=%f\n", it, eps);
					thread_number = omp_get_thread_num();
					printf("Working_thread_number: %d\n\n", thread_number);
					if (eps < maxeps) break;
				}
			}

			//time1 = clock();
			time1 = omp_get_wtime (); 
			SYSTEMTIME curTimeEnd = wtimeMiliSeconds();
			wtime(&time1);

			printf("Time in seconds = %d\t", time1 - time0/* / CLOCKS_PER_SEC*/);
			printf("Time in miliseconds = %d\t", curTimeEnd.wMilliseconds = curTimeStart.wMilliseconds);
			verify();
			free(A);
			free(B);
		}
	}

	_getche();
	system("pause");
	return 0;

}

void init() {
	A = (double*)malloc(N*N*N * sizeof(double));
	B = (double*)malloc(N*N*N * sizeof(double));

	printf("Memory need:%d Bytes\n", 2 * N*N*N * sizeof(double));

	if ((A != NULL) && (B != NULL)) {
		#pragma omp parallel num_threads(processor_amount) 
		{
		#pragma omp for //
		for (i = 0; i < N; i++)
			for (j = 0; j < N; j++)
				for (k = 0; k < N; k++) {
					double element = A[i*N*N + j*N + k];
					if (i == 0 || i == N - 1 || j == 0 || j == N - 1 || k == 0 || k == N - 1) A[i*N*N + j*N + k] = 0;
					else A[i*N*N + j*N + k] = (4. + i + j + k);
				}
			}
		}
	else {
		printf("There is not enough memory for arrays A and B");
		do_it = 0;
	}
}

void relax()
{
	#pragma omp parallel num_threads(processor_amount)
	{
		#pragma omp for schedule(static, processor_amount)
		for (i = 1; i <= N - 2; i++)
			for (j = 1; j <= N - 2; j++)
				for (k = 1; k <= N - 2; k++) {
					B[i*N*N + j*N + k] = (A[(i - 1)*N*N + j*N + k] + A[(i + 1)*N*N + j*N + k] + A[i*N*N + (j - 1)*N + k] + A[i*N*N + (j + 1)*N + k] + A[i*N*N + j*N + k - 1] + A[i*N*N + j*N + k + 1]) / 6.;
				}
	}
}

void resid() {
	int thread;
	#pragma omp parallel num_threads(1)
	{
		#pragma omp for schedule(static, processor_amount)
		for (i = 1; i <= N - 2; i++)
			for (j = 1; j <= N - 2; j++)
				for (k = 1; k <= N - 2; k++) {
					double e;
					e = fabs(A[i*N*N + j*N + k] - B[i*N*N + j*N + k]);
					A[i*N*N + j*N + k] = B[i*N*N + j*N + k];
					eps = Max(eps, e);
				}
	}
}

void verify() {
	double s;
	s = 0.;
	#pragma omp parallel reduction(+:s) num_threads(processor_amount)
	{
		#pragma omp for
		for (i = 0; i <= N - 1; i++)
			for (j = 0; j <= N - 1; j++)
				for (k = 0; k <= N - 1; k++)
				{
					s = s + A[i*N*N + j*N + k] * (i + 1)*(j + 1)*(k + 1) / (N*N*N);
				}
	}
	printf(" S = %f\n", s);
}



