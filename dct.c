/*
e.g., cl /O2 dct.c /link /stack:2097152
*/

/* The Fourier transform (K. F. Gauss discovered it earlier) is beautiful. The discrete, bounded problem is solved once one notices that a regular polygon's radii' sum, considered as vectors, is 0. DCT is a special case. */

/*
DFT:
N complex samples zn=xn + iyn, n=└-N/2┘, ... 0, ... └(N-1)/2┘
   └(N-1)/2┘
wn = Σ zj*e^(i*j*n*2*π/N), "wholes"
   j=-└N/2┘
       └(N-1)/2┘
zn = 1/N*Σ wj*e^(-i*j*n*2*π/N) =
       j=-└N/2┘
    └(N-1)/2┘                                      └(N-1)/2┘
= 1/N*Σ[Σ zk*e^(i*k*j*2*π/N)]*e^(-i*j*n*2*π/N) = 1/N*ΣΣ zk*e^[i*j*(k - n)*2*π/N)]
    j,k=-└N/2┘                                     j,k=-└N/2┘

imagine this structure rotating (wi clocks i*k units, the sum is N*zk)
       4*z0 
    w0  w3
    /\  .\w2
   /  .  /
  / .  \/
 /.    w1
0

DCT:
N real samples xn => 2*N-1 virtual samples, even function of n=-N+1, ... N-1
    N-1                          N-1                                                    N-1
yn = Σ xj*e^[i*j*n*2*π/(2*N-1)] = Σ xj*cos[j*n*2*π/(2*N-1)] (sines cancel out.) = x0 + 2*Σ xj*cos[j*n*2*π/(2*N-1)]
   j=-N+1                       j=-N+1                                                  j=1
xn real => im(y-n) = -im(yn) & r(yn) = r(y-n)
              N-1                                             N-1
xn = 1/(2*N-1)*Σ yj*e^[-i*j*n*2*π/(2*N-1)] = 1/(2*N-1)*{y0 + 2*Σ yj*cos[j*n*2*π/(2*N-1)]} (sines cancel out & r(yj) = r(y-j))
             j=-N+1                                           j=1
                   N-1             N-1
= 1/(2*N-1)*{x0 + 2*Σ [xj + [x0 + 2*Σ xk*cos(k*j*2*π/(2*N-1))]*cos[j*n*2*π/(2*N-1)]]}
                   j=1             j=1
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define PI	3.1415926535897932384626433832795

int main(int argc, char *argv[])
{
	if (argc > 1) {
#define SAMPLE_COUNT	65536
		double a[SAMPLE_COUNT];
		int i, j, k;
		for (i = 1; i < argc; i++)
			a[i-1] = atof(argv[i]);
		{
			double b[SAMPLE_COUNT] = { 0 };
			double x;
			for (i--; i--;) {
				for (j = 1; j < argc-1; j++)
					b[i] += a[j]*cos(j*i*2*PI/(2*(argc-1)-1));
				b[i] = a[0] + 2*b[i];
			}
#define X	0
			for (i++; i < argc-1; i++) {
				x = 0;
				for (j = 1; j < argc-1-X; j++)
					x += b[j]*cos(j*i*2*PI/(2*(argc-1)-1));
				x = (b[0] + 2*x)/(2*(argc-1)-1);
				printf("%.6f ", x);
			}
		}
/*
		{
			double y, z;
			for (i = 0; i < argc-1; i++) {
				for (y = 0, j = 1; j < argc-1; j++) {
					for (z = 0, k = 1; k < argc-1; k++)
						z += a[k]*cos(k*j*2*PI/(2*(argc-1)-1));
					y += (a[0] + 2*z)*cos(j*i*2*PI/(2*(argc-1)-1)) + a[j];
				}
				y = (a[0] + 2*y)/(2*(argc-1)-1);
//				if (fabs(a[i] - (a[0] + 2*y)/(2*(argc-1)-1)) > .000000001)
//					printf("\n%d: %.16f/%.16f ", i, a[i], (a[0] + 2*y)/(2*(argc-1)-1));
			}
		}
*/
	}

	return 0;
}

