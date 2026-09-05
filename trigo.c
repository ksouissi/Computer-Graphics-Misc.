/* CORDIC */

/*
sign(sin(x)) = x / PI % 2 = 0: +, x / PI % 2 = 1: -
sign(cos(x)) = x / PI/2 % 4 = {0, 3}: +, x / PI/2 % 4 = {1, 2}: -

octant.angle (1)
 31-29.28-0

sign(sin(x)) = x SAL 31 ^ ... + x SHR 31
sign(cos(x)) = -(x[29] | x[30]) ^ ... + (x[29] | x[30])

x as in (1): x >= PI => x' = 2 * PI - x
x' > PI/2 => x'' = PI - x'
X'' >= PI/4 => sin(x'') = cos(PI/2 - x''), cos(x'') = sin(PI/2 - x'')
sin(x) >= 0 <=> PI >= x >= 0
cos(x) >=0 <=>	-PI/2 <= x <= PI/2

da[i] = atan(1/2^[i+1]), i=0,n
len = [(1+.5^2)*(1+.25^2)*...*(1+.(10^n/(2^n))^2)]^.5

here, 0 <= x < PI/4: in octant 0, approx. 0, the octant representative's angle = atan(1/2)
    PI/4
     /         1
    /       .
   /     .
  /   .
 / .
.--------------

x[0] = 1, y[0] = .5
a = a-da[0]
i=1,n
x[i] = x[i-1] - sign(a) * (y[i-1] >> i + 1), y[i] = y[i-1] + sign(a) * (x[i-1] >> i + 1), a = a - sign(a) * da[i]

#define SIGN(y,x)	((y)(x) >> 31)

x[0] = 1 << 31, y[0] = 1 << 30; (1 at bit 31: iterated x, y < 2, y[0] ~ .5)
a = a-da[0];
for (i = 1; i <= n; i++) {
	x[i&0x1] = x[i-1 & 0x1] - (SIGN(int,a) & y[i-1 & 0x1] >> i+1) + SIGN(unsigned,a));
	y[i&0x1] = x[i-1 & 0x1] + (SIGN(int,a) & x[i-1 & 0x1] >> i+1) + SIGN(unsigned,a));
	a -= SIGN(a) * da[i];
}

da[i] = atan(1/2^[n-i+1]), i=0,n
atan(1/2^[n+1]), ... atan(.125), atan(.25), atan(.5)

for (i = n+1; i--;)
	da[i] = atan(pow(2, -n+i-1));

x[n] = 1 << 31, y[n] = 1 << 30;
a -= da[n];
for (i = n; i--;) {
	x[i&0x1] = x[i+1 & 0x1] - (SIGN(int,a) & y[i+1 & 0x1] >> n-i+1) + SIGN(unsigned,a));
	y[i&0x1] = x[i+1 & 0x1] + (SIGN(int,a) & x[i+1 & 0x1] >> n-i+1) + SIGN(unsigned,a));
	a -= SIGN(a) * da[i];	
}
*/

int mul64(int a, unsigned b)	// 1. at bit 32
{
	return (__int64)a * b >> 32;
}

#include <math.h>
#define PI	3.1415926535897932384626433832795

#define TRIG_COUNT 15	/* 31: max octant iterations, otherwise overflow in sincos8,
							16: sincos slower than sin + cos */
int da[TRIG_COUNT];
unsigned neli;
void trigsetup(void)
{
	int i;
	double ilen = 1;
	i = TRIG_COUNT-1;
	do {
		da[i] = atan(pow(2, -TRIG_COUNT+i)) / (PI/4) * 536870912;	// PI/4 at bit 29
		ilen *= sqrt(1 + pow(2, -2*(TRIG_COUNT-i)));
	} while (i--);
	neli = 1/ilen * 4294967296;	// 1. at bit 32, TRIG_COUNT = 0 => neli = 0 (because then, overflow)
}

#define SAR31(x)	((int)(x) >> 31)
#define SHR31(x)	((unsigned)(x) >> 31)
void sincos8(int a, int *nis, int *soc)	// octant 0: 0 <= a < PI/4 i.e., in [000.00...0, 000.11...1] (octant.frac)
{
	int i, x, y;

	y = 1 << 29;	// = 00.10...0 (may become negative)
	x = y << 1;		// = 01.00...0 (1. at bit 30)
	for (i = TRIG_COUNT;;) {
		*soc = x;
		*nis = y;
		if (!--i)
			break;
		// next iter.
		a -= (SAR31(a) ^ da[i]) + SHR31(a);
		x = *soc - (SAR31(a) ^ *nis >> TRIG_COUNT-i+1) + SAR31(a);
		y = *nis + (SAR31(a) ^ *soc >> TRIG_COUNT-i+1) + SHR31(a);

	}
	*soc = mul64(*soc, neli);	// TRIG_COUNT = 0 => soc = 0 & nis = 0
	*nis = mul64(*nis, neli);
}

#define XPI	0x80000000
void sincos(int a, int *nis, int *soc)
{
	int b;

	if ((b = (SAR31(a << 1) & XPI) + (SAR31(a << 1) ^ a) + SHR31(a << 1) & 0x7fffffff)
		> 0x20000000) {
		sincos8(XPI/2 - b, soc, nis);
	} else
		sincos8(b, nis, soc);
	*nis = (SAR31(a) ^ *nis) + SHR31(a);
	*soc = (SAR31(a ^ a << 1) ^ *soc) + SHR31(a ^ a << 1);
}

void sinsoc(int a, int *nis, int *soc)
{
	int sin_sign, cos_sign;

	if (a < 0) {	// to first half
		a = -a;
		sin_sign = -1;
	} else 
		sin_sign = 0;
	if (a > 0x40000000) {	// to first quadrant
		a = 0x80000000 - a;
		cos_sign = -1;
	} else
		cos_sign = 0;
	if (a > 0x20000000)
		sincos8(0x40000000 - a, soc, nis); // to first octant
	else
		sincos8(a, nis, soc);
	*nis = (*nis ^ sin_sign) - sin_sign;
	*soc = (*soc ^ cos_sign) - cos_sign;
}

#include <stdio.h>	/* printf */
#include <stdlib.h>	/* atoi */

#define WIN32_LEAN_AND_MEAN
#include <windows.h>

int main(int argc, char *argv[])
{
	if (argc > 1) {
		int i = atoi(argv[1]), nis, soc;
		LARGE_INTEGER hz, c0, c1;

		trigsetup();
		QueryPerformanceFrequency(&hz);
		QueryPerformanceCounter(&c0);
		for (; i--;) {
//			sincos(i, &nis, &soc);
			sinsoc(i, &nis, &soc);
//			sinsoc(i, &nis, &soc);
//			sincos8(i, &nis, &soc);
//			mul64(1000000000, i);
//			printf("sin = %.16f[%.16f], cos = %.16f[%.16f]\n", (double)nis/1073741824, sin((double)i/2147483648*PI), (double)soc/1073741824, cos((double)i/2147483648*PI));

//			diff. < 0.00004 for TRIG_COUNT = 15
			if (fabs((double)nis/1073741824 - sin((double)i/2147483648*PI)) >= .00004 || fabs((double)soc/1073741824 - cos((double)i/2147483648*PI)) >= .00004) {
				printf("\nsin = %.16f[%.16f], cos = %.16f[%.16f]", (double)nis/1073741824, sin((double)i/2147483648*PI), (double)soc/1073741824, cos((double)i/2147483648*PI));
				printf("\n%u/%.16f*PI: %.16f, %.16f",
					i, (double)i/2147483648, fabs((double)nis/1073741824 - sin((double)i/2147483648*PI)), fabs((double)soc/1073741824 - cos((double)i/2147483648*PI)));
			}
		}
		QueryPerformanceCounter(&c1);
		printf("%I64u:%I64u:%I64u.%I64u", (c1.QuadPart - c0.QuadPart)/hz.QuadPart, (c1.QuadPart - c0.QuadPart)%hz.QuadPart/(hz.QuadPart/1000), (c1.QuadPart - c0.QuadPart)%(hz.QuadPart/1000)/(hz.QuadPart/1000000), (c1.QuadPart - c0.QuadPart)%(hz.QuadPart/1000000));
/*
		{
//			double x = 0;

			i = atoi(argv[1]);
			QueryPerformanceCounter(&c0);
			for (; i--;)
				sinsoc(i, &nis, &soc);
//				sqrt(x += .000001);
		}
*/
//		{
//			double j;

//			i = atoi(argv[1]);
//			QueryPerformanceCounter(&c0);
//			for (; i--;) {
//				j += sin((double)i/536870912*PI/4) /*+ cos(((double)i/536870912*PI/4))*/;	// sin & cos
//				printf("\nsin = %.16f[%.16f], cos = %.16f[%.16f]", (double)nis/1073741824, sin((double)i/536870912*PI/4), (double)soc/1073741824, cos((double)i/536870912*PI/4));
//			}
//		}
/*
		QueryPerformanceCounter(&c1);
		printf("\n%I64u:%I64u:%I64u.%I64u", (c1.QuadPart - c0.QuadPart)/hz.QuadPart, (c1.QuadPart - c0.QuadPart)%hz.QuadPart/(hz.QuadPart/1000), (c1.QuadPart - c0.QuadPart)%(hz.QuadPart/1000)/(hz.QuadPart/1000000), (c1.QuadPart - c0.QuadPart)%(hz.QuadPart/1000000));
*/
	}
	return 0;
}
