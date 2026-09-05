/* Determine the extrema of a box wrt. the proj. point, reducing the number of divisions. */

/*
sort sorted-in-cclk-order-ijk wrt. u (vec. to cube's top-left-front corner=cube org.)
sort i, j, k in cclk. order:
cclk_ijk 0: ij 1:ik 2:jk
0	0	000	210
9	16	+00 
18	32	-00
3	4	0+0 20		1=org. \2
12	20	++0	210		   /2	.-0,1
21	36	-+0	120		1-.-0
6	8	0-0 02		1=org.
15	24	+-0	021		1-.-0
24	36	--0	012		 /2
1	1	00+	10		2=org.   /1
10	17	+0+	210			  2-.-0		 /1
19	33	-0+	120		 /1,2			.-0,2	2 1
4	5	0++	210		.-0 					|/
13	21	+++	210								.-0
22	37	-++	120
7	9	0-+	102			\1
16	25	+-+	210			 .-0
25	41	--+	102			  \2
2	2	00- 01		2=org.
11	18	+0-	201		.-0,2
20	34	-0-	012		 \1
5	6	0+- 201
14	22	++- 201
23	38	-+-	012
8	10	0-- 012
17	26	+-- 021
26	42	--- 012

id3	id4		order
0	0	000			// i, j, k not all colinear (orthonormal system)
1	1	00+	120*	// * need the one=0 in the middle for case ijk_u0 = 000
2	2	00- 021*
3	4	0+0 210*
4	5	0++	210
5	6	0+- 201
6	8	0-0 012*
7	9	0-+	102
8	10	0-- 012
9	16	+00 201*		// ijk0=0
10	17	+0+	210
11	18	+0-	201
12	20	++0	210
13	21	+++	210
14	22	++- 201
15	24	+-0	021
16	25	+-+	210
17	26	+-- 021
18	32	-00 102*		// ijk0=0
19	33	-0+	120
20	34	-0-	012
21	36	-+0	120
22	37	-++	120
23	38	-+-	012
24	36	--0	012
25	41	--+	102
26	42	--- 012

\|/ <-cclk-, from extremum to extremum
 .

0: col.	0: end	bit 0: +? bit 1: -? 0=both + and - (=3?)
+:1 cclk	1: 0
-:2 clk		2: 1
			3: 2
	210	210+	210-
0	000	03		01	// case u=0, note: case one of ijk=0, the same must be in the middle/non-extremum
1	001	0123	023
2	002	0		01
3	010	023		03
4	011	0123	03
5	012	023		013
6	020	01		021
7	021	01		02
8	022	0		021
9	100	03		0
10	101	013		0
11	102	03		01
12	110	023		0
13	111	0123	0
14	112	023		01
15	120	031		021
16	121	031		02
17	122	03		021
18	200	021		0321
19	201	012		032
20	202 0		031
21	210	02		03
22	211	012		03
23	212	02		013
24	220	01		0321
25	221	01		032
26	222	0		0321

x > 0 = (unsigned)~(x-1) >> 31	// 1 =(unsigned)-x >> 31
x < 0 = x >> 31					// -1
#define SIGN(x)	((x) >> 31 | (unsigned)~((x)-1) >> 31)

code = (x <= 0) << 1 | x >= 0
#define SIGN2(x)	((unsigned)((x)-1) >> 30 & 0x2 | (unsigned)~(x) >> 31)
int sign(int x)
{
	return (unsigned)((x)-1) >> 30 & 0x2 | (unsigned)~(x) >> 31;
}
*/

#define B42(x1, x0)	((x1) << 2 | (x0))
#define B43(x2, x1, x0)	((x2) << 4 | (x1) << 2 | (x0))	/* 6-bit code */

char ijk_order[] = { 0, B43(1,2*1,0), B43(0,2*1,1), B43(2,1*1,0), B43(2,1,0), B43(2,0,1), B43(0,1*1,2), B43(1,0,2), B43(0,1,2), B43(2,0*1,1), B43(2,1,0), B43(2,0,1), B43(2,1,0), B43(2,1,0), B43(2,0,1), B43(0,2,1), B43(2,1,0), B43(0,2,1), B43(1,0*1,2), B43(1,2,0), B43(0,1,2), B43(1,2,0), B43(1,2,0), B43(0,1,2), B43(0,1,2), B43(1,0,2), B43(0,1,2) };

/* ijk/u0 +:1 -:2 0:0 */
char ijk_u0[][2] = { /* 0 */ { 3, 1 }, /* 1 */ { B43(1,2,3), B42(2,3) }, /* 2 */ { 0, 1 }, /* 3 */ { B42(2,3), 3 }, /* 4 */ { B43(1,2,3), 3 }, /* 5 */ { B42(2,3), B42(1,3) }, /* 6 */ { 1, B42(2,1) }, /* 7 */ { 1, 2 }, /* 8 */ { 0, B42(2,1) }, /* 9 */ { 3, 0 }, /* 10 */ { B42(1,3), 0 }, /* 11 */ { 3, 1 }, /* 12 */ { B42(2,3), 0 }, /* 13 */ { B43(1,2,3), 0 }, /* 14 */ { B42(2,3), 1 }, /* 15 */ { B42(3,1), B42(2,1) }, /* 16 */ { B42(3,1), 2 }, /* 17 */ { 3, B42(2,1) }, /* 18 */ { B42(2,1), B43(3,2,1) }, /* 19 */ { B42(1,2), B42(3,2) }, /* 20 */ { 0, B42(3,1) }, /* 21 */ { 2, 3 }, /* 22 */ { B42(1,2), 3 }, /* 23 */ { 2, B42(1,3) }, /* 24 */ { 1, B43(3,2,1) }, /* 25 */ { 1, B42(3,2) }, /* 26 */ { 0, B43(3,2,1) } };

#define VECMOV(u, v, i)	do { u[i] = v[i]; u[2] = v[2]; } while (0)
#define VECADD(u, v, i)	do { u[i] += v[i]; u[2] += v[2]; } while (0)

#define SIGN3(x)	((unsigned)(x) >> 30 & 0x2 | (unsigned)-(x) >> 31)

int sign3(int x)
{
	return (unsigned)x >> 30 & 0x2 | (unsigned)-x >> 31;	/* -x = ~(x-1) */
}

/* u->v cclk */
#define cclk(u, v, i)	(v[2]*u[i] - v[i]*u[2])	/* a function would be perhaps better */

/* test */
int ijk[3][3] = { { 0, -2, 0 }, { 2, 0, 0 }, { 0, 0, 2 } };

#include <stdio.h>

int main()
{
	int i;
	int u0[3] = { -2, -1, 1 } /* test */, lrtb[4][3];	/* need be [][3] because of cclk */

	int ijk_seq[2] = { ijk_order[sign3(cclk(ijk[1], ijk[2], 0)) * 9 +
								 sign3(cclk(ijk[0], ijk[2], 0)) * 3 +
								 sign3(cclk(ijk[0], ijk[1], 0))],
					   ijk_order[sign3(cclk(ijk[1], ijk[2], 1)) * 9 +
								 sign3(cclk(ijk[0], ijk[2], 1)) * 3 +
								 sign3(cclk(ijk[0], ijk[1], 1))] },
		/* proj. vecs. 0, 1, 2 (proj. i, j, k in cclk. order) wrt. u0 */
		k[] = { sign3(cclk(u0, ijk[ijk_seq[0] >> 4], 0)) * 9 +
				sign3(cclk(u0, ijk[ijk_seq[0] >> 2 & 0x3], 0)) * 3 +
				sign3(cclk(u0, ijk[ijk_seq[0] & 0x3], 0)),
				sign3(cclk(u0, ijk[ijk_seq[1] >> 4], 1)) * 9 +
				sign3(cclk(u0, ijk[ijk_seq[1] >> 2 & 0x3], 1)) * 3 +
				sign3(cclk(u0, ijk[ijk_seq[1] & 0x3], 1)) };
	{
		int j, n;

		i = 3;	/* bottom=3, top=2, right=1, left=0 */
		do {
			VECMOV(lrtb[i], u0, i >> 1);
			j = ijk_u0[k[i >> 1]][i&0x1];

			if (j) {
				n = ijk_seq[i >> 1] >> ((j & 0x3) - 1)*2 & 0x3;
				do {
					VECADD(lrtb[i], ijk[n], i >> 1);
				} while ((j >>= 2) &&
						(-(i&0x1) ^
						cclk(lrtb[i], ijk[n = ijk_seq[i >> 1] >> ((j & 0x3) - 1)*2 & 0x3], i >> 1)) +
						(i&0x1) >= 0);	/* =0: for case proj. i/j/k colinear & opposite to u0 */
			}
		} while (i--);
	}
	printf("l: %d,%d\nr: %d,%d\nt: %d,%d\nb: %d,%d", lrtb[0][0], lrtb[0][2],
													 lrtb[1][0], lrtb[1][2],
													 lrtb[2][1], lrtb[2][2],
													 lrtb[3][1], lrtb[3][2]);
	return 0;
}
