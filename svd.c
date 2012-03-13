/*
 * svdcomp - SVD decomposition routine.
 * Takes an mxn matrix a and decomposes it into udv, where u,v are
 * left and right orthogonal transformation matrices, and d is a
 * diagonal matrix of singular values.
 *
 * This routine is adapted from svdecomp.c in XLISP-STAT 2.1 which is
 * code from Numerical Recipes adapted by Luke Tierney and David Betz.
 *
 * Input to svd_decompose is as follows:
 *   a = mxn matrix to be decomposed, gets overwritten with u
 *   m = row dimension of a
 *   n = column dimension of a
 *   w = returns the vector of singular values of a
 *   v = returns the right orthogonal transformation matrix
*/

#include <math.h>
#include <stdio.h>
#include <malloc.h>

/* returns ABS(A)*s', where s is +1 if `B >= 0', -1 otherwise. */
static float svd_sign(float a, float b)
{
	const float fabsa = fabsf(a);
	if (b < 0.0f)
	{
		return -fabsa;
	}
	return fabsa;
}

static float svd_max(const float a, const float b)
{
	if (a > b)
	{
		return a;
	}
	return b;
}

static float svd_pythag(float a, float b)
{
	float at = fabsf(a);
	float bt = fabsf(b);
	float ct;
	float result;

	if (at > bt)
	{
		ct = bt / at;
		result = at * sqrtf(1.0f + ct * ct);
	}
	else if (bt > 0.0f)
	{
		ct = at / bt;
		result = bt * sqrtf(1.0f + ct * ct);
	}
	else
	{
		result = 0.0f;
	}
	return(result);
}

int svd_decompose(float** a, unsigned int M, unsigned int N, float* w, float** v)
{
	const int m = (int)M;
	const int n = (int)N;
	int flag, i, its, j, jj, k, l, nm;
	float c, f, h, s, x, y, z;
	float anorm = 0.0f, g = 0.0f, scale = 0.0f;
	float *rv1;

	if (m < n)
	{
		fprintf(stderr, "#rows must be > #cols \n");
		return(0);
	}

	rv1 = (float *)malloc((unsigned int) n*sizeof(float));

	/* Householder reduction to bidiagonal form */
	for (i = 0; i < n; i++)
	{
		/* left-hand reduction */
		l = i + 1;
		rv1[i] = scale * g;
		g = s = scale = 0.0f;
		if (i < m)
		{
			for (k = i; k < m; k++)
			{
				scale += fabsf((float)a[k][i]);
			}
			if (scale > 0.0f)
			{
				for (k = i; k < m; k++)
				{
					a[k][i] = (float)((float)a[k][i]/scale);
					s += ((float)a[k][i] * (float)a[k][i]);
				}
				f = (float)a[i][i];
				g = -svd_sign(sqrtf(s), f);
				h = f * g - s;
				a[i][i] = (float)(f - g);
				if (i != n - 1)
				{
					for (j = l; j < n; j++)
					{
						for (s = 0.0f, k = i; k < m; k++)
						{
							s += ((float)a[k][i] * (float)a[k][j]);
						}
						f = s / h;
						for (k = i; k < m; k++)
						{
							a[k][j] += (float)(f * (float)a[k][i]);
						}
					}
				}
				for (k = i; k < m; k++)
				{
					a[k][i] = (float)((float)a[k][i]*scale);
				}
			}
		}
		w[i] = (float)(scale * g);

		/* right-hand reduction */
		g = s = scale = 0.0f;
		if (i < m && i != n - 1)
		{
			for (k = l; k < n; k++)
			{
				scale += fabsf((float)a[i][k]);
			}
			if (scale > 0.0f)
			{
				for (k = l; k < n; k++)
				{
					a[i][k] = (float)((float)a[i][k]/scale);
					s += ((float)a[i][k] * (float)a[i][k]);
				}
				f = (float)a[i][l];
				g = -svd_sign(sqrtf(s), f);
				h = f * g - s;
				a[i][l] = (float)(f - g);
				for (k = l; k < n; k++)
				{
					rv1[k] = (float)a[i][k] / h;
				}
				if (i != m - 1)
				{
					for (j = l; j < m; j++)
					{
						for (s = 0.0f, k = l; k < n; k++)
						{
							s += ((float)a[j][k] * (float)a[i][k]);
						}
						for (k = l; k < n; k++)
						{
							a[j][k] += (float)(s * rv1[k]);
						}
					}
				}
				for (k = l; k < n; k++)
				{
					a[i][k] = (float)((float)a[i][k]*scale);
				}
			}
		}
		anorm = svd_max(anorm, (fabsf((float)w[i]) + fabsf(rv1[i])));
	}

	/* accumulate the right-hand transformation */
	for (i = n - 1; i >= 0; i--)
	{
		if (i < n - 1)
		{
			if (g > 0.0f)
			{
				for (j = l; j < n; j++)
				{
					v[j][i] = (float)(((float)a[i][j] / (float)a[i][l]) / g);
				}
				/* double division to avoid underflow */
				for (j = l; j < n; j++)
				{
					for (s = 0.0f, k = l; k < n; k++)
					{
						s += ((float)a[i][k] * (float)v[k][j]);
					}
					for (k = l; k < n; k++)
					{
						v[k][j] += (float)(s * (float)v[k][i]);
					}
				}
			}
			for (j = l; j < n; j++)
			{
				v[i][j] = v[j][i] = 0.0f;
			}
		}
		v[i][i] = 1.0f;
		g = rv1[i];
		l = i;
	}

	/* accumulate the left-hand transformation */
	for (i = n - 1; i >= 0; i--)
	{
		l = i + 1;
		g = (float)w[i];
		if (i < n - 1)
		{
			for (j = l; j < n; j++)
			{
				a[i][j] = 0.0f;
			}
		}
		if (g > 0.0f)
		{
			g = 1.0f / g;
			if (i != n - 1)
			{
				for (j = l; j < n; j++)
				{
					for (s = 0.0f, k = l; k < m; k++)
					{
						s += ((float)a[k][i] * (float)a[k][j]);
					}
					f = (s / (float)a[i][i]) * g;
					for (k = i; k < m; k++)
					{
						a[k][j] += (float)(f * (float)a[k][i]);
					}
				}
			}
			for (j = i; j < m; j++)
			{
				a[j][i] = (float)((float)a[j][i]*g);
			}
		}
		else
		{
			for (j = i; j < m; j++)
			{
				a[j][i] = 0.0f;
			}
		}
		++a[i][i];
	}

	/* diagonalize the bidiagonal form */
	for (k = n - 1; k >= 0; k--)
	{                             /* loop over singular values */
		for (its = 0; its < 30; its++)
		{                         /* loop over allowed iterations */
			flag = 1;
			for (l = k; l >= 0; l--)
			{                     /* test for splitting */
				float temp;
				nm = l - 1;
				temp = fabsf(rv1[l]) + anorm;
				temp -= anorm;
				if (((temp > 0.0f) || (temp < 0.0f)) == 0)
				{
					flag = 0;
					break;
				}
				temp = (fabsf((float)w[nm]) + anorm);
				temp -= anorm;
				if (((temp > 0.0f) || (temp < 0.0f)) == 0)
				{
					break;
				}
			}
			if (flag)
			{
				c = 0.0f;
				s = 1.0f;
				for (i = l; i <= k; i++)
				{
					float temp;
					f = s * rv1[i];
					temp = (fabsf(f) + anorm);
				 	temp -=	anorm;
					if ((temp < 0.0f) || (temp > 0.0f))
					{
						g = (float)w[i];
						h = svd_pythag(f, g);
						w[i] = (float)h;
						h = 1.0f / h;
						c = g * h;
						s = (- f * h);
						for (j = 0; j < m; j++)
						{
							y = (float)a[j][nm];
							z = (float)a[j][i];
							a[j][nm] = (float)(y * c + z * s);
							a[j][i] = (float)(z * c - y * s);
						}
					}
				}
			}
			z = (float)w[k];
			if (l == k)
			{                  /* convergence */
				if (z < 0.0f)
				{              /* make singular value nonnegative */
					w[k] = (float)(-z);
					for (j = 0; j < n; j++)
					{
						v[j][k] = (-v[j][k]);
					}
				}
				break;
			}
			if (its >= 30)
			{
				free((void*) rv1);
				fprintf(stderr, "No convergence after 30,000! iterations \n");
				return(0);
			}

			/* shift from bottom 2 x 2 minor */
			x = (float)w[l];
			nm = k - 1;
			y = (float)w[nm];
			g = rv1[nm];
			h = rv1[k];
			f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0f * h * y);
			g = svd_pythag(f, 1.0f);
			f = ((x - z) * (x + z) + h * ((y / (f + svd_sign(g, f))) - h)) / x;

			/* next QR transformation */
			c = s = 1.0f;
			for (j = l; j <= nm; j++)
			{
				i = j + 1;
				g = rv1[i];
				y = (float)w[i];
				h = s * g;
				g = c * g;
				z = svd_pythag(f, h);
				rv1[j] = z;
				c = f / z;
				s = h / z;
				f = x * c + g * s;
				g = g * c - x * s;
				h = y * s;
				y = y * c;
				for (jj = 0; jj < n; jj++)
				{
					x = (float)v[jj][j];
					z = (float)v[jj][i];
					v[jj][j] = (float)(x * c + z * s);
					v[jj][i] = (float)(z * c - x * s);
				}
				z = svd_pythag(f, h);
				w[j] = (float)z;
				if (z > 0.0f)
				{
					z = 1.0f / z;
					c = f * z;
					s = h * z;
				}
				f = (c * g) + (s * y);
				x = (c * y) - (s * g);
				for (jj = 0; jj < m; jj++)
				{
					y = (float)a[jj][j];
					z = (float)a[jj][i];
					a[jj][j] = (float)(y * c + z * s);
					a[jj][i] = (float)(z * c - y * s);
				}
			}
			rv1[l] = 0.0f;
			rv1[k] = f;
			w[k] = (float)x;
		}
	}
	free((void*) rv1);
	return(1);
}
