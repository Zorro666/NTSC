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

/*  svd.c -- Singular value decomposition. Translated to 'C' from the
 *           original Algol code in "Handbook for Automatic Computation,
 *           vol. II, Linear Algebra", Springer-Verlag.
 *
 *  (C) 2000, C. Bond. All rights reserved.
 *
 *  This is almost an exact translation from the original, except that
 *  an iteration counter is added to prevent stalls. This corresponds
 *  to similar changes in other translations.
 *
 *  Returns an error code = 0, if no errors and 'k' if a failure to
 *  converge at the 'kth' singular value.
 * 
 */

int svd(int m,int n,int withu,int withv,float eps,float tol,
	float **a,float *q,float **u,float **v)
{
	int i,j,k,l,l1,iter,retval;
	float c,f,g,h,s,x,y,z;
	float *e;

	e = (float *)calloc((size_t)n,sizeof(float));
	retval = 0;

	 l = 0;

/* Copy 'a' to 'u' */    
	for (i=0;i<m;i++) {
		for (j=0;j<n;j++)
			u[i][j] = a[i][j];
	}
/* Householder's reduction to bidiagonal form. */
	g = x = 0.0;    
	for (i=0;i<n;i++) {
		e[i] = g;
		s = 0.0;
		l = i+1;
		for (j=i;j<m;j++)
			s += (u[j][i]*u[j][i]);
		if (s < tol)
			g = 0.0;
		else {
			f = u[i][i];
			g = (f < 0) ? sqrtf(s) : -sqrtf(s);
			h = f * g - s;
			u[i][i] = f - g;
			for (j=l;j<n;j++) {
				s = 0.0;
				for (k=i;k<m;k++)
					s += (u[k][i] * u[k][j]);
				f = s / h;
				for (k=i;k<m;k++)
					u[k][j] += (f * u[k][i]);
			} /* end j */
		} /* end s */
		q[i] = g;
		s = 0.0;
		for (j=l;j<n;j++)
			s += (u[i][j] * u[i][j]);
		if (s < tol)
			g = 0.0;
		else {
			f = u[i][i+1];
			g = (f < 0) ? sqrtf(s) : -sqrtf(s);
			h = f * g - s;
			u[i][i+1] = f - g;
			for (j=l;j<n;j++) 
				e[j] = u[i][j]/h;
			for (j=l;j<m;j++) {
				s = 0.0;
				for (k=l;k<n;k++) 
					s += (u[j][k] * u[i][k]);
				for (k=l;k<n;k++)
					u[j][k] += (s * e[k]);
			} /* end j */
		} /* end s */
		y = fabsf(q[i]) + fabsf(e[i]);                         
		if (y > x)
			x = y;
	} /* end i */

/* accumulation of right-hand transformations */
	if (withv) {
		for (i=n-1;i>=0;i--) {
			if ((g > 0.0f) || (g < 0.0f)){
				h = u[i][i+1] * g;
				for (j=l;j<n;j++)
					v[j][i] = u[i][j]/h;
				for (j=l;j<n;j++) {
					s = 0.0;
					for (k=l;k<n;k++) 
						s += (u[i][k] * v[k][j]);
					for (k=l;k<n;k++)
						v[k][j] += (s * v[k][i]);

				} /* end j */
			} /* end g */
			for (j=l;j<n;j++)
				v[i][j] = v[j][i] = 0.0;
			v[i][i] = 1.0;
			g = e[i];
			l = i;
		} /* end i */
 
	} /* end withv, parens added for clarity */

/* accumulation of left-hand transformations */
	if (withu) {
		for (i=n;i<m;i++) {
			for (j=n;j<m;j++)
				u[i][j] = 0.0;
			u[i][i] = 1.0;
		}
	}
	if (withu) {
		for (i=n-1;i>=0;i--) {
			l = i + 1;
			g = q[i];
			for (j=l;j<m;j++)  /* upper limit was 'n' */
				u[i][j] = 0.0f;
			if ((g > 0.0f) || (g < 0.0f)){
				h = u[i][i] * g;
				for (j=l;j<m;j++) { /* upper limit was 'n' */
					s = 0.0f;
					for (k=l;k<m;k++)
						s += (u[k][i] * u[k][j]);
					f = s / h;
					for (k=i;k<m;k++) 
						u[k][j] += (f * u[k][i]);
				} /* end j */
				for (j=i;j<m;j++) 
					u[j][i] /= g;
			} /* end g */
			else {
				for (j=i;j<m;j++)
					u[j][i] = 0.0f;
			}
			u[i][i] += 1.0f;
		} /* end i*/
	} /* end withu, parens added for clarity */

/* diagonalization of the bidiagonal form */
	eps *= x;
	for (k=n-1;k>=0;k--) {
		iter = 0;
test_f_splitting:
		for (l=k;l>=0;l--) {
			if (fabsf(e[l]) <= eps) goto test_f_convergence;
			if (fabsf(q[l-1]) <= eps) goto cancellation;
		} /* end l */

/* cancellation of e[l] if l > 0 */
cancellation:
		c = 0.0;
		s = 1.0;
		l1 = l - 1;
		for (i=l;i<=k;i++) {
			f = s * e[i];
			e[i] *= c;
			if (fabsf(f) <= eps) goto test_f_convergence;
			g = q[i];
			h = q[i] = sqrtf(f*f + g*g);
			c = g / h;
			s = -f / h;
			if (withu) {
				for (j=0;j<m;j++) {
					y = u[j][l1];
					z = u[j][i];
					u[j][l1] = y * c + z * s;
					u[j][i] = -y * s + z * c;
				} /* end j */
			} /* end withu, parens added for clarity */
		} /* end i */
test_f_convergence:
		z = q[k];
		if (l == k) goto convergence;

/* shift from bottom 2x2 minor */
		iter++;
		if (iter > 30) {
			retval = k;
			break;
		}
		x = q[l];
		y = q[k-1];
		g = e[k-1];
		h = e[k];
		f = ((y-z)*(y+z) + (g-h)*(g+h)) / (2*h*y);
		g = sqrtf(f*f + 1.0f);
		f = ((x-z)*(x+z) + h*(y/((f<0)?(f-g):(f+g))-h))/x;
/* next QR transformation */
		c = s = 1.0;
		for (i=l+1;i<=k;i++) {
			g = e[i];
			y = q[i];
			h = s * g;
			g *= c;
			e[i-1] = z = sqrtf(f*f+h*h);
			c = f / z;
			s = h / z;
			f = x * c + g * s;
			g = -x * s + g * c;
			h = y * s;
			y *= c;
			if (withv) {
				for (j=0;j<n;j++) {
					x = v[j][i-1];
					z = v[j][i];
					v[j][i-1] = x * c + z * s;
					v[j][i] = -x * s + z * c;
				} /* end j */
			} /* end withv, parens added for clarity */
			q[i-1] = z = sqrtf(f*f + h*h);
			c = f/z;
			s = h/z;
			f = c * g + s * y;
			x = -s * g + c * y;
			if (withu) {
				for (j=0;j<m;j++) {
					y = u[j][i-1];
					z = u[j][i];
					u[j][i-1] = y * c + z * s;
					u[j][i] = -y * s + z * c;
				} /* end j */
			} /* end withu, parens added for clarity */
		} /* end i */
		e[l] = 0.0;
		e[k] = f;
		q[k] = x;
		goto test_f_splitting;
convergence:
		if (z < 0.0) {
/* q[k] is made non-negative */
			q[k] = - z;
			if (withv) {
				for (j=0;j<n;j++)
					v[j][k] = -v[j][k];
			} /* end withv, parens added for clarity */
		} /* end z */
	} /* end k */
	
	free(e);
	return retval;
}
