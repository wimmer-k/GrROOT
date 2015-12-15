////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////
////////////                       GrROOT
////////////
////////////          Purpose:
////////////                   To assist in the analysis of data from
////////////                 the gretina/S800 experimental setup.
////////////                          
////////////          Current Maintainers:
////////////                 Kathrin Wimmer  (wimmer@phys.s.u-tokyo.ac.jp)
////////////                 Eric Lunderberg (lunderberg@nscl.msu.edu)
////////////
////////////          Distribution:
////////////                   Please do not redistribute this software directly.
////////////                   If someone new wants a copy of this software,
////////////                 email one of the maintainers for the download link.
////////////                   This allows us to keep track of who has the software,
////////////                 and properly distribute updates and bug fixes.
////////////                 
////////////          Suggestions:
////////////                   We view the development of the software as a collaborative
////////////                 effort, and as such, welcome and appreciate any suggestions
////////////                 for bug fixes and improvements.
////////////
////////////          Disclaimer:
////////////                 This software is provided as-is, with no warranty.
////////////                 No current or future support is guaranteed for this software.
////////////
////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
#include "Fit.hh"
short gError;
void gaussian(double x, double a[], double *y, double dyda[], int na)
{
  // Gaussian with offset
  double arg, ex, fac;
  arg = (x - a[2]) / a[3];
  ex = exp(-arg * arg);
  fac = a[1] * ex * 2 * arg;
  *y = a[1] * ex + a[4];
  dyda[1] = ex;
  dyda[2] = fac / a[3];
  dyda[3] = fac * arg / a[3];
  dyda[4] = 1;
}

int Fit(double *x, double *y, double *sig, double *a, int MA, int Npoints, double *chisq)
{
  int i, *ia, itst, k, v1=MA;
  double alamda, ochisq, **covar, **alpha;

  ia = f6(1, MA);
  covar = f8(1, MA, 1, MA);
  alpha = f8(1, MA, 1, MA);

  for (i=1; i<=v1; i++) ia[i]=1;
  gError = 0;

  // initialize fit
  alamda = -1;
  f1(x, y, sig, Npoints, a, ia, MA, covar, alpha, chisq, gaussian, &alamda);
  k = 1;
  itst = 0;

  // fitting loop
  for (;;) {
    k++;
    ochisq = *chisq;
    f1(x, y, sig, Npoints, a, ia, MA, covar, alpha, chisq, gaussian, &alamda);
    // if we have encountered some kind of error during the fit
    // we forget about this event and keep going
    if (gError == 1) {
      ff8(alpha, 1, MA, 1, MA);
      ff8(covar, 1, MA, 1, MA);
      ff6(ia, 1, MA);
      return 0;
    }
    if (*chisq > ochisq)
      itst = 0;
    else if (fabs(ochisq - *chisq) < 0.1)
      itst++;
    if (itst < 4) continue;
    alamda = 0;
    f1(x, y, sig, Npoints, a, ia, MA, covar, alpha, chisq, gaussian, &alamda);
    break;
  }
  // success!
  ff8(alpha, 1, MA, 1, MA);
  ff8(covar, 1, MA, 1, MA);
  ff6(ia, 1, MA);
  return 1;
}

void f1(double x[], double y[], double sig[], int ndata, double a[], int ia[],
	int ma, double **covar, double **alpha, double *chisq,
	void (*funcs)(double, double [], double *, double [], int), double *alamda)
{
	void f4(double **covar, int ma, int ia[], int v1);
	void f3(double **a, int n, double **b, int m);
	void f2(double x[], double y[], double sig[], int ndata, double a[],
		int ia[], int ma, double **alpha, double v5[], double *chisq,
		void (*funcs)(double, double [], double *, double [], int));
	int j,k,l,m;
	static int v1;
	static double v3,*v4,*v5,*v6,**v2;

	if (*alamda < 0.0) {
		v4=f7(1,ma);
		v5=f7(1,ma);
		v6=f7(1,ma);
		for (v1=0,j=1;j<=ma;j++)
			if (ia[j]) v1++;
		v2=f9(1,v1,1,1);
		*alamda=0.001;
		f2(x,y,sig,ndata,a,ia,ma,alpha,v5,chisq,funcs);
		v3=(*chisq);
		for (j=1;j<=ma;j++) v4[j]=a[j];
	}
	for (j=0,l=1;l<=ma;l++) {
		if (ia[l]) {
			for (j++,k=0,m=1;m<=ma;m++) {
				if (ia[m]) {
					k++;
					covar[j][k]=alpha[j][k];
				}
			}
			covar[j][j]=alpha[j][j]*(1.0+(*alamda));
			v2[j][1]=v5[j];
		}
	}
	f3(covar,v1,v2,1);
	for (j=1;j<=v1;j++) v6[j]=v2[j][1];
	if (*alamda == 0.0) {
		f4(covar,ma,ia,v1);
		ff9(v2,1,v1,1,1);
		ff7(v6,1,ma);
		ff7(v5,1,ma);
		ff7(v4,1,ma);
		return;
	}
	for (j=0,l=1;l<=ma;l++)
		if (ia[l]) v4[l]=a[l]+v6[++j];
	f2(x,y,sig,ndata,v4,ia,ma,covar,v6,chisq,funcs);
	if (*chisq < v3) {
		*alamda *= 0.1;
		v3=(*chisq);
		for (j=0,l=1;l<=ma;l++) {
			if (ia[l]) {
				for (j++,k=0,m=1;m<=ma;m++) {
					if (ia[m]) {
						k++;
						alpha[j][k]=covar[j][k];
					}
				}
				v5[j]=v6[j];
				a[l]=v4[l];
			}
		}
	} else {
		*alamda *= 10.0;
		*chisq=v3;
	}
}

void f2(double x[], double y[], double sig[], int ndata, double a[], int ia[],
	int ma, double **alpha, double beta[], double *chisq,
	void (*funcs)(double, double [], double *, double [], int))
{
	int i,j,k,l,m,mfit=0;
	double ymod,wt,sig2i,dy,*dyda;

	dyda=f7(1,ma);
	for (j=1;j<=ma;j++)
		if (ia[j]) mfit++;
	for (j=1;j<=mfit;j++) {
		for (k=1;k<=j;k++) alpha[j][k]=0.0;
		beta[j]=0.0;
	}
	*chisq=0.0;
	for (i=1;i<=ndata;i++) {
		(*funcs)(x[i],a,&ymod,dyda,ma);
		sig2i=1.0/(sig[i]*sig[i]);
		dy=y[i]-ymod;
		for (j=0,l=1;l<=ma;l++) {
			if (ia[l]) {
				wt=dyda[l]*sig2i;
				for (j++,k=0,m=1;m<=l;m++)
					if (ia[m]) alpha[j][++k] += wt*dyda[m];
				beta[j] += dy*wt;
			}
		}
		*chisq += dy*dy*sig2i;
	}
	for (j=2;j<=mfit;j++)
		for (k=1;k<j;k++) alpha[k][j]=alpha[j][k];
	ff7(dyda,1,ma);
}

#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

void f3(double **a, int n, double **b, int m)
{
	int *indxc,*indxr,*ipiv;
	int i,j,k,l,ll;
	double big,dum,pivinv,temp;
	int icol=-1, irow = -1; //Prevent compiler from complaining about lack of initialization.

	indxc=f6(1,n);
	indxr=f6(1,n);
	ipiv=f6(1,n);
	for (j=1;j<=n;j++) ipiv[j]=0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if (ipiv[j] != 1)
				for (k=1;k<=n;k++) {
					if (ipiv[k] == 0) {
						if (fabs(a[j][k]) >= big) {
							big=fabs(a[j][k]);
							irow=j;
							icol=k;
						}
					} else if (ipiv[k] > 1) f5("Error: singular matrix");
				}
		++(ipiv[icol]);
		if (irow != icol) {
			for (l=1;l<=n;l++) SWAP(a[irow][l],a[icol][l])
			for (l=1;l<=m;l++) SWAP(b[irow][l],b[icol][l])
		}
		indxr[i]=irow;
		indxc[i]=icol;
		if (a[icol][icol] == 0.0) f5("Error: singular matrix");
		pivinv=1.0/a[icol][icol];
		a[icol][icol]=1.0;
		for (l=1;l<=n;l++) a[icol][l] *= pivinv;
		for (l=1;l<=m;l++) b[icol][l] *= pivinv;
		for (ll=1;ll<=n;ll++)
			if (ll != icol) {
				dum=a[ll][icol];
				a[ll][icol]=0.0;
				for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
				for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
			}
	}
	for (l=n;l>=1;l--) {
		if (indxr[l] != indxc[l])
			for (k=1;k<=n;k++)
				SWAP(a[k][indxr[l]],a[k][indxc[l]]);
	}
	ff6(ipiv,1,n);
	ff6(indxr,1,n);
	ff6(indxc,1,n);
}
#undef SWAP

#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}

void f4(double **covar, int ma, int ia[], int mfit)
{
	int i,j,k;
	double swap;

	for (i=mfit+1;i<=ma;i++)
		for (j=1;j<=i;j++) covar[i][j]=covar[j][i]=0.0;
	k=mfit;
	for (j=ma;j>=1;j--) {
		if (ia[j]) {
			for (i=1;i<=ma;i++) SWAP(covar[i][k],covar[i][j])
			for (i=1;i<=ma;i++) SWAP(covar[k][i],covar[j][i])
			k--;
		}
	}
}
#undef SWAP


#define NR_END 1
#define FREE_ARG char*

void f5(const char* error_text)
/* Numerical Recipes standard error handler */
{
//	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	//	fprintf(stderr,"...now exiting to system...\n");
	//	exit(1);
	gError = 1;
}

int *f6(long nl, long nh)
/* allocate an int Vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) f5("allocation failure!");
	return v-nl+NR_END;
}

double *f7(long nl, long nh)
/* allocate a double Vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) f5("allocation failure!");
	return v-nl+NR_END;
}

double **f8(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) f5("allocation failure!");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) f5("allocation failure!");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

double **f9(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) f5("allocation failure!");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) f5("allocation failure!");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

void ff6(int *v, long nl, long nh)
/* free an int Vector allocated with iVector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void ff7(double *v, long nl, long nh)
/* free a double Vector allocated with dVector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void ff8(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by matrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void ff9(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}
