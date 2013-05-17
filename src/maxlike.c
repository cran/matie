
#include <R.h>
#include <Rmath.h>


double kerFunc(double x, double kw){
    double kerVal = exp(-R_pow_di(x,2)/(2.0*R_pow_di(kw,2)));
	return(kerVal);
}

/* likelihoodVec is the meat of the MATIE package.
 * It computes a vector of likelioods for the joint model, one for each data point.
 * it is also used to comute likelihoods for the null model by computing a
 * series of joint likelihoods (one for each dependency grouping).
 *
 * The following allows a significant speedup in likelihood computations
 * The routine assumes that the rankData is sorted on the first variable and uses this
 * assumption to compute the liklihood of each data point in rankData.
 * The trick is to progress down rankData, keeping left and right pointers at the
 * bounds of the kernel for the first variable.  All data points in this strip
 * are then considered for accumulating likelihoods.
 */

void likelihoodVec(int *nIn, int *mIn, double *rankData, 
		double *kvec, double *loob, double *liVec )
{

	int n = nIn[0];
    int m = mIn[0];
    double lb = loob[0];
	
	double kw = kvec[0];
	
	double temp, sum;

	int k, i, j, jj;
	int kr, kn;
	kr = (int)(3.0 * kw + 1.0); // kernel radius
	kn = 2 * kr +1; // kernel width

    /* compute the 1D kernel */
	double* ker;
	ker = (double *) R_alloc(kn, sizeof(double));
	sum=0.0;
	for (i=-kr; i<kr+1; i++){
		temp = kerFunc((double)i,kw);
		sum += temp;
		ker[i+kr] = temp;
	}
    
	/* normalize the 1D kernel */ 
    for (i=0; i<kn; i++){
        ker[i] /= sum;
    }
    
	/* compute the height of the 1D kernel at it's centre */
	double kh = ker[kr];
    
    /* compute accumulated 1D kernel values 
    double* accKer;
	accKer = (double *) R_alloc(kn, sizeof(double));
	sum=0.0;
	for (i=-kr; i<kr+1; i++){
		sum += ker[i+kr];
		accKer[i+kr] = sum;
	} */
    
    // Rprintf(" accKer[%i] = %f",kr,accKer[kr]);
    
    /* compute kernel scale factors for each data point 
    double* scaleFactor;
    scaleFactor = (double *) R_alloc(n, sizeof(double));
    for (i=0; i<n; i++){
        scaleFactor[i]=1.0;
        for (j=0; j<m; j++){
            int p = rankData[i+j*n]-1;
            if(p<kr)scaleFactor[i] /= accKer[kr+p]; 
            if(p>n-1-kr)scaleFactor[i] /= accKer[kr+n-(p+1)];
        }
        // Rprintf(" scaleFactor[%i] = %f",i,scaleFactor[i]);
    } */

	/* assuming the rankData is sorted on the first coord */
	int pl=0;
	int pm=0;
	int pr=0;
	while (pm < n) {
        liVec[pm] = 0.0;
		int xpm = rankData[pm]-1;
		int xpl = rankData[pl]-1;
		int xpr = rankData[pr]-1;
		while(xpm-xpl>kr && pl<pm){
			pl++;
			xpl = rankData[pl]-1;
		}
    /* bug fix: changing && pr<n to && pr<(n-1) */
		while(xpr-xpm<=kr && pr<(n-1)){
			pr++;
			xpr = rankData[pr]-1;
		}
		for (j=pl; j<pr; j++){
            temp = 1.0;
            for(jj = 0; jj < m; jj++){
                int jyp = rankData[j+jj*n]-1;
                int ypm = rankData[pm+jj*n]-1;
                int yp = jyp-ypm+kr;
                if(yp>=0 && yp<kn){
                     temp = temp * ker[yp];
                } else {
                    temp = 0.0;
                }
            }
            liVec[pm] += temp;  // *scaleFactor[pm];
		}
		pm++;
	}

	/* compute likelihood on a leave one out basis
	   at each data point subract the kernel height and 
	   scale by total distribution - one kernel contribution */
	for (k=0; k<n; k++){
		liVec[k] = (liVec[k]-lb*R_pow_di(kh,m) /* *scaleFactor[k] */ )/((double)(n)-lb);
	}
		
}

/* The next routine computes a distribution on a 2D grid 
 * using a weighted sum of marginal and full kernel heights
 * The marginal kernel width, full kernel width and weight
 * are passed as parameters.
 *
 * rankData is assumed to be of dimensions n x 2.
 *
 * The routine makes no assumptions on rankData order
 * The distribution is computed on an n x n grid 
 */

void getDistribution(int *nIn, int *rankData, 
                   double *kvec, double *weight, double *dist)
{
    
	int n = nIn[0];
    
	double mkw = kvec[0];
    double kw = kvec[1];
    
    double w = weight[0];
	
	double temp, sum;
    
	int k, i, j;

    /* set up the 1D kernel for computing marginals */
	int mkr, mkn;
	mkr = (int)(3.0 * mkw + 1.0); // marginal kernel radius
	mkn = 2 * mkr + 1; // marginal kernel length

	/* compute the 1D kernel for estimating marginals */
	double* sker;
	sker = (double *) R_alloc(mkn, sizeof(double));
	sum=0.0;
	for (i=-mkr; i<mkr+1; i++){
		temp = kerFunc((double)i,mkw);
		sum += temp;
		sker[i+mkr] = temp;
	}
    
	/* normalize the 1D kernel */ 
    for (i=0; i<mkn; i++){
        sker[i] /= sum;
    }
    
    /* compute the marginals by adding 1D kernels at each interior point */
    double *margDist;
    margDist = (double *) malloc(n * sizeof(double));
	for (k=0; k<n; k++){
		margDist[k] = 0.0;
	}

	for (k=0; k<n; k++){
		for (i=-mkr; i<mkr+1; i++){
			temp = sker[(i+mkr)];
			if( (k+i)>=0 && (k+i)<n)margDist[k+i] += temp;
		}
	}
    
    
    /* now do the computations for the full model */
    int kr, kn;
	kr = (int)(3.0 * kw + 1.0); // full kernel radius
	kn = 2 * kr +1; // full kernel length
	
    /* compute a 1D kernel for building the full distribution */
	double* ker;
	ker = (double *) R_alloc(kn, sizeof(double));
	sum=0.0;
	for (i=-kr; i<kr+1; i++){
		temp = kerFunc((double)i,kw);
		sum += temp;
		ker[i+kr] = temp;
	}
    
	/* normalize the 1D kernel */ 
    for (i=0; i<kn; i++){
        ker[i] /= sum;
    }
    
    /* compute the full distribution by adding kernel products at each data point */
    double *fullDist;
    fullDist = (double *) malloc(n * n * sizeof(double));

	int n2 = n * n;
	for (k=0; k<n2; k++){
        fullDist[k] = 0.0;
	}
    
	for (k=0; k<n; k++){
		int xp = rankData[k]-1;
		int yp = rankData[k+n]-1;
		for (i=-kr; i<kr+1; i++){
			for(j=-kr; j<kr+1; j++){
				// temp = ker[(j+kr)*kn+(i+kr)];
                temp = ker[i+kr] * ker[j+kr];
				int ip = xp+i;
				int jp = yp+j;
				if(ip>=0 && ip<n && jp>=0 && jp<n)fullDist[jp*n+ip] += temp;
			}
		}
        
	}
    
    
	/* return the weighted distribution in dist */
	for (i=0; i<n; i++){
		for(j=0; j<n; j++){
			dist[j*n+i] = (w * margDist[i] * margDist[j]) + 
                (1.0-w) * fullDist[j*n+i] * (double) n;	
		}
	}
    
    /* free up the array space */
    free(margDist);
    free(fullDist);
}

