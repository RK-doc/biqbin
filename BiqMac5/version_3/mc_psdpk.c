/* 
 * solver for 
 * (P) max <L,X>              (D) min e'y
 *     s.t. diag(X) = e           s.t. Diag(y) - L=Z
 *          X psd                      Z psd
 * 
 * (predictor/corrector)
 *
 * input:  L, n
 * output: X, y
 * 
 * translation from mc_psdpk.m (f. rendl, 4/03)
 * dec 2004, a.wiegele
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "declarations.h"

#define MAX_ITER 30

int mc_psdpk(L,n,silent,pX,py,pphi)
    /* input */
    double *L;
    int n;
    int silent;
    /* output */
    double **pX; /* primal matrix */
    double **py; /* dual variables */
    double *pphi; /* dual function value of the sdp */
{
    /* vars for call to blas- and lapack-routines */
    int inc=1;
    double alpha;
    char up='U';
    int info,*ipiv;
    //double *work;
    double *ap,*ap2;

    /* and some other vars */
    int i,j,k,nn,nsq;
    double *Z; /* dual matrix */
    double phi,psi;
    double delta;
    double mu;
    double *Zi,*dy,*dZ,*dX,*Xnew;
    double *dy1, *dX1;
    double alphap=1.0,alphad=1.0;
    double *tmp, *tmp1;  /* need non-symm mat when computing Zi*diag(dy)*X  */
    int index1,index2;
    double *p,*p2;

    /**********************************
     *                                *
     * initialization:                *
     *                                *
     * starting triplet               *
     * primal/dual cost, gap          *
     **********************************/

    nn=n*(n+1)/2;
    nsq=n*n;

    /* allocate space for Z */
    Z=(double *)malloc(nn*sizeof(double));
    if (Z == NULL)
        {
            printf("Storage allocation failed!\n");
            exit(10);
        };

    /* X = eye(n); y=ones(n,1); */ /*    Diag(pX,*py,n); */
    for (i=1,p=(*pX),p2=(*py)+1;i<=n;i++,p2++)
      {         
		(*p2)=1.0;
		for (j=1;j<i;j++,p++)
		  (*p)=0.0;
		(*p)=1.0;
		p++;
      };


    /* y = sum( abs(L) )' + 1.; */
    for (i=1,p=(*py)+1;i<=n;i++,p++)
      {
		//(*py)[i]+=Abs(L[ijtokp(j,i,n)]);

		for (j=1,p2=L+i*(i-1)/2;j<=i;j++,p2++)
		  (*p)+=Abs(*p2);
		for (j=i+1,p2+=(i-1);j<=n;j++)
		  {
			(*p)+=Abs(*p2);
			p2+=j;
		  };
      };

    /* Z = diag(y) - L; */
    Diag(&Z,*py,n);                     /* Z=diag(y) */
    alpha=-1.0;
    daxpy_(&nn,&alpha,L,&inc,Z,&inc); /* Z=Z-L */        

    /* phi = ones(n,1)'*y; psi = L(:)'*X(:); */
    /* psi=traceprod(L,*pX,n); */    /* remember that X=I */
    phi=0.0;
    psi=0.0;
    for (i=1,p=(*py)+1,p2=L;i<=n;i++,p++)
        {
		  phi+=*p;  //phi+=(*py)[i];
		  psi+=*p2; //psi+=L[ijtokp(i,i,n)];
		  p2+=(i+1);
        };

    delta=phi-psi;

    /* mu = Z(:)'*X(:) / (2*n); */
    mu=traceprod(Z,*pX,n)/(2*n);

    /*****************************
     *                           *
     * main loop                 *
     *                           *
     *****************************/

    /* 
     * allocate space for dX, dy (which will also be used as rhs), 
     * dZ, Xnew, Zi=inv(Z), ap, tmp 
     */
    dX=(double *)malloc(nn*sizeof(double));
    dy=(double *)malloc((n+1)*sizeof(double));
    dZ=(double *)malloc(nn*sizeof(double));
    if ((dy == NULL)||(dX == NULL)||(dZ == NULL))
        {
            printf("Storage allocation failed!\n");
            exit(10);
        };
    dy1=(double *)malloc((n+1)*sizeof(double));
    dX1=(double *)malloc(nn*sizeof(double));
    Xnew=(double *)malloc(nn*sizeof(double));
    Zi=(double *)malloc(nn*sizeof(double));
    ap=(double *)malloc(nn*sizeof(double));
    ap2=(double *)malloc(nn*sizeof(double));
    if ((Zi == NULL)||(dy1 == NULL)||(dX1 == NULL)||(Xnew == NULL)||(ap == NULL)||(ap2 == NULL))
        {
            printf("Storage allocation failed!\n");
            exit(10);
        };

    tmp=(double *)malloc(nsq*sizeof(double));  
    tmp1=(double *)malloc(nsq*sizeof(double));  /* general matrices! */
    if ((tmp == NULL)||(tmp1 == NULL))
        {
            printf("Storage allocation failed!\n");
            exit(10);
        };

    /* 
     * allocate space for int array ipiv which is used in 
     * the call to dspsv. 
     * also, allocate work space for dsptri.
     */
    ipiv=(int *)malloc((n+1)*sizeof(int));
    //work=(double *)malloc(n*sizeof(double));
    if ((ipiv == NULL))//||(work == NULL))
        {
            printf("Storage allocation failed!\n");
            exit(10);
        };

    if (silent == 0)
        {
            printf("\n iter alphap alphad  log(gap)    lower     upper \n");
            printf("                    ");
            printf(" %7.5f  %8.3f  %8.3f\n",log10(delta),psi,phi); 
        };

    for (i=1;(delta > EPSILON6)&&(delta > Abs(phi)*EPSILON6);i++) /* while duality gap too large */
        { 
            /*
             * compute Zi=inv(Z);
             */
            dcopy_(&nn,Z,&inc,Zi,&inc);        /* copy Z to Zi */
            //dsptrf_(&up,&n,Zi,ipiv+1,&info);   /* computes factorization used for dsptri */
			dpptrf_(&up,&n,Zi,&info);          /* computes factoraization used for dpptri */
            if (info > 0)
                {
                    fprintf(stderr,"mc_psd: Z is singular! info = %d\n",info);
				//fwrite(L,sizeof(double),n*(n+1)/2,stderr);
                    exit(10);
                };

            //dsptri_(&up,&n,Zi,ipiv+1,work,&info);  /* compute Zi=inv(Z) */
			dpptri_(&up,&n,Zi,&info);              /* compute Zi=inv(Z) */

            if (info > 0)
                {
                    fprintf(stderr,"mc_psd: Z is singular! info = %d\n",info);
				//fwrite(L,sizeof(double),n*(n+1)/2,stderr);
                    exit(10);
                };

            /* 
             * predictor step 
             * solves:  Z * X + diag(dy1) * X + Z * dX1 = 0 
             *
             */

            /* dy1 = (Zi .* X) \ (-e) */
            for (j=1,p=dy1+1;j<=n;j++,p++)
				*p=-1.0; //dy1[j]=-1.0;
            /* ap = Zi .* X */
            for (j=0,p=ap;j<nn;j++,p++)
				*p=Zi[j]*(*pX)[j]; //ap[j]=Zi[j]*(*pX)[j];
            /*  dy1 = ap \ dy1; */
            dspsv_(&up,&n,&inc,ap,ipiv+1,dy1+1,&n,&info); /* NOTE: ap is changed on exit! */
            if (info != 0)
                {
                    fprintf(stderr,"mc_psdpk (predictor step): ");
                    fprintf(stderr,"problems in solving a system of linear equations! %d\n",info);
		    //fwrite(L,sizeof(double),n*(n+1)/2,stderr);
                    exit(10);
                };
            /* dX1 = -Zi*diag(dy1)*X - X */
            for (j=1;j<=n;j++)            /* tmp = -Zi*diag(dy1) */
                for (k=1;k<=n;k++)        
                    if (k <= j)
                        tmp[ijtok(j,k,n)]=-dy1[k]*Zi[ijtokp(k,j,n)];
                    else
                        tmp[ijtok(j,k,n)]=-dy1[k]*Zi[ijtokp(j,k,n)];

            /* multiply general matrix and symm packed matrix */
            my_dgsmm(&n,tmp,*pX);              /* tmp = tmp*X */

            /* symmetrize tmp and store it in dX1 */
            for (k=1;k<=n;k++)
                for (j=1;j<=k;j++)
                    dX1[ijtokp(j,k,n)]=(tmp[ijtok(j,k,n)]+tmp[ijtok(k,j,n)])/2.0;

            alpha=-1.0;
            daxpy_(&nn,&alpha,*pX,&inc,dX1,&inc); /* dX1 = dX1 - X */


           /*
             * corrector step 
             * solves: diag(dy2)*X + Z*dX2-muI + diag(dy1)*dX1 = 0
             *
             */

            /* rhs = mu*diag(Zi) - (Zi .* dX1)*dy1; */
            diag(Zi,&dy,n);            /* dy = diag(Zi) */
            for (j=0;j<nn;j++)         /* ap2 = Zi .* dX1 */
                ap2[j]=Zi[j]*dX1[j];
            /* dy = -(ap2)*dy1 + mu*dy */
            alpha=-1.0;
            dspmv_(&up,&n,&alpha,ap2,dy1+1,&inc,&mu,dy+1,&inc);

            /*  dy = (Zi .* X) \ (mu*diag(Zi) - (Zi .* dX1)*dy1); */
            /* ap = Zi .* X */
            for (j=0,p=ap;j<nn;j++,p++)
			*p=Zi[j]*(*pX)[j]; //ap[j]=Zi[j]*(*pX)[j];
            dspsv_(&up,&n,&inc,ap,ipiv+1,dy+1,&n,&info);
            if (info != 0)
                {
                    printf("mc_psd (corrector step): ");
                    printf("problems in solving a system of linear equations! %d\n",info);
                    exit(10);
                };


            /* 
             * dX2 = mu*Zi - Zi*diag(dy2)*X - Zi*diag(dy1)*dX1); 
             */

            for (j=1;j<=n;j++)            /* tmp = -Zi*diag(dy) */
                for (k=1;k<=n;k++)        /* tmp1 = -Zi*diag(dy1) */
				  {
					index1=ijtok(j,k,n);
					index2=( k<=j ? ijtokp(k,j,n) : ijtokp(j,k,n));
					tmp[index1]=-dy[k]*Zi[index2];
					tmp1[index1]=-dy1[k]*Zi[index2];
				  };


            /* multiply general matrix and symm packed matrix */
            my_dgsmm(&n,tmp,*pX);              /* tmp = tmp*X */
            my_dgsmm(&n,tmp1,dX1);             /* tmp1 = tmp1*dX1 */
            alpha=1.0;                         /* tmp = tmp+tmp1 */
            daxpy_(&nsq,&alpha,tmp1,&inc,tmp,&inc);
            /* symmetrize  tmp -> dX */
            for (k=1;k<=n;k++)
                for (j=1;j<=k;j++)
                    dX[ijtokp(j,k,n)]=(tmp[ijtok(j,k,n)]+tmp[ijtok(k,j,n)])/2.0;

            daxpy_(&nn,&mu,Zi,&inc,dX,&inc);  /* dX = mu*Zi + dX  */


            /* final steps */
            alpha=1.0;
            daxpy_(&n,&alpha,dy1+1,&inc,dy+1,&inc);     /* dy=dy1+dy */
            daxpy_(&nn,&alpha,dX1,&inc,dX,&inc);        /* dX=dX1+dX */


            /* 
             * find steplength alphap and alphad
             */
            for (alphap=1.0,info=1;info > 0;)
                {
                    dcopy_(&nn,*pX,&inc,Xnew,&inc);        /* Xnew=X */
                    daxpy_(&nn,&alphap,dX,&inc,Xnew,&inc); /* Xnew=alphap*dX+Xnew */
                    dpptrf_(&up,&n,Xnew,&info);
                    if (info > 0)
                        alphap*=0.8;
                }
            if (alphap < 1)                       /* stay away from boundary */
                alphap*=0.95;
            daxpy_(&nn,&alphap,dX,&inc,*pX,&inc); /* X=alphap*dX+X */

            Diag(&dZ,dy,n);                       /* dZ=diag(dy) */
            for (alphad=1.0,info=1;info > 0;)
                {
                    dcopy_(&nn,Z,&inc,Xnew,&inc);          /* Xnew=Z */
                    daxpy_(&nn,&alphad,dZ,&inc,Xnew,&inc); /* Xnew=alphad*dZ+Xnew */
                    dpptrf_(&up,&n,Xnew,&info);
                    if (info > 0)
                        alphad*=0.8;
                }
            if (alphad < 1)
                alphad*=0.95;

            /*
             * update
             */
            daxpy_(&n,&alphad,dy+1,&inc,*py+1,&inc);  /* y=y+alphad*dy */
            daxpy_(&nn,&alphad,dZ,&inc,Z,&inc);       /* Z=Z+alphad*dZ */


            /* mu = Z(:)'*X(:) / (2*n); */
            mu=traceprod(Z,*pX,n)/(2*n);

            if (alphap + alphad > 1.6)
                mu*=0.5;
            if (alphap + alphad > 1.9)
                mu/=5;  /* reduce mu if stepsize good */

            /* phi = ones(n,1)'*y; psi = L(:)'*X(:); */
            phi=0.0;
            for (j=1,p=(*py)+1;j<=n;j++,p++)
				phi+=(*p); // phi+=(*py)[j];
            psi=traceprod(L,*pX,n);  
            delta=phi-psi;


            if (silent == 0)
                printf("   %2d %6.3f %6.3f  %7.5f  %8.3f  %8.3f \n",i,alphap,alphad,log10(delta),psi,phi); 

        } 

    /* free and return */
    *pphi=phi;
    //free(work);
    free(ipiv);
    free(tmp1);free(tmp);
    free(ap2);free(ap);free(Zi);free(Xnew);free(dX1);free(dy1);
    free(dZ);free(dy);free(dX);
    free(Z);
    return(0);
} 


void my_dgsmm(pn,pA,pB)
     int *pn;
     double *pA,*pB;
{
  int i,j,nsq;
  double *C,*B;
  double dmy;
  double *p;
  /* vars for call to blas- and lapack-routines */
  double alpha,beta;
  char up='U',ri='R';
  
  B=(double *)malloc((*pn)*(*pn)*sizeof(double));
  C=(double *)calloc((*pn)*(*pn),sizeof(double));
  if ((B == NULL)||(C == NULL))
    {
      printf("Storage allocation failed!\n");
      exit(10);
    };
  
                       
  for (i=1,p=pB;i<=*pn;i++)
    for (j=1;j<=i;j++,p++)
      {
	dmy=*p; //pB[ijtokp(j,i,(*pn))];
	B[ijtok(i,j,(*pn))]=dmy;
	B[ijtok(j,i,(*pn))]=dmy;
      };


  /* C = (*pA)*B + 0.0*C */

  nsq=(*pn)*(*pn);
  alpha=1.0; beta=0.0;
  dsymm_(&ri,&up,pn,pn,&alpha,B,pn,pA,pn,&beta,C,pn);    
  
  j=1;
  dcopy_(&nsq,C,&j,pA,&j);
  free(C);free(B);
}
