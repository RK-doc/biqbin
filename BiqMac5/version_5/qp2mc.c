/*
  Read in the problem formulated as
  min y'*Q*y s.t. y in {0,1}^n

  and transfrom it to the graph underlying the max cut formulation
  max x'*L*x s.t. x in {-1,1}^(n+1)
  
  Return 0 if ok, 1 if failure.

  first line is n m (n=dim of Q, m=number of given values)
  followed by m lines of the form: i j value

  you may put some comment lines at the beginning of the file. 
  comment lines have to start with the character # 

  output is the graph in rudy-format.

  a. wiegele, iasi, sept 7, 2007

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "declarations.h"

int main(int argc,char *argv[])
{
  char *fname;
  FILE *fid;
  int i,j;
  char c;
  int indexi;
  int indexj;
  double cost;
  int intcost;
  int ret;
  int n,m,m2;
  int n01,m01;
  double *btilde;
  double offset=0.0;
  int integrality=1;

  /*
   * Open file for reading.
   */
  
  if(argc != 2) { 
    printf("qp2mc: wrong call to qp2mc\n");
    exit(1);
  }
  
  fname=argv[1];
  
  fid=fopen(fname,"r");
  if (fid == (FILE *) NULL) {
    printf("qp2mc: Couldn't open problem file for reading! \n");
    printf("Please check your input file or send an email to biqmac@uni-klu.ac.at\n");
    exit(1);
  };
  

  /*
   * Ignore the comment lines.
   */
 
  c=getc(fid);
  while (c == '#'){  
    c=getc(fid);
    while (c!='\n')
      c=getc(fid);
    c=getc(fid);
  };	
  ungetc(c,fid);
  
  /*
   *  Read in dimension n and number of nonzero entries.
   */
  
  ret=fscanf(fid,"%d %d",&n01,&m01);
  if (ret != 2) {
    printf("qp2mc: Incorrect input file. First line has to give n and m values.\n");
    printf("Please check your input file or send an email to biqmac@uni-klu.ac.at\n");
    fclose(fid);
    exit(1);
  };
  
  if (n01<=0 || m01<=0) {
    printf("qp2mc: Incorrect input file. n and m must be positive.\n");
    printf("Please check your input file or send an email to biqmac@uni-klu.ac.at\n");
    fclose(fid);
    exit(1);
  };
 

  /*
   *  Loop through the entries and calculate btilde.
   */

  btilde=(double *)calloc(n01+1,sizeof(double));
  if (btilde == NULL) {
    printf("qp2mc: Storage allocation failed!\n");
    printf("Please send an email to biqmac@uni-klu.ac.at\n");
    exit(1);
  };


  for (i=0,m=0;i<m01;i++) {

    ret=fscanf(fid,"%d %d %lf",&indexi,&indexj,&cost);

    if (ret != 3)
      {
	printf("qp2mc: Incorrect input file. Can't read values.\n");
	printf("Please check your input file or send an email to biqmac@uni-klu.ac.at\n");
	fclose(fid);
	exit(1);
      };
    
    if (!Zero(cost)) {
      if(indexi < 1 || indexi > n01) {
	printf("qp2mc: ERROR: node %d of edge %d is out of range.\n",indexi,1);
	printf("Please check your input file or send an email to biqmac@uni-klu.ac.at\n");
	exit(1);
      }
      if(indexj < 1 || indexj > n01) {
	printf("qp2mc: ERROR: node %d of edge %d is out of range.\n",indexj,1);
	printf("Please check your input file or send an email to biqmac@uni-klu.ac.at\n");
	exit(1);
      }
      if(integrality && !Integer(cost))
	integrality=0;
      if (indexi != indexj) {
	m++;
	btilde[indexi]+=cost;
	btilde[indexj]+=cost;
	offset+=2*cost;
      }
      else {
	btilde[indexi]+=cost;
	offset+=cost;
      };
    };

  };

  fclose(fid);

  /* Print the offset (i.e., e'Qe) as comment. */

  /*
  if (integrality)
    printf("# %d = offset\n",(int) offset);
  else
    printf("# %lf = offset\n",offset);
  */

  for (i=1,m2=0;i<=n01;i++)
    if (!Zero(btilde[i]))
      m2++;

  n=n01+1;
  printf("%d %d\n",n,m+m2); /* entries + btilde-entries */

  
  /*
   * Read the actual entries.
   */

  fid=fopen(fname,"r");
  if (fid == (FILE *) NULL) {
    printf("qp2mc: Couldn't open problem file for reading! \n");
    printf("Please check your input file or send an email to biqmac@uni-klu.ac.at\n");
    exit(1);
  };


  /*
   * Ignore the comment lines and the line specifying n and m.
   */
 
  c=getc(fid);
  while (c == '#'){  
    c=getc(fid);
    while (c!='\n')
      c=getc(fid);
    c=getc(fid);
  };	
  ungetc(c,fid);
  
  ret=fscanf(fid,"%d %d",&indexi,&indexj);

  /*
   * 
   * now there is twice almost the same code: one for integer weights,
   * one for rational weigths.  
   *
   */

  if (integrality) { /* all weights are integer */

    for (i=1;i<=n01;i++) {
      if (!Zero(btilde[i]))
	printf("%d %d %d\n",1,i+1,-((int) btilde[i]));
    };
    
    /* read and write edges */
    
    ret=fscanf(fid,"%d %d %d",&indexi,&indexj,&intcost);
    do {
      if  ((indexi != indexj) && (intcost != 0)) {
	printf("%d %d %d\n",indexi+1,indexj+1,intcost);
      };
      ret=fscanf(fid,"%d %d %d",&indexi,&indexj,&intcost);
    } while (ret == 3);
  
  }
  else { /* there exist rational weights */

    for (i=1;i<=n01;i++) {
      if (!Zero(btilde[i]))
	printf("%d %d %lf\n",1,i+1,-btilde[i]);
    };
    
    /* read and write edges */
    
    ret=fscanf(fid,"%d %d %lf",&indexi,&indexj,&cost);
    do {
      if  ((indexi != indexj) && (!Zero(cost))) {
	printf("%d %d %lf\n",indexi+1,indexj+1,cost);
      };
      ret=fscanf(fid,"%d %d %lf",&indexi,&indexj,&cost);
    } while (ret == 3);
    };
  
  /*
   * Free memory, close file and return.
   */
  
  free(btilde);
  fclose(fid);
  return(0);
}

