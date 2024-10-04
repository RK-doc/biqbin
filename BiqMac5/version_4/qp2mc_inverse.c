/*
  Transform back the Biq Mac output of a Max-Cut into a 
  binary quadratic problem solution. 
  min y'*Q*y s.t. y in {0,1}^n

  Return 0 if ok, 1 if failure.

  Read through the lines until a line start with 'one side of the cut'.
  The next line states one side of the shore.

  NOTE: the program works for graphs with size at most 1000! 
  (if you need larger sizes, adjust M.)

  a. wiegele, iasi, sept 7, 2007

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "declarations.h"

#define M 1000

int main(int argc,char *argv[])
{
  int ret,i,m,n;
  int timelimit=0;
  int first,c;
  int cut;
  double cutvalue;


  /*
   * Read until you find the line 'The input graph has...'
   */
  for(first=getchar();first != EOF && first != 'T';first=getchar()) {
    putchar(first);
    if (first != '\n') {
      do {
	c=getchar();putchar(c);
      }while (c!='\n');
    };
  }
    
  if (first == EOF) {
    printf("qp2mc: Problem with converting Biq Mac output!\n");
    printf("Please send an email to biqmac@uni-klu.ac.at\n");
    exit(1);
  }
  if (first == 'T')
    do {
      c=getchar();
    }while (c!='s');
  scanf("%d",&n);

  for (i=1;i<=10;i++) 
    c=getchar();

  scanf("%d",&m);
  printf("The problem has %d variables,",n-1);
  for (i=1;i<=12;i++) c=getchar();
  do {
    c=getchar();putchar(c);
  }while (c!='\n');



  /* 
   * Ignore lines until you find a line 'cut value:' or 'Upper bound:'
   */ 

  for(first=getchar();first != EOF && first != 'c' && first != 'U';first=getchar()) {
    putchar(first);
    do {
      c=getchar();putchar(c);
    }while (c!='\n');
  }
    
  if (first == EOF) {
    printf("qp2mc: Problem with converting Biq Mac output!\n");
    printf("Please send an email to biqmac@uni-klu.ac.at\n");
    exit(1);
  }
  if (first == 'U') { 
    do {
      c=getchar();
    }while (c!=':');
    scanf("%lf",&cutvalue);
    printf("Lower bound: %lf\n",(-1.0)*cutvalue);
    c=getchar(); //should be \n
    timelimit=1;

    for(first=getchar();first != EOF && first != 'c';first=getchar()) {
      putchar(first);
      do {
	c=getchar();putchar(c);
      }while (c!='\n');
    }
    
    if (first == EOF) {
      printf("qp2mc: Problem with converting Biq Mac output!\n");
      printf("Please send an email to biqmac@uni-klu.ac.at\n");
      exit(1);
    }
  };


  if (first == 'c')
    do {
      c=getchar();
    }while (c!=':');
  scanf("%lf",&cutvalue);
  if (timelimit) {
    if (Zero(cutvalue))
      printf("current best found value: %lf\n",cutvalue);
    else
      printf("current best found value: %lf\n",(-1.0)*cutvalue);
  }
  else {
    if (Zero(cutvalue))
      printf("optimal value: %lf\n",cutvalue);
    else
      printf("optimal value: %lf\n",(-1.0)*cutvalue);
  };
  c=getchar(); //should be \n
  


  /* 
   * Ignore lines until you find a line 'one side of the cut:'
   * (note that it is expected that no other line starts with 
   * the character 'o'!!!)
   */ 

  for(first=getchar();first != EOF && first != 'o';first=getchar()) {
    putchar(first);
    do {
      c=getchar();putchar(c);
    }while (c!='\n');
  }
    
  if (first == EOF) {
    printf("qp2mc: Problem with converting Biq Mac output!\n");
    printf("Please send an email to biqmac@uni-klu.ac.at\n");
    exit(1);
  }
  if (first == 'o')
    do {
      c=getchar();
    }while (c!='\n');

  /*
   * Read in the nodes of one side of the shore.
   * If 1 is in the shore, then we have the values being equal to 1, otherwise 0. ????????
   */
  
  c=getchar();
  if ((c=='\n')||(c==EOF))
    printf("trivial solution. x_i=0 for i=1,..,%d\n",n-1);
  else {
    ret=scanf("%d",&cut);
    if (cut==1) {
      printf("variables equal to 0 are those indexed by:\n");
      for(ret=scanf("%d",&cut);(cut != EOF)&&(ret==1);ret=scanf("%d",&cut))
	printf(" %d",cut-1);
    }
    else { 
      if ((cut>1)&&(cut<n)){
	printf("variables equal to 1 are those indexed by:\n");
	do {
	  printf(" %d",cut-1);
	  ret=scanf("%d",&cut);
	} while ((cut != EOF)&&(ret==1));
      }
      else  {
	printf("qp2mc: Problem with converting Biq Mac output!\n");
	printf("Please send an email to biqmac@uni-klu.ac.at\n");
	exit(1);
      };
    };
    printf("\n");
  };

  /*
   * Free memory, close file and return.
   */
  
  return(0);
}

