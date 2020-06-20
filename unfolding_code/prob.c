#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <string.h>
#include <malloc.h>
#include "lib.h"


int main (int argc, char **argv){

  int i,j,k,kk,n,jj,atoms,holes,r,so,DN,bin;
  double eval[100000],gamma,xmin,xmax,h,L,ge;
  FILE *lg,*lh;
  char nv[100];
  int atm[]={5,9};
  int hol[]={9,4};

  r=2; so=5; bin=50; ge=1.0;

  for(kk=1;kk<2;kk++){
  atoms=(int)atm[kk]; holes=(int)hol[kk];
  sprintf(nv,"../Add_Paper/energies_s%dNL%d%d_g%0.3lfFr%d.dat",so,atoms,holes,ge,r); 
  lg=fopen(nv,"r");
  i=0; do{ fscanf(lg,"%lf",&eval[i]);  i++;}while(!feof(lg)); fclose(lg);
  int dim=i-1;
  double evals[dim],auxevals[dim],mev,s=0.0,mls[dim-1],ms,err,emax,err0=1;;
  stupidsort(eval,dim);

  // unfolding :D
  for(i=0;i<dim;i++)evals[i]=eval[i]-eval[0];
  for(i=0;i<dim-1;i++)mls[i]=(evals[i+1]-evals[i])/evals[dim-1];
  if(dim<1200){
    ms=mtwert(mls,dim-1); n=(int)(0.7/ms); bin=90;
  }
  else{
    n=1401; bin=120;
  }
  k=0;
  gamma=0.0009; // initial gamma!
  do{
    for(i=0;i<dim;i++)auxevals[i]=evals[i]/evals[dim-1];
    printf("(n = %d, %02d) gamma = %lf\t",n,k+1,gamma);
    ms=0;err=0;
    // STUPID ERROR BECAUSE OF NOT CHANGING EVALS 
    Loren_unfold(gamma,auxevals,dim,n,&ms,&err,&emax);
    //easyunfold(gamma,auxevals,dim,&ms,&err,&emax); 
    printf("ms = %lf, err = %0.11lf\n",ms,err);
    
    if(fabs(ms-1.0)<4.0e-3 && err<1.0e-4)
      break;
    else if(err>err0 && ms>0.98){
      printf("err>err0 ...!");
      break;
    }
    else
      gamma=gamma/1.05;
    err0=err;
    k++;
  }while(k<400);
  printf("k-value = %d\n",k);
  
  for(i=0;i<dim;i++)evals[i]=emax*auxevals[i];
  sprintf(nv,"datafile/ufenergies_so%d_NL%d%d_g%0.3lfFr%d.dat",so,atoms,holes,ge,r); lg=fopen(nv,"w");
  for(i=0;i<dim;i++)
    fprintf(lg,"%0.22lf\n",evals[i]/evals[dim-1]);
  fclose(lg);
  
  double ds[dim-1];
  
  for(i=0;i<dim-1;i++)ds[i]=(evals[i+1]-evals[i])/ms;
  
  sprintf(nv,"datafile/histogram_so%d_NL%d%d_g%0.3lfFr%d.dat",so,atoms,holes,ge,r);
  hist(ds,dim-1,0.0,10.0,bin,nv); // Histogram
  
  // Plotting... 
  lg=fopen("h.gpl","w");
  fprintf(lg,"set term postscript enhanced color\n");
  //fprintf(lg,"set term png\n");
  fprintf(lg,"set output 'datafile/Ps_ge%0.3lf_so%d_r%d_NL%d%d.eps'\n",ge,so,r,atoms,holes);
  fprintf(lg,"set grid\n");
  fprintf(lg,"set key sample 3\n");
  fprintf(lg,"set xrange [0:5]\n");
  fprintf(lg,"set yrange [0:1.05]\n");
  fprintf(lg,"POISSON(x)=exp(-x);\n");
  fprintf(lg,"GOE(x)=(%0.22lf/2.0)*x*exp(-(%0.22lf/4.0)*x**2)\n",M_PI,M_PI);
  fprintf(lg,"set rmargin 10\n");
  fprintf(lg,"set lmargin 10\n");
  fprintf(lg,"set bmargin 5\n");
  fprintf(lg,"set xlabel 's' font 'Helvetica,39'\n");
  fprintf(lg,"set ylabel 'P(s)' font 'Helvetica,30'\n");
  fprintf(lg,"set xtics font 'Helvetica,18'\n");
  fprintf(lg,"set ytics font 'Helvetica,18'\n");
  fprintf(lg,"plot\tGOE(x) title 'GOE  ' w l lw 3 lt 2 lc rgb '#000000',POISSON(x) title 'POISSON' w l lw 3 lt 2 lc rgb '#ff0000','datafile/histogram_so%d_NL%d%d_g%0.3lfFr%d.dat' title 'Numeric' w l lt 1 lw 2 lc rgb '#000000'\n",so,atoms,holes,ge,r);
  fclose(lg);
  system("gnuplot h.gpl");
  system("rm h.gpl");
  
  // Number variance
  FILE *nl;
  sprintf(nv,"datafile/Numvars_so%d_NL%d%d_g%0.3lfFr%d.dat",so,atoms,holes,ge,r);
  nl=fopen(nv,"w");
  for(L=0.2;L<=25;L+=0.25){
    //    ll=(int)((rint)L);
    DN=(rint)(evals[dim-1]/L);
    double ne[DN];
    for(i=0;i<DN;i++)ne[i]=0.0;
    
    for(i=0;i<DN;i++){
      xmin=i*L; xmax=(i+1)*L;
      /*
      h=(xmax-xmin)/(1.0*n);
      for(j=0;j<n;j++)xi[j]=j*h+xmin;
      Loren_dens(xi,yi,gamma,evals,dim,n);
      simp(n,h,yi,&s);
      ne[i]=pow(s-L,2);
      */
      j=0;
      for(k=0;k<dim;k++){
	if(evals[k]<=xmax && evals[k]>xmin){
	  j++;
	}
      }
      ne[i]=pow(1.0*j-L,2);
    }
    fprintf(nl,"%lf\t%lf\n",L,mtwert(ne,DN));
    //    printf("%lf\t%lf\n",L,mtwert(ne,DN));
  }
  fclose(nl);

  // Plotting... 
  lg=fopen("h.gpl","w");
  fprintf(lg,"set term postscript enhanced color\n");
  //fprintf(lg,"set term png\n");
  fprintf(lg,"set output 'datafile/Nv_ge%0.3lf_so%d_r%d_NL%d%d.eps'\n",ge,so,r,atoms,holes);
  fprintf(lg,"set grid\n");
  fprintf(lg,"set key sample 3\n");
  fprintf(lg,"set xrange [0:25]\n");
  fprintf(lg,"set yrange [0:5]\n");
  fprintf(lg,"POISSON(x)=x;\n");
  fprintf(lg,"GOE(x)=(2.0/%0.22lf**2)*(log(2.0*%0.22lf*x)+0.57722+1-%0.22lf**2/8.0)\n",M_PI,M_PI,M_PI);
  fprintf(lg,"set rmargin 10\n");
  fprintf(lg,"set lmargin 10\n");
  fprintf(lg,"set bmargin 5\n");
  fprintf(lg,"set xlabel 'L' font 'Helvetica,39'\n");
  fprintf(lg,"set ylabel '{/Symbol S}^2(L)' font 'Helvetica,30'\n");
  fprintf(lg,"set xtics font 'Helvetica,18'\n");
  fprintf(lg,"set ytics font 'Helvetica,18'\n");
  fprintf(lg,"plot\tGOE(x) title 'GOE  ' w l lw 3 lt 2 lc rgb '#000000',POISSON(x) title 'POISSON' w l lw 3 lt 2 lc rgb '#ff0000','datafile/Numvars_so%d_NL%d%d_g%0.3lfFr%d.dat' title 'Numeric' w p ps 1.5 pt 6 lc rgb '#000000'\n",so,atoms,holes,ge,r);
  fclose(lg);
  system("gnuplot h.gpl");
  
  }  

  return 0;
}
