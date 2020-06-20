#include "lib.h"

double mtwert(double *data,int length){
  double y;
  int i;

  y=0.0;
  for(i=0;i<length;i++){
    y=y+data[i];
  }
  return y/(1.0*length);
}

void stupidsort(double *a, int length){
  int i, j;
  double temp;

  for(i=1; i<length; i++){
    for(j=0; j<i; j++){
      if(a[j] >= a[i]){
	temp = a[i];
	a[i] = a[j];
	a[j] = temp;
      }
    }
  }

}

void simp(int n,double h,double f[],double *s){
  int i;
  double s0,s1,s2;
  
  s0 = 0;
  s1 = 0;
  s2 = 0;
  for (i = 1; i < n-1; i = i+2)
    {
      s0 = s0+f[i];
      s1 = s1+f[i-1];
      s2 = s2+f[i+1];
    }
  *s = h*(s1+4*s0+s2)/3;
  
  /* If n is even, add the last slice separately */
  
  if (n%2 == 0) *s = *s+h*(5*f[n-1]+8*f[n-2]-f[n-3])/12;
}

void hist(double *array,int size,double xi,double xf,int bin,char *file){
  double delta,a[bin+1],b[bin],x;
  int i,n,k,m;
  FILE *lg;

  delta=(xf-xi)/(1.0*bin);
  for(i=0;i<bin+1;i++)a[i]=i*delta+xi;
  
  m=0;
  for(i=0;i<bin;i++){
    n=0;
    for(k=0;k<size;k++){
      if(array[k]<=a[i+1] && array[k]>=a[i])n++;
    }
    b[i]=1.0*n;
    m+=n;
  }
  
  for(i=0;i<bin;i++)b[i]=b[i]/(delta*m);
  
  //  lg=fopen("datafile/histogram.dat","w");
  lg=fopen(file,"w");
  for(i=0;i<bin;i++){
    for(x=a[i];x<=a[i+1];x+=delta/(50.0)){
      fprintf(lg,"%0.17lf\t%0.17lf\n",x,b[i]);
    }
  }
  fclose(lg);
}


void Loren_dens(double *x,double *y,double gamma,double *eval,int dim,int size){
   
  double y0;
  int i,j,k;

  stupidsort(eval,dim); // sorting
  
  for(i=0;i<size;i++){
    y0=0.0;
    for(j=0;j<dim;j++){
      y0+=(gamma/M_PI)*(1.0/(pow(x[i]-eval[j],2)+pow(gamma,2)));
    }
    y[i]=y0;
  }
}


void Loren_unfold(double gamma,double *eval,int dim,int size,double *mls,double *err,double *emax){
   
  double s,y0,evals[dim],xmin,xmax,ms,h,x[size],y[size],dx[dim];
  int i,j,k,p;

  for(i=0;i<size;i++){x[i]=0.0;y[i]=0.0;}
  for(i=0;i<dim;i++)
    evals[i]=eval[i]-eval[0];

  // unfolding
  for(i=0;i<dim;i++)eval[i]=0;
  for(p=0;p<dim;p++){
    s=0;
    xmin=1.0*evals[0];  xmax=1.0*evals[p];
    h=(xmax-xmin)/(1.0*size);
    
    for(i=0;i<size;i++){
      x[i]=i*h-xmin;
    }
    
    for(i=0;i<size;i++){
      y0=0.0;
      for(j=0;j<dim;j++){
	y0+=(gamma/M_PI)*(1.0/(pow(x[i]-evals[j],2)+pow(gamma,2)));
      }
      y[i]=y0;
    }
    simp(size,h,y,&s);
    eval[p]=s;
  }
  ms=0;
  for(i=0;i<dim-1;i++)
    ms+=eval[i+1]-eval[i];

  *mls=ms/(1.0*dim); // mean level spacing
  *emax=eval[dim-1];
  for(i=0;i<dim;i++){
    eval[i]=eval[i]/eval[dim-1];
    dx[i]=pow(eval[i]-i/(1.0*dim),2);
  }
  *err = mtwert(dx,dim);
}

void easyunfold(double gamma,double *eval,int dim,double *mls,double *err,double *emax){
  int i,j;
  double auxeval[dim],dx[dim-1],y,ms;
  
  for(i=0;i<dim;i++){
    y=0;
    for(j=0;j<dim;j++){
      y+=atan((eval[i]-eval[j])/gamma)-atan(-eval[j]/gamma);
    }
    auxeval[i]=y/M_PI;
  }
  ms=0;
  for(i=0;i<dim-1;i++)
    ms+=auxeval[i+1]-auxeval[i];

  *mls=ms/(1.0*dim); // mean level spacing
  *emax=auxeval[dim-1];
  for(i=0;i<dim;i++){
    auxeval[i]=auxeval[i]/auxeval[dim-1];
    dx[i]=pow(auxeval[i]-i/(1.0*dim),2);
  }
  *err = mtwert(dx,dim-1);
  for(i=0;i<dim;i++)eval[i]=auxeval[i];
}

