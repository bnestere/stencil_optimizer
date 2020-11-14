#include <omp.h>
#include <math.h>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <stencil_config.h>
#if defined(_OPENMP)
#endif
double seconds() {
   struct timeval tv;
   gettimeofday(&tv,NULL);
   return (double)(tv.tv_sec)+1e-6*tv.tv_usec;
}
int validate_results(int x_max,int y_max,int z_max,int T_MAX,void** _opt_res) {
   double tim1;double tim2;double nFlops;
   int i;int j;int k;int t;
   int x;int y;int z;
   double alpha;double beta;
   double(*u_0_0)[x_max] = (double*)malloc(2*x_max*sizeof (double));
   alpha = 1.f/(double)x_max;
   beta = 2.f/(double)y_max;
   for (i=0; i<x_max; i+=1)
     {
        u_0_0[0][i] = 1.+i*0.1;
     }
   tim1 = seconds();
   for (t=0; t<T_MAX; t+=1)
     {
        for (x=1; x<x_max-1; x+=1)
          {
             u_0_0[1-t%2][x] = u_0_0[t%2][x]+0.125*(u_0_0[t%2][x+1]-2*u_0_0[t%2][x]+u_0_0[t%2][x-1]);
          }
     }
   tim2 = seconds();
   double(*opt_res)[x_max] = _opt_res;
   for (x=1; x<x_max-1; x+=1)
     {
        if (u_0_0[1][x]!=opt_res[1][x]) 
          {
             fprintf(stderr,"Index %d differs by %.16lf",x,u_0_0[1][x]-opt_res[1][x]);
             exit(EXIT_FAILURE);
          }
     }
   printf("%s\n%s\n%s\n","---------------------------------------------------------","Opt stencil successfully validated against normal version","---------------------------------------------------------");
   nFlops = (double)((x_max-2)*T_MAX*5.0);
   printf("ORIGINAL FLOPs in stencil code:      %e\n",nFlops);
   printf("ORIGINAL Time spent in stencil code: %f\n",tim2-tim1);
   printf("ORIGINAL Performance in GFlop/s:     %f\n",nFlops/(1e9*(tim2-tim1)));
   free(u_0_0);
   return 0;
}
int main(int argc,char** argv) {
   int x_max;int y_max;int z_max;
   int i;int j;int k;int t;
   int x;int y;int z;
   int T_MAX;
   double tim1;double tim2;double nFlops;
   double alpha;double beta;
   if (argc!=5) 
     {
        printf("Wrong number of parameters, <xmax> <ymax> <zmax> <timesteps>.\n",argv[0]);
        exit(-1);
     }
   x_max = atoi(argv[1]);
   y_max = atoi(argv[2]);
   z_max = atoi(argv[3]);
   T_MAX = atoi(argv[4]);
   double(*u_0_0)[x_max] = (double*)malloc(2*x_max*sizeof (double));
   alpha = 1.f/(double)x_max;
   beta = 2.f/(double)y_max;
   for (i=0; i<x_max; i+=1)
     {
        u_0_0[0][i] = 1.+i*0.1;
     }
   tim1 = seconds();
  int t1, t2, t3, t4;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
if ((T_MAX >= 1) && (x_max >= 3)) {
  for (t1=-1;t1<=floord(T_MAX-1,16);t1++) {
    lbp=max(ceild(t1,2),ceild(32*t1-T_MAX+2,32));
    ubp=min(floord(16*t1+x_max+13,32),floord(x_max+T_MAX-3,32));
#pragma omp parallel for private(lbv,ubv,t3,t4)
    for (t2=lbp;t2<=ubp;t2++) {
      for (t3=max(max(max(0,16*t1),32*t1-32*t2+1),32*t2-x_max+2);t3<=min(min(min(T_MAX-1,16*t1+31),32*t2+30),32*t1-32*t2+x_max+29);t3++) {
        lbv=max(max(32*t2,t3+1),-32*t1+32*t2+2*t3-31);
        ubv=min(min(32*t2+31,-32*t1+32*t2+2*t3),t3+x_max-2);
#pragma ivdep
#pragma vector always
        for (t4=lbv;t4<=ubv;t4++) {
          u_0_0[-(t3 % 2) + 1][(-t3+t4)] = (u_0_0[t3 % 2][(-t3+t4)] + (0.125 * ((u_0_0[t3 % 2][(-t3+t4) + 1] - (2 * u_0_0[t3 % 2][(-t3+t4)])) + u_0_0[t3 % 2][(-t3+t4) - 1])));;
        }
      }
    }
  }
}
   tim2 = seconds();
   validate_results(x_max,y_max,z_max,T_MAX,u_0_0);
   nFlops = (double)((x_max-2)*T_MAX*5.0);
   printf("\nOPTIMIZED FLOPs in stencil code:      %e\n",nFlops);
   printf("OPTIMIZED Time spent in stencil code: %f\n",tim2-tim1);
   printf("OPTIMIZED Performance in GFlop/s:     %f\n",nFlops/(1e9*(tim2-tim1)));
   free(u_0_0);
   return EXIT_SUCCESS;
}
