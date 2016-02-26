#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "fluid.h"

//----------------------------

#include "loutron.h"

extern gimage gimp_image;

//----------------------------

void _debug_isnan(double* d)
{
  int i;
  for(i=0;i<SIZE;i++)
    if(isnan(d[i]))
    {
      fprintf(stderr,"NAN %d %3.3f\n",i,d[i]);
      exit(42);
    }

  return;
}

void FLUID_init()
{
  memset(u,0,SIZE*sizeof(double));
  memset(u_p,0,SIZE*sizeof(double));
  memset(u_src,0,SIZE*sizeof(double));

  memset(v,0,SIZE*sizeof(double));
  memset(v_p,0,SIZE*sizeof(double));
  memset(v_src,0,SIZE*sizeof(double));

  memset(dens,0,SIZE*sizeof(double));
  memset(dens_p,0,SIZE*sizeof(double));

  memset(dens_T,0,SIZE*sizeof(double));
  memset(dens_Tp,0,SIZE*sizeof(double));

  memset(dens_src,0,SIZE*sizeof(double));
  memset(dens_Tsrc,0,SIZE*sizeof(double));
  
  memset(bnd,0,SIZE*sizeof(char));

  memset(w,0,SIZE*sizeof(double));
  memset(w_p,0,SIZE*sizeof(double));

  memset(tm,0,SIZE*sizeof(double));
  memset(tm_p,0,SIZE*sizeof(double));
  memset(tm_src,0,SIZE*sizeof(double));

	int i,j;
  double d;
  
  for(i=0;i<=N;i++)
  for(j=0;j<=N;j++)
  {
    tm[IX(i,j)] = +0.0001;
  }

  double R;
  int px,py;

  for(i=0;i<=N;i++)
  for(j=0.9*N;j<=N;j++)
  {
    dens[IX(i,j)] = 300;
    dens_T[IX(i,j)] = 900;
  }

  px = 0.5*N;
  py = 0.85*N;

  R=0.05*N;
  for(j=-R; j<R; j++)
  for(i=-R; i<R; i++)
  {
    dens_T[IX(i+px,j+py)] = 20000;
    dens[IX(i+px,j+py)] = 2000;
    //double _t = tm[IX(i+px,j+py)] - DRAND(0.00001);
    //tm[IX(i+px,j+py)] = (_t>0)?_t:0.000001;
  }
/*
  px = 0.5*N;
  py = 0.5*N;

  R=0.1*N;
  for(j=-R; j<R; j++)
  for(i=-R; i<R; i++)
  {
    bnd[IX(i+px,j+py)] = 1;
  }
*/  
}

void FLUID_update()
{
  static int t=0; t++;
/*
  if(i>(0.01/DT(i,j))) memset(dens_src,0,SIZE*sizeof(double));
  if(i>(0.01/DT(i,j))) memset(dens_Tsrc,0,SIZE*sizeof(double));
  if(i>(0.01/DT(i,j))) memset(u_src,0,SIZE*sizeof(double));
  if(i>(0.01/DT(i,j))) memset(v_src,0,SIZE*sizeof(double));
*/ 
 /*
  if(1)
  {
    int d,c,x1,x2,ds,dst,vs;
    int i,j;

    d = rand()%(N/8);
    c = rand()%N;
    ds = DRAND(600.0);
    dst = DRAND(1000.0);
    vs = DRAND(1000.0) - 500.0;

    x1 = c-(d/2);
    x2 = c+(d/2);
    
    if(x1<0) x1=0;
    if(x2>N) x2=N;
    
    for(i=x1;i<x2;i++)
    for(j=N/2-5;j<N/2+5;j++)
    {
      v[IX(i,j)] = vs;
      dens[IX(i,j)] = ds;
      dens_T[IX(i,j)] = dst;
    }
  }
*/

  //FLUID_density_step(tm,tm_p, tm_src, u, v, DIFFUSION_DM, ADVECT_DM);
  FLUID_density_step(dens,dens_p, dens_src, u, v, DIFFUSION_D, 1.0);
  FLUID_density_step(dens_T,dens_Tp, dens_Tsrc, u, v, DIFFUSION_DT, 1.0);
  FLUID_gravity_step(u,v,dens,GF);

  FLUID_temp_step(u,v,dens,dens_T,TU,TE,CF);
  //FLUID_cool_step(dens_T, CF);

  //FLUID_random_step(u,v,dens_T);
  
  FLUID_vorticity_step(w,w_p,u,v,VORT,WMAX);
  
  FLUID_velocity_step(u,v,u_p,v_p, DIFFUSION_V);

  FLUID_friction(u,FFLUID,FVISC);
  FLUID_friction(v,FFLUID,FVISC);

}

void FLUID_friction(double* x,
                      double ff, double fv)
{
  int i,j;
  double v,h;
  for(j=1;j<=N;j++)
  for(i=1;i<=N;i++)
  {
    if(x[IX(i,j)] != 0.0)
    {
      v = fabs(x[IX(i,j)]);
      h = (ff*x[IX(i,j)] + fv*x[IX(i,j)]*v)*DT(i,j);
      x[IX(i,j)] -= h;
    }
  }

  return;
}

void FLUID_vorticity_step(double* w, double* w_p,
                          double* u, double* v,
                          double vorticity, double wmax)
{
  int i,j;
  double Nx,Ny,Nm;
  double lw=0;
  
  //_______________________________
  // Compute vorticity grid
  for(j=1;j<=N;j++)
  for(i=1;i<=N;i++)
  {
    lw = ( (v[IX(i+1,j)] - v[IX(i-1,j)])
                  -(u[IX(i,j+1)] - u[IX(i,j-1)]));
    
    if(lw > wmax) lw = wmax;
    if(lw < -wmax) lw = -wmax;

 
    w[IX(i,j)] = lw;
  }
    
  //________________________________
  // Compute velocity grid
  for(j=1;j<=N;j++)
  for(i=1;i<=N;i++)
  {
    Nx = w[IX(i+1,j)] - w[IX(i-1,j)];
    Ny = w[IX(i,j+1)] - w[IX(i,j-1)];
    
      
    Nm = sqrt(Nx*Nx+Ny*Ny);

    if(Nm == 0)
      continue;

    Nx /= Nm;
    Ny /= Nm;

    lw = vorticity
      *(0.01*0.01*N*N)*fabs(w[IX(i,j)])
      *DT(i,j)*10000
      *DT(i,j)*10000;

    u[IX(i,j)] += Ny*lw;
    v[IX(i,j)] += -Nx*lw;
  }
}

void FLUID_random_step(double* u, double* v, double* d)
{
  int i,j;
  int dm;
  for(j=1;j<=N;j++)
  for(i=1;i<=N;i++)
  {
    dm = (int)(d[IX(i,j)]+1);
    u[IX(i,j)] +=  RE*((rand()%dm) - dm/2);
    v[IX(i,j)] +=  RE*((rand()%dm) - dm/2);
  }
}

void FLUID_cool_step( double* denst, double cf)
{
  int i,j;
  for(j=1;j<=N;j++)
  for(i=1;i<=N;i++)
  {
    denst[IX(i,j)] -= cf*DQ(DT(i,j));
    if(denst[IX(i,j)] <= 0) denst[IX(i,j)] = 0;
  }
}

void FLUID_temp_step( double* u, double* v,
                        double* dens, double* dens_T,
                        double fu, double fe, double cf)
{
  double f;
  double Nx,Ny,Nm;
  int i,j;

  for(j=1;j<=N;j++)
  for(i=1;i<=N;i++)
  {
    f = DT(i,j)*fe;

    Nx = (dens_T[IX(i+1,j)] - dens_T[IX(i-1,j)]);
    Ny = (dens_T[IX(i,j+1)] - dens_T[IX(i,j-1)]);

    u[IX(i,j)] += -1.0*Nx*f;
    v[IX(i,j)] += -1.0*Ny*f;

    v[IX(i,j)] += -1.0*dens[IX(i,j)]*dens_T[IX(i,j)]*DT(i,j)*fu;
  }
}

void FLUID_gravity_step( double* u, double* v, double* dens,
                          double g)
{
  int i,j;
  for(j=1;j<=N;j++)
  for(i=1;i<=N;i++)
  {
    v[IX(i,j)] += 1.0*dens[IX(i,j)]*DT(i,j)*g;
  }

}

void FLUID_velocity_step( double* u, double* v, 
                double* u0, double* v0,
                double visc)
{
  FLUID_addSource( u, u_src );
  FLUID_addSource( v, v_src );
 
  SWAP( u0, u ); FLUID_diffuse( 1, u, u0, visc);
  SWAP( v0, v ); FLUID_diffuse( 2, v, v0, visc);

  FLUID_project( u, v, u0, v0 );

  SWAP( u0, u ); SWAP( v0, v );

  FLUID_advect( 1, u, u0, u0, v0, 1.0);
  FLUID_advect( 2, v, v0, u0, v0, 1.0);
  FLUID_project( u, v, u0, v0 );
}


void FLUID_density_step(double* x, double* x0, 
                        double* src,
                        double* u, double* v,
                        double diff,double advf)
{

  FLUID_addSource(x, src);

  SWAP(x0,x); FLUID_diffuse( 0, x, x0, diff);
  SWAP(x0,x); FLUID_advect( 0, x, x0, u, v, advf);
}

void FLUID_diffuse ( int b, double * x, double * x0, double diff )
{
  int i, j, k;
  double a;

  for ( k=0 ; k<20 ; k++ ) {
   for ( j=1 ; j<=N ; j++ ) 
     for ( i=1 ; i<=N ; i++ ) {

        a=DT(i,j)*diff*N*N;

        x[IX(i,j)] = (x0[IX(i,j)] 
                        + a*(x[IX(i-1,j)]+x[IX(i+1,j)]
                        +    x[IX(i,j-1)]+x[IX(i,j+1)]))/(1+4*a);
      }
    FLUID_set_bnd ( b, x );
  }
}

void FLUID_advect ( int b, 
              double * d, double * d0,
              double * u, double * v,
              double f)
{
  int i=0, j=0, i0=0, j0=0, i1=0, j1=0;
  double x=0, y=0, s0=0, t0=0, s1=0, t1=0, dt0=0;
  for ( j=1 ; j<=N ; j++ ) {
    for ( i=1 ; i<=N ; i++ ) {

      dt0 = DT(i,j)*N*f;
      x = i-dt0*u[IX(i,j)]; y = j-dt0*v[IX(i,j)];
      if (x<0.5) x=0.5; if (x>N+0.5) x=N+ 0.5; i0=floorf(x); i1=i0+1;
      if (y<0.5) y=0.5; if (y>N+0.5) y=N+ 0.5; j0=floorf(y); j1=j0+1;
      s1 = x-i0; s0 = 1-s1; t1 = y-j0; t0 = 1-t1;
      d[IX(i,j)] = s0*(t0*d0[IX(i0,j0)]
                  +t1*d0[IX(i0,j1)])
                  +s1*(t0*d0[IX(i1,j0)]
                  +t1*d0[IX(i1,j1)]);
    }
  }
  FLUID_set_bnd ( b, d );
}


void FLUID_addSource(double* x, double* s)
{
  int i,j;

  for(j=1;j<=N;j++)
  for(i=1;i<=N;i++)
    x[IX(i,j)] += DQ(DT(i,j))*s[IX(i,j)];

}

void FLUID_project( double* u, double* v, double* p, double* div)
{
  int i, j, k;
  int b=0;
  double h;
  h = 1.0/N;

  for ( i=1 ; i<=N ; i++ ) {
    for ( j=1 ; j<=N ; j++ ) {
      div[IX(i,j)] = -0.5*h*(u[IX(i+1,j)]-u[IX(i-1,j)]+
                             v[IX(i,j+1)]-v[IX(i,j-1)]);
      p[IX(i,j)] = 0;
    }
  }

  FLUID_set_bnd( 0, div ); FLUID_set_bnd( 0, p );
  for ( k=0 ; k<20 ; k++ ) {
    for ( i=1 ; i<=N ; i++ ) {
      for ( j=1 ; j<=N ; j++ ) {
        p[IX(i,j)] = (div[IX( i,j)]+p[IX(i-1,j)]+p[IX(i+1,j)]+
                      p[IX(i,j-1)]+p[IX(i,j+1)])/4;
      }
    }
    FLUID_set_bnd( 0, p );
  }
  for ( i=1 ; i<=N ; i++ ) {
    for ( j=1 ; j<=N ; j++ ) {
      u[IX(i,j)] -= 0.5*(p[IX(i+1,j)]-p[IX(i-1,j)])/h;
      v[IX(i,j)] -= 0.5*(p[IX(i,j+1)]-p[IX(i,j-1)])/h;
    }
  }
  FLUID_set_bnd( 1, u ); FLUID_set_bnd( 2, v );
}

void FLUID_set_bnd ( int b, double * x )
{
  int i,j;
  for ( i=1 ; i<=N ; i++ ) {
    x[IX(0 ,i)]  = (b==1)?(-x[IX(1,i)]):(x[IX(1,i)]);
    x[IX(N+1,i)] = (b==1)?(-x[IX(N,i)]):(x[IX(N,i)]);
    x[IX(i,0 )]  = (b==2)?(-x[IX(i,1)]):(x[IX(i,1)]);
    x[IX(i,N+1)] = (b==2)?(-x[IX(i,N)]):(x[IX(i,N)]);
  }

  x[IX(0 ,0 )] =   0.5*(x[IX(1,0 )]+x[IX(0 ,1)]);
  x[IX(0 ,N+1)] =  0.5*(x[IX(1,N+1)]+x[IX(0 ,N )]);
  x[IX(N+1,0 )] =  0.5*(x[IX(N,0 )]+x[IX(N+1,1)]);
  x[IX(N+1,N+1)] = 0.5*(x[IX(N,N+1)]+x[IX(N+1,N )]);

  for(j=1;j<=N;j++)
  for(i=1;i<=N;i++)
  {
    if(bnd[IX(i,j)])
    {
      if( bnd[IX(i-1,j)] ) x[IX(i-1,j)] = -x[IX(i-1,j)];
      if( bnd[IX(i+1,j)] ) x[IX(i+1,j)] = -x[IX(i+1,j)];
      if( bnd[IX(i,j-1)] ) x[IX(i,j-1)] = -x[IX(i,j-1)];
      if( bnd[IX(i,j+1)] ) x[IX(i,j+1)] = -x[IX(i,j+1)];
    }
  }
}
