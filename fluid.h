
#define N 200

#define DRAND(x) (1.0*x*rand()/(double)RAND_MAX)

// (2*Pi)^(3/2)
#define M_2PI3_2 ((double)(15.749610295))

#define IX(i,j) ((i)+(N+2)*(j))
#define SIZE ((N+2)*(N+2))

#define SWAP(x0,x) {double* tmp=(x0);(x0)=(x);(x)=tmp;}

#define DQ(dt) (dt)

#define LIM(x,h) (x)=((x)>(h))?(h):(x)
#define BILIN_INT(x,k,l,kf,lf) \
((x)[IX((k),(l))]*(1-(kf))*(1-(lf))\
+(x)[IX((k)+1,(l))]*(kf)*(1-(lf))\
+(x)[IX((k),(l)+1)]*(1-(kf))*(lf)\
+(x)[IX((k)+1,(l)+1)]*(kf)*(lf))

#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))

#define ABS(a) ((a)>0?(a):(-(a)))

#define  DT(i,j) (0.0001)// (tm[IX(i,j)])

#define  DIFFUSION_D  0.0001
#define  DIFFUSION_DT  2.0
#define  DIFFUSION_V  0.1

#define  DIFFUSION_DM  1.0
#define  ADVECT_DM  0.9

#define  GF  100.0

#define  TU  0.1
#define  TE  0
#define  CF  1000.0

#define  RE 0.02

#define  VORT  0.1
#define  WMAX  50.0

#define  FFLUID  1.0
#define  FVISC  0.5

double u[SIZE], u_p[SIZE];
double u_src[SIZE];

double v[SIZE], v_p[SIZE];
double v_src[SIZE];

double dens[SIZE], dens_p[SIZE];
double dens_src[SIZE];
double dens_T[SIZE], dens_Tp[SIZE];
double dens_Tsrc[SIZE];

double w[SIZE];
double w_p[SIZE];

double tm[SIZE];
double tm_p[SIZE];
double tm_src[SIZE];

char bnd[SIZE];

//--------------------------------

void FLUID_init();
void FLUID_update();

//-------

void FLUID_friction(double* x,
                      double ff, double fv);

void FLUID_random_step(double* u, double* v, double* d);

void FLUID_vorticity_step(double* w, double* wp,
                          double* u, double* v,
                          double vorticity, double wmax);

void FLUID_density_step(double* x, double* x0, double* xsrc,
                          double* u, double* v,
                          double diff, double advf);

void FLUID_velocity_step( double* u, double* v, 
                          double* u0, double* v0,
                          double visc);

void FLUID_diffuse(int b, double* x, double* x0,
                          double diff);

void FLUID_advect(int b, double* d, double* d0,
                         double* u, double* v,
                         double f);

void FLUID_gravity_step( double* u, double* v, double* dens,
                          double g);

void FLUID_temp_step( double* u, double* v,
                        double* dens, double* dens_T,
                        double fu, double fe, double cf);

void FLUID_cool_step( double* dens, double cf);

void FLUID_project( double* u, double* v, double* p, double* div);

void FLUID_addSource(double* x, double* s);

void FLUID_set_bnd( int b, double* x);
