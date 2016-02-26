#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <SDL/SDL.h>
#include <time.h>

#include "sdl_utils.h"

#include "version.h"

#include "fluid.h"

//#define BILINEAR

//#define SAVETOFILE

int mode=0;
int Nmodes=7;

char updateScreen = 1;

//-------------------------------

void SDL_window_resize(int width,int height);

void SDL_draw_line(SDL_Surface* s, 
                    int x, int y, int xh, int yh, Uint16);

void SDL_draw_pixel(SDL_Surface* s,int x, int y, Uint16);

void SDL_draw_bilinear();

void SDL_draw_blocs();

void SDL_write_surface_toBMP(SDL_Surface* s, char* filename);

//-------------------------------

SDL_Surface* SDL_screen = NULL;

int SDL_screen_x = 500;
int SDL_screen_y = 500;
int SDL_screen_bpp = 32;
Uint32 SDL_screen_flags = SDL_HWSURFACE;

int t = 0;

int main(int argc, char** argv)
{
  srand( time(NULL) );

  printf("\nSulphate version %d.%d.%d %s %s\n", VERSION_H, VERSION_L, VERSION_N, __TIME__,__DATE__);
  
  // Init SDL
  if( SDL_Init(SDL_INIT_VIDEO) < 0 )
  {
    fprintf(stderr, "Unable to initialize SDL: %s\n",SDL_GetError());
    exit(1);
  }
  
  // trap exit()
  atexit(SDL_Quit);
 
  // Init video
  SDL_screen = SDL_SetVideoMode(SDL_screen_x, SDL_screen_y, 
                                SDL_screen_bpp,
                                SDL_screen_flags);
  
  if(!SDL_screen)
  {
    fprintf(stderr, "Unable to set video mode: %s\n",SDL_GetError());
    exit(2);
  }

  FLUID_init();
  
  int fps_c=0;
  int fps_t = time(NULL);
  while(1)
  {
    t++;
    fps_c++;
    if(!(fps_c%5))
    {
      int t=time(NULL);
      fprintf(stderr,"\rmode %2.2d frame #%3.3d %.2f fps  ",
                mode,fps_c,(double)fps_c/(double)(t-fps_t));
    }

    FLUID_update();

#ifdef SAVETOFILE
    static int _frame=0; _frame++;
    char file[1024];
    sprintf(file,"./out/f%6.6d.bmp",_frame);

    SDL_write_surface_toBMP(SDL_screen,file);
#endif

    if(updateScreen)
    {  
      SDL_LockSurface(SDL_screen);
 
#ifdef BILINEAR
      SDL_draw_bilinear();
#else
      SDL_draw_blocs();
#endif 
      SDL_UnlockSurface(SDL_screen);

      SDL_UpdateRect(SDL_screen,0,0,0,0);
    }

    int R,px,py;
    int i,j,k;
  

    static int moux=-1;
    static int mouy=-1;

    SDL_Event event;
    while(SDL_PollEvent(&event))
    {
      switch(event.type)
      {

        case SDL_MOUSEMOTION:

          px = event.motion.x;
          py = event.motion.y;
 
          px = (px*N)/SDL_screen_x;
          py = (py*N)/SDL_screen_y;


          R=0.01*N;
          for(j=-R; j<R; j++)
          for(i=-R; i<R; i++)
          {
            u[IX(i+px,j+py)] += -(50.0+DRAND(25.0))*(moux - px)
																-(25-DRAND(50.0))*(mouy-py);
            v[IX(i+px,j+py)] += -(50.0+DRAND(25.0))*(mouy - py)
																-(25-DRAND(50.0))*(moux-px);
          }
					
          R=0.02*N;
          for(j=-R; j<R; j++)
          for(i=-R; i<R; i++)
          {
            dens_T[IX(i+px,j+py)] += DRAND(300.0);
            dens[IX(i+px,j+py)] += DRAND(100.0);
            //double _t = tm[IX(i+px,j+py)] - DRAND(0.00001);
            //tm[IX(i+px,j+py)] = (_t>0)?_t:0.000001;
          }
					
          moux = px;
          mouy = py;

          break;

        case SDL_MOUSEBUTTONDOWN:

          px = event.button.x;
          py = event.button.y;
          
          px = (px*N)/SDL_screen_x;
          py = (py*N)/SDL_screen_y;

          R=0.05*N;
          for(j=-R+1;j<R;j++)
          {
            i=(int)((double)R*cos(asin((double)j/R)));
            for(k=-i;k<i;k++)
            {
              tm[IX(k+px,j+py)] = 0.00001;
            }
          }

          break;

        case SDL_QUIT:
          exit(0);
          break;
        case SDL_VIDEORESIZE:
          SDL_window_resize(event.resize.w,event.resize.h);
          break;
        case SDL_KEYDOWN:
          switch(event.key.keysym.sym)
          {
            case SDLK_ESCAPE:
              exit(0);
              break;

            case SDLK_z:
              mode++; mode %= Nmodes;
              break;

            case SDLK_a:
              mode--;
              if(mode<0)
                mode=Nmodes-1;
              break;
            
            case SDLK_s:
              SDL_write_surface_toBMP(SDL_screen,"./out.bmp");
              exit(0);
              break;

            case SDLK_q:
              updateScreen = !updateScreen;
              break;

            default:break;
          }

          break;
      }
    }
  }


  return 0;
}

inline Uint32 getColorFromParam(double e, double d, double w,
                                double uf, double vf, double t)
{
  double h,s,v;

  e = e/2000;
  if(e>1.0) e=1.0;
  if(e<-1.0) e=-1.0;

  d = d/505;
  if(d>1.0) d=1.0;
  if(d<-1.0) d=-1.0;

  w = 0.02*w;
  if(w>1.0) w=1.0;
  if(w<-1.0) w=-1.0;

  uf = 0.01*uf;
  if(uf>1.0) uf=1.0;
  if(uf<-1.0) uf=-1.0;

  vf = 0.01*vf;
  if(vf>1.0) vf=1.0;
  if(vf<-1.0) vf=-1.0;

  t = 2000*t;
  if(t>1.0) t=1.0;
  if(t<-1.0) t=-1.0;

  switch(mode)
  {
    case 0:
      h=0.5+d*0.10 + e*0.10;
      s=1.1 - d*1.0;
      v=d*0.8+e*0.2;
      break;
    case 1:
      h=0.5-0.45*w;
      s=1.0;
      v=1.0;
      break;
    case 2:
      h=(atan2f(vf,uf)+M_PI)*0.5/M_PI;
      s=1.0;
      v=sqrt(uf*uf+vf*vf);
      break;
    case 3:
      h=(4.0/6.0)*(1.0-d);
      s=1.0;
      v=1.0;
      break;
    case 4:
      h=(4.0/6.0)*(1.0-e);
      s=1.0;
      v=1.0;
      break;
    case 5:
      h=0.3+d*0.20 + e*0.20;
      s=1.0 - e*1.0;
      v=d*0.5+e*0.5;
      break;

    case 6:
      h=0.5+t*0.8;
      s=1.0-fabs(t);
      v=1.0;
      break;

    default:
      h=0;s=0;v=0;
      break;
  }

  LIM(h,1.0);
  LIM(s,1.0);
  LIM(v,1.0);

  return colorFromHSV(SDL_screen, h, s, v);
}

void SDL_draw_blocs()
{
  int i,j,pi,pj,x,y;

  double Ex = N/(double)SDL_screen_x;
  double Ey = N/(double)SDL_screen_y;

  Uint32 tmp[N*N];
  Uint32 color;
  
  for(j=0;j<N;j++)
  for(i=0;i<N;i++)
  {
      color = getColorFromParam(dens_T[IX(i,j)],
                        dens[IX(i,j)],
                        w[IX(i,j)],
                        u[IX(i,j)], v[IX(i,j)],
                        tm[IX(i,j)]);

    tmp[i+j*N] = color;
  }


  for(j=0;j<SDL_screen_y;j++)
  {
    pj = j*Ey;
    for(i=0;i<SDL_screen_x;i++)
    {
      pi = i*Ex;
      SDL_PIXEL(i,j) = tmp[pi+pj*N];
    }
  }

}


void SDL_draw_bilinear()
{
  int i,j;
  int k,l;
  double kf,lf,d;

  for(i=0;i<SDL_screen_x;i++)
  {
    kf=((i*(N-2))/(double)SDL_screen_x);
    k = kf;
    kf = modf(kf,&d);
    
    for(j=0;j<SDL_screen_y;j++)
    {
      lf=((j*(N-2))/(double)SDL_screen_y);
      l = lf;
      lf = modf(lf,&d);

      if(bnd[IX(k,l)])
      {
        SDL_PIXEL(i,j) =
          SDL_MapRGB(SDL_screen->format, 100, 100, 100);
      }
      else
      {
          double e = BILIN_INT(dens_T,k,l,kf,lf);
          double d = BILIN_INT(dens,k,l,kf,lf);
          double wf = BILIN_INT(w,k,l,kf,lf);
          double uf = BILIN_INT(u,k,l,kf,lf);
          double vf = BILIN_INT(v,k,l,kf,lf);
          double t = BILIN_INT(tm,k,l,kf,lf);

          SDL_PIXEL(i,j) = getColorFromParam(e,d,wf,uf,vf,t);
      }
    }
  }
}

void SDL_window_resize(int width,int height)
{

  SDL_screen_x = width;
  SDL_screen_y = height;

  /* changer le mode video */
  if( (SDL_screen = SDL_SetVideoMode(width, height,
                                     SDL_screen_bpp,
                                     SDL_screen_flags)) == 0)
  {
    fprintf(stderr,"Unable to switch video mode: %s\n",SDL_GetError());
    exit(3);
  }

  return;
}

void SDL_draw_pixel(SDL_Surface* s, int x, int y, Uint16 color)
{
  if(x<0) return;
  if(y<0) return;
  if(x>=(s->w)) return;
  if(y>=(s->h)) return;

  SDL_PIXEL_S(s,x,y) = color;
}

void SDL_draw_line(SDL_Surface* s, int x, int y, int xh, int yh,
                    Uint16 color)
{
  //fprintf(stderr,"%d %d\n",xh,yh);
  int i,j;

  if(xh==0) return;
  if(yh==0) return;

  if(xh<0) { x+=xh; xh=-xh; }
  if(yh<0) { y+=yh; yh=-yh; }
  
  double a = (double)yh/(double)xh;


  for(i=0;i<xh;i++)
  for(j=0;j<yh;j++)
  {
 //   fprintf(stderr,"i*a=%d j=%d\n",(int)(a*i),(int)j);
    if( abs((int)(a*i)-(int)j) <= 1)
    SDL_draw_pixel(s,x+i,y+j,color);
  }
}

void SDL_write_surface_toBMP(SDL_Surface* s, char* filename)
{
  FILE* fout = fopen(filename, "w");

  int padding;

  Uint32 tmp32;
  Uint16 tmp16;

  Uint32 image_size = (s->w)*(s->h)*(SDL_screen_bpp/8);
  // FILE HEADER

  fprintf(fout,"BM");

  //___________
  // File size

  tmp32 = 0x36 + image_size;

  tmp32 += (padding = 4 - (tmp32%4));

  fwrite(&tmp32,1,sizeof(Uint32),fout);
    
  //__________
  // Reserved
  fwrite("\0\0\0\0",4,sizeof(char),fout);

  //______________
  // Image offset
  tmp32 = 0x36;
  fwrite(&tmp32,1,sizeof(Uint32),fout);
  
  // IMAGE HEADER
  tmp32 = 0x28;
  fwrite(&tmp32,1,sizeof(Uint32),fout);
  
  //_____________
  // Image width
  tmp32 = s->w;
  fwrite(&tmp32,1,sizeof(Uint32),fout);
 
  //______________
  // Image height
  tmp32 = s->h;
  fwrite(&tmp32,1,sizeof(Uint32),fout);

  //_______
  // Plans
  tmp16 = 1;
  fwrite(&tmp16,1,sizeof(Uint16),fout);

  //_____________
  // Color depth
  tmp16 = SDL_screen_bpp;
  fwrite(&tmp16,1,sizeof(Uint16),fout);

  //_____________
  // Compression
  tmp32 = 0;
  fwrite(&tmp32,1,sizeof(Uint32),fout);

  //____________
  // Image size
  fwrite(&image_size,1,sizeof(Uint32),fout);

  //_________
  // H pix/m
  tmp32 = 0x0B13;
  fwrite(&tmp32,1,sizeof(Uint32),fout);

  //_________
  // v pix/m
  tmp32 = 0x0B13;
  fwrite(&tmp32,1,sizeof(Uint32),fout);

  //________________
  // Palette colors
  tmp32 = 0;
  fwrite(&tmp32,1,sizeof(Uint32),fout);

  //__________________
  // Palette I colors
  tmp32 = 0;
  fwrite(&tmp32,1,sizeof(Uint32),fout);

  //__________
  // Pixels !
  int i,j;
  for(j=(s->h)-1;j>=0;j--)
  for(i=0;i<(s->w);i++)
  {
    fwrite( &(SDL_PIXEL_S(s,i,j)),
            1, sizeof(Uint32),
            fout);
  }

  for(i=0;i<padding;i++)
    fwrite("\0",1,sizeof(char),fout);

  fclose(fout);

  return;
}
