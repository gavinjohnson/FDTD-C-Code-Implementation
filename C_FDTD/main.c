//#include "params.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

// read in the parameters for the computation

// Global Domain Size
float Dx = 1;
float Dy = 1;
float Dz = 1;

// array indexing macros
#define Ex(i,j,k) (Ex[(i) + NX*(j) + NX*NY*(k)])
#define Ey(i,j,k) (Ey[(i) + NX*(j) + NX*NY*(k)])
#define Ez(i,j,k) (Ez[(i) + NX*(j) + NX*NY*(k)])
#define Hx(i,j,k) (Hx[(i) + NX*(j) + NX*NY*(k)])
#define Hy(i,j,k) (Hy[(i) + NX*(j) + NX*NY*(k)])
#define Hz(i,j,k) (Hz[(i) + NX*(j) + NX*NY*(k)])

// Number of nodes in each direction
#ifndef NX
#define NX 6
#endif
#ifndef NY
#define NY 6
#endif
#ifndef NZ
#define NZ 6
#endif

// Courant Number
float CFLN = 0.99;

// Simulation time
float simtime = 4e-7;

// sample location informatoin struct
typedef struct sampleLoc{
    char type[1];
    char dir[1];
    int i;
    int j;
    int k;
}sample;

// source information struct
typedef struct source{
    char type[1];
    char dir[1];
    int i1;
    int i2;
    int j1;
    int j2;
    int k1;
    int k2;
}source;




// source function
void sig(float *t, float *val)
{
    static float tw,to;
    tw = 1e-9;
    to = 5*tw;
    static float mag = 1;
    *val = -(mag*exp(-pow(((*t) - to),2)/pow(tw,2))*(2.0*(*t) - 2.0*to))/tw;
}




void Eupdate(float * Ex, float * Ey, float * Ez, float * Hx, float * Hy, float * Hz,
             float * cExy, float * cExz, float * cEyz, float * cEyx, float * cEzx, float * cEzy){
    int i,j,k;
    const float c_cExy = *cExy;
    const float c_cExz = *cExy;
    const float c_cEyz = *cEyz;
    const float c_cEyx = *cEyx;
    const float c_cEzx = *cEzx;
    const float c_cEzy = *cEzy;
    // Ex update
    // constant_cExy = *    <----------  set to const for speed, *pointer
    //    % Ex-update  (excludes boundaries)
    //    for k = 2:nz-1
    //        for j = 2:ny-1
    //            for i = 1:nx-1
    //                Ex(i,j,k) = Ex(i,j,k)+((Hz(i,j,k)-Hz(i,j-1,k))*cexz(j)-...
    //                                       (Hy(i,j,k)-Hy(i,j,k-1))*cexy(k))*epsxm1(i,j,k);
    //    end
    //    end
    //    end
    //
    //    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //    % Ey-update  (excludes boundaries)
    //    for k = 2:nz-1
    //        for j = 1:ny-1
    //            for i = 2:nx-1
    //                Ey(i,j,k) = Ey(i,j,k)+((Hx(i,j,k)-Hx(i,j,k-1))*ceyx(k)-...
    //                                       (Hz(i,j,k)-Hz(i-1,j,k))*ceyz(i))*epsym1(i,j,k);
    //    end
    //    end
    //    end
    //
    //    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //    % Ez-update  (excludes boundaries)
    //    for k = 1:nz-1
    //        for j = 2:ny-1
    //            for i = 2:nx-1
    //                Ez(i,j,k) = Ez(i,j,k)+((Hy(i,j,k)-Hy(i-1,j,k))*cezy(i)-...
    //                                       (Hx(i,j,k)-Hx(i,j-1,k))*cezx(j))*epszm1(i,j,k);
    //    end
    //    end
    //    end
    for (k=1;k<NX-2;k++){
        for (j=1;j<NY-2;j++){
            for (i=0;i<NZ-2;i++){
                Ex(i,j,k) = Ex(i,j,k) + c_cExy * (Hz(i,j,k) - Hz(i,j-1,k)) -
                c_cExz * (Hy(i,j,k) - Hy(i,j,k-1));
            }
        }
    }
    // Ey update
    for (k=1;k<NX-2;k++){
        for (j=0;j<NY-2;j++){
            for (i=1;i<NZ-2;i++){
                Ey(i,j,k) = Ey(i,j,k) + c_cEyz * (Hx(i,j,k) - Hx(i,j,k-1)) -
                c_cEyx * (Hz(i,j,k) - Hz(i-1,j,k));
            }
        }
    }
    // Ez update
    for (k=0;k<NX-2;k++){
        for (j=1;j<NY-2;j++){
            for (i=1;i<NZ-2;i++){
                Ez(i,j,k) = Ez(i,j,k) + c_cEzx * (Hy(i,j,k) - Hy(i-1,j,k)) -
                c_cEzy * (Hx(i,j,k) - Hx(i,j-1,k));
            }
        }
    }
    
}




void Hupdate(float * Ex, float * Ey, float * Ez, float * Hx, float * Hy, float * Hz,
             float * cHxy, float * cHxz, float * cHyz, float * cHyx, float * cHzx, float * cHzy){
    int i,j,k;
    const float c_cHxy = *cHxy;
    const float c_cHxz = *cHxy;
    const float c_cHyz = *cHyz;
    const float c_cHyx = *cHyx;
    const float c_cHzx = *cHzx;
    const float c_cHzy = *cHzy;

    
    // Hx update
    for (k=1;k<NZ-2;k++){
        for (j=1;j<NY-2;j++){
            for (i=0;i<NX-2;i++){
                Hx(i,j,k) = Hx(i,j,k) - c_cHxy * (Ez(i,j+1,k) - Ez(i,j,k)) +
                c_cHxz * (Ey(i,j,k+1) - Ey(i,j,k));
            }
        }
    }
    // Hy update
    for (k=1;k<NZ-2;k++){
        for (j=0;j<NY-2;j++){
            for (i=1;i<NX-2;i++){
                Hy(i,j,k) = Hy(i,j,k) - c_cHyz * (Ex(i,j,k+1) - Ex(i,j,k)) +
                c_cHyx * (Ez(i+1,j,k) - Ez(i,j,k));
            }
        }
    }
    // Hz update
    for (k=0;k<NZ-2;k++){
        for (j=1;j<NY-2;j++){
            for (i=1;i<NX-2;i++){
                Hz(i,j,k) = Hz(i,j,k) - c_cHzx * (Ey(i+1,j,k) - Ey(i,j,k)) +
                c_cHzy * (Ex(i,j+1,k) - Ex(i,j,k));
            }
        }
    }
    
}










int main(){
    //struct source src;
    
    // Set up the sampling
    //struct sampleLoc sample;
    struct sampleLoc sloc;
    strcpy(sloc.type, "E");
    strcpy(sloc.dir, "z");
    sloc.i = 3;
    sloc.j = 3;
    sloc.k = 3;
    
    // set up the source
    struct source src;
    strcpy(src.type, "J");
    strcpy(src.dir, "z");
    src.i1 = 2;
    src.i2 = 2;
    src.j1 = 2;
    src.j2 = 2;
    src.k1 = 2;
    src.k2 = 3;
    
    // material parameters and speed of light
    float mu_o = 1.2566370614e-6;
    float eps_o = 8.854187817e-12;
    float mu_r = 1;
    float eps_r = 1;
    float mu= mu_r * mu_o;
    float eps = eps_r * eps_o;
    float c = 1.0/sqrt(mu*eps);
    
    // discretization size computation
    float dx=Dx/(NX-1);
    float dy=Dy/(NY-1);
    float dz=Dz/(NZ-1);
    float dt=CFLN/(c*sqrt(pow(dx,-2)+pow(dy,-2)+pow(dz,-2)));
    
    int nt = floor(simtime/dt);
    
    // build the E and H space
    // E-field space
    float * Ex = malloc((NX-1)*NY*NZ*sizeof(float));
    float * Ey = malloc(NX*(NY-1)*NZ*sizeof(float));
    float * Ez = malloc(NX*NY*(NZ-1)*sizeof(float));
    // H-field space
    float * Hx = malloc(NX*(NY-1)*(NZ-1)*sizeof(float));
    float * Hy = malloc((NX-1)*NY*(NZ-1)*sizeof(float));
    float * Hz = malloc((NX-1)*(NY-1)*NZ*sizeof(float));
    
    // build coefficient matrices
    // E-field coefficients
    float cExy = (dt/(eps*dy));
    float cExz = (dt/(eps*dz));
    float cEyx = (dt/(eps*dx));
    float cEyz = (dt/(eps*dz));
    float cEzx = (dt/(eps*dx));
    float cEzy = (dt/(eps*dy));
    // H-field coefficients
    float cHxy = (dt/(mu*dy));
    float cHxz = (dt/(mu*dz));
    float cHyx = (dt/(mu*dx));
    float cHyz = (dt/(mu*dz));
    float cHzx = (dt/(mu*dx));
    float cHzy = (dt/(mu*dy));
    
    // source coefficients
    float cJ = dt/eps;
    float cM = dt/mu;
    
    int i;
    // build an init sampling vector
    int N = pow(2, ceil(log(nt)/log(2)))*2;
    float * sample = malloc(N*sizeof(float));
    for (i=0;i<N;i++){
        sample[i] = 0;
    }
    
    // time (make first update be at t=0)
    float t = -0.5*dt;
    
    int n;
    
    for(n=0 ; n<=nt*2 ; n++){
        // increment time by dt/2
        t=t+0.5*dt;
        // update the E-field
        Eupdate(Ex,Ey,Ez,Hx,Hy,Hz,&cExy,&cExz,&cEyz,&cEyx,&cEzx,&cEzy);
        // increment time by dt/2
        t=t+0.5*dt;
        // update the H-field
        Hupdate(Ex,Ey,Ez,Hx,Hy,Hz,&cHxy,&cHxz,&cHyz,&cHyx,&cHzx,&cHzy);
        // update the sample
        sample[n] = Ez(sloc.i, sloc.j, sloc.k);
        
        if(n%100==0){printf("%f\n",sample[n]); }
        
        //Flat[x + WIDTH * (y + DEPTH * z)] = Original[x, y, z]
    }
    
    
    
    // Free the E and H field memory
    free(Ex);
    free(Ey);
    free(Ez);
    free(Hx);
    free(Hy);
    free(Hz);
    free(sample);
}




