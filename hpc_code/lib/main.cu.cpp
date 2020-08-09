//
//  main.cpp
//  diffusion_3D_RMS
//
//  Update Journal:
//  -- 3/11/2019: equal step length random leap, one fiber, IAS, 3D, cuda version
//  -- 4/25/2019: implement mitochondria, high permeability, short T2, same diffusivity as IAS
//  -- 1/29/2020: implement generalized realistic microstructure simulator (RMS): elastic reflection, water exchange (no permeability) and T2 relaxation
//
//  Created by Hong-Hsi Lee in January, 2020.
//

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <iomanip>
#include <time.h>
#include <cstdlib>
#include <algorithm>
#include <string>
#include <complex>
#include <string>

#include <cuda.h>
#include <curand_kernel.h>
#include <cuda_runtime.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>

using namespace std;
    
#define Pi 3.14159265
#define timepoints 1000
#define nite 4
#define Nc_max 3

// ********** cuda kernel **********
__device__ double atomAdd(double* address, double val)
{
    unsigned long long int* address_as_ull =
    (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    
    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val +
                                             __longlong_as_double(assumed)));
        
        // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
    } while (assumed != old);
    
    return __longlong_as_double(old);
}

__global__ void setup_kernel(curandStatePhilox4_32_10_t *state, unsigned long seed){
    int idx = threadIdx.x+blockDim.x*blockIdx.x;
    curand_init(seed, idx, 0, &state[idx]);
}

__global__ void propagate(curandStatePhilox4_32_10_t *state, double *sig0, double *dx2, double *dx4, double *sigRe, double *NPar_count, const double *TD, const int TN, const int NPar, const int Nbvec, const int Nc, const double res, const double *step, const int NPix1, const int NPix2, const int NPix3, const double *T2, const double *Pij, const double *bval, const double *bvec, const int *APix){
    int idx = threadIdx.x + blockDim.x * blockIdx.x;
    int stride = blockDim.x * gridDim.x;
    curandStatePhilox4_32_10_t localstate=state[idx];
    
    int Tstep=TN/timepoints;

    for (int k=idx; k<NPar; k+=stride){
        // Random number
        double vRand=0;
        
        // Particle position on a grid
        int xParGi[3]={0}, xParGj[3]={0};
        
        // xi: initial particle position, xt: particle position at the i-th step
        // tt: distance between particle and y-z, x-z, x-y box wall
        // vt: a unit vector indicating the hopping direction
        // xTmp: a temporary variable to save the position
        double xi[3]={0}, xt[3]={0}, tt[3]={0}, vt[3]={0};
//        double xTmp[3]={0};
        
        double cos_theta=0, sin_theta=0;
        int tidx=0, nidx=0;
        int NPix[3]={0}; NPix[0]=NPix1; NPix[1]=NPix2; NPix[2]=NPix3;
        
        // Signal weighted by T2 relaxation
        double s0=0;
        
        // The flags of hitting the medium boundary
        bool flip0[3]={false}, flip1[3]={false};
//        bool flip0_tmp[3]={false}, flip1_tmp[3]={false};
//        bool flipx0=false, flipx1=false, flipy0=false, flipy1=false, flipz0=false, flipz1=false;
        
        // q=\gamma * g * \delta (1/µm)
        double qx=0;
        double dx=0, dy=0, dz=0;
        
        // fstep: remaining fraction of one step
        // tmp: temporary variable
        // tmin: the shortest distance between the particle and y-z, x-z, x-y plane
        double fstep=0, tmp=0, tmin=0;
        
        // Elements of lookup table APix
        int ai=0, aj=0;
        
        // The box wall hit by the particle. 1:y-z plane, 2:x-z plane, 3:x-y plane
        int ii_hit = 0;
        
        // The time staying in compartments
        double t[Nc_max]={0};
//        printf("step size=%.4f\n",step[0]);
//        printf("step size=%.4f\n",step[1]);
        
        //********** Initialize Particle Positions in IAS *********
        while (1){
            xi[0]=curand_uniform_double(&localstate)*static_cast<double>(NPix1);
            xi[1]=curand_uniform_double(&localstate)*static_cast<double>(NPix2);
            xi[2]=curand_uniform_double(&localstate)*static_cast<double>(NPix3);
            
            // Whether the particle is inside compartments
            xParGi[0]=floor(xi[0]); xParGi[1]=floor(xi[1]); xParGi[2]=floor(xi[2]);
            if ( APix[ NPix2*NPix3*xParGi[0] + NPix3*xParGi[1] + xParGi[2] ] !=0 ){ break; }
        }
//        printf("position=%.2f, %.2f, %.2f\n",xi[0],xi[1],xi[2]);

        // ********** Simulate diffusion **********
        xt[0]=xi[0]; xt[1]=xi[1]; xt[2]=xi[2];
        
        xParGi[0]=floor(xt[0]); xParGi[1]=floor(xt[1]); xParGi[2]=floor(xt[2]);
        ai=APix[ NPix2*NPix3*xParGi[0] + NPix3*xParGi[1] + xParGi[2] ];
        
        for (int i=0; i<TN; i++){
            fstep = 1.0;
            for (int jj=0; jj<3; jj++) {
                flip0[jj]=false; flip1[jj]=false;
            }
//            flipx0=false; flipx1=false; flipy0=false; flipy1=false; flipz0=false; flipz1=false;
            vRand=curand_uniform_double(&localstate);
            cos_theta=1.0-2.0*vRand;
            sin_theta=2.0*sqrt(vRand*(1.0-vRand)); //sin(acos(cos_theta));
            
            vRand=curand_uniform_double(&localstate);
            
//            xParGi[0]=floor(xt[0]); xParGi[1]=floor(xt[1]); xParGi[2]=floor(xt[2]);
//            ai=APix[ NPix2*NPix3*xParGi[0] + NPix3*xParGi[1] + xParGi[2] ];
            
            vt[0]=sin_theta*cos(2.0*Pi*vRand);
            vt[1]=sin_theta*sin(2.0*Pi*vRand);
            vt[2]=cos_theta;
            
            for (int j=0; (j<nite) && (fstep>0) ; j++){
//                xTmp[0]=xt[0]+step[ai]*fstep*vt[0];
//                xTmp[1]=xt[1]+step[ai]*fstep*vt[1];
//                xTmp[2]=xt[2]+step[ai]*fstep*vt[2];
                
//                for (int jj=0; jj<3; jj++) {
//                    flip0_tmp[jj]=flip0[jj];
//                    flip1_tmp[jj]=flip1[jj];
//                    if (xTmp[jj]<0) {
//                        xTmp[jj]=-xTmp[jj];
//                        flip0_tmp[jj]=true;
//                    }
//                    if (xTmp[jj]>NPix[jj]) {
//                        xTmp[jj]=2.0*static_cast<double>(NPix[jj])-xTmp[jj];
//                        flip1_tmp[jj]=true;
//                    }
//                }
                
//                if (xTmp[0]<0){
//                    xTmp[0]=-xTmp[0]; flipx0=true;
//                }
//                if (xTmp[1]<0){
//                    xTmp[1]=-xTmp[1]; flipy0=true;
//                }
//                if (xTmp[2]<0){
//                    xTmp[2]=-xTmp[2]; flipz0=true;
//                }
//                
//                if (xTmp[0]>NPix1){
//                    xTmp[0]=2.0*static_cast<double>(NPix1)-xTmp[0]; flipx1=true;
//                }
//                if (xTmp[1]>NPix2){
//                    xTmp[1]=2.0*static_cast<double>(NPix2)-xTmp[1]; flipy1=true;
//                }
//                if (xTmp[2]>NPix3){
//                    xTmp[2]=2.0*static_cast<double>(NPix3)-xTmp[2]; flipz1=true;
//                }
                
//                xParGj[0]=floor(xTmp[0]); xParGj[1]=floor(xTmp[1]); xParGj[2]=floor(xTmp[2]);
//                if ( (xParGi[0]==xParGj[0]) && (xParGi[1]==xParGj[1]) && (xParGi[2]==xParGj[2]) ) {
//                    xt[0]=xTmp[0]; xt[1]=xTmp[1]; xt[2]=xTmp[2];
//                    t[ai]=t[ai]+fstep;
//                    for (int jj=0; jj<3; jj++) { flip0[jj]=flip0_tmp[jj]; flip1[jj]=flip1_tmp[jj]; }
//                    break;
//                }
                
                tmin=2.0; ii_hit=-1;
                for (int ii=0;ii<3;ii++) {
                    if (vt[ii]>0.0) {
                        tmp = static_cast<double>(xParGi[ii]) + 1.0 - xt[ii];
                    } else if (vt[ii]<0.0) {
                        tmp = static_cast<double>(xParGi[ii]) - xt[ii];
                    } else {
                        tmp = 2.0;
                    }
                    
//                    tmp=static_cast<double>(xParGi[ii]) + fmax(0.0,static_cast<double>(vt[ii]>0.0)) - xt[ii];
                    if ( fabs(tmp) > fabs(vt[ii]) ) {
                        tt[ii] = 2.0;
                    } else {
                        tt[ii] = max(0.0,tmp/vt[ii]);
                    }
                    if (tt[ii]<tmin) {
                        tmin=tt[ii];
                        ii_hit=ii;
                    }
                }
                
//                if ( ii_hit<0 ) {
//                    printf("Error1: walker does not encounter the box wall, tmin=%.4f, xGi=%i,%i,%i, xGj=%i,%i,%i\n",tmin,xParGi[0],xParGi[1],xParGi[2],xParGj[0],xParGj[1],xParGj[2]);
//                    break;
//                }
                
                if ( tmin<0.0 ) {
                    printf("Error: walker jumps into wrong direction, tmin=%.4f\n",tmin);
                    break;
                }

                if ( (fstep*step[ai]) >= tmin ) {
                    fstep=fstep-tmin/step[ai];
                    xt[0]=xt[0]+tmin*vt[0];
                    xt[1]=xt[1]+tmin*vt[1];
                    xt[2]=xt[2]+tmin*vt[2];
                    t[ai]=t[ai]+tmin/step[ai];
                } else {
                    xt[0]=xt[0]+step[ai]*fstep*vt[0];
                    xt[1]=xt[1]+step[ai]*fstep*vt[1];
                    xt[2]=xt[2]+step[ai]*fstep*vt[2];
                    t[ai]=t[ai]+fstep;
//                    printf("Error3: walker does not encounter the box wall, tmin=%.4f, xGi=%i,%i,%i, xGj=%i,%i,%i, vt=%.4f, ii_hit=%i\n",tmin,xParGi[0],xParGi[1],xParGi[2],xParGj[0],xParGj[1],xParGj[2],vt[ii_hit],ii_hit);
                    break;
                }
                
                xParGj[0]=xParGi[0]; xParGj[1]=xParGi[1]; xParGj[2]=xParGi[2];
//                if ( (~flipx0) && (~flipx1) && (~flipy0) && (~flipy1) && (~flipz0) && (~flipz1) ) {
//                if ( (~flip0[0]) && (~flip1[0]) && (~flip0[1]) && (~flip1[1]) && (~flip0[2]) && (~flip1[2]) ) {
                if (vt[ii_hit]>=0.0) {
                    xParGj[ii_hit]=xParGj[ii_hit]+1;
                    if (xParGj[ii_hit]>=NPix[ii_hit]){
                        xParGj[ii_hit]=NPix[ii_hit]-1;
                        flip1[ii_hit]=true;
                        vt[ii_hit]=-vt[ii_hit];
                        continue;
                    }
                } else {
                    xParGj[ii_hit]=xParGj[ii_hit]-1;
                    if (xParGj[ii_hit]<0) {
                        xParGj[ii_hit]=0;
                        flip0[ii_hit]=true;
                        vt[ii_hit]=-vt[ii_hit];
                        continue;
                    }
                }
                
                if ( ii_hit<0 ) {
                    printf("Error: walker does not encounter the box wall, tmin=%.4f, xGi=%i,%i,%i, xGj=%i,%i,%i\n",tmin,xParGi[0],xParGi[1],xParGi[2],xParGj[0],xParGj[1],xParGj[2]);
                    break;
                }
                
                aj=APix[ NPix2*NPix3*xParGj[0] + NPix3*xParGj[1] + xParGj[2] ];
                
                if (Pij[ai*Nc+aj]>0.999999) {
                    ai=aj;
                    xParGi[ii_hit]=xParGj[ii_hit];
                } else if (Pij[ai*Nc+aj]<0.000001) {
                    vt[ii_hit]=-vt[ii_hit];
//                    printf("ai=%i, aj=%i, xGj=%i,%i,%i\n",ai,aj,xParGj[0],xParGj[1],xParGj[2]);
                } else {
//                    printf("ai=%i, aj=%i, Pij=%.4f\n",ai,aj,Pij[ai*Nc+aj]);
                    vRand=curand_uniform_double(&localstate);
                    if (vRand<Pij[ai*Nc+aj]) {
                        ai=aj;
                        xParGi[ii_hit]=xParGj[ii_hit];
                    } else {
                        vt[ii_hit]=-vt[ii_hit];
                    }
                }
//                }
                
//                for (int jj=0; jj<3; jj++) {
//                    if (flip0[jj] || flip1[jj]) { vt[jj]=-vt[jj]; }
//                }
                
//                if ( flipx0 || flipx1 ) { vt[0]=-vt[0]; }
//                if ( flipy0 || flipy1 ) { vt[1]=-vt[1]; }
//                if ( flipz0 || flipz1 ) { vt[2]=-vt[2]; }
                
            }
            
            for (int jj=0; jj<3; jj++) {
                if (flip0[jj]) {
                    xi[jj]=-xi[jj];
                }
                if (flip1[jj]) {
                    xi[jj]=2.0*static_cast<double>(NPix[jj])-xi[jj];
                }
                flip0[jj]=false; flip1[jj]=false;
            }

//            if (flipx0){
//                xi[0]=-xi[0];
//            }
//            if (flipy0){
//                xi[1]=-xi[1];
//            }
//            if (flipz0){
//                xi[2]=-xi[2];
//            }
//            
//            if (flipx1){
//                xi[0]=2.0*static_cast<double>(NPix1)-xi[0];
//            }
//            if (flipy1){
//                xi[1]=2.0*static_cast<double>(NPix2)-xi[1];
//            }
//            if (flipz1){
//                xi[2]=2.0*static_cast<double>(NPix3)-xi[2];
//            }
//            flipx0=false; flipx1=false; flipy0=false; flipy1=false; flipz0=false; flipz1=false;
            
            if ( (i%Tstep)==0 ) { // Save moment tensor for dx^2 and dx^4, and signal for the b-table
                s0=0.0;
                for (int j=0; j<Nc; j++) {
                    s0=s0+(t[j]/T2[j]);
                }
                s0=exp(-1.0*s0);
//                s0=1.0;
                
                tidx=i/Tstep;
                nidx=Nc*tidx+ai;
                
                atomAdd(&sig0[tidx],s0);
                atomAdd(&NPar_count[nidx],1);
                
                dx=(xt[0]-xi[0])*res;
                dy=(xt[1]-xi[1])*res;
                dz=(xt[2]-xi[2])*res;
                
                atomAdd(&dx2[6*tidx+0],s0*dx*dx);
                atomAdd(&dx2[6*tidx+1],s0*dx*dy);
                atomAdd(&dx2[6*tidx+2],s0*dx*dz);
                atomAdd(&dx2[6*tidx+3],s0*dy*dy);
                atomAdd(&dx2[6*tidx+4],s0*dy*dz);
                atomAdd(&dx2[6*tidx+5],s0*dz*dz);
                
                atomAdd(&dx4[15*tidx+0],s0*dx*dx*dx*dx);
                atomAdd(&dx4[15*tidx+1],s0*dx*dx*dx*dy);
                atomAdd(&dx4[15*tidx+2],s0*dx*dx*dx*dz);
                atomAdd(&dx4[15*tidx+3],s0*dx*dx*dy*dy);
                atomAdd(&dx4[15*tidx+4],s0*dx*dx*dy*dz);
                atomAdd(&dx4[15*tidx+5],s0*dx*dx*dz*dz);
                atomAdd(&dx4[15*tidx+6],s0*dx*dy*dy*dy);
                atomAdd(&dx4[15*tidx+7],s0*dx*dy*dy*dz);
                atomAdd(&dx4[15*tidx+8],s0*dx*dy*dz*dz);
                atomAdd(&dx4[15*tidx+9],s0*dx*dz*dz*dz);
                atomAdd(&dx4[15*tidx+10],s0*dy*dy*dy*dy);
                atomAdd(&dx4[15*tidx+11],s0*dy*dy*dy*dz);
                atomAdd(&dx4[15*tidx+12],s0*dy*dy*dz*dz);
                atomAdd(&dx4[15*tidx+13],s0*dy*dz*dz*dz);
                atomAdd(&dx4[15*tidx+14],s0*dz*dz*dz*dz);
                
                for (int j=0; j<Nbvec; j++) {
                    qx= sqrt( bval[j] / TD[tidx] ) *(dx*bvec[j*3+0]+dy*bvec[j*3+1]+dz*bvec[j*3+2]);
                    atomAdd(&sigRe[Nbvec*tidx+j],s0*cos(qx));
//                    atomAdd(&sigIm[Nbvec*tidx+j],-s0*sin(qx));
                }
            }

        }
    }
    state[idx]=localstate;
}

    //********** Define tissue parameters **********

    int main(int argc, char *argv[]) {
        
        clock_t begin=clock();
        clock_t end=clock();
        
        // Define index number
        int i=0, j=0, k=0;
        
        //********** Load mictostructure **********
        
        double dt=0.0;              // Time step in ms
        int TN=0;                   // # time steps
        int NPar=0;                 // # particles
        int Nbvec=0;                // # gradient directions
        int Nc=0;                   // # compartments
        double res=0.0;             // voxel size
        int NPix1=0, NPix2=0, NPix3=0;  // medium matrix dimension
        
//        string root, target;        // root for loading packing, target for saving files

//        // root
//        ifstream myfile_root ("root.txt", ios::in);
//        myfile_root>>root;
//        myfile_root.close();
        
        // simulation parameter
        ifstream myfile0 ("simParamInput.txt", ios::in);
        myfile0>>dt; myfile0>>TN; myfile0>>NPar;
        myfile0>>Nbvec; myfile0>>Nc;
        myfile0>>res;
        myfile0>>NPix1; myfile0>>NPix2; myfile0>>NPix3;
        myfile0.close();
        
        thrust::host_vector<double> D(Nc);
        ifstream myfile0D ("diffusivity.txt", ios::in);
        for (i=0; i<Nc; i++){
            myfile0D>>D[i];
            cout<<"diffusivity="<<D[i]<<endl;
        }
        myfile0D.close();
        
        thrust::host_vector<double> T2(Nc);
        ifstream myfile0T2 ("T2.txt", ios::in);
        for (i=0; i<Nc; i++){
            myfile0T2>>T2[i];
            T2[i]/=dt;
//            cout<<"T2="<<T2[i]<<endl;
        }
        myfile0T2.close();
        
        thrust::host_vector<double> step(Nc); // Step size in IAS in µm
        for (i=0; i<Nc; i++){
            step[i]=sqrt(6.0*dt*D[i])/res;
//            cout<<"step size="<<step[i]<<endl;
        }
        
        thrust::host_vector<double> Pij(Nc*Nc);
        k=0;
        for (i=0; i<Nc; i++){
            for (j=0; j<Nc; j++){
                if ( (i==0) || (j==0) ) {
                    Pij[k]=0.0;
                } else if ( i==j ) {
                    Pij[k]=1.0;
                } else {
                    Pij[k]=fmin(1.0,sqrt(D[j]/D[i]));
                }
                cout<<"permeation probability="<<Pij[k]<<endl;
                k++;
            }
        }
        
        // Pixelized matrix A for the fiber
        thrust::host_vector<int> APix(NPix1*NPix2*NPix3);
        ifstream myfile1 ("fiber.txt", ios::in);
        for (i=0; i<NPix1*NPix2*NPix3; i++){
            myfile1>>APix[i];
        }
        myfile1.close();
        
        // b-value
        thrust::host_vector<double> bval(Nbvec);
        ifstream myfile2 ("bval.txt", ios::in);
        for (i=0; i<Nbvec; i++){
            myfile2>>bval[i];
        }
        myfile2.close();
        
        // b-vector
        thrust::host_vector<double> bvec(Nbvec*3);
        ifstream myfile3 ("bvec.txt", ios::in);
        for (i=0; i<Nbvec; i++){
            myfile3>>bvec[i*3+0]; myfile3>>bvec[i*3+1]; myfile3>>bvec[i*3+2];
        }
        myfile3.close();
        
        thrust::host_vector<double> TD(timepoints);
        for (i=0; i<timepoints; i++){
        TD[i]=(i*(TN/timepoints)+1)*dt;
        }

        // ********** Simulate diffusion **********
        
        // Initialize seed
        unsigned long seed=0;
        FILE *urandom;
        urandom = fopen("/dev/random", "r");
        fread(&seed, sizeof (seed), 1, urandom);
        fclose(urandom);
        
        // Initialize state of RNG
        int blockSize = 64;
        int numBlocks = (NPar + blockSize - 1) / blockSize;
        cout<<numBlocks<<endl<<blockSize<<endl;
        
        thrust::device_vector<curandStatePhilox4_32_10_t> devState(numBlocks*blockSize);
        setup_kernel<<<numBlocks, blockSize>>>(devState.data().get(),seed);
        
        // Initialize output
        thrust::host_vector<double> sig0(timepoints);
        thrust::host_vector<double> dx2(timepoints*6);
        thrust::host_vector<double> dx4(timepoints*15);
        thrust::host_vector<double> sigRe(timepoints*Nbvec);
//        thrust::host_vector<double> sigIm(timepoints*Nbvec);
        thrust::host_vector<double> NPar_count(timepoints*Nc);
        for (i=0;i<timepoints;i++){ sig0[i]=0; }
        for (i=0;i<timepoints*6;i++){ dx2[i]=0; }
        for (i=0;i<timepoints*15;i++){ dx4[i]=0; }
        for (i=0;i<timepoints*Nbvec;i++){ sigRe[i]=0; }
        for (i=0;i<timepoints*Nc;i++){ NPar_count[i]=0; }
        
        // Move data from host to device
        thrust::device_vector<double> d_sig0=sig0;
        thrust::device_vector<double> d_dx2=dx2;
        thrust::device_vector<double> d_dx4=dx4;
        thrust::device_vector<double> d_sigRe=sigRe;
//        thrust::device_vector<double> d_sigIm=sigIm;
        thrust::device_vector<double> d_NPar_count = NPar_count;
        thrust::device_vector<double> d_TD=TD;
        thrust::device_vector<double> d_step=step;
        thrust::device_vector<double> d_T2=T2;
        thrust::device_vector<double> d_Pij=Pij;
        thrust::device_vector<double> d_bval=bval;
        thrust::device_vector<double> d_bvec=bvec;
        thrust::device_vector<int> d_APix=APix;
        
//        double *NPar_count; cudaMallocManaged(&NPar_count,sizeof(double)); NPar_count[0] = 0;
        // Parallel computation
        begin=clock();
        propagate<<<numBlocks, blockSize>>>(devState.data().get(), d_sig0.data().get(), d_dx2.data().get(), d_dx4.data().get(), d_sigRe.data().get(), d_NPar_count.data().get(), d_TD.data().get(), TN, NPar, Nbvec, Nc, res, d_step.data().get(), NPix1, NPix2, NPix3, d_T2.data().get(), d_Pij.data().get(), d_bval.data().get(), d_bvec.data().get(), d_APix.data().get());
        cudaDeviceSynchronize();
        end=clock();
        cout << "Done! Elpased time "<<double((end-begin)/CLOCKS_PER_SEC) << " s"<< endl;
        
        thrust::copy(d_sig0.begin(), d_sig0.end(), sig0.begin());
        thrust::copy(d_dx2.begin(), d_dx2.end(), dx2.begin());
        thrust::copy(d_dx4.begin(), d_dx4.end(), dx4.begin());
        thrust::copy(d_sigRe.begin(), d_sigRe.end(), sigRe.begin());
//        thrust::copy(d_sigIm.begin(), d_sigIm.end(), sigIm.begin());
        thrust::copy(d_NPar_count.begin(), d_NPar_count.end(), NPar_count.begin());
        
        // Save results
        ofstream fs0out("sig0.txt");
        ofstream fdx2out("dx2.txt");
        ofstream fdx4out("dx4.txt");
        ofstream fsRout("sigRe.txt");
//        ofstream fsIout("sigIm.txt");
        ofstream fNpout("NPar_count.txt");
        fs0out.precision(15);
        fdx2out.precision(15);
        fdx4out.precision(15);
        fsRout.precision(15);
//        fsIout.precision(15);
        for (i=0; i<timepoints; i++) {
            fs0out<<sig0[i]<<endl;
            for (j=0;j<6;j++){
                if (j==5){
                    fdx2out<<dx2[i*6+j]<<endl;
                }
                else {
                    fdx2out<<dx2[i*6+j]<<"\t";
                }
            }
            for (j=0;j<15;j++){
                if (j==14){
                    fdx4out<<dx4[i*15+j]<<endl;
                }
                else {
                    fdx4out<<dx4[i*15+j]<<"\t";
                }
            }
            
            for (j=0; j<Nbvec; j++) {
                if (j==Nbvec-1) {
                    fsRout<<sigRe[i*Nbvec+j]<<endl;
//                    fsIout<<sigIm[i*Nbvec+j]<<endl;
                }
                else {
                    fsRout<<sigRe[i*Nbvec+j]<<"\t";
//                    fsIout<<sigIm[i*Nbvec+j]<<"\t";
                }
            }
            
            for (j=0; j<Nc; j++) {
                if (j==Nc-1) {
                    fNpout<<NPar_count[i*Nc+j]<<endl;
                }
                else {
                    fNpout<<NPar_count[i*Nc+j]<<"\t";
                }
            }
        }
        fdx2out.close();
        fdx4out.close();
        fsRout.close();
//        fsIout.close();
        fNpout.close();
        
        ofstream paraout("sim_para.txt");
        paraout<<dt<<endl<<TN<<endl<<NPar<<endl;
        paraout<<res<<endl;
        paraout.close();
        
        ofstream TDout("diff_time.txt");
        for (i=0; i<timepoints; i++){
            TDout<<TD[i]<<endl;
        }
        TDout.close();
        
        
//        ofstream NParout("NPar_count.txt");
//        NParout<<NPar_count[0];
//        NParout.close();
}

