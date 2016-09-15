#include "gpufiltering.h"
#define SAMPLE_STEP 8
// k_p(p,q)
__global__ void filtering_kernel_p(float * X, float * values, float * sigmas, float * output, size_t N) 
{
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    float p_c[3];
    float sigmasquared;
    if(idx>=N)
        return;
    p_c[0] = X[idx*3];
    p_c[1] = X[idx*3+1];
    p_c[2] = X[idx*3+2];
    sigmasquared = sigmas[idx]*sigmas[idx]*(-2.0);
    output[idx] = 1;
    float diff = 0;
    float out = 0;
    for(int n=0;n<N;n+=SAMPLE_STEP){
        diff=(p_c[0]-X[n*3])*(p_c[0]-X[n*3])+(p_c[1]-X[n*3+1])*(p_c[1]-X[n*3+1])+(p_c[2]-X[n*3+2])*(p_c[2]-X[n*3+2]);
        out += __expf(diff/sigmasquared)*values[n];
    }
    output[idx] = out*SAMPLE_STEP;
}

// k_q(p,q)
__global__ void filtering_kernel_q(float * X, float * values, float * sigmas, float * output, size_t N) 
{
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    float p_c[3];
    if(idx>=N)
        return;
    p_c[0] = X[idx*3];
    p_c[1] = X[idx*3+1];
    p_c[2] = X[idx*3+2];
    output[idx] = 1;
    float diff = 0;
    float out = 0;
    for(int n=0;n<N;n+=SAMPLE_STEP){
        diff=(p_c[0]-X[n*3])*(p_c[0]-X[n*3])+(p_c[1]-X[n*3+1])*(p_c[1]-X[n*3+1])+(p_c[2]-X[n*3+2])*(p_c[2]-X[n*3+2]);
        out += __expf(diff/sigmas[n]/sigmas[n]/(-2.0))*values[n];
    }
    output[idx] = out*SAMPLE_STEP;
}

// k_pq(p,q)
__global__ void filtering_kernel_pq(float * X, float * values, float * sigmas, float * output, size_t N) 
{
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    float p_c[3];
    if(idx>=N)
        return;
    p_c[0] = X[idx*3];
    p_c[1] = X[idx*3+1];
    p_c[2] = X[idx*3+2];
    output[idx] = 1;
    float diff = 0;
    float out = 0;
    float sigma_p = sigmas[idx];
    for(int n=0;n<N;n+=SAMPLE_STEP){
        diff=(p_c[0]-X[n*3])*(p_c[0]-X[n*3])+(p_c[1]-X[n*3+1])*(p_c[1]-X[n*3+1])+(p_c[2]-X[n*3+2])*(p_c[2]-X[n*3+2]);
        out += __expf(diff/sigmas[n]/sigma_p/(-2.0))*values[n];
    }
    output[idx] = out*SAMPLE_STEP;
}

int gpufilteringdemo(int N) {
    cout << "filtering for N = " <<N<< endl; // prints !!!Hello World!!!
    srand(time(NULL));
    clock_t startTime,endTime;
    // intialization
    startTime = clock();
    float * X = new float[N*3];
    float * values = new float[N];
    float * sigmas = new float[N];
    for(int n=0;n<N;n++){
        X[n*3] = (rand()%1000)/1000.0;
        X[n*3+1] = (rand()%1000)/1000.0;
        X[n*3+2] = (rand()%1000)/1000.0;
        values[n] = (rand()%1000)/1000.0;
        sigmas[n] = 0.01+(rand()%100)/100000.0;
    }
    endTime = clock();
    cout << "Initialization: "<<float( endTime- startTime ) / (float)CLOCKS_PER_SEC<< " seconds." << endl;
    // parallel dense filtering
    startTime = clock();
    float *Xd, *valuesd, *sigmasd, * outputd;
    cudaMalloc((void **)&Xd, sizeof(float)*N*3);
    cudaMalloc((void **)&valuesd, sizeof(float)*N);
    cudaMalloc((void **)&sigmasd, sizeof(float)*N);
    cudaMalloc((void **)&outputd, sizeof(float)*N);
    cudaMemcpy(Xd, X, sizeof(float)*N*3, cudaMemcpyHostToDevice);
    cudaMemcpy(valuesd, values, sizeof(float)*N, cudaMemcpyHostToDevice);
    cudaMemcpy(sigmasd, sigmas, sizeof(float)*N, cudaMemcpyHostToDevice);
    endTime = clock();
    cout << "memory copy: "<<float( endTime- startTime ) / (float)CLOCKS_PER_SEC<< " seconds." << endl;
    int num_block = ((N%BLOCK_SIZE)==0)?(N/BLOCK_SIZE):(N/BLOCK_SIZE+1);
    startTime = clock();
    cout<<"number of blocks: "<<num_block<<endl;
    // launch the kernel
    filtering_kernel_q<<<num_block, BLOCK_SIZE>>>(Xd, valuesd, sigmasd, outputd, N);
    cudaThreadSynchronize();
    cudaError_t code = cudaGetLastError();
    if(code!=cudaSuccess){
        printf("Cuda error -- %s\n",cudaGetErrorString(code));
    } else printf("Cuda success\n");
    float * para_output = new float[N];
    cudaMemcpy(para_output, outputd, sizeof(float)*N, cudaMemcpyDeviceToHost);
    endTime = clock();
    float kernel1_time = float( endTime- startTime ) / (float)CLOCKS_PER_SEC;
    cout << "kernel time : "<<kernel1_time<< " seconds." << endl;

    delete [] X;
    delete [] values;
    delete [] sigmas;
    delete [] para_output;
    cudaFree(Xd);
    cudaFree(valuesd);
    cudaFree(sigmasd);
    cudaFree(outputd);
    printf("End of program.\n");
    return 0;
}

float * gpufiltering(float * X, float * values, float * sigmas, size_t N, char flag)
{
	//int devicecount=0;
	//cudaGetDeviceCount(&devicecount);
	//printf("devicec count %d\n",devicecount);
    float * output = new float[N];
    float *Xd, *valuesd, *sigmasd, * outputd;
    cudaMalloc((void **)&Xd, sizeof(float)*N*3);
    cudaMalloc((void **)&valuesd, sizeof(float)*N);
    cudaMalloc((void **)&sigmasd, sizeof(float)*N);
    cudaMalloc((void **)&outputd, sizeof(float)*N);
    cudaMemcpy(Xd, X, sizeof(float)*N*3, cudaMemcpyHostToDevice);
    cudaMemcpy(valuesd, values, sizeof(float)*N, cudaMemcpyHostToDevice);
    cudaMemcpy(sigmasd, sigmas, sizeof(float)*N, cudaMemcpyHostToDevice);
    int num_block = ((N%BLOCK_SIZE)==0)?(N/BLOCK_SIZE):(N/BLOCK_SIZE+1);
    // launch the kernel
    if(flag=='p')
        filtering_kernel_p<<<num_block, BLOCK_SIZE>>>(Xd, valuesd, sigmasd, outputd, N);
    else if(flag=='q')
        filtering_kernel_q<<<num_block, BLOCK_SIZE>>>(Xd, valuesd, sigmasd, outputd, N);
    else if(flag=='m')
        filtering_kernel_pq<<<num_block, BLOCK_SIZE>>>(Xd, valuesd, sigmasd, outputd, N);
    cudaThreadSynchronize();
    cudaError_t code = cudaGetLastError();
    if(code!=cudaSuccess){
        printf("Cuda error -- %s\n",cudaGetErrorString(code));
    } //else printf("Cuda success\n");
    cudaMemcpy(output, outputd, sizeof(float)*N, cudaMemcpyDeviceToHost);
    // deallocation
    cudaFree(Xd);
    cudaFree(valuesd);
    cudaFree(sigmasd);
    cudaFree(outputd);
    return output;
}
