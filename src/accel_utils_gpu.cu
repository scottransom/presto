//GPU functions for accelsearch
//by Jintao Luo, NRAO

// includes, CUDA
#include <cuda_runtime.h>
#include <cufft.h>
#include <helper_functions.h>
#include <helper_cuda.h>
#include <stdio.h>

#include "device_functions.h"

#include "accel_utils_gpu.h"

typedef float2 Complex;

//define a texture memory
texture<Complex> tex_d_kernel;
texture<Complex> tex_d_data;

texture<Complex> tex_d_result;

texture<Complex, 2> tex_d_result_2D;
cudaChannelFormatDesc channelDesc ;

texture<unsigned short> tex_d_zinds;
texture<unsigned short> tex_d_rinds;

texture<float> tex_d_fundamental;

static __device__ __host__ inline fcomplex ComplexScale(fcomplex, float);
static __device__ __host__ inline fcomplex ComplexMul(fcomplex, fcomplex);
static __device__ __host__ inline fcomplex ComplexMul_02(fcomplex, fcomplex);

static __global__ void ComplexPointwiseMulAndScale_one_loop(fcomplex *c, fcomplex *a, const fcomplex *b, int numkern_in_array, int data_size, float scale, int kernel_array_offset);

static __global__ void ComplexPointwiseMulAndScale_one_loop_02(fcomplex *c, fcomplex *a, const fcomplex *b, int numkern_in_array, int data_size, float scale, int kernel_array_offset);

static __global__ void Complex_Pow_and_Chop(float *d_pow, fcomplex *d_result, int fftlen, int numkern_in_array, int chopbins, int numtocopy);

static __global__ void add_ffdotpows_on_gpu(float *d_fundamental, fcomplex *d_result, int numzs_full, int numrs_full, unsigned short *zinds, unsigned short *rinds, int fftlen, int chopbins, double obs_zlo, double fullrlo, double harm_fract, int zlo, int rlo);


extern "C"
extern int  search_ffdotpows_gpu(float powcut, float *d_fundamental, accel_cand_gpu * cand_array_search_gpu, accel_cand_gpu * cand_array_sort_gpu, int numzs, int numrs, accel_cand_gpu *cand_gpu_cpu);

static __global__ void  search_ffdotpows_kernel(float powcut, float *d_fundamental, accel_cand_gpu * cand_array_search_gpu, int numzs, int numrs, int *d_addr_jluo);

extern "C"
fcomplex * prep_data_on_gpu(subharminfo **subharminfs, int numharmstages);

extern "C"
fcomplex * prep_result_on_gpu(subharminfo **subharminfs, int numharmstages);

extern "C"
fcomplex * cp_kernel_array_to_gpu(subharminfo **subharminfs, int numharmstages, int **offset_array);

extern "C"
fcomplex * cp_input_to_gpu(fcomplex *input_vect_on_cpu, long long numbins, long long N);

extern "C"
unsigned short *prep_rz_inds_on_gpu(int size_inds);

extern "C"
float * prep_float_vect_on_gpu( int size);

extern "C"
accel_cand_gpu *prep_cand_array(int size);

extern "C"
void cudaFree_kernel_vect(fcomplex *in);

extern "C"
void select_cuda_dev(int cuda_inds);

extern "C"
void complex_corr_conv_gpu(fcomplex * data, fcomplex * kernel_vect_on_gpu,
                       int numdata, 
                       int numkern_in_array, 
                       int stage, int harmtosum, int harmnum, 
                       int ** offset_array, fcomplex *d_data, fcomplex *d_result,
                       int chopbins, int numtocopy,
                       unsigned short *zinds, unsigned short *rinds,
                       int numzs_full, int numrs_full,
                       float *d_fundamental,
                       unsigned short *d_zinds, unsigned short *d_rinds,
                       int datainf_flag,
                       presto_ffts ffts, presto_optype type,
                       double obs_zlo, double fullrlo, double harm_fract, int zlo, int rlo);


extern "C"
void init_cuFFT_plans(subharminfo **subharminfs, int numharmstages);

extern "C"
void destroy_cuFFT_plans(subharminfo **subharminfs, int numharmstages);

cufftHandle plan_data_array[16][16];		
cufftHandle plan_result_array[16][16];		

/******************************************** complex_corr_conv ********************************************************************/
void complex_corr_conv_gpu(fcomplex * data, fcomplex * kernel_vect_on_gpu,
                       int numdata, 
                       int numkern_in_array, 
                       int stage, int harmtosum, int harmnum, 
                       int ** offset_array, fcomplex *d_data, fcomplex *d_result,
                       int chopbins, int numtocopy,
                       unsigned short *zinds, unsigned short *rinds,
                       int numzs_full, int numrs_full,
                       float *d_fundamental,
                       unsigned short *d_zinds, unsigned short *d_rinds,
                       int datainf_flag,
                       presto_ffts ffts, presto_optype type,
                       double obs_zlo, double fullrlo, double harm_fract, int zlo, int rlo)
{


	int fftlen = numdata;
	int kernel_array_offset;

   if (ffts > 3) {
      printf("\nIllegal 'ffts' option (%d) in complex_corr_conv().\n", ffts);
      printf("Exiting.\n\n");
      exit(1);
   }
   if (type > 3) {
      printf("\nIllegal 'type' option (%d) in complex_corr_conv().\n", type);
      printf("Exiting.\n\n");
      exit(1);
   }
	
	if(harmtosum==1 && harmnum==1){
		kernel_array_offset = offset_array[0][0];
	}
	if(harmtosum > 1){
    kernel_array_offset = offset_array[stage][harmnum-1];
	}
	 
	//copy data to GPU memory
	checkCudaErrors(cudaMemcpy(d_data, data, sizeof(fcomplex) * fftlen, cudaMemcpyHostToDevice));
	//FFT data on GPU
	if(datainf_flag == 1)
	{					
		if(harmtosum==1 && harmnum==1){
			checkCudaErrors(cufftExecC2C(plan_data_array[0][0], (cufftComplex *)d_data, (cufftComplex *)d_data, CUFFT_FORWARD));
		}
		if(harmtosum > 1){
  	  checkCudaErrors(cufftExecC2C(plan_data_array[stage][harmnum-1], (cufftComplex *)d_data, (cufftComplex *)d_data, CUFFT_FORWARD));
		}

	}

	//Bind data and kernel to Texture Memory
	cudaBindTexture(NULL, tex_d_kernel, kernel_vect_on_gpu, sizeof(fcomplex) * ( kernel_array_offset + 1 + numkern_in_array * fftlen ) );
	cudaBindTexture(NULL, tex_d_data, d_data, sizeof(fcomplex) * fftlen);
	//Mul the FFTed data with Kernels
	if (type == CORR || type == INPLACE_CORR) {
		ComplexPointwiseMulAndScale_one_loop<<<512, 512>>>(d_result, d_data, kernel_vect_on_gpu, numkern_in_array, fftlen, 1.0/fftlen, kernel_array_offset);
	}
	else {
		ComplexPointwiseMulAndScale_one_loop_02<<<512, 512>>>(d_result, d_data, kernel_vect_on_gpu, numkern_in_array, fftlen, 1.0/fftlen, kernel_array_offset);
	}	
	//unbind the data and kenerl from Texture memory
	cudaUnbindTexture(tex_d_kernel);
	cudaUnbindTexture(tex_d_data);  	

   
  //Inverse FFT   
	if(harmtosum==1 && harmnum==1){
		checkCudaErrors(cufftExecC2C(plan_result_array[0][0], (cufftComplex *)d_result, (cufftComplex *)d_result, CUFFT_INVERSE));	
	}
	if(harmtosum > 1){
 	  checkCudaErrors(cufftExecC2C(plan_result_array[stage][harmnum-1], (cufftComplex *)d_result, (cufftComplex *)d_result, CUFFT_INVERSE));	
	} 

	//bind the FFTed result
		//cudaBindTexture(NULL, tex_d_result, d_result, sizeof(fcomplex) * fftlen * numkern_in_array );	
		//bind d_result to its 2D texture memory
		checkCudaErrors(cudaBindTexture2D(NULL, tex_d_result_2D, d_result, channelDesc, fftlen, numkern_in_array, sizeof(fcomplex) * fftlen ));
	//sum harmonics
	if(harmtosum==1 && harmnum==1){//if fundamental
		Complex_Pow_and_Chop<<<512, 512>>>(d_fundamental, d_result, fftlen, numkern_in_array, chopbins, numtocopy);
	}
	if(harmtosum > 1){	//if harmonics
		//move zinds and rinds to GPU
		checkCudaErrors(cudaMemcpy(d_zinds, zinds, sizeof(unsigned short) * numzs_full, cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemcpy(d_rinds, rinds, sizeof(unsigned short) * numrs_full, cudaMemcpyHostToDevice));		
		//bind zinds and rinds to Texture Memory
		cudaBindTexture(NULL, tex_d_zinds, d_zinds, sizeof(unsigned short) * numzs_full );
		cudaBindTexture(NULL, tex_d_rinds, d_rinds, sizeof(unsigned short) * numrs_full );		
		//add_ffdotpows_on_gpu
		add_ffdotpows_on_gpu<<<512, 512>>>(d_fundamental, d_result, numzs_full, numrs_full, d_zinds, d_rinds, fftlen, chopbins, obs_zlo, fullrlo, harm_fract, zlo, rlo);
		//Unbind zinds and rinds from Texture Memory		
		cudaUnbindTexture(tex_d_zinds);    
		cudaUnbindTexture(tex_d_rinds);    
	}
 
	//cudaUnbindTexture(tex_d_result);
	checkCudaErrors(cudaUnbindTexture(tex_d_result_2D));
  
}                   

/******************************************** add fftdot pows on GPU, choping included ************************************************/
static __global__ void add_ffdotpows_on_gpu(float *d_fundamental, fcomplex *d_result, int numzs_full, int numrs_full, unsigned short *zinds, unsigned short *rinds, int fftlen, int chopbins, double obs_zlo, double fullrlo, double harm_fract, int zlo, int rlo)
{

    const int numThreads = blockDim.x * gridDim.x;
    const int threadID = blockIdx.x * blockDim.x + threadIdx.x;		
    
    //int addr_z, addr_r, addr_result, addr_fundamental ;
    int addr_z, addr_r ;				
    int addr_result ;
		
		Complex buf;
		
		int z_index, r_index;
		
		for (int i = threadID; i < numrs_full * numzs_full; i += numThreads){    		
						
				z_index = i/numrs_full ;
				addr_z = tex1Dfetch(tex_d_zinds, z_index) ;	
				
				r_index = i -  z_index * numrs_full ;				
				addr_r = tex1Dfetch(tex_d_rinds, r_index) ;
				
				//addr_result = addr_z * fftlen + chopbins + addr_r ;
				//buf = tex1Dfetch(tex_d_result, addr_result) ;    		
				
				buf = tex2D(tex_d_result_2D, chopbins + addr_r, addr_z);
				
				d_fundamental[i] += buf.x * buf.x + buf.y * buf.y ;
				
		}

}




/********************************************* Calculate complex ffdot pows, choping included *****************************************/
static __global__ void Complex_Pow_and_Chop(float *d_pow, fcomplex *d_result, int fftlen, int numkern_in_array, int chopbins, int numtocopy)
/*
	calaulate POW of data in d_result
	d_pw:				where the result stored
	d_result:		the array contains data to Pow_and_chop
	fftlen:			length of FFT used in the complex_corr_conv_gpu
	numkern_in_array:	num of kernels in array in complex_corr_conv_gpu
	chopbins:		num of bins that should be discarded in the head and tail of the result vectors
	numtocopy:	num of data should be kept in each vector of d_result, after the first chopbins points	
*/
{
    const int numThreads = blockDim.x * gridDim.x;
    const int threadID = blockIdx.x * blockDim.x + threadIdx.x;		
    
    int addr_result;
    
    Complex buf;
    
    int z_ind, r_ind ;
    
    for (int i = threadID; i < numkern_in_array*numtocopy; i += numThreads){
    	
    	z_ind = i/numtocopy;
    	r_ind = i - z_ind * numtocopy;
    	
    	//addr_result = z_ind * fftlen + chopbins +  r_ind;
    	
    	//buf = tex1Dfetch(tex_d_result, addr_result) ;
    
    	buf = tex2D(tex_d_result_2D, chopbins + r_ind, z_ind);
    
    	d_pow[i] = buf.x * buf.x + buf.y * buf.y ;
    }   

}


/**********************************************************************************************************************************************************/
//complex operations

// Complex scale
static __device__ __host__ inline fcomplex ComplexScale(fcomplex a, float s)
{
    fcomplex c;
    c.r = s * a.r;
    c.i = s * a.i;
    return c;
}

// Complex multiplication
static __device__ __host__ inline fcomplex ComplexMul(fcomplex a, fcomplex b)
{
    fcomplex c;
    c.r = a.r * b.r + a.i * b.i;
    c.i = a.i * b.r - a.r * b.i;
    return c;
}

static __device__ __host__ inline fcomplex ComplexMul_02(fcomplex a, fcomplex b)
{
    fcomplex c;
    c.r = a.r * b.r - a.i * b.i;
    c.i = a.i * b.r + a.r * b.i;
    return c;
}


//Use one loop to realize Complex Point Mul and Scale 
static __global__ void ComplexPointwiseMulAndScale_one_loop(fcomplex *c, fcomplex *a, const fcomplex *b, int numkern_in_array, int data_size, float scale, int kernel_array_offset)
{
    const int numThreads = blockDim.x * gridDim.x;
    const int threadID = blockIdx.x * blockDim.x + threadIdx.x;		
		
		fcomplex a_buf ;
		fcomplex b_buf ;
		Complex buf_buf_a, buf_buf_b;
		int i;
		
		int a_index, b_index;
		
		int total_num = data_size*numkern_in_array;
		
		for (i = threadID; i < total_num; i += numThreads)
		{
				a_index = i/data_size;
				a_index = i - a_index * data_size ;
				buf_buf_a = tex1Dfetch(tex_d_data, a_index);
    		a_buf.r = buf_buf_a.x;
    		a_buf.i = buf_buf_a.y;
    		
    		b_index = i + kernel_array_offset;
    		buf_buf_b = tex1Dfetch(tex_d_kernel, b_index) ;
				
				b_buf.r = buf_buf_b.x;
	    	b_buf.i = buf_buf_b.y;
        c[i] = ComplexScale(ComplexMul(a_buf, b_buf), scale);                                
        	
        
		}

}

static __global__ void ComplexPointwiseMulAndScale_one_loop_02(fcomplex *c, fcomplex *a, const fcomplex *b, int numkern_in_array, int data_size, float scale, int kernel_array_offset)
{
    const int numThreads = blockDim.x * gridDim.x;
    const int threadID = blockIdx.x * blockDim.x + threadIdx.x;		
		
		fcomplex a_buf ;
		fcomplex b_buf ;
		Complex buf_buf_a, buf_buf_b;
		int i;

		int a_index, b_index;
		
		int total_num = data_size*numkern_in_array;
		
		for (i = threadID; i < total_num; i += numThreads)
		{
				a_index = i/data_size;
				a_index = i - a_index * data_size ;
				buf_buf_a = tex1Dfetch(tex_d_data, a_index);
    		a_buf.r = buf_buf_a.x;
    		a_buf.i = buf_buf_a.y;
    		
    		b_index = i + kernel_array_offset;
    		buf_buf_b = tex1Dfetch(tex_d_kernel, b_index) ;
				
				b_buf.r = buf_buf_b.x;
	    	b_buf.i = buf_buf_b.y;
        c[i] = ComplexScale(ComplexMul_02(a_buf, b_buf), scale);                                
        	
        
		}		

}
/**********************************************************************************************************************************************************/

//-------------------------Prepare search_ffdotpows cand array
accel_cand_gpu *prep_cand_array(int size)
{
	
	accel_cand_gpu * cand_on_gpu;
	checkCudaErrors(cudaMalloc((void **)&cand_on_gpu, size * sizeof(accel_cand_gpu)));
	
	return cand_on_gpu ;
	
}

int  search_ffdotpows_gpu(float powcut, float *d_fundamental, accel_cand_gpu * cand_array_search_gpu, accel_cand_gpu * cand_array_sort_gpu, int numzs, int numrs, accel_cand_gpu *cand_gpu_cpu)
{
	
	int *d_addr;
	int h_addr;
	
	checkCudaErrors(cudaMalloc((void **)&d_addr, sizeof(int) * 1));	
	checkCudaErrors(cudaMemset(d_addr, 0, sizeof(int))); // set d_addr to 0

	//bind d_fundamental to texture
	cudaBindTexture(NULL, tex_d_fundamental, d_fundamental, sizeof(float) * numzs * numrs );

	//search ffdot_pow
	search_ffdotpows_kernel<<<512, 512>>>(powcut, d_fundamental, cand_array_search_gpu, numzs, numrs, d_addr);

	//get nof_cand
	checkCudaErrors(cudaMemcpy(&h_addr, d_addr, sizeof(int) * 1, cudaMemcpyDeviceToHost));	
	int nof_cand ;		
	nof_cand = h_addr;		

	//get the candicates
	checkCudaErrors(cudaMemcpy(cand_gpu_cpu, cand_array_search_gpu, sizeof(accel_cand_gpu) * nof_cand, cudaMemcpyDeviceToHost));
	
	cudaFree(d_addr);
	
	cudaUnbindTexture(tex_d_fundamental);  
	
	return nof_cand ;
}

static __global__ void  search_ffdotpows_kernel(float powcut, float *d_fundamental, accel_cand_gpu * cand_array_search_gpu, int numzs, int numrs, int *d_addr_jluo)
{

    const int numThreads = blockDim.x * gridDim.x;
    const int threadID = blockIdx.x * blockDim.x + threadIdx.x;		

		int i ;		
		int nof_cand = 0;

		float pow ;
				
		int addr_search=0;	

		accel_cand_gpu cand_tmp ;
		
		int z_ind, r_ind;

		for (i = threadID; i < numzs*numrs; i += numThreads)
    {    	
    
    	nof_cand = 0;   	    	

   		pow = tex1Dfetch(tex_d_fundamental, i);

    		if(pow > powcut)
    		{    			
    			cand_tmp.pow = pow ;
    			nof_cand += 1 ;
    			
    			cand_tmp.nof_cand = nof_cand ;
    			
    			z_ind = (int)(i/numrs);
    			cand_tmp.z_ind = z_ind;
    			
    			r_ind = i - z_ind * numrs ;
    			cand_tmp.r_ind = r_ind;  			

    			addr_search = atomicAdd(&d_addr_jluo[0], 1);    			
    			
    			cand_array_search_gpu[ addr_search ] = cand_tmp ;    			   			    			
    		}
    }
}

//-------------------------Prepare vectors for zinds and rinds
unsigned short *prep_rz_inds_on_gpu(int size_inds)
{

	unsigned short *inds;
	checkCudaErrors(cudaMalloc((void **)&inds, size_inds * sizeof(unsigned short)));

	return inds;

}

//-------------------------Prepare a float-point vectors 
float * prep_float_vect_on_gpu( int size)
{

	float * f_vect;
	checkCudaErrors(cudaMalloc((void **)&f_vect, size * sizeof(float)));
	
	return f_vect;

}

//----------Copy input data to GPU---------------------------------
fcomplex * cp_input_to_gpu(fcomplex *input_vect_on_cpu, long long numbins, long long N)
{
	fcomplex *input_vect_on_gpu;		
	
  checkCudaErrors(cudaMalloc((void **)&input_vect_on_gpu, sizeof(fcomplex) * numbins));       		
	
	//check some parameters
	printf("\nNum of data points: %d\n", N);
	printf("\nNum of bins : %d\n", numbins);	
	
	checkCudaErrors(cudaMemcpy( input_vect_on_gpu, input_vect_on_cpu, sizeof(fcomplex) * numbins, cudaMemcpyHostToDevice) );
	
	return input_vect_on_gpu;
}

//----------Prepare d_data on GPU---------------------
fcomplex * prep_data_on_gpu(subharminfo **subharminfs, int numharmstages)
//cudaMalloc d_data for the whole search process
{

	int harm, harmtosum, stage;		
	int fftlen;
	int size_d_data;
	
	fcomplex *d_data;

	size_d_data = 0;
		
	fftlen = subharminfs[0][0].kern[0].fftlen;
	
	size_d_data = fftlen ;

if (numharmstages > 1) {

	for(stage=1; stage<numharmstages; stage++){	
		harmtosum = 1 << stage;		
		for (harm = 1; harm < harmtosum; harm += 2) {               					
			
			fftlen = subharminfs[stage][harm-1].kern[0].fftlen;
			
			if( size_d_data < fftlen )
			{
				size_d_data = fftlen;
			}
			
		}	
	}

}
	
	//alloc memory for device data
	 checkCudaErrors(cudaMalloc((void **)&d_data, sizeof(fcomplex) * size_d_data));
	
	return d_data ;	
}                            

//----------Prepare d_result on GPU---------------------
fcomplex * prep_result_on_gpu(subharminfo **subharminfs, int numharmstages)
//cudaMalloc d_result for the whole search process
{

	int harm, harmtosum, stage;		
	int numkern, fftlen;
	int size_d_result;
	
	fcomplex *d_result;	
	size_d_result = 0;

	numkern = subharminfs[0][0].numkern;
	fftlen = subharminfs[0][0].kern[0].fftlen;

	size_d_result = numkern * fftlen;

if (numharmstages > 1) {

	for(stage=1; stage<numharmstages; stage++){	
		harmtosum = 1 << stage;		
		for (harm = 1; harm < harmtosum; harm += 2) {               		
						
			numkern = subharminfs[stage][harm-1].numkern;
			fftlen = subharminfs[stage][harm-1].kern[0].fftlen;
								
			if(size_d_result < numkern * fftlen )
			{
				size_d_result = numkern * fftlen;
			}			
		}	
	}	
}
	//Alloc mem for result on GPU  
  checkCudaErrors(cudaMalloc((void **)&d_result, size_d_result * sizeof(fcomplex)));  
	
	//also prepare its 2D texture memory
	channelDesc = cudaCreateChannelDesc(32, 32, 0, 0, cudaChannelFormatKindFloat);
	tex_d_result_2D.addressMode[0] = cudaAddressModeClamp;
  tex_d_result_2D.addressMode[1] = cudaAddressModeClamp;
  //texRef2D.filterMode = cudaFilterModePoint;
  tex_d_result_2D.filterMode = cudaFilterModePoint;
  tex_d_result_2D.normalized = false;
	
	return d_result ;
	
}                            

//----------Copy kernek array to GPU-------------------------------
fcomplex * cp_kernel_array_to_gpu(subharminfo **subharminfs, int numharmstages, int **offset_array)
//cp kernel arrays contained in **subharminfs to GPU memory *kernel_array_on_gpu
//and store the offset within *kernel_array_on_gpu for every array, offsets stored in **offset_array
//subharminfs: input
//numharmstages: input
//offset_array: output
//kernel_array_on_gpu: output
//
{

	int harm, harmtosum, stage;		
	int ii, jj;
	int kernel_total_size =0 ;
	int numkern, fftlen;
	unsigned int offset_base = 0, offset_tmp;
	fcomplex *kernel_vect_host;
	fcomplex *kernel_vect_on_gpu;
	
	//printf("cp_kernel_array_to_gpu: Check CUDA devices\n");
	
	numkern = subharminfs[0][0].numkern;
	fftlen = subharminfs[0][0].kern[0].fftlen;
	offset_array[0][0] = kernel_total_size ; 
	kernel_total_size = kernel_total_size + numkern * fftlen;
	
	for(stage=1; stage<numharmstages; stage++){	
		harmtosum = 1 << stage;		
		for (harm = 1; harm < harmtosum; harm += 2) {               					
			numkern = subharminfs[stage][harm-1].numkern;
			fftlen = subharminfs[stage][harm-1].kern[0].fftlen;
			offset_array[stage][harm-1] = kernel_total_size ; 
			kernel_total_size = kernel_total_size + numkern * fftlen;
		}	
	}	
	//printf("    ----Total size of kernel arrays is : %d\n", kernel_total_size);
	checkCudaErrors(cudaMalloc((void **)&kernel_vect_on_gpu, sizeof(fcomplex) * kernel_total_size));       
	
	kernel_vect_host = (fcomplex *)malloc( sizeof(fcomplex) * kernel_total_size);
	//take kernel_vect_host as buffer
	 //stage : 0
	 numkern = subharminfs[0][0].numkern;
   fftlen = subharminfs[0][0].kern[0].fftlen;
   offset_base = offset_array[0][0];
   for(ii=0; ii<numkern; ii++)
   	{   		
   		for(jj=0; jj<fftlen; jj++){   		
   			offset_tmp = offset_base + fftlen*ii + jj;
   			kernel_vect_host[offset_tmp].r = subharminfs[0][0].kern[ii].data[jj].r;
   			kernel_vect_host[offset_tmp].i = subharminfs[0][0].kern[ii].data[jj].i;
   		}
   	}

   //other stages
if (numharmstages > 1) {   
   for(stage=1; stage<numharmstages; stage++){	
			harmtosum = 1 << stage;		
			for(harm = 1; harm < harmtosum; harm += 2) {         

				offset_base = offset_array[stage][harm-1];

				fftlen = subharminfs[stage][harm-1].kern[0].fftlen;
				numkern = subharminfs[stage][harm-1].numkern;

				for(ii=0; ii<numkern; ii++){
					for(jj=0; jj<fftlen; jj++){
						offset_tmp = offset_base + fftlen*ii + jj ;
						kernel_vect_host[offset_tmp].r = subharminfs[stage][harm-1].kern[ii].data[jj].r;
		   			kernel_vect_host[offset_tmp].i = subharminfs[stage][harm-1].kern[ii].data[jj].i;
						}
					}
					
				}
			}
}   
   checkCudaErrors(cudaMemcpy( kernel_vect_on_gpu, kernel_vect_host, sizeof(fcomplex) * kernel_total_size, cudaMemcpyHostToDevice) );

   free(kernel_vect_host);
   
   return kernel_vect_on_gpu;
	
}

//----------------------initialize cuFFT plans ----------------------
void init_cuFFT_plans(subharminfo **subharminfs, int numharmstages)
{
	int harm, harmtosum, stage;		
	int numkern, fftlen;
	
	printf("\ninit_cuFFT_plans\n");
	
	numkern = subharminfs[0][0].numkern;
	fftlen = subharminfs[0][0].kern[0].fftlen;

	checkCudaErrors(cufftPlan1d(&plan_data_array[0][0], fftlen, CUFFT_C2C, 1));		
	checkCudaErrors(cufftPlan1d(&plan_result_array[0][0], fftlen, CUFFT_C2C, numkern));		

	printf("stage: 0, fftlen : %d, numkern: %d\n", fftlen, numkern);

if (numharmstages > 1) {

	for(stage=1; stage<numharmstages; stage++){	
		harmtosum = 1 << stage;		
		for (harm = 1; harm < harmtosum; harm += 2) {               					
			numkern = subharminfs[stage][harm-1].numkern;
			fftlen = subharminfs[stage][harm-1].kern[0].fftlen;
			
			printf("stage: %d, fftlen : %d, numkern: %d\n", stage, fftlen, numkern);
				
			checkCudaErrors(cufftPlan1d(&plan_data_array[stage][harm-1], fftlen, CUFFT_C2C, 1));		
			checkCudaErrors(cufftPlan1d(&plan_result_array[stage][harm-1], fftlen, CUFFT_C2C, numkern));		
						
		}	
	}	

}
	
}

//----------------------destroy cuFFT plans ----------------------
void destroy_cuFFT_plans(subharminfo **subharminfs, int numharmstages)
{
	int harm, harmtosum, stage;		
	int numkern, fftlen;
	
	printf("\ndestroy_cuFFT_plans\n");
	
	numkern = subharminfs[0][0].numkern;
	fftlen = subharminfs[0][0].kern[0].fftlen;

	checkCudaErrors(cufftDestroy(plan_data_array[0][0]));
	checkCudaErrors(cufftDestroy(plan_result_array[0][0]));

	printf("stage: 0, fftlen : %d, numkern: %d\n", fftlen, numkern);

if (numharmstages > 1) {

	for(stage=1; stage<numharmstages; stage++){	
		harmtosum = 1 << stage;		
		for (harm = 1; harm < harmtosum; harm += 2) {               					
			numkern = subharminfs[stage][harm-1].numkern;
			fftlen = subharminfs[stage][harm-1].kern[0].fftlen;
			
			printf("stage: %d, fftlen : %d, numkern: %d\n", stage, fftlen, numkern);
	
			checkCudaErrors(cufftDestroy(plan_data_array[stage][harm-1]));
			checkCudaErrors(cufftDestroy(plan_result_array[stage][harm-1]));
					
		}	
	}	
}
	
}


//----------------------select a cpu to play with --------------------
void select_cuda_dev(int cuda_inds)
{
	
	cudaSetDevice(cuda_inds);
  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, cuda_inds);

  printf("\nGPU Device %d: \"%s\" with Capability: %d.%d\n", cuda_inds, deviceProp.name, deviceProp.major, deviceProp.minor);
	
	cudaDeviceReset();
}

//-------------------------clean up on GPU ----------------------------
void cudaFree_kernel_vect(fcomplex *in)
{
	cudaFree( in );
}
