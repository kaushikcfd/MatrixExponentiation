#include <stdio.h>
#include <stdio.h>



__global__ void matrixExponentiation(float *A,float *C, float *D, float *E, float *Coeff, unsigned *N, unsigned *rowsPerThread)
{
	unsigned i,j,k,row,order;

	/*Initializing by the identity matrix. And initialized the matrix C by matrix A.*/
	for(i=0;i<(*rowsPerThread);i++)
	{
		row	=	blockIdx.x*blockDim.x	+	threadIdx.x*(*rowsPerThread)	+	i;
		for(j=0;j<(*N);j++)
		{
			E[row*(*N)+j]	=	Coeff[1]*A[row*(*N)+j];
			D[row*(*N)+j]	=	A[row*(*N)+j];
		}
		E[row*(*N)+row]	+=	Coeff[0];
	}
	__syncthreads();

	for(order=2;order<12;order++)
	{

		for(i=0;i<(*rowsPerThread);i++)
		{
			row	=	blockIdx.x*blockDim.x	+	threadIdx.x*(*rowsPerThread)	+	i;
			for(j=0;j<(*N);j++)
			{
				C[row*(*N)+j]	=	0.0;
				for(k=0;k<(*N);k++)
					C[row*(*N)+j]+=(D[row*(*N)+k]*A[k*(*N)+j]);
			}
		}

		__syncthreads();
	
	
		for(i=0;i<(*rowsPerThread);i++)
		{
			row	=	blockIdx.x*blockDim.x	+	threadIdx.x*(*rowsPerThread)	+	i;
			for(j=0;j<(*N);j++)
			{
				E[row*(*N)+j]+=(Coeff[order]*C[row*(*N)+j]);
				D[row*(*N)+j] =(C[row*(*N)+j]);	
			}
		}
		__syncthreads();
	}
	
}

void makeA(float *A,unsigned N)
{
	unsigned i,j;
	for(i=0;i<N;i++)
		for(j=0;j<N;j++)
			A[i*N+j]=1e-3;	
	return ;
}

void printMatrix(float *A, unsigned N)
{
	unsigned i,j;
	for(i=0;i<N;i++)
	{
		for(j=0;j<N;j++)
			printf("%6.4f\t",A[i*N+j]);
		printf("\n");
	}
	return ;
}

int main()
{
	unsigned N,blocks, threads,i;
	unsigned *dev_N;
	float *A,*Exp,*Coeff;
	float *dev_A,*dev_B,*dev_C,*dev_Exp,*dev_Coeff;
	unsigned size;
	unsigned rowsPerThread, *dev_rowsPerThread;	

	printf("The order of matrix to be used\n");
	scanf("%d",&N);	
	printf("Enter the number of blocks.\n");
	scanf("%d",&blocks);
	printf("Enter the number of threads per block.\n");
	scanf("%d",&threads);
	
	if((N%(threads*blocks))!=0)
	{
		printf("The order of the matrix `N` must be divisible by the product `threads*blocks`\n. Aborting the program!\n");
	}
	else
	{
		cudaEvent_t	start, stop;
		cudaEventCreate(&start);
		cudaEventCreate(&stop);
		size		=	N*N*sizeof(float);	
		rowsPerThread	=	(N/(threads*blocks));
				
		A		=	(float *)malloc(size);
		Exp		=	(float *)malloc(size);
		Coeff		=	(float *)malloc(12*sizeof(float));
		Coeff[0]	=	1.0;
		for(i=1;i<12;i++)
			Coeff[i]	=	(Coeff[i-1]/(1.0*i));
		cudaEventRecord(start);
		cudaMalloc((void**)&dev_A,size);
		cudaMalloc((void**)&dev_B,size);
		cudaMalloc((void**)&dev_C,size);
		cudaMalloc((void**)&dev_Exp,size);
		cudaMalloc((void**)&dev_Coeff,12*sizeof(float));
		cudaMalloc((void**)&dev_rowsPerThread,sizeof(unsigned));
		cudaMalloc((void**)&dev_N,sizeof(unsigned));

		makeA(A,N);

		cudaMemcpy(dev_A,	A,	size,cudaMemcpyHostToDevice);	
		cudaMemcpy(dev_Coeff,	Coeff,	12*sizeof(float),cudaMemcpyHostToDevice);
		cudaMemcpy(dev_rowsPerThread,&rowsPerThread,sizeof(unsigned),cudaMemcpyHostToDevice);
		cudaMemcpy(dev_N,&N,sizeof(unsigned),cudaMemcpyHostToDevice);		

		matrixExponentiation<<<blocks,threads>>>(dev_A, dev_B, dev_C, dev_Exp, dev_Coeff ,dev_N,dev_rowsPerThread);		

		cudaMemcpy(Exp,dev_Exp,size,cudaMemcpyDeviceToHost);

		cudaFree(dev_A);
		cudaFree(dev_B);
		cudaFree(dev_C);
		cudaFree(dev_rowsPerThread);
		cudaFree(dev_N);

		cudaEventRecord(stop);
		cudaEventSynchronize(stop);
		float	milliseconds	=	0.0;
		cudaEventElapsedTime(&milliseconds,start,stop);
		fprintf (stderr,"Time for the Matrix Multiplication of order %d : %f s using blocks %d and threads per block %d.\n\n",N ,0.001*milliseconds,blocks,threads);
		/*freopen("A.dat","w",stdout);
		printMatrix(A,N);
		fclose(stdout);
		freopen("B.dat","w",stdout);
		printMatrix(B,N);
		fclose(stdout);*/
		freopen("C.dat","w",stdout);
		printMatrix(Exp,N);
		fclose(stdout);
	}
	return 0;
}

