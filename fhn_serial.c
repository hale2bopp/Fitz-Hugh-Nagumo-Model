#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>


float f;
float C;
int N;

int stopper;
float vth = 0.3;
float vsl = 0.001;
int A = 3;
int B = 3;
float e = 0.01;
//float e = 10;
float c = -0.1;
float gsyn = 0.05;
float h = 0.001;

int* mat = NULL;
float* I=NULL;
float* E=NULL;
float** g = NULL;


void array_add(int len_array_in,float array_in[N*(N+1)], float array_add[N*(N+1)],float array_out_add[N*(N+1)]);
void integrator_rk4(float dt,float t, float* p1,  float yout[N*(N+1)]);
float ff(float v, float w);
float gf(float v, float w);
float Nf(float v);
void isin(int* array, int len_arr, int array_out[N]);
void oscnetwork_opt(float t, float y[N*(N+1)],float dydt[N*(N+1)]);
void array_mul(int len_array_in,float array_in[N*(N+1)], float num,float array_out[N*(N+1)]);


int main(int argc, char *argv[]){

	N = atoi(argv[1]);
	f = atof(argv[2]);
	C = atof(argv[3]);

	
	srand(1);
	/* initializations*/

	int i,j,iter;
	
	
	/* allocating space for all the global arrays of variable size*/
//	stimulation array
	if (I != 0) {
	    I = (float*) realloc(I, N * sizeof(float));
	} else {
	    I = (float*) malloc(N * sizeof(float));
	}
//	matrix giving the locations of voltage in the overall array	
	if (mat != 0) {
	    mat = (int*) realloc(mat, N * sizeof(int));
	} else {
	    mat = (int*) malloc(N * sizeof(int));
	}

// 	reversal potential array ( -ve  : inhibitory, +ve : excitatory)
	if (E != 0) {
	    E = (float*) realloc(E, N * sizeof(float));
	} else {
	    E = (float*) malloc(N * sizeof(float));
	}

//	weights for the adjacency matrix : g synaptic
	g = malloc(N * sizeof(float *));
	for(i = 0; i < N; i++){
		g[i] = malloc(N * sizeof(float));
		if(g[i] == NULL){
	        	fprintf(stderr, "out of memory\n");
        		break;
		}
	}
	
// 	fraction of neurons that are stimulated
	float stim = 0.5;
	
	int sum_adj;
//	time span for which to run simulation
	int tspan = 100;
	
//	initialize a pointer		
	float *q;	

	
	float randomnum;
	
// 	total number of time iterations = tspan*step_size	
	int tot_time = (int) ceil(tspan/h);
	
// 	initialize the fraction variables
	float count_exc, count_inh;

//	Time array
	float T[tot_time];


// 	Time variable	
	float t;


//	Initialize the stimulation variable with magnitude
	float stimag = 1.0;	 
	for (i = 0;i<(int) floor(stim*N);i++){
		I[i]= stimag;	

	}
	
	for (i = (int) floor(stim*N);i<N;i++){
		I[i]= 0; 	

	}

//	printf("rand %d %d ", rand(), RAND_MAX);
//	randomize the stimulation array	
	for (i = 0; i < N; i++) 
        {
		j = i + rand() / ((RAND_MAX / (N - i)) + 1);
		t = I[j];
		//	  printf(" j: I[%d] = %f  ",j,t);	
		I[j] = I[i];
		I[i] = t;
		//	  printf("j: %d,I[%d]= %f ",j,i,I[i]);
}
	printf("\n");

// 	initialize the reversal potential array
	for (i = 0;i<(int) floor(f*N);i++){
		E[i]= -5;	

	}
	
	for (i =(int) floor(f*N);i<N;i++){
		E[i]= 5;	


	}

// 	randomize the reversal potential array
	for (i = 0; i < N; i++) 
        {
		j = i + rand() / ((RAND_MAX / (N - i)) + 1);
		t = E[j];
		E[j] = E[i];
		E[i] = t;
        }
// 	initialize the adjacency matrix
	for (i=0;i<N;i++){	
		sum_adj = 0;
		for(j=0;j<N;j++){
			randomnum = (float) rand()/(float) RAND_MAX;	
			if (randomnum<C){
				g[i][j] = gsyn;
				sum_adj += 1;
			}
			else{ 	
				g[i][j] = 0;
			}
			if (i==j){
				g[i][j] = 0;
			} 	
		}
		for (j=0;j<N;j++){
			g[i][j] = g[i][j]/(sum_adj+1);
//			printf("%f ",g[i][j]);
		}
//		printf("\n");
	}
//	printf("\n\n\n");
		
		
// 	vector to hold values for each differential variable for all time iterations
	float Y[N*(N+1)];

//	initialize the vector of indices showing the locations of the voltage variables
//	The first location is 0 because the first variable in the differential vector will
//	always be v1, having the index 0.
	mat[0] = 0;	
//	printf("mat[0] = 0  ");
//	To get to the location of the next voltage variable, add 2 for vi and wi, 
//	add the sum of the adjacency matrix for the previous row	
	for (i=1;i<N;i++){	
		mat[i] = mat[i-1] + N+1;	
	}
	
// 	initial conditions vector for time = 0
	for (i = 0;i<N*(N+1);i++)
		Y[i] = 0;
	
//	set the time array
	T[0] = 0;

//	flag to stop current after a certain time	
	stopper = 1;
	
//	receive the result into this array
	float Y1[N*(N+1)];

// 	This loop calls the RK4 code
	for (i=0;i<tot_time-1;i++){
		q = Y;
//		call the RK4 integrator with current time value, and current 
//		values of voltage			
		integrator_rk4(h,T[i],Y,Y1);  
		
//		Return the time output of integrator into the next iteration of time
		T[i+1] = T[i]+h;	
		
//		copy the output of the integrator into the next iteration of voltage		
		q = memcpy(q, Y1, N*(N+1) * sizeof(float));
//		for (i=0;i<N*(N+1);i++)
//			Y[i] = t_y.y[i];

		printf("%f ",T[i+1]);
		count_exc = 0;
		count_inh = 0;
		for (iter = 0;iter<N;iter++){
//			printf("%f ", t_y.y[mat[iter]]);//*(p+mat[iter]));
			printf("%f ",Y[mat[iter]]); 
			
			if ((Y[mat[iter]]>(0.1)) && (E[iter] > 0 )){
				count_exc +=1;

			}
			else if ((Y[mat[iter]]>(0.1)) && (E[iter] < 0 )){
				count_inh +=1;
			}	
		}
                printf(" %f %f", ((float) count_exc)/(((float) N)*(1-f)), ((float) count_inh)/((float )(N)*(f))); 

//		printf(" %f %f", count_exc/(N*(1-f)),count_inh/(N*(f)));		
//		printf(" %f %f", count_exc,count_inh);		
		printf("\n");
	}
// 	show the actual value of the Connection density
//	printf("\nactual C: %f %d",(float) sum_adj/(N*N),N*(N+1));

// 	free all the memory that was allocated for the arrays		
	free(E);
	free(I);

	for(i = 0; i < N; i++)
    		free(g[i]);
	free(g);
	free(mat);

//	free(t_y.y);



	return 0;
//	exit(0);
}

void integrator_rk4(float dt,float t, float y[N*(N+1)], float yout[N*(N+1)])
{	
//	initialize all the pointers
	float y1[N*(N+1)],y2[N*(N+1)],y3[N*(N+1)];
	float tout,dt_half;
	float k1[N*(N+1)],k2[N*(N+1)],k3[N*(N+1)],k4[N*(N+1)];
// 	initialize iterator
	int i;

	tout = t+dt;
	dt_half = 0.5*dt;
	float addition[N*(N+1)];

//	return the differential array into k1
	oscnetwork_opt(t,y,k1);
//	#pragma omp parallel for 
// 	multiply the array k1 by dt_half
	for(i=0;i<N*(N+1);i++)
		y1[i]=y[i]+(k1[i])*dt_half;	


// 	do the same thing 3 times
	oscnetwork_opt(t+dt_half,y1,k2);
//	#pragma omp parallel for 
	for(i=0;i<N*(N+1);i++)
		y2[i]=y[i]+(k2[i])*dt_half;	
	
	oscnetwork_opt(t+dt_half,y2,k3);
//	#pragma omp parallel for 
	for(i=0;i<N*(N+1);i++)
		y3[i]=y[i]+(k3[i])*dt;	

	oscnetwork_opt(tout,y3,k4);
//	Make the final additions with k1,k2,k3 and k4 according to the RK4 code
//	#pragma omp parallel for 
	for (i=0;i<N*(N+1);i++){
		addition[i] = ((k1[i]) + (k2[i])*2 + (k3[i])*2 + (k4[i])) *dt/6;
	}
//	add this to the original array
//	#pragma omp parallel for 
	for(i=0;i<N*(N+1);i++)
		yout[i]=y[i]+addition[i];		
}

 
// function to set the voltage variable v
float ff(float v, float w){
	float ret_val;
	ret_val =  v*(v+c)*(1-v)-w;
	return ret_val;
}

// function to set the recovery variable w
float gf(float v, float w){
	float ret_val;
	ret_val =  v-0.5*w;
	return ret_val;	
}

// function to set the synaptic gate fraction variables sji
float Nf(float v){
	float ret_val;
	ret_val = 0.5*( (1 + atan((v-vth)/vsl) ) );
	return ret_val;
}

int sign(float t){
	int flag;
	if (t<0)
		flag = -1;
	else if (t>0)
		flag = 1;
	else
		flag =  0;
	return flag;
}	



// function to return the vector with coupled differential variables for each time iteration
void oscnetwork_opt(float t, float y[N*(N+1)],float dydt[N*(N+1)]){
//	initialize iterators	
	int i,j;


	float sum_syn[N];


	int count1[N], count2[N];



//	this little piece of code stops the stimulation after time = 0.05 seconds
// 	if stopper is set to 0 inside the loop

	if (stopper==1)
		if (t>7){
			stopper = 0;
		}


//	d_iter over each value in the differential vector: 
//	d_iter : 0...2*N + sum_adj -1
	

// 	i is reserved to iterate over the indices of voltage as mat[i].
//	i : 0...N-1

	for (i=0;i<N;i++){
		
//		find the sum of contributions from different neurons		
		sum_syn[i] = 0;

//		count : 0... sum_column
		
		count1[i] = 0;

//		only consider those columns which have non-zero values for the ith row		
		for (j = 0;j<N;j++){
			if (j!=i){
				count1[i]+=1;
				sum_syn[i] += g[i][j]*(y[mat[i]]-E[j])*y[mat[i]+count1[i]+1];		
			}	
//			printf(" i: %d,   i*(N+1): %f ,   i*(N+1)+1 :%F, sum_syn[i]: %f , count1[i] : %d \n",i,dydt[i*(N+1)], dydt[i*(N+1)+1], sum_syn[i], count1[i] );
//                        printf("i : %d, j = %d, count1[i] : %d, sym_syn[i]: %f\n",i,j,count1[i],sum_syn[i]);
		}
		
		
//		set up the differential vector for voltage variable v
//		dydt[d_iter] = ( ff(y[mat[i]],y[mat[i]+1]) + (stopper*I[i])-sum_syn)/e;
		dydt[i*(N+1)] = ( ff(y[mat[i]],y[mat[i]+1]) + ((0.5+0.5*sign(t-1))*I[i])-sum_syn[i])/e;
// 		recovery variable w
		dydt[i*(N+1)+1] = gf(y[mat[i]],y[mat[i]+1]);
//		printf(" i: %d,   i*(N+1): %f ,   i*(N+1)+1 :%F, sum_syn[i]: %f , count1[i] : %d \n",i,dydt[i*(N+1)], dydt[i*(N+1)+1], sum_syn[i], count1[i] );
		count2[i] = 0;
//		printf(" d_iter: %d, d_iter+1 : %d ", d_iter, d_iter+1);
//		for each of the synaptic variables
		for (j = 0;j<N;j++){
			if (j!=i){
				count2[i] += 1;
				dydt[i*(N+1)+1+count2[i]] = A*Nf(y[mat[j]])*(1-y[mat[i]+count2[i]+1])-B*y[mat[i]+count2[i]+1];
//				printf("d_iter+1+count: %d ", d_iter+1+count );
			}
		}

//		update the variable by 2 to get to the next voltage variable
//		d_iter += count+2;
}
}
