#include <iostream>
#include <random>
#include <vector>
#include "mpi.h"
#include <ctime>
#include <iomanip>


using namespace std;

// Size of the side of the grid
const int N=10;

// Function to initialize the grid
void initGrid(int,double,double**);
// Function for calculating the black dots (red-black algorithm) 
void blackDots(double**,double**,double ,int,int,int,int,int );
// Function for calculating the red dots ( red-black algorithm)
void redDots(double**,double**,double,int,int,int,int,int);
// Function for printing the grid
void printGrid(double**);
// Function for allocating 2d array
double **allocate( int,int);
// Function g
double g(int,int);
      


int main(int argc,char ** argv){

  cout << std::fixed<< std::setprecision(3);
  
  int id,ntasks;
  
  clock_t t1,t2;
  
  t1=clock();
  
  // Initializing MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&id);
  MPI_Comm_size(MPI_COMM_WORLD,&ntasks);
  MPI_Request send_request, recv_request;
  MPI_Status status;
  
  // Allocating the grid
  double **grid=0;
  grid=allocate(N,N);

  
  // Boundary condition
  double bound;
  bound = 3;
  
  initGrid(N,bound,grid);
  
  if(id==0){
    cout<<"Initial grid: "<<endl;
    printGrid(grid);

  }
  
  int iterations = 100;

  // Relaxation parameter
  double relPara = 0.5;
  
  int redOrBlack=0;

  
  int j=0;

  // Allocating the grid, where are the components of 'grid' from previous round
  double **oldGrid=0;
  oldGrid=allocate(N,N);
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      oldGrid[i][j]=grid[i][j];
    }
  }

  // How many rows every process calculates
  int rows= N/ntasks;
  // The highest row of the process 
  int highRow= rows*id;
  // The lowest row of the process
  int lowRow=highRow+rows-1; 
 

  // Solving the Poisson's equation using SOR algorithm and red-black algorithm
  for(int t=0;t<iterations;t++){
   
    redDots(grid,oldGrid,relPara, N,lowRow,highRow,id,ntasks);
    blackDots(grid,oldGrid,relPara, N,lowRow,highRow,id,ntasks);
 
    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
	oldGrid[i][j]=grid[i][j];
      }
    }

  }
  // All processes must be here to proceed
  MPI_Barrier(MPI_COMM_WORLD);

  // Array for receiving
  double **gridPart;
  gridPart=allocate(N,N);

  // All processes except 0 send their array to process 0
  if(id>0){
    if(id==ntasks-1){
      cout<<"The grid of the last process before sending it to process 0"<<endl;
      printGrid(grid);
    }
    MPI_Isend(&(grid[0][0]),N*N,MPI_DOUBLE,0,0,MPI_COMM_WORLD,&send_request);
     
    MPI_Wait(&send_request,&status);
    
  }
   
  MPI_Barrier(MPI_COMM_WORLD);

  // Process 0 receives every process's array and replaces the right part of it's array with the received array
  if(id==0){
 
    for(int idn=1;idn<ntasks;idn++){
      
      MPI_Irecv(&(gridPart[0][0]),N*N,MPI_DOUBLE,idn,0,MPI_COMM_WORLD,&recv_request);
   
     
      for(int j= rows*idn;j<rows*idn+rows;j++){
	for(int i=0;i<N;i++){
	  grid[j][i]=gridPart[j][i];
	   
	}
      }
    }
  }

        
     
  // Printing the solved grid
  if(id==0){
    
    cout<<"Solved grid: "<<endl;
    printGrid(grid);

  }
  MPI_Barrier(MPI_COMM_WORLD);
  t2=clock();
  if(id==0){
    t2=clock();
    cout<<"The CPU time used: "<<(double)(t2-t1)/CLOCKS_PER_SEC<<" seconds"<<endl;
  }
  MPI_Finalize();
  return 0;
}

// Initializing the grid and setting boundary conditions
void initGrid(int s,double bound,double ** grid){
  mt19937 rnd(453);
 
  for(int i=1;i<s-1;i++){
    for(int j=1;j<s-1;j++){
      grid[i][j]=(double)rnd()/rnd.max();
    }
  }
  
  for(int i=0;i<s;i++){
    grid[i][0]=bound;
    grid[i][s-1]=bound;
  }
  for(int i=1;i<s-1;i++){
    grid[0][i]=bound;
    grid[s-1][i]=bound;
  }
  
}

// Calculating the red dots
void redDots(double ** grid,double ** oldGrid,double relPara,int s,int lowRow, int highRow,int id,int ntasks){
  
  int j=highRow;

  // Array for receiving
  double *upper=new double[1];
  double *lower=new double[1];
  
  MPI_Request send_request, recv_request;
  MPI_Status status;

  // Setting the limits right
  int limit = lowRow;
  if(lowRow == s-1){
    limit = lowRow-1;
  }
  if(highRow == 0){
    highRow=1;
  }
  j=highRow;
  // Every process calculates different part of the array
  while(j<=limit){

    // Every other row
    if(j %2 == 0  && j != 0){
       
      for(int i=1;i<s-1;i+=2){
	if(j==lowRow && id != ntasks-1){
          // Processes send their lowest row coordinates to the process that has the lower part of the grid
	  MPI_Isend(&(grid[j][i]),1,MPI_DOUBLE,id+1,0,MPI_COMM_WORLD,&send_request);
	  MPI_Wait(&send_request,&status);
	}
	if(j== highRow && id != ntasks-1){
	  // Processes in their highest row receiving from above
	  MPI_Irecv(&(upper[0]),1,MPI_DOUBLE,id-1,0,MPI_COMM_WORLD,&recv_request);

	}
      
	if(j == highRow && id != 0){
	  // Processes sending their highest row coordinates to the process that has the upper part of the grid
	  MPI_Isend(&(oldGrid[j][i]),1,MPI_DOUBLE,id-1,0,MPI_COMM_WORLD,&send_request);
	  MPI_Wait(&send_request,&status);
	}
	if(j == lowRow && id != 0){
	  // Processes in their lowest row receiving from below
	  MPI_Irecv(&lower[0],1,MPI_DOUBLE,id+1,0,MPI_COMM_WORLD,&recv_request);

	}
	// Calculating the Poisson's equation
	if(j==lowRow & id != ntasks-1){
	  grid[j][i]=(1-relPara)*oldGrid[j][i]+(relPara/4)*(oldGrid[j+1][i]+grid[j-1][i]+oldGrid[j][i+1]+lower[0])-relPara/(4*s*s)*g(i,j);
	}else if( j== highRow && id != 0){
	  grid[j][i]=(1-relPara)*oldGrid[j][i]+(relPara/4)*(oldGrid[j+1][i]+grid[j-1][i]+upper[0]+grid[j][i-1])-relPara/(4*s*s)*g(i,j);
	}
	if( j != lowRow && j != highRow){

	  grid[j][i]=(1-relPara)*oldGrid[j][i]+(relPara/4)*(oldGrid[j+1][i]+grid[j-1][i]+oldGrid[j][i+1]+grid[j][i-1])-relPara/(4*s*s)*g(i,j);
	}
      }
   
    }else if(j != 0 && j%2!= 0){
      
      for(int i=2;i<s-1;i+=2){
	if(j==lowRow && id != ntasks-1){
          // Processes send their lowest row coordinates to the process that has the lower part of the grid
	  MPI_Isend(&(grid[j][i]),1,MPI_DOUBLE,id+1,0,MPI_COMM_WORLD,&send_request);
	  MPI_Wait(&send_request,&status);
	}
	if(j== highRow && id != ntasks-1){
	  // Processes in their highest row receiving from above
	  MPI_Irecv(&upper[0],1,MPI_DOUBLE,id-1,0,MPI_COMM_WORLD,&recv_request);
	}
      
	if(j == highRow && id != 0){
	  // Processes sending their highest row coordinates to the process that has the upper part of the grid
	  MPI_Isend(&(oldGrid[j][i]),1,MPI_DOUBLE,id-1,0,MPI_COMM_WORLD,&send_request);
	  MPI_Wait(&send_request,&status);
	}
	if(j == lowRow && id != 0){
	  // Processes in their lowest row receiving from below
	  MPI_Irecv(&lower[0],1,MPI_DOUBLE,id+1,0,MPI_COMM_WORLD,&recv_request);
	}
	// Calculating the Poisson's equation
	if(j==lowRow & id != ntasks-1){
	  grid[j][i]=(1-relPara)*oldGrid[j][i]+(relPara/4)*(oldGrid[j+1][i]+grid[j-1][i]+oldGrid[j][i+1]+lower[0])-relPara/(4*s*s)*g(i,j);
	}else if( j== highRow && id != 0){
	  grid[j][i]=(1-relPara)*oldGrid[j][i]+(relPara/4)*(oldGrid[j+1][i]+grid[j-1][i]+upper[0]+grid[j][i-1])-relPara/(4*s*s)*g(i,j);
	}
	if( j != lowRow && j != highRow){
   	  grid[j][i]=(1-relPara)*oldGrid[j][i]+(relPara/4)*(oldGrid[j+1][i]+grid[j-1][i]+oldGrid[j][i+1]+grid[j][i-1])-relPara/(4*s*s)*g(i,j);
	}
  
      }
     
    }
    j++;
   
  }
  
}



  
void blackDots(double **  grid, double ** oldGrid,double relPara,int s,int lowRow,int highRow,int id,int ntasks){
  
  int j=highRow;
  int row = 1;
  // Arrays for receiving
  double *upper = new double[1];
  double *lower = new double[1];

  MPI_Request send_request, recv_request;
  MPI_Status status;

  // Setting the limits right
  int limit = lowRow;
  if(lowRow == s-1){
    limit = lowRow-1;
  }
  if(highRow == 0){
    highRow=1;
  }
  j=highRow;

  
  while(j<=limit){
    
    
    if(j %2!=0 && j != 0){
    
      for(int i=1;i<s-1;i+=2){
	if(j==lowRow && id != ntasks-1){
	  // Processes send their lowest row coordinates to the process that has the lower part of the grid
	  MPI_Isend(&(grid[j][i]),1,MPI_DOUBLE,id+1,0,MPI_COMM_WORLD,&send_request);
	  MPI_Wait(&send_request,&status);
	}
	if(j== highRow && id != ntasks-1){
	  // Processes in their highest row receiving from above
	  MPI_Irecv(&upper[0],1,MPI_DOUBLE,id-1,0,MPI_COMM_WORLD,&recv_request);

	}
      
	if(j == highRow && id != 0){
	  // Processes sending their highest row coordinates to the process that has the upper part of the grid
	  MPI_Isend(&(oldGrid[j][i]),1,MPI_DOUBLE,id-1,0,MPI_COMM_WORLD,&send_request);
	  MPI_Wait(&send_request,&status);
	}
	if(j == lowRow && id != 0){
	  // Processes in their lowest row receiving from below
	  MPI_Irecv(&lower[0],1,MPI_DOUBLE,id+1,0,MPI_COMM_WORLD,&recv_request);

	}
	// Calculating the Poisson's equation
	if(j==lowRow & id != ntasks-1){
	  grid[j][i]=(1-relPara)*oldGrid[j][i]+(relPara/4)*(oldGrid[j+1][i]+grid[j-1][i]+oldGrid[j][i+1]+lower[0])-relPara/(4*s*s)*g(i,j);
	}else if( j== highRow && id != 0){
	  grid[j][i]=(1-relPara)*oldGrid[j][i]+(relPara/4)*(oldGrid[j+1][i]+grid[j-1][i]+upper[0]+grid[j][i-1])-relPara/(4*s*s)*g(i,j);
	}
	if( j != lowRow && j != highRow){

	  
	  grid[j][i]=(1-relPara)*oldGrid[j][i]+(relPara/4)*(oldGrid[j+1][i]+grid[j-1][i]+oldGrid[j][i+1]+grid[j][i-1])-relPara/(4*s*s)*g(i,j);
	}
      }
      
    }else if( j%2==0 && j!= 0){
   
      for(int i=2;i<s-1;i+=2){
      
	if(j==lowRow && id != ntasks-1){
	  // Processes send their lowest row coordinates to the process that has the lower part of the grid
	  MPI_Isend(&(grid[j][i]),1,MPI_DOUBLE,id+1,0,MPI_COMM_WORLD,&send_request);
	  MPI_Wait(&send_request,&status);
	}
	if(j== highRow && id != ntasks-1){
     	  // Processes in their highest row receiving from above
	  MPI_Irecv(&upper[0],1,MPI_DOUBLE,id-1,0,MPI_COMM_WORLD,&recv_request);

	}
      
	if(j == highRow && id != 0){
	  // Processes sending their highest row coordinates to the process that has the upper part of the grid
	  MPI_Isend(&oldGrid[j][i],1,MPI_DOUBLE,id-1,0,MPI_COMM_WORLD,&send_request);
	  MPI_Wait(&send_request,&status);
	}
	if(j == lowRow && id != 0){
	  // Processes in their lowest row receiving from below

	  MPI_Irecv(&lower[0],1,MPI_DOUBLE,id+1,0,MPI_COMM_WORLD,&recv_request);

	}
	// Calculating the Poisson's equation
	if(j==lowRow & id != ntasks-1){
	  grid[j][i]=(1-relPara)*oldGrid[j][i]+(relPara/4)*(oldGrid[j+1][i]+grid[j-1][i]+oldGrid[j][i+1]+lower[0])-relPara/(4*s*s)*g(i,j);
	}else if( j== highRow && id != 0){
	  grid[j][i]=(1-relPara)*oldGrid[j][i]+(relPara/4)*(oldGrid[j+1][i]+grid[j-1][i]+upper[0]+grid[j][i-1])-relPara/(4*s*s)*g(i,j);
	}
	if( j != lowRow && j != highRow){

	  
	  grid[j][i]=(1-relPara)*oldGrid[j][i]+(relPara/4)*(oldGrid[j+1][i]+grid[j-1][i]+oldGrid[j][i+1]+grid[j][i-1])-relPara/(4*s*s)*g(i,j);
	}
      
      }
    }
   
     
  
    j++;
  }

  
}


void printGrid(double ** grid){
  
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      cout<<grid[i][j]<<"   ";
    }
    cout<<"\n";
  }

}

double **allocate( int Rows, int Cols){
  double **Array;

  Array = new double*[Rows];
  for( int i = 0 ; i < Rows ; i++ ){
    Array[i] = new double [Cols];

  }
  return Array;
}

double g(int i,int j){
  double value = i/N+j/N;
  return value;
}

