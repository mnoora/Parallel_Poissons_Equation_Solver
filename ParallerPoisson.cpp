#include <iostream>
#include <random>
#include <vector>
#include "mpi.h"


using namespace std;

const int N=10;
void initGrid(int,double,double**);
void blackDots(double**,double**,double,double ,int,int,int,int,int );
void redDots(double**,double**,double,double,int,int,int,int,int);
void printGrid(double**);

double **allocate( int nRows, int nCols){
      double **dynamicArray;

      dynamicArray = new double*[nRows];
      for( int i = 0 ; i < nRows ; i++ )
      dynamicArray[i] = new double [nCols];

      return dynamicArray;
}



int main(int argc,char ** argv){
  int id,ntasks;

// Initializing MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&id);
  MPI_Comm_size(MPI_COMM_WORLD,&ntasks);
  MPI_Request send_request, recv_request;
  MPI_Status status;
  
  
  double **grid=0;
  grid=allocate(N,N);

  
 
  double bound;
  bound = 3;
  initGrid(N,bound,grid);
  if(id==0){
      cout<<"Initial grid: "<<endl;
     printGrid(grid);

  }
  int iterations = 400000;
  double relPara = 0.5;
  int redOrBlack=0;
  double g = 0;
  int j=0;
  double **oldGrid=0;
  oldGrid=allocate(N,N);
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      oldGrid[i][j]=grid[i][j];
    }
  }

 int rows= N/ntasks;
  int highRow= rows*id;
  int lowRow=highRow+rows-1; 
 
 
  for(int t=0;t<iterations;t++){

    redDots(grid,oldGrid,relPara,g, N,lowRow,highRow,id,ntasks);
    blackDots(grid,oldGrid,relPara,g, N,lowRow,highRow,id,ntasks);
  
 
  
    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
	oldGrid[i][j]=grid[i][j];
      }
    }

  }
  MPI_Barrier(MPI_COMM_WORLD);
  
  double **gridPart;
  gridPart=allocate(N,N);
   if(id>0){
 
      //MPI_Sendrecv(grid,N*N,MPI_DOUBLE,0,0,gridPart,N*N,MPI_INT,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
     if(id==4){
       cout<<"id="<<id<<" count="<<N*N<<endl;
       cout<<"lähetettävä:"<<endl;
       printGrid(grid);
     }
     MPI_Isend(&(grid[0][0]),N*N,MPI_DOUBLE,0,0,MPI_COMM_WORLD,&send_request);
     
     MPI_Wait(&send_request,&status);
    
    }
   MPI_Barrier(MPI_COMM_WORLD);
   if(id==0){
 
     for(int idn=1;idn<ntasks;idn++){
       //printGrid(gridPart);
       if(idn==4){
       cout<<"idn="<<idn<<" count="<<N*N<<endl;
       }
       MPI_Irecv(&(gridPart[0][0]),N*N,MPI_DOUBLE,idn,0,MPI_COMM_WORLD,&recv_request);
       if(idn==4){
       cout<<"Vastaanotto"<<endl;
       printGrid(gridPart);
       }
       cout<<"\n";
       for(int j= rows*idn;j<rows*idn+rows;j++){
	 for(int i=0;i<N;i++){
	   grid[j][i]=gridPart[j][i];
	   
	 }
       }
     }
   }

        
     
    
  if(id==0){
    
      cout<<"Solved grid: "<<endl;
     printGrid(grid);

  }
  MPI_Finalize();
  return 0;
}


void initGrid(int s,double bound,double ** grid){
  mt19937 rnd(453);
 
  for(int i=1;i<s-1;i++){
    for(int j=1;j<s-1;j++){
      grid[i][j]=0;//(double)rnd()/rnd.max();
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

void redDots(double ** grid,double ** oldGrid,double relPara, double g,int s,int lowRow, int highRow,int id,int ntasks){
  int j=highRow;
  int row = 1;
  double *upper=new double[1];
  double *lower=new double[1];
  MPI_Request send_request, recv_request;
  MPI_Status status;
  int limit = lowRow;
  if(lowRow == s-1){
    limit = lowRow-1;
  }
   if(highRow == 0){
    highRow=1;
  }
   j=highRow;
  while(j<=limit){
    
    if(j %2 == 0  && j != 0){
       
      for(int i=1;i<s-1;i+=2){
	if(j==lowRow && id != ntasks-1){
	  //MPI_Sendrecv(&oldGrid[i][j],1,MPI_DOUBLE,id-1,0,&upper[0],1,MPI_DOUBLE,id+1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	  MPI_Isend(&(grid[j][i]),1,MPI_DOUBLE,id+1,0,MPI_COMM_WORLD,&send_request);
	  MPI_Wait(&send_request,&status);
	}
	if(j== highRow && id != ntasks-1){
	  MPI_Irecv(&(upper[0]),1,MPI_DOUBLE,id-1,0,MPI_COMM_WORLD,&recv_request);

	}
      
	if(j == highRow && id != 0){
	  //  MPI_Sendrecv(&grid[i][j],1,MPI_DOUBLE,id+1,0,&lower[0],1,MPI_DOUBLE,id-1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	  MPI_Isend(&(oldGrid[j][i]),1,MPI_DOUBLE,id-1,0,MPI_COMM_WORLD,&send_request);
	  MPI_Wait(&send_request,&status);
	}
	if(j == lowRow && id != 0){
	  MPI_Irecv(&lower[0],1,MPI_DOUBLE,id+1,0,MPI_COMM_WORLD,&recv_request);

	  }
	if(j==lowRow & id != ntasks-1){
	  grid[j][i]=(1-relPara)*oldGrid[j][i]+(relPara/4)*(oldGrid[j+1][i]+grid[j-1][i]+oldGrid[j][i+1]+lower[0])-relPara/(4*s*s)*g;
	}else if( j== highRow && id != 0){
	  grid[j][i]=(1-relPara)*oldGrid[j][i]+(relPara/4)*(oldGrid[j+1][i]+grid[j-1][i]+upper[0]+grid[j][i-1])-relPara/(4*s*s)*g;
	}
	if( j != lowRow && j != highRow){

	  
	  grid[j][i]=(1-relPara)*oldGrid[j][i]+(relPara/4)*(oldGrid[j+1][i]+grid[j-1][i]+oldGrid[j][i+1]+grid[j][i-1])-relPara/(4*s*s)*g;
	}
      }
      row=0;
    }else if(j != 0 && j%2!= 0){
      for(int i=2;i<s-1;i+=2){
	if(j==lowRow && id != ntasks-1){
	  //MPI_Sendrecv(&oldGrid[i][j],1,MPI_DOUBLE,id-1,0,&upper[0],1,MPI_DOUBLE,id+1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	  MPI_Isend(&(grid[j][i]),1,MPI_DOUBLE,id+1,0,MPI_COMM_WORLD,&send_request);
	  MPI_Wait(&send_request,&status);
	}
	if(j== highRow && id != ntasks-1){
	  MPI_Irecv(&upper[0],1,MPI_DOUBLE,id-1,0,MPI_COMM_WORLD,&recv_request);

	}
      
	if(j == highRow && id != 0){
	  //  MPI_Sendrecv(&grid[i][j],1,MPI_DOUBLE,id+1,0,&lower[0],1,MPI_DOUBLE,id-1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	  MPI_Isend(&(oldGrid[j][i]),1,MPI_DOUBLE,id-1,0,MPI_COMM_WORLD,&send_request);
	  MPI_Wait(&send_request,&status);
	}
	if(j == lowRow && id != 0){
	  MPI_Irecv(&lower[0],1,MPI_DOUBLE,id+1,0,MPI_COMM_WORLD,&recv_request);

	  }
	if(j==lowRow & id != ntasks-1){
	  grid[j][i]=(1-relPara)*oldGrid[j][i]+(relPara/4)*(oldGrid[j+1][i]+grid[j-1][i]+oldGrid[j][i+1]+lower[0])-relPara/(4*s*s)*g;
	}else if( j== highRow && id != 0){
	  grid[j][i]=(1-relPara)*oldGrid[j][i]+(relPara/4)*(oldGrid[j+1][i]+grid[j-1][i]+upper[0]+grid[j][i-1])-relPara/(4*s*s)*g;
	}
	if( j != lowRow && j != highRow){

	  
	  grid[j][i]=(1-relPara)*oldGrid[j][i]+(relPara/4)*(oldGrid[j+1][i]+grid[j-1][i]+oldGrid[j][i+1]+grid[j][i-1])-relPara/(4*s*s)*g;
	}
      row=1;
      }
     
    }
  j++;
   
  }
  
}



  
void blackDots(double **  grid, double ** oldGrid,double relPara, double g,int s,int lowRow,int highRow,int id,int ntasks){
  int j=highRow;
  int row = 1;
  double *upper = new double[1];
  double *lower = new double[1];

  MPI_Request send_request, recv_request;
  MPI_Status status;
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
	  //MPI_Sendrecv(&oldGrid[i][j],1,MPI_DOUBLE,id-1,0,&upper[0],1,MPI_DOUBLE,id+1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	  MPI_Isend(&(grid[j][i]),1,MPI_DOUBLE,id+1,0,MPI_COMM_WORLD,&send_request);
	  MPI_Wait(&send_request,&status);
	}
	if(j== highRow && id != ntasks-1){
	  MPI_Irecv(&upper[0],1,MPI_DOUBLE,id-1,0,MPI_COMM_WORLD,&recv_request);

	}
      
	if(j == highRow && id != 0){
	  //  MPI_Sendrecv(&grid[i][j],1,MPI_DOUBLE,id+1,0,&lower[0],1,MPI_DOUBLE,id-1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	  MPI_Isend(&(oldGrid[j][i]),1,MPI_DOUBLE,id-1,0,MPI_COMM_WORLD,&send_request);
	  MPI_Wait(&send_request,&status);
	}
	if(j == lowRow && id != 0){
	  MPI_Irecv(&lower[0],1,MPI_DOUBLE,id+1,0,MPI_COMM_WORLD,&recv_request);

	  }
	if(j==lowRow & id != ntasks-1){
	  grid[j][i]=(1-relPara)*oldGrid[j][i]+(relPara/4)*(oldGrid[j+1][i]+grid[j-1][i]+oldGrid[j][i+1]+lower[0])-relPara/(4*s*s)*g;
	}else if( j== highRow && id != 0){
	  grid[j][i]=(1-relPara)*oldGrid[j][i]+(relPara/4)*(oldGrid[j+1][i]+grid[j-1][i]+upper[0]+grid[j][i-1])-relPara/(4*s*s)*g;
	}
	if( j != lowRow && j != highRow){

	  
	  grid[j][i]=(1-relPara)*oldGrid[j][i]+(relPara/4)*(oldGrid[j+1][i]+grid[j-1][i]+oldGrid[j][i+1]+grid[j][i-1])-relPara/(4*s*s)*g;
	}
      }
      
    }else if( j%2==0 && j!= 0){
   
      for(int i=2;i<s-1;i+=2){
      
      if(j==lowRow && id != ntasks-1){
	  //MPI_Sendrecv(&oldGrid[i][j],1,MPI_DOUBLE,id-1,0,&upper[0],1,MPI_DOUBLE,id+1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	MPI_Isend(&(grid[j][i]),1,MPI_DOUBLE,id+1,0,MPI_COMM_WORLD,&send_request);
	  MPI_Wait(&send_request,&status);
	}
	if(j== highRow && id != ntasks-1){
	  MPI_Irecv(&upper[0],1,MPI_DOUBLE,id-1,0,MPI_COMM_WORLD,&recv_request);

	}
      
	if(j == highRow && id != 0){
	  //  MPI_Sendrecv(&grid[i][j],1,MPI_DOUBLE,id+1,0,&lower[0],1,MPI_DOUBLE,id-1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	   MPI_Isend(&oldGrid[j][i],1,MPI_DOUBLE,id-1,0,MPI_COMM_WORLD,&send_request);
	  MPI_Wait(&send_request,&status);
	}
	if(j == lowRow && id != 0){
	  MPI_Irecv(&lower[0],1,MPI_DOUBLE,id+1,0,MPI_COMM_WORLD,&recv_request);

	  }
	if(j==lowRow & id != ntasks-1){
	  grid[j][i]=(1-relPara)*oldGrid[j][i]+(relPara/4)*(oldGrid[j+1][i]+grid[j-1][i]+oldGrid[j][i+1]+lower[0])-relPara/(4*s*s)*g;
	}else if( j== highRow && id != 0){
	  grid[j][i]=(1-relPara)*oldGrid[j][i]+(relPara/4)*(oldGrid[j+1][i]+grid[j-1][i]+upper[0]+grid[j][i-1])-relPara/(4*s*s)*g;
	}
	if( j != lowRow && j != highRow){

	  
	  grid[j][i]=(1-relPara)*oldGrid[j][i]+(relPara/4)*(oldGrid[j+1][i]+grid[j-1][i]+oldGrid[j][i+1]+grid[j][i-1])-relPara/(4*s*s)*g;
	}
      
      }
    }
   
     
  
      j++;
  }

  
}

void printGrid(double ** grid){
  
  for(int i=0;i<N;i++){
    for(int j=0;j<N;j++){
      cout<<grid[i][j]<<" ";
    }
    cout<<"\n";
  }

}



