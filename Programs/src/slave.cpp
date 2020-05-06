//This file contains the code that the master process will execute.

#include <iostream>
#include <mpi.h>
#include <math.h>
#include <stdlib.h>

#include "RayTrace.h"
#include "slave.h"

void slaveMain(ConfigData* data){
     //Depending on the partitioning scheme, different things will happen.
    //You should have a different function for each of the required
    //schemes that returns some values that you need to handle.
    switch (data->partitioningMode)
    {
        case PART_MODE_NONE:
            //The slave will do nothing since this means sequential operation.
            break;
        case PART_MODE_STATIC_STRIPS_HORIZONTAL:
            slaveStripsHorizontal(data);
            break;
        case PART_MODE_STATIC_BLOCKS:
            slaveStaticBlocks(data);
            break;
        case PART_MODE_DYNAMIC:
            slaveDynamicBlocks(data);
            break;
        case PART_MODE_STATIC_CYCLES_VERTICAL:
            slaveCyclesVertical(data);
            break;
        default:
            std::cout << "This mode (" << data->partitioningMode;
            std::cout << ") is not currently implemented." << std::endl;
            break;
    }
}
void slaveStripsHorizontal(ConfigData* data){
  float* myPixels = new float[3* data->width];
    //Start the computation time timer.
    double computationStart,computationStop,computationTime;
    MPI_Status status;
    int barrier,row,column;
    //Render the scene.
    int iters = data->height/data->mpi_procs;
    int myiters = iters;
    int offset = data->height % data->mpi_procs;
    if (data->mpi_rank < offset){
      myiters++;
      offset = data->mpi_rank;
    }
    for( int i = 0; i < myiters; i++ )
    {
        for( int j = 0; j < data->width; ++j )
        {
            row = i + (offset) + (iters * data->mpi_rank);
            column = j;
            int baseIndex = 3 * ( column );
            computationStart = MPI_Wtime();
            shadePixel(&(myPixels[baseIndex]),row,column,data);
            computationStop = MPI_Wtime();
            computationTime += (computationStop - computationStart);
        }
        MPI_Send(myPixels,3*data->width,MPI_FLOAT,0,row + 1,MPI_COMM_WORLD);

    }
    //Stop the comp. timer
    MPI_Recv(&barrier,1,MPI_INT,0,0,MPI_COMM_WORLD,&status);
    MPI_Send(&computationTime,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
    delete[] myPixels;

}
void slaveStaticBlocks(ConfigData* data){
  //Start the computation time timer.
  double computationStart,computationStop,computationTime;
  MPI_Status status;
  int barrier,row,column;
  //Render the scene.
  int iters = (data->height * data->width) /(data->mpi_procs);
  iters = (int) sqrt(iters);
  float* myPixels = new float[3* iters * iters];

  computationStart = MPI_Wtime();
  for( int i = 0; i < iters; i++ ){
      for( int j = 0; j < iters; j++ )
      {
        row = ((((data->mpi_rank - 1)/(data->width/iters) ) % (data->width/iters)) * iters) + i ;
        column = ((data->mpi_rank - 1) % (data->width/iters)) * iters + j;
        int baseIndex = 3 * (((i) * iters) + j );
        shadePixel(&(myPixels[baseIndex]),row,column,data);
      }
  }
  computationStop = MPI_Wtime();
  computationTime += (computationStop - computationStart);
  MPI_Send(myPixels,3*iters*iters,MPI_FLOAT,0,1,MPI_COMM_WORLD);
  //Stop the comp. timer
  MPI_Recv(&barrier,1,MPI_INT,0,0,MPI_COMM_WORLD,&status);
  MPI_Send(&computationTime,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
  delete[] myPixels;
}
void slaveDynamicBlocks(ConfigData *data){
  //Start the computation time timer.
  double computationStart,computationStop,computationTime;
  MPI_Status status;
  int row,column;
  //Render the scene.
  int iters = 0;
  int block_size = 3* data->dynamicBlockWidth*data->dynamicBlockHeight;
  float* myPixels = NULL;
  int*   myblocks = NULL;
  int* block = new int[2];
  int ready;

  while( 1 ){
     MPI_Send(&ready,1,MPI_INT,0,0,MPI_COMM_WORLD);
     MPI_Recv(block,2,MPI_INT,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
     if( status.MPI_TAG == 0){ break;}
     computationStart = MPI_Wtime();
     myPixels = (float*) realloc(myPixels,++iters * block_size * sizeof(float));
     myblocks = (int*) realloc(myblocks,iters * 2 * sizeof(int));
     myblocks[(iters-1)*2] = block[0];
     myblocks[(iters-1)*2 + 1] = block[1];
     for( int i = 0; i < data->dynamicBlockHeight; i++ ){
         for( int j = 0; j < data->dynamicBlockWidth; j++ )
         {
           row = block[0] + i;
           column = block[1] + j;
           int baseIndex = (3 * ((i * data->dynamicBlockWidth) + (j) ))
                              + ((iters-1) * block_size);
           if(row >= data->height || column >= data->width){
             myPixels[baseIndex] = 0;
           }else {
             shadePixel(&(myPixels[baseIndex]),row,column,data);
           }
         }
     }
     computationStop = MPI_Wtime();
     computationTime += (computationStop - computationStart);
  }
  MPI_Send(&iters,1,MPI_INT,0,1,MPI_COMM_WORLD);
  MPI_Send(myblocks,iters*2,MPI_INT,0,1,MPI_COMM_WORLD);
  MPI_Send(myPixels,iters*block_size,MPI_FLOAT,0,1,MPI_COMM_WORLD);
  MPI_Send(&computationTime,1,MPI_DOUBLE,0,1,MPI_COMM_WORLD);

  delete[] myPixels;
}
void slaveCyclesVertical(ConfigData* data){
    float* myPixels = new float[3* data->height];

    //Start the computation time timer.
    double computationStart,computationTime;
    MPI_Status status;
    int barrier,row,column;
    //Render the scene.
    int iters = (data->width/(data->mpi_procs * data->cycleSize));
    if( data->mpi_rank * data->cycleSize < (data->width % (data->mpi_procs * data->cycleSize))){
      iters++;
    }

    for( int i = 0; i < iters; i++ ){
      for(int k = 0; k < data->cycleSize; ++k){
        column = (((i*data->mpi_procs)+data->mpi_rank)* data->cycleSize) + k;
        if ( column >= data->width) { break;}
        for( int j = 0; j < data->height; ++j )
        {
            row = j;
            int baseIndex = 3 * ( row );
            computationStart = MPI_Wtime();
            shadePixel(&(myPixels[baseIndex]),row,column,data);
            computationTime += (MPI_Wtime() - computationStart);
        }
        MPI_Send(myPixels,3*data->height,MPI_FLOAT,0,column,MPI_COMM_WORLD);
      }
    }
    //Stop the comp. timer
    MPI_Recv(&barrier,1,MPI_INT,0,0,MPI_COMM_WORLD,&status);
    MPI_Send(&computationTime,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
    delete[] myPixels;
}
