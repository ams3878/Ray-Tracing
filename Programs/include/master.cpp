//This file contains the code that the master process will execute.

#include <iostream>
#include <mpi.h>
#include <math.h>

#include "RayTrace.h"
#include "master.h"

void masterMain(ConfigData* data)
{
    //Depending on the partitioning scheme, different things will happen.
    //You should have a different function for each of the required
    //schemes that returns some values that you need to handle.

    //Allocate space for the image on the master.
    float* pixels = new float[3 * data->width * data->height];

    //Execution time will be defined as how long it takes
    //for the given function to execute based on partitioning
    //type.
    double renderTime = 0.0, startTime, stopTime;
	//Add the required partitioning methods here in the case statement.
	//You do not need to handle all cases; the default will catch any
	//statements that are not specified. This switch/case statement is the
	//only place that you should be adding code in this function. Make sure
	//that you update the header files with the new functions that will be
	//called.
	//It is suggested that you use the same parameters to your functions as shown
	//in the sequential example below.
    switch (data->partitioningMode)
    {
        case PART_MODE_NONE:
          //Call the function that will handle this.
          startTime = MPI_Wtime();
          masterSequential(data, pixels);
          stopTime = MPI_Wtime();
          break;
        case PART_MODE_STATIC_STRIPS_HORIZONTAL:
          startTime = MPI_Wtime();
          masterStripsHorizontal(data, pixels);
          stopTime = MPI_Wtime();
          break;
//        case PART_MODE_STATIC_STRIPS_VERTICAL:
//            startTime = MPI_Wtime();
//            masterCyclesHorizontal(data, pixels);
//            stopTime = MPI_Wtime();
//            break;
        case PART_MODE_STATIC_BLOCKS:
          startTime = MPI_Wtime();
          masterStaticBlocks(data, pixels);
          stopTime = MPI_Wtime();
          break;
//        case PART_MODE_STATIC_CYCLES_HORIZONTAL:
//            startTime = MPI_Wtime();
//            masterCyclesHorizontal(data, pixels);
//            stopTime = MPI_Wtime();
//            break;
        case PART_MODE_STATIC_CYCLES_VERTICAL:
          startTime = MPI_Wtime();
          masterCyclesVertical(data, pixels);
          stopTime = MPI_Wtime();
          break;
        default:
          std::cout << "This mode (" << data->partitioningMode;
          std::cout << ") is not currently implemented." << std::endl;
          break;
    }

    renderTime = stopTime - startTime;
    std::cout << "Execution Time: " << renderTime << " seconds" << std::endl << std::endl;

    //After this gets done, save the image.
    std::cout << "Image will be save to: ";
    std::string file = generateFileName(data);
    std::cout << file << std::endl;
    savePixels(file, pixels, data);

    //Delete the pixel data.
    delete[] pixels;
}

void masterStripsHorizontal(ConfigData* data, float* pixels)
{
    MPI_Status status;
    int barrier;
    //Start the computation time timer.
    double computationStart = MPI_Wtime();
    //Render the scene.
    int iters = data->height/data->mpi_procs;
    if (data->mpi_rank < (data->height % data->mpi_procs)){
      iters++;
    }
    for( int i = 0; i < iters; i++ ){
        for( int j = 0; j < data->width; ++j ){
            int row = i;
            int column = j;
            //Calculate the index into the array.
            int baseIndex = 3 * ( row * data->width + column );
            //Call the function to shade the pixel.
            shadePixel(&(pixels[baseIndex]),row,column,data);
          }
    }
    //Stop the comp. timer
    double computationStop = MPI_Wtime();
    double communicationStart = MPI_Wtime();
    double computationTime = computationStop - computationStart;
    double tempTime = 0.0;
    float* tempPixel = new float[3* data->width];
    for(int i = 0; i < data->height - iters ;i++ ){
      MPI_Recv(tempPixel,3*data->width,MPI_FLOAT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
      memcpy(&(pixels[3*(status.MPI_TAG-1)*data->width]), tempPixel,3*data->width * sizeof(float));
    }
    for(int i = 1; i < data->mpi_procs;i++ ){
      MPI_Send(&barrier,1,MPI_INT,i,0,MPI_COMM_WORLD);
      MPI_Recv(&tempTime,1,MPI_DOUBLE,i,0,MPI_COMM_WORLD,&status);
      computationTime += tempTime;
    }
    //After receiving from all processes, the communication time will
    //be obtained.
    double communicationStop = MPI_Wtime();
    double communicationTime = communicationStop - communicationStart;
    //Print the times and the c-to-c ratio
	//This section of printing, IN THIS ORDER, needs to be included in all of the
	//functions that you write at the end of the function.
    std::cout << "Total Computation Time: " << computationTime << " seconds" << std::endl;
    std::cout << "Total Communication Time: " << communicationTime << " seconds" << std::endl;
    double c2cRatio = communicationTime / computationTime;
    std::cout << "C-to-C Ratio: " << c2cRatio << std::endl;
}
//NOT DONE
void masterStripsVertical(ConfigData* data, float* pixels)
{
    MPI_Status status;
    int barrier;
    //Start the computation time timer.
    double communicationStart = MPI_Wtime();
    double computationStart = MPI_Wtime();
    //Render the scene.
    int iters = data->width/data->mpi_procs;
    if (data->mpi_rank < (data->width % data->mpi_procs)){
      iters++;
    }
    for( int i = 0; i < iters; i++ ){
        for( int j = 0; j < data->height; ++j ){
            int row = j;
            int column = i;
            //Calculate the index into the array.
            int baseIndex = 3 * ( row * data->width + column );
            //Call the function to shade the pixel.
            shadePixel(&(pixels[baseIndex]),row,column,data);
          }
    }
    //Stop the comp. timer
    double computationStop = MPI_Wtime();
    double computationTime = computationStop - computationStart;
    double tempTime = 0.0;
    float* tempPixel = new float[3* data->width];
    for(int i = 0; i <= data->height - iters-1 ;i++ ){
      MPI_Recv(tempPixel,3*data->width,MPI_FLOAT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
      memcpy(&(pixels[3*((status.MPI_TAG-1) * data->mpi_procs + status.MPI_SOURCE)*data->width]), tempPixel,3*data->width * sizeof(float));
    }
    for(int i = 1; i < data->mpi_procs;i++ ){
      MPI_Send(&barrier,1,MPI_INT,i,0,MPI_COMM_WORLD);
      MPI_Recv(&tempTime,1,MPI_DOUBLE,i,0,MPI_COMM_WORLD,&status);
      computationTime += tempTime;
    }
    //After receiving from all processes, the communication time will
    //be obtained.
    double communicationStop = MPI_Wtime();
    double communicationTime = communicationStop - communicationStart - computationTime;
    //Print the times and the c-to-c ratio
	//This section of printing, IN THIS ORDER, needs to be included in all of the
	//functions that you write at the end of the function.
    std::cout << "Total Computation Time: " << computationTime << " seconds" << std::endl;
    std::cout << "Total Communication Time: " << communicationTime << " seconds" << std::endl;
    double c2cRatio = communicationTime / computationTime;
    std::cout << "C-to-C Ratio: " << c2cRatio << std::endl;
}
void masterStaticBlocks(ConfigData* data, float* pixels)
{
  MPI_Status status;
  int barrier;
  //Start the computation time timer.
  double computationStart = MPI_Wtime();
  //Render the scene.
  int iters = (data->width * data->width) /(data->mpi_procs);
  iters = (int) sqrt(iters);
  if (data->mpi_rank < (data->height % data->mpi_procs)){
    iters++;
  }
  for( int i = 0; i < iters; i++ ){
      for( int j = 0; j < data->width; ++j ){
          int row = i;
          int column = j;
          //Calculate the index into the array.
          int baseIndex = 3 * ( row * data->width + column );
          //Call the function to shade the pixel.
          shadePixel(&(pixels[baseIndex]),row,column,data);
        }
  }
  //Stop the comp. timer
  double computationStop = MPI_Wtime();
  double communicationStart = MPI_Wtime();
  double computationTime = computationStop - computationStart;
  double tempTime = 0.0;
  float* tempPixel = new float[3* data->width];
  for(int i = 0; i < data->mpi_procs;i++ ){
    MPI_Recv(tempPixel,3*iters*iters,MPI_FLOAT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
    memcpy(&(pixels[3*(((status.MPI_SOURCE - 1)/iters) % (data->height/iters)) * iters
                    + ((status.MPI_SOURCE - 1) % (data->height/iters)) * iters)]), tempPixel[0],3 *iters* sizeof(float));
    memcpy(&(pixels[3*(((status.MPI_SOURCE - 1)/iters) % (data->height/iters)) * iters
                    + ((status.MPI_SOURCE - 1) % (data->height/iters)) * iters)]), tempPixel[3 *iters],3 *iters* sizeof(float));
    memcpy(&(pixels[3*(((status.MPI_SOURCE - 1)/iters) % (data->height/iters)) * iters
                    + ((status.MPI_SOURCE - 1) % (data->height/iters)) * iters)]), tempPixel[6 *iters],3 *iters* sizeof(float));
  }
  for(int i = 1; i < data->mpi_procs;i++ ){
    MPI_Send(&barrier,1,MPI_INT,i,0,MPI_COMM_WORLD);
    MPI_Recv(&tempTime,1,MPI_DOUBLE,i,0,MPI_COMM_WORLD,&status);
    computationTime += tempTime;
  }
  //After receiving from all processes, the communication time will
  //be obtained.
  double communicationStop = MPI_Wtime();
  double communicationTime = communicationStop - communicationStart;
  //Print the times and the c-to-c ratio
//This section of printing, IN THIS ORDER, needs to be included in all of the
//functions that you write at the end of the function.
  std::cout << "Total Computation Time: " << computationTime << " seconds" << std::endl;
  std::cout << "Total Communication Time: " << communicationTime << " seconds" << std::endl;
  double c2cRatio = communicationTime / computationTime;
  std::cout << "C-to-C Ratio: " << c2cRatio << std::endl;
}
void masterDynamicBlocks(ConfigData *data, float* pixels){

}

//NOT TESTED
void masterCyclesHorizontal(ConfigData* data, float* pixels)
{
    MPI_Status status;
    int barrier;
    //Start the computation time timer.
    double communicationStart = MPI_Wtime();
    double computationStart = MPI_Wtime();
    //Render the scene.
    int iters = data->height/data->mpi_procs;
    if (data->mpi_rank < (data->height % data->mpi_procs)){
      iters++;
    }
    for( int i = 0; i < iters; i++ ){
        for( int j = 0; j < data->width; ++j ){
            int row = (i*data->mpi_procs);
            int column = j;
            //Calculate the index into the array.
            int baseIndex = 3 * ( row * data->width + column );
            //Call the function to shade the pixel.
            shadePixel(&(pixels[baseIndex]),row,column,data);
          }
    }
    //Stop the comp. timer
    double computationStop = MPI_Wtime();
    double computationTime = computationStop - computationStart;
    double tempTime = 0.0;
    float* tempPixel = new float[3* data->width];
    for(int i = 0; i <= data->height - iters ;i++ ){
      MPI_Recv(tempPixel,3*data->width,MPI_FLOAT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
      memcpy(&(pixels[3*((status.MPI_TAG-1) * data->mpi_procs + status.MPI_SOURCE)*data->width]), tempPixel,3*data->width * sizeof(float));
    }
    for(int i = 1; i < data->mpi_procs;i++ ){
      MPI_Send(&barrier,1,MPI_INT,i,0,MPI_COMM_WORLD);
      MPI_Recv(&tempTime,1,MPI_DOUBLE,i,0,MPI_COMM_WORLD,&status);
      computationTime += tempTime;
    }
    //After receiving from all processes, the communication time will
    //be obtained.
    double communicationStop = MPI_Wtime();
    double communicationTime = communicationStop - communicationStart - computationTime;
    //Print the times and the c-to-c ratio
	//This section of printing, IN THIS ORDER, needs to be included in all of the
	//functions that you write at the end of the function.
    std::cout << "Total Computation Time: " << computationTime << " seconds" << std::endl;
    std::cout << "Total Communication Time: " << communicationTime << " seconds" << std::endl;
    double c2cRatio = communicationTime / computationTime;
    std::cout << "C-to-C Ratio: " << c2cRatio << std::endl;
}
void masterCyclesVertical(ConfigData* data, float* pixels)
{
    float* pixelsT = new float[3 * data->width * data->height];
    MPI_Status status;
    int barrier = 0;
    int row;
    int column;
    //Start the computation time timer.
    double computationStart = MPI_Wtime();
    //Render the scene.
    int iters = (data->width/(data->mpi_procs * data->cycleSize));
    if( data->mpi_rank * data->cycleSize < (data->width % (data->mpi_procs * data->cycleSize))){
      iters++;
    }
    for( int i = 0; i < iters; i++ ){
      for(int k = 0; k < data->cycleSize; ++k){
        column = i*(data->mpi_procs * data->cycleSize) + k ;
        barrier++;
        if ( column >= data->width) { break;}
        for( int j = 0; j < data->height; ++j )
        {
            row = j;
            //Calculate the index into the array.
            int baseIndex = 3 * ( row * data->width + column );
            //Call the function to shade the pixel.
            shadePixel(&(pixels[baseIndex]),row,column,data);
        }
      }
    }
    for( int i=0; i < data->width; ++i){
       for(int j=0; j < data->height; ++j)
       {
           pixelsT[3*(j * data->height + i)] = pixels[3*(i *data->width + j)];
           pixelsT[3*(j * data->height + i) + 1] = pixels[3*(i *data->width + j) + 1];
           pixelsT[3*(j * data->height + i) + 2] = pixels[3*(i *data->width + j) + 2];
       }
     }
    //Stop the comp. timer
    double computationStop = MPI_Wtime();
    double computationTime = computationStop - computationStart;
    double tempTime = 0.0;
    double communicationStart = MPI_Wtime();
    float* tempPixel = new float[3* data->width];
    for(int i = 0; i < data->width - barrier ;i++ ){
      MPI_Recv(tempPixel,3*data->height,MPI_FLOAT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
      memcpy(&(pixelsT[3*(status.MPI_TAG)*data->height]), tempPixel,3*data->height * sizeof(float));
    }
    for(int i=0; i < data->width; ++i){
       for(int j=0; j < data->height; ++j)
       {
         pixels[3*(j * data->width + i)] = pixelsT[3*(i *data->height + j)];
         pixels[3*(j * data->width + i) + 1] = pixelsT[3*(i *data->height + j) + 1];
         pixels[3*(j * data->width + i) + 2] = pixelsT[3*(i *data->height + j) + 2];
       }
    }
    for(int i = 1; i < data->mpi_procs;i++ ){
      MPI_Send(&barrier,1,MPI_INT,i,0,MPI_COMM_WORLD);
      MPI_Recv(&tempTime,1,MPI_DOUBLE,i,0,MPI_COMM_WORLD,&status);
      computationTime += tempTime;
    }

    //After receiving from all processes, the communication time will
    //be obtained.
    double communicationStop = MPI_Wtime();
    double communicationTime = communicationStop - communicationStart;
    //Print the times and the c-to-c ratio
	//This section of printing, IN THIS ORDER, needs to be included in all of the
	//functions that you write at the end of the function.
    std::cout << "Total Computation Time: " << computationTime << " seconds" << std::endl;
    std::cout << "Total Communication Time: " << communicationTime << " seconds" << std::endl;
    double c2cRatio = communicationTime / computationTime;
    std::cout << "C-to-C Ratio: " << c2cRatio << std::endl;
}
void masterSequential(ConfigData* data, float* pixels)
{
    //Start the computation time timer.
    double computationStart = MPI_Wtime();

    //Render the scene.
    for( int i = 0; i < data->height; ++i )
    {
        for( int j = 0; j < data->width; ++j )
        {
            int row = i;
            int column = j;

            //Calculate the index into the array.
            int baseIndex = 3 * ( row * data->width + column );

            //Call the function to shade the pixel.
            shadePixel(&(pixels[baseIndex]),row,column,data);
        }
    }

    //Stop the comp. timer
    double computationStop = MPI_Wtime();
    double computationTime = computationStop - computationStart;

    //After receiving from all processes, the communication time will
    //be obtained.
    double communicationTime = 0.0;

    //Print the times and the c-to-c ratio
	//This section of printing, IN THIS ORDER, needs to be included in all of the
	//functions that you write at the end of the function.
    std::cout << "Total Computation Time: " << computationTime << " seconds" << std::endl;
    std::cout << "Total Communication Time: " << communicationTime << " seconds" << std::endl;
    double c2cRatio = communicationTime / computationTime;
    std::cout << "C-to-C Ratio: " << c2cRatio << std::endl;
}
