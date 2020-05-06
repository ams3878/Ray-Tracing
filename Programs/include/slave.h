#ifndef __SLAVE_PROCESS_H__
#define __SLAVE_PROCESS_H__

#include "RayTrace.h"

void slaveMain( ConfigData *data );
//This function will perform ray tracing when no MPI use was
//given.
//
//Inputs:
//    data - the ConfigData that holds the scene information.
//
//Outputs: None
void slaveStripsHorizontal(ConfigData *data);
void slaveStripsVertical(ConfigData *data);
void slaveStaticBlocks(ConfigData *data);
void slaveDynamicBlocks(ConfigData *data);
void slaveCyclesHorizontal(ConfigData *data);
void slaveCyclesVertical(ConfigData *data);
#endif
