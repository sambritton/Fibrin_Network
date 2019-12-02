//look up forward declaration and inheritance if this doesn't make sense. 
#include "PlateletForceDevice.h" 
#include "NodeSystemDevice.h"


typedef thrust::tuple<unsigned, bool> Tub; //tuple holding id of a node and if that node is within reach of the platelet

__host__ __device__

void PlateletForceOnDevice(
    NodeInfoVecs& nodeInfoVecs,
	WLCInfoVecs& wlcInfoVecs,
	GeneralParams& generalParams,
    PlateletParams& plateletParams,
    PlateletInfoVecs& plateletInfoVecs) {

        thrust::counting_iterator<unsigned> indexBegin(0);
        thrust::counting_iterator<unsigned> indexEnd(plateletParams.maxPlateletCount);
        
        //Each of these 
        thrust::device_vector<double> tempNodeForceX;
        thrust::device_vector<double> tempNodeForceY;
        thrust::device_vector<double> tempNodeForceZ;

        //Call the platelet force functor
        thrust::for_each(indexBegin, indexEnd,
            PlateletForceFunctor(
                thrust::raw_pointer_cast(plateletInfoVecs.plateletLocX.data()),
                thrust::raw_pointer_cast(plateletInfoVecs.plateletLocY.data()),
                thrust::raw_pointer_cast(plateletInfoVecs.plateletLocZ.data()),
                thrust::raw_pointer_cast(plateletInfoVecs.plateletForceX.data()),
                thrust::raw_pointer_cast(plateletInfoVecs.plateletForceY.data()),
                thrust::raw_pointer_cast(plateletInfoVecs.plateletForceZ.data()),
                plateletParams.numberOfArms,
                plateletParams.pullingForce,
                thrust::raw_pointer_cast(nodeInfoVecs.nodeLocX.data()),
                thrust::raw_pointer_cast(nodeInfoVecs.nodeLocY.data()),
                thrust::raw_pointer_cast(nodeInfoVecs.nodeLocZ.data()),
                thrust::raw_pointer_cast(nodeInfoVecs.nodeForceX.data()),
                thrust::raw_pointer_cast(nodeInfoVecs.nodeForceY.data()),
                thrust::raw_pointer_cast(nodeInfoVecs.nodeForceZ.data())
                ))
    
        //now call a sort by key followed by a reduce by key to figure out which nodes are have force applied.
        //then make a functor that takes the id and force (4 tuple) and takes that force and adds it to the id'th entry in nodeInfoVecs.nodeForceX,Y,Z 

    };

    void AdvancePlateletPosition(//stuff here) {
        //stuff here
    };