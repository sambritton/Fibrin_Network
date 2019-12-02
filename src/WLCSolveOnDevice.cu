/*
* WLCSolveOnDevice.cu
*
* Created on 8/7/2017 
* 		Author: SRB
*/
#include "NodeSystemDeviceFunctors.h"
#include "NodeSystemDevice.h"
#include "WLCSolveOnDevice.h" 
#include <functional>
#include <algorithm>    // std::transform
#include <vector>     
#include <math.h>  

inline double CVec3_dot(CVec3 v1, CVec3 v2) {
	return (thrust::get<0>(v1)*thrust::get<0>(v2) +
		thrust::get<1>(v1)*thrust::get<1>(v2) +
		thrust::get<2>(v1)*thrust::get<2>(v2));
};
void WLCSolveOnDevice(
	NodeInfoVecs& nodeInfoVecs,
	WLCInfoVecs& wlcInfoVecs,  
	GeneralParams& generalParams) {
 
 
	thrust::counting_iterator<unsigned> startEdgeIter(0);
			  
	//
	thrust::for_each( 
		thrust::make_zip_iterator( 
			thrust::make_tuple(startEdgeIter,
								nodeInfoVecs.springDivisionCount.begin(),
								nodeInfoVecs.isNodeFixed.begin() )),
		thrust::make_zip_iterator(
			thrust::make_tuple(startEdgeIter,
								nodeInfoVecs.springDivisionCount.begin(),
								nodeInfoVecs.isNodeFixed.begin() )) + generalParams.maxNodeCount,
		WLCfunctor(
			thrust::raw_pointer_cast(nodeInfoVecs.nodeLocX.data()),
			thrust::raw_pointer_cast(nodeInfoVecs.nodeLocY.data()),
			thrust::raw_pointer_cast(nodeInfoVecs.nodeLocZ.data()),
			thrust::raw_pointer_cast(nodeInfoVecs.nodeForceX.data()),
			thrust::raw_pointer_cast(nodeInfoVecs.nodeForceY.data()),
			thrust::raw_pointer_cast(nodeInfoVecs.nodeForceZ.data()),
 
			generalParams.kB,
			generalParams.persistenceLengthMon,
			generalParams.CLM,
			generalParams.temperature,
			generalParams.maxNeighborCount,
			generalParams.maxNodeCount,

			thrust::raw_pointer_cast(wlcInfoVecs.lengthZero.data()),
			thrust::raw_pointer_cast(nodeInfoVecs.springDivisionCount.data()),
			thrust::raw_pointer_cast(wlcInfoVecs.globalNeighbors.data()),
			thrust::raw_pointer_cast(wlcInfoVecs.currentNodeEdgeCountVector.data()),
			thrust::raw_pointer_cast(wlcInfoVecs.numOriginalNeighborsNodeVector.data()) ) );
};

void GetStrainParameters(
	NodeInfoVecs& nodeInfoVecs,
	WLCInfoVecs& wlcInfoVecs,  
	GeneralParams& generalParams,
	DomainParams& domainParams) {
		


		//count positive and negative strains for edges that are not added. If an edge is added, a zero is placed on that strain.
		//notice that each thread will count edges twice, so divide by two at the end
		thrust::counting_iterator<unsigned> startStrainIter(0);

		thrust::fill(nodeInfoVecs.discretizedEdgeStrain.begin(), nodeInfoVecs.discretizedEdgeStrain.end(),0.0);
		thrust::fill(nodeInfoVecs.discretizedEdgeAlignment.begin(), nodeInfoVecs.discretizedEdgeAlignment.end(),0.0);	

		thrust::transform(
			thrust::make_zip_iterator(
				thrust::make_tuple(
					nodeInfoVecs.deviceEdgeLeft.begin(),
					nodeInfoVecs.deviceEdgeRight.begin())),
					 
			thrust::make_zip_iterator(
				thrust::make_tuple(
					nodeInfoVecs.deviceEdgeLeft.begin(),
					nodeInfoVecs.deviceEdgeRight.begin())) + generalParams.currentEdgeCount,
					
			//outputs discretized strain etc			
			thrust::make_zip_iterator(
				thrust::make_tuple(
					nodeInfoVecs.discretizedEdgeStrain.begin(),
					nodeInfoVecs.discretizedEdgeAlignment.begin())),
					
			CalculateStrainParamsFunctor(	
				generalParams.originLinkCount,
				generalParams.originEdgeCount,
				generalParams.originNodeCount,
				generalParams.maxNodeCount,
				generalParams.maxNeighborCount,
				thrust::raw_pointer_cast(nodeInfoVecs.nodeLocX.data()),
				thrust::raw_pointer_cast(nodeInfoVecs.nodeLocY.data()),
				thrust::raw_pointer_cast(nodeInfoVecs.nodeLocZ.data()),
				thrust::raw_pointer_cast(wlcInfoVecs.numOriginalNeighborsNodeVector.data()),
				thrust::raw_pointer_cast(wlcInfoVecs.currentNodeEdgeCountVector.data()),
				thrust::raw_pointer_cast(wlcInfoVecs.globalNeighbors.data()),
				thrust::raw_pointer_cast(wlcInfoVecs.lengthZero.data()) ));
			
}; 

  