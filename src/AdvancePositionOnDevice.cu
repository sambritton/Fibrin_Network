/*
* AdvancePositionOnDevice.cu
*
* Created on 11/5/2017 
* 		Author: SRB
*/
#include "NodeSystemDevice.h" 
#include "AdvancePositionOnDevice.h"
 
double AdvancePositionOnDevice(
	NodeInfoVecs& nodeInfoVecs,
	GeneralParams& generalParams,
	DPDParticleVariables& dpdParticleVariables) {

   
		//At this point, the previous node location is the same as the current node, 
		//we can therefore use previous node locations to update nodeLoc. 
		
		//set high for not use 
		 unsigned _seed = rand();
    	thrust::device_vector<double> gaussianData;
    	gaussianData.resize(generalParams.maxNodeCount); // 
		thrust::counting_iterator<unsigned> index_sequence_begin(_seed);
		
    	thrust::transform(thrust::device, index_sequence_begin, index_sequence_begin + (generalParams.maxNodeCount), 
        gaussianData.begin(), psrunifgen(-1.0, 1.0));//unif [-1,1];
		  
		
		thrust::counting_iterator<unsigned> nodeIndexBegin(0);
		  
		thrust::transform( 
			thrust::make_zip_iterator(
				thrust::make_tuple(
					nodeIndexBegin,
					nodeInfoVecs.nodeLocX.begin(),
					nodeInfoVecs.nodeLocY.begin(),
					nodeInfoVecs.nodeLocZ.begin())),
			thrust::make_zip_iterator(
				thrust::make_tuple(
					nodeIndexBegin,
					nodeInfoVecs.nodeLocX.begin(),
					nodeInfoVecs.nodeLocY.begin(),
					nodeInfoVecs.nodeLocZ.begin())) + generalParams.maxNodeCount,
			//second vector begin
			thrust::make_zip_iterator(
				thrust::make_tuple(
					gaussianData.begin(),
				/*	nodeInfoVecs.prevNodeForceX.begin(),
					nodeInfoVecs.prevNodeForceY.begin(),
					nodeInfoVecs.prevNodeForceZ.begin(),*/
					nodeInfoVecs.nodeForceX.begin(),
					nodeInfoVecs.nodeForceY.begin(),
					nodeInfoVecs.nodeForceZ.begin())),
			//save result in third vector to test values
			thrust::make_zip_iterator(
				thrust::make_tuple(
					nodeInfoVecs.nodeLocX.begin(),
					nodeInfoVecs.nodeLocY.begin(),
					nodeInfoVecs.nodeLocZ.begin(),
					nodeInfoVecs.nodeVelocity.begin())),
			SaxpyFunctorDim3(generalParams.dtTemp,
				generalParams.viscousDamp, 
				generalParams.temperature,
				generalParams.kB,
				generalParams.nodeMass,
				generalParams.maxNodeCount,
				thrust::raw_pointer_cast(nodeInfoVecs.isNodeFixed.data())));
		cudaThreadSynchronize();
			/*	thrust::copy(nodeInfoVecs.nodeForceX.begin() + generalParams.maxNodeCount, 
					nodeInfoVecs.nodeForceX.end(), 
					nodeInfoVecs.prevNodeForceX.begin() + generalParams.maxNodeCount);
				thrust::copy(nodeInfoVecs.nodeForceY.begin() + generalParams.maxNodeCount, 
					nodeInfoVecs.nodeForceY.end(), 
					nodeInfoVecs.prevNodeForceY.begin() + generalParams.maxNodeCount);
				thrust::copy(nodeInfoVecs.nodeForceZ.begin() + generalParams.maxNodeCount, 
					nodeInfoVecs.nodeForceZ.end(), 
					nodeInfoVecs.prevNodeForceZ.begin() + generalParams.maxNodeCount);*/
			
		//finally, clear the random data.
        gaussianData.clear();
        gaussianData.shrink_to_fit();

	return generalParams.dtTemp;
		//now that nodeLoc is different, we can calculate change and then set previous location
		//to the current location. 
	
}
