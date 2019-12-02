/*
* TorsionSolveOnDevice.cu
*
* Created on 11/7/2017
* 		Author: SRB
*/ 

#include "NodeSystemDevice.h"
#include "TorsionSolveOnDevice.h"

void TorsionSolveOnDevice(
	NodeInfoVecs& nodeInfoVecs,
	TorsionInfoVecs& torsionInfoVecs,
	GeneralParams& generalParams)  {
	
const double PI = 3.14159265358979323846;  
if (generalParams.totalTorsionCount>0) { 

		thrust::counting_iterator<unsigned> startTorsionIter(0);
		thrust::counting_iterator<unsigned> endTorsionIter(generalParams.totalTorsionCount);
 
		//for_each guarrantees order. This is needed for iter count and saving to torsion force vectors.
		//forces are filled using 3 counters left = counter, center = counter + totalTorsionCount etc. 
		//Thus, in the force vector, only the first 3*totalTorsionCount entries are filled. 
		thrust::for_each(
			thrust::make_zip_iterator(
				thrust::make_tuple(
					startTorsionIter,
					torsionInfoVecs.leftIndex.begin(),
					torsionInfoVecs.centerIndex.begin(),
					torsionInfoVecs.rightIndex.begin(),
					torsionInfoVecs.angleZero.begin())),
			thrust::make_zip_iterator(
				thrust::make_tuple(
					startTorsionIter,
					torsionInfoVecs.leftIndex.begin(),
					torsionInfoVecs.centerIndex.begin(),
					torsionInfoVecs.rightIndex.begin(),
					torsionInfoVecs.angleZero.begin())) + generalParams.totalTorsionCount,
			TorsionFunctor(
				thrust::raw_pointer_cast(nodeInfoVecs.nodeLocX.data()),
				thrust::raw_pointer_cast(nodeInfoVecs.nodeLocY.data()),
				thrust::raw_pointer_cast(nodeInfoVecs.nodeLocZ.data()),
				thrust::raw_pointer_cast(torsionInfoVecs.forceX.data()),
				thrust::raw_pointer_cast(torsionInfoVecs.forceY.data()),
				thrust::raw_pointer_cast(torsionInfoVecs.forceZ.data()),
				thrust::raw_pointer_cast(nodeInfoVecs.isNodeFixed.data()),
				generalParams.torsionStiffness,
				generalParams.maxNodeCount,
				generalParams.totalTorsionCount,
				PI));  

		cudaThreadSynchronize();
		//reduce by key to get forces.Notice leftIndex is 1/3rd the length of torsion.forceX
		//this vector will be sorted each iteration, so it needs to be recopied unfortunately.
		//fill must end before non-set id's
		thrust::copy(torsionInfoVecs.leftIndex.begin(), torsionInfoVecs.leftIndex.begin() + generalParams.totalTorsionCount,
			torsionInfoVecs.tempTorIndices.begin());

		thrust::copy(torsionInfoVecs.centerIndex.begin(), torsionInfoVecs.centerIndex.begin() + generalParams.totalTorsionCount, 
			torsionInfoVecs.tempTorIndices.begin() + generalParams.totalTorsionCount);

		thrust::copy(torsionInfoVecs.rightIndex.begin(), torsionInfoVecs.rightIndex.begin() + generalParams.totalTorsionCount, 
			torsionInfoVecs.tempTorIndices.begin() + 2 * generalParams.totalTorsionCount);


		//key, then value. Each vector returns sorted		
		thrust::sort_by_key(torsionInfoVecs.tempTorIndices.begin(), torsionInfoVecs.tempTorIndices.begin() + 3 * generalParams.totalTorsionCount,
			thrust::make_zip_iterator(
				thrust::make_tuple(
					torsionInfoVecs.forceX.begin(),
					torsionInfoVecs.forceY.begin(),
					torsionInfoVecs.forceZ.begin())), thrust::less<unsigned>());


		thrust::fill(torsionInfoVecs.tempForceX.begin(), torsionInfoVecs.tempForceX.end(), 0);
		thrust::fill(torsionInfoVecs.tempForceY.begin(), torsionInfoVecs.tempForceY.end(), 0);
		thrust::fill(torsionInfoVecs.tempForceZ.begin(), torsionInfoVecs.tempForceZ.end(), 0);
		thrust::fill(torsionInfoVecs.reducedIds.begin(), torsionInfoVecs.reducedIds.end(), 0);

		unsigned endKey = thrust::get<0>(
			thrust::reduce_by_key(
				torsionInfoVecs.tempTorIndices.begin(), 
				torsionInfoVecs.tempTorIndices.begin() + 3*generalParams.totalTorsionCount,
			thrust::make_zip_iterator(
				thrust::make_tuple(
					torsionInfoVecs.forceX.begin(),
					torsionInfoVecs.forceY.begin(),
					torsionInfoVecs.forceZ.begin())),
			torsionInfoVecs.reducedIds.begin(),
			thrust::make_zip_iterator(
				thrust::make_tuple(
					torsionInfoVecs.tempForceX.begin(),
					torsionInfoVecs.tempForceY.begin(),
					torsionInfoVecs.tempForceZ.begin())),
			thrust::equal_to<unsigned>(), CVec3Add())) - torsionInfoVecs.reducedIds.begin();//binary_pred, binary_op

		cudaThreadSynchronize();
		
		//std::cout<<"endkey: "<< endKey << std::endl;
		//std::cout<<"totalTorsion: "<< generalParams.totalTorsionCount << std::endl;

		thrust::for_each(
			thrust::make_zip_iterator(//1st begin
				thrust::make_tuple(
					torsionInfoVecs.reducedIds.begin(),
					torsionInfoVecs.tempForceX.begin(),
					torsionInfoVecs.tempForceY.begin(),
					torsionInfoVecs.tempForceZ.begin())),
			thrust::make_zip_iterator(//1st end
				thrust::make_tuple(
					torsionInfoVecs.reducedIds.begin(),
					torsionInfoVecs.tempForceX.begin(),
					torsionInfoVecs.tempForceY.begin(),
					torsionInfoVecs.tempForceZ.begin())) + endKey,
			AddTorsionForceFunctor(
				generalParams.maxNodeCount,
				thrust::raw_pointer_cast(nodeInfoVecs.nodeForceX.data()),
				thrust::raw_pointer_cast(nodeInfoVecs.nodeForceY.data()),
				thrust::raw_pointer_cast(nodeInfoVecs.nodeForceZ.data())));

	}

	
}