#include "NodeSystemDevice.h"
#include "IncrementExternalForceOnDevice.h"


void IncrementExternalForceOnDevice(
	NodeInfoVecs& nodeInfoVecs,
	GeneralParams& generalParams,
	CompressionParams& compressionParams,
	DomainParams& domainParams) {

		//WARNING: index only goes up to original nodes.
		//calculate averages

		
		thrust::counting_iterator<unsigned> indexBeginB(0);
		thrust::counting_iterator<unsigned> indexBeginC(0);

		compressionParams.averageUpperStrain = (thrust::transform_reduce(
			thrust::make_zip_iterator(
				thrust::make_tuple(
					indexBeginB,
					nodeInfoVecs.nodeUpperChoiceForStrain.begin(),
					nodeInfoVecs.nodeLocZ.begin())),
			thrust::make_zip_iterator(
				thrust::make_tuple(
					indexBeginB,
					nodeInfoVecs.nodeUpperChoiceForStrain.begin(),
					nodeInfoVecs.nodeLocZ.begin())) + generalParams.maxNodeCount,
			StrainFunctor(generalParams.maxNodeCount, compressionParams.originalNetworkLength),
				0.0,
			thrust::plus<double>())) / generalParams.numUpperStrainNodes;	
			
			compressionParams.averageLowerStrain = (thrust::transform_reduce(
				thrust::make_zip_iterator(
					thrust::make_tuple(
						indexBeginC,
						nodeInfoVecs.nodeLowerChoiceForStrain.begin(),
						nodeInfoVecs.nodeLocZ.begin())),
				thrust::make_zip_iterator(
					thrust::make_tuple(
						indexBeginC,
						nodeInfoVecs.nodeLowerChoiceForStrain.begin(),
						nodeInfoVecs.nodeLocZ.begin())) + generalParams.maxNodeCount,
				StrainFunctor(generalParams.maxNodeCount, compressionParams.originalNetworkLength),
					0.0,
				thrust::plus<double>())) / generalParams.numLowerStrainNodes;
	if (generalParams.iterationCounter == 1) {
		compressionParams.originAverageUpperStrain = compressionParams.averageUpperStrain;
		compressionParams.originAverageLowerStrain = compressionParams.averageLowerStrain;
	}
						
	

	//Apply External Force
	thrust::counting_iterator<unsigned> indexBeginA(0);

	thrust::for_each(
		thrust::make_zip_iterator( 
			thrust::make_tuple(
				indexBeginA,
				nodeInfoVecs.nodeLocZ.begin(),
				nodeInfoVecs.isNodeFixed.begin(),
				nodeInfoVecs.nodeUpperChoiceForStrain.begin(),
				nodeInfoVecs.nodeLowerChoiceForStrain.begin())),
		thrust::make_zip_iterator(
			thrust::make_tuple(
				indexBeginA,
				nodeInfoVecs.nodeLocZ.begin(),
				nodeInfoVecs.isNodeFixed.begin(),
				nodeInfoVecs.nodeUpperChoiceForStrain.begin(),
				nodeInfoVecs.nodeLowerChoiceForStrain.begin())) + generalParams.maxNodeCount,
		IncrementFunctor(
			thrust::raw_pointer_cast(nodeInfoVecs.isNodeFixed.data()),
			thrust::raw_pointer_cast(nodeInfoVecs.nodeForceX.data()),
			thrust::raw_pointer_cast(nodeInfoVecs.nodeForceY.data()),
			thrust::raw_pointer_cast(nodeInfoVecs.nodeForceZ.data()),

			generalParams.magnitudeForce,
			compressionParams.originalNetworkLength,
			compressionParams.strainProportion,
			compressionParams.averageLowerStrain,
			compressionParams.averageUpperStrain));

};

