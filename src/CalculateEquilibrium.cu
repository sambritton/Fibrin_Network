/*
* CalculateEquilibrium.cu
*
* Created on 11/5/2017
* 		Author: SRB
*/
#include "NodeSystemDevice.h"
#include "CalculateEquilibrium.h"


typedef thrust::tuple<double, double, double, double, double, double> CVec6;

double CalculateEquilibrium(
	NodeInfoVecs& nodeInfoVecs,
	GeneralParams& generalParams) {


	thrust::counting_iterator<uint32_t> indexBegin(0);
	thrust::counting_iterator<uint32_t> indexEnd(generalParams.maxNodeCount);

	return thrust::transform_reduce(indexBegin, indexEnd,
			EquilibriumFunctor(
				thrust::raw_pointer_cast(nodeInfoVecs.prevNodeLocX.data()),
				thrust::raw_pointer_cast(nodeInfoVecs.prevNodeLocY.data()),
				thrust::raw_pointer_cast(nodeInfoVecs.prevNodeLocZ.data()),
				thrust::raw_pointer_cast(nodeInfoVecs.prevNodeVelX.data()),
				thrust::raw_pointer_cast(nodeInfoVecs.prevNodeVelZ.data()),
				thrust::raw_pointer_cast(nodeInfoVecs.prevNodeVelY.data()),
				thrust::raw_pointer_cast(nodeInfoVecs.nodeLocX.data()),
				thrust::raw_pointer_cast(nodeInfoVecs.nodeLocY.data()),
				thrust::raw_pointer_cast(nodeInfoVecs.nodeLocZ.data())),
			0.0, thrust::plus<double>());
	
		thrust::copy(nodeInfoVecs.nodeVelX.begin(),nodeInfoVecs.nodeVelX.begin() + generalParams.maxNodeCount, nodeInfoVecs.prevNodeVelX.begin());
		thrust::copy(nodeInfoVecs.nodeVelY.begin(),nodeInfoVecs.nodeVelY.begin() + generalParams.maxNodeCount, nodeInfoVecs.prevNodeVelY.begin());
		thrust::copy(nodeInfoVecs.nodeVelZ.begin(),nodeInfoVecs.nodeVelZ.begin() + generalParams.maxNodeCount, nodeInfoVecs.prevNodeVelZ.begin());
}