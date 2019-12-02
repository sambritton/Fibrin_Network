
#ifndef CALCULATEEQUILIBRIUM_H_
#define CALCULATEEQUILIBRIUM_H_

#include <thrust/sort.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <thrust/for_each.h>
#include <thrust/transform.h>

#include <thrust/transform_reduce.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/tuple.h>
#include <thrust/device_vector.h>
#include <thrust/device_ptr.h>
#include <thrust/remove.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/copy.h>
#include <thrust/pair.h>
#include <stdint.h>

typedef thrust::tuple<double, double, double, double, double, double> CVec6;
struct NodeInfoVecs;
struct GeneralParams;


double CalculateEquilibrium(
	NodeInfoVecs&,
	GeneralParams&);
 



struct EquilibriumFunctor { 
		double* prevNodeLocXAddr;
		double* prevNodeLocYAddr;
		double* prevNodeLocZAddr;
		double* prevNodeVelXAddr;
		double* prevNodeVelYAddr;
		double* prevNodeVelZAddr;
		double* nodeLocXAddr;
		double* nodeLocYAddr;
		double* nodeLocZAddr;

	__host__ __device__
	
	EquilibriumFunctor(
		double* _prevNodeLocXAddr,
		double* _prevNodeLocYAddr, 
		double* _prevNodeLocZAddr,
		double* _prevNodeVelXAddr,
		double* _prevNodeVelYAddr,
		double* _prevNodeVelZAddr,
		double* _nodeLocXAddr, 
		double* _nodeLocYAddr, 
		double* _nodeLocZAddr):
		prevNodeLocXAddr(_prevNodeLocXAddr),
		prevNodeLocYAddr(_prevNodeLocYAddr),
		prevNodeLocZAddr(_prevNodeLocZAddr),
		prevNodeVelXAddr(_prevNodeVelXAddr),
		prevNodeVelYAddr(_prevNodeVelYAddr),
		prevNodeVelZAddr(_prevNodeVelZAddr),
		nodeLocXAddr(_nodeLocXAddr),
		nodeLocYAddr(_nodeLocYAddr),
		nodeLocZAddr(_nodeLocZAddr) {}


	__device__
		
	double operator() (const unsigned& id) {

		double result = sqrt(
			(prevNodeLocXAddr[id] - nodeLocXAddr[id]) * (prevNodeLocXAddr[id] - nodeLocXAddr[id]) +
			(prevNodeLocYAddr[id] - nodeLocYAddr[id]) * (prevNodeLocYAddr[id] - nodeLocYAddr[id]) +
			(prevNodeLocZAddr[id] - nodeLocZAddr[id]) * (prevNodeLocZAddr[id] - nodeLocZAddr[id]));

		//now reset the previous location to be the new one.
		prevNodeLocXAddr[id] = nodeLocXAddr[id];
		prevNodeLocYAddr[id] = nodeLocYAddr[id];
		prevNodeLocZAddr[id] = nodeLocZAddr[id];
		if (isnan(result))
			return 0;
		else
			return result;
	}
	
};







#endif /*CALCULATEEQUILIBRIUM_H_*/