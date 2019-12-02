/*
* TorsionSolveOnDevice.h
*
*	Created 11/4/2017 SRB
*/
#pragma once



#ifndef TORSIONSOLVEONDEVICE_H_
#define TORSIONSOLVEONDEVICE_H_

#include <thrust/version.h>
#include <thrust/sort.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <thrust/for_each.h>
#include <thrust/transform.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/tuple.h>
#include <thrust/device_vector.h>
#include <thrust/device_ptr.h>
#include <thrust/remove.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/copy.h>
#include <thrust/pair.h>
#include <stdint.h>
#include <cmath>
//for atomic add
//#include "cuda_runtime.h"
#include <stdlib.h>

//#include "support/common.h"


#include "glm/glm.hpp"
 /*
#if !defined(__CUDA_ARCH__) || __CUDA_ARCH__ >= 600
#else
__device__ double atomicAdd(double* a, double b) { return b; }
#endif	*/


typedef thrust::tuple<unsigned, unsigned, unsigned, unsigned, double> Tuuuud;
typedef thrust::tuple<unsigned, double, double, double> Tuddd;
typedef thrust::tuple<double, double, double> CVec3;
typedef thrust::tuple<unsigned, bool> Tub;
struct NodeInfoVecs;
struct TorsionInfoVecs;
struct GeneralParams;

void TorsionSolveOnDevice(
	NodeInfoVecs& nodeInfoVecs,
	TorsionInfoVecs& torsionInfoVecs,
	GeneralParams& generalParams);

__device__
inline CVec3 cross_product(CVec3 v1, CVec3 v2) {
	return thrust::make_tuple(thrust::get<1>(v1)*thrust::get<2>(v2) - thrust::get<2>(v1)*thrust::get<1>(v2),
		-(thrust::get<0>(v1)*thrust::get<2>(v2) - thrust::get<2>(v1)*thrust::get<0>(v2)),
		thrust::get<0>(v1)*thrust::get<1>(v2) - thrust::get<1>(v1)*thrust::get<0>(v2));
};

__device__
inline double dot_product(CVec3 v1, CVec3 v2) {
	return (thrust::get<0>(v1)*thrust::get<0>(v2) +
		thrust::get<1>(v1)*thrust::get<1>(v2) +
		thrust::get<2>(v1)*thrust::get<2>(v2));
}

struct CVec3Add : public thrust::binary_function<CVec3, CVec3, CVec3> {
	__host__ __device__ 
		CVec3 operator()(const CVec3 &vec1, const CVec3 &vec2) {
		return thrust::make_tuple(thrust::get<0>(vec1) + thrust::get<0>(vec2),
			thrust::get<1>(vec1) + thrust::get<1>(vec2),
			thrust::get<2>(vec1) + thrust::get<2>(vec2));
	}
};

struct AddTorsionForceFunctor {
	unsigned maxNodeCount;
	double* forceXAddr;
	double* forceYAddr;
	double* forceZAddr;

	__host__ __device__
	//
		AddTorsionForceFunctor(
				unsigned& _maxNodeCount,
				double* _forceXAddr,
				double* _forceYAddr,
				double* _forceZAddr) :
			maxNodeCount(_maxNodeCount),
			forceXAddr(_forceXAddr),
			forceYAddr(_forceYAddr),
			forceZAddr(_forceZAddr) {}

	__device__
	void operator() (const Tuddd& u1d3) {
			unsigned idToAssign = thrust::get<0>(u1d3);
			
			if (idToAssign < maxNodeCount) {
				if (!isnan(thrust::get<1>(u1d3)) && !isnan(thrust::get<2>(u1d3)) && !isnan(thrust::get<3>(u1d3))) {

					forceXAddr[idToAssign] += thrust::get<1>(u1d3);
					forceYAddr[idToAssign] += thrust::get<2>(u1d3);
					forceZAddr[idToAssign] += thrust::get<3>(u1d3);
				}
			}

	}

};

struct TorsionFunctor {
	double * locXAddr;
	double * locYAddr;
	double * locZAddr;
	double * forceXAddr;
	double * forceYAddr;
	double * forceZAddr;
	bool* isFixedVector;
	double torsionStiffness;
	unsigned maxNodeCount;
	unsigned totalTorsionCount;
	const double PI;

	__host__ __device__
		//
		TorsionFunctor(
			double* _locXAddr,
			double* _locYAddr,
			double* _locZAddr,
			double* _forceXAddr,
			double* _forceYAddr,
			double* _forceZAddr,
			bool* _isFixedVector,
			double& _torsionStiffness,
			unsigned& _maxNodeCount,
			unsigned& _totalTorsionCount,
			const double& _PI) :

		locXAddr(_locXAddr),
		locYAddr(_locYAddr),
		locZAddr(_locZAddr),
		forceXAddr(_forceXAddr),
		forceYAddr(_forceYAddr),
		forceZAddr(_forceZAddr),
		isFixedVector(_isFixedVector),
		torsionStiffness(_torsionStiffness),
		maxNodeCount(_maxNodeCount),
		totalTorsionCount(_totalTorsionCount),
		PI(_PI) {}

	__device__
		//maybe hand in id and vector of neighbors???
		void operator() (const Tuuuud& u4d1) {
		unsigned counter = thrust::get<0>(u4d1);
		unsigned indexLeft = thrust::get<1>(u4d1);
		unsigned indexCenter = thrust::get<2>(u4d1);
		unsigned indexRight = thrust::get<3>(u4d1);
		double angleZero = thrust::get<4>(u4d1);


		unsigned indexLeftApply = counter;
		unsigned indexCenterApply = counter + totalTorsionCount;
		unsigned indexRightApply = counter + 2*totalTorsionCount;

		double distLCX = locXAddr[indexLeft] - locXAddr[indexCenter];
		double distLCY = locYAddr[indexLeft] - locYAddr[indexCenter];
		double distLCZ = locZAddr[indexLeft] - locZAddr[indexCenter];
		double distRCX = locXAddr[indexRight] - locXAddr[indexCenter];
		double distRCY = locYAddr[indexRight] - locYAddr[indexCenter];
		double distRCZ = locZAddr[indexRight] - locZAddr[indexCenter];

		//lengths between left & center, right & center
		double lenLC = sqrt(distLCX*distLCX + distLCY*distLCY + distLCZ*distLCZ);
		double lenRC = sqrt(distRCX*distRCX + distRCY*distRCY + distRCZ*distRCZ);


		//normalized dot product
		double cosTheta = (distLCX*distRCX + distLCY*distRCY + distLCZ*distRCZ) / (lenLC * lenRC);
		
		//rounding errors
		if (cosTheta < -1.0) {
			cosTheta = -1.0;	
		}
		else if (cosTheta > 1.0) {
			cosTheta = 1.0;
		}

		//should be zero to pi.
		double currentAngle = acos(cosTheta);


		double magForce = -1 * torsionStiffness * (angleZero - currentAngle);

		//calculate force vectors
		//careful here. If the vectors are parallel, axb = 0;
		CVec3 axb = cross_product(thrust::make_tuple(distLCX, distLCY, distLCZ),
			thrust::make_tuple(distRCX, distRCY, distRCZ));


		CVec3 c1 = cross_product(axb, thrust::make_tuple(distLCX, distLCY, distLCZ));
		CVec3 c2 = cross_product(thrust::make_tuple(distRCX, distRCY, distRCZ), axb);

		double divisorA = sqrt(dot_product(c1, c1));
		double divisorB = sqrt(dot_product(c2, c2));

			CVec3 a = thrust::make_tuple(magForce * thrust::get<0>(c1) / divisorA,
				magForce * thrust::get<1>(c1) / divisorA,
				magForce * thrust::get<2>(c1) / divisorA);

			//calculate force on right node
			CVec3 b = thrust::make_tuple(magForce*thrust::get<0>(c2) / divisorB,
				magForce*thrust::get<1>(c2) / divisorB,
				magForce*thrust::get<2>(c2) / divisorB);

		if (isfinite(thrust::get<0>(a))
			&& isfinite(thrust::get<1>(a))
			&& isfinite(thrust::get<2>(a)) 
			&& isfinite(thrust::get<0>(b))
			&& isfinite(thrust::get<1>(b)) 
			&& isfinite(thrust::get<2>(b))) {

			//calculate force on left node


			//apply force to nodes left, right, and center
			//node id indexCenter must be moved up neighborsCount multiples of itself.

			//notice we write force before applying it. We thus rewrite and do not add to the current value.
		if (!isFixedVector[indexCenter]){
			forceXAddr[indexCenterApply] = -thrust::get<0>(a) - thrust::get<0>(b);
			forceYAddr[indexCenterApply] = -thrust::get<1>(a) - thrust::get<1>(b);
			forceZAddr[indexCenterApply] = -thrust::get<2>(a) - thrust::get<2>(b);
		}
		else {
			forceXAddr[indexCenterApply] = 0;
			forceYAddr[indexCenterApply] = 0;
			forceZAddr[indexCenterApply] = 0;
		}

		if (!isFixedVector[indexLeft]) {
			forceXAddr[indexLeftApply] = thrust::get<0>(a);
			forceYAddr[indexLeftApply] = thrust::get<1>(a);
			forceZAddr[indexLeftApply] = thrust::get<2>(a);
		}
		else {
			forceXAddr[indexLeftApply] = 0;
			forceYAddr[indexLeftApply] = 0;
			forceZAddr[indexLeftApply] = 0;
		}

		if (!isFixedVector[indexRight]){
			forceXAddr[indexRightApply] = thrust::get<0>(b);
			forceYAddr[indexRightApply] = thrust::get<1>(b);
			forceZAddr[indexRightApply] = thrust::get<2>(b);
		}
		else {
			forceXAddr[indexRightApply] = 0;
			forceYAddr[indexRightApply] = 0;
			forceZAddr[indexRightApply] = 0;
		}



		}
		else  {
			forceXAddr[indexRightApply] = 0;
			forceYAddr[indexRightApply] = 0;
			forceZAddr[indexRightApply] = 0;
			forceXAddr[indexCenterApply] = 0;
			forceYAddr[indexCenterApply] = 0;
			forceZAddr[indexCenterApply] = 0;
			forceXAddr[indexLeftApply] = 0;
			forceYAddr[indexLeftApply] = 0;
			forceZAddr[indexLeftApply] = 0;
		}

	}
};


	

								
#endif /*TORSIONSOLVEONDEVICE_H_*/