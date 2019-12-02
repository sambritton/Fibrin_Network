#ifndef NODESYSTEMIMPLDEVICEFUNCTORS_H_
#define NODESYSTEMIMPLDEVICEFUNCTORS_H_

#include <memory>
#include <thrust/extrema.h>
#include <thrust/sort.h>
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
#include <thrust/execution_policy.h>
#include <thrust/pair.h>
#include <thrust/unique.h>
#include <stdint.h>
 
#pragma once
typedef thrust::tuple<unsigned, unsigned> Tuu;
typedef thrust::tuple<unsigned, double> Tud;
typedef thrust::tuple<bool, double> Tbd;
typedef thrust::tuple<unsigned, unsigned, unsigned> Tuuu;
typedef thrust::tuple<unsigned, unsigned, unsigned,unsigned> Tuuuu;

typedef thrust::tuple<double, double, double, unsigned> Tdddu;

typedef thrust::tuple<bool, double, double, double, double, double, double> BoolCVec6;
typedef thrust::tuple<double, double, double, double, double, double> CVec6;
typedef thrust::tuple<double, double, double, double> CVec4;
typedef thrust::tuple<double, double, double> CVec3;
typedef thrust::tuple<double, double> CVec2;
 
struct TorsionAngleFunctor {
	double* locXAddr;
	double* locYAddr;
	double* locZAddr;
	__host__ __device__
		TorsionAngleFunctor(
			double* _locXAddr,
			double* _locYAddr,
			double* _locZAddr) :
			locXAddr(_locXAddr),
			locYAddr(_locYAddr),
			locZAddr(_locZAddr) {}

	__device__
	double operator() (const Tuuu &u3) {
		unsigned indexLeft = thrust::get<0>(u3);
		unsigned indexCenter = thrust::get<1>(u3);
		unsigned indexRight = thrust::get<2>(u3);

		double distLCX = locXAddr[indexLeft] - locXAddr[indexCenter];
		double distLCY = locYAddr[indexLeft] - locYAddr[indexCenter];
		double distLCZ = locZAddr[indexLeft] - locZAddr[indexCenter];
		double distRCX = locXAddr[indexRight] - locXAddr[indexCenter];
		double distRCY = locYAddr[indexRight] - locYAddr[indexCenter];
		double distRCZ = locZAddr[indexRight] - locZAddr[indexCenter];
//		CVec3 distLC = 

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

		//double currentAngle = acos(cosTheta);

		return acos(cosTheta);
	}
};



//used to calculate strain
struct AveStrainFunctor {
  __host__ __device__ 
  double operator() (const Tbd &b1d1) {
	  bool isStrainNode = thrust::get<0>(b1d1);
	  if (isStrainNode) 
		  return thrust::get<1>(b1d1);
	  else 
    	return 0.0; 
  }
};




struct NormFunctor {
	
		__host__ __device__
		double operator() (const CVec3& vec) {
		//divide force by fiber cross section to get stress
		double result = sqrt(thrust::get<0>(vec) * thrust::get<0>(vec) +
			thrust::get<1>(vec) * thrust::get<1>(vec) +
			thrust::get<2>(vec) * thrust::get<2>(vec));
		if (isfinite(result))
			return result;
		else
			return 0.0; 

	}
};
 
//returns true if is greater than level
struct IsGreaterThanLevel {
	double limit;

	__host__ __device__ 
		IsGreaterThanLevel(
			double& _limit) : 
			limit(_limit) {}
		
	__device__
		bool operator() (double val) {

			return (val > limit);
		}
}; 

//returns true if is greater than level
struct IsGreaterOrLessThanLevel {
	double upperLimit;
	double lowerLimit;

	__host__ __device__ 
		IsGreaterOrLessThanLevel(
			double& _upperLimit,
			double& _lowerLimit) : 
			upperLimit(_upperLimit),
			lowerLimit(_lowerLimit) {}
		
	__device__
	//replaces value with 1 if returns true 
		bool operator() (double val) {
			if (val > upperLimit) {
				return true;
			}
			else if (val < lowerLimit) {
				return true;
			}
			else {
				return false;
			}
			
		}
}; 

//returns true if less than llevel.
struct IsLessThanLevel {
		double limit;
 
	__host__ __device__ 
		IsLessThanLevel(
			double& _limit) : 
			limit(_limit) {}
		
	__device__
	//replaces value with 1 if returns true 
	bool operator() (double val) {
		return (val <  abs( limit) );//((1-percentPull) * networkLength));
	}
};

struct IsLessThanLevelAndNonZero {
	double limit;
 
	__host__ __device__ 
		IsLessThanLevelAndNonZero(
			double& _limit) : 
			limit(_limit) {}
		
	__device__
	//replaces value with 1 if returns true 
	bool operator() (double val) {
		if (val != 0.0) {
			return (val < abs(limit) );
		}
		else 
			return false;
	}
};

struct IsLessGreaterThanLevelAndNonZero {
	double limit;
 
	__host__ __device__ 
		IsLessGreaterThanLevelAndNonZero (
			double& _limit) : 
			limit(_limit) {}
		
	__device__
	//replaces value with 1 if returns true 
	bool operator() (double val) {
		if (val != 0.0) {
			return ( (val < abs(limit)) && ( val > -abs(limit)) );
		}
		else 
			return false;
	}
};

struct IsEqualToOne {
	 __host__ __device__

	bool operator()(const unsigned& x) {
		if ( x == 1 ) {
			return true;
		}
		else {
			return false;
		}
	}
};

#endif /* NODESYSTEMIMPLDEVICEFUNCTORS_H_*/