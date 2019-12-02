#ifndef INCREMENTEXTERNALFORCEONDEVICE_H_
#define INCREMENTEXTERNALFORCEONDEVICE_H_

#pragma once

#include <thrust/extrema.h>
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
#include <thrust/iterator/counting_iterator.h>
#include <thrust/copy.h>
#include <thrust/pair.h>
#include <stdint.h>

struct NodeInfoVecs;
struct GeneralParams;
struct CompressionParams;
struct DomainParams;

typedef thrust::tuple<unsigned, double, double, double> Tuddd;

typedef thrust::tuple<unsigned, double, bool, bool, bool> Tudbbb;
typedef thrust::tuple<unsigned, bool, double> Tubd;
typedef thrust::tuple<unsigned, bool, bool> Tubb;
typedef thrust::tuple<double, double, double> CVec3;
typedef thrust::tuple<unsigned, bool> Tub;

void IncrementExternalForceOnDevice(
	NodeInfoVecs& nodeInfoVecs,
	GeneralParams& generalParams,
	CompressionParams& compressionParams,
	DomainParams& domainParams);

struct StrainFunctor {
	unsigned maxNodeCount;
	double maxNetworkLength;
	
	__host__ __device__
	StrainFunctor(
		unsigned& _maxNodeCount,
		double& _maxNetworkLength):
		maxNodeCount(_maxNodeCount),
		maxNetworkLength(_maxNetworkLength) {}
		
	__device__
	double operator() (const Tubd& u1b1d1) {
		unsigned id = thrust::get<0>(u1b1d1);
		bool isStrainNode = thrust::get<1>(u1b1d1);

		if (isStrainNode) {
			double zpos = thrust::get<2>(u1b1d1);
			return zpos;
		}
		else {
			return (0.0);
		}
	}
	
};


struct IncrementFunctor {
	bool* isNodeFixedAddr;
	double* forceXAddr;
	double* forceYAddr;
	double* forceZAddr;
	double magForce;
	double originalNetworkLength;
	double strainProportion;
	double averageLowerStrain;
	double averageUpperStrain;

	__host__ __device__
		IncrementFunctor(
			bool* _isNodeFixedAddr,
			double*	_forceXAddr,
			double*	_forceYAddr,
			double*	_forceZAddr,
			double& _magForce,
			double& _originalNetworkLength,
			double& _strainProportion,
			double& _averageLowerStrain,
			double& _averageUpperStrain) :
		isNodeFixedAddr(_isNodeFixedAddr),
		forceXAddr(_forceXAddr),
		forceYAddr(_forceYAddr),
		forceZAddr(_forceZAddr),
		magForce(_magForce),
		originalNetworkLength(_originalNetworkLength),
		strainProportion(_strainProportion),
		averageLowerStrain(_averageLowerStrain),
		averageUpperStrain(_averageUpperStrain) {}

	__device__
	void operator()(const Tudbbb& u1d1b3) {
		
		unsigned id = thrust::get<0>(u1d1b3);
		double locZ = thrust::get<1>(u1d1b3);
		bool isFixed = thrust::get<2>(u1d1b3);
		bool isUpperStrainNode = thrust::get<3>(u1d1b3);
		bool isLowerStrainNode = thrust::get<4>(u1d1b3);
		
		//pull top
		if ((!isFixed) && (isUpperStrainNode)) {
			//if not fixed, we can apply force unless it is too far away, then we fix it. 
			if (locZ > originalNetworkLength * (0.1 + strainProportion) ) {
				isNodeFixedAddr[id] = true;
			}

			double upperDiff = abs(locZ - averageUpperStrain);

			//only apply force if within 10 of the average. 
			if (upperDiff < 5) {
				double dirX = 0.0;//tForceX / normForce;
				double dirY = 0.0;//tForceY / normForce;
				double dirZ = 1.0;
				forceXAddr[id] = dirX * (magForce);
				forceYAddr[id] = dirY * (magForce);
				forceZAddr[id] = dirZ * (magForce);
			}
		}

		
		//pull bottom	
		else if ((!isFixed) && (isLowerStrainNode) ) {
			//safety fix
			if (locZ < -(originalNetworkLength * (0.1 + strainProportion))) {
				isNodeFixedAddr[id] = true;
			}

			
			double lowerDiff = abs(locZ - averageLowerStrain);

			if (lowerDiff < 5) {
				double dirX = 0.0;//tForceX / normForce;
				double dirY = 0.0;//tForceY / normForce;
				double dirZ = -1.0;
				forceXAddr[id] = dirX * (magForce);
				forceYAddr[id] = dirY * (magForce);
				forceZAddr[id] = dirZ * (magForce);
			}

		}
	}
};


#endif /*INCREMENTEXTERNALFORCEONDEVICE_H_*/