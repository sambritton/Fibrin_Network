
#ifndef WLCSOLVEONDEVICE_H_
#define WLCSOLVEONDEVICE_H_


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



/*
the structure of lengthZero_index is 
0  1  2  3 
4  5  6  7 
8  9  10 11
12 13 14 15 for a 4 node system. 
index/4 = row,
index%4 = col. If you apply force to column node always or row node always then 
each thread will apply opposing forces to springs. 
if you decide to apply force to column instead of rows, you'll need sign change
LengthZero_value is symmetric, so values line up correctly.
*/
//
struct NodeInfoVecs;
struct WLCInfoVecs;
struct GeneralParams; 
struct DomainParams;
struct CompressionParams;

typedef thrust::tuple<double, double, double> CVec3;
typedef thrust::tuple<unsigned, unsigned, double> Tuud;
typedef thrust::tuple<unsigned, unsigned, bool> Tuub;
typedef thrust::tuple<unsigned, unsigned, CVec3> TuuCV;
typedef thrust::tuple<unsigned, unsigned, double, CVec3> TuudCV;
typedef thrust::tuple<unsigned, double> Tud;
typedef thrust::tuple<unsigned, unsigned> Tuu;
typedef thrust::tuple<unsigned, unsigned,unsigned> Tuuu;
typedef thrust::tuple<double, double> Tdd;

typedef thrust::tuple<CVec3, CVec3, CVec3> Mat_3x3;

void WLCSolveOnDevice(
	NodeInfoVecs& nodeInfoVecs,
	WLCInfoVecs& wlcInfoVecs,
	GeneralParams& generalParams);
void GetStrainParameters(NodeInfoVecs& nodeInfoVecs,
	WLCInfoVecs& wlcInfoVecs,  
	GeneralParams& generalParams,
	DomainParams& domainParams);

	
struct WLCfunctor {
	double* locXAddr;
	double* locYAddr;
	double* locZAddr;
	double* forceXAddr;
	double* forceYAddr;
	double* forceZAddr;

	double Kb; //convert to nN and microns
	double PLengthMon;
	double CLM;
	double Temp;
	unsigned maxNeighborCount;
	unsigned maxNodeCount;

	double* lenZero;
	unsigned* edgeCountVec;
	unsigned* globalNeighbors;
	unsigned* springDivisionCount;
	unsigned* numOriginalNeighborsVec;

	__host__ __device__

		WLCfunctor(
			double* _locXAddr, 
			double* _locYAddr, 
			double* _locZAddr,
			double* _forceXAddr, 
			double* _forceYAddr, 
			double* _forceZAddr, 

			double& _Kb, 
			double& _PLengthMon, 
			double& _CLM, 
			double& _Temp,
			unsigned& _maxNeighborCount,
			unsigned& _maxNodeCount,

			double* _lenZero,
			unsigned* _springDivisionCount,
			unsigned* _globalNeighbors,
			unsigned* _edgeCountVec,
			unsigned* _numOriginalNeighborsVec) :

		locXAddr(_locXAddr),
		locYAddr(_locYAddr),
		locZAddr(_locZAddr),
		forceXAddr(_forceXAddr),
		forceYAddr(_forceYAddr),
		forceZAddr(_forceZAddr),

		Kb(_Kb), 
		PLengthMon(_PLengthMon), 
		CLM(_CLM), 
		Temp(_Temp),
		maxNeighborCount(_maxNeighborCount),
		maxNodeCount(_maxNodeCount),

		lenZero(_lenZero),
		springDivisionCount(_springDivisionCount),
		globalNeighbors(_globalNeighbors),
		edgeCountVec(_edgeCountVec),
		numOriginalNeighborsVec(_numOriginalNeighborsVec) {}
 
	__device__
	void operator()(const Tuub& u2b1) {
		//idA represents row.
		unsigned idA = thrust::get<0>(u2b1);
		unsigned numOriginalNeighbors = numOriginalNeighborsVec[idA];

		//variable nbrcountMultiplier == 1 if main node, or number of dividing edges if subnode.
		unsigned nbrCountMultiplierA = thrust::get<1>(u2b1);
		unsigned nbrCountMultiplierB;

		bool isFixed = thrust::get<2>(u2b1);
		double sumForceX = 0;
		double sumForceY = 0;
		double sumForceZ = 0;
		//double numberOfExtendedEdges = 0.0;

		if (!isFixed) {
			//only apply force if not fixed. 
			//id ranges from 0 to maxNodeCount, to apply force to node id, we find those connected 
			//to it. If lenZero_index[i]/maxNodeCount == id, then it is in the same row and id 
			//is linked to the node with id lenZero_index[i]%maxNodeCount since idMat ranges from 0
			//to maxNodeCount^2.

			
			unsigned beginIndex = idA * maxNeighborCount;
			unsigned endIndex = beginIndex + maxNeighborCount;
			bool isAddedLink = false;

			for (unsigned i = beginIndex; i < endIndex; i++) {//currentSpringCount is the length of index and value vectors
				unsigned idB = globalNeighbors[i];//look through possible neighbors. May contain ULONG_MAX
				if (idB < maxNodeCount){
					//now idB is not ULONG_MAX
					if (numOriginalNeighbors < (i % maxNeighborCount) ) {
						//then the index is past the origial edges and the edge is linked. 
						isAddedLink = true;
					}
							

					if (nbrCountMultiplierA == 1) {
						//then the node idA is a main node, and we choose nbCount from idB. 
						//if idB is also a main, node, nbrCount stays at 1, else it, was not 1, 
						//and we do not want to change it. 
						//nbrCountMultiplierA = springDivisionCount[idB];

					}
					nbrCountMultiplierB = springDivisionCount[idB];
					//in the case of main node, the two multipliers are equal, but in the case of a crosslinker
					//we have difference. To maintain stability, we average. 
					//the unsigned will be 
					double aveCount = max(nbrCountMultiplierA, nbrCountMultiplierB);
					if ((i % maxNeighborCount) >= numOriginalNeighbors) {
						//then the neighbor is linked, so we do not scale for series of springs.
						//e.g. Say a node k has two initial neighbors, the third neighbor is located in slot k*maxNeighborCount+2, if we modulo, we are left with 2. So we must use >= and not >
						aveCount = 1.0;
					}

					double lengthZero = lenZero[i];//
					if (lengthZero > 0) {

						double posXA_XB = locXAddr[idB] - locXAddr[idA];
						double posYA_YB = locYAddr[idB] - locYAddr[idA];
						double posZA_ZB = locZAddr[idB] - locZAddr[idA];
		
						double currentLength = sqrt(
							(posXA_XB) * (posXA_XB)+
							(posYA_YB) * (posYA_YB)+
							(posZA_ZB) * (posZA_ZB));
						double strain = ((currentLength - lengthZero) / lengthZero);

						if ( strain > 0 ) {
							
							//numberOfExtendedEdges += 0.5;
							//Warning: No ave count used!
							//This decision is based on data. 
							double dL_norm = strain / ( CLM);//CLM is unitless since it was already normalized. 

							double magForce = (1100*(Kb*Temp) / PLengthMon) * ( 0.25 * pow(1.0 - dL_norm, -2.0) - 0.25 + dL_norm);
							//double magForce = ((Kb*Temp) / PLength)*(0.25*(1 / ((1-dL_norm)*(1-dL_norm))) - 0.25 + dL_norm);
							//double magForce = 1* (lengthZero - currentLength);
							if ((isAddedLink) && (strain > 0.0) ) {
								magForce = 2*magForce;
							}
							//if I update magForceX then memory errors
							double magForceX = (posXA_XB / currentLength) * magForce;
							double magForceY = (posYA_YB / currentLength) * magForce;
							double magForceZ = (posZA_ZB / currentLength) * magForce;
							
							sumForceX += magForceX;
							sumForceY += magForceY;
							sumForceZ += magForceZ;
						}
						else {

							//only difference here:
							double normalizerForNegStrain = 1;
							double dL_norm = strain / ( CLM);//CLM is unitless since it was already normalized. 
							if (isAddedLink ) {
								dL_norm = abs(dL_norm);
							}

							double magForce = (1100*(Kb*Temp) / PLengthMon)*(0.25*pow(1.0 - dL_norm, -2.0) - 0.25 + dL_norm);
							
							if ((isAddedLink) && (strain < 0.0) ) {
								magForce = -2*magForce;
							}

							//double magForce = ((Kb*Temp) / PLength)*(0.25*(1 / ((1-dL_norm)*(1-dL_norm))) - 0.25 + dL_norm);
							//double magForce = 1* (lengthZero - currentLength);
			

							double magForceX = (posXA_XB / currentLength) * magForce;
							double magForceY = (posYA_YB / currentLength) * magForce;
							double magForceZ = (posZA_ZB / currentLength) * magForce;
							
							sumForceX += normalizerForNegStrain * magForceX;
							sumForceY += normalizerForNegStrain * magForceY;
							sumForceZ += normalizerForNegStrain * magForceZ;
							
						}
					}
				}
			}
			
			
			if (isfinite(sumForceX))
				forceXAddr[idA] += sumForceX;
			
			if (isfinite(sumForceY))
				forceYAddr[idA] += sumForceY;
			
			if (isfinite(sumForceY))
				forceZAddr[idA] += sumForceZ;	
		}

	}
};

struct CalculateStrainParamsFunctor {
	unsigned originLinkCount;
	unsigned originEdgeCount;
	unsigned originNodeCount;
	unsigned maxNodeCount;
	unsigned maxNeighborCount;

	double* locXAddr;
	double* locYAddr;
	double* locZAddr;

	unsigned* originNbrVec;
	unsigned* currentNbrVec;
	unsigned* globalNeighbors;
	double* lenZero;

	__host__ __device__
	CalculateStrainParamsFunctor(
		unsigned& _originLinkCount,
		unsigned& _originEdgeCount,
		unsigned& _originNodeCount,
		unsigned& _maxNodeCount,
		unsigned& _maxNeighborCount,	

		double* _locXAddr, 
		double* _locYAddr, 
		double* _locZAddr,

		unsigned* _originNbrVec,
		unsigned* _currentNbrVec,
		unsigned* _globalNeighbors,
		double* _lenZero) :
		originLinkCount(_originLinkCount),
		originEdgeCount(_originEdgeCount),
		originNodeCount(_originNodeCount),
		maxNodeCount(_maxNodeCount),
		maxNeighborCount(_maxNeighborCount),

		locXAddr(_locXAddr), 
		locYAddr(_locYAddr), 
		locZAddr(_locZAddr),
		
		originNbrVec(_originNbrVec),
		currentNbrVec(_currentNbrVec),
		globalNeighbors(_globalNeighbors),
		lenZero(_lenZero){}
 
	__device__
	Tdd operator() (const Tuu& u2) {

		unsigned idL = thrust::get<0>(u2);//use IDL as row
		unsigned idR = thrust::get<1>(u2);//use IDR as col

		//identify matrix location
		unsigned idMat;
		for (unsigned i = idL; i < maxNeighborCount; i++) {
			unsigned idCol = globalNeighbors[i];//represents nbr id
			if (idR == idCol) {
				idMat = idCol;
				break;
			}
		}		

		double lengthZero = lenZero[idMat];

		unsigned originNbrCount = originNbrVec[idL];//number of original neighbors on idRow
		unsigned currentNbrCount = currentNbrVec[idL];

		double edgeStrain = 0.0;
		double edgeAlignment = 0.0;
		
		if ((lengthZero > 0.0) ) {
			double posXA_XB = locXAddr[idL] - locXAddr[idR];
			double posYA_YB = locYAddr[idL] - locYAddr[idR];
			double posZA_ZB = locZAddr[idL] - locZAddr[idR];
			double currentLength = sqrt(
				(posXA_XB) * (posXA_XB) +
				(posYA_YB) * (posYA_YB) +
				(posZA_ZB) * (posZA_ZB));
			edgeStrain = ((currentLength - lengthZero) / lengthZero);
			edgeAlignment = abs(1.0 * posZA_ZB / currentLength);
		}


		
		
		return thrust::make_tuple(edgeStrain, edgeAlignment);
	}
};




#endif /*WLCSOLVEONDEVICE*/
