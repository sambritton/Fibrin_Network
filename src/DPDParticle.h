#ifndef DPDPARTICLE_H_
#define DPDPARTICLE_H_

#include <thrust/random.h>
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

typedef thrust::tuple<bool, double, double, double, double, double, double> BoolCVec6;
typedef thrust::tuple<double, double, double, double, double, double> CVec6;
typedef thrust::tuple<double, double, double> CVec3;

//#include <thrust/random/default_random_engine.h>

struct NodeInfoVecs;
struct DPDParticleVariables;
struct AuxVecs;
struct GeneralParams;

typedef thrust::tuple<unsigned, unsigned> Tuu;

void DPDParticleForce(
    NodeInfoVecs& nodeInfoVecs,
    DPDParticleVariables& dpdParticleVariables,
    AuxVecs& auxVecs,
    GeneralParams& generalParams);

void DPDParticleAdvance(
    NodeInfoVecs& nodeInfoVecs,
    DPDParticleVariables& dpdParticleVariables,
    GeneralParams& generalParams);

struct psrngen {
    
    double a, b;

    __host__ __device__ psrngen(double _a, double _b) : a(_a), b(_b) {;}
 
    __host__ __device__ double operator()(const unsigned n) const
    {
        thrust::default_random_engine rng(n);
        thrust::uniform_real_distribution<float> dist(a, b);
        rng.discard(n);
        return dist(rng);
    }
 
};

struct ParticleForceFunctor {
    double dt;
    double R_c;
    double alpha;
    unsigned particleCount;
    double kB;
    double temperature;
    double viscosity;
    unsigned maxNodeCount;
    double* gaussianData;
    double* nodeLocXAddr;
	double* nodeLocYAddr;
	double* nodeLocZAddr;
    double* nodeVelXAddr;
    double* nodeVelYAddr;
    double* nodeVelZAddr;
    double* nodeForceXAddr;
    double* nodeForceYAddr;
    double* nodeForceZAddr;
    double* prevNodeForceXAddr;
    double* prevNodeForceYAddr;
    double* prevNodeForceZAddr;
    unsigned* bucketKeys;
	unsigned* bucketValues;
	unsigned* bucketNbrsExp;
	unsigned* keyBegin;
	unsigned* keyEnd;

    __host__ __device__
    ParticleForceFunctor(
        double _dt,
        double _R_c,
        double _alpha,
        unsigned _particleCount,
        double _kB,
        double _temperature,
        double _viscosity,
        unsigned _maxNodeCount,
        double* _gaussianData,
        double* _nodeLocXAddr,
        double* _nodeLocYAddr,
        double* _nodeLocZAddr,
        double* _nodeVelXAddr,
        double* _nodeVelYAddr,
        double* _nodeVelZAddr,
        double* _nodeForceXAddr,
        double* _nodeForceYAddr,
        double* _nodeForceZAddr,
        double* _prevNodeForceXAddr,
        double* _prevNodeForceYAddr,
        double* _prevNodeForceZAddr,
        unsigned* _bucketKeys,
		unsigned* _bucketValues,
		unsigned* _bucketNbrsExp,
		unsigned* _keyBegin,
		unsigned* _keyEnd) :
        dt(_dt),
        R_c(_R_c),
        alpha(_alpha),
        particleCount(_particleCount),
        kB(_kB),
        temperature(_temperature),
        viscosity(_viscosity),
        maxNodeCount(_maxNodeCount),
        gaussianData(_gaussianData),
        nodeLocXAddr(_nodeLocXAddr),
        nodeLocYAddr(_nodeLocYAddr),
        nodeLocZAddr(_nodeLocZAddr),
        nodeVelXAddr(_nodeVelXAddr),
        nodeVelYAddr(_nodeVelYAddr),
        nodeVelZAddr(_nodeVelZAddr),
        nodeForceXAddr(_nodeForceXAddr),
        nodeForceYAddr(_nodeForceYAddr),
        nodeForceZAddr(_nodeForceZAddr),
        prevNodeForceXAddr(_prevNodeForceXAddr),
        prevNodeForceYAddr(_prevNodeForceYAddr),
        prevNodeForceZAddr(_prevNodeForceZAddr),
        bucketKeys(_bucketKeys),
        bucketValues(_bucketValues),
        bucketNbrsExp(_bucketNbrsExp),
        keyBegin(_keyBegin),
        keyEnd(_keyEnd) {}

    __device__
    void operator() (const Tuu& u2) {

        unsigned bucketId = thrust::get<0>(u2);//bucket containing nodeId
		unsigned nodeId = thrust::get<1>(u2);//node to attempt link from.
        if (nodeId < (maxNodeCount + particleCount)) {
           
            //beginning and end of attempted interaction end id's in bucketNbrsExp
		    unsigned beginIndex = keyBegin[bucketId];
		    unsigned endIndex = keyEnd[bucketId];

            //positions
            double xPos = nodeLocXAddr[nodeId];
            double yPos = nodeLocYAddr[nodeId];
            double zPos = nodeLocZAddr[nodeId]; 

            double forceXHolder = 0;
            double forceYHolder = 0;
            double forceZHolder = 0;

            for (unsigned iter = beginIndex; iter < endIndex; iter++ ) {
		    	unsigned candidateId = bucketNbrsExp[iter];//test id

				
                if (candidateId < (maxNodeCount + particleCount) ) {
                //if nodeId and candidateId < maxNodeCount, no interaction.
                //only interaction if one is larger 
				
				//trial for dpd to test only self interaction.
                    if ((maxNodeCount <= nodeId) || (maxNodeCount <= candidateId)) {
                        double gaussNoise = 0.0;
                        //in order to maintain conservation of momentum, the noise must by symmetric
                        //we look up the bucket


                        double xDiff = xPos - nodeLocXAddr[candidateId];
                        double yDiff = yPos - nodeLocYAddr[candidateId];
                        double zDiff = zPos - nodeLocZAddr[candidateId]; 
                        double dist = sqrt(xDiff*xDiff + yDiff*yDiff + zDiff*zDiff);

                        if ((dist < R_c) && (dist > 0)) {

                            //look up gaussian noise if interaciton
                            if (candidateId > nodeId) {
                                //then we use the noise data in bucket of nodeId in the location of 
                                //candidateId. 
                                gaussNoise = gaussianData[iter];
                            }
                            else {
                                //we look up the bucket of candidateId and use the data at the point where nodeId is stored 
                                unsigned bucketOfCandidate = bucketKeys[candidateId];
                                unsigned beginIndexCandidate = keyBegin[bucketOfCandidate];
                                unsigned endIndexCandidate = keyEnd[bucketOfCandidate];
                                for (unsigned iterCand = beginIndexCandidate; iterCand < endIndexCandidate; iterCand++) {
                                    unsigned attemptNodeId = bucketNbrsExp[iterCand];
                                    if (attemptNodeId == nodeId) {
                                        //then we found our data
                                        gaussNoise = gaussianData[iterCand];
                                        break;
                                    }
                                }
                            }

                            double velXDiff = nodeVelXAddr[nodeId] - nodeVelXAddr[candidateId] ;
                            double velYDiff = nodeVelYAddr[nodeId] - nodeVelYAddr[candidateId] ;
                            double velZDiff = nodeVelZAddr[nodeId] - nodeVelZAddr[candidateId] ;

                            double xDir = xDiff / dist;
                            double yDir = yDiff / dist;                
                            double zDir = zDiff / dist;

                            //calculate conservative force
                            double W_C = (1-dist/R_c);

                            forceXHolder += alpha * W_C * xDir;
                            forceYHolder += alpha * W_C * yDir;
                            forceZHolder += alpha * W_C * zDir;

                            //calculate dissipative/frictional force
                            double W_D = (1-dist/R_c)*(1-dist/R_c);
                            //dot produce of normalized direction and velocity
                            double magnitude = velXDiff * xDir + velYDiff * yDir + velZDiff * zDir;

                            forceXHolder += -viscosity * W_D * magnitude * xDir; 
                            forceYHolder += -viscosity * W_D * magnitude * yDir;
                            forceZHolder += -viscosity * W_D * magnitude * zDir;

                            //calculate 
                            double W_R = sqrt(1-dist/R_c);
                            double sigma = sqrt(2*viscosity*kB*temperature);
                            //I need a different noise for each iteration. 
                            //use bucketNbrsExpanded size vector filled with random numbers
//TEST NOT NOI
                            forceXHolder += sigma * W_R*gaussNoise*(1/sqrt(dt))*xDir;
                            forceYHolder += sigma * W_R*gaussNoise*(1/sqrt(dt))*yDir;
                            forceZHolder += sigma * W_R*gaussNoise*(1/sqrt(dt))*zDir;
					
                        }

                    }
                }

            }
          //add forces to nodes. 
            if (isfinite(forceXHolder)) 
                nodeForceXAddr[nodeId] += forceXHolder;
            
            if (isfinite(forceYHolder)) 
                nodeForceYAddr[nodeId] += forceYHolder;
            
            if (isfinite(forceZHolder)) 
                nodeForceZAddr[nodeId] += forceZHolder;
            
            //Now that we have force for step n+1 and n, we can update velocity.
            //for the numerical scheme on the particles. 
            if (nodeId > maxNodeCount) {
                nodeVelXAddr[nodeId] += 0.5*dt*(nodeForceXAddr[nodeId] + prevNodeForceXAddr[nodeId]);
                nodeVelYAddr[nodeId] += 0.5*dt*(nodeForceYAddr[nodeId] + prevNodeForceYAddr[nodeId]);
                nodeVelZAddr[nodeId] += 0.5*dt*(nodeForceZAddr[nodeId] + prevNodeForceZAddr[nodeId]);
            }
        }
    }
};


struct SaxpyFunctorDPD : public thrust::binary_function<CVec6, CVec3, CVec6> {
	double dt;
	double lampda=0.5;
	double mass;
	__host__ __device__
		//
		SaxpyFunctorDPD(double _dt, double _mass) :
		dt(_dt),
		mass(_mass) {}

	__device__
		CVec6 operator()(const CVec6 &p3v3, const CVec3 &f3) {
		//Solving  
		/*
		x' = v
		v' = f/m
		so k1 = [v, a], k2 = [v+adt/2, a], k2 = same etc
		*/

		//only update position and velocity of non fixed nodes. 
		//same as k1Pos
		double velX = thrust::get<3>(p3v3);
		double velY = thrust::get<4>(p3v3);
		double velZ = thrust::get<5>(p3v3);
		
		//same as k1Vel (don't forget v' = F/m)
		//We subtract off a viscous term
		double accX = thrust::get<0>(f3);
		double accY = thrust::get<1>(f3);
		double accZ = thrust::get<2>(f3);		
		
        //result new pos = oldpos + dt*vel + .5*dt^2*acc;
		double xPosResult = thrust::get<0>(p3v3) + dt * velX + 0.5 * dt * dt * accX;
		double yPosResult = thrust::get<1>(p3v3) + dt * velY + 0.5 * dt * dt * accY;
		double zPosResult = thrust::get<2>(p3v3) + dt * velZ + 0.5 * dt * dt * accZ;
		
        //second calculation for velocity update
		//result
		//corrector
		double velXHat = velX + lampda * dt * accX;
		double velYHat = velY + lampda * dt * accY;
		double velZHat = velZ + lampda * dt * accZ;
				

			
		return thrust::make_tuple(xPosResult, yPosResult, zPosResult, velXHat, velYHat, velZHat);
                 
	}                                 

};






#endif /* DPDPARTICLE_H_ */