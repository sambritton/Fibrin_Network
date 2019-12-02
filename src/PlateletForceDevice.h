
#ifndef PLATELETFORCEDEVICE_H_
#define PLATELETFORCEDEVICE_H_

#include <vector>


typedef thrust::tuple<unsigned, bool> Tub; //tuple holding id of a node and if that node is within reach of the platelet

//look up forward declaration if this doesn't make sense.
struct NodeInfoVecs;
struct WLCInfoVecs;
struct GeneralParams;
struct PlateletParams;
struct PlateletInfoVecs;

__host__ __device__

void PlateletForceOnDevice(
    NodeInfoVecs& nodeInfoVecs,
	WLCInfoVecs& wlcInfoVecs,
	GeneralParams& generalParams,
    PlateletParams& plateletParams,
    PlateletInfoVecs& plateletInfoVecs);

void AdvancePlateletPosition(stuff here);

struct AdvancePosition {

}

//go through and add appropriate entries for input
struct PlateletForceFunctor {
    double* plateletLocXAddr;
	double* plateletLocYAddr;
	double* plateletLocZAddr;
	double* plateletForceXAddr;
	double* plateletForceYAddr;
	double* plateletForceZAddr;
    bool* plateletCanPullAddr;
    double bodyRadius;
    double armDist;
    unsigned numberOfArms;
    unsigned maxPlateletCount;
    double pullingForce; 
    double* nodeLocXAddr;
	double* nodeLocYAddr;
	double* nodeLocZAddr;
	double* nodeForceXAddr;
	double* nodeForceYAddr;
	double* nodeForceZAddr;
    unsigned* bucketValues;
	unsigned* bucketNbrsExp;
	unsigned* keyBegin;
	unsigned* keyEnd;


   __host__ __device__

       PlateletForceFunctor(
            double* _plateletLocXAddr,
            double* _plateletLocYAddr,
            double* _plateletLocZAddr,
            double* _plateletForceXAddr,
            double* _plateletForceYAddr,
            double* _plateletForceZAddr,
            unsigned& _numberOfArms, 
            double& pullingForce,
            double* _nodeLocXAddr,
            double* _nodeLocYAddr,
            double* _nodeLocZAddr,
            double* _nodeForceXAddr,
            double* _nodeForceYAddr,
            double* _nodeForceZAddr) :
        plateletLocXAddr(_plateletLocXAddr),
		plateletLocYAddr(_plateletLocYAddr),
		plateletLocZAddr(_plateletLocZAddr),
		plateletForceXAddr(_plateletForceXAddr),
		plateletForceYAddr(_plateletForceYAddr),
		plateletForceZAddr(_plateletForceZAddr),
        numberOfArms(_numberOfArms), 
        pullingForce(_pullingForce),

        nodeLocXAddr(_nodeLocXAddr),
		nodeLocYAddr(_nodeLocYAddr),
		nodeLocZAddr(_nodeLocZAddr),
		nodeForceXAddr(_nodeForceXAddr),
		nodeForceYAddr(_nodeForceYAddr),
		nodeForceZAddr(_nodeForceZAddr),
        bucketValues(_bucketValues),
		bucketNbrsExp(_bucketNbrsExp),
		keyBegin(_keyBegin),
		keyEnd(_keyEnd), {}


   __device__
   void operator() (const Tuu& u2) {

       //choice 1: functor pulls 4 times
       //choice 2: 
        unsigned id = thrust::get<0>(u1);
        unsigned bucketId = thrust::get<0>(u1);
        
        //beginning and end of attempted interaction network nodes. 
		unsigned beginIndex = keyBegin[bucketId];
		unsigned endIndex = keyEnd[bucketId];
   

        unsigned storageLocation = id * maxPlateletCount;

        double pLocX = plateletLocXAddr[id];
        double pLocY = plateletLocYAddr[id];
        double pLocZ = plateletLocZAddr[id];
        		

        double sumPlateletForceX = 0;   
        double sumPlateletForceY = 0;
        double sumPlateletForceZ = 0;

        double forcePerArm = pullingforce / numberOfArms;

    //Loop through the number of available neighbors for each platelet.
    for(unsigned i=beginIndex; i < endIndex; i++) {

    
        
            
         //Choose a neighbor.
         unsigned pullNode_id = bucketNbrsExp[i];
         bool plateletCanPull = plateletCanPullAddr[pullNode_id];
         //Get position of node
         double vecN_PX = pLocX - nodeLocXAddr[pullNode_id];
         double vecN_PY = pLocY - nodeLocYAddr[pullNode_id];
         double vecN_PZ = pLocZ - nodeLocZAddr[pullNode_id];
         //Calculate distance from platelet to node.
         double dist = sqrt(
             (vecN_PX) * (vecN_PX) +
             (vecN_PY) * (vecN_PY) +
             (vecN_PZ) * (vecN_PZ));
         
         if (dist < bodyRadius) && (plateletCanPull) {
             //then we detach the node from the network and make sure it is never pulled again. 
             plateletCanPullAddr[pullNode_id] = false;

             //loop through global neighbors and set all edges of pullNode_id equal to ULONG_MAX.
             //you'll need to add global neighbors
         }

        //only pull as many as are arms. 
        if (pullCounter < numberOfArms) {

            if ((dist < armDist) && dist > (bodyRadius)) {
                //node only affects platelet position if it is pulled. 
                //Determine direction of force based on positions and multiply magnitude force
                double forceNodeX = (vecP_NX / dist) * (forcePerArm);        
                double forceNodeY = (vecP_NY / dist) * (forcePerArm);
                double forceNodeZ = (vecP_NZ / dist) * (forcePerArm); 

                //count force for platelet.
                sumPlateletForceX += (-1)*forceNodeX;        
                sumPlateletForceY += (-1)*forceNodeY;
                sumPlateletForceZ += (-1)*forceNodeZ;

                //store force in temporary vector. Call reduction later. 
                forceAppliedXAddr[storageLocation + pullCounter] = forceNodeX;
                forceAppliedYAddr[storageLocation + pullCounter] = forceNodeY;
                forceAppliedZAddr[storageLocation + pullCounter] = forceNodeZ;
                idForceAppliedAddr[storageLocation + pullCounter] = pullNode_id;

                pullCounter++

            }
        }


    }
        

    

    plateletForceX[id] = sumPlateletForceX;
    plateletForceY[id] = sumPlateletForceY;
    plateletForceZ[id] = sumPlateletForceZ;


   }
}