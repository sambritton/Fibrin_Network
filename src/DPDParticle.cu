#include "DPDParticle.h"
#include "NodeSystemDevice.h"
  

//updates force and velocity in that order
void DPDParticleForce( 
    NodeInfoVecs& nodeInfoVecs,
    DPDParticleVariables& dpdParticleVariables,
    AuxVecs& auxVecs,
    GeneralParams& generalParams) {  
 
    unsigned interactionCount = auxVecs.bucketValuesIncludingNeighbor.size();
      
		
    unsigned _seed = rand();
    thrust::device_vector<double> gaussianData;
    gaussianData.resize(interactionCount); // 
    thrust::counting_iterator<unsigned> index_sequence_begin(_seed);
    thrust::transform(thrust::device, index_sequence_begin, index_sequence_begin + (interactionCount), 
        gaussianData.begin(), psrngen(-1.0f, 1.0f));
     
    std::cout<<"gaussianData size: " << gaussianData.size() << std::endl;
    //std::cout<< "adding dpd forces" << std::endl;
   
     
    //adds force to dpd particles and 
    thrust::for_each(  
    		thrust::make_zip_iterator(
    			thrust::make_tuple(
    				auxVecs.bucketKeys.begin(),  
    				auxVecs.bucketValues.begin())),
    		thrust::make_zip_iterator(
    			thrust::make_tuple(
    				auxVecs.bucketKeys.begin(),
    				auxVecs.bucketValues.begin())) + (generalParams.maxNodeCount + dpdParticleVariables.particleCount),
            ParticleForceFunctor( 
                generalParams.dtTemp, 
                dpdParticleVariables.R_c,
                dpdParticleVariables.alpha,
                dpdParticleVariables.particleCount,
                generalParams.kB, 
                generalParams.temperature,
                generalParams.viscousDamp,
                generalParams.maxNodeCount,
                thrust::raw_pointer_cast(gaussianData.data()), 
                thrust::raw_pointer_cast(nodeInfoVecs.nodeLocX.data()),
    			thrust::raw_pointer_cast(nodeInfoVecs.nodeLocY.data()),
    			thrust::raw_pointer_cast(nodeInfoVecs.nodeLocZ.data()),
                thrust::raw_pointer_cast(nodeInfoVecs.nodeVelX.data()),
    			thrust::raw_pointer_cast(nodeInfoVecs.nodeVelY.data()),
    			thrust::raw_pointer_cast(nodeInfoVecs.nodeVelZ.data()),
                thrust::raw_pointer_cast(nodeInfoVecs.nodeForceX.data()),
    			thrust::raw_pointer_cast(nodeInfoVecs.nodeForceY.data()),
    			thrust::raw_pointer_cast(nodeInfoVecs.nodeForceZ.data()),
                thrust::raw_pointer_cast(nodeInfoVecs.prevNodeForceX.data()),
    			thrust::raw_pointer_cast(nodeInfoVecs.prevNodeForceY.data()),
    			thrust::raw_pointer_cast(nodeInfoVecs.prevNodeForceZ.data()),
    			thrust::raw_pointer_cast(auxVecs.bucketKeys.data()),
    			thrust::raw_pointer_cast(auxVecs.bucketValues.data()),
    			thrust::raw_pointer_cast(auxVecs.bucketValuesIncludingNeighbor.data()),
    			thrust::raw_pointer_cast(auxVecs.keyBegin.data()),
    			thrust::raw_pointer_cast(auxVecs.keyEnd.data())));


        //force swap to clear
        gaussianData.clear();
        
        gaussianData.shrink_to_fit();
        //thrust::device_vector<double>().swap(gaussianData);
        std::cout<< " in dpd force count : " <<nodeInfoVecs.nodeForceX.size()<< std::endl;
 

};

//updates position and velocity corrector.
void DPDParticleAdvance(
    NodeInfoVecs& nodeInfoVecs,
    DPDParticleVariables& dpdParticleVariables,
    GeneralParams& generalParams) {

 
		thrust::transform(
				thrust::make_zip_iterator(
					thrust::make_tuple(
						nodeInfoVecs.nodeLocX.begin(),
						nodeInfoVecs.nodeLocY.begin(),
						nodeInfoVecs.nodeLocZ.begin(),
						nodeInfoVecs.nodeVelX.begin(),
						nodeInfoVecs.nodeVelY.begin(),
						nodeInfoVecs.nodeVelZ.begin())) + generalParams.maxNodeCount,
				thrust::make_zip_iterator(
					thrust::make_tuple(
						nodeInfoVecs.nodeLocX.begin(),
						nodeInfoVecs.nodeLocY.begin(),
						nodeInfoVecs.nodeLocZ.begin(),
						nodeInfoVecs.nodeVelX.begin(),
						nodeInfoVecs.nodeVelY.begin(),
						nodeInfoVecs.nodeVelZ.begin())) + (generalParams.maxNodeCount + dpdParticleVariables.particleCount),
				//second vector begin
				thrust::make_zip_iterator(
					thrust::make_tuple(
						nodeInfoVecs.nodeForceX.begin(),
						nodeInfoVecs.nodeForceY.begin(),
						nodeInfoVecs.nodeForceZ.begin())) + generalParams.maxNodeCount,
				//save result in third vector 
				thrust::make_zip_iterator(
					thrust::make_tuple(
						nodeInfoVecs.nodeLocX.begin(),
						nodeInfoVecs.nodeLocY.begin(), 
						nodeInfoVecs.nodeLocZ.begin(),								 
						nodeInfoVecs.nodeVelX.begin(),
						nodeInfoVecs.nodeVelY.begin(),
						nodeInfoVecs.nodeVelZ.begin())) + generalParams.maxNodeCount,
				SaxpyFunctorDPD(generalParams.dtTemp,
					dpdParticleVariables.particleMass));
    

    //after advancing positions, we have used the force from step n, so copy the force from 
    //current step for use next iteration.
	thrust::copy(nodeInfoVecs.nodeForceX.begin() + generalParams.maxNodeCount, 
        nodeInfoVecs.nodeForceX.end(), 
        nodeInfoVecs.prevNodeForceX.begin() + generalParams.maxNodeCount);
	thrust::copy(nodeInfoVecs.nodeForceY.begin() + generalParams.maxNodeCount, 
        nodeInfoVecs.nodeForceY.end(), 
        nodeInfoVecs.prevNodeForceY.begin() + generalParams.maxNodeCount);
	thrust::copy(nodeInfoVecs.nodeForceZ.begin() + generalParams.maxNodeCount, 
        nodeInfoVecs.nodeForceZ.end(), 
        nodeInfoVecs.prevNodeForceZ.begin() + generalParams.maxNodeCount);


};

