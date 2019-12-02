
#include "LinkNodesOnDevice.h"
#include "NodeSystemDevice.h"
#include <thrust/execution_policy.h>


void LinkNodesOnDevice(
	NodeInfoVecs& nodeInfoVecs,
	WLCInfoVecs& wlcInfoVecs,
	AuxVecs& auxVecs,
	TorsionInfoVecs & torsionInfoVecs,
	GeneralParams& generalParams) {
 
			//std::cout<<"begin"<<std::endl;	
			/*if (generalParams.magnitudeForce == 52) {
				
				std::cout<<"pre link"<<std::endl;
				unsigned maxLinked = *( thrust::max_element(wlcInfoVecs.currentEdgeCountVector.begin(),wlcInfoVecs.currentEdgeCountVector.end()) );
				std::cout<<"maximum neighbors currently: "<< maxLinked<<std::endl;
				std::cout<<"maximum neighbors: "<< generalParams.maxNeighborCount <<std::endl;	
			}
			if (generalParams.magnitudeForce == 30) {
				unsigned start = 880 * generalParams.maxNeighborCount;
							unsigned end = start + generalParams.maxNeighborCount;
							std::cout<<"id's"<<std::endl;
							for (unsigned i = start; i < end; i++) {
								std::cout<<" "<< wlcInfoVecs.globalNeighbors[i]<< " ";
							}
							std::cout<<" "<<std::endl;
							std::cout<<"lengths's"<<std::endl;
							for (unsigned i = start; i < end; i++) {
								std::cout<<" "<< wlcInfoVecs.lengthZero[i] << " ";
							}
							std::cout<<" "<<std::endl;
			}*/   
	try { thrust::for_each(  
				thrust::make_zip_iterator(
					thrust::make_tuple(
						auxVecs.bucketKeys.begin(), 
						auxVecs.bucketValues.begin())),
				thrust::make_zip_iterator(
					thrust::make_tuple(
						auxVecs.bucketKeys.begin(),
						auxVecs.bucketValues.begin())) + generalParams.maxNodeCount,
				
				LinkNodesFunctor(
					thrust::raw_pointer_cast(nodeInfoVecs.nodeLocX.data()),
					thrust::raw_pointer_cast(nodeInfoVecs.nodeLocY.data()),
					thrust::raw_pointer_cast(nodeInfoVecs.nodeLocZ.data()),
					thrust::raw_pointer_cast(wlcInfoVecs.globalNeighbors.data()),
					thrust::raw_pointer_cast(wlcInfoVecs.currentNodeEdgeCountVector.data()),
					thrust::raw_pointer_cast(wlcInfoVecs.lengthZero.data()),
					thrust::raw_pointer_cast(auxVecs.bucketKeys.data()),
					thrust::raw_pointer_cast(auxVecs.bucketValues.data()),
					thrust::raw_pointer_cast(auxVecs.bucketValuesIncludingNeighbor.data()),
					thrust::raw_pointer_cast(auxVecs.keyBegin.data()),
					thrust::raw_pointer_cast(auxVecs.keyEnd.data()), 
					generalParams.fiberDiameter,
					generalParams.maxNeighborCount,
					generalParams.maxNodeCount, 
					generalParams.maxLinksPerIteration, 
					thrust::raw_pointer_cast(nodeInfoVecs.idEdgesMadeTemp.data()) ) );
  
			cudaThreadSynchronize();
 
			/*thrust::counting_iterator<unsigned> startDeLinkIter(0);
			thrust::counting_iterator<unsigned> endDeLinkIter(generalParams.maxNodeCount);
		 
			thrust::for_each( 
				startDeLinkIter, endDeLinkIter, 
				DeLinkCopiesFunctor(
					thrust::raw_pointer_cast(wlcInfoVecs.globalNeighbors.data()),
					thrust::raw_pointer_cast(wlcInfoVecs.lengthZero.data()),
					thrust::raw_pointer_cast(wlcInfoVecs.currentEdgeCountVector.data()),
					generalParams.maxNeighborCount, 
					generalParams.maxNodeCount));
				
			cudaThreadSynchronize();*/
		
		}
		catch(thrust::system_error &e) { std::cerr << "Error linking functor: " << e.what() << std::endl; exit(-1); }


 
		unsigned numAdded = thrust::count_if( nodeInfoVecs.idEdgesMadeTemp.begin(), nodeInfoVecs.idEdgesMadeTemp.end(), isNotEqualZero() );
		//std::cout<<"totalToraion: " << generalParams.totalTorsionCount << std::endl;
		if (numAdded != 0) {
			//be aware that if an edge was made, it will appear twice, but if it was made once and removed, it will appear once.
			//sort in increasing order. Then when we hit zero, we are done.  
			try { thrust::sort(nodeInfoVecs.idEdgesMadeTemp.begin(), nodeInfoVecs.idEdgesMadeTemp.end(),thrust::greater<unsigned>() ); }
			catch(thrust::system_error &e){ std::cerr << "Error link sorting: " << e.what() << std::endl; exit(-1); }
	
			try {	
				unsigned idLast = nodeInfoVecs.idEdgesMadeTemp[0];
				unsigned count = 0;
				for (unsigned i = 1; i<nodeInfoVecs.idEdgesMadeTemp.size(); i++) {
					//add extra edges and preferred lengths. Notice the lower and upper must be added since each imparts force to one single node and 
					//not the neighboring node to the edge. This is b/c edges are solved per node and not per edge
					unsigned id = nodeInfoVecs.idEdgesMadeTemp[i];
					
					if (id == 0){
						break; 
					}
		
					if (id == idLast) { 
						//then id has shown up twice and can be added
						count +=1; 
					}  
					else { 
						count = 0;
					}
		
					if ((id !=0) && (count > 0) ) {
						//count edges 
						//std::cout<<"placing id: "<< id<<std::endl;
					
						nodeInfoVecs.idEdgesMadeHost.push_back(id);
						unsigned idLeft = id % generalParams.maxNodeCount;
						unsigned idRight = id / generalParams.maxNodeCount;

						//std::cout<< "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
						//std::cout<< " new Id linked: " << idLeft << " " << idRight << std::endl;
						//std::cout<< "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;

						nodeInfoVecs.deviceEdgeLeft[generalParams.currentEdgeCount] = (idLeft);
						nodeInfoVecs.deviceEdgeRight[generalParams.currentEdgeCount] = (idRight);
						generalParams.currentEdgeCount += 1;
/*
						//now that we added the id's that got linked. let's add the torsion springs. 
						unsigned id_left_1=0; //iterates through neighbors of id_center_1
						unsigned id_center_1 = idLeft;
						unsigned id_right_1 = idRight;

						unsigned begin_center_1 = generalParams.maxNeighborCount * id_center_1;
						unsigned end_center_1 = generalParams.maxNeighborCount + begin_center_1;
						for (unsigned nbr = begin_center_1; nbr < end_center_1; nbr++) {
							unsigned possible_nbr = wlcInfoVecs.globalNeighbors[nbr];
							if ((possible_nbr < generalParams.maxNodeCount) && 
								(possible_nbr != id_center_1) && 
								(possible_nbr != id_right_1))  {

								id_left_1 = possible_nbr;
								//std::cout<< "Adding torsion: " << id_left_1 << " " << id_center_1 << " " << id_right_1 << " " << std::endl;
								
								//then id_left_1, id_center_1 and id_right_1 are the id's. 
								//calculate cosTheta for the spring. 
								double distLCX = nodeInfoVecs.nodeLocX[id_left_1] - nodeInfoVecs.nodeLocX[id_center_1];
								double distLCY = nodeInfoVecs.nodeLocY[id_left_1] - nodeInfoVecs.nodeLocY[id_center_1];
								double distLCZ = nodeInfoVecs.nodeLocZ[id_left_1] - nodeInfoVecs.nodeLocZ[id_center_1];
								double distRCX = nodeInfoVecs.nodeLocX[id_right_1] - nodeInfoVecs.nodeLocX[id_center_1];
								double distRCY = nodeInfoVecs.nodeLocY[id_right_1] - nodeInfoVecs.nodeLocY[id_center_1];
								double distRCZ = nodeInfoVecs.nodeLocZ[id_right_1] - nodeInfoVecs.nodeLocZ[id_center_1];
						
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

								// now store everything and add to totolTorsionCount
								//std::cout<< "Adding torsion: " << id_left_1 << " " << id_center_1 << " " << id_right_1 << " " << cosTheta << std::endl;

								torsionInfoVecs.leftIndex[generalParams.totalTorsionCount] = id_left_1;
								torsionInfoVecs.centerIndex[generalParams.totalTorsionCount] = id_center_1;
								torsionInfoVecs.rightIndex[generalParams.totalTorsionCount] = id_right_1;
								
								torsionInfoVecs.angleZero[generalParams.totalTorsionCount] = cosTheta;
								generalParams.totalTorsionCount += 1;


							}
						}


						//now we finished torsion for the left point, add the torsion for the right. 
						id_left_1 = idLeft;//iterates through neighbors of id_center_1
						id_center_1 = idRight;
						id_right_1=0;
						//std::cout<< " changing center to idRight: " << id_center_1 << std::endl;
						
						begin_center_1 = generalParams.maxNeighborCount * id_center_1;
						end_center_1 = generalParams.maxNeighborCount + begin_center_1;
						for (unsigned nbr = begin_center_1; nbr < end_center_1; nbr++) {
							unsigned possible_nbr = wlcInfoVecs.globalNeighbors[nbr];
							if ((possible_nbr < generalParams.maxNodeCount) && 
								(possible_nbr != id_center_1) && 
								(possible_nbr != id_left_1))  {

								id_right_1 = possible_nbr;
								//std::cout<< "Adding torsion: " << id_left_1 << " " << id_center_1 << " " << id_right_1 << " " << std::endl;
								//then id_left_1, id_center_1 and id_right_1 are the id's. 
								//calculate cosTheta for the spring. 
								double distLCX = nodeInfoVecs.nodeLocX[id_left_1] - nodeInfoVecs.nodeLocX[id_center_1];
								double distLCY = nodeInfoVecs.nodeLocY[id_left_1] - nodeInfoVecs.nodeLocY[id_center_1];
								double distLCZ = nodeInfoVecs.nodeLocZ[id_left_1] - nodeInfoVecs.nodeLocZ[id_center_1];
								double distRCX = nodeInfoVecs.nodeLocX[id_right_1] - nodeInfoVecs.nodeLocX[id_center_1];
								double distRCY = nodeInfoVecs.nodeLocY[id_right_1] - nodeInfoVecs.nodeLocY[id_center_1];
								double distRCZ = nodeInfoVecs.nodeLocZ[id_right_1] - nodeInfoVecs.nodeLocZ[id_center_1];
						
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

								// now store everything and add to totolTorsionCount
								//std::cout<< "Adding torsion: " << id_left_1 << " " << id_center_1 << " " << id_right_1 << " " << cosTheta << std::endl;
								torsionInfoVecs.leftIndex[generalParams.totalTorsionCount] = id_left_1;
								torsionInfoVecs.centerIndex[generalParams.totalTorsionCount] = id_center_1;
								torsionInfoVecs.rightIndex[generalParams.totalTorsionCount] = id_right_1;
								
								torsionInfoVecs.angleZero[generalParams.totalTorsionCount] = cosTheta;
								generalParams.totalTorsionCount += 1;
								
								//std::cout<<"space left: " << torsionInfoVecs.angleZero.size() << std::endl; 
								//std::cout<<"totalToraion: " << generalParams.totalTorsionCount << std::endl;

							}
						}*/
					} 
					
					idLast = id;//set last id to current
				} 
				//std::cout<<"post link: " << std::endl;
				//std::cout<<"totalToraion: " << generalParams.totalTorsionCount << std::endl;
			}
			catch(thrust::system_error &e){ std::cerr << "Error link loop: " << e.what() << std::endl; exit(-1); }
			catch (std::exception& e) { std::cout <<"error link loop "<<  e.what(); }
		}																					
		

		thrust::fill(nodeInfoVecs.idEdgesMadeTemp.begin(), 
				nodeInfoVecs.idEdgesMadeTemp.end(), 0);
		
};
