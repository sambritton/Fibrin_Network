	
#include "BucketSchemeOnDevice.h"
#include "NodeSystemDevice.h"

	
typedef thrust::device_vector< unsigned >                IntVector;
typedef IntVector::iterator                         IntIterator;
typedef thrust::tuple< IntIterator, IntIterator >   IntIteratorTuple;
typedef thrust::zip_iterator< IntIteratorTuple >    ZipIterator;
 
void initDimensionBucketScheme(
	NodeInfoVecs& nodeInfoVecs,
	DomainParams& domainParams,
	AuxVecs& auxVecs,
	GeneralParams& generalParams,
	DPDParticleVariables& dpdParticleVariables,
	CompressionParams& compressionParams) {
	
	domainParams.minX = (*(thrust::min_element(nodeInfoVecs.nodeLocX.begin(), nodeInfoVecs.nodeLocX.end())));// - 1.0;
	domainParams.maxX = (*(thrust::max_element(nodeInfoVecs.nodeLocX.begin(), nodeInfoVecs.nodeLocX.end())));// + 1.0;
	domainParams.minY = (*(thrust::min_element(nodeInfoVecs.nodeLocY.begin(), nodeInfoVecs.nodeLocY.end())));// - 1.0;
	domainParams.maxY = (*(thrust::max_element(nodeInfoVecs.nodeLocY.begin(), nodeInfoVecs.nodeLocY.end())));// + 1.0;
	domainParams.minZ = (*(thrust::min_element(nodeInfoVecs.nodeLocZ.begin(), nodeInfoVecs.nodeLocZ.end())));// - 1.0;
	domainParams.maxZ = (*(thrust::max_element(nodeInfoVecs.nodeLocZ.begin(), nodeInfoVecs.nodeLocZ.end())));// + 1.0;
	

	if (generalParams.iterationCounter == 0) {
		domainParams.originMinX = domainParams.minX;
		domainParams.originMaxX = domainParams.maxX;		
		domainParams.originMinY = domainParams.minY;
		domainParams.originMaxY = domainParams.maxY;
		domainParams.originMinZ = domainParams.minZ;
		domainParams.originMaxZ = domainParams.maxZ;
	}
	//all takes place in quadrant 1 
	
	//std::cout<<"BSOD Zminmax: "<<domainParams.minZ <<" "<< domainParams.maxZ	<< std::endl;
 	//std::cout<<"BSOD Xminmax: "<<domainParams.minX <<" "<< domainParams.maxX	<< std::endl;
 
 	if (generalParams.iterationCounter == 0) {
		domainParams.gridSpacing = (dpdParticleVariables.R_c);
		domainParams.XBucketCount = 1.5 * ceil(domainParams.maxX - domainParams.minX) / domainParams.gridSpacing + 1;
		domainParams.YBucketCount = 1.5 * ceil(domainParams.maxY - domainParams.minY) / domainParams.gridSpacing + 1;
		domainParams.ZBucketCount = ceil(domainParams.originMaxZ * (compressionParams.strainProportion + 1.0) ) / domainParams.gridSpacing + 1;// - domainParams.minZ) / domainParams.gridSpacing + 1;

		domainParams.totalBucketCount = domainParams.XBucketCount * domainParams.YBucketCount * domainParams.ZBucketCount;
		
		if (generalParams.iterationCounter == 0 ) 
			std::cout<<"total bucket count: "<< domainParams.totalBucketCount<<std::endl;


	/*	std::cout<<"maxX: "<< domainParams.maxX<<  std::endl;
		std::cout<<"minZ: "<< domainParams.minZ<<  std::endl;
		std::cout<<"maxZ: "<< domainParams.maxZ<<  std::endl;
		std::cout << "totalBucketsX" << domainParams.XBucketCount << std::endl;
		std::cout << "totalBucketsY" << domainParams.YBucketCount << std::endl;
		std::cout << "totalBucketsZ" << domainParams.ZBucketCount << std::endl;
		*/
		auxVecs.keyBegin.resize(domainParams.totalBucketCount);
		auxVecs.keyEnd.resize(domainParams.totalBucketCount);
	}
	thrust::fill(auxVecs.keyBegin.begin(),auxVecs.keyBegin.end(),0);
	thrust::fill(auxVecs.keyEnd.begin(),auxVecs.keyEnd.end(),0);

}

void extendBucketScheme(
	NodeInfoVecs& nodeInfoVecs,
	DomainParams& domainParams,
	AuxVecs& auxVecs) {
		
	//memory is already allocated. 
	unsigned endIndexExpanded = (auxVecs.endIndexBucketKeys) * 27;

	//test for removing copies. 
	unsigned valuesCount = auxVecs.bucketValues.size();
	thrust::fill(auxVecs.bucketKeysExpanded.begin(),auxVecs.bucketKeysExpanded.end(),0);
	thrust::fill(auxVecs.bucketValuesIncludingNeighbor.begin(),auxVecs.bucketValuesIncludingNeighbor.end(),0);

	

	/*
	* beginning of constant iterator
	*/
	thrust::constant_iterator<unsigned> first(27);
	/** 
	* end of constant iterator.
	* the plus sign only indicate movement of position, not value.
	* e.g. movement is 5 and first iterator is initialized as 9
	* result array is [9,9,9,9,9];
	*/
	thrust::constant_iterator<unsigned> last = first + (auxVecs.endIndexBucketKeys); // this is NOT numerical addition!

	expand(first, last,
		thrust::make_zip_iterator(
			thrust::make_tuple(
				auxVecs.bucketKeys.begin(),
				auxVecs.bucketValues.begin())),
		thrust::make_zip_iterator(
			thrust::make_tuple(
				auxVecs.bucketKeysExpanded.begin(),
				auxVecs.bucketValuesIncludingNeighbor.begin())));
		

	thrust::counting_iterator<unsigned> countingBegin(0);

	thrust::transform(
		thrust::make_zip_iterator(
			thrust::make_tuple(
				auxVecs.bucketKeysExpanded.begin(),
				countingBegin)),												  
		thrust::make_zip_iterator(
			thrust::make_tuple(
				auxVecs.bucketKeysExpanded.begin(),
				countingBegin)) + endIndexExpanded,
		thrust::make_zip_iterator(
			thrust::make_tuple(
				auxVecs.bucketKeysExpanded.begin(),
				countingBegin)),
		NeighborFunctor(
			domainParams.XBucketCount,
			domainParams.YBucketCount,
			domainParams.ZBucketCount));

/*		std::cout<<"after transform " << std::endl;
		for (unsigned i = 0; i<auxVecs.bucketValuesIncludingNeighbor.size(); i++ ) {
			std::cout<< " " << auxVecs.bucketValuesIncludingNeighbor[i] <<" ";			
			
		}
		std::cout<<" " << std::endl;
		for (unsigned i = 0; i<auxVecs.bucketValuesIncludingNeighbor.size(); i++ ) {
			std::cout<< " " << auxVecs.bucketKeysExpanded[i] <<" ";			
			
		}
		
		std::cout<<" " << std::endl;*/

	unsigned numberOfOutOfRange = thrust::count_if(auxVecs.bucketKeysExpanded.begin(),
		auxVecs.bucketKeysExpanded.end(), is_greater_than(domainParams.totalBucketCount) );
	unsigned numberInsideRange = endIndexExpanded - numberOfOutOfRange;

	//unsigned endIndexSearch = endIndexExpanded - numberOfOutOfRange;

	thrust::sort_by_key(auxVecs.bucketKeysExpanded.begin(),
		auxVecs.bucketKeysExpanded.begin() + endIndexExpanded,
		auxVecs.bucketValuesIncludingNeighbor.begin());

//	std::cout<<"out of range: " << numberOfOutOfRange << std::endl;
	
//	std::cout<<"in of range: " << numberInsideRange << std::endl;

	auxVecs.bucketKeysExpanded.erase(
			auxVecs.bucketKeysExpanded.begin() + numberInsideRange,
			auxVecs.bucketKeysExpanded.end());

	auxVecs.bucketValuesIncludingNeighbor.erase(
			auxVecs.bucketValuesIncludingNeighbor.begin() + numberInsideRange,
			auxVecs.bucketValuesIncludingNeighbor.end());

	/*std::cout<<"after erase " << std::endl;
		for (unsigned i = 0; i<auxVecs.bucketValuesIncludingNeighbor.size(); i++ ) {
			std::cout<< " " << auxVecs.bucketValuesIncludingNeighbor[i] <<" ";			
			
		}
		std::cout<<" " << std::endl;
		for (unsigned i = 0; i<auxVecs.bucketValuesIncludingNeighbor.size(); i++ ) {
			std::cout<< " " << auxVecs.bucketKeysExpanded[i] <<" ";			
			
		}
		
		std::cout<<" " << std::endl; */
	//now we only want to keep unique tuples
/*	ZipIterator newEnd = thrust::unique( 
		thrust::make_zip_iterator(
			thrust::make_tuple(
				auxVecs.bucketKeysExpanded.begin(),
				auxVecs.bucketValuesIncludingNeighbor.begin())),
		thrust::make_zip_iterator(
			thrust::make_tuple(
				auxVecs.bucketKeysExpanded.end(),
				auxVecs.bucketValuesIncludingNeighbor.end())),
		tupleEqual());	

		// Trim the vectors
		IntIteratorTuple endTuple = newEnd.get_iterator_tuple();
		auxVecs.bucketKeysExpanded.erase( thrust::get<0>( endTuple ), auxVecs.bucketKeysExpanded.end() );
		auxVecs.bucketValuesIncludingNeighbor.erase( thrust::get<1>( endTuple ), auxVecs.bucketValuesIncludingNeighbor.end() );
		std::cout<<"new size " << auxVecs.bucketKeysExpanded.size()<< std::endl;
*/

//another Trim
/*	ZipIterator newEnd2 = thrust::unique(
		thrust::make_zip_iterator(
			thrust::make_tuple(
				auxVecs.bucketValuesIncludingNeighbor.begin(),
				
				auxVecs.bucketKeysExpanded.begin())),
		thrust::make_zip_iterator(
			thrust::make_tuple(
				auxVecs.bucketValuesIncludingNeighbor.end(),
				auxVecs.bucketKeysExpanded.end())),
		tupleEqual());	

		// Trim the vectors
		IntIteratorTuple endTuple2 = newEnd2.get_iterator_tuple();
		auxVecs.bucketKeysExpanded.erase( thrust::get<0>( endTuple2 ), auxVecs.bucketKeysExpanded.end() );
		auxVecs.bucketValuesIncludingNeighbor.erase( thrust::get<1>( endTuple2 ), auxVecs.bucketValuesIncludingNeighbor.end() );
		



	thrust::sort_by_key(auxVecs.bucketKeysExpanded.begin(),
			auxVecs.bucketKeysExpanded.end(),
			auxVecs.bucketValuesIncludingNeighbor.begin());


*/








/*		std::cout<<" " << std::endl;
		for (unsigned i = 0; i<auxVecs.bucketValuesIncludingNeighbor.size(); i++ ) {
			std::cout<< " " << auxVecs.bucketValuesIncludingNeighbor[i] <<" ";			
			
		}
		std::cout<<" " << std::endl;
		for (unsigned i = 0; i<auxVecs.bucketValuesIncludingNeighbor.size(); i++ ) {
			std::cout<< " " << auxVecs.bucketKeysExpanded[i] <<" ";			
			
		}
		
		std::cout<<" " << std::endl;*/

	
	thrust::counting_iterator<unsigned> search_begin(0);

	thrust::lower_bound(auxVecs.bucketKeysExpanded.begin(),
		auxVecs.bucketKeysExpanded.end(), search_begin,
		search_begin + domainParams.totalBucketCount,
		auxVecs.keyBegin.begin());

	thrust::upper_bound(auxVecs.bucketKeysExpanded.begin(),
		auxVecs.bucketKeysExpanded.end(),search_begin,
		search_begin + domainParams.totalBucketCount, 
		auxVecs.keyEnd.begin()); 

		unsigned idTemp = auxVecs.bucketValues[1];
		unsigned bucket = auxVecs.bucketKeys[1];

      /*  for (unsigned i = auxVecs.keyBegin[bucket]; i < auxVecs.keyEnd[bucket]; i++ ) {
			unsigned id = auxVecs.bucketValuesIncludingNeighbor[i];
            std::cout<< "nbr " << id << std::endl;
			std::cout<< "force idTemp " << nodeInfoVecs.nodeForceX[idTemp]<< " " <<nodeInfoVecs.nodeForceY[idTemp]  << " "<< nodeInfoVecs.nodeForceZ[idTemp]<<std::endl;
			std::cout<< "force idnbr " << nodeInfoVecs.nodeForceX[id]<< " " <<nodeInfoVecs.nodeForceY[id]  << " " <<nodeInfoVecs.nodeForceZ[id]<<std::endl;
			double dist = sqrt(
				(nodeInfoVecs.nodeLocX[idTemp]-nodeInfoVecs.nodeLocX[id])*(nodeInfoVecs.nodeLocX[idTemp]-nodeInfoVecs.nodeLocX[id])
				+(nodeInfoVecs.nodeLocY[idTemp]-nodeInfoVecs.nodeLocY[id])*(nodeInfoVecs.nodeLocY[idTemp]-nodeInfoVecs.nodeLocY[id])
				+(nodeInfoVecs.nodeLocZ[idTemp]-nodeInfoVecs.nodeLocZ[id])*(nodeInfoVecs.nodeLocZ[idTemp]-nodeInfoVecs.nodeLocZ[id]));
			std::cout<<"dist between: " << idTemp<< " and " << id<< " " << dist<< std::endl;
			std::cout<< "loc main " << nodeInfoVecs.nodeLocX[idTemp]<< " " <<nodeInfoVecs.nodeLocY[idTemp]  << " "<< nodeInfoVecs.nodeLocZ[idTemp]<<std::endl;
			std::cout<< "loc idnbr " << nodeInfoVecs.nodeLocX[id]<< " " <<nodeInfoVecs.nodeLocY[id]  << " " <<nodeInfoVecs.nodeLocZ[id]<<std::endl;
		}*/
		
	/*if (nodeInfoVecs.nodeForceX.size() > 1600) {
		idTemp = auxVecs.bucketValues[1600];
		bucket = auxVecs.bucketKeys[1600];

        for (unsigned i = auxVecs.keyBegin[bucket]; i < auxVecs.keyEnd[bucket]; i++ ) {
			unsigned id = auxVecs.bucketValuesIncludingNeighbor[i];
            std::cout<< "nbr " << id <<" in bucket " << auxVecs.bucketKeysExpanded[i]<< std::endl;
			//std::cout<< "force idTemp " << nodeInfoVecs.nodeForceX[idTemp]<< " " <<nodeInfoVecs.nodeForceY[idTemp]  << " "<< nodeInfoVecs.nodeForceZ[idTemp]<<std::endl;
			//std::cout<< "force idnbr " << nodeInfoVecs.nodeForceX[id]<< " " <<nodeInfoVecs.nodeForceY[id]  << " " <<nodeInfoVecs.nodeForceZ[id]<<std::endl;
			double dist = sqrt(
				(nodeInfoVecs.nodeLocX[idTemp]-nodeInfoVecs.nodeLocX[id])*(nodeInfoVecs.nodeLocX[idTemp]-nodeInfoVecs.nodeLocX[id])
				+(nodeInfoVecs.nodeLocY[idTemp]-nodeInfoVecs.nodeLocY[id])*(nodeInfoVecs.nodeLocY[idTemp]-nodeInfoVecs.nodeLocY[id])
				+(nodeInfoVecs.nodeLocZ[idTemp]-nodeInfoVecs.nodeLocZ[id])*(nodeInfoVecs.nodeLocZ[idTemp]-nodeInfoVecs.nodeLocZ[id]));
			std::cout<<"dist between: " << idTemp<< " and " << id<< " " << dist<< std::endl;
			std::cout<< "loc main " << nodeInfoVecs.nodeLocX[idTemp]<< " " <<nodeInfoVecs.nodeLocY[idTemp]  << " "<< nodeInfoVecs.nodeLocZ[idTemp]<<std::endl;
			std::cout<< "loc idnbr " << nodeInfoVecs.nodeLocX[id]<< " " <<nodeInfoVecs.nodeLocY[id]  << " " <<nodeInfoVecs.nodeLocZ[id]<<std::endl;
		}
	}*/
}	   

 
void buildBucketScheme(
	NodeInfoVecs& nodeInfoVecs,
	DomainParams& domainParams,
	AuxVecs& auxVecs,
	GeneralParams& generalParams,
	DPDParticleVariables& dpdParticleVariables) {


	thrust::counting_iterator<unsigned> indexBucketBegin(0);
	// takes counting iterator and coordinates
	// return tuple of keys and values
	// transform the points to their bucket indices
	thrust::for_each(
		thrust::make_zip_iterator(
			thrust::make_tuple(
				nodeInfoVecs.nodeLocX.begin(),
				nodeInfoVecs.nodeLocY.begin(),
				nodeInfoVecs.nodeLocZ.begin(),
				indexBucketBegin)),
		thrust::make_zip_iterator(
			thrust::make_tuple(
				nodeInfoVecs.nodeLocX.begin(),
				nodeInfoVecs.nodeLocY.begin(),
				nodeInfoVecs.nodeLocZ.begin(),
				indexBucketBegin)) + generalParams.maxNodeCount + dpdParticleVariables.particleCount,			    
		BucketIndexer(domainParams.minX, domainParams.maxX, domainParams.minY,
			domainParams.maxY, domainParams.minZ, domainParams.maxZ,
			domainParams.gridSpacing,
			thrust::raw_pointer_cast(auxVecs.bucketKeys.data()),
			thrust::raw_pointer_cast(auxVecs.bucketValues.data())));

	// sort the points by their bucket index
	
/*	thrust::sort_by_key(auxVecs.bucketKeys.begin(),
		auxVecs.bucketKeys.begin() + generalParams.maxNodeCount + dpdParticleVariables.particleCount,
		auxVecs.bucketValues.begin());

	unsigned numberOutOfRange = thrust::count(auxVecs.bucketKeys.begin(),
			auxVecs.bucketKeys.begin() + dpdParticleVariables.particleCount + generalParams.maxNodeCount, ULONG_MAX);*/

//test sorting by node instaed of bucket index
thrust::sort_by_key(auxVecs.bucketValues.begin(),
		auxVecs.bucketValues.begin() + generalParams.maxNodeCount + dpdParticleVariables.particleCount,
		auxVecs.bucketKeys.begin());
	unsigned numberOutOfRange = thrust::count(auxVecs.bucketKeys.begin(),
			auxVecs.bucketKeys.begin() + dpdParticleVariables.particleCount + generalParams.maxNodeCount, ULONG_MAX);

	auxVecs.endIndexBucketKeys =  dpdParticleVariables.particleCount + generalParams.maxNodeCount - numberOutOfRange;
 
										
	// for those nodes that are inactive, key value of UINT_MAX will be returned.
	// we need to removed those keys along with their values.
	//int numberOfOutOfRange = thrust::count(auxVecs.bucketKeys.begin(),
	//auxVecs.bucketKeys.begin() + totalActiveNodes, UINT_MAX);
	//endIndx_M = totalActiveNodes - numberOfOutOfRange;

};


