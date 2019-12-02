#ifndef ADVANCEPOSITIONONDEVICE_H_
#define ADVANCEPOSITIONONDEVICE_H_

#include <thrust/random.h>

typedef thrust::tuple<double, double, double, double, double, double, double> CVec7;
typedef thrust::tuple<unsigned, double, double, double, double, double, double> UCVec6;
typedef thrust::tuple<unsigned, double, double, double> UCVec3;
typedef thrust::tuple<double, double, double, double, double, double> CVec6;

typedef thrust::tuple<double, double, double, double> CVec4;
typedef thrust::tuple<double, double, double> CVec3;

struct NodeInfoVecs;
struct GeneralParams;
struct DPDParticleVariables;

double AdvancePositionOnDevice(
	NodeInfoVecs& nodeInfoVecs,
	GeneralParams& generalParams,
	DPDParticleVariables& dpdParticleVariables);

struct DiffCVec3NormFunctor {
	__host__ __device__

		double operator() (const CVec6& vec) {
		//divide force by fiber cross section to get stress
		return fabs(
			((thrust::get<0>(vec) - thrust::get<3>(vec))) +
			((thrust::get<1>(vec) - thrust::get<4>(vec))) +
			((thrust::get<2>(vec) - thrust::get<5>(vec))));
	}
};

struct psrnormgen {
    
    double a, b;

    __host__ __device__ 
	psrnormgen(
		double _a, 
		double _b) : 
		a(_a), 
		b(_b) {}
 
    __device__ double operator()(const unsigned n) const
    {
        thrust::default_random_engine rng(n);
        thrust::normal_distribution<double> dist(a, b);
        rng.discard(n);
        return dist(rng);
    }
 
};

struct psrunifgen {
    
    double a, b;

    __host__ __device__ 
	psrunifgen(
		double _a, 
		double _b) : 
		a(_a), 
		b(_b) {}
 
    __device__ double operator()(const unsigned n) const
    {
        thrust::default_random_engine rng(n);
        thrust::uniform_real_distribution<double> dist(a, b);
        rng.discard(n);
        return dist(rng);
    }
 
};
//used to advance position and velocity from known force
//Velocity Verlet Like algorithm used : https://arxiv.org/pdf/1212.1244.pdf
struct SaxpyFunctorDim3 : public thrust::binary_function<UCVec3, CVec4, CVec4> {
	double dt;
	double viscosity;
	double temperature;
	double kB;
	double mass;
	bool* isNodeFixedAddr; 
	unsigned maxNodeCount;

	__host__ __device__
		//
		SaxpyFunctorDim3(
			double& _dt, 
			double& _viscosity, 
			double& _temperature,
			double& _kB,
			double& _mass,
			unsigned& _maxNodeCount,
			bool* _isNodeFixedAddr) :
		dt(_dt),
		viscosity(_viscosity),
		temperature(_temperature),
		kB(_kB),
		mass(_mass),
		maxNodeCount(_maxNodeCount),
		isNodeFixedAddr(_isNodeFixedAddr) {}

	__device__
		CVec4 operator()(const UCVec3 &p3, const CVec4 &g1f3) {
		// 
		/*
		x' = v 
		v' = f/m
		so k1 = [v, a], k2 = [v+adt/2, a], k2 = same etc
		*/
		unsigned id = thrust::get<0>(p3);
		bool isFixed =false;//true if fixed, false if movable. 
		if (id < maxNodeCount) {
			isFixed = isNodeFixedAddr[id];
		}

		if (!isFixed) {
			double locX = thrust::get<1>(p3);
			double locY = thrust::get<2>(p3);
			double locZ = thrust::get<3>(p3);

		/*	double velX = thrust::get<4>(p3v3);
			double velY = thrust::get<5>(p3v3);
			double velZ = thrust::get<6>(p3v3);*/
			
			//not using: https://www2.ph.ed.ac.uk/~dmarendu/MVP/MVP03.pdf numeric scheme
			//using: https://arxiv.org/pdf/1212.1244.pdf

/*			double prevAccX = thrust::get<1>(g1pf3f3);
			double prevAccY = thrust::get<2>(g1pf3f3);
			double prevAccZ = thrust::get<3>(g1pf3f3);*/

			//random data
			double gaussianData = thrust::get<0>(g1f3);
			//double b = 1/ (1 + (viscosity * dt / 2));

			//normally you would have x_(n+1) - x(n) / dt = F / eta + F_b / eta, and F_b = sqrt(2*kb*T*eta/dt) after multiplication, additive noise 
			//becomes sqrt(2*kb*t*dt/eta) * N(0,1)
			double noise = sqrt(2.0 * kB* temperature * dt  / viscosity) * gaussianData;



			//same as k1Vel (don't forget v' = F/m)
			//We subtract off a viscous term 
			double accX = (thrust::get<1>(g1f3));//  - viscosity * velX) / mass;
			double accY = (thrust::get<2>(g1f3));//  - viscosity * velY) / mass;
			double accZ = (thrust::get<3>(g1f3));//  - viscosity * velZ) / mass;


			//update positions
			double xLocRes = locX + (dt/viscosity) * (accX) + noise;/*b * dt * velX +
							b * dt * dt * accX * (1 / (2 * mass)) + 
							((b * dt) / (2 * mass)) * noise;*/

			double yLocRes = locY + (dt/viscosity) * (accY) + noise;/*b * dt * velY + 
							b * dt * dt * accY * (1 / (2 * mass)) + 
							((b * dt) / (2 * mass)) * noise;*/

			double zLocRes = locZ + (dt/viscosity) * (accZ) + noise;/*b * dt * velZ + 
							b * dt * dt * accZ * (1 / (2 * mass)) + 
							((b * dt) / (2 * mass)) * noise;*/
			
			double velocity = sqrt((xLocRes - locX) * (xLocRes - locX) + (yLocRes - locY) * (yLocRes - locY) + (zLocRes - locZ) * (zLocRes - locZ));
			//update velocity using updated and old positions and forces. 
			/*double xVelRes = abs(xLocRes - locX); /* velX + ((dt / 2 * mass)) * (prevAccX + accX) - 
				viscosity * (xLocRes - locX) + (1 / mass) * noise;
			//double yVelRes = abs(yLocRes - locY); /*velY + ((dt / 2 * mass)) * (prevAccY + accY) - 
				viscosity * (yLocRes - locY) + (1 / mass) * noise;
			//double zVelRes = abs(zLocRes - locZ)/dt; /*velZ + ((dt / 2 * mass)) * (prevAccZ + accZ) - 
				viscosity * (zLocRes - locZ) + (1 / mass) * noise;*/

			
			return thrust::make_tuple(xLocRes, yLocRes, zLocRes, velocity);//, xVelRes, yVelRes, zVelRes);
		}
		else {
			//do not move position or velocity. 
			return thrust::make_tuple(thrust::get<1>(p3), 
			thrust::get<2>(p3), thrust::get<3>(p3),0.0);/*
			thrust::get<4>(p3v3), thrust::get<5>(p3v3), 
			thrust::get<6>(p3v3));*/
		}                     
	}                                 

};

//old rk stuff
			/*
			//first calculation for position update
			double k2X = velX + dt / 2 * accX;
			double k2Y = velY + dt / 2 * accY;
			double k2Z = velZ + dt / 2 * accZ;
			
			double k3X = velX + dt / 2 * accX;
			double k3Y = velY + dt / 2 * accY;
			double k3Z = velZ + dt / 2 * accZ;
			
			double k4X = velX + dt * accX;
			double k4Y = velY + dt * accY;
			double k4Z = velZ + dt * accZ;
			
			//result new pos = oldpos + (rk4stuff);
			double xPosRes = thrust::get<1>(p3v3) + dt / 6 * (velX + 2 * k2X + 2 * k3X + k4X);
			double yPosRes = thrust::get<2>(p3v3) + dt / 6 * (velY + 2 * k2Y + 2 * k3Y + k4Y);
			double zPosRes = thrust::get<3>(p3v3) + dt / 6 * (velZ + 2 * k2Z + 2 * k3Z + k4Z);
			
			//second calculation for velocity update
			//result
			double xVelRes = velX + dt * accX;
			double yVelRes = velY + dt * accY;
			double zVelRes = velZ + dt * accZ;
			*/

//other numeric method
			/*double locX = thrust::get<1>(p3v3);
			double locY = thrust::get<2>(p3v3);
			double locZ = thrust::get<3>(p3v3);

			double velX = thrust::get<4>(p3v3);
			double velY = thrust::get<5>(p3v3);
			double velZ = thrust::get<6>(p3v3);
			
			//same as k1Vel (don't forget v' = F/m)
			//We subtract off a viscous term
			double accX = (thrust::get<4>(g1pf3f3)  - viscosity * velX) / mass;
			double accY = (thrust::get<5>(g1pf3f3)  - viscosity * velY) / mass;
			double accZ = (thrust::get<6>(g1pf3f3)  - viscosity * velZ) / mass;

			double prevAccX = thrust::get<1>(g1pf3f3);
			double prevAccY = thrust::get<2>(g1pf3f3);
			double prevAccZ = thrust::get<3>(g1pf3f3);

			//random data
			double gaussianData = thrust::get<0>(g1pf3f3);

			double b = 1/ (1 + (viscosity*dt / 2));
			double noise = sqrt(2*viscosity*kB*temperature) * gaussianData;

			//update positions
			double xLocRes = xLocRes = locX + b * dt * velX + (b*dt*dt/2)*accX + (b*dt/2) * noise;
			double yLocRes = yLocRes = locY + b * dt * velY + (b*dt*dt/2)*accY + (b*dt/2) * noise;
			double zLocRes = zLocRes = locZ + b * dt * velZ + (b*dt*dt/2)*accZ + (b*dt/2) * noise;

			//update velocity using updated and old positions and forces. 
			double xVelRes = velX + (dt/2)*(prevAccX + accX) - viscosity*(xLocRes - locX) + noise;
			double yVelRes = velY + (dt/2)*(prevAccY + accY) - viscosity*(yLocRes - locY) + noise;
			double zVelRes = velZ + (dt/2)*(prevAccZ + accZ) - viscosity*(zLocRes - locZ) + noise;*/

////////////////////////////////////////////////
/*double prevAccX = thrust::get<1>(g1pf3f3);
			double prevAccY = thrust::get<2>(g1pf3f3);
			double prevAccZ = thrust::get<3>(g1pf3f3);

			//random data
			double gaussianData = thrust::get<0>(g1pf3f3);

			double noise = sqrt(2 * dt * viscosity * kB * temperature / mass) * gaussianData;
			
			double accTestX = prevAccX/mass - viscosity * velX;
			double accTestY = prevAccY/mass - viscosity * velY;
			double accTestZ = prevAccZ/mass - viscosity * velZ;

			double vel_tildeX = velX + dt/2*(accTestX) + noise;
			double vel_tildeY = velY + dt/2*(accTestY) + noise;
			double vel_tildeZ = velZ + dt/2*(accTestZ) + noise;
			
			double xLocRes = xLocRes = locX + dt * vel_tildeX;
			double yLocRes = yLocRes = locY + dt * vel_tildeY;
			double zLocRes = zLocRes = locZ + dt * vel_tildeZ;



			//update positions
			double accX = (thrust::get<4>(g1pf3f3)  - viscosity * vel_tildeX) / mass;
			double accY = (thrust::get<5>(g1pf3f3)  - viscosity * vel_tildeY) / mass;
			double accZ = (thrust::get<6>(g1pf3f3)  - viscosity * vel_tildeZ) / mass;
		

			//update velocity using updated and old positions and forces. 
			double xVelRes = vel_tildeX + (dt/2)*(accX) + noise;
			double yVelRes = vel_tildeY + (dt/2)*(accY) + noise;
			double zVelRes = vel_tildeZ + (dt/2)*(accZ) + noise;
*/

#endif /*ADVANCEPOSITIONONDEVICE_H_*/
