/**
 * Handles the data structure and force spreading calculations for the 
 * bending force connections.
 *
 * Boundary Conditions: For PointId < 0, the point is assumed to be fictitious. 
 *
 * NOTE: When adding new force connection, you need to update ../../solvers/cythonlib/IBHelper.pxi
 * to handle conversion from Python Dictionary to C IB Data structure.
 *
 * Author: Jeffrey Wiens
 **/


#ifndef GUARD_BENDINGFORCE
#define GUARD_BENDINGFORCE

#include "IForceConnection.h"
#include "IBPoint.h"
#include "../utils/PetscUtil.h"
#include "../utils/Util.h"
#include "../io/IO.h"
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <petscdmda.h>
#include <stdexcept>
#include <sstream>
#include <unordered_map>
#include <map>
#include <vector>
#include <math.h>
#include <string>
#include <algorithm>
#include <functional>

namespace IB
{

   template <int dim> class ImmersedBoundary;

   namespace ForceConnection
   {
   
      template <int dim>
      struct bendingforce_data_template
      {
         int ForceID;
         int FCID;
         int L2PointID;
         int LPointID;
         int RPointID;
         int R2PointID;
         int PointID;
         double sigma;
         //double D2X_TARGET[dim];
         //double D2X_LTARGET[dim];   
         //double D2X_RTARGET[dim];
	double A;
	double K;
	double Omega;
         double dv;
         double ds;
      } ;      
       
      template <int dim>
      class BendingForce : public IForceConnection<dim>
      {
         private:
            typedef IB::IBPointStruct<dim> Point;
            typedef NetworkIBPointStruct<dim> NetworkPoint;
            typedef IB::ImmersedBoundary<dim> IBPoints;
            typedef bendingforce_data_template<dim> bendingforce_data;
            PetscErrorCode ierr;
            PetscMPIInt size,rank;
                
         public:
            // Constructor
            BendingForce(); 

            // Deconstructor
            ~BendingForce();
    
            // Distribute Immersed Boundary Points and Force Connections
            PetscErrorCode DistributeIBPoints(IBPoints *ib, std::unordered_map<int, Point > *pointDict, const int dir, std::unordered_map<int, int > *distribute_points, std::vector< int > *output_buffer, std::vector< int > *force_output_buffer);

            // Calculate Force Density
            PetscErrorCode CalculateForceDensity(std::unordered_map<int, Point > *pointDict, const double *domain_length,const double *spatial_stepsize, const double dt, const double time);

            // Clean Immersed Boundary Points and Force Connections Not in the Domain
            PetscErrorCode CleanDataNotInDomain(std::unordered_map<int, Point > *pointDict);
            static bool CleanupForceConnections(const bendingforce_data_template<dim> fc, void **params);
            
            // Add Force Connections
            PetscErrorCode AddLocalIBForceConnections(std::unordered_map<int, Point > *pointDict, const int *buffer, int *i, const int exclude_not_in_domain  );

            // Public data structure
            std::vector<bendingforce_data> force_connections;
            std::unordered_map<int, int > force_lookup;
      };
   }

}


#endif
