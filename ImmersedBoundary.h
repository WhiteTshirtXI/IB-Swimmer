/**
 * Evolves immersed boundary, handles force spreading, and distributes data across processes.
 * Author: Jeffrey Wiens
 * Updated by Saeed June 30, 2016
 **/


#ifndef GUARD_IMMERSEDBOUNARY
#define GUARD_IMMERSEDBOUNARY

#include "IBPoint.h"
#include "IForceConnection.h"
#include "ElasticForce.h"
#include "ElasticForce_2Pt.h"
#include "BendingForce.h"
#include "ElasticForce_OscStiffness.h"
#include "ElasticForce_OscStiffness_2Pt.h"
#include "ElasticForce_OscLength_2Pt.h"
#include "ElasticForce_OscLength.h"
#include "ElasticForce_Tether_2Pt.h"
#include "ElasticForce_Tether_OscStiffness_2Pt.h"
#include "Force_Pulse_2Pt.h"
#include "PenaltyMass.h"
//Added by Saeed
#include "BendingForceCurve.h"
#include "ElasticForceEnergy.h"

#include "../utils/PetscUtil.h"
#include "../utils/Util.h"
#include "../io/IO.h"
#include "../io/RestartSim.h"
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

namespace IB
{
  // Identifies different Force Connections
  const int FORCE_CONNECTION_ELASTIC = 1;
  const int FORCE_CONNECTION_ELASTIC_TWOPOINT = 2;
  const int FORCE_CONNECTION_ELASTIC_OSCILLATINGSTIFFNESS = 3;
  const int FORCE_CONNECTION_ELASTIC_OSCILLATINGSTIFFNESS_TWOPOINT = 4;
  const int FORCE_CONNECTION_ELASTIC_OSCILLATINGLENGTH = 5;
  const int FORCE_CONNECTION_ELASTIC_OSCILLATINGLENGTH_TWOPOINT = 6;
  const int FORCE_CONNECTION_ELASTIC_TETHER_TWOPOINT = 7;
  const int FORCE_CONNECTION_ELASTIC_TETHER_OSCILLATINGSTIFFNESS_TWOPOINT = 8;
  const int FORCE_CONNECTION_BENDING = 9;
  const int FORCE_CONNECTION_PENALTYMASS = 10;
  const int FORCE_CONNECTION_FORCEPULSE_TWOPOINT = 11;
  //Added by Saeed
  const int FORCE_CONNECTION_BENDING_CURVE = 12;
  const int FORCE_CONNECTION_ELASTIC_ENERGY = 13;

  // Identifies different Force Connections
  const int IBPOINT_AND_FORCECONNECTION_DIVIDER = -1;  

  const int IBDISTRIBUTE_NOTHING = 0;
  const int IBDISTRIBUTE_SENDIB = 1;
  const int IBDISTRIBUTE_SENDFORCEANDIB = 2;

  double inline phi( const double r);

  template <int dim>
  class ImmersedBoundary
  {
     protected:
       typedef IBPointStruct<dim> Point;
       typedef NetworkIBPointStruct<dim> NetworkPoint;
       typedef std::unordered_map< int, IB::ForceConnection::IForceConnection<dim>* > FC_Map;
       typedef typename std::unordered_map< int, Point>::iterator IB_Iterator;
       typedef typename FC_Map::iterator FC_Iterator;
       typedef PetscErrorCode (*CallUpdateIBData)(double *X, double *Xh, double *Uold, double *Unew, const double dt, const int initialize);

     public:
       double dt;

       // Constructor
       ImmersedBoundary(); 

       // Deconstructor
       ~ImmersedBoundary();
    
       // Initialize
       PetscErrorCode Initialize( std::vector<int> point_buffer, std::vector<int> force_buffer, DM *da_vector, double *h, double *domain_length, PetscInt _delta_span, double _dt );

       // Update Immersed Boundary Points Position
       PetscErrorCode UpdateIBPosition(const DM *da_vector, const Vec *Umac, const int initialize);
       PetscErrorCode UpdateIBPosition(const DM *da_vector, const Vec *Umac, const int initialize, const int interpolate_nextstep, CallUpdateIBData update);

       // Distribute Immersed Boundary Points
       PetscErrorCode DistributeIBPoints(DM *da_vector);

       // Calculate Force Density and Spread Force
       PetscErrorCode SpreadForce(DM *da_vector, Vec *F, const double time);
       PetscErrorCode SpreadForce(DM *da_vector, Vec *F, const double time, const int zero_force, const double magnify_force);

       // Clean Immersed Boundary Points and Force Connections Not in the Domain
       PetscErrorCode CleanDataNotInDomain(unsigned int keep_nextstep_points);

       // Save Immersed Boundary
       PetscErrorCode OutputIBVar(IO::OutputHDF5 *logger, DM *da);

       // Project IB onto Cartesian grid for output
       PetscErrorCode ConstructOmega(DM *da, Vec *Omega);

       // Determine if IB node should be sent to particular neighbour
       int DisributeToNeighbour(const int dir, const double *X);

     protected:
       PetscBool reduced_fluid_output;

       // Add Immersed Boundary Points
       PetscErrorCode AddLocalIBPoints( const int *buffer, const int buffer_length, int *last_index, const int exclude_not_in_domain );

       // Add Force Connections
       PetscErrorCode AddLocalIBForceConnections( const int *buffer, const int buffer_length, const int exclude_not_in_domain );

       // Set inside/outside domain flags
       PetscErrorCode SetIBPointDomainFlags(Point *ptr_point);

       PetscMPIInt size,rank;

       FC_Map fcDict;
       std::unordered_map<int, Point > pointDict;
       
       std::vector< int > output_buffer;
       std::vector< std::vector<int> > input_buffer;
       std::vector< int > force_output_buffer;

       double xstart[3]; double xend[3];
       double xghoststart[3]; double xghostend[3];
       double extended_xinsideghoststart[3]; double extended_xinsideghostend[3];
       double xinsideghoststart[3]; double xinsideghostend[3];
       double domain_length[3];
       double spatial_stepsize[3];
       PetscInt spatial_gridsize[3];
       unsigned int delta_span;
       PetscErrorCode ierr;

       PetscLogEvent IB_INIT, IB_UPDATEIBPOSITION, IB_FORCEDENSITY, IB_SPREADFORCE, IB_DISTRIBUTE, IB_OUTPUT, IB_CLEANUP;
  };

}


#endif
