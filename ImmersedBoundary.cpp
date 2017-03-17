/**
 * Evolves immersed boundary, handles force spreading, and distributes data across processes.
 * Author: Jeffrey Wiens
 *Updated by Saeed June30, 2016
 **/

#include "ImmersedBoundary.h"

namespace IB {

/**
 * Constructor and Deconstructor for Immersed Boundary
 */
template <int dim>
ImmersedBoundary<dim>::ImmersedBoundary()
{ 
   output_buffer.reserve( 400*sizeof(NetworkPoint)/sizeof(int) );
   force_output_buffer.reserve( 200*sizeof(IB::ForceConnection::elasticforce_data)/sizeof(int) );

   int max_dir = (dim==2 ? 9 : 27);
   input_buffer.resize(max_dir);

   for(int i = 0; i < max_dir; i++ )
       input_buffer[i].reserve( 400*sizeof(NetworkPoint)/sizeof(int) );

   // Setup log
   PetscLogEventRegister("IB:Init", 0 , &IB_INIT);
   PetscLogEventRegister("IB:EvolveIB", 0 , &IB_UPDATEIBPOSITION);
   PetscLogEventRegister("IB:ForceDensity", 0 , &IB_FORCEDENSITY);
   PetscLogEventRegister("IB:SpreadForce", 0 , &IB_SPREADFORCE);
   PetscLogEventRegister("IB:DistributeIB", 0 , &IB_DISTRIBUTE);
   PetscLogEventRegister("IB:Cleanup", 0 , &IB_CLEANUP);
   PetscLogEventRegister("IB:OutputData", 0 , &IB_OUTPUT);
   
   // Add Force Connections
   fcDict[FORCE_CONNECTION_ELASTIC] = new IB::ForceConnection::ElasticForce<dim>();
   fcDict[FORCE_CONNECTION_ELASTIC_TWOPOINT] = new IB::ForceConnection::ElasticForceTwoPoint<dim>(); 
   fcDict[FORCE_CONNECTION_BENDING] = new IB::ForceConnection::BendingForce<dim>();
   fcDict[FORCE_CONNECTION_ELASTIC_OSCILLATINGSTIFFNESS] = new IB::ForceConnection::ElasticForceOscillatingStiffness<dim>();
   fcDict[FORCE_CONNECTION_ELASTIC_OSCILLATINGLENGTH] = new IB::ForceConnection::ElasticForceOscillatingLength<dim>();
   fcDict[FORCE_CONNECTION_ELASTIC_OSCILLATINGSTIFFNESS_TWOPOINT] = new IB::ForceConnection::ElasticForceOscillatingStiffnessTwoPoint<dim>();
   fcDict[FORCE_CONNECTION_ELASTIC_OSCILLATINGLENGTH_TWOPOINT] = new IB::ForceConnection::ElasticForceOscillatingLengthTwoPoint<dim>();   
   fcDict[FORCE_CONNECTION_ELASTIC_TETHER_TWOPOINT] = new IB::ForceConnection::ElasticForceTetherTwoPoint<dim>(); 
   fcDict[FORCE_CONNECTION_ELASTIC_TETHER_OSCILLATINGSTIFFNESS_TWOPOINT] = new IB::ForceConnection::ElasticForceTetherOscillatingStiffnessTwoPoint<dim>(); 
   fcDict[FORCE_CONNECTION_PENALTYMASS] = new IB::ForceConnection::PenaltyMass<dim>();
   fcDict[FORCE_CONNECTION_FORCEPULSE_TWOPOINT] = new IB::ForceConnection::ForcePulseTwoPoint<dim>();
   //added by Saeed
   fcDict[FORCE_CONNECTION_BENDING_CURVE] = new IB::ForceConnection::BendingForceCurve<dim>();
   fcDict[FORCE_CONNECTION_ELASTIC_ENERGY] = new IB::ForceConnection::ElasticForceEnergy<dim>();
}

template <int dim>
ImmersedBoundary<dim>::~ImmersedBoundary()
{
   // delete Force Connections
   delete fcDict[FORCE_CONNECTION_ELASTIC] ;
   delete fcDict[FORCE_CONNECTION_ELASTIC_TWOPOINT] ;
   delete fcDict[FORCE_CONNECTION_BENDING] ;
   delete fcDict[FORCE_CONNECTION_ELASTIC_OSCILLATINGSTIFFNESS] ;
   delete fcDict[FORCE_CONNECTION_ELASTIC_OSCILLATINGLENGTH] ;
   delete fcDict[FORCE_CONNECTION_ELASTIC_OSCILLATINGSTIFFNESS_TWOPOINT] ;
   delete fcDict[FORCE_CONNECTION_ELASTIC_OSCILLATINGLENGTH_TWOPOINT] ;
   delete fcDict[FORCE_CONNECTION_ELASTIC_TETHER_TWOPOINT] ;
   delete fcDict[FORCE_CONNECTION_ELASTIC_TETHER_OSCILLATINGSTIFFNESS_TWOPOINT] ;
   delete fcDict[FORCE_CONNECTION_PENALTYMASS] ;
   delete fcDict[FORCE_CONNECTION_FORCEPULSE_TWOPOINT] ;
   delete fcDict[FORCE_CONNECTION_BENDING_CURVE] ;
   delete fcDict[FORCE_CONNECTION_ELASTIC_ENERGY] ;
}

/**
 * Initialize Immersed Boundary
 * 
 * Parameters:
 *        point_buffer: Buffer filled with all the immersed boundary points. The buffer is of
 *                      type int, so they need to be cast to correct struct type. Also, points
 *                      not in the domain need to be thrown out.
 *        force_buffer: Buffer filled with all the force connections. The buffer is of
 *                      type int, so they need to be cast to correct struct type. Also, points
 *                      not in the domain need to be thrown out.
 *        da_vector:    The DA used for fluid variables, used to obtain boundary of cartesian 
 *                      grid owned by process.
 *        h:            Spatial step size in each dimension.
 *        domain_length: Domain length
 *        delta_span:    Span of the discrete delta function.
 *        dt:            The time-step
 */
template <int dim>
PetscErrorCode ImmersedBoundary<dim>::Initialize( std::vector<int> point_buffer, std::vector<int> force_buffer, DM *da_vector, double *h, double *_domain_length, PetscInt _delta_span, double _dt )
{
   PetscLogEventBegin(IB_INIT, 0, 0, 0, 0);

   // Get Rank and Size
   ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);
   ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &rank); CHKERRQ(ierr);

   // Time-step, Span of discrete delta function
   delta_span = _delta_span;
   dt = _dt;

   // Retrieve domain length and step size
   for(int i = 0; i < dim; i++) domain_length[i] = _domain_length[i];
   for(int i = 0; i < dim; i++) spatial_stepsize[i] = h[i];

   // Get number of grid points in each spatial dimension
   ierr = DMDAGetInfo(*da_vector, PETSC_NULL, &spatial_gridsize[0], &spatial_gridsize[1], &spatial_gridsize[2], PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL); CHKERRQ(ierr);

   // Retrieve the start and end of the domain
   PetscInt is, im, js, jm, ks, km;
   ierr = DMDAGetCorners(*da_vector, &is, &js, &ks, &im, &jm, &km); CHKERRQ(ierr);

   // Domain owned by Processor (xstart <= x < xend)
   xstart[0] = is*h[0]; 
   xstart[1] = js*h[1];
   if( dim == 3) xstart[2] = ks*h[2];

   xend[0] = (is+im)*h[0]; 
   xend[1] = (js+jm)*h[1];
   if( dim == 3) xend[2] = (ks+km)*h[2];

   // Ghost cell region outside the processor's domain (xghoststart < xstart <= x < xend < xghostend)
   // IB points AND Force Connections in this region get copied to adjacent processors .
   // (xghoststart <= x <= xinsideghoststart OR xinsideghostend <= x <= xghostend)
   ierr = DMDAGetGhostCorners(*da_vector, &is, &js, &ks, &im, &jm, &km); CHKERRQ(ierr);

   xghoststart[0] = is*h[0]; 
   xghoststart[1] = js*h[1];
   if( dim == 3) xghoststart[2] = ks*h[2];

   xghostend[0] = (is+im)*h[0]; 
   xghostend[1] = (js+jm)*h[1];
   if( dim == 3) xghostend[2] = (ks+km)*h[2];

   // Ghost cell region inside the processor's domain (xstart < xinsideghoststart <= x < xinsideghostend < xend)
   // IB points AND Force Connections in this region get copied to adjacent processors.
   // (xghoststart <= x <= xinsideghoststart OR xinsideghostend <= x <= xghostend)
   xinsideghoststart[0] = xghostend[0] - xend[0] + xstart[0];
   xinsideghoststart[1] = xghostend[1] - xend[1] + xstart[1];
   if( dim == 3 ) xinsideghoststart[2] = xghostend[2] - xend[2] + xstart[2];

   xinsideghostend[0] = xend[0] - (xstart[0] - xghoststart[0]);
   xinsideghostend[1] = xend[1] - (xstart[1] - xghoststart[1]);
   if( dim == 3 ) xinsideghostend[2] = xend[2] - (xstart[2] - xghoststart[2]);

   // Extended Ghost cell region inside the processor's domain 
   // (xstart < xinsideghoststart < extended_xinsideghoststart <= x < extended_xinsideghostend < xinsideghostend < xend)
   // IB points (NO FORCE CONNECTIONS) in this region get copied to adjacent processors.
   // This ensures that IB points aren't missing when calculating the force density.
   extended_xinsideghoststart[0] = xstart[0] + 1.5*(xghostend[0] - xend[0]);
   extended_xinsideghoststart[1] = xstart[1] + 1.5*(xghostend[1] - xend[1]);
   if( dim == 3) extended_xinsideghoststart[2] = xstart[2] + 1.5*(xghostend[2] - xend[2]);

   extended_xinsideghostend[0] = xend[0] - 1.5*(xstart[0] - xghoststart[0]);
   extended_xinsideghostend[1] = xend[1] - 1.5*(xstart[1] - xghoststart[1]);
   if( dim == 3) extended_xinsideghostend[2] = xend[2] - 1.5*(xstart[2] - xghoststart[2]);

   // Replace buffer with data from checkpoint, if we are restarting
   ierr = IO::ReplaceIBPointBuffer<dim>(&point_buffer); CHKERRQ(ierr);

   // Add Force Connections and IB Points from buffer into data structure
   ierr = AddLocalIBPoints( &point_buffer[0], point_buffer.size(), PETSC_NULL, 1 ); CHKERRQ(ierr);

   ierr = AddLocalIBForceConnections( &force_buffer[0], force_buffer.size(), 1 ); CHKERRQ(ierr);

   // Clear any of the IB Points and Force Connections that don't reside in the processors domain
   ierr = CleanDataNotInDomain(0); CHKERRQ(ierr);

   // Check - Output less fluid data to reduce size
   ierr = PetscOptionsHasName(PETSC_NULL, "-ReducedFluidOutput", &reduced_fluid_output);

   PetscLogEventEnd(IB_INIT, 0, 0, 0, 0);

   return 0;
}

/**
 * Update the position of all immersed boundary points.
 * 
 * Parameters:
 *        da_vector:    The distributed array that contains info on fluid velocity field
 *        Umac:         The fluid velocity field
 *        initialize:   Is this the first time-step? (initialize > 0 => yes)
 *        interpolate_nextstep: Interpolate position of IB point at the next step (1 = interpolate Xh, 0 = interpolate X)
 *        update:       A function which updates the IB data structure (allows for customizations)
 */
template <int dim>
PetscErrorCode ImmersedBoundary<dim>::UpdateIBPosition(const DM *da_vector, const Vec *Umac, const int initialize, const int interpolate_nextstep, CallUpdateIBData update)
{
   PetscScalar ****array_Umac;
   PetscScalar ***array_Umac_2d;

   double Unew[dim];
   double X[dim];
   int index_start[3];
   int index_end[3];
   double delta_val[6];
   double inv_dx[dim];
   double delta_param[dim];
   
   PetscLogEventBegin(IB_UPDATEIBPOSITION, 0, 0, 0, 0);

   index_start[2] = 0;
   index_end[2] = 0;

   for(unsigned int i = 0; i < dim;i++) {
      inv_dx[i] = 1.0/spatial_stepsize[i];
   }

   if( dim == 2 ) {
      ierr = DMDAVecGetArrayDOF(*da_vector, *Umac, &array_Umac_2d); CHKERRQ(ierr);
      array_Umac = &array_Umac_2d;
   } else {
      ierr = DMDAVecGetArrayDOF(*da_vector, *Umac, &array_Umac); CHKERRQ(ierr);
   }

   for (IB_Iterator it = pointDict.begin(); it != pointDict.end();it++)
   {
       // initialize data for immersed boundary point
       Point *pt = &(it->second);
       for(unsigned int i = 0; i < dim;i++) Unew[i] = 0;
       if(interpolate_nextstep) { for(unsigned int i = 0; i < dim; i++){ X[i] = pt->Xh[i]; } }
       else { for(unsigned int i = 0; i < dim; i++){ X[i] = pt->X[i]; } }

       for( unsigned int i = 0; i < dim; i++ )
       {
       	  // Compute portion of parameter used in discrete delta function
          delta_param[i] = (X[i])*inv_dx[i];
          
          // Interpolate the fluid velocity onto the immersed boundary point
          index_start[i] = floor( inv_dx[i]*X[i] - .5*delta_span );
          index_end[i] = index_start[i] + delta_span;
       }


       for(int k = index_start[2]; k <= index_end[2]; k++)
       {
           delta_val[2] = (dim == 2 ? 1.0 : phi( k - delta_param[2] ));
           delta_val[5] = (dim == 2 ? 1.0 : phi( k+.5 - delta_param[2] ));

           for(int j = index_start[1]; j <= index_end[1]; j++)
           {
               delta_val[1] = phi( j - delta_param[1] );
               delta_val[4] = phi( j+.5 - delta_param[1] );

               for(int i = index_start[0]; i <= index_end[0]; i++)
               {
                   delta_val[0] = phi( i - delta_param[0] );
                   delta_val[3] = phi( i+.5 - delta_param[0] );
                   
                   Unew[0] += array_Umac[k][j][i][0]*delta_val[0]*delta_val[4]*delta_val[5];
                   Unew[1] += array_Umac[k][j][i][1]*delta_val[1]*delta_val[3]*delta_val[5];
                   if(dim==3) Unew[2] += array_Umac[k][j][i][2]*delta_val[2]*delta_val[3]*delta_val[4];
               }
           }
       }

       // Evolve IB forward in time
       ierr = update(pt->X, pt->Xh, pt->U, Unew, dt, initialize); CHKERRQ(ierr);

       for(unsigned int i = 0; i < dim;i++)
       {
          // Set back to zero if too close
          pt->X[i] = ( fabs(pt->X[i]) < 1e-10 ? 0.0 : pt->X[i] );
          pt->Xh[i] = ( fabs(pt->Xh[i]) < 1e-10 ? 0.0 : pt->Xh[i] );

          // Zero Force Density
          pt->ForceDensity[i] = 0.0;
       }

       ierr = SetIBPointDomainFlags(pt); CHKERRQ(ierr);
   }

   if( dim == 2 ){
      ierr = DMDAVecRestoreArrayDOF(*da_vector, *Umac, &array_Umac_2d); CHKERRQ(ierr);
   } else {
      ierr = DMDAVecRestoreArrayDOF(*da_vector, *Umac, &array_Umac); CHKERRQ(ierr);
   }

   PetscLogEventEnd(IB_UPDATEIBPOSITION, 0, 0, 0, 0);

   return 0;
}

template <int dim>
PetscErrorCode ImmersedBoundary<dim>::UpdateIBPosition(const DM *da_vector, const Vec *Umac, const int initialize)
{
   SETERRQ(PETSC_COMM_WORLD,1,"UpdateIBPosition function is not implemented.");
   return 1;
}


/**
 * Calculate the force density and Spread onto Cartesian Grid
 * 
 * Parameters:
 *        da_vector:    The distributed array that contains info on external force field
 *        F:            The external force field
 *        time:		The time at the current time-step
 *        zero_force:	When 1, the force vector F will be set to zero before spreading force density
 *        magnify_force: Multiple the force density by this value.
 */
template <int dim>
PetscErrorCode ImmersedBoundary<dim>::SpreadForce(DM *da_vector, Vec *F, const double time, const int zero_force, const double magnify_force)
{
   PetscScalar ****array_Fmac;
   PetscScalar ***array_Fmac_2d;
   PetscInt is, im, js, jm, ks, km, i, j, k;
   double iDist, jDist, kDist;

   int index_start[3];
   double delta_val[6];
   double inv_dx[dim];
   double inv_dv = 1.0;
   
   for(int i = 0; i < dim; i++)
   {
      inv_dx[i] = 1.0/spatial_stepsize[i];
      inv_dv = inv_dv*inv_dx[i];
   }
   
   // Calculate Force Density
   PetscLogEventBegin(IB_FORCEDENSITY, 0, 0, 0, 0);   
   for (FC_Iterator it = fcDict.begin(); it != fcDict.end();it++)
   {
       ierr = (it->second)->CalculateForceDensity(&pointDict, &domain_length[0], &spatial_stepsize[0], dt, time); CHKERRQ(ierr);
   }
   PetscLogEventEnd(IB_FORCEDENSITY, 0, 0, 0, 0);


   PetscLogEventBegin(IB_SPREADFORCE, 0, 0, 0, 0);

   // Set External Force to zero
   if(zero_force) { ierr = VecSet(*F, 0.0); CHKERRQ(ierr);}

   ierr = DMDAGetCorners(*da_vector, &is, &js, &ks, &im, &jm, &km); CHKERRQ(ierr);

   if( dim == 2 ) {
      ierr = DMDAVecGetArrayDOF(*da_vector, *F, &array_Fmac_2d); CHKERRQ(ierr);
      array_Fmac = &array_Fmac_2d;
   } else {
      ierr = DMDAVecGetArrayDOF(*da_vector, *F, &array_Fmac); CHKERRQ(ierr);
   }

   for (IB_Iterator it = pointDict.begin(); it != pointDict.end();it++)
   {
       // Get Point
       Point *pt =  &(it->second);

       // Find Span of delta function for IB point
       for( int d = 0; d < dim; d++) {
          index_start[d] = Mod( (int)floor(pt->Xh[d]*inv_dx[d] -  .5*delta_span), spatial_gridsize[d] );
       }

       // Spread Force
       for(int kk = 0; kk <= (dim==2 ? 1 : delta_span); kk++)
       {
           if( dim == 2){k = kk; delta_val[2] = 1.0; delta_val[5] = 1.0;}
           else {
              k = (index_start[2] + kk)%spatial_gridsize[2]; 
              kDist = ibdistance( k*spatial_stepsize[2], pt->Xh[2], domain_length[2]);
              delta_val[2] = phi( kDist*inv_dx[2] );
              delta_val[5] = phi( (.5*spatial_stepsize[2] + kDist)*inv_dx[2] );
           }


           for(int jj = 0; jj <= delta_span; jj++)
           {
              j = (index_start[1] + jj)%spatial_gridsize[1]; 
              jDist = ibdistance( j*spatial_stepsize[1], pt->Xh[1], domain_length[1]);
              delta_val[1] = phi( jDist*inv_dx[1] );
              delta_val[4] = phi( (.5*spatial_stepsize[1] + jDist)*inv_dx[1] ); 
              
              for(int ii = 0; ii <= delta_span; ii++)
              {
                   i = (index_start[0] + ii)%spatial_gridsize[0]; 
                   
                   if( (k >= ks && k < ks+km) && (j >= js && j < js+jm) && (i >= is && i < is+im)) 
                   {
                      iDist = ibdistance( i*spatial_stepsize[0], pt->Xh[0], domain_length[0]);
                      delta_val[0] = phi( iDist*inv_dx[0] );
                      delta_val[3] = phi( (.5*spatial_stepsize[0] + iDist)*inv_dx[0] );     
                        
                      array_Fmac[k][j][i][0] += pt->ForceDensity[0]*delta_val[0]*delta_val[4]*delta_val[5]*inv_dv*magnify_force;
                      array_Fmac[k][j][i][1] += pt->ForceDensity[1]*delta_val[1]*delta_val[3]*delta_val[5]*inv_dv*magnify_force;
                      if( dim==3 ) array_Fmac[k][j][i][2] += pt->ForceDensity[2]*delta_val[2]*delta_val[3]*delta_val[4]*inv_dv*magnify_force;
                   }
              }
           }
       }
   }

   if( dim == 2 ){
      ierr = DMDAVecRestoreArrayDOF(*da_vector, *F, &array_Fmac_2d); CHKERRQ(ierr);
   } else {
      ierr = DMDAVecRestoreArrayDOF(*da_vector, *F, &array_Fmac); CHKERRQ(ierr);
   }

   PetscLogEventEnd(IB_SPREADFORCE, 0, 0, 0, 0);

   return 0;
}

template <int dim>
PetscErrorCode ImmersedBoundary<dim>::SpreadForce(DM *da_vector, Vec *F, const double time)
{
   return SpreadForce(da_vector, F, time, 1, 1.0);
}


/**
 * Add immersed boundary points to data structure. Update flags indicating where on the domain 
 * the point lies.
 * 
 * Parameters:
 *        buffer: Buffer filled with all the immersed boundary points. The buffer is of
 *                      type int, so they need to be cast to correct struct type. 
 *        buffer_length: Length of the buffer.
 *        last_index: The index of the last data point that was read. The buffer may contain force connections.
 *                    The function will stop when the force connections are reached.
 *        exclude_not_in_domain: Flag determining if the point should be excluded if it
 *                                isn't in the proc's domain.
 */
template <int dim>
PetscErrorCode ImmersedBoundary<dim>::AddLocalIBPoints( const int *buffer, const int buffer_length, int *last_index, const int exclude_not_in_domain )
{
   int i = 0;
   while(i < buffer_length && buffer[i] != IBPOINT_AND_FORCECONNECTION_DIVIDER)
   {
      // Retrieve Point
      NetworkPoint network_pt = *((NetworkPoint *) &(buffer[i]) );
     
      // Load Point
      Point pt;
      pt.PointID = network_pt.PointID;
      for(int index = 0; index < dim; index++)
      {
          pt.X[index] = network_pt.X[index];
          pt.Xh[index] = network_pt.Xh[index];
          pt.U[index] = network_pt.U[index];
          pt.ForceDensity[index] = 0.0;
      }

      // Set to zero if too close
      pt.X[0] = ( fabs(pt.X[0]) < 1e-10 ? 0.0 : pt.X[0] );
      pt.X[1] = ( fabs(pt.X[1]) < 1e-10 ? 0.0 : pt.X[1] );
      if(dim == 3) pt.X[2] = ( fabs(pt.X[2]) < 1e-10 ? 0.0 : pt.X[2] );

      pt.Xh[0] = ( fabs(pt.Xh[0]) < 1e-10 ? 0.0 : pt.Xh[0] );
      pt.Xh[1] = ( fabs(pt.Xh[1]) < 1e-10 ? 0.0 : pt.Xh[1] );
      if(dim == 3) pt.Xh[2] = ( fabs(pt.Xh[2]) < 1e-10 ? 0.0 : pt.Xh[2] );

      // Enforce that IB point resides in the domain
      pt.X[0] = ibposition( pt.X[0], (xstart[0]+xend[0])/2, domain_length[0]); 
      pt.X[1] = ibposition( pt.X[1], (xstart[1]+xend[1])/2, domain_length[1]); 
      if(dim == 3) pt.X[2] =ibposition( pt.X[2], (xstart[2]+xend[2])/2, domain_length[2]); 

      pt.Xh[0] = ibposition( pt.Xh[0], (xstart[0]+xend[0])/2, domain_length[0]); 
      pt.Xh[1] = ibposition( pt.Xh[1], (xstart[1]+xend[1])/2, domain_length[1]); 
      if(dim == 3) pt.Xh[2] =ibposition( pt.Xh[2], (xstart[2]+xend[2])/2, domain_length[2]); 

      // Update IB flags
      ierr = SetIBPointDomainFlags(&pt); CHKERRQ(ierr);

      // Add to data structure
      if( exclude_not_in_domain == 0 || pt.I == IBPOINT_INDOMIAN ) pointDict[ pt.PointID ] = pt;

      i += sizeof(NetworkPoint)/sizeof(int);
   }

   if( last_index != PETSC_NULL )
     *last_index = ( i < buffer_length ? i+1 : i );

   return 0;
}

/**
 * Add force connection to appriopriate data structures.
 * 
 * Parameters:
 *        buffer: Buffer filled with all the force connections. The buffer is of
 *                      type int, so they need to be cast to correct struct type. 
 *        buffer_length: Length of the buffer.
 *        exclude_not_in_domain: Flag determining if the point should be excluded if it
 *                                isn't in the proc's domain.
 */
template <int dim>
PetscErrorCode ImmersedBoundary<dim>::AddLocalIBForceConnections( const int *buffer, const int buffer_length, const int exclude_not_in_domain  )
{
   int i = 0;

   while(i < buffer_length )
   {
      if( fcDict.count( buffer[i] ) > 0 ) {
          ierr = fcDict[buffer[i]]->AddLocalIBForceConnections(&pointDict, buffer, &i, exclude_not_in_domain  ); CHKERRQ(ierr);
      } else {
          SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_CORRUPT,"Invalid Force connection type. %D \n",buffer[i]);
          return 1;
      }
   }

   return 0;
}

/**
 * Sets flags indicating if the point resides in our outside the processors domain
 * 
 * Parameters:
 *        pt: Immersed boundary point
 */
template <int dim>
PetscErrorCode ImmersedBoundary<dim>::SetIBPointDomainFlags(Point *pt)
{
   // Update Points Position to account for periodicity of the domain
   // Note, these "updated" points are only used in a few locations
   double X[3];
   double Xh[3];
   for(unsigned int i=0; i<dim;i++) 
   {
       X[i] = ibposition( pt->X[i], .5*(xstart[i]+xend[i]), domain_length[i]); 
       Xh[i] = ibposition( pt->Xh[i], .5*(xstart[i]+xend[i]), domain_length[i]); 
   }

   pt->I = IBPOINT_OUTSIDEDOMIAN;
   if(xstart[0] <= pt->X[0] && xend[0] > pt->X[0] && xstart[1] <= pt->X[1] && xend[1] > pt->X[1] && 
        (dim==2 || (xstart[2] <= pt->X[2] && xend[2] > pt->X[2])) ) 
   {
      pt->I = IBPOINT_INDOMIAN;
   }
   else if(xghoststart[0] <= X[0] && xghostend[0] > X[0] && xghoststart[1] <= X[1] && xghostend[1] > X[1] &&
             (dim==2 || (xghoststart[2] <= X[2] && xghostend[2] > X[2])) ) 
   {
      pt->I = IBPOINT_INGHOSTDOMIAN;
   }

   pt->Ih = IBPOINT_OUTSIDEDOMIAN;
   if(xstart[0] <= pt->Xh[0] && xend[0] > pt->Xh[0] && xstart[1] <= pt->Xh[1] && xend[1] > pt->Xh[1] &&
         (dim==2 || (xstart[2] <= pt->Xh[2] && xend[2] > pt->Xh[2])) )
   {
      pt->Ih = IBPOINT_INDOMIAN;
   }
   else if(xghoststart[0] <= Xh[0] && xghostend[0] > Xh[0] && xghoststart[1] <= Xh[1] && xghostend[1] > Xh[1] &&
              (dim==2 || (xghoststart[2] <= Xh[2] && xghostend[2] > Xh[2])) ) 
   {
      pt->Ih = IBPOINT_INGHOSTDOMIAN;
   }

   return 0;
}

/**
 * Distribute IB points and force density to neighbouring processes
 */
template <int dim>
PetscErrorCode ImmersedBoundary<dim>::DistributeIBPoints(DM *da_vector)
{
    int max_dir = (dim==2 ? 9 : 27);
    std::unordered_map<int, int > distribute_points;
    int input_size[max_dir];
    int senders_buffer_size;
    int recvers_buffer_size;

    PetscLogEventBegin(IB_DISTRIBUTE, 0, 0, 0, 0);

    for(int dir = 0; dir < max_dir; dir++)
    {
        distribute_points.clear();
        output_buffer.clear();
        force_output_buffer.clear();

        for (typename FC_Map::iterator it = fcDict.begin(); it != fcDict.end();it++)
        {
              ierr = (it->second)->DistributeIBPoints(this, &pointDict, dir, &distribute_points, &output_buffer, &force_output_buffer);
        }

        // Merge Force buffer into point buffer (include separator)
        output_buffer.reserve( output_buffer.size() + force_output_buffer.size() + 1 );
        if(output_buffer.size() > 0) output_buffer.push_back( IBPOINT_AND_FORCECONNECTION_DIVIDER );
        output_buffer.insert( output_buffer.end(), force_output_buffer.begin(), force_output_buffer.end() );

        // Send and Receive Data in Buffer
        senders_buffer_size = output_buffer.size();

        ierr = PetscUtil::SendToNeighbour(dir, *da_vector, &senders_buffer_size, 1, &recvers_buffer_size, 2, &input_size[dir]); CHKERRQ(ierr);
        if( recvers_buffer_size >=  input_buffer[dir].capacity()) input_buffer[dir].reserve( recvers_buffer_size+1 );
        
        ierr = PetscUtil::SendToNeighbour(dir, *da_vector, &output_buffer[0], senders_buffer_size, &input_buffer[dir][0], input_buffer[dir].capacity(), &input_size[dir]); CHKERRQ(ierr);
    }

    // Move input buffers into IB data structures
    for(int dir = 0; dir < max_dir; dir++)
    {
        if( input_size[dir] > 0 ) 
        {
            int last_index;

            // Add Immersed Boundary Points
            ierr = AddLocalIBPoints( &input_buffer[dir][0], input_size[dir], &last_index, 0 ); CHKERRQ(ierr);

            // Add Force Connections
            if( last_index < input_size[dir] )
            {
                ierr = AddLocalIBForceConnections( &input_buffer[dir][last_index], input_size[dir]-last_index, 0 ); CHKERRQ(ierr);
            }
        }
    }

   PetscLogEventEnd(IB_DISTRIBUTE, 0, 0, 0, 0);

    return 0;
}


/**
 * Clean points and force connections not in the domain
 * 
 * Parameters:
 *        keep_nextstep_points: Keep the IB points pt->Xh that are in the domain (keep_nextstep_points=1 uses pt->Ih, keep_nextstep_points=0 uses pt->I )
 */
template <int dim>
PetscErrorCode ImmersedBoundary<dim>::CleanDataNotInDomain(unsigned int keep_nextstep_points)
{
   PetscLogEventBegin(IB_CLEANUP, 0, 0, 0, 0);

   // Erase Points not in the domain
   for (typename std::unordered_map<int,Point>::iterator it = pointDict.begin(); it != pointDict.end();)
   {
        if( !keep_nextstep_points && (it->second).I != IBPOINT_INDOMIAN ) it = pointDict.erase(it);
        else if( keep_nextstep_points && (it->second).Ih != IBPOINT_INDOMIAN ) it = pointDict.erase(it);
        else it++;
   }

   // Erase Force Connections not in the domain
   for (typename FC_Map::iterator it = fcDict.begin(); it != fcDict.end();it++)
   {
       ierr = (it->second)->CleanDataNotInDomain(&pointDict);
   }

   PetscLogEventEnd(IB_CLEANUP, 0, 0, 0, 0);

   return 0;
}


/**
 * Save position of IB to HDF5
 */
template <int dim>
PetscErrorCode ImmersedBoundary<dim>::OutputIBVar(IO::OutputHDF5 *logger, DM *da)
{
   PetscLogEventBegin(IB_OUTPUT, 0, 0, 0, 0);

   // Output Lagrangian IB Values
   Vec PointIds, IBPoints, ForceDensity;
   PetscInt i;
   DM ib_da1,ib_da2;
   double id;
   PetscMPIInt size;
   ierr = MPI_Comm_size(PETSC_COMM_WORLD, &size); CHKERRQ(ierr);
   PetscInt lx[size];
   PetscInt lvector[size];

   int ib_count = pointDict.size();
   ierr = MPI_Allgather(&ib_count, 1, MPIU_INT, lvector, 1, MPIU_INT, PETSC_COMM_WORLD);CHKERRQ(ierr);
   
   int move_values = 0;
   for(int i = 0; i< size; i++) {
      lx[i] = lvector[i];
      if (lx[i] == 0) {
         move_values++;
         lvector[i] = 1;
      }
   }
   
   if( move_values > 0 ) {
      for(int i = 0; i< size; i++) {
         if( lvector[i] > move_values) {
           lvector[i] -= move_values;
           move_values = 0;
           break;
         }
      }
   }
   
   ierr = DMDACreate1d(PETSC_COMM_WORLD, DMDA_BOUNDARY_NONE, Util::Sum(lvector, size), dim, 0, lvector, &ib_da1);
   ierr = DMDACreate1d(PETSC_COMM_WORLD, DMDA_BOUNDARY_NONE, Util::Sum(lvector, size), 1, 0, lvector, &ib_da2);

   ierr = DMCreateGlobalVector(ib_da1, &IBPoints); CHKERRQ(ierr);
   ierr = PetscObjectSetName((PetscObject) IBPoints, "IB Position" );
   ierr = DMCreateGlobalVector(ib_da1, &ForceDensity); CHKERRQ(ierr);
   ierr = PetscObjectSetName((PetscObject) ForceDensity, "IB Force Density" );
   ierr = DMCreateGlobalVector(ib_da2, &PointIds); CHKERRQ(ierr);
   ierr = PetscObjectSetName((PetscObject) PointIds, "IB PointIds");

   i = (rank > 0 ? Util::Sum(lx, rank) : 0);
   for (typename std::unordered_map<int,Point>::iterator it = pointDict.begin(); it != pointDict.end();it++)
   {
       id = ((it->second).PointID);
       ierr = VecSetValues(PointIds,1,&i,&id,INSERT_VALUES); CHKERRQ(ierr);
       i++;
   }

   i = (rank > 0 ? dim*(Util::Sum(lx, rank)) : 0);
   for (typename std::unordered_map<int,Point>::iterator it = pointDict.begin(); it != pointDict.end();it++)
   {
       for(int d = 0; d < dim; d++)
       {
           ierr = VecSetValues(IBPoints,1,&i,&((it->second).X[d]),INSERT_VALUES); CHKERRQ(ierr);
           ierr = VecSetValues(ForceDensity,1,&i,&((it->second).ForceDensity[d]),INSERT_VALUES); CHKERRQ(ierr);
           i++;
       }
   }

   ierr = VecAssemblyBegin(PointIds); CHKERRQ(ierr);
   ierr = VecAssemblyBegin(IBPoints); CHKERRQ(ierr);
   ierr = VecAssemblyBegin(ForceDensity); CHKERRQ(ierr);
   ierr = VecAssemblyEnd(PointIds);
   ierr = VecAssemblyEnd(IBPoints); CHKERRQ(ierr);
   ierr = VecAssemblyEnd(ForceDensity); CHKERRQ(ierr);

   ierr = logger->SaveIBPoints( PointIds, IBPoints, ForceDensity  ); CHKERRQ(ierr);

   VecDestroy(&PointIds);
   VecDestroy(&IBPoints);
   VecDestroy(&ForceDensity);
   DMDestroy(&ib_da1);
   DMDestroy(&ib_da2);

   // Project IB onto Eulerian Grid and output it
   if( !reduced_fluid_output ){
      Vec Omega;
      ierr = DMCreateGlobalVector(*da, &Omega); CHKERRQ(ierr);
      ierr = PetscObjectSetName((PetscObject) Omega, "Omega");
      ierr = ConstructOmega(da, &Omega); CHKERRQ(ierr);
      ierr = logger->SaveScalarField( Omega ); CHKERRQ(ierr);
      VecDestroy(&Omega);
   }
   PetscLogEventEnd(IB_OUTPUT, 0, 0, 0, 0);


   return 0;
}

/**
 * Construct Eulerian Vector containing position of IB
 */
template <>
PetscErrorCode ImmersedBoundary<2>::ConstructOmega(DM *da, Vec *Omega)
{
   PetscScalar ***array_Omega;
   PetscInt is, im, js, jm, ks, km;
   int index_start[2];
   int index_end[2];
   ierr = DMDAVecGetArrayDOF(*da, *Omega, &array_Omega); CHKERRQ(ierr);
   ierr = DMDAGetCorners(*da, &is, &js, &ks, &im, &jm, &km); CHKERRQ(ierr);

   for (std::unordered_map<int,Point>::iterator it = pointDict.begin(); it != pointDict.end();it++)
   {
       Point *pt = &(it->second);

       // Find Span of delta function for IB point
       index_start[0] = floor( (*pt).Xh[0]/spatial_stepsize[0] - 1 );
       index_end[0] = index_start[0] + 1;
       if( index_start[0] < is) index_start[0] = is;
       if( index_end[0] > is+im-1) index_end[0] = is+im-1;

       index_start[1] = floor( (*pt).Xh[1]/spatial_stepsize[1] - 1 );
       index_end[1] = index_start[1] + 1;
       if( index_start[1] < js) index_start[1] = js;
       if( index_end[1] > js+jm-1) index_end[1] = js+jm-1;
   
       for(int i = index_start[0]; i <= index_end[0]; i++)
       {
           for(int j = index_start[1]; j <= index_end[1]; j++)
           {
               array_Omega[j][i][0] = 1.0;
           }
       }
   }

   ierr = DMDAVecRestoreArrayDOF(*da, *Omega, &array_Omega); CHKERRQ(ierr);
   return 0;
}

template <>
PetscErrorCode ImmersedBoundary<3>::ConstructOmega(DM *da, Vec *Omega)
{
   PetscScalar ****array_Omega;
   PetscInt is, im, js, jm, ks, km;
   int index_start[3];
   int index_end[3];
   ierr = DMDAVecGetArrayDOF(*da, *Omega, &array_Omega); CHKERRQ(ierr);
   ierr = DMDAGetCorners(*da, &is, &js, &ks, &im, &jm, &km); CHKERRQ(ierr);

   for (std::unordered_map<int,Point>::iterator it = pointDict.begin(); it != pointDict.end();it++)
   {
       Point *pt = &(it->second);

       // Find Span of delta function for IB point
       index_start[0] = floor( (*pt).Xh[0]/spatial_stepsize[0] - 1 );
       index_end[0] = index_start[0] + 1;
       if( index_start[0] < is) index_start[0] = is;
       if( index_end[0] > is+im-1) index_end[0] = is+im-1;

       index_start[1] = floor( (*pt).Xh[1]/spatial_stepsize[1] - 1 );
       index_end[1] = index_start[1] + 1;
       if( index_start[1] < js) index_start[1] = js;
       if( index_end[1] > js+jm-1) index_end[1] = js+jm-1;

       index_start[2] = floor( (*pt).Xh[2]/spatial_stepsize[2] - 1 );
       index_end[2] = index_start[2] + 1;   
       if( index_start[2] < ks) index_start[2] = ks;
       if( index_end[2] > ks+km-1) index_end[2] = ks+km-1;
   
       for(int i = index_start[0]; i <= index_end[0]; i++)
       {
           for(int j = index_start[1]; j <= index_end[1]; j++)
           {
              for(int k = index_start[2]; k <= index_end[2]; k++)
              {
                   array_Omega[k][j][i][0] = 1.0;
              }
           }
       }
   }

   ierr = DMDAVecRestoreArrayDOF(*da, *Omega, &array_Omega); CHKERRQ(ierr);
   return 0;
}

template ImmersedBoundary<2>::ImmersedBoundary();
template ImmersedBoundary<2>::~ImmersedBoundary();
template PetscErrorCode ImmersedBoundary<2>::Initialize( std::vector<int> point_buffer, std::vector<int> force_buffer, DM *da_vector, double *h, double *domain_length, PetscInt _delta_span, double _dt );
template PetscErrorCode ImmersedBoundary<2>::AddLocalIBPoints( const int *buffer, const int buffer_length, int *last_index, const int exclude_not_in_domain );
template PetscErrorCode ImmersedBoundary<2>::AddLocalIBForceConnections( const int *buffer, const int buffer_length, const int exclude_not_in_domain  );
template PetscErrorCode ImmersedBoundary<2>::SetIBPointDomainFlags(Point *ptr_point);
template PetscErrorCode ImmersedBoundary<2>::DistributeIBPoints(DM *da_vector);
template PetscErrorCode ImmersedBoundary<2>::CleanDataNotInDomain(unsigned int keep_nextstep_points);
template PetscErrorCode ImmersedBoundary<2>::OutputIBVar(IO::OutputHDF5 *logger, DM *da);
template PetscErrorCode ImmersedBoundary<2>::UpdateIBPosition(const DM *da_vector, const Vec *Umac, const int initialize);
template PetscErrorCode ImmersedBoundary<2>::UpdateIBPosition(const DM *da_vector, const Vec *Umac, const int initialize, const int interpolate_nextstep, CallUpdateIBData update);
template PetscErrorCode ImmersedBoundary<2>::SpreadForce(DM *da_vector, Vec *F, const double time);
template PetscErrorCode ImmersedBoundary<2>::SpreadForce(DM *da_vector, Vec *F, const double time, const int zero_force, const double magnify_force);
template struct IBPointStruct<2>;

template ImmersedBoundary<3>::ImmersedBoundary();
template ImmersedBoundary<3>::~ImmersedBoundary();
template PetscErrorCode ImmersedBoundary<3>::Initialize( std::vector<int> point_buffer, std::vector<int> force_buffer, DM *da_vector, double *h, double *domain_length, PetscInt _delta_span, double _dt );
template PetscErrorCode ImmersedBoundary<3>::AddLocalIBPoints( const int *buffer, const int buffer_length, int *last_index, const int exclude_not_in_domain );
template PetscErrorCode ImmersedBoundary<3>::AddLocalIBForceConnections( const int *buffer, const int buffer_length, const int exclude_not_in_domain  );
template PetscErrorCode ImmersedBoundary<3>::SetIBPointDomainFlags(Point *ptr_point);
template PetscErrorCode ImmersedBoundary<3>::DistributeIBPoints(DM *da_vector);
template PetscErrorCode ImmersedBoundary<3>::CleanDataNotInDomain(unsigned int keep_nextstep_points);
template PetscErrorCode ImmersedBoundary<3>::OutputIBVar(IO::OutputHDF5 *logger, DM *da);
template PetscErrorCode ImmersedBoundary<3>::UpdateIBPosition(const DM *da_vector, const Vec *Umac, const int initialize);
template PetscErrorCode ImmersedBoundary<3>::UpdateIBPosition(const DM *da_vector, const Vec *Umac, const int initialize, const int interpolate_nextstep, CallUpdateIBData update);
template PetscErrorCode ImmersedBoundary<3>::SpreadForce(DM *da_vector, Vec *F, const double time);
template PetscErrorCode ImmersedBoundary<3>::SpreadForce(DM *da_vector, Vec *F, const double time, const int zero_force, const double magnify_force);
template struct IBPointStruct<3>;

/**
 * The 1d discrete delta function
 **/
double inline phi( const double r)
{
    double abs_r = fabs(r);

    if(abs_r < 1) {
        return 0.125*( 3.0 - 2.0*abs_r + sqrt(1.0 + 4.0*abs_r - 4.0*abs_r*abs_r) );
    }
    else if(abs_r < 2) {
        return 0.125*( 5.0 - 2.0*abs_r - sqrt(-7.0 + 12.0*abs_r - 4.0*abs_r*abs_r) );
    }
    else {
        return 0.0;
    }

    return 0.0;
}

/**
 * Determine whether IB point and/or force point should be distributed to a 
 * particular neighbouring node.
 *
 * Parameters:
 *        dir: Integer representing the direction of neighbour. 
 *             The number corresponds to the output of petsc's DMDAGetNeighbors.
 *        X: The location of the immersed boundary point.
 */
template <>
int ImmersedBoundary<2>::DisributeToNeighbour(const int dir, const double *X)
{
    switch ( dir )
    {
        case 0:
            if( X[0] <= xinsideghoststart[0] &&  X[1] <= xinsideghoststart[1] ) return IBDISTRIBUTE_SENDFORCEANDIB;
            else if( X[0] <= extended_xinsideghoststart[0] &&  X[1] <= extended_xinsideghoststart[1] ) return IBDISTRIBUTE_SENDIB;
            else return IBDISTRIBUTE_NOTHING;
        case 1:
            if( xghoststart[0] <= X[0] && X[0] <= xghostend[0] && X[1] <= xinsideghoststart[1] ) return IBDISTRIBUTE_SENDFORCEANDIB;
            else if( xghoststart[0] <= X[0] && X[0] <= xghostend[0] && X[1] <= extended_xinsideghoststart[1] ) return IBDISTRIBUTE_SENDIB;
            else return IBDISTRIBUTE_NOTHING;
        case 2:
            if( xinsideghostend[0] <= X[0] &&  X[1] <= xinsideghoststart[1] ) return IBDISTRIBUTE_SENDFORCEANDIB;
            else if( extended_xinsideghostend[0] <= X[0] &&  X[1] <= extended_xinsideghoststart[1] ) return IBDISTRIBUTE_SENDIB;
            else return IBDISTRIBUTE_NOTHING;
        case 3:
            if( X[0] <= xinsideghoststart[0] &&  xghoststart[1] <= X[1] && X[1] <= xghostend[1] ) return IBDISTRIBUTE_SENDFORCEANDIB;
            else if( X[0] <= extended_xinsideghoststart[0] &&  xghoststart[1] <= X[1] && X[1] <= xghostend[1] ) return IBDISTRIBUTE_SENDIB;
        case 4:
            return IBDISTRIBUTE_NOTHING;
        case 5:
            if( xinsideghostend[0] <= X[0] &&  xghoststart[1] <= X[1] &&  X[1] <= xghostend[1] ) return IBDISTRIBUTE_SENDFORCEANDIB;
            else if( extended_xinsideghostend[0] <= X[0] &&  xghoststart[1] <= X[1] &&  X[1] <= xghostend[1] ) return IBDISTRIBUTE_SENDIB;
            else return IBDISTRIBUTE_NOTHING;
        case 6:
            if( X[0] <= xinsideghoststart[0] &&  xinsideghostend[1] <= X[1] ) return IBDISTRIBUTE_SENDFORCEANDIB;
            else if( X[0] <= extended_xinsideghoststart[0] &&  extended_xinsideghostend[1] <= X[1] ) return IBDISTRIBUTE_SENDIB;
            else return IBDISTRIBUTE_NOTHING;
        case 7:
            if( xghoststart[0] <= X[0] && X[0] <= xghostend[0] &&  xinsideghostend[1] <= X[1] ) return IBDISTRIBUTE_SENDFORCEANDIB;
            else if( xghoststart[0] <= X[0] && X[0] <= xghostend[0] &&  extended_xinsideghostend[1] <= X[1] ) return IBDISTRIBUTE_SENDIB;
            else return IBDISTRIBUTE_NOTHING;
        case 8:
            if( xinsideghostend[0] <= X[0]  &&  xinsideghostend[1] <= X[1] ) return IBDISTRIBUTE_SENDFORCEANDIB;
            else if( extended_xinsideghostend[0] <= X[0]  &&  extended_xinsideghostend[1] <= X[1] ) return IBDISTRIBUTE_SENDIB;
            else return IBDISTRIBUTE_NOTHING;
        default:
            return IBDISTRIBUTE_NOTHING;
    }
}


template <>
int ImmersedBoundary<3>::DisributeToNeighbour(const int dir, const double *X)
{
    switch ( dir )
    {
        case 0:
            if( X[0] <= xinsideghoststart[0] &&  X[1] <= xinsideghoststart[1] &&  X[2] <= xinsideghoststart[2] ) return IBDISTRIBUTE_SENDFORCEANDIB;
            else if( X[0] <= extended_xinsideghoststart[0] &&  X[1] <= extended_xinsideghoststart[1] &&  X[2] <= extended_xinsideghoststart[2] ) return IBDISTRIBUTE_SENDIB;
            else return IBDISTRIBUTE_NOTHING;
        case 1:
            if( xghoststart[0] <= X[0] && X[0] <= xghostend[0] && X[1] <= xinsideghoststart[1] && X[2] <= xinsideghoststart[2] ) return IBDISTRIBUTE_SENDFORCEANDIB;
            else if( xghoststart[0] <= X[0] && X[0] <= xghostend[0] && X[1] <= extended_xinsideghoststart[1] && X[2] <= extended_xinsideghoststart[2] ) return IBDISTRIBUTE_SENDIB;
            else return IBDISTRIBUTE_NOTHING;
        case 2:
            if( xinsideghostend[0] <= X[0] &&  X[1] <= xinsideghoststart[1] &&  X[2] <= xinsideghoststart[2] ) return IBDISTRIBUTE_SENDFORCEANDIB;
            else if( extended_xinsideghostend[0] <= X[0] &&  X[1] <= extended_xinsideghoststart[1] &&  X[2] <= extended_xinsideghoststart[2] ) return IBDISTRIBUTE_SENDIB;
            else return IBDISTRIBUTE_NOTHING;
        case 3:
            if( X[0] <= xinsideghoststart[0] &&  xghoststart[1] <= X[1] && X[1] <= xghostend[1] && X[2] <= xinsideghoststart[2] ) return IBDISTRIBUTE_SENDFORCEANDIB;
            else if( X[0] <= extended_xinsideghoststart[0] &&  xghoststart[1] <= X[1] && X[1] <= xghostend[1] && X[2] <= extended_xinsideghoststart[2] ) return IBDISTRIBUTE_SENDIB;
            else return IBDISTRIBUTE_NOTHING;
        case 4:
            if( xghoststart[0] <= X[0] && X[0] <= xghostend[0] && xghoststart[1] <= X[1] && X[1] <= xghostend[1] && X[2] <= xinsideghoststart[2] ) return IBDISTRIBUTE_SENDFORCEANDIB;
            else if( xghoststart[0] <= X[0] && X[0] <= xghostend[0] && xghoststart[1] <= X[1] && X[1] <= xghostend[1] && X[2] <= extended_xinsideghoststart[2] ) return IBDISTRIBUTE_SENDIB;
            else return IBDISTRIBUTE_NOTHING;
        case 5:
            if( xinsideghostend[0] <= X[0] &&  xghoststart[1] <= X[1] &&  X[1] <= xghostend[1] && X[2] <= xinsideghoststart[2] ) return IBDISTRIBUTE_SENDFORCEANDIB;
            else if( extended_xinsideghostend[0] <= X[0] &&  xghoststart[1] <= X[1] &&  X[1] <= xghostend[1] && X[2] <= extended_xinsideghoststart[2] ) return IBDISTRIBUTE_SENDIB;
            else return IBDISTRIBUTE_NOTHING;
        case 6:
            if( X[0] <= xinsideghoststart[0] &&  xinsideghostend[1] <= X[1] && X[2] <= xinsideghoststart[2] ) return IBDISTRIBUTE_SENDFORCEANDIB;
            else if( X[0] <= extended_xinsideghoststart[0] &&  extended_xinsideghostend[1] <= X[1] && X[2] <= extended_xinsideghoststart[2] ) return IBDISTRIBUTE_SENDIB;
            else return IBDISTRIBUTE_NOTHING;
        case 7:
            if( xghoststart[0] <= X[0] && X[0] <= xghostend[0] &&  xinsideghostend[1] <= X[1] && X[2] <= xinsideghoststart[2] ) return IBDISTRIBUTE_SENDFORCEANDIB;
            else if( xghoststart[0] <= X[0] && X[0] <= xghostend[0] &&  extended_xinsideghostend[1] <= X[1] && X[2] <= extended_xinsideghoststart[2] ) return IBDISTRIBUTE_SENDIB;
            else return IBDISTRIBUTE_NOTHING;
        case 8:
            if( xinsideghostend[0] <= X[0]  &&  xinsideghostend[1] <= X[1] && X[2] <= xinsideghoststart[2] ) return IBDISTRIBUTE_SENDFORCEANDIB;
            else if( extended_xinsideghostend[0] <= X[0]  &&  extended_xinsideghostend[1] <= X[1] && X[2] <= extended_xinsideghoststart[2] ) return IBDISTRIBUTE_SENDIB;
            else return IBDISTRIBUTE_NOTHING;

        case 9:
            if( X[0] <= xinsideghoststart[0] &&  X[1] <= xinsideghoststart[1] && xghoststart[2] <= X[2] && X[2] <= xghostend[2] ) return IBDISTRIBUTE_SENDFORCEANDIB;
            else if( X[0] <= extended_xinsideghoststart[0] &&  X[1] <= extended_xinsideghoststart[1] && xghoststart[2] <= X[2] && X[2] <= xghostend[2] ) return IBDISTRIBUTE_SENDIB;
            else return IBDISTRIBUTE_NOTHING;
        case 10:
            if( xghoststart[0] <= X[0] && X[0] <= xghostend[0] && X[1] <= xinsideghoststart[1] && xghoststart[2] <= X[2] && X[2] <= xghostend[2] ) return IBDISTRIBUTE_SENDFORCEANDIB;
            else if( xghoststart[0] <= X[0] && X[0] <= xghostend[0] && X[1] <= extended_xinsideghoststart[1] && xghoststart[2] <= X[2] && X[2] <= xghostend[2] ) return IBDISTRIBUTE_SENDIB;
            else return IBDISTRIBUTE_NOTHING;
        case 11:
            if( xinsideghostend[0] <= X[0] &&  X[1] <= xinsideghoststart[1] && xghoststart[2] <= X[2] && X[2] <= xghostend[2] ) return IBDISTRIBUTE_SENDFORCEANDIB;
            else if( extended_xinsideghostend[0] <= X[0] &&  X[1] <= extended_xinsideghoststart[1] && xghoststart[2] <= X[2] && X[2] <= xghostend[2] ) return IBDISTRIBUTE_SENDIB;
            else return IBDISTRIBUTE_NOTHING;
        case 12:
            if( X[0] <= xinsideghoststart[0] &&  xghoststart[1] <= X[1] && X[1] <= xghostend[1] && xghoststart[2] <= X[2] && X[2] <= xghostend[2] ) return IBDISTRIBUTE_SENDFORCEANDIB;
            else if( X[0] <= extended_xinsideghoststart[0] &&  xghoststart[1] <= X[1] && X[1] <= xghostend[1] && xghoststart[2] <= X[2] && X[2] <= xghostend[2] ) return IBDISTRIBUTE_SENDIB;
            else return IBDISTRIBUTE_NOTHING;
        case 13:
            return IBDISTRIBUTE_NOTHING;
        case 14:
            if( xinsideghostend[0] <= X[0] &&  xghoststart[1] <= X[1] &&  X[1] <= xghostend[1] && xghoststart[2] <= X[2] && X[2] <= xghostend[2] ) return IBDISTRIBUTE_SENDFORCEANDIB;
            else if( extended_xinsideghostend[0] <= X[0] &&  xghoststart[1] <= X[1] &&  X[1] <= xghostend[1] && xghoststart[2] <= X[2] && X[2] <= xghostend[2] ) return IBDISTRIBUTE_SENDIB;
            else return IBDISTRIBUTE_NOTHING;
        case 15:
            if( X[0] <= xinsideghoststart[0] &&  xinsideghostend[1] <= X[1] && xghoststart[2] <= X[2] && X[2] <= xghostend[2] ) return IBDISTRIBUTE_SENDFORCEANDIB;
            else if( X[0] <= extended_xinsideghoststart[0] &&  extended_xinsideghostend[1] <= X[1] && xghoststart[2] <= X[2] && X[2] <= xghostend[2] ) return IBDISTRIBUTE_SENDIB;
            else return IBDISTRIBUTE_NOTHING;
        case 16:
            if( xghoststart[0] <= X[0] && X[0] <= xghostend[0] &&  xinsideghostend[1] <= X[1] && xghoststart[2] <= X[2] && X[2] <= xghostend[2] ) return IBDISTRIBUTE_SENDFORCEANDIB;
            else if( xghoststart[0] <= X[0] && X[0] <= xghostend[0] &&  extended_xinsideghostend[1] <= X[1] && xghoststart[2] <= X[2] && X[2] <= xghostend[2] ) return IBDISTRIBUTE_SENDIB;
            else return IBDISTRIBUTE_NOTHING;
        case 17:
            if( xinsideghostend[0] <= X[0]  &&  xinsideghostend[1] <= X[1] && xghoststart[2] <= X[2] && X[2] <= xghostend[2] ) return IBDISTRIBUTE_SENDFORCEANDIB;
            else if( extended_xinsideghostend[0] <= X[0]  &&  extended_xinsideghostend[1] <= X[1] && xghoststart[2] <= X[2] && X[2] <= xghostend[2] ) return IBDISTRIBUTE_SENDIB;
            else return IBDISTRIBUTE_NOTHING;

        case 18:
            if( X[0] <= xinsideghoststart[0] &&  X[1] <= xinsideghoststart[1] && xinsideghostend[2] <= X[2] ) return IBDISTRIBUTE_SENDFORCEANDIB;
            else if( X[0] <= extended_xinsideghoststart[0] &&  X[1] <= extended_xinsideghoststart[1] && extended_xinsideghostend[2] <= X[2] ) return IBDISTRIBUTE_SENDIB;
            else return IBDISTRIBUTE_NOTHING;
        case 19:
            if( xghoststart[0] <= X[0] && X[0] <= xghostend[0] && X[1] <= xinsideghoststart[1] && xinsideghostend[2] <= X[2] ) return IBDISTRIBUTE_SENDFORCEANDIB;
            else if( xghoststart[0] <= X[0] && X[0] <= xghostend[0] && X[1] <= extended_xinsideghoststart[1] && extended_xinsideghostend[2] <= X[2] ) return IBDISTRIBUTE_SENDIB;
            else return IBDISTRIBUTE_NOTHING;
        case 20:
            if( xinsideghostend[0] <= X[0] &&  X[1] <= xinsideghoststart[1] && xinsideghostend[2] <= X[2] ) return IBDISTRIBUTE_SENDFORCEANDIB;
            else if( extended_xinsideghostend[0] <= X[0] &&  X[1] <= extended_xinsideghoststart[1] && extended_xinsideghostend[2] <= X[2] ) return IBDISTRIBUTE_SENDIB;
            else return IBDISTRIBUTE_NOTHING;
        case 21:
            if( X[0] <= xinsideghoststart[0] &&  xghoststart[1] <= X[1] && X[1] <= xghostend[1] && xinsideghostend[2] <= X[2] ) return IBDISTRIBUTE_SENDFORCEANDIB;
            else if( X[0] <= extended_xinsideghoststart[0] &&  xghoststart[1] <= X[1] && X[1] <= xghostend[1] && extended_xinsideghostend[2] <= X[2] ) return IBDISTRIBUTE_SENDIB;
            else return IBDISTRIBUTE_NOTHING;
        case 22:
            if( xghoststart[0] <= X[0] && X[0] <= xghostend[0] && xghoststart[1] <= X[1] && X[1] <= xghostend[1] && xinsideghostend[2] <= X[2] ) return IBDISTRIBUTE_SENDFORCEANDIB;
            else if( xghoststart[0] <= X[0] && X[0] <= xghostend[0] && xghoststart[1] <= X[1] && X[1] <= xghostend[1] && extended_xinsideghostend[2] <= X[2] ) return IBDISTRIBUTE_SENDIB;
            else return IBDISTRIBUTE_NOTHING;
        case 23:
            if( xinsideghostend[0] <= X[0] &&  xghoststart[1] <= X[1] &&  X[1] <= xghostend[1] && xinsideghostend[2] <= X[2] ) return IBDISTRIBUTE_SENDFORCEANDIB;
            else if( extended_xinsideghostend[0] <= X[0] &&  xghoststart[1] <= X[1] &&  X[1] <= xghostend[1] && extended_xinsideghostend[2] <= X[2] ) return IBDISTRIBUTE_SENDIB;
            else return IBDISTRIBUTE_NOTHING;
        case 24:
            if( X[0] <= xinsideghoststart[0] &&  xinsideghostend[1] <= X[1] && xinsideghostend[2] <= X[2] ) return IBDISTRIBUTE_SENDFORCEANDIB;
            else if( X[0] <= extended_xinsideghoststart[0] &&  extended_xinsideghostend[1] <= X[1] && extended_xinsideghostend[2] <= X[2] ) return IBDISTRIBUTE_SENDIB;
            else return IBDISTRIBUTE_NOTHING;
        case 25:
            if( xghoststart[0] <= X[0] && X[0] <= xghostend[0] &&  xinsideghostend[1] <= X[1] && xinsideghostend[2] <= X[2] ) return IBDISTRIBUTE_SENDFORCEANDIB;
            else if( xghoststart[0] <= X[0] && X[0] <= xghostend[0] &&  extended_xinsideghostend[1] <= X[1] && extended_xinsideghostend[2] <= X[2] ) return IBDISTRIBUTE_SENDIB;
            else return IBDISTRIBUTE_NOTHING;
        case 26:
            if( xinsideghostend[0] <= X[0]  &&  xinsideghostend[1] <= X[1] && xinsideghostend[2] <= X[2] ) return IBDISTRIBUTE_SENDFORCEANDIB;
            else if( extended_xinsideghostend[0] <= X[0]  &&  extended_xinsideghostend[1] <= X[1] && extended_xinsideghostend[2] <= X[2] ) return IBDISTRIBUTE_SENDIB;
            else return IBDISTRIBUTE_NOTHING;

        default:
            return IBDISTRIBUTE_NOTHING;
    }
}

}

