/**
 * Handles the data structure and force spreading calculations for the
 * elastic force connections.
 *
 * Boundary Conditions: For PointId < 0, the point is assumed to be fictitious.
 *
 * NOTE: When adding new force connection, you need to update ../../solvers/cythonlib/IBHelper.pxi
 * to handle conversion from Python Dictionary to C IB Data structure.
 *
 * Modified: Saeed Mirazimi, last: Jan 17
 **/

#include "BendingForce.h"
#include "ImmersedBoundary.h"

namespace IB {
namespace ForceConnection {

/**
 * Constructor and Deconstructor for Force Connections
 */
template <int dim>
BendingForce<dim>::BendingForce()
{
   // Get Rank and Size
   MPI_Comm_size(PETSC_COMM_WORLD, &size);
   MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
}

template <int dim>
BendingForce<dim>::~BendingForce()
{

}


/**
 * Calculate the force density
 *
 * Parameters:
 *        pointDict:    	Associate array containing all immersed boundary points.
 *        domain_length:	Length of the entire domain
 *        time:		The time at the current time-step
 */
template <int dim>
PetscErrorCode BendingForce<dim>::CalculateForceDensity(std::unordered_map<int, Point > *pointDict, const double *domain_length, const double *spatial_stepsize, const double dt, const double time)
{
   double lX[dim], rX[dim], l2X[dim], r2X[dim], Cl[dim], Clm1[dim], Clp1[dim], ForceDensity[dim];

   for (typename std::vector<bendingforce_data>::iterator it=force_connections.begin() ; it < force_connections.end(); it++)
   {

      // initialize data for immersed boundary points
      if( (*pointDict).find(it->PointID) == (*pointDict).end() ) {
         SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_ARG_CORRUPT,"Immersed Boundary Point doesn't exist. PointID: %D, Rank: %D \n",it->PointID, rank);
         return 1;
      }
      Point *pt = &( (*pointDict)[it->PointID] );

      if( it->LPointID >= 0 && (*pointDict).find(it->LPointID) == (*pointDict).end() ){
         SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_ARG_CORRUPT,"Immersed Boundary Point doesn't exist. LPointID: %D, Rank: %D \n",it->LPointID, rank);
         return 1;
      }
      Point *lpt = ( it->LPointID >= 0 ? &( (*pointDict)[it->LPointID] ) : NULL);

      if( it->L2PointID >= 0 && (*pointDict).find(it->L2PointID) == (*pointDict).end() ){
         SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_ARG_CORRUPT,"Immersed Boundary Point doesn't exist. L2PointID: %D, Rank: %D \n",it->L2PointID, rank);
         return 1;
      }
      Point *l2pt = ( it->L2PointID >= 0 ? &( (*pointDict)[it->L2PointID] ) : NULL);

      if( it->RPointID >= 0 && (*pointDict).find(it->RPointID) == (*pointDict).end() ){
         SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_ARG_CORRUPT,"Immersed Boundary Point doesn't exist. RPointID: %D, Rank: %D \n",it->RPointID, rank);
         return 1;
      }
      Point *rpt = ( it->RPointID >= 0 ? &( (*pointDict)[it->RPointID] ) : NULL);

      if( it->R2PointID >= 0 && (*pointDict).find(it->R2PointID) == (*pointDict).end() ){
         SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_ARG_CORRUPT,"Immersed Boundary Point doesn't exist. R2PointID: %D, Rank: %D \n",it->R2PointID, rank);
         return 1;
      }
      Point *r2pt = ( it->R2PointID >= 0 ? &( (*pointDict)[it->R2PointID] ) : NULL);

      // Align IB position correctly
      for(int i = 0; i < dim; i++)
      {
         if( lpt!=NULL ) lX[i] = ibposition( (*lpt).Xh[i], (*pt).Xh[i], domain_length[i]);
         if( rpt!=NULL ) rX[i] = ibposition( (*rpt).Xh[i], (*pt).Xh[i], domain_length[i]);
         if( l2pt!=NULL )l2X[i] = ibposition( (*l2pt).Xh[i], (*pt).Xh[i], domain_length[i]);
         if( r2pt!=NULL ) r2X[i] = ibposition( (*r2pt).Xh[i], (*pt).Xh[i], domain_length[i]);
      }

      // Calculate Fictious End Points using boundary condition
      if( lpt==NULL || rpt==NULL || l2pt==NULL || r2pt==NULL )
      {
         // Error handling
         if( lpt==NULL && l2pt!=NULL ){
            SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_CORRUPT,"If LPointID is a fictitious point, then L2PointID must also be a fictitious point.\n");
            return 1;
         }
         if( rpt==NULL && r2pt!=NULL ){
            SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_CORRUPT,"If RPointID is a fictitious point, then R2PointID must also be a fictitious point.\n");
            return 1;
         }
         if( rpt==NULL && lpt==NULL ){
            SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_CORRUPT,"LPointID and RPointID can not both be fictitious points.\n");
            return 1;
         }

         // Fill End Points
         /**
         for(int i = 0; i < dim; i++)
         {
             if( l2pt==NULL ) l2X[i] = 3*(*pt).Xh[i] - 2*rX[i];
             if( lpt==NULL ) lX[i] = 2*(*pt).Xh[i] - rX[i];
             if( r2pt==NULL ) r2X[i] = 3*(*pt).Xh[i] - 2*lX[i];
             if( rpt==NULL ) rX[i] = 2*(*pt).Xh[i] - lX[i];
         }
         **/
      }


      // Calculate Force Density
      double ds = it->ds;
      double ds_sqr = ds*ds;
      double r_ds = 1.0/ds;
      double val = it->sigma * r_ds * r_ds;

      // Kludge!!! Bending stiffness is REALLY unstable. Code breaks when
      // immersed boundary moves through periodic domain. This is "fixed"
      // by using a smaller sigma for points at the edge of the boundary/
      /**
      for(int i = 0; i < dim; i++) {
          double up_crit = ((domain_length[i])) - ((spatial_stepsize[i]));
          double down_crit = ((spatial_stepsize[i]));
          if( pt->Xh[i] > up_crit || pt->Xh[i] < down_crit ) {
              val = 0.0;
              break;
          }
      }**/

      // if( lpt!=NULL && rpt!=NULL ) {
      //      for(int i = 0; i < dim; i++) Cl[i] = (lX[i] - 2*pt->Xh[i] + rX[i])*r_ds*r_ds - (it->D2X_TARGET[i]);
      // } else {
      //      for(int i = 0; i < dim; i++) Cl[i] = 0.0;
      // }

      // if( lpt!=NULL && l2pt!=NULL ) {
      //      for(int i = 0; i < dim; i++) Clm1[i] = (l2X[i] - 2*lX[i] + pt->Xh[i])*r_ds*r_ds - (it->D2X_LTARGET[i]);
      // } else {
      //      for(int i = 0; i < dim; i++) Clm1[i] = 0.0;
      // }

      // if( rpt!=NULL && r2pt!=NULL ) {
      //      for(int i = 0; i < dim; i++) Clp1[i] = (r2X[i] - 2*rX[i] + pt->Xh[i])*r_ds*r_ds - (it->D2X_RTARGET[i]);
      // } else {
      //      for(int i = 0; i < dim; i++) Clp1[i] = 0.0;
      // }
      double Rdlpt[dim], Ldlpt[dim], Rdpt[dim], Ldpt[dim], Rdrpt[dim], Ldrpt[dim];


      if(lpt == NULL && l2pt == NULL){

	Point *ptFict = &((*pointDict)[IBsize-1]);
	Point *lptFict = &((*pointDict)[IBsize-2]);
	Point *l2ptFict = &((*pointDict)[IBsize-3]);

	double lXFict[dim], l2XFict[dim];

	for(int i = 0; i < dim; i++){
	  lXFict[i] = ibposition( (*lptFict).Xh[i], (*ptFict).Xh[i], domain_length[i]);
	  l2XFict[i] = ibposition( (*l2ptFict).Xh[i], (*ptFict).Xh[i], domain_length[i]);
	}

	for(int i = 0; i < dim; i++){
	  Rdlpt[i] = (*ptFict).Xh[i] - lXFict[i];
	  Ldlpt[i] = l2XFict[i] - lXFict[i];

	  Rdpt[i] = rX[i] - (*pt).Xh[i];
	  Ldpt[i] = lXFict[i] - (*ptFict).Xh[i];

	  Rdrpt[i] = r2X[i] - rX[i];
	  Ldrpt[i] = (*pt).Xh[i] - rX[i];
	}

      }else if(rpt == NULL && r2pt == NULL){

	Point *ptFict = &((*pointDict)[0]);
	Point *rptFict = &((*pointDict)[1]);
	Point *r2ptFict = &((*pointDict)[2]);

	double rXFict[dim], r2XFict[dim];

	for(int i = 0; i < dim; i++){
	  rXFict[i] = ibposition( (*rptFict).Xh[i], (*ptFict).Xh[i], domain_length[i]);
	  r2XFict[i] = ibposition( (*r2ptFict).Xh[i], (*ptFict).Xh[i], domain_length[i]);
	}

	for(int i = 0; i < dim; i++){
	  Rdlpt[i] = (*pt).Xh[i] - lX[i];
	  Ldlpt[i] = l2X[i] - lX[i];

	  Rdpt[i] = rXFict[i] - (*ptFict).Xh[i];
	  Ldpt[i] = lX[i] - (*pt).Xh[i];

	  Rdrpt[i] = r2XFict[i] - rXFict[i];
	  Ldrpt[i] = (*ptFict).Xh[i] - rXFict[i];
	}
      }else if(lpt != NULL && l2pt == NULL){

	Point *lptFict = &((*pointDict)[IBsize - 1]);
	Point *l2ptFict = &((*pointDict)[IBsize - 2]);

	double l2XFict[dim];

	for(int i = 0; i < dim; i++){
	  l2XFict[i] = ibposition( (*l2ptFict).Xh[i], (*lptFict).Xh[i], domain_length[i]);
	}

	for(int i = 0; i < dim; i++){
	  Rdlpt[i] = (*pt).Xh[i] - lX[i];
	  Ldlpt[i] = l2XFict[i] - (*lptFict).Xh[i];

	  Rdpt[i] = rX[i] - (*pt).Xh[i];
	  Ldpt[i] = lX[i] - (*pt).Xh[i];

	  Rdrpt[i] = r2X[i] - rX[i];
	  Ldrpt[i] = (*pt).Xh[i] - rX[i];
	}

      }else if(rpt != NULL && r2pt == NULL){

	Point *rptFict = &((*pointDict)[0]);
	Point *r2ptFict = &((*pointDict)[1]);

	double r2XFict[dim];

	for(int i = 0; i < dim; i++){
	  Rdlpt[i] = (*pt).Xh[i] - lX[i];
	  Ldlpt[i] = l2X[i] - lX[i];

	  Rdpt[i] = rX[i] - (*pt).Xh[i];
	  Ldpt[i] = lX[i] - (*pt).Xh[i];

	  Rdrpt[i] = (*r2ptFict).Xh[i] - (*rptFict).Xh[i];
	  Ldrpt[i] = (*pt).Xh[i] - rX[i];
	}

      }else{

	for(int i = 0; i < dim; i++){
	  Rdlpt[i] = (*pt).Xh[i] - lX[i];
	  Ldlpt[i] = l2X[i] - lX[i];

	  Rdpt[i] = rX[i] - (*pt).Xh[i];
	  Ldpt[i] = lX[i] - (*pt).Xh[i];

	  Rdrpt[i] = r2X[i] - rX[i];
	  Ldrpt[i] = (*pt).Xh[i] - rX[i];
	}

      }
       for(int i = 0; i < dim; i++) Cl[i] = (Rdpt[i] + Ldpt[i])*r_ds*r_ds - (it->D2X_TARGET[i]);
       for(int i = 0; i < dim; i++) Clm1[i] = (Rdlpt[i] + Ldlpt[i])*r_ds*r_ds - (it->D2X_LTARGET[i]);
       for(int i = 0; i < dim; i++) Clp1[i] = (Rdrpt[i] + Ldrpt[i])*r_ds*r_ds - (it->D2X_RTARGET[i]);

       for(int i = 0; i < dim; i++) ForceDensity[i] = -1.0*val* ( Clm1[i] - 2.0*Cl[i] + Clp1[i] );

       printf("Index: %d, ",(*it)->PointID );
       printf("Rdpt : (%f,%f), Ldpt : (%f,%f), Rdlpt : (%f,%f), Ldlpt : (%f,%f), Rdrpt: (%f,%f), Cl : (%f,%f), Clm: (%f,%f), Clp : (%f,%f), ForceDensity : (%f,%f) \n ", Rdpt[0], Rdpt[1], Ldpt[0],Ldpt[1], Rdlpt[0] ,Rdlpt[1], Ldlpt[0],Ldlpt[1], Rdrpt[0], Rdrpt[1], Cl[0], Cl[1], Clm[0], Clm[1], Clp[0], Clp[1], ForceDensity[0], ForceDensity[1] );




      // if( (*it).PointID == 100) printf("Bending: %e %e (%e %e %e) (%e %e %e %e %e)\n",ForceDensity[1], time, Clm1[1], Clp1[1], Cl[1], (*pt).Xh[1], lX[1], rX[1], l2X[1], r2X[1]);
      // if( (*it).PointID == 0) printf("%d (%e %e) (%e %e) (%e %e) (%e %e) \n", (*it).PointID, ForceDensity[0], ForceDensity[1], Clm1[0], Clm1[1], Cl[0], Cl[1], Clp1[0], Clp1[1]);
      // if( (*it).PointID == 1) printf("%d (%e %e) (%e %e) (%e %e) (%e %e) \n", (*it).PointID, ForceDensity[0], ForceDensity[1], Clm1[0], Clm1[1], Cl[0], Cl[1], Clp1[0], Clp1[1]);

      for(int i = 0; i < dim; i++) pt->ForceDensity[i] += ForceDensity[i];
   }

   return 0;
}

/**
 * Add force connection to appriopriate data structures.
 *
 * Parameters:
 *        pointDict: List of immersed boundary points.
 *        buffer: Buffer filled with all the force connections. The buffer is of
 *                      type int, so they need to be cast to correct struct type.
 *        i: current index of the buffer. move index after pulling data.
 *        exclude_not_in_domain: Flag determining if the point should be excluded if it
 *                                isn't in the proc's domain.
 */
template <int dim>
PetscErrorCode BendingForce<dim>::AddLocalIBForceConnections(std::unordered_map<int, Point > *pointDict, const int *buffer, int *i, const int exclude_not_in_domain  )
{
   bendingforce_data force = *((bendingforce_data *) &(buffer[*i]) );

   if( force_lookup.find(force.FCID) == force_lookup.end() &&
         ( exclude_not_in_domain == 0 ||
               (pointDict->find(force.PointID) != pointDict->end() && (*pointDict)[force.PointID].I == IB::IBPOINT_INDOMIAN)
         )
      )
   {
      force_connections.push_back(force);
      force_lookup[force.FCID] = 1;
   }

   *i += sizeof(bendingforce_data)/sizeof(int);

   return 0;
}

/**
 * Clean force connections not in the domain
 */
template <int dim>
PetscErrorCode BendingForce<dim>::CleanDataNotInDomain(std::unordered_map<int, Point > *pointDict)
{
   void* params[2] = {(void*)pointDict, (void*) &force_lookup};
   force_connections.erase( remove_if(force_connections.begin(), force_connections.end(), std::bind2nd(std::ptr_fun(&BendingForce<dim>::CleanupForceConnections), &params[0])) , force_connections.end());
   return 0;
}

template <int dim>
bool BendingForce<dim>::CleanupForceConnections(const bendingforce_data_template<dim> fc, void **params)
{
   std::unordered_map<int, Point > *pointDict = (std::unordered_map<int, Point > *) params[0];
   std::unordered_map<int, int > *force_lookup =  (std::unordered_map<int, int > *) params[1];

   bool missing_points = pointDict->count(fc.PointID) == 0;
   if( missing_points ) force_lookup->erase( fc.FCID );

   return missing_points;
}

/**
 * Distributed Force Connections and IB Points
 */
template <int dim>
PetscErrorCode BendingForce<dim>::DistributeIBPoints(IBPoints *ib, std::unordered_map<int, Point > *pointDict, const int dir, std::unordered_map<int, int > *distribute_points, std::vector< int > *output_buffer, std::vector< int > *force_output_buffer)
{
   Point pt;
   NetworkPoint network_pt;
   int *ptr;

   for (typename std::vector<bendingforce_data>::iterator it=force_connections.begin() ; it < force_connections.end(); it++)
   {
      // Check if node and/or force connection needs to be sent to Neighbour in direction dir
      int send = ib->DisributeToNeighbour(dir, ((*pointDict)[(*it).PointID]).X );
      int send2 = ib->DisributeToNeighbour(dir, ((*pointDict)[(*it).PointID]).Xh );
      if ( send == IB::IBDISTRIBUTE_SENDFORCEANDIB || send2 == IB::IBDISTRIBUTE_SENDFORCEANDIB) send = IB::IBDISTRIBUTE_SENDFORCEANDIB;
      else if ( send == IB::IBDISTRIBUTE_SENDIB || send2 == IB::IBDISTRIBUTE_SENDIB) send = IB::IBDISTRIBUTE_SENDIB;

      // Add IB point to buffer
      if ( send == IB::IBDISTRIBUTE_SENDIB ||  send == IB::IBDISTRIBUTE_SENDFORCEANDIB)
      {
         // Add Points to buffer
         int ids[5] = {(*it).PointID, (*it).LPointID, (*it).RPointID, (*it).L2PointID, (*it).R2PointID};

         for( int index = 0; index < 5; index++)
         {
            if( distribute_points->count( ids[index] ) == 0 &&  pointDict->count( ids[index] ) > 0)
            {
               // Load Network Point
               ptr = (int *) &network_pt;
               pt = ((*pointDict)[ids[index]]);
               IB::LoadNetworkPoint<dim>(&pt, &network_pt);

               // Insert into buffer
               for(unsigned  int i=0; i < sizeof(NetworkPoint)/sizeof(int); i++)  output_buffer->push_back( ptr[i] );
               (*distribute_points)[ ids[index] ] = 1;
            }
         }
      }

      // Add Force connection to buffer
      if ( send == IB::IBDISTRIBUTE_SENDFORCEANDIB)
      {
         // Add Force to temporary buffer
         ptr = (int *) &(*it);
         for(unsigned int i=0; i < sizeof(bendingforce_data)/sizeof(int); i++)  force_output_buffer->push_back( ptr[i] );
      }

   }

   return 0;
}

template BendingForce<2>::BendingForce();
template BendingForce<2>::~BendingForce();
template PetscErrorCode BendingForce<2>::CalculateForceDensity(std::unordered_map<int, Point > *pointDict, const double *domain_length, const double *spatial_stepsize, const double dt, const double time);
template PetscErrorCode BendingForce<2>::AddLocalIBForceConnections(std::unordered_map<int, Point > *pointDict, const int *buffer, int *i, const int exclude_not_in_domain  );
template PetscErrorCode BendingForce<2>::CleanDataNotInDomain(std::unordered_map<int, Point > *pointDict);
template PetscErrorCode BendingForce<2>::DistributeIBPoints(IBPoints *ib, std::unordered_map<int, Point > *pointDict, const int dir, std::unordered_map<int, int > *distribute_points, std::vector< int > *output_buffer, std::vector< int > *force_output_buffer);
template bool BendingForce<2>::CleanupForceConnections(const bendingforce_data_template<2>  fc, void **params);

template BendingForce<3>::BendingForce();
template BendingForce<3>::~BendingForce();
template PetscErrorCode BendingForce<3>::CalculateForceDensity(std::unordered_map<int, Point > *pointDict, const double *domain_length, const double *spatial_stepsize, const double dt, const double time);
template PetscErrorCode BendingForce<3>::AddLocalIBForceConnections(std::unordered_map<int, Point > *pointDict, const int *buffer, int *i, const int exclude_not_in_domain  );
template PetscErrorCode BendingForce<3>::CleanDataNotInDomain(std::unordered_map<int, Point > *pointDict);
template PetscErrorCode BendingForce<3>::DistributeIBPoints(IBPoints *ib, std::unordered_map<int, Point > *pointDict, const int dir, std::unordered_map<int, int > *distribute_points, std::vector< int > *output_buffer, std::vector< int > *force_output_buffer);
template bool BendingForce<3>::CleanupForceConnections(const bendingforce_data_template<3> fc, void **params);
}
}
