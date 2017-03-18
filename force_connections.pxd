"""
Define all the force connection structs used in C++ portion of code/
Note: The struct variable order needs to be identical to the C++ definition.
"""
# Updated by Saeed Mirazimi
# June Jan 17
# To handle the data structure of the the new force class, BendingForceCurve
cdef struct elasticforceenergy_data:
    int ForceID
    int FCID
    int LPointID
    int RPointID
    int PointID
    double sigma
    double L
    double dv
    double ds

cdef struct bendingforcecurve_data_3d:
    int ForceID
    int FCID
    int L2PointID
    int LPointID
    int RPointID
    int R2PointID
    int PointID
    double C2
    double epsil
    double length
    double Omega
    double Kappa
    double dv
    double ds

cdef struct bendingforcecurve_data_2d:
    int ForceID
    int FCID
    int L2PointID
    int LPointID
    int RPointID
    int R2PointID
    int PointID
    double C2
    double epsil
    double length
    double Omega
    double Kappa
    double dv
    double ds

cdef struct bendingforce_data_3d:
    int ForceID
    int FCID
    int L2PointID
    int LPointID
    int RPointID
    int R2PointID
    int PointID
    double sigma
    #double D2X_TARGET[3]
    #double D2X_LTARGET[3]
    #double D2X_RTARGET[3]
    double A
    double K
    double Omega

    double dv
    double ds

cdef struct bendingforce_data_2d:
    int ForceID
    int FCID
    int L2PointID
    int LPointID
    int RPointID
    int R2PointID
    int PointID
    double sigma
    #double D2X_TARGET[2]
    #double D2X_LTARGET[2]
    #double D2X_RTARGET[2]
    double A
    double K
    double Omega

    double dv
    double ds

cdef struct elasticforce_tether_2pt_data_2d:
    int ForceID
    int FCID
    int PointID
    double sigma
    double L
    double PtInitLoc[2]
    double PtVel[2]

cdef struct elasticforce_tether_2pt_data_3d:
    int ForceID
    int FCID
    int PointID
    double sigma
    double L
    double PtInitLoc[3]
    double PtVel[3]

cdef struct elasticforce_tether_osc_stiffness_2pt_data_2d:
    int ForceID
    int FCID
    int PointID
    double sigma
    double L
    int OscillatingTypeID
    double OscillatingFrequency
    double OscillatingAmplitude
    double onoff_ratio
    double time_offset
    double PtInitLoc[2]
    double PtVel[2]

cdef struct elasticforce_tether_osc_stiffness_2pt_data_3d:
    int ForceID
    int FCID
    int PointID
    double sigma
    double L
    int OscillatingTypeID
    double OscillatingFrequency
    double OscillatingAmplitude
    double onoff_ratio
    double time_offset
    double PtInitLoc[3]
    double PtVel[3]

cdef struct elasticforce_data:
    int ForceID
    int FCID
    int LPointID
    int RPointID
    int PointID
    double sigma
    double L
    double dv
    double ds

cdef struct elasticforce_2pt_data:
    int ForceID
    int FCID
    int PointID1
    int PointID2
    double sigma
    double L

cdef struct elasticforce_osc_stiffness_data:
    int ForceID
    int FCID
    int LPointID
    int RPointID
    int PointID
    double sigma
    double L
    double dv
    double ds
    int OscillatingTypeID
    double OscillatingFrequency
    double OscillatingAmplitude
    double onoff_ratio
    double time_offset

cdef struct elasticforce_osc_length_data:
    int ForceID
    int FCID
    int LPointID
    int RPointID
    int PointID
    double sigma
    double L
    double dv
    double ds
    int OscillatingTypeID
    double OscillatingFrequency
    double OscillatingAmplitude
    double onoff_ratio
    double time_offset

cdef struct elasticforce_osc_stiffness_2pt_data:
    int ForceID
    int FCID
    int PointID1
    int PointID2
    double sigma
    double L
    int OscillatingTypeID
    double OscillatingFrequency
    double OscillatingAmplitude
    double onoff_ratio
    double time_offset

cdef struct elasticforce_osc_length_2pt_data:
    int ForceID
    int FCID
    int PointID1
    int PointID2
    double sigma
    double L
    int OscillatingTypeID
    double OscillatingFrequency
    double OscillatingAmplitude
    double onoff_ratio
    double time_offset

cdef struct penaltymass_data_2d:
    int ForceID
    int FCID
    int PointID
    double sigma
    double mass
    double PtX[2]
    double PtVel[2]
    double PtVelOld[2]
    double dv

cdef struct penaltymass_data_3d:
    int ForceID
    int FCID
    int PointID
    double sigma
    double mass
    double PtX[3]
    double PtVel[3]
    double PtVelOld[3]
    double dv

cdef struct force_pulse_2pt_data:
    int ForceID
    int FCID
    int PointID1
    int PointID2
    double sigma
    int OscillatingTypeID
    double OscillatingFrequency
    double OscillatingAmplitude
    double onoff_ratio
    double time_offset

