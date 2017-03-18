# Updated by Saeed Mirazimi
# Jan 17
# To handle the data structure of the the new force classes
cimport cythonlib.ibpoint as ibpoint
cimport cythonlib.force_connections as force_connections

cdef int FORCE_CONNECTION_ELASTIC = 1
cdef int FORCE_CONNECTION_ELASTIC_TWOPOINT = 2
cdef int FORCE_CONNECTION_ELASTIC_OSCILLATINGSTIFFNESS = 3
cdef int FORCE_CONNECTION_ELASTIC_OSCILLATINGSTIFFNESS_TWOPOINT = 4
cdef int FORCE_CONNECTION_ELASTIC_OSCILLATINGLENGTH = 5
cdef int FORCE_CONNECTION_ELASTIC_OSCILLATINGLENGTH_TWOPOINT = 6
cdef int FORCE_CONNECTION_ELASTIC_TETHER_TWOPOINT = 7
cdef int FORCE_CONNECTION_ELASTIC_TETHER_OSCILLATINGSTIFFNESS_TWOPOINT = 8
cdef int FORCE_CONNECTION_BENDING = 9
cdef int FORCE_CONNECTION_PENALTYMASS = 10
cdef int FORCE_CONNECTION_FORCE_PULSE_TWOPOINT = 11
#Added by Saeed
cdef int FORCE_CONNECTION_BENDING_CURVE = 12
cdef int FORCE_CONNECTION_ELASTIC_ENERGY = 13


cdef ConstructIB_C_DataStructure(ibData, vector[int] *points, vector[int] *forces):
    """
    Construct C data structure that stores all the IB data.

    Parameters:
       ibData: Dictionary containing IB data.
    """

    cdef int dim = ibData["Dimensions"]

    if dim == 2:
        ConstructIB_C_DataStructure2d(ibData, points, forces)
    elif dim == 3:
        ConstructIB_C_DataStructure3d(ibData, points, forces)
    else:
        raise Exception("Invalid Dimension")


cdef ConstructIB_C_DataStructure2d(ibData, vector[int] *points, vector[int] *forces   ):
    """
    Construct C data structure that stores all the 2d IB data.

    Parameters:
       ibData: Dictionary containing IB data.
    """

    cdef int dim = ibData["Dimensions"]
    cdef ibpoint.IBPointStruct2d c_pt
    cdef int *c_intpt = <int *> &c_pt

    # Import immersed boundary points
    for pt in ibData["Points"]:
         c_pt.PointID = pt["PointID"]
         c_pt.X[0] = pt["X1"]
         c_pt.X[1] = pt["X2"]
         c_pt.Xh[0] = pt["X1"]
         c_pt.Xh[1] = pt["X2"]

         for i in range(sizeof(ibpoint.IBPointStruct2d)/sizeof(int)):
              points.push_back( c_intpt[i] )

    # Import Force Connection
    ConstructIB_C_ForceConnection(ibData, forces )


cdef ConstructIB_C_DataStructure3d(ibData, vector[int] *points, vector[int] *forces   ):
    """
    Construct C data structure that stores all the 2d IB data.

    Parameters:
       ibData: Dictionary containing IB data.
    """

    cdef int dim = ibData["Dimensions"]
    cdef ibpoint.IBPointStruct3d c_pt
    cdef int *c_intpt = <int *> &c_pt

    # Import immersed boundary points
    for pt in ibData["Points"]:
         c_pt.PointID = pt["PointID"]
         c_pt.X[0] = pt["X1"]
         c_pt.X[1] = pt["X2"]
         c_pt.X[2] = pt["X3"]
         c_pt.Xh[0] = pt["X1"]
         c_pt.Xh[1] = pt["X2"]
         c_pt.Xh[2] = pt["X3"]

         for i in range(sizeof(ibpoint.IBPointStruct3d)/sizeof(int)):
              points.push_back( c_intpt[i] )

    # Import Force Connection
    ConstructIB_C_ForceConnection(ibData, forces )


cdef ConstructIB_C_ForceConnection(ibData, vector[int] *forces   ):
    """
    Construct C data structure that stores all the 2d IB data.

    Parameters:
       ibData: Dictionary containing IB data.
    """

    # Import Force Connection
    cdef force_connections.elasticforce_data c_elasticforce
    cdef force_connections.elasticforce_2pt_data c_elasticforce_2pt
    cdef force_connections.bendingforce_data_2d c_bendingforce_2d
    cdef force_connections.bendingforce_data_3d c_bendingforce_3d
    cdef force_connections.elasticforce_osc_stiffness_data c_elasticforce_osc_stiffness
    cdef force_connections.elasticforce_osc_length_data c_elasticforce_osc_length
    cdef force_connections.elasticforce_osc_stiffness_2pt_data c_elasticforce_osc_stiffness_2pt
    cdef force_connections.elasticforce_osc_length_2pt_data c_elasticforce_osc_length_2pt
    cdef force_connections.elasticforce_tether_2pt_data_2d c_elasticforce_tether_2pt_2d
    cdef force_connections.elasticforce_tether_2pt_data_3d c_elasticforce_tether_2pt_3d
    cdef force_connections.elasticforce_tether_osc_stiffness_2pt_data_2d c_elasticforce_tether_osc_stiffness_2pt_2d
    cdef force_connections.elasticforce_tether_osc_stiffness_2pt_data_3d c_elasticforce_tether_osc_stiffness_2pt_3d
    cdef force_connections.penaltymass_data_2d c_penaltymass_data_2d
    cdef force_connections.penaltymass_data_3d c_penaltymass_data_3d
    cdef force_connections.force_pulse_2pt_data c_force_pulse_2pt
#Added by Saeed
    cdef force_connections.bendingforcecurve_data_2d c_bendingforcecurve_2d
    cdef force_connections.bendingforcecurve_data_3d c_bendingforcecurve_3d
    cdef force_connections.elasticforceenergy_data c_elasticforceenergy

    cdef int *c_int_elasticforce = <int *> &c_elasticforce
    cdef int *c_int_elasticforce_2pt = <int *> &c_elasticforce_2pt
    cdef int *c_int_bendingforce_2d = <int *> &c_bendingforce_2d
    cdef int *c_int_bendingforce_3d = <int *> &c_bendingforce_3d
    cdef int *c_int_elasticforce_osc_stiffness = <int *> &c_elasticforce_osc_stiffness
    cdef int *c_int_elasticforce_osc_length = <int *> &c_elasticforce_osc_length
    cdef int *c_int_elasticforce_osc_stiffness_2pt = <int *> &c_elasticforce_osc_stiffness_2pt
    cdef int *c_int_elasticforce_osc_length_2pt = <int *> &c_elasticforce_osc_length_2pt
    cdef int *c_int_elasticforce_tether_2pt_2d = <int *> &c_elasticforce_tether_2pt_2d
    cdef int *c_int_elasticforce_tether_2pt_3d = <int *> &c_elasticforce_tether_2pt_3d
    cdef int *c_int_elasticforce_tether_osc_stiffness_2pt_2d = <int *> &c_elasticforce_tether_osc_stiffness_2pt_2d
    cdef int *c_int_elasticforce_tether_osc_stiffness_2pt_3d = <int *> &c_elasticforce_tether_osc_stiffness_2pt_3d
    cdef int *c_int_penaltymass_data_2d = <int *> &c_penaltymass_data_2d
    cdef int *c_int_penaltymass_data_3d = <int *> &c_penaltymass_data_3d
    cdef int *c_int_force_pulse_2pt = <int *> &c_force_pulse_2pt
#Added by Saeed
    cdef int *c_int_bendingforcecurve_2d = <int *> &c_bendingforcecurve_2d
    cdef int *c_int_bendingforcecurve_3d = <int *> &c_bendingforcecurve_3d
    cdef int *c_int_elasticforceenergy = <int *> &c_elasticforceenergy

    # PointID Map
    PointIDMap = {}
    for (index, d) in enumerate(ibData["Points"]):
       PointIDMap[d["PointID"]] = index

    for fc in ibData["ForceConnections"]:
        if fc["ForceID"] == "ELASTIC":
           c_elasticforce.ForceID = FORCE_CONNECTION_ELASTIC
           c_elasticforce.FCID = fc["FCID"]
           c_elasticforce.LPointID = fc["LPointID"]
           c_elasticforce.RPointID = fc["RPointID"]
           c_elasticforce.PointID = fc["PointID"]
           c_elasticforce.sigma = fc["sigma"]
           c_elasticforce.L = fc["L"]
           c_elasticforce.ds = fc["ds"]
           c_elasticforce.dv = fc["dv"]

           for i in range(sizeof(force_connections.elasticforce_data)/sizeof(int)):
              forces.push_back( c_int_elasticforce[i] )
        elif fc["ForceID"] == "ELASTIC_ENERGY":
           c_elasticforce.ForceID = FORCE_CONNECTION_ELASTIC_ENERGY
           c_elasticforce.FCID = fc["FCID"]
           c_elasticforce.LPointID = fc["LPointID"]
           c_elasticforce.RPointID = fc["RPointID"]
           c_elasticforce.PointID = fc["PointID"]
           c_elasticforce.sigma = fc["sigma"]
           c_elasticforce.L = fc["L"]
           c_elasticforce.ds = fc["ds"]
           c_elasticforce.dv = fc["dv"]

           for i in range(sizeof(force_connections.elasticforce_data)/sizeof(int)):
              forces.push_back( c_int_elasticforce[i] )

        elif fc["ForceID"] == "BENDING" and ibData["Dimensions"] == 2:
           c_bendingforce_2d.ForceID = FORCE_CONNECTION_BENDING
           c_bendingforce_2d.FCID = fc["FCID"]
           c_bendingforce_2d.L2PointID = fc["L2PointID"]
           c_bendingforce_2d.LPointID = fc["LPointID"]
           c_bendingforce_2d.R2PointID = fc["R2PointID"]
           c_bendingforce_2d.RPointID = fc["RPointID"]
           c_bendingforce_2d.PointID = fc["PointID"]
           c_bendingforce_2d.sigma = fc["sigma"]

           #for i in range(1,ibData["Dimensions"]+1): c_bendingforce_2d.D2X_TARGET[i-1] = fc["D2X_TARGET_"+str(i)]
           #for i in range(1,ibData["Dimensions"]+1): c_bendingforce_2d.D2X_LTARGET[i-1] = fc["D2X_LTARGET_"+str(i)]
           #for i in range(1,ibData["Dimensions"]+1): c_bendingforce_2d.D2X_RTARGET[i-1] = fc["D2X_RTARGET_"+str(i)]
           c_bendingforce_2d.A = fc["A"]
           c_bendingforce_2d.K = fc["K"]
           c_bendingforce_2d.Omega = fc["Omega"]

           c_bendingforce_2d.ds = fc["ds"]
           c_bendingforce_2d.dv = fc["dv"]

           for i in range(sizeof(force_connections.bendingforce_data_2d)/sizeof(int)):
              forces.push_back( c_int_bendingforce_2d[i] )

        elif fc["ForceID"] == "BENDING" and ibData["Dimensions"] == 3:
           c_bendingforce_3d.ForceID = FORCE_CONNECTION_BENDING
           c_bendingforce_3d.FCID = fc["FCID"]
           c_bendingforce_3d.L2PointID = fc["L2PointID"]
           c_bendingforce_3d.LPointID = fc["LPointID"]
           c_bendingforce_3d.R2PointID = fc["R2PointID"]
           c_bendingforce_3d.RPointID = fc["RPointID"]
           c_bendingforce_3d.PointID = fc["PointID"]
           c_bendingforce_3d.sigma = fc["sigma"]

           #for i in range(1,ibData["Dimensions"]+1): c_bendingforce_3d.D2X_TARGET[i-1] = fc["D2X_TARGET_"+str(i)]
           #for i in range(1,ibData["Dimensions"]+1): c_bendingforce_3d.D2X_LTARGET[i-1] = fc["D2X_LTARGET_"+str(i)]
           #for i in range(1,ibData["Dimensions"]+1): c_bendingforce_3d.D2X_RTARGET[i-1] = fc["D2X_RTARGET_"+str(i)]
           c_bendingforce_2d.A = fc["A"]
           c_bendingforce_2d.K = fc["K"]
           c_bendingforce_2d.Omega = fc["Omega"]

           c_bendingforce_3d.ds = fc["ds"]
           c_bendingforce_3d.dv = fc["dv"]

           for i in range(sizeof(force_connections.bendingforce_data_3d)/sizeof(int)):
              forces.push_back( c_int_bendingforce_3d[i] )

        elif fc["ForceID"] == "ELASTIC_TWOPOINT":
           c_elasticforce_2pt.ForceID = FORCE_CONNECTION_ELASTIC_TWOPOINT
           c_elasticforce_2pt.FCID = fc["FCID"]
           c_elasticforce_2pt.PointID1 = fc["PointID1"]
           c_elasticforce_2pt.PointID2 = fc["PointID2"]
           c_elasticforce_2pt.sigma = fc["sigma"]
           c_elasticforce_2pt.L = fc["L"]

           for i in range(sizeof(force_connections.elasticforce_2pt_data)/sizeof(int)):
              forces.push_back( c_int_elasticforce_2pt[i] )

        elif fc["ForceID"] == "ELASTIC_OSCILLATINGSTIFFNESS":
           c_elasticforce_osc_stiffness.ForceID = FORCE_CONNECTION_ELASTIC_OSCILLATINGSTIFFNESS
           c_elasticforce_osc_stiffness.FCID = fc["FCID"]
           c_elasticforce_osc_stiffness.PointID = fc["PointID"]
           c_elasticforce_osc_stiffness.LPointID = fc["LPointID"]
           c_elasticforce_osc_stiffness.RPointID = fc["RPointID"]
           c_elasticforce_osc_stiffness.sigma = fc["sigma"]
           c_elasticforce_osc_stiffness.L = fc["L"]
           c_elasticforce_osc_stiffness.OscillatingTypeID = fc["OscillatingTypeID"]
           c_elasticforce_osc_stiffness.OscillatingFrequency = fc["OscillatingFrequency"]
           c_elasticforce_osc_stiffness.OscillatingAmplitude = fc["OscillatingAmplitude"]
           c_elasticforce_osc_stiffness.onoff_ratio = fc["onoff_ratio"]

           if "time_offset" in fc: c_elasticforce_osc_stiffness.time_offset = fc["time_offset"]
           else: c_elasticforce_osc_stiffness.time_offset = 0.0

           c_elasticforce_osc_stiffness.ds = fc["ds"]
           c_elasticforce_osc_stiffness.dv = fc["dv"]

           for i in range(sizeof(force_connections.elasticforce_osc_stiffness_data)/sizeof(int)):
              forces.push_back( c_int_elasticforce_osc_stiffness[i] )

        elif fc["ForceID"] == "ELASTIC_OSCILLATINGLENGTH":

           c_elasticforce_osc_length.ForceID = FORCE_CONNECTION_ELASTIC_OSCILLATINGLENGTH
           c_elasticforce_osc_length.FCID = fc["FCID"]
           c_elasticforce_osc_length.PointID = fc["PointID"]
           c_elasticforce_osc_length.LPointID = fc["LPointID"]
           c_elasticforce_osc_length.RPointID = fc["RPointID"]
           c_elasticforce_osc_length.sigma = fc["sigma"]
           c_elasticforce_osc_length.L = fc["L"]
           c_elasticforce_osc_length.OscillatingTypeID = fc["OscillatingTypeID"]
           c_elasticforce_osc_length.OscillatingFrequency = fc["OscillatingFrequency"]
           c_elasticforce_osc_length.OscillatingAmplitude = fc["OscillatingAmplitude"]
           c_elasticforce_osc_length.onoff_ratio = fc["onoff_ratio"]

           if "time_offset" in fc: c_elasticforce_osc_length.time_offset = fc["time_offset"]
           else: c_elasticforce_osc_length.time_offset = 0.0

           c_elasticforce_osc_length.ds = fc["ds"]
           c_elasticforce_osc_length.dv = fc["dv"]

           for i in range(sizeof(force_connections.elasticforce_osc_length_data)/sizeof(int)):
              forces.push_back( c_int_elasticforce_osc_length[i] )


        elif fc["ForceID"] == "ELASTIC_OSCILLATINGSTIFFNESS_TWOPOINT":
           c_elasticforce_osc_stiffness_2pt.ForceID = FORCE_CONNECTION_ELASTIC_OSCILLATINGSTIFFNESS_TWOPOINT
           c_elasticforce_osc_stiffness_2pt.FCID = fc["FCID"]
           c_elasticforce_osc_stiffness_2pt.PointID1 = fc["PointID1"]
           c_elasticforce_osc_stiffness_2pt.PointID2 = fc["PointID2"]
           c_elasticforce_osc_stiffness_2pt.sigma = fc["sigma"]
           c_elasticforce_osc_stiffness_2pt.L = fc["L"]
           c_elasticforce_osc_stiffness_2pt.OscillatingTypeID = fc["OscillatingTypeID"]
           c_elasticforce_osc_stiffness_2pt.OscillatingFrequency = fc["OscillatingFrequency"]
           c_elasticforce_osc_stiffness_2pt.OscillatingAmplitude = fc["OscillatingAmplitude"]
           c_elasticforce_osc_stiffness_2pt.onoff_ratio = fc["onoff_ratio"]

           if "time_offset" in fc: c_elasticforce_osc_stiffness_2pt.time_offset = fc["time_offset"]
           else: c_elasticforce_osc_stiffness_2pt.time_offset = 0.0

           for i in range(sizeof(force_connections.elasticforce_osc_stiffness_2pt_data)/sizeof(int)):
              forces.push_back( c_int_elasticforce_osc_stiffness_2pt[i] )

        elif fc["ForceID"] == "ELASTIC_OSCILLATINGLENGTH_TWOPOINT":
           c_elasticforce_osc_length_2pt.ForceID = FORCE_CONNECTION_ELASTIC_OSCILLATINGLENGTH_TWOPOINT
           c_elasticforce_osc_length_2pt.FCID = fc["FCID"]
           c_elasticforce_osc_length_2pt.PointID1 = fc["PointID1"]
           c_elasticforce_osc_length_2pt.PointID2 = fc["PointID2"]
           c_elasticforce_osc_length_2pt.sigma = fc["sigma"]
           c_elasticforce_osc_length_2pt.L = fc["L"]
           c_elasticforce_osc_length_2pt.OscillatingTypeID = fc["OscillatingTypeID"]
           c_elasticforce_osc_length_2pt.OscillatingFrequency = fc["OscillatingFrequency"]
           c_elasticforce_osc_length_2pt.OscillatingAmplitude = fc["OscillatingAmplitude"]
           c_elasticforce_osc_length_2pt.onoff_ratio = fc["onoff_ratio"]

           if "time_offset" in fc: c_elasticforce_osc_length_2pt.time_offset = fc["time_offset"]
           else: c_elasticforce_osc_length_2pt.time_offset = 0.0

           for i in range(sizeof(force_connections.elasticforce_osc_length_2pt_data)/sizeof(int)):
              forces.push_back( c_int_elasticforce_osc_length_2pt[i] )

        elif (fc["ForceID"] == "ELASTIC_TETHER" or fc["ForceID"] == "ELASTIC_TETHER_TWOPOINT") and ibData["Dimensions"] == 2:
           c_elasticforce_tether_2pt_2d.ForceID = FORCE_CONNECTION_ELASTIC_TETHER_TWOPOINT
           c_elasticforce_tether_2pt_2d.FCID = fc["FCID"]
           c_elasticforce_tether_2pt_2d.PointID = fc["PointID"]
           c_elasticforce_tether_2pt_2d.sigma = fc["sigma"]
           c_elasticforce_tether_2pt_2d.L = fc["L"]

           for i in range(1,ibData["Dimensions"]+1): c_elasticforce_tether_2pt_2d.PtInitLoc[i-1] = fc["PtInitLoc_"+str(i)]
           for i in range(1,ibData["Dimensions"]+1): c_elasticforce_tether_2pt_2d.PtVel[i-1] = fc["PtVel_"+str(i)]

           for i in range(sizeof(force_connections.elasticforce_tether_2pt_data_2d)/sizeof(int)):
              forces.push_back( c_int_elasticforce_tether_2pt_2d[i] )

        elif (fc["ForceID"] == "ELASTIC_TETHER" or fc["ForceID"] == "ELASTIC_TETHER_TWOPOINT") and ibData["Dimensions"] == 3:
           c_elasticforce_tether_2pt_3d.ForceID = FORCE_CONNECTION_ELASTIC_TETHER_TWOPOINT
           c_elasticforce_tether_2pt_3d.FCID = fc["FCID"]
           c_elasticforce_tether_2pt_3d.PointID = fc["PointID"]
           c_elasticforce_tether_2pt_3d.sigma = fc["sigma"]
           c_elasticforce_tether_2pt_3d.L = fc["L"]

           for i in range(1,ibData["Dimensions"]+1): c_elasticforce_tether_2pt_3d.PtInitLoc[i-1] = fc["PtInitLoc_"+str(i)]
           for i in range(1,ibData["Dimensions"]+1): c_elasticforce_tether_2pt_3d.PtVel[i-1] = fc["PtVel_"+str(i)]

           for i in range(sizeof(force_connections.elasticforce_tether_2pt_data_3d)/sizeof(int)):
              forces.push_back( c_int_elasticforce_tether_2pt_3d[i] )

        elif (fc["ForceID"] == "ELASTIC_TETHER_OSCILLATINGSTIFFNESS" or fc["ForceID"] == "ELASTIC_TETHER_OSCILLATINGSTIFFNESS_TWOPOINT") and ibData["Dimensions"] == 2:
           c_elasticforce_tether_osc_stiffness_2pt_2d.ForceID = FORCE_CONNECTION_ELASTIC_TETHER_OSCILLATINGSTIFFNESS_TWOPOINT
           c_elasticforce_tether_osc_stiffness_2pt_2d.FCID = fc["FCID"]
           c_elasticforce_tether_osc_stiffness_2pt_2d.PointID = fc["PointID"]
           c_elasticforce_tether_osc_stiffness_2pt_2d.sigma = fc["sigma"]
           c_elasticforce_tether_osc_stiffness_2pt_2d.L = fc["L"]
           c_elasticforce_tether_osc_stiffness_2pt_2d.OscillatingTypeID = fc["OscillatingTypeID"]
           c_elasticforce_tether_osc_stiffness_2pt_2d.OscillatingFrequency = fc["OscillatingFrequency"]
           c_elasticforce_tether_osc_stiffness_2pt_2d.OscillatingAmplitude = fc["OscillatingAmplitude"]
           c_elasticforce_tether_osc_stiffness_2pt_2d.onoff_ratio = fc["onoff_ratio"]

           for i in range(1,ibData["Dimensions"]+1): c_elasticforce_tether_osc_stiffness_2pt_2d.PtInitLoc[i-1] = fc["PtInitLoc_"+str(i)]
           for i in range(1,ibData["Dimensions"]+1): c_elasticforce_tether_osc_stiffness_2pt_2d.PtVel[i-1] = fc["PtVel_"+str(i)]

           for i in range(sizeof(force_connections.elasticforce_tether_osc_stiffness_2pt_data_2d)/sizeof(int)):
              forces.push_back( c_int_elasticforce_tether_osc_stiffness_2pt_2d[i] )

        elif (fc["ForceID"] == "ELASTIC_TETHER_OSCILLATINGSTIFFNESS" or fc["ForceID"] == "ELASTIC_TETHER_OSCILLATINGSTIFFNESS_TWOPOINT") and ibData["Dimensions"] == 3:
           c_elasticforce_tether_osc_stiffness_2pt_3d.ForceID = FORCE_CONNECTION_ELASTIC_TETHER_OSCILLATINGSTIFFNESS_TWOPOINT
           c_elasticforce_tether_osc_stiffness_2pt_3d.FCID = fc["FCID"]
           c_elasticforce_tether_osc_stiffness_2pt_3d.PointID = fc["PointID"]
           c_elasticforce_tether_osc_stiffness_2pt_3d.sigma = fc["sigma"]
           c_elasticforce_tether_osc_stiffness_2pt_3d.L = fc["L"]
           c_elasticforce_tether_osc_stiffness_2pt_3d.OscillatingTypeID = fc["OscillatingTypeID"]
           c_elasticforce_tether_osc_stiffness_2pt_3d.OscillatingFrequency = fc["OscillatingFrequency"]
           c_elasticforce_tether_osc_stiffness_2pt_3d.OscillatingAmplitude = fc["OscillatingAmplitude"]
           c_elasticforce_tether_osc_stiffness_2pt_3d.onoff_ratio = fc["onoff_ratio"]

           for i in range(1,ibData["Dimensions"]+1): c_elasticforce_tether_osc_stiffness_2pt_3d.PtInitLoc[i-1] = fc["PtInitLoc_"+str(i)]
           for i in range(1,ibData["Dimensions"]+1): c_elasticforce_tether_osc_stiffness_2pt_3d.PtVel[i-1] = fc["PtVel_"+str(i)]

           for i in range(sizeof(force_connections.elasticforce_tether_osc_stiffness_2pt_data_3d)/sizeof(int)):
              forces.push_back( c_int_elasticforce_tether_osc_stiffness_2pt_3d[i] )

        elif fc["ForceID"] == "PENALTY_MASS" and ibData["Dimensions"] == 2:
           c_penaltymass_data_2d.ForceID = FORCE_CONNECTION_PENALTYMASS
           c_penaltymass_data_2d.FCID = fc["FCID"]
           c_penaltymass_data_2d.PointID = fc["PointID"]
           c_penaltymass_data_2d.sigma = fc["sigma"]
           c_penaltymass_data_2d.mass = fc["mass"]
           c_penaltymass_data_2d.dv = fc["dv"]

           for i in range(1,ibData["Dimensions"]+1): c_penaltymass_data_2d.PtX[i-1] = ibData["Points"][PointIDMap[fc["PointID"]]]["X"+str(i)]
           for i in range(1,ibData["Dimensions"]+1): c_penaltymass_data_2d.PtVel[i-1] = 0.0
           for i in range(1,ibData["Dimensions"]+1): c_penaltymass_data_2d.PtVelOld[i-1] = 0.0

           for i in range(sizeof(force_connections.penaltymass_data_2d)/sizeof(int)):
              forces.push_back( c_int_penaltymass_data_2d[i] )

        elif fc["ForceID"] == "PENALTY_MASS" and ibData["Dimensions"] == 3:
           c_penaltymass_data_3d.ForceID = FORCE_CONNECTION_PENALTYMASS
           c_penaltymass_data_3d.FCID = fc["FCID"]
           c_penaltymass_data_3d.PointID = fc["PointID"]
           c_penaltymass_data_3d.sigma = fc["sigma"]
           c_penaltymass_data_3d.mass = fc["mass"]
           c_penaltymass_data_3d.dv = fc["dv"]

           for i in range(1,ibData["Dimensions"]+1): c_penaltymass_data_3d.PtX[i-1] = ibData["Points"][PointIDMap[fc["PointID"]]]["X"+str(i)]
           for i in range(1,ibData["Dimensions"]+1): c_penaltymass_data_3d.PtVel[i-1] = 0.0
           for i in range(1,ibData["Dimensions"]+1): c_penaltymass_data_3d.PtVelOld[i-1] = 0.0

           for i in range(sizeof(force_connections.penaltymass_data_3d)/sizeof(int)):
              forces.push_back( c_int_penaltymass_data_3d[i] )

        elif fc["ForceID"] == "FORCE_PULSE_TWOPOINT":
           c_force_pulse_2pt.ForceID = FORCE_CONNECTION_FORCE_PULSE_TWOPOINT
           c_force_pulse_2pt.FCID = fc["FCID"]
           c_force_pulse_2pt.PointID1 = fc["PointID1"]
           c_force_pulse_2pt.PointID2 = fc["PointID2"]
           c_force_pulse_2pt.sigma = fc["sigma"]
           c_force_pulse_2pt.OscillatingTypeID = fc["OscillatingTypeID"]
           c_force_pulse_2pt.OscillatingFrequency = fc["OscillatingFrequency"]
           c_force_pulse_2pt.OscillatingAmplitude = fc["OscillatingAmplitude"]
           c_force_pulse_2pt.onoff_ratio = fc["onoff_ratio"]

           if "time_offset" in fc: c_force_pulse_2pt.time_offset = fc["time_offset"]
           else: c_force_pulse_2pt.time_offset = 0.0

           for i in range(sizeof(force_connections.force_pulse_2pt_data)/sizeof(int)):
              forces.push_back( c_int_force_pulse_2pt[i] )
#Added by Saeed to handle the data structure for the BendingForceCurve class
        elif fc["ForceID"] == "BENDING_CURVE" and ibData["Dimensions"] == 2:
           c_bendingforcecurve_2d.ForceID = FORCE_CONNECTION_BENDING_CURVE
           c_bendingforcecurve_2d.FCID = fc["FCID"]
           c_bendingforcecurve_2d.L2PointID = fc["L2PointID"]
           c_bendingforcecurve_2d.LPointID = fc["LPointID"]
           c_bendingforcecurve_2d.R2PointID = fc["R2PointID"]
           c_bendingforcecurve_2d.RPointID = fc["RPointID"]
           c_bendingforcecurve_2d.PointID = fc["PointID"]
           c_bendingforcecurve_2d.C2 = fc["C2"]
           c_bendingforcecurve_2d.epsil = fc["epsil"]
           c_bendingforcecurve_2d.length = fc["length"]
           c_bendingforcecurve_2d.Omega = fc["Omega"]
           c_bendingforcecurve_2d.Kappa = fc["Kappa"]
           c_bendingforcecurve_2d.ds = fc["ds"]
           c_bendingforcecurve_2d.dv = fc["dv"]

           for i in range(sizeof(force_connections.bendingforcecurve_data_2d)/sizeof(int)):
              forces.push_back( c_int_bendingforcecurve_2d[i] )

        elif fc["ForceID"] == "BENDING_CURVE" and ibData["Dimensions"] == 3:
           c_bendingforcecurve_3d.ForceID = FORCE_CONNECTION_BENDING_CURVE
           c_bendingforcecurve_3d.FCID = fc["FCID"]
           c_bendingforcecurve_3d.L2PointID = fc["L2PointID"]
           c_bendingforcecurve_3d.LPointID = fc["LPointID"]
           c_bendingforcecurve_3d.R2PointID = fc["R2PointID"]
           c_bendingforcecurve_3d.RPointID = fc["RPointID"]
           c_bendingforcecurve_3d.PointID = fc["PointID"]
           c_bendingforcecurve_3d.C2 = fc["C2"]
           c_bendingforcecurve_3d.epsil = fc["epsil"]
           c_bendingforcecurve_3d.length = fc["length"]
           c_bendingforcecurve_3d.Omega = fc["Omega"]
           c_bendingforcecurve_3d.Kappa = fc["Kappa"]
           c_bendingforcecurve_3d.ds = fc["ds"]
           c_bendingforcecurve_3d.dv = fc["dv"]
           for i in range(sizeof(force_connections.bendingforcecurve_data_3d)/sizeof(int)):
              forces.push_back( c_int_bendingforcecurve_3d[i] )


        else:
           raise Exception("NEWHELPER: Invalid Immersed Boundary Type")

