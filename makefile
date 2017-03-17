include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules
include ../../conf/ibrules

LOCDIR = src/ib/

SOURCES = ImmersedBoundary.cpp IForceConnection.cpp ElasticForce.cpp ElasticForce_Tether_2Pt.cpp ElasticForce_Tether_OscStiffness_2Pt.cpp ElasticForce_OscStiffness_2Pt.cpp ElasticForce_OscStiffness.cpp ElasticForce_OscLength_2Pt.cpp ElasticForce_OscLength.cpp ElasticForce_2Pt.cpp BendingForce.cpp PenaltyMass.cpp Force_Pulse_2Pt.cpp BendingForceCurve.cpp ElasticForceEnergy.cpp
OBJECTS=$(SOURCES:.cpp=.o)

CLEANFILES = $(OBJECTS) ${LIBDIR}/ib.a

compile_src: ${OBJECTS}
	ar rcs ${LIBDIR}/ib.a $(OBJECTS)

include ${PETSC_DIR}/conf/test

