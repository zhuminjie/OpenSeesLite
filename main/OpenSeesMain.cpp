#include <iostream>
#include <OPS_Globals.h>
#include <elementAPI.h>
#include <StandardStream.h>
#include <UniaxialMaterial.h>
#include <NDMaterial.h>
#include <SectionForceDeformation.h>
#include <SectionRepres.h>
#include <TimeSeries.h>
#include <CrdTransf.h>
#include <BeamIntegration.h>
#include <NodalLoad.h>
#include <AnalysisModel.h>
#include <PlainHandler.h>
#include <RCM.h>
#include <AMDNumberer.h>
#include <LimitCurve.h>
#include <DamageModel.h>
#include <FrictionModel.h>
#include <HystereticBackbone.h>
#include <YieldSurface_BC.h>
#include <CyclicModel.h>
#include <FileStream.h>
#include <CTestNormUnbalance.h>
#include <NewtonRaphson.h>
#include <TransformationConstraintHandler.h>
#include <Newmark.h>
#include <ProfileSPDLinSolver.h>
#include <ProfileSPDLinDirectSolver.h>
#include <ProfileSPDLinSOE.h>
#include <SymBandEigenSolver.h>
#include <SymBandEigenSOE.h>
#include <FullGenEigenSolver.h>
#include <FullGenEigenSOE.h>
#include <ArpackSOE.h>
#include <LoadControl.h>
#include <CTestPFEM.h>
#include <PFEMIntegrator.h>
#include <TransientIntegrator.h>
#include <PFEMSolver.h>
#include <PFEMLinSOE.h>
#include <Accelerator.h>
#include <KrylovAccelerator.h>
#include <AcceleratedNewton.h>
#include <RaphsonAccelerator.h>
#include <SecantAccelerator2.h>
#include <PeriodicAccelerator.h>
#include <LineSearch.h>
#include <InitialInterpolatedLineSearch.h>
#include <BisectionLineSearch.h>
#include <SecantLineSearch.h>
#include <RegulaFalsiLineSearch.h>
#include <NewtonLineSearch.h>
#include <FileDatastore.h>
#include <UniaxialMaterial.h>
#include <Domain.h>
#include <StaticAnalysis.h>
#include <DirectIntegrationAnalysis.h>
#include <LinearSOE.h>
#include <EigenSOE.h>
#include <DOF_Numberer.h>
#include <ConstraintHandler.h>
#include <StaticIntegrator.h>
#include <TransientIntegrator.h>
#include <EquiSolnAlgo.h>
#include <FE_Datastore.h>
#include <FEM_ObjectBrokerAllClasses.h>
#include <PFEMAnalysis.h>
#include <VariableTimeStepDirectIntegrationAnalysis.h>
#ifdef _RELIABILITY
#include <ReliabilityStaticAnalysis.h>
#include <ReliabilityDirectIntegrationAnalysis.h>
#endif
#include <Timer.h>
#include <SimulationInformation.h>
#include <SP_Constraint.h>
#include <LinearCrdTransf2d.h>
#include <ElasticBeam2d.h>
#include <LinearSeries.h>
#include <LoadPattern.h>
#include <DOF_Numberer.h>

// Global variables
StandardStream sserr;
OPS_Stream *opserrPtr = &sserr;

Domain __theDomain;
int __ndm = 0;
int __ndf = 0;

// few APIs
Domain* OPS_GetDomain(void)
{
    return &__theDomain;
}

int OPS_GetNDF()
{
    return __ndf;
}

int OPS_GetNDM()
{
    return __ndm;
}

void OPS_SetNDF(int ndf)
{
    __ndf = ndf;
}

void OPS_SetNDM(int ndm)
{
    __ndm = ndm;
}

int OPS_wipe()
{
    Domain* domain = OPS_GetDomain();
    
    // wipe domain
    domain->clearAll();

    // time set to zero
    ops_Dt = 0.0;

    // wipe uniaxial material
    OPS_clearAllUniaxialMaterial();
    OPS_clearAllNDMaterial();

    // wipe sections
    OPS_clearAllSectionForceDeformation();
    OPS_clearAllSectionRepres();

    // wipe time series
    OPS_clearAllTimeSeries();

    // wipe GeomTransf
    OPS_clearAllCrdTransf();

    // wipe BeamIntegration
    OPS_clearAllBeamIntegrationRule();

    // wipe limit state curve
    OPS_clearAllLimitCurve();

    // wipe damages model
    OPS_clearAllDamageModel();

    // wipe friction model
    OPS_clearAllFrictionModel();

    // wipe HystereticBackbone
    OPS_clearAllHystereticBackbone();

    // wipe YieldSurface_BC
    OPS_clearAllYieldSurface_BC();

    // wipe CyclicModel
    OPS_clearAllCyclicModel();
    return 0;
}

int OPS_GetNumRemainingInputArgs()
{
    return 0;
}

int OPS_GetIntInput(int *numData, int*data)
{
    return 0;
}

int OPS_SetIntOutput(int *numData, int*data)
{
    return 0;
}

int OPS_GetDoubleInput(int *numData, double *data)
{
    return 0;
}

int OPS_SetDoubleOutput(int *numData, double *data)
{
    return 0;
}

const char * OPS_GetString(void)
{
    return 0;
}

int OPS_SetString(const char* str)
{
    return 0;
}

int OPS_ResetCurrentInputArg(int cArg)
{
    return 0;
}



// main function
int main(int argc, char** argv)
{
    Domain* domain = OPS_GetDomain();

    int ndm = 2;
    int ndf = 3;

    int numsteps = 1000;
    
    OPS_SetNDM(ndm);
    OPS_SetNDF(ndf);

    // create two end nodes
    int ndtag1 = 1, ndtag2 = 2;
    Node* node1 = new Node(ndtag1, ndf, 0.0, 0.0);
    Node* node2 = new Node(ndtag2, ndf, 1.0, 1.0);

    // add to domain
    domain->addNode(node1);
    domain->addNode(node2);

    // fully fix node 1
    for (int i=0; i<ndf; ++i) {
	SP_Constraint* sp = new SP_Constraint(ndtag1, i, 0.0, true);
	domain->addSP_Constraint(sp);
    }

    // create a geometric transformation
    Vector jntOffsetI(2), jntOffsetJ(2);
    int transfTag = 1;
    LinearCrdTransf2d* transf = new LinearCrdTransf2d(transfTag,jntOffsetI,jntOffsetJ);

    // create an elastic beam element
    double A = 1.0;
    double E = 1e8;
    double I = 1.0;
    int eleTag = 1;
    ElasticBeam2d* beam = new ElasticBeam2d(eleTag,A,E,I,ndtag1,ndtag2,*transf);
    domain->addElement(beam);

    // create a time series
    int tsTag = 1;
    TimeSeries* ts = new LinearSeries(tsTag, 1.0);

    // create a load pattern
    int patternTag = 1;
    LoadPattern* pattern = new LoadPattern(patternTag, 1.0);
    pattern->setTimeSeries(ts);
    domain->addLoadPattern(pattern);

    // create a nodal load at node 2
    int loadTag = 1;
    Vector forces(ndf);
    forces(1) = -1.0;
    NodalLoad* load = new NodalLoad(loadTag, ndtag2, forces, false);
    domain->addNodalLoad(load,patternTag);

    // print domain information
    opserr << *domain;

    // create static analysis
    AnalysisModel* theAnalysisModel = new AnalysisModel;
    CTestNormUnbalance* theTest = new CTestNormUnbalance(1.0e-6,25,0);
    NewtonRaphson* theAlgorithm = new NewtonRaphson(*theTest);
    PlainHandler* theHandler = new PlainHandler;
    RCM* theRCM = new RCM(false);
    DOF_Numberer* theNumberer = new DOF_Numberer(*theRCM);
    LoadControl* theStaticIntegrator = new LoadControl(1, 1, 1, 1);

    ProfileSPDLinDirectSolver* theSolver = new ProfileSPDLinDirectSolver;
    ProfileSPDLinSOE* theSOE = new ProfileSPDLinSOE(*theSolver);

    StaticAnalysis* theStaticAnalysis = new StaticAnalysis(*domain,
							   *theHandler,
							   *theNumberer,
							   *theAnalysisModel,
							   *theAlgorithm,
							   *theSOE,
							   *theStaticIntegrator,
							   theTest);

    // perform analysis
    if (theStaticAnalysis->analyze(numsteps) < 0) {
	opserr << "WARNING: failed to do analysis\n";
	return -1;
    }

    // print node 2 displacements
    const Vector& disp = node2->getTrialDisp();
    opserr<<"Node 2 Disp = "<<disp;

    // clean up
    theStaticAnalysis->clearAll();
    delete theStaticAnalysis;
    OPS_wipe();
    
    return 0;
}
