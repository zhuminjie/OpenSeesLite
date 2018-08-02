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
#include <ElasticSection2d.h>
#include <LegendreBeamIntegration.h>
#include <ElementParameter.h>
#include <DispBeamColumn2d.h>
#include <DisplacementControl.h>
#include <UmfpackGenLinSOE.h>
#include <UmfpackGenLinSolver.h>
#include <Linear.h>
#include <PlainNumberer.h>
#include <ElasticPPMaterial.h>
#include <WideFlangeSectionIntegration.h>
#include <FiberSection2d.h>
#include <CorotCrdTransf2d.h>
#include <Steel01.h>
#include <vector>

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

    OPS_SetNDM(ndm);
    OPS_SetNDF(ndf);

    // variables
    double L = 100.0;
    double E = 29000.0;
    double A = 20.0;
    double I = 1400.0;
    double P = 30.0;
    double Fy = 36000.0;
    double hardeningRatio = 0.02;
    int numEles = 10;
    double hsize = L/numEles;

    // W 24 x 62
    double d = 23.7;
    double tw = 0.43;
    double bf = 7.04;
    double tf = 0.59;
    int Nfw = 10;
    int Nff = 10;

    // create nodes
    int numNodes = numEles+1;
    std::vector<Node*> nodes(numNodes);
    for (int i=0; i<numNodes; ++i) {
	nodes[i] = new Node(i+1, ndf, i*hsize, 0.0);
	domain->addNode(nodes[i]);
    }

    // fully fix node 1
    for (int i=0; i<ndf; ++i) {
	SP_Constraint* sp = new SP_Constraint(nodes[0]->getTag(), i, 0.0, true);
	domain->addSP_Constraint(sp);
    }

    // create a geometric transformation
    Vector jntOffsetI(2), jntOffsetJ(2);
    int transfTag = 1;
    //LinearCrdTransf2d* transf = new LinearCrdTransf2d(transfTag,jntOffsetI,jntOffsetJ);
    CorotCrdTransf2d* transf = new CorotCrdTransf2d(transfTag,jntOffsetI,jntOffsetJ);

    // create steel01 material
    int matTag = 1;
    Steel01* mat = new Steel01(matTag, Fy, E, hardeningRatio);

    // create a wide flange section
    WideFlangeSectionIntegration wfsect(d, tw, bf, tf, Nfw, Nff);
    int numFibers = wfsect.getNumFibers();
    UniaxialMaterial** mats = new UniaxialMaterial*[numFibers];
    wfsect.arrangeFibers(mats, mat);

    int secTag = 1;
    FiberSection2d* section = new FiberSection2d(secTag, numFibers, mats, wfsect);
    delete [] mats;

    // create an elastic section
    // int secTag = 1;
    // ElasticSection2d* section = new ElasticSection2d(secTag,E,A,I);
    

    // create a beam integration
    BeamIntegration* bi = new LegendreBeamIntegration;

    // create a displacement based beam element
    int eleTag = 1;
    int numSec = 2;
    double mass = 0.0;
    int cmass = 0;
    SectionForceDeformation* s[] = {section,section};

    std::vector<Element*> eles(numEles);
    for (int i=0; i<numEles; ++i) {
	eles[i] = new DispBeamColumn2d(i+1,i+1,i+2,numSec,s,*bi,*transf,mass,cmass);
	domain->addElement(eles[i]);
    }

    // create a parameter
    int paramTag = 1;
    const char* argvp[] = {"E"};
    Parameter* param = new ElementParameter(paramTag,eleTag,argvp,1);
    domain->addParameter(param);

    for (int i=1; i<numEles; ++i) {
	param->addComponent(eles[i], argvp, 1);
    }

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
    forces(1) = P;
    NodalLoad* load = new NodalLoad(loadTag, nodes.back()->getTag(), forces, false);
    domain->addNodalLoad(load,patternTag);

    // print domain information
    opserr << *domain;

    // create static analysis
    AnalysisModel* theAnalysisModel = new AnalysisModel;
    CTestNormUnbalance* theTest = new CTestNormUnbalance(1.0e-8,10,1);
    NewtonRaphson* theAlgorithm = new NewtonRaphson(*theTest);
    //Linear* theAlgorithm = new Linear;
    PlainHandler* theHandler = new PlainHandler;
    DOF_Numberer* theNumberer = new PlainNumberer;

    int dof = 0;
    double loadfactor = 44.35;
    int numsteps = 10;
    int numIter = 1;
    double incr = loadfactor/numsteps;
    LoadControl* theStaticIntegrator = new LoadControl(incr,numIter,incr,incr);


    // int dof = 1;
    // double incr = 0.1*100;
    // int numIter = 1;
    // int numsteps = 10;
    // DisplacementControl* theStaticIntegrator =
    // 	new DisplacementControl(nodes.back()->getTag()
    // 				,dof,incr,domain,numIter,
    // 				incr,incr);
    
    UmfpackGenLinSolver* theSolver = new UmfpackGenLinSolver;
    UmfpackGenLinSOE* theSOE = new UmfpackGenLinSOE(*theSolver);

    StaticAnalysis* theStaticAnalysis = new StaticAnalysis(*domain,
							   *theHandler,
							   *theNumberer,
							   *theAnalysisModel,
							   *theAlgorithm,
							   *theSOE,
							   *theStaticIntegrator,
							   theTest);

    

    // set sensitivity algorithm
    int analysisTypeTag = 1;
    theStaticIntegrator->setComputeType(analysisTypeTag);
    theStaticIntegrator->activateSensitivityKey();

    // perform analysis
    if (theStaticAnalysis->analyze(numsteps) < 0) {
	opserr << "WARNING: failed to do analysis\n";
	return -1;
    }

    // print node 2 displacements
    opserr<<"Node 2 Disp = "<<nodes.back()->getTrialDisp();
    int gradIndex = param->getGradIndex();
    dof = 3;
    double value = nodes.back()->getDispSensitivity(dof,gradIndex);
    opserr<<"Node 2 Disp sensitivity = "<<value<<"\n";

    double factor = pattern->getLoadFactor();
    opserr<<"Elastic Linear soln = "<<-P*factor*L*L/(2*E*E*I)<<"\n";

    opserr<<"Load factor = "<<factor<<"\n";

    // clean up
    theStaticAnalysis->clearAll();
    delete theStaticAnalysis;
    OPS_wipe();
    
    return 0;
}
