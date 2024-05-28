#include "NKM_SurfPack.hpp"
#include "NKM_KrigingModel.hpp"
//#include "Accel.hpp"
//#include "NKM_LinearRegressionModel.hpp"
#include <math.h>
#include <iostream>
#include <cfloat>


#ifdef SURFPACK_HAVE_BOOST_SERIALIZATION
BOOST_CLASS_EXPORT(nkm::KrigingModel)
#endif


namespace nkm {

using std::cout;
using std::cerr;
using std::endl;
using std::ostringstream;


//#define __KRIG_ERR_CHECK__
#define __NKM_UNBIASED_LIKE__



// typical constructor
KrigingModel::KrigingModel(const SurfData& sd, const ParamMap& params)
  : SurfPackModel(sd,sd.getIOut()), numVarsr(sd.getNVarsr()),
    numTheta(numVarsr), numPoints(sdBuild.getNPts()), XR(sdBuild.xr)
{
  //printf("calling the right KrigingModel constructor\n"); fflush(stdout);

  //if the SurfDataScaler class does what it's supposed to (the only private content in sdBuild that a Model can access are the scaling bits and then only through SurfDataScaler, and only the model can see inside the scaler) the next line will cause an error when you try to compile with it uncommented, that is intentional
  //printf("scaler->mySd.iout=%d\n",scaler.mySd.iout);

  // OPTIONS PARSING
  ParamMap::const_iterator param_it;

  // *************************************************************
  // control verbosity outputLevel
  // *************************************************************
  param_it = params.find("verbosity");
  if (param_it != params.end() && param_it->second.size() > 0)
    outputLevel=static_cast<short>(std::atoi(param_it->second.c_str()));
  // outputLevel is a member of nkm::SurfPackModel which nkm::KrigingModel
  // is derived from

  // ********************************************************************
  // does the user want to use derivative information to build the
  // Kriging model (e.g. Gradient Enhanced Kriging)
  // ********************************************************************

  buildDerOrder=0; //default is regular Kriging (i.e. not GEK)
  param_it = params.find("derivative_order");
  if(param_it != params.end() && param_it->second.size() > 0)
    buildDerOrder = std::atoi(param_it->second.c_str());
  if(buildDerOrder==0) {//Kriging
    numEqnAvail=numPoints;
    nDer=1;
    Der.newSize(numVarsr,nDer); Der.zero();
  } else if(buildDerOrder==1) { //Gradient Enhanced Kriging (GEK)
    numEqnAvail=(1+numVarsr)*numPoints;
    multi_dim_poly_power(Der, numVarsr, 1);  //use all mixed partial
    //derivatives, up to first order, of the basis functions
    nDer=Der.getNCols(); //for GradKrigingModel nDer=(1+numVarsr);
    //printf("nDer=%d\n",nDer);
    int data_der_order=sdBuild.getDerOrder();
    if(data_der_order<1) {
      std::cerr << "the order of derivative information available in the "
		<< "build data is " << data_der_order << "\n"
		<< "You need to supply gradients of the output in order to "
		<< "construct a\nGradient Enhanced Kriging (GEK) Model."
		<< std::endl;
      assert(false);
    }
  }
  else{
    std::cerr << "derivative_order=" << buildDerOrder
	      << " in the nkm::KrigingModel constructor.\n"
	      << "For Kriging you must use derivative_order=0.\n"
	      << "For Gradient Enhanced Kriging (GEK) you must use "
	      << "derivative_order=1.\nHigher order derivative "
	      << "enhanced Kriging (e.g. Hessian Enhanced Kriging)\n"
	      << "has not been implemented." << std::endl;
    assert(false);
  }

  // *************************************************************
  // detect an anchor point if present this is the one point that
  // we make sure that the equationSelectingCholR does not discard
  // *************************************************************
  iAnchorPoint=0;
  ifHaveAnchorPoint=false;
  param_it = params.find("anchor_index");
  if (param_it != params.end() && param_it->second.size() > 0) {
    ifHaveAnchorPoint=true;
    //printf("nkm::KrigingModel() sees an anchor point\n");
    //fflush(stdout);
    iAnchorPoint=std::atoi(param_it->second.c_str());
    //printf("iAnchorPoint=%d\n",iAnchorPoint);
    //fflush(stdout);
    if(!((0<=iAnchorPoint)&&(iAnchorPoint<numPoints))) {
      std::cerr << "You can't specify an anchor point that isn't one of "
		<< "the build points" << std::endl;
      assert(false);
    }
  }

  // *************************************************************
  // this starts the input section about scaling the data
  // *************************************************************

  MtxDbl min_max_xr(numVarsr, 2);
  bool if_user_specified_lower_bounds=false;
  param_it = params.find("lower_bounds");
  if (param_it != params.end() && param_it->second.size() > 0) {
    if_user_specified_lower_bounds=true;
    if(min_max_xr.putCols(param_it->second,0)) {
      std::cerr << "You didn't enter the right number of lower bounds"
		<< std::endl;
      assert(false);
    }
  }

  bool if_user_specified_upper_bounds=false;
  param_it = params.find("upper_bounds");
  if (param_it != params.end() && param_it->second.size() > 0) {
    if_user_specified_upper_bounds=true;
    if(min_max_xr.putCols(param_it->second,1)) {
      std::cerr << "You didn't enter the right number of upper bounds"
		<< std::endl;
      assert(false);
    }
  }

  if(!(if_user_specified_lower_bounds==if_user_specified_upper_bounds)) {
    std::cerr << "Your options are to\n(A) specify both the upper and lower, or\n(B) specify neither the upper nor lower,\nbounds of the domain of the Kriging Model\n";
    assert(false);
  }

  if(if_user_specified_lower_bounds==true) {
    for(int ixr=0; ixr<numVarsr; ++ixr)
      if(!(min_max_xr(ixr,0)<=min_max_xr(ixr,1))) {
	std::cerr << "The lower bound of the domain of the Kriging Model must be less than or equal to the upper bound of the domain of the Kriging Model\n";
	assert(min_max_xr(ixr,0)<=min_max_xr(ixr,1));
      }
    //printf("lower_bounds = [%g",min_max_xr(0,0));
    //for(int ixr=1; ixr<numVarsr; ++ixr)
    //printf(", %g",min_max_xr(ixr,0));
    //printf("]^T, upper_bounds = [%g",min_max_xr(0,1));
    //for(int ixr=1; ixr<numVarsr; ++ixr)
    //printf(", %g",min_max_xr(ixr,1));
    //printf("]^T\n");
    sdBuild.setUnscaledDomainSize(min_max_xr);
  }

  //printf("KrigingModel constructor should have just written out domain bounds\n");

  param_it = params.find("dimension_groups");
  if (param_it != params.end() && param_it->second.size() > 0) {
    MtxInt dim_groups(numVarsr,1);
    if(dim_groups.putCols(param_it->second,0)) {
      std::cerr << "If you specify dimension_groups for any dimensions, "
		<< "you must specify groups\nfor all dimensions. If you "
		<< "don't want some of the dimensions to be grouped\n"
		<< "with other dimensions during scaling, give each of "
		<< "them their own group." << std::endl;
      assert(false);
    }
    sdBuild.setDimGroups(dim_groups);
  }

  scaler.scaleToDefault(); //scale outputs to -0.5<=Y<=0.5 and scale
  //real inputs to volume 1 hyper-rectangle centered at 0 if real
  //input dimensions are locked or the unit hypercube centered at 0 if
  //no dimensions are locked.  The scaling is done to let us define
  //the feasible region simply (region is defined in create);

  if(buildDerOrder==0) {
    sdBuild.getY(Yall);
    Yall.reshape(numPoints,1);
  }
  else if(buildDerOrder==1) {
    sdBuild.getUpToDerY(Yall,1);
    Yall.reshape(numEqnAvail,1);
    //Yall is now a column vector that contains
    //[y0, dy_0/dxr_0, ..., dy_0/dxr_{numVarsr-1}, y1, dy_1/dxr_0, ..., y2, ...
    // y_{numPoints-1}, dy_{numPoints-1}/dxr_0, ...,
    // dy_{numPoints-1}/dx_{numVarsr-1}]^T
  }

  // *************************************************************
  // this starts the input section about optimizing or directly
  // specifying correlation lengths, it must come after the
  // scaling section
  // *************************************************************

  // current options are none (fixed correl) | sampling (guess) | local |
  // global | global_local
  optimizationMethod = "global"; //the default
  //optimizationMethod = "none"; //the default
  param_it = params.find("optimization_method");
  if (param_it != params.end() && param_it->second.size() > 0)
    optimizationMethod = param_it->second;

  if(optimizationMethod.compare("none")==0)
    maxTrials=1;
  else if(optimizationMethod.compare("local")==0)
    maxTrials=20;
  else if(optimizationMethod.compare("sampling")==0)
    maxTrials=2*numVarsr+1;
  else if(optimizationMethod.compare("global")==0)
    maxTrials = 10000;
  else if(optimizationMethod.compare("global_local")==0) {
    maxTrials = 10000; //ensure it has non-zero as a fail safe but this
    //shouldn't be used
    maxTrialsGlobal = 500;
    maxTrialsLocal = 20;
  }
  else{ //error checking the input
    std::cerr << "KrigingModel() unknown optimization_method [" << optimizationMethod << "]  aborting\n";
    assert(false);
  }

  //std::cout << "optimization_method=\"" << optimizationMethod << "\"\n";

  //numStarts is the number of starting locations in a multi-start local search
  numStarts=1; //default is a single starting location
  param_it = params.find("num_starts");
  if (param_it != params.end() && param_it->second.size() > 0) {
    numStarts = std::atoi(param_it->second.c_str());
    if(numStarts<1) {
      std::cerr << "You can't specify fewer than one starting location "
		<< "for the optimization\nof correlation lenghts"
		<< std::endl;
      assert(false);
    }
  }

  if(!((numStarts==1)||(optimizationMethod.compare("local")==0))) {
    std::cerr << "Local optimization is the only optimization method for Kriging that uses the \"num_starts\" key word. Check your input file for errors.\n";
    assert(false);
  }

  //std::cout << "num_starts=" << numStarts << "\n";


  // does the user want to specify correlation lengths directly?
  ifUserSpecifiedCorrLengths=false; //the default is no
  param_it = params.find("correlation_lengths");
  if (param_it != params.end() && param_it->second.size() > 0) {
    ifUserSpecifiedCorrLengths=true;
    //printf("User specifying correlation lengths\n"); fflush(stdout);

    // make sure that the user didn't
    // * say they want to global optimize __AND__
    // * specify correlation lengths
    if(optimizationMethod.compare("global")==0) {
      std::cerr << "You can't both \n (A) use the global optimization method to choose, and \n (B) directly specify \n correlation lengths for the Kriging model.\n";
      assert(false);
    }
    else if(optimizationMethod.compare("global_local")==0) {
      //they can't coarse global followed by local either
      std::cerr << "You can't both \n (A) use the coarse global polished by local optimization method to choose, and \n (B) directly specify \n correlation lengths for the Kriging model.\n";
      assert(false);
    }
    else if(optimizationMethod.compare("sampling")==0) {
      // this is only the default number of samples/maxTrials; the user can
      // still overide this below
      maxTrials+=1;
    }

    natLogCorrLen.newSize(numVarsr,1); //allocate space

    //read the correlation lengths in from the string
    if(natLogCorrLen.putCols(param_it->second,0)) {
      std::cerr << "The specified correlation lengths had the wrong "
		<< "number of input dimensions." << std::endl;
      assert(false);
    }
    // "natLogCorrLen" currently holds the unscaled correlation LENGTHS, not
    // the natural log of the scaled correlation length, we need to fix that
    // but first we need to check the input for errors
    for(int ixr=0; ixr<numVarsr; ++ixr)
      if(!(natLogCorrLen(ixr,0)>0.0)) {
	std::cerr << "For the Kriging Model, correlation lengths must be strictly positive\n.";
	assert(false);
      }

    //printf("unscaled corr lens = [%12.6g",natLogCorrLen(0,0));
    //for(int ixr=1; ixr<numVarsr; ++ixr)
    //printf(", %12.6g",natLogCorrLen(ixr,0));
    //printf("]\n");

    scaler.scaleXrDist(natLogCorrLen); //scale the lengths
    //scaler.scaleXrOther(natLogCorrLen); //error
    //printf("scaled corr lens = [%12.6g",natLogCorrLen(0,0));
    //for(int ixr=1; ixr<numVarsr; ++ixr)
    // printf(", %12.6g",natLogCorrLen(ixr,0));
    //printf("]\n");
    //fflush(stdout);

    //compute the natural log of the correlation lengths
    for(int ixr=0; ixr<numVarsr; ++ixr)
      natLogCorrLen(ixr,0)=std::log(natLogCorrLen(ixr,0));

    natLogCorrLen.reshape(numVarsr,1);
    //natLogCorrLen will be the first of the initial iterates (guesses), this happens in the create() function below
  }
  //printf("If user specified correlationlengths we should have just printed them\n");

  // maximum objective evals for optimization or guess
  param_it = params.find("max_trials");
  if (param_it != params.end() && param_it->second.size() > 0) {
    maxTrials = std::atoi(param_it->second.c_str());
  }

  if(!(maxTrials > 0)) {
    std::cerr << "You can't specify a maximum number of trials that is "
	      << "less than or equal\nto zero." << std::endl;
    assert(false);
  }

  //printf("maxTrials=%d\n",maxTrials);


  // *************************************************************
  // this starts the input section about the trend function
  // *************************************************************
  polyOrderRequested = 2;
  ifReducedPoly=false;
  param_it = params.find("order");
  if (param_it != params.end() && param_it->second.size() > 0) {
    polyOrderRequested = std::atoi(param_it->second.c_str());
    //ssstd::cerr << "polyOrderRequested=" << polyOrderRequested << std::endl;
    if(!(polyOrderRequested >= 0)) {
      std::cerr << "You can't use a trend function with a polynomial "
		<< "order less than zero." << std::endl;
      assert(false);
    }
  }
  else{
    //if they don't specify a polynomial order use a main effects
    //polynomial with order 2 for the trend function, (if they do
    //specify a polynomial order assume they mean a full polynomial
    //order unless they specify that it's a reduced_polynomial)
    ifReducedPoly=true;
  }
  numTrend.newSize(polyOrderRequested+1,1);

  //cout << "order=" << polyOrder << "\n";

  //polyOrder = 2; //for debug
  //main_effects_poly_power(Poly, numVarsr, polyOrder); //for debug
  //commented out for debug

  param_it = params.find("reduced_polynomial");
  if (param_it != params.end() && param_it->second.size() > 0)
    if((std::atoi(param_it->second.c_str()))!=0)
      ifReducedPoly=true;

  //cout << "ifReducedPoly=" << ifReducedPoly << "\n";

  if(ifReducedPoly) {
    for(polyOrder=0; polyOrder<=polyOrderRequested; ++polyOrder)
      numTrend(polyOrder,0)=polyOrder*numVarsr+1;
    main_effects_poly_power(Poly, numVarsr, polyOrderRequested);
  }
  else{
    for(polyOrder=0; polyOrder<=polyOrderRequested; ++polyOrder)
      numTrend(polyOrder,0)=num_multi_dim_poly_coef(numVarsr, polyOrder);
    multi_dim_poly_power(Poly, numVarsr, polyOrderRequested);
  }


  // ********************************************************************
  // this starts the section about the choice of correlation functions
  // need to do build derivative order before this
  // ********************************************************************
  corrFunc=DEFAULT_CORR_FUNC;

  //POW_EXP_CORR_FUNC
  powExpCorrFuncPow=0.0; //only 1.0<=powExpCorrFunc<=2.0 are allowed
  //later if corrFunc==POW_EXP_CORR_FUNC and powExpCorrFuncPow==0.0 we know
  //we have an error
  param_it = params.find("powered_exponential");
  if(param_it != params.end() && param_it->second.size() > 0) {
    if(corrFunc!=DEFAULT_CORR_FUNC) {
      std::cerr << "You can only specify one correlation function\n";
      assert(false);
    }
    corrFunc=POW_EXP_CORR_FUNC;
    powExpCorrFuncPow=std::atof(param_it->second.c_str());
    if(!((1.0<=powExpCorrFuncPow)&&(powExpCorrFuncPow<=2.0))){
      std::cerr << "The powered exponential correlation function must have 1.0<=power<=2.0\n";
      assert(false);
    }
    //need to require 1<powExpCorrFuncPow if first derivatives are used
    //(otherwise no derivative is continuous at build points
    //will need to require powExpCorrFuncPow==2 of 2nd or higher order
    //derivatives are used
    if(powExpCorrFuncPow==1.0)
      corrFunc=EXP_CORR_FUNC;
    else if(powExpCorrFuncPow==2.0)
      corrFunc=GAUSSIAN_CORR_FUNC;
  }

  //MATERN_CORR_FUNC
  maternCorrFuncNu=0.0; //only 0.5, 1.5, 2.5, and infinity will be allowed
  //later if corrFunc==MATERN_CORR_FUNC and maternCorrFuncNu=0.0 we know
  //we have an error
  param_it = params.find("matern");
  if(param_it != params.end() && param_it->second.size() > 0) {
    if(corrFunc!=DEFAULT_CORR_FUNC) {
      std::cerr << "You can only specify one correlation function\n";
      assert(false);
    }
    if(param_it->second.compare("infinity")==0) {
      corrFunc=GAUSSIAN_CORR_FUNC;
      //matern nu=infinty is the Gaussian correlation function
    }
    else{
      corrFunc=MATERN_CORR_FUNC;
      maternCorrFuncNu=std::atof(param_it->second.c_str());
      if(!((maternCorrFuncNu==0.5)||(maternCorrFuncNu==1.5)||
	   (maternCorrFuncNu==2.5))) {
	//could allow more later if 3rd+ order derivatives are enabled later
	std::cerr << "For the Matern correlation function the only allowed values for nu are 0.5, 1.5, 2.5, and infinity\n";
	assert(false);
      }
      if(maternCorrFuncNu==0.5) {
	corrFunc=EXP_CORR_FUNC; //matern nu=0.5 is the exponential correlation function
	//need to disallow maternCorrFuncNu=0.5 if gradients or higher order derivatives are used to construct the Kriging model
      }
      //need to disallow maternCorrFuncNu=1.5 it hessians or higher order derivatives are used to construct the Kriging model
    }
  }

  if(corrFunc==DEFAULT_CORR_FUNC)
    corrFunc=GAUSSIAN_CORR_FUNC;

  // *************************************************************
  // this starts the input section HOW to bound the condition
  // number, this determines which derivatives of the constraint
  // function can be computed analytically so handle that here too
  // *************************************************************
  //constraintType="rcond"; //rcond is now the only option for type of
  //constraint against ill conditioning
  numConFunc=1;

  //convert to the Dakota bitflag convention for derivative orders
  int num_analytic_obj_ders_in=0; //analytical derivatives have been removed
  int num_analytic_con_ders_in=0; //analytical derivatives have been removed
  maxObjDerMode=(static_cast<int>(std::pow(2.0,num_analytic_obj_ders_in+1)))-1; //analytical gradients of objective function
  maxConDerMode=(static_cast<int> (std::pow(2.0,num_analytic_con_ders_in+1)))-1; //analytical gradients of constraint function(s)

  maxCondNum=std::pow(1024.0,4);

  // *************************************************************
  // this starts the input section about the nugget which can be
  // used to smooth the data and also decrease the condition
  // number
  // *************************************************************

  ifChooseNug = false;
  //ifChooseNug = true;
  ifAssumeRcondZero=false;
  param_it = params.find("find_nugget");
  if (param_it != params.end() && param_it->second.size() > 0) {
    ifChooseNug = true;
    int zero_or_one = std::atoi(param_it->second.c_str());
    if(zero_or_one==0)
      ifAssumeRcondZero=true;
  }
  //ifChooseNug = true ; std::cout << "ifChooseNug=" << ifChooseNug << "\n";

  // fixed value for now
  nug = 0.0; //default
  ifPrescribedNug=false;
  param_it = params.find("nugget");
  if (param_it != params.end() && param_it->second.size() > 0) {
    if(!(ifChooseNug==false)) {
      std::cerr << "You can do at most 1 of the following (A) auto-select "
		<< "the nugget\n(approximately the minimum needed to "
		<< "satisfy the condition number bound)\n(B) directly "
		<< "specify a nugget.  The default is not to use a nugget "
		<< "at all\n(i.e. use a nugget of zero)." << std::endl;
      assert(false);
    }
    nug = std::atof(param_it->second.c_str());
    if(!(nug >= 0.0)) {
      std::cerr << "The nugget must be greater than or equal to zero."
		<< std::endl;
      assert(false);
    }
    ifPrescribedNug=true;
  }

  // *************************************************************
  // this ends the input parsing now finish up the prep work
  // *************************************************************
  preAllocateMaxMemory(); //so we don't have to constantly dynamically
  //allocate, the SurfMat class can use a subset of the allocated memory
  //without using dynamic reallocation

  // precompute and store the trend function for all build points
  if(buildDerOrder==0) //for Kriging
    eval_trend_fn(Gall, XR);
  else if(buildDerOrder>=1) { //for GEK
    //actually this is generic to higher order derivative enhanced Kriging
    //(e.g. Hessian Enhanced Kriging) provided that Der is appropriately
    //defined
    eval_der_trend_fn(Gall, Der, XR);
  } else{
    std::cerr << "bogus buildDerOrder=" << buildDerOrder
	      << " in the constructor when evaluating Gall" << std::endl;
    assert(false);
  }

  if((ifChooseNug==true)||(ifPrescribedNug==true)) {
    //if we're using a nugget then we aren't using pivoted cholesky to
    //select an optimal subset of points, that means that the order of
    //points aren't going to change so we'll set Y and Gtran to what
    //we know they need to be
    iPtsKeep.newSize(numPoints,1);
    for(int ipt=0; ipt<numPoints; ++ipt)
      iPtsKeep(ipt,0)=ipt;

    Y.copy(Yall);

    //polyOrder=polyOrderRequested;
    nTrend=numTrend(polyOrderRequested,0);
    Gtran.newSize(numEqnAvail,nTrend);
    for(int itrend=0; itrend<nTrend; ++itrend)
      for(int i=0; i<numEqnAvail; ++i)
	Gtran(i,itrend)=Gall(itrend,i);

    if(buildDerOrder==0)
      numExtraDerKeep=0;
    else if(buildDerOrder==1)
      numExtraDerKeep=numVarsr;
    else{
      std::cerr << "buildDerOrder=" << buildDerOrder
		<< " in void KrigingModel::nuggetSelectingCholR(); "
		<< "for Kriging buildDerOrder must be 0; "
		<< "for Gradient Enhanced Kriging buildDerOrder must be 1; "
		<< "Higher order derivative enhanced Kriging "
		<< "(e.g Hessian Enhanced Kriging) has not been implemented"
		<< std::endl;
      assert(false);
    }
    numPointsKeep=numPoints;
    numRowsR=numEqnAvail;
  }

  gen_Z_matrix();  //initializes deltaXR and Z matrices (what the Z
  // matrix contains depends on the choose of correlation function but
  // as of 2012.06.25 all correlation functions involve an
  // coefficient*exp(Z^T*theta) (reshaped to the lower triangular part of R)
  // for the powered Exponential family of correlation functions that
  // coefficient is 1.0, for Matern 1.5 and Matern 2.5 correlation functions
  // it's not 1.0 and we only do the multiplcation when the coefficient isn't
  // 1.0 to save computation.

  //printf("completed the KrigingModel constructor\n"); fflush(stdout);
}

void KrigingModel::create()
{
  //printf("entered create()\n"); fflush(stdout);

  prevObjDerMode=prevConDerMode=0; //tells us not to reuse previous work used
  //to calculate the objective, constraints and their derivatives the first
  //time they are asked for
  prevTheta.newSize(numTheta,1);
  prevTheta.zero(); //not necessary just useful to debug

  //printf("KM.create():1: nug=%g\n",nug);

  // -
  // solve the optimization problem for the correlations
  // -

  //printf("numVarsr=%d\n",numVarsr); fflush(stdout);
  OptimizationProblem opt(*this, numVarsr, numConFunc);


  // set the bounds for the plausible region for correlation lengths
  // (assumes input space has a volume of 1, and data points are
  // uniformly distributed)
  aveDistBetweenPts=std::pow(numPoints,-1.0/numVarsr);

  // note: we should explore different bounds on the correlation lengths
  // for different choices of the correlation function, but this has not
  // been done yet, the "procedure" for determining how far the lengths
  // should extend is to say that X% probability mass is located a
  // certain distance away (5% at 16 neighbors for the upper bound, and
  // for the lower bound you want the nearest neighbor to be
  // "essentially uncorrelated" while halfway between nearest neigbors
  // is slightly correlated)

  // For the maximum correlation length = aveDistBetweenPts*8.0
  // the Gaussian Correlation function (the original one this GP a.k.a.
  // Kriging model was developed for) has about ~5% confidence (2 std
  // devs away) in what points 16 neighbors away have to say. If points
  // are correlated well at even greater distances then either
  // * that same information will be contained in nearby points OR
  // * you shouldn't be using a Gaussian process error model
  double max_corr_length = aveDistBetweenPts*8.0;
  maxNatLogCorrLen=std::log(max_corr_length);

  // For the minimum correlation length = aveDistBetweenPts/4.0
  // the Gaussian Correlation function (the original one this GP a.k.a.
  // Kriging model was developed for) has about ~5% confidence (2 std
  // midway between neighboring points... i.e. you're 4 std devs away
  // from your nearest neighbor so all sample points are treated as being
  // essentially uncorrelated
  double min_corr_length = aveDistBetweenPts/4.0;
  minNatLogCorrLen=std::log(min_corr_length);

  //Choose dead center (in log(correlation length)) of the feasible region
  //as the default initial guess for the Gaussian Process error model
  double init_guess=0.5*(maxNatLogCorrLen+minNatLogCorrLen);

  ///set the bounds and the initial iterates
  if(ifUserSpecifiedCorrLengths==true) {
    // the first guess is what the user told us he/she wanted to use
    for(int ixr=0; ixr<numVarsr; ++ixr) {
      opt.lower_bound(ixr, minNatLogCorrLen);
      opt.upper_bound(ixr, maxNatLogCorrLen);
      opt.initial_iterate(ixr, natLogCorrLen(ixr,0));
    }
    // the second guess is the center of the small feasible region
    MtxDbl the_second_guess(numVarsr,1);
    for(int ixr=0; ixr<numVarsr; ++ixr)
      the_second_guess(ixr,0)=init_guess;
    opt.add_initial_iterates(the_second_guess);
  } else {

    // since the user didn't specify an initial guess we will use the center
    // of the small feasible region as our first initial guess
    for(int ixr=0; ixr< numVarsr; ++ixr) {
      opt.lower_bound(ixr, minNatLogCorrLen);
      opt.upper_bound(ixr, maxNatLogCorrLen);
      opt.initial_iterate(ixr, init_guess);
    }
  }

  // add a binning optimal (read as space filling) random design with
  // 2*numVars more guesses
  // bins are the endpoints of a randomly rotated axis
  MtxDbl axes_of_guesses(numVarsr,2*numVarsr);
  gen_rand_axis_bin_opt_samples_0to1(axes_of_guesses,numVarsr);
  for(int iguess=0; iguess<2*numVarsr; ++iguess) {
    for(int ixr=0; ixr<numVarsr; ++ixr) {
      axes_of_guesses(ixr,iguess)=(maxNatLogCorrLen-minNatLogCorrLen)*axes_of_guesses(ixr,iguess)+minNatLogCorrLen;
    }
  }
  opt.add_initial_iterates(axes_of_guesses);

  //choose the optimizer you want to use
  if(optimizationMethod.compare("none")==0) {
    natLogCorrLen.resize(numVarsr,1);
    opt.retrieve_initial_iterate(0,natLogCorrLen);
  }
  else{
    if(optimizationMethod.compare("local")==0) {
      //local optimization
      if(numStarts==1)
	//from a single starting location
	opt.conmin_optimize();
      else{
	//doing multi-start local optimization
	opt.multistart_conmin_optimize(numStarts);
      }
    }
    else if(optimizationMethod.compare("global")==0)
      //global optimization via the "DIvision of RECTangles" method
      opt.direct_optimize();
    else if(optimizationMethod.compare("sampling")==0)
      //randomly generate candidates and pick the best guess
      opt.best_guess_optimize(maxTrials);
    else if(optimizationMethod.compare("global_local")==0){
      //a coarse global optimization that is polished by
      //local optimization
      maxTrials=maxTrialsGlobal;
      opt.direct_optimize();
      natLogCorrLen = opt.best_point();
      maxTrials=maxTrialsLocal;
      opt.conmin_optimize();
    }
    else{
      std::cerr << "KrigingModel:create() unknown optimization_method [" << optimizationMethod << "]\naborting" << std::endl;
      assert(false);
    }
    natLogCorrLen = opt.best_point();
  }

  MtxDbl corr_len(numVarsr,1);
  for(int ixr=0; ixr<numVarsr; ++ixr)
    corr_len(ixr,0)=std::exp(natLogCorrLen(ixr,0));
  correlations.newSize(numVarsr,1);
  get_theta_from_corr_len(correlations,corr_len);

  //printf("}\n");

  //printf("scaled correlations=[%12.6g",correlations(0,0));
  //for(int ixr=1; ixr<numVarsr; ++ixr)
  //printf(", %12.6g",correlations(ixr,0));
  //printf("]\n");

  masterObjectiveAndConstraints(correlations, 1, 0);

  //keep only the "optimal" subset of trend basis function in Poly that was
  //selected by the pivoted Cholesky factorization of G*R^-1*G^
  if(nTrend<numTrend(polyOrderRequested,0)) {
    //for this to work, the basis function indicices in iTrendKeep must
    //be in monotonically increasing order

    //we are guaranteed to keep the constant term of the trend function so
    //start loop from 1 not zero
    for(int itrend=1; itrend<nTrend; ++itrend) {
      int isrc=iTrendKeep(itrend,0);
      if(itrend<isrc)
	for(int ixr=0; ixr<numVarsr; ++ixr)
	   Poly(ixr,itrend)=Poly(ixr,isrc);
    }

    //now reduce the size of Poly
    Poly.resize(numVarsr,nTrend);
  }

  //determine the maximum total order of any term in the part of the
  //trend that was retained
  polyOrder=Poly(0,nTrend-1);
  for(int ixr=1; ixr<numVarsr; ++ixr)
    polyOrder+=Poly(ixr,nTrend-1);



  //make a reordered copy of (the retained portion of) XR for evaluation speed
  XRreorder.newSize(numVarsr,numPointsKeep);
  for(int ipt=0; ipt<numPointsKeep; ++ipt) {
    int isrc=iPtsKeep(ipt,0);
    for(int ixr=0; ixr<numVarsr; ++ixr)
      XRreorder(ixr,ipt)=XR(ixr,isrc);
  }

  if(outputLevel >= NORMAL_OUTPUT) {
    std::cout << model_summary_string();
    //std::cout << std::endl;
  }


  //variables whose values needed to be retained between sequential call to masterObjectiveAndConstraints for precompute and store strategy to work
  prevObjDerMode=prevConDerMode=0;

  //deallocate matrices we no longer need after emulator has been created
  //these were made member variables (instead of local variables) to avoid
  //the cost of dynamic allocation and deallocation each cycle of the
  //optimization of the correlation parameters
  scaleRChol.clear(); //matrix
  sumAbsColR.clear(); //vector
  oneNormR.clear(); //vector
  lapackRcondR.clear(); //vector
  rcondDblWork.clear();  //vector
  rcondIntWork.clear(); //vector
  Yall.clear(); //vector
  Gall.clear(); //matrix
  Gtran.clear(); //matrix
  iTrendKeep.clear(); //vector
  Z.clear(); //matrix
  Ztran_theta.clear(); //vector
  deltaXR.clear(); //matrix
  R.clear(); //matrix
  G_Rinv_Gtran.clear(); //matrix
  G_Rinv_Gtran_Chol_Scale.clear(); //vector
  G_Rinv_Gtran_Chol_DblWork.clear(); //vector
  G_Rinv_Gtran_Chol_IntWork.clear(); //vector
  G_Rinv_Y.clear(); //vector
  eps.clear(); //vector
  prevTheta.clear(); //vector
  con.clear(); //vector
}

std::string KrigingModel::get_corr_func() const {
  std::ostringstream oss;

  switch(corrFunc) {
  case GAUSSIAN_CORR_FUNC:
    oss << "Gaussian";
    break;
  case EXP_CORR_FUNC:
    oss << "exponential";
    break;
  case POW_EXP_CORR_FUNC:
    oss << "powered exponential with power=" << powExpCorrFuncPow;
    break;
  case MATERN_CORR_FUNC:
    oss << "Matern " << static_cast<int>(maternCorrFuncNu*2.0) << "/2";
    break;
  default:
    std::cerr << "unknown correlation function enumerated as " << corrFunc
	      << std::endl;
    assert(false);
  }
  return (oss.str());
}


std::string KrigingModel::model_summary_string() const {
  MtxDbl temp_out_corr_lengths(numVarsr,1);
  get_corr_len_from_theta(temp_out_corr_lengths,correlations);
  scaler.unScaleXrDist(temp_out_corr_lengths);

  //printf("numPoints=%d numTrend=%d numPointsKeep=%d numWholePointsKeep=%d numExtraDerKeep=%d\n",numPoints,numTrend(polyOrder,0),numPointsKeep,numWholePointsKeep,numExtraDerKeep);

  std::ostringstream oss;
  oss << "--- Surfpack Kriging Diagnostics ---\n";
  if(buildDerOrder==0)
    oss << "KM: #real inputs=" << numVarsr << "; #pts=" << numPoints
	<< "; used " << numPointsKeep << "/" << numPoints << " pts;\n";
  else if(buildDerOrder==1)
    oss << "GEK: #real inputs=" << numVarsr << "; #pts=" << numPoints
	<< "; #eqns=" << numEqnAvail << "; used "
	<< numRowsR << "/" << numEqnAvail << " eqns;\n";
  else{
    oss << "error std::string KrigingModel::model_summary_string() const\n"
	<< "buildDerOrder=" << buildDerOrder << "; it should be 0 for Kriging"
	<< " or 1 for Gradient Enhanced Kriging (GEK);"
	<< " the model_summary_string() function will need to be modified "
	<< "to handle other build derivative orders.\n";
  }
  oss << "using the ";
  if(corrFunc==GAUSSIAN_CORR_FUNC)
    oss << "Gaussian";
  else if(corrFunc==EXP_CORR_FUNC)
    oss << "exponential";
  else if(corrFunc==POW_EXP_CORR_FUNC)
    oss << "powered exponential (with power = " << powExpCorrFuncPow << ")";
  else if(corrFunc==MATERN_CORR_FUNC)
    oss << "Matern " << maternCorrFuncNu;
  else{
    std::cerr << "unknown corr func in model_summary_string()" << std::endl;
    assert(false);
  }
  oss << " correlation function with (unscaled)\n"
      << "Correlation lengths=[" << temp_out_corr_lengths(0,0);
  for(int ixr=1; ixr<numVarsr; ++ixr)
    oss << ", " << temp_out_corr_lengths(ixr,0);
  oss << "]^T\nfound by the \"" << optimizationMethod
      << "\" optimization_method;\nunadjusted variance="
      << estVarianceMLE * scaler.unScaleFactorVarY()
      << "; \"per equation\" log(likelihood)=" << likelihood << ";\n"
      << "rcond(R)=" << rcondR << "; rcond(G_Rinv_Gtran)="
      << rcond_G_Rinv_Gtran << "; [if either rcond is less\n"
      << "than 2^-40 (approx 9.095*10^-13) then the matrix is ill-conditioned "
      << "and\nthat \"voids the warranty\" of the Kriging Model]; nugget="
      << nug << ".  A ";
  if(polyOrder>1) {
    if(ifReducedPoly==true)
      oss << "reduced_";
    else oss <<"full ";
  }
  oss << "polynomial\nof order " << polyOrderRequested << " (with "
      << numTrend(polyOrderRequested,0) << " terms) was requested "
      << "for the trend function; the build\ndata was ";
  if(nTrend<numTrend(polyOrderRequested,0) )
    oss << "NOT ";
  oss << "sufficient to use the requested trend function; "
      << "the highest total\npolynomial order of any term in the "
      << "utlized trend function is " << polyOrder << ";\n"
      << "for SCALED inputs and outputs the utilized trend function is\n"
      << "betaHat^T*g(x)=";
  int nterm_on_this_line=0;
  for(int itrend=0; itrend<nTrend; ++itrend) {
    ++nterm_on_this_line;
    oss << betaHat(itrend,0);
    for(int ixr=0; ixr<numVarsr; ++ixr)
      if(Poly(ixr,itrend)>0) {
	oss << "*x" << ixr;
	if(Poly(ixr,itrend)>1)
	  oss << "^" << Poly(ixr,itrend);
      }
    if(itrend<nTrend-1) {
      oss << " ";
      if(betaHat(itrend+1,0)>=0.0)
	oss << "+ ";
      if(nterm_on_this_line==3) {
	oss << "...\n               ";
	nterm_on_this_line=0;
      }
    }
  }
  oss << "\n------------------------------------\n";
  return (oss.str());
}

void KrigingModel::preAllocateMaxMemory() {
  //this preallocates the maximum sizce of arrays whose size depends on how many equations were discarded by pivoted Cholesky and they could possibly be allocated to a different size than their maximum the first time they are allocated.

  nTrend=numTrend(polyOrderRequested,0);
  Y.newSize(numEqnAvail,1);
  Gtran.newSize(numEqnAvail,nTrend);
  Rinv_Gtran.newSize(numEqnAvail,nTrend);
  G_Rinv_Gtran.newSize(nTrend,nTrend);
  G_Rinv_Gtran_Chol.newSize(nTrend,nTrend);
  rhs.newSize(numEqnAvail,1);
  betaHat.newSize(nTrend,1);
  G_Rinv_Y.newSize(nTrend,1);
  eps.newSize(numEqnAvail,1);
  iPtsKeep.newSize(numPoints,1);
  RChol.newSize(numEqnAvail,numEqnAvail);
  int nrows = nTrend;
  if((ifChooseNug==false)&&(ifPrescribedNug==false)&&(numEqnAvail>nTrend))
    nrows=numEqnAvail;
  scaleRChol.newSize(nrows,3);
  lapackRcondR.newSize(nrows,1);

  return;
}

// BMA TODO: combine these two functions?

/// evaluate (y) the Kriging Model at a single point (xr)
double KrigingModel::evaluate(const MtxDbl& xr)
{
  if(buildDerOrder==0) {
    //you wouldn't want to do this for Gradient Enhanced Kriging
    //i.e. if gradients of y were used as inputs
    double singular_y;
    if(scaler.isYSingular(0,singular_y))
      return singular_y;
  }

  //assert( (numVarsr == xr.getNCols()) && (xr.getNRows() == 1) );
  MtxDbl g(nTrend,1), r(numRowsR,1);

  /*
  printf("double evaluate()\n");
  printf("xr=[%20.14g", xr(0,0));
  for(int ixr=1; ixr<numVarsr; ++ixr)
    printf(", %20.14g",xr(ixr,0));
    printf("]^T\n");
  */

  if(scaler.isUnScaled()) {
    eval_trend_fn(g, xr);
    correlation_matrix(r, xr);
  }
  else{
    MtxDbl xr_scaled(xr);
    scaler.scaleXrOther(xr_scaled);
    eval_trend_fn(g, xr_scaled);
    correlation_matrix(r, xr_scaled);
  }

  double y = dot_product(g, betaHat) + dot_product(r, rhs);

  //double yus=scaler.unScaleYOther(y);
  //printf("y=%g yunscaled=%g\n",y,yus);
  //return yus;

  return (scaler.unScaleYOther(y));
}


/// evaluate (y) the Kriging Model at a collection of points (xr)
MtxDbl& KrigingModel::evaluate(MtxDbl& y, const MtxDbl& xr)
{
  int nptsxr=xr.getNCols();
  //printf("nptsxr=%d nvarsrxr=%d",nptsxr,xr.getNCols());

  y.newSize(1,nptsxr);
  if(buildDerOrder==0) {
    //you wouldn't want to do this for Gradient Enhanced Kriging
    //i.e. if gradients of y were used as inputs
    double singular_y;
    if(scaler.isYSingular(0,singular_y)) {
      for(int ipt=0; ipt<nptsxr; ++ipt)
	y(0,ipt)=singular_y;
      return y;
    }
  }
  //assert(numVarsr == xr.getNRows());
  MtxDbl g(nTrend, nptsxr), r(numRowsR, nptsxr);

  if(scaler.isUnScaled()) {
    eval_trend_fn(g, xr);
    correlation_matrix(r, xr);
  }
  else{
    MtxDbl xr_scaled(xr);
    scaler.scaleXrOther(xr_scaled);
    eval_trend_fn(g, xr_scaled);
    correlation_matrix(r, xr_scaled);
  }

  //y=0.0*y+1.0*betaHat^T*g => y = betaHat^T*g
  matrix_mult(y, betaHat, g, 0.0, 1.0,'T','N');

  //y=1.0*y+1.0*r*rhs where rhs=R^-1*(Y-G(XR)^T*betaHat), initial y=betaHat^T*g => y=betaHat^T*g+rhs^T*r
  matrix_mult(y, rhs    , r, 1.0, 1.0,'T','N');

  scaler.unScaleYOther(y);

  //printf("y is correct for ValidateMain because it isn't being unscaled\n");

  return y;
}

MtxDbl& KrigingModel::evaluate_d1y(MtxDbl& d1y, const MtxDbl& xr)
{
  int nptsxr=xr.getNCols();
#ifdef __KRIG_ERR_CHECK__
  assert((numVarsr == xr.getNRows())&&(0<nptsxr));
#endif
  d1y.newSize(numVarsr, nptsxr);
  if(buildDerOrder==0) {
    //you wouldn't want to do this for Gradient Enhanced Kriging
    //i.e. if gradients of y were used as inputs
    double singular_y;
    if(scaler.isYSingular(0,singular_y)) {
      d1y.zero();
      return d1y;
    }
  }

  /*
  printf("evaluate_d1y()\n");
  for(int ipt=0; ipt<numPoints; ++ipt) {
    printf("XR(:,%3d)=[%12.6g",ipt,XR(0,ipt));
    for(int ixr=1; ixr<numVarsr; ++ixr)
      printf(", %12.6g",XR(ixr,ipt));
    printf("]^T Y(%3d)=%12.6g\n",ipt,Y(0,ipt));
  }
  */

  MtxDbl xr_scaled(xr);
  if(~(scaler.isUnScaled())) {
    //printf("scaling xr_scaled\n");
    scaler.scaleXrOther(xr_scaled);
  }

  /*
  printf("xr       =[%12.6g, %12.6g]\n",xr(0,0),xr(1,0));
  printf("xr_scaled=[%12.6g, %12.6g]\n",xr_scaled(0,0),xr_scaled(1,0));
  */

  int nder=num_multi_dim_poly_coef(numVarsr,-1);
  MtxInt der(numVarsr,nder);
  multi_dim_poly_power(der,numVarsr,-1); //equivalent to der.identity();

  evaluate_poly_der(d1y,flyPoly,derivBetaHat,Poly,der,betaHat,xr_scaled);

  MtxDbl r(numRowsR,nptsxr);
  correlation_matrix(r, xr_scaled);
  //apply_nugget_eval(r);
  MtxDbl d1r(numRowsR,nptsxr);
  MtxDbl temp_vec(1,nptsxr);

  for(int ider=0; ider<nder; ++ider) {

    //find the single dimension we are taking the first derviative of
    int ixr;
    for(ixr=0; ixr<numVarsr; ++ixr)
      if(der(ixr,ider)>0)
	break;
    //printf("ixr=%d ",ixr);
#ifdef __KRIG_ERR_CHECK__
    assert(ixr==ider);
#endif

    double d1y_unscale_factor=scaler.unScaleFactorDerY(ixr);
    //printf("d1y_usf=%g\n",d1y_unscale_factor);

    dcorrelation_matrix_dxI(d1r, r, xr_scaled, ixr);
    matrix_mult(temp_vec,rhs,d1r,0.0,1.0,'T');

    for(int ipt=0; ipt<nptsxr; ++ipt)
      d1y(ider,ipt)=(d1y(ider,ipt)+temp_vec(0,ipt))*d1y_unscale_factor;
  }
  /*
  printf("d1y(:,0)=[%g",d1y(0,0));
  for(int ider=1; ider<numVarsr; ++ider)
    printf(", %g",d1y(ider,0));
  printf("]\n");
  */
  return d1y;
}

MtxDbl& KrigingModel::evaluate_d2y(MtxDbl& d2y, const MtxDbl& xr)
{
  int nptsxr=xr.getNCols();
  int nder=num_multi_dim_poly_coef(numVarsr,-2);
  d2y.newSize(nder,nptsxr);
  if(buildDerOrder==0) {
    double singular_y;
    if(scaler.isYSingular(0,singular_y)) {
      //you wouldn't want to do this for gradient based Kriging
      //if gradients of y were used as inputs
      d2y.zero();
      return d2y;
    }
  }

  MtxDbl xr_scaled(xr);
  if(~(scaler.isUnScaled()))
    scaler.scaleXrOther(xr_scaled);
  //assert(numVarsr == xr.getNCols());

  MtxInt der(numVarsr,nder);
  MtxInt thisder(numVarsr,1);
  multi_dim_poly_power(der,numVarsr,-2);

  evaluate_poly_der(d2y,flyPoly,derivBetaHat,Poly,der,betaHat,xr_scaled);

  MtxDbl r(numRowsR,nptsxr);
  correlation_matrix(r, xr);
  //apply_nugget_eval(r);
  MtxDbl d1r(numRowsR,nptsxr);
  MtxDbl d2r(numRowsR,nptsxr);
  MtxDbl temp_vec(1,nptsxr);

  for(int ider=0; ider<nder; ++ider) {
    int ixr, jxr, ixrold=-1;

    der.getCols(thisder,ider);
    double d2y_unscale_factor=scaler.unScaleFactorDerY(thisder);
    //std::cout << "thisder=[" << thisder(0,0) << ", " << thisder(1,0)
    //<< "]^T unscalefactor=" << d2y_unscale_factor << std::endl;

    //find the first dimension we are taking a first derviative of
    for(ixr=0; ixr<numVarsr; ++ixr)
      if(der(ixr,ider)>0)
	break;

    if(ixr!=ixrold) {
      ixrold=ixr;
      dcorrelation_matrix_dxI(d1r, r, xr_scaled, ixr);
    }

    //find the second dimension we are taking a first derivative of
    if(der(ixr,ider)==2)
      jxr=ixr;
    else
      for(jxr=ixr+1; jxr<numVarsr; ++jxr)
	if(der(jxr,ider)>0)
	  break;
#ifdef __KRIG_ERR_CHECK__
    assert(jxr<numVarsr);
#endif

    //dcorrelation_matrix_dxI(d2r, d1r, xr_scaled, jvar);
    d2correlation_matrix_dxIdxJ(d2r, d1r, r, xr_scaled, ixr, jxr);
    //std::cout << "ider=" << ider << " size(d2r)=[" << d2r.getNRows()
    //	      << ", " << d2r.getNCols() << "]" << std::endl;
    matrix_mult(temp_vec,rhs,d2r,0.0,1.0,'T','N');

    for(int ipt=0; ipt<nptsxr; ++ipt)
      d2y(ider,ipt)=(d2y(ider,ipt)+temp_vec(0,ipt))*d2y_unscale_factor;
  }

  return d2y;
}



/** matrix Ops evaluation of adjusted variance at a single point
    adj_var=unadjvar*
            (1-r^T*R^-1*r+(g-G*R^-1*r)^T*(G*R^-1*G^T)^-1*(g-G*R^-1*r))
    on a point by point basis */
double KrigingModel::eval_variance(const MtxDbl& xr)
{
#ifdef __KRIG_ERR_CHECK__
  assert( (numVarsr==xr.getNRows()) && (xr.getNCols()==1) );
#endif
  MtxDbl g_minus_G_Rinv_r(nTrend,1), r(numRowsR,1);

  double unscaled_unadj_var=estVarianceMLE;
  if(scaler.isUnScaled()) {
    eval_trend_fn(g_minus_G_Rinv_r, xr);
    correlation_matrix(r, xr);
  }
  else{
    unscaled_unadj_var*=scaler.unScaleFactorVarY();
    MtxDbl xr_scaled(xr);
    scaler.scaleXrOther(xr_scaled);
    eval_trend_fn(g_minus_G_Rinv_r, xr_scaled);
    correlation_matrix(r, xr_scaled);
  }
  //at this point g_minus_G_Rinv_r holds g

  MtxDbl Rinv_r(numRowsR,1);
  MtxDbl G_Rinv_Gtran_inv_g_minus_G_Rinv_r(nTrend,1);


  solve_after_Chol_fact(Rinv_r,RChol,r);

  matrix_mult(g_minus_G_Rinv_r,Rinv_Gtran,r,1.0,-1.0,'T','N');
  //at this point g_minus_G_Rinv_r holds g-G*R^-*r (i.e. the name is correct)

  solve_after_Chol_fact(G_Rinv_Gtran_inv_g_minus_G_Rinv_r,
			G_Rinv_Gtran_Chol,g_minus_G_Rinv_r);


  double adj_var=unscaled_unadj_var*
    (1.0-dot_product(Rinv_r,r)+
     dot_product(G_Rinv_Gtran_inv_g_minus_G_Rinv_r,g_minus_G_Rinv_r));
  //if(!(adj_var>0.0)) {
  //printf("adj_var=%g unscaled_unadj_var=%g rcondR=%g\n",adj_var,unscaled_unadj_var,rcondR);
  //fflush(stdout);
  //}
  adj_var=std::fabs(adj_var); //hack to handle "negative zero" variance (numerical precision round off error)
  if(adj_var<0.0) {
    printf("NKM setting adj_var to zero adj_var=%g unadj_var=%g rcondR=%g\n",adj_var,unscaled_unadj_var,rcondR);
    adj_var=0.0;
  }
  else if(adj_var==0.0)
    printf("NKM adj_var is zero =%g\n",adj_var);
  else if(!(adj_var>=0.0))
    printf("double NKM_KrigingModel::eval_variance(...) adj_var=nan rcondR=%g\n",rcondR);

  return adj_var;
}

/** Evaluate the adjusted variance for a collection of points using matrix
    ops (i.e. BLAS and LAPACK) as much as possible)
    adj_var=unadjvar*
            (1-r^T*R^-1*r+(g-G*R^-1*r)^T*(G*R^-1*G^T)^-1*(g-G*R^-1*r))
    on a point by point basis */
MtxDbl& KrigingModel:: eval_variance(MtxDbl& adj_var, const MtxDbl& xr)
{
#ifdef __KRIG_ERR_CHECK__
  assert(numVarsr==xr.getNRows());
#endif
  int nptsxr=xr.getNCols();
  adj_var.newSize(1,nptsxr);
  MtxDbl g_minus_G_Rinv_r(nTrend,nptsxr), r(numRowsR,nptsxr);

  double unscaled_unadj_var=estVarianceMLE;
  if(scaler.isUnScaled()) {
    eval_trend_fn(g_minus_G_Rinv_r, xr);
    correlation_matrix(r, xr);
  }
  else{
    unscaled_unadj_var*=scaler.unScaleFactorVarY();
    MtxDbl xr_scaled(xr);
    scaler.scaleXrOther(xr_scaled);
    eval_trend_fn(g_minus_G_Rinv_r, xr_scaled);
    correlation_matrix(r, xr_scaled);
  }
  //right now g_minus_G_Rinv_r actually holds g

  MtxDbl Rinv_r(numRowsR,nptsxr);
  MtxDbl G_Rinv_Gtran_inv_g_minus_G_Rinv_r(nTrend,nptsxr);

  solve_after_Chol_fact(Rinv_r,RChol,r);
  matrix_mult(g_minus_G_Rinv_r,Rinv_Gtran,r,1.0,-1.0,'T','N');
  //g_minus_G_Rinv_r now holds g-G*R^-1*r (i.e. it's name is correct)

  solve_after_Chol_fact(G_Rinv_Gtran_inv_g_minus_G_Rinv_r,
			G_Rinv_Gtran_Chol,g_minus_G_Rinv_r);

  for(int ipt=0; ipt<nptsxr; ++ipt) {
    //saved 2*nptsxr loops
    adj_var(0,ipt)=1.0-r(0,ipt)*Rinv_r(0,ipt)+
      g_minus_G_Rinv_r(0,ipt)*G_Rinv_Gtran_inv_g_minus_G_Rinv_r(0,ipt);

    //looks a lot like matrix mult but only N^2 ops... it's the diagonal
    //of a matrix matrix multiply
    for(int iR=1; iR<numRowsR; ++iR)
      adj_var(0,ipt)-=r(iR,ipt)*Rinv_r(iR,ipt);

    //looks a lot like matrix mult but only N^2 ops ... it's the diagonal
    //of a matrix matrix multiply
    for(int itrend=1; itrend<nTrend; ++itrend)
      adj_var(0,ipt)+=g_minus_G_Rinv_r(itrend,ipt)*
	G_Rinv_Gtran_inv_g_minus_G_Rinv_r(itrend,ipt);

    adj_var(0,ipt)*=unscaled_unadj_var;

    if(adj_var(0,ipt)<0.0)
      adj_var(0,ipt)=std::fabs(adj_var(0,ipt)); //zero to within round off and the magnitude of the negative value will give us an idea of how big round off is
    else if(!(adj_var(0,ipt)>=0.0))
      printf("MtxDbl& NKM_KrigingModel::eval_variance(...) adj_var(%d)=nan rcondR=%g\n",ipt,rcondR);
  }

  return adj_var;
}

/** set R=(R+nug*I), where the original R is the correlation matrix for the
    data that the model is built from.  For GEK this generalizes to
    R(i,i)=R(i,i)*(1+nug); Modifying the correlation matrix by the inclusion
    of a nugget causes the KrigingModel to smooth the data, i.e. approximate
    it rather than interpolate it (which is good if you know how big your
    measurement noise is), it can also be used to fix an ill conditioned
    correlation matrix.  The convention is that capital matrices are for
    data the model is built from, lower case matrices are for arbitrary
    points to evaluate the model at */
void KrigingModel::apply_nugget_build() {
  if(!(nug>0.0)) return;
  //printf("applying nugget=%22.16g\n",nug);

  int nrowsR=R.getNRows();
#ifdef __KRIG_ERR_CHECK__
  assert(nrowsR==R.getNCols());
#endif

  double one_plus_nug=1.0+nug;
  for(int i=0; i<nrowsR; ++i)
    R(i,i)*=one_plus_nug;

  return;
}

// convert from correlation lengths to theta (a.k.a. correlation parameters)
MtxDbl& KrigingModel::get_theta_from_corr_len(MtxDbl& theta,
					      const MtxDbl& corr_len) const{
  theta.newSize(numVarsr,1);
  if(corrFunc==GAUSSIAN_CORR_FUNC)
    for(int ixr=0; ixr<numVarsr; ++ixr)
      theta(ixr,0)=0.5/(corr_len(ixr,0)*corr_len(ixr,0));
  else if(corrFunc==EXP_CORR_FUNC) {
#ifdef __KRIG_ERR_CHECK__
    assert(buildDerOrder==0);
#endif
    for(int ixr=0; ixr<numVarsr; ++ixr)
      theta(ixr,0)=1.0/corr_len(ixr,0);
  }
  else if(corrFunc==POW_EXP_CORR_FUNC) {
#ifdef __KRIG_ERR_CHECK__
    assert(buildDerOrder==0);
#endif
    for(int ixr=0; ixr<numVarsr; ++ixr)
      theta(ixr,0)=1.0/
	(powExpCorrFuncPow*std::pow(corr_len(ixr,0),powExpCorrFuncPow));
  }
  else if(corrFunc==MATERN_CORR_FUNC)
    for(int ixr=0; ixr<numVarsr; ++ixr)
      theta(ixr,0)=std::sqrt(2.0*maternCorrFuncNu)/corr_len(ixr,0);
  else{
    std::cerr << "unknown corrFunc in get_theta_from_corr_len()\n";
    assert(false);
  }
  return theta;
}

// convert from theta (a.k.a. correlation parameters) to correlation lengths
MtxDbl& KrigingModel::get_corr_len_from_theta(MtxDbl& corr_len,
					      const MtxDbl& theta) const{
  corr_len.newSize(numVarsr,1);
  if(corrFunc==GAUSSIAN_CORR_FUNC)
    for(int ixr=0; ixr<numVarsr; ++ixr)
      corr_len(ixr,0)=std::sqrt(0.5/theta(ixr,0));
  else if(corrFunc==EXP_CORR_FUNC) {
#ifdef __KRIG_ERR_CHECK__
    assert(buildDerOrder==0);
#endif
    for(int ixr=0; ixr<numVarsr; ++ixr)
      corr_len(ixr,0)=1.0/theta(ixr,0);
  }
  else if(corrFunc==POW_EXP_CORR_FUNC) {
#ifdef __KRIG_ERR_CHECK__
    assert(buildDerOrder==0);
#endif
    for(int ixr=0; ixr<numVarsr; ++ixr)
      corr_len(ixr,0)=
	std::pow(powExpCorrFuncPow*theta(ixr,0),-1.0/powExpCorrFuncPow);
  }
  else if(corrFunc==MATERN_CORR_FUNC)
    for(int ixr=0; ixr<numVarsr; ++ixr)
      corr_len(ixr,0)=std::sqrt(2.0*maternCorrFuncNu)/theta(ixr,0);
  else{
    std::cerr << "unknown corrFunc in get_theta_from_corr_len()\n";
    assert(false);
  }
  return corr_len;
}




/** the inline function
    MtxDbl& KrigingModel::correlation_matrix(MtxDbl& r, const MtxDbl& xr) const
    calls either
    MtxDbl& KrigingModel::eval_kriging_correlation_matrix(MtxDbl& r, const MtxDbl& xr) const (i.e. this function)
    OR
    MtxDbl& KrigingModel::eval_gek_correlation_matrix(MtxDbl& r, const MtxDbl& xr) const

    r (lower case r) is the Kriging correlation matrix between the
    interpolation points and build data points, it used to EVALUATE but
    not construct the emulator's Gaussian process error model
    i.e. E(y(xr)|Y(XR))=betaHat^T*g(xr)+rhs^T*r  where
    rhs=R^-1*(Y-G(XR)^T*betaHat)
    choices for correlation function are gaussian, exponential,
    powered exponential with 1<power<2, matern with nu=1.5 or 2.5
    KRD wrote this */
MtxDbl& KrigingModel::eval_kriging_correlation_matrix(MtxDbl& r, const MtxDbl& xr) const
{
  if(buildDerOrder!=0) {
    std::cerr << "You should only call eval_kriging_correlation_matrix when you want to evaluate regular Kriging (not GEK)\n";
    assert(buildDerOrder==0);
  }

  int nptsxr=xr.getNCols(); //points at which we are evalutating the model
#ifdef __KRIG_ERR_CHECK__
  //  std::cerr<< "xr.getNRows()=" << xr.getNRows()
  //	   << " numVarsr=" << numVarsr
  //	   << " nptsxr=" << nptsxr << std::endl;

  assert((xr.getNRows()==numVarsr)&&(0<nptsxr));
#endif
  r.newSize(numRowsR,nptsxr);
  int i; //row index of the Kriging r matrix (also reorderd XR point index)
  int j; //column index of the Kriging r matrix (also xr point index)
  int k; //dimension index
  double deltax;

  if(corrFunc==GAUSSIAN_CORR_FUNC) {
    // ******************************************************************
    // the Gaussian Correlation function
    // ******************************************************************
    if(numVarsr==1) {
      //special case for when there is only 1 input variable
      double theta=correlations(0,0);
      for(j=0; j<nptsxr; ++j)
	for(i=0; i<numPointsKeep; ++i) {
	  deltax=xr(0,j)-XRreorder(0,i);
	  r(i,j)=std::exp(-theta*deltax*deltax);
	}
    } else {
      //general case there is more than 1 input variable
      //even if nptsxr==1 outer looping once isn't a big performance hit
      //so don't duplicate the code; smallest, i.e. k, loop is inside but
      //that enables a single writing pass through the output array "r"
      double sum_neg_theta_dx_squared;
      for(j=0; j<nptsxr; ++j)
	for(i=0; i<numPointsKeep; ++i) {
	  deltax=xr(0,j)-XRreorder(0,i);
	  sum_neg_theta_dx_squared=-correlations(0,0)* //=- is correct
	    deltax*deltax;
	  for(k=1; k<numVarsr-1; ++k) {
	    deltax=xr(k,j)-XRreorder(k,i);
	    sum_neg_theta_dx_squared-=correlations(k,0)* //-= is correct
	      deltax*deltax;
	  }
	  k=numVarsr-1;
	  deltax=xr(k,j)-XRreorder(k,i);
	  r(i,j)=std::exp(sum_neg_theta_dx_squared
			  -correlations(k,0)*deltax*deltax);

	}
    }
  } else if(corrFunc==EXP_CORR_FUNC) {
    // ******************************************************************
    // the exponential correlation function
    // ******************************************************************
    if(numVarsr==1) {
      //special case for when there is only 1 input variable
      double theta=correlations(0,0);
      for(j=0; j<nptsxr; ++j)
	for(i=0; i<numPointsKeep; ++i)
	  r(i,j)=std::exp(-theta*std::fabs(xr(0,j)-XRreorder(0,i)));
    }
    else {
      //general case there is more than 1 input variable
      //even if nptsxr==1 outer looping once isn't a big performance hit
      //so don't duplicate the code; smallest, i.e. k, loop is inside but
      //that enables a single writing pass through the output array "r"
      double sum_neg_theta_abs_dx;
      for(j=0; j<nptsxr; ++j)
	for(i=0; i<numPointsKeep; ++i) {
	  sum_neg_theta_abs_dx=-correlations(0,0)* //=- is correct
	    std::fabs(xr(0,j)-XRreorder(0,i));
	  for(k=1; k<numVarsr-1; ++k)
	    sum_neg_theta_abs_dx-=correlations(k,0)* //-= is correct
	      std::fabs(xr(k,j)-XRreorder(k,i));
	  k=numVarsr-1;
	  r(i,j)=std::exp(sum_neg_theta_abs_dx
			  -correlations(k,0)*std::fabs(xr(k,j)-XRreorder(k,i)));
	}
    }
  } else if(corrFunc==POW_EXP_CORR_FUNC) {
    // ******************************************************************
    // the powered exponential correlation function 1<powExpCorrFuncPow<2
    // because exponention and Gaussian (a.k.a. squared exponential) were
    // pulled out
    // ******************************************************************
    if(numVarsr==1) {
      //special case for when there is only 1 input variable
      double theta=correlations(0,0);
      for(i=0; i<numPointsKeep; ++i)
	for(j=0; j<nptsxr; ++j)
	  r(i,j)=std::exp(-theta*std::pow(std::fabs(xr(0,j)-XRreorder(0,i)),
					  powExpCorrFuncPow));
    } else {
      //general case there is more than 1 input variable
      //even if nptsxr==1 outer looping once isn't a big performance hit
      //so don't duplicate the code; smallest, i.e. k, loop is inside but
      //that enables a single writing pass through the output array "r"
      double sum_neg_theta_abs_dx_pow;
      for(j=0; j<nptsxr; ++j)
	for(i=0; i<numPointsKeep; ++i) {
	  sum_neg_theta_abs_dx_pow=-correlations(0,0)* //=- is correct
	    std::pow(std::fabs(xr(0,j)-XRreorder(0,i)),powExpCorrFuncPow);
	  for(k=1; k<numVarsr-1; ++k)
	    sum_neg_theta_abs_dx_pow-=correlations(k,0)* //-= is correct
	      std::pow(std::fabs(xr(k,j)-XRreorder(k,i)),powExpCorrFuncPow);
	  k=numVarsr-1;
	  r(i,j)=std::exp(sum_neg_theta_abs_dx_pow-correlations(k,0)*
			  std::pow(std::fabs(xr(k,j)-XRreorder(k,i)),
				   powExpCorrFuncPow));
	}
    }
  } else if((corrFunc==MATERN_CORR_FUNC)&&(maternCorrFuncNu==1.5)) {
    // ******************************************************************
    // the Matern 3/2 Correlation function
    // ******************************************************************
    double theta_abs_dx;
    if(numVarsr==1) {
      //special case for when there is only 1 input variable
      double theta=correlations(0,0);
      for(i=0; i<numPointsKeep; ++i)
	for(j=0; j<nptsxr; ++j) {
	  theta_abs_dx=theta*std::fabs(xr(0,j)-XRreorder(0,i));
	  r(i,j)=(1.0+theta_abs_dx)*std::exp(-theta_abs_dx);
	}
    } else {
      //general case there is more than 1 input variable
      //even if nptsxr==1 outer looping once isn't a big performance hit
      //so don't duplicate the code; smallest, i.e. k, loop is inside but
      //that enables a single writing pass through the output array "r"
      double sum_neg_theta_abs_dx;
      double matern_coef_prod;
      for(j=0; j<nptsxr; ++j)
	for(i=0; i<numPointsKeep; ++i) {
	  theta_abs_dx=correlations(0,0)*std::fabs(xr(0,j)-XRreorder(0,i));
	  matern_coef_prod=1.0+theta_abs_dx;
	  sum_neg_theta_abs_dx=-theta_abs_dx; //=- is correct
	  for(k=1; k<numVarsr-1; ++k) {
	    theta_abs_dx=correlations(k,0)*std::fabs(xr(k,j)-XRreorder(k,i));
	    matern_coef_prod*=(1.0+theta_abs_dx);
	    sum_neg_theta_abs_dx-=theta_abs_dx; //-= is correct
	  }
	  k=numVarsr-1;
	  theta_abs_dx=correlations(k,0)*std::fabs(xr(k,j)-XRreorder(k,i));
	  r(i,j)=matern_coef_prod*(1.0+theta_abs_dx)*
	    std::exp(sum_neg_theta_abs_dx-theta_abs_dx);
	}
    }
  } else if((corrFunc==MATERN_CORR_FUNC)&&(maternCorrFuncNu==2.5)) {
    // ******************************************************************
    // the Matern 5/2 Correlation function
    // ******************************************************************
    double theta_abs_dx;
    const double one_third=1.0/3.0;
    if(numVarsr==1) {
      //special case for when there is only 1 input variable
      double theta=correlations(0,0);
      for(i=0; i<numPointsKeep; ++i)
	for(j=0; j<nptsxr; ++j) {
	  theta_abs_dx=theta*std::fabs(xr(0,j)-XRreorder(0,i));
	  r(i,j)=(1.0+theta_abs_dx+theta_abs_dx*theta_abs_dx*one_third)*
	    std::exp(-theta_abs_dx);
	}
    } else {
      //general case there is more than 1 input variable
      //even if nptsxr==1 outer looping once isn't a big performance hit
      //so don't duplicate the code; smallest, i.e. k, loop is inside but
      //that enables a single writing pass through the output array "r"
      double sum_neg_theta_abs_dx;
      double matern_coef_prod;
      for(j=0; j<nptsxr; ++j)
	for(i=0; i<numPointsKeep; ++i) {
	  theta_abs_dx=correlations(0,0)*std::fabs(xr(0,j)-XRreorder(0,i));
	  matern_coef_prod=1.0+theta_abs_dx+theta_abs_dx*theta_abs_dx*one_third;
	  sum_neg_theta_abs_dx=-theta_abs_dx; //=- is correct
	  for(k=1; k<numVarsr-1; ++k) {
	    theta_abs_dx=correlations(k,0)*std::fabs(xr(k,j)-XRreorder(k,i));
	    matern_coef_prod*=
	      (1.0+theta_abs_dx+theta_abs_dx*theta_abs_dx*one_third);
	    sum_neg_theta_abs_dx-=theta_abs_dx; //-= is correct
	  }
	  k=numVarsr-1;
	  theta_abs_dx=correlations(k,0)*std::fabs(xr(k,j)-XRreorder(k,i));
	  r(i,j)=matern_coef_prod*
	    (1.0+theta_abs_dx+theta_abs_dx*theta_abs_dx*one_third)*
	    std::exp(sum_neg_theta_abs_dx-theta_abs_dx);
	}
    }
  } else{
      std::cerr << "unknown corrFunc in MtxDbl& eval_kriging_correlation_matrix(MtxDbl& r, const MtxDbl& xr) const\n";
      assert(false);
  }

  return r;
}


/** the inline function
    MtxDbl& KrigingModel::correlation_matrix(MtxDbl& r, const MtxDbl& xr) const
    calls either
    MtxDbl& KrigingModel::eval_kriging_correlation_matrix(MtxDbl& r, const MtxDbl& xr) const
    OR
    MtxDbl& KrigingModel::eval_gek_correlation_matrix(MtxDbl& r, const MtxDbl& xr) const (i.e. this function)

    r (lower case r) is the GEK correlation matrix between the
    interpolation points and build data points, it used to EVALUATE but
    not construct the emulator's Gaussian process error model
    i.e. E(y(xr)|Y(XR))=g(xr)^T*betaHat+r^T*R^-1*eps where
    eps=(Y-G(XR)^T*betaHat)
    choices for correlation function are gaussian, and matern with
    nu=1.5 or 2.5 KRD wrote this */
MtxDbl& KrigingModel::eval_gek_correlation_matrix(MtxDbl& r, const MtxDbl& xr) const
{
  if(buildDerOrder!=1) {
    std::cerr << "You should only call eval_gek_correlation_matrix when you want to evaluate Gradient Enhanced Kriging\n";
    assert(buildDerOrder==1);
  }

  int nptsxr=xr.getNCols(); //points at which we are evalutating the model
#ifdef __KRIG_ERR_CHECK__
  assert((xr.getNCols()==numVarsr)&&(0<nptsxr));
#endif

  r.newSize(numRowsR,nptsxr);
  int i; //row index of the GEK r matrix 0<=i<numRowsR
  int j; //column index of the GEK r matrix 0<=j<nptsxr also the xr point index
  int k; //dimension index 0<=k<numVarsr (num real variables)
  int ipt; //XRreorder point index
  double deltax;  //xr(k,j)-XRreorder(k,ipt)
  double krig_r;

  //note to future developers on a point that may be confusing otherwise (you
  //might otherwise mistake this for a bug) all of the derivatives in this
  //file are with respect to XR not xr, the matern_1pt5_d1_mult_r and
  //matern_2pt5_d1_mult_r functions in the .hpp file would be for derivatives
  //with respect to xr (assuming deltax=xr-XR), the difference is here I've
  //absorbed the negative into deltax going from -deltax (where deltax=XR-xr)
  //to deltax=xr-XR;
  int neqn_per_pt=numVarsr+1;

  if(corrFunc==GAUSSIAN_CORR_FUNC) {
    if(numVarsr==1) {
      double theta=correlations(0,0); //save matrix access lookup
      double two_theta = 2.0*theta;
      for(j=0; j<nptsxr; ++j) {
	for(ipt=0, i=0; ipt<numWholePointsKeep; ++ipt, i+=2) {
	  deltax=(xr(0,j)-XRreorder(0,ipt));
	  krig_r=std::exp(-theta*deltax*deltax);
	  r(i  ,j)=krig_r;
	  r(i+1,j)=two_theta*deltax*krig_r; //this is a first
	  //derivative with respect to XR not xr
	}
	if(numPointsKeep>numWholePointsKeep) {
	  //since there's part of another point left and we know that
	  //there is only one derivative it means that were missing that
	  //derivative and only have the function value
#ifdef __KRIG_ERR_CHECK__
	  assert((ipt==numWholePointsKeep)&&
		 (ipt==numPointsKeep-1)&&
		 (i==neqn_per_pt*numWholePointsKeep)&&
		 (i==numRowsR-1));
#endif
	  deltax=(xr(0,j)-XRreorder(0,ipt));
	  r(i,j)=std::exp(-theta*deltax*deltax);
	}
      }
    } else{ //there is more than 1 dimensions and more than one
      //evaluation point
      for(j=0; j<nptsxr; ++j) {
	for(ipt=0, i=0; ipt<numWholePointsKeep; ++ipt, i+=neqn_per_pt) {
	  deltax=xr(0,j)-XRreorder(0,ipt);
	  r(i+1,j)=correlations(0,0)*deltax; //dr_dXR/(2*r)
	  krig_r=-correlations(0,0)*deltax*deltax; //=- is correct
	  for(k=1; k<numVarsr-1; ++k) {
	    deltax=xr(k,j)-XRreorder(k,ipt);
	    r(i+1+k,j)=correlations(k,0)*deltax; //dr_dXR/(2*r)
	    krig_r-=correlations(k,0)*deltax*deltax; //-= is correct
	  }
	  k=numVarsr-1;
	  deltax=xr(k,j)-XRreorder(k,ipt);
	  krig_r=std::exp(krig_r-correlations(k,0)*deltax*deltax);
	  r(i,j)=krig_r; //r(XR(i,:),xr(j,:)) (the correlation function)
	  krig_r*=2.0; //now it's 2*kriging's correlation function to save
	  //some ops
	  //dr_dXR_k=2*theta(k)*(xr(k,j)-XRreorder(k,ipt))*r(xr(k,j),XRreorder(k,ipt))
	  r(i+1+k,j)=correlations(k,0)*deltax*krig_r; //dr_dXR
	  for(k=0; k<numVarsr-1; ++k)
	    r(i+1+k,j)*=krig_r; //dr_dXR
	}
	if(numPointsKeep>numWholePointsKeep) {
	  //the last XR point isn't a "whole point" we dropped some derivatives
	  //out of its gradient to meet the bound on rcond, numExtraDerKeep
	  //is the number of derivatives kept for the last point.  The
	  //derivatives of the last point have NOT been reordered, they appear
	  //in the same order as the input variables
#ifdef __KRIG_ERR_CHECK__
	  assert((ipt==numWholePointsKeep)&&
		 (ipt==numPointsKeep-1)&&
		 (i==neqn_per_pt*numWholePointsKeep));
#endif
	  //printf("deltax=xr(0,%d)-XRreorder(0,%d);\n",j,ipt);
	  deltax=xr(0,j)-XRreorder(0,ipt);
	  krig_r=-correlations(0,0)*deltax*deltax; //=- is correct
	  for(k=1; k<numVarsr-1; ++k) {
	    deltax=xr(k,j)-XRreorder(k,ipt);
	    krig_r-=correlations(k,0)*deltax*deltax; //-= is correct
	  }
	  k=numVarsr-1;
	  deltax=xr(k,j)-XRreorder(k,ipt);
	  krig_r=std::exp(krig_r-correlations(k,0)*deltax*deltax);
	  r(i,j)=krig_r; //r(XR(i,:),xr) (the correlation function)
	  krig_r*=2.0; //now it's 2*kriging's correlation function to save
	  //some ops
	  //dr_dXR_k=2*theta(k)*(xr(k,j)-XRreorder(k,ipt))*r(xr(k,j),XRreorder(k,ipt))
	  for(k=0; k<numExtraDerKeep; ++k)
	    r(i+1+k,j)=correlations(k,0)*(xr(k,j)-XRreorder(k,ipt))*krig_r; //dr_dXR
	}
      }
    }
  } else if((corrFunc==MATERN_CORR_FUNC)&&(maternCorrFuncNu==1.5)) {
    //this starts the section for the Matern 3/2 correlation function

    double theta_abs_dx;
    if(numVarsr==1) {
      double theta=correlations(0,0); //save array access lookup
      double theta_squared= theta*theta;
      double exp_neg_theta_abs_dx;
      for(j=0; j<nptsxr; ++j) {
	for(ipt=0, i=0; ipt<numWholePointsKeep; ++ipt, i+=2) {
	  deltax=(xr(0,j)-XRreorder(0,ipt));
	  theta_abs_dx=theta*std::fabs(deltax);
	  exp_neg_theta_abs_dx=std::exp(-theta_abs_dx);
	  r(i  ,j)=(1.0+theta_abs_dx)*exp_neg_theta_abs_dx; //1D correlation
	  //function
	  r(i+1,j)=theta_squared*deltax*exp_neg_theta_abs_dx; //this is a first
	  //derivative with respect to XR not xr
	}
	if(numPointsKeep>numWholePointsKeep) {
	  //since there's part of another point left and we know that
	  //there is only one derivative it means that were missing that
	  //derivative and only have the function value
#ifdef __KRIG_ERR_CHECK__
	  assert((ipt==numWholePointsKeep)&&
		 (ipt==numPointsKeep-1)&&
		 (i==neqn_per_pt*numWholePointsKeep)&&
		 (i==numRowsR-1));
#endif
	  theta_abs_dx=theta*std::fabs(xr(0,j)-XRreorder(0,ipt));
	  r(i,j)=(1.0+theta_abs_dx)*std::exp(-theta_abs_dx); //1D correlation
	  //function
	}
      }
    }
    else{ //there is more than 1 dimension
      double matern_coef, matern_coef_prod, sum_neg_theta_abs_dx;
      for(j=0; j<nptsxr; ++j) {
	for(ipt=0, i=0; ipt<numWholePointsKeep; ++ipt, i+=neqn_per_pt) {
	  deltax=xr(0,j)-XRreorder(0,ipt);
	  theta_abs_dx=correlations(0,0)*std::fabs(deltax);
	  matern_coef=1.0+theta_abs_dx;
	  matern_coef_prod=matern_coef;
	  r(i+1,j)= //dr_dXR/r
	    correlations(0,0)*correlations(0,0)*deltax/matern_coef;
	  sum_neg_theta_abs_dx=-theta_abs_dx; //=- is correct
	  for(k=1; k<numVarsr-1; ++k) {
	    deltax=xr(k,j)-XRreorder(k,ipt);
	    theta_abs_dx=correlations(k,0)*std::fabs(deltax);
	    matern_coef=1.0+theta_abs_dx;
	    matern_coef_prod*=matern_coef;
	    r(i+1+k,j)= //dr_dXR/r
	      correlations(k,0)*correlations(k,0)*deltax/matern_coef;
	    sum_neg_theta_abs_dx-=theta_abs_dx; //-= is correct
	  }
	  k=numVarsr-1;
	  deltax=xr(k,j)-XRreorder(k,ipt);
	  theta_abs_dx=correlations(k,0)*std::fabs(deltax);
	  matern_coef=1.0+theta_abs_dx;
	  krig_r=matern_coef_prod*matern_coef*
	    std::exp(sum_neg_theta_abs_dx-theta_abs_dx);
	  r(i,j)=krig_r; //r(XR(i,:),xr) (the correlation function)
	  r(i+1+k,j)= //dr_dXR
	    correlations(k,0)*correlations(k,0)*deltax/matern_coef*krig_r;
	  for(k=0; k<numVarsr-1; ++k)
	    r(i+1+k,j)*=krig_r; //dr_dXR
	}
	if(numPointsKeep>numWholePointsKeep) {
	  //the last XR point isn't a "whole point" we dropped some derivatives
	  //out of its gradient to meet the bound on rcond, numExtraDerKeep
	  //is the number of derivatives kept for the last point.  The
	  //derivatives of the last point have NOT been reordered, they appear
	  //in the same order as the input variables
#ifdef __KRIG_ERR_CHECK__
	  assert((ipt==numWholePointsKeep)&&
		 (ipt==numPointsKeep-1)&&
		 (i==neqn_per_pt*numWholePointsKeep));
#endif
	  theta_abs_dx=correlations(0,0)*std::fabs(xr(0,j)-XRreorder(0,ipt));
	  matern_coef_prod=1.0+theta_abs_dx;
	  sum_neg_theta_abs_dx=-theta_abs_dx; //=- is correct
	  for(k=1; k<numVarsr-1; ++k) {
	    theta_abs_dx=correlations(k,0)*std::fabs(xr(k,j)-XRreorder(k,ipt));
	    matern_coef_prod*=(1.0+theta_abs_dx);
	    sum_neg_theta_abs_dx-=theta_abs_dx; //-= is correct
	  }
	  k=numVarsr-1;
	  theta_abs_dx=correlations(k,0)*std::fabs(xr(k,j)-XRreorder(k,ipt));
	  krig_r=matern_coef_prod*(1.0+theta_abs_dx)*
	    std::exp(sum_neg_theta_abs_dx-theta_abs_dx);
	  r(i,j)=krig_r; //r(XR(i,:),xr) (the correlation function)
	  for(k=0; k<numExtraDerKeep; ++k) {
	    deltax=xr(k,j)-XRreorder(k,ipt);
	    r(i+1+k,j)=krig_r * //r(i+1+k,j)=dr_dXR
	      correlations(k,0)*correlations(k,0)*deltax/
	      (1.0+correlations(k,0)*std::fabs(deltax));
	  }
	}
      }
    }
  } else if((corrFunc==MATERN_CORR_FUNC)&&(maternCorrFuncNu==2.5)) {
    //this starts the section for the Matern 5/2 correlation function

    const double one_third=1.0/3.0;
    double theta_abs_dx;
    if(numVarsr==1) {
      double theta=correlations(0,0); //save array access lookup
      double theta_squared= theta*theta;
      double exp_neg_theta_abs_dx;
      for(j=0; j<nptsxr; ++j) {
	for(ipt=0, i=0; ipt<numWholePointsKeep; ++ipt, i+=2) {
	  deltax=(xr(0,j)-XRreorder(0,ipt));
	  theta_abs_dx=theta*std::fabs(deltax);
	  exp_neg_theta_abs_dx=std::exp(-theta_abs_dx);
	  r(i  ,j)=(1.0+theta_abs_dx+theta_abs_dx*theta_abs_dx*one_third)*
	    exp_neg_theta_abs_dx; //1D correlation function
	  r(i+1,j)=theta_squared*deltax*(1.0+theta_abs_dx)*one_third*
	    exp_neg_theta_abs_dx; //this is a first derivative with respect
	  //to XR not xr
	}
	if(numPointsKeep>numWholePointsKeep) {
	  //since there's part of another point left and we know that
	  //there is only one derivative it means that were missing that
	  //derivative and only have the function value
#ifdef __KRIG_ERR_CHECK__
	  assert((ipt==numWholePointsKeep)&&
		 (ipt==numPointsKeep-1)&&
		 (i==neqn_per_pt*numWholePointsKeep)&&
		 (i==numRowsR-1));
#endif
	  theta_abs_dx=theta*std::fabs(xr(0,j)-XRreorder(0,ipt));
	  r(i  ,j)=(1.0+theta_abs_dx+theta_abs_dx*theta_abs_dx*one_third)*
	    std::exp(-theta_abs_dx); //1D correlation function
	}
      }
    }
    else{ //there is more than 1 dimension
      double matern_coef, matern_coef_prod, sum_neg_theta_abs_dx;
      for(j=0; j<nptsxr; ++j) {
	for(ipt=0, i=0; ipt<numWholePointsKeep; ++ipt, i+=neqn_per_pt) {
	  deltax=xr(0,j)-XRreorder(0,ipt);
	  theta_abs_dx=correlations(0,0)*std::fabs(deltax);
	  matern_coef=1.0+theta_abs_dx+theta_abs_dx*theta_abs_dx*one_third;
	  matern_coef_prod=matern_coef;
	  r(i+1,j)= //dr_dXR/r
	    correlations(0,0)*correlations(0,0)*deltax*(1.0+theta_abs_dx)*
	    one_third/matern_coef;
	  sum_neg_theta_abs_dx=-theta_abs_dx; //=- is correct
	  for(k=1; k<numVarsr-1; ++k) {
	    deltax=xr(k,j)-XRreorder(k,ipt);
	    theta_abs_dx=correlations(k,0)*std::fabs(deltax);
	    matern_coef=1.0+theta_abs_dx+theta_abs_dx*theta_abs_dx*one_third;
	    matern_coef_prod*=matern_coef;
	    r(i+1+k,j)= //dr_dXR/r
	      correlations(k,0)*correlations(k,0)*deltax*(1.0+theta_abs_dx)*
	      one_third/matern_coef;
	    sum_neg_theta_abs_dx-=theta_abs_dx; //-= is correct
	  }
	  k=numVarsr-1;
	  deltax=xr(k,j)-XRreorder(k,ipt);
	  theta_abs_dx=correlations(k,0)*std::fabs(deltax);
	  matern_coef=1.0+theta_abs_dx+theta_abs_dx*theta_abs_dx*one_third;
	  krig_r=matern_coef_prod*matern_coef*
	    std::exp(sum_neg_theta_abs_dx-theta_abs_dx);
	  r(i,j)=krig_r; //r(XR(i,:),xr) (the correlation function)
	  r(i+1+k,j)= //dr_dXR
	    correlations(k,0)*correlations(k,0)*deltax*(1.0+theta_abs_dx)*
	    one_third/matern_coef*krig_r;
	  for(k=0; k<numVarsr-1; ++k)
	    r(i+1+k,j)*=krig_r; //dr_dXR
	}
	if(numPointsKeep>numWholePointsKeep) {
	  //the last XR point isn't a "whole point" we dropped some derivatives
	  //out of its gradient to meet the bound on rcond, numExtraDerKeep
	  //is the number of derivatives kept for the last point.  The
	  //derivatives of the last point have NOT been reordered, they appear
	  //in the same order as the input variables
#ifdef __KRIG_ERR_CHECK__
	  assert((ipt==numWholePointsKeep)&&
		 (ipt==numPointsKeep-1)&&
		 (i==neqn_per_pt*numWholePointsKeep));
#endif
	  theta_abs_dx=correlations(0,0)*std::fabs(xr(0,j)-XRreorder(0,ipt));
	  matern_coef_prod=1.0+theta_abs_dx+theta_abs_dx*theta_abs_dx*one_third;
	  sum_neg_theta_abs_dx=-theta_abs_dx; //=- is correct
	  for(k=1; k<numVarsr-1; ++k) {
	    theta_abs_dx=correlations(k,0)*std::fabs(xr(k,j)-XRreorder(k,ipt));
	    matern_coef_prod*=
	      (1.0+theta_abs_dx+theta_abs_dx*theta_abs_dx*one_third);
	    sum_neg_theta_abs_dx-=theta_abs_dx; //-= is correct
	  }
	  k=numVarsr-1;
	  theta_abs_dx=correlations(k,0)*std::fabs(xr(k,j)-XRreorder(k,ipt));
	  krig_r=matern_coef_prod*
	    (1.0+theta_abs_dx+theta_abs_dx*theta_abs_dx*one_third)*
	    std::exp(sum_neg_theta_abs_dx-theta_abs_dx);
	  r(i,j)=krig_r; //r(XR(i,:),xr) (the correlation function)
	  for(k=0; k<numExtraDerKeep; ++k) {
	    deltax=xr(k,j)-XRreorder(k,ipt);
	    theta_abs_dx=correlations(k,0)*std::fabs(deltax);
	    r(i+1+k,j)=krig_r * //r(i+1+k,j)=dr_dXR
	      correlations(k,0)*correlations(k,0)*deltax*(1.0+theta_abs_dx)/
	      (3.0*(1.0+theta_abs_dx)+theta_abs_dx*theta_abs_dx);

	  }
	}
      }
    }
  } else{
    std::cerr << "Unknown or Invalid Correlation function for Gradient Enhanced Kriging in MtxDbl& KrigingModel::eval_gek_correlation_matrix(MtxDbl& r, const MtxDbl& xr) const\n";
    assert(false);
  }


  return r;
}


///Ider is the variable/dimension not the point
MtxDbl& KrigingModel::eval_kriging_dcorrelation_matrix_dxI(MtxDbl& dr, const MtxDbl& r, const MtxDbl& xr, int Ider) const
{
  if(buildDerOrder!=0) {
    std::cerr << "You should only call eval_kriging_dcorrelation_matrix_dxI when you want to evaluate regular Kriging's (not GEK's) first derivative.\n";
    assert(buildDerOrder==0);
  }
  int nptsxr=xr.getNCols();
#ifdef __KRIG_ERR_CHECK__
  assert((r.getNCols()==nptsxr)&&(r.getNRows()==numRowsR)&&
	 (xr.getNRows()==numVarsr)&&(0<=Ider)&&(Ider<numVarsr));
#endif
  dr.newSize(numRowsR,nptsxr);
  int i; //row index of r & dr, also the point index of reordered XR
  int j; //column index of r & dr, also the point index of xr

  if(corrFunc==GAUSSIAN_CORR_FUNC) {
    // *******************************************************************
    // Gaussian Correlation Function
    // GAUSSIAN_CORR_FUNC is infinitely differentiable
    // *******************************************************************
    double neg_two_theta=-2.0*correlations(Ider,0); //save matrix dereference
    //for speed
    for(j=0; j<nptsxr; ++j)
      for(i=0; i<numPointsKeep; ++i)
	dr(i,j)=r(i,j)*neg_two_theta*(xr(Ider,j)-XRreorder(Ider,i));
  } else if(corrFunc==EXP_CORR_FUNC) {
    // *******************************************************************
    // Exponential Correlation Function
    // 1D EXP_CORR_FUNC r(x1,x2) is differentiable except where x1==x2
    // this is correct for x1!=x2
    // *******************************************************************
    double neg_theta=-correlations(Ider,0); //save matrix dereference for
    //speed
    for(j=0; j<nptsxr; ++j)
      for(i=0; i<numPointsKeep; ++i)
	dr(i,j)=r(i,j)*neg_theta*dsign(xr(Ider,j)-XRreorder(Ider,i));
  } else if(corrFunc==POW_EXP_CORR_FUNC) {
    // *******************************************************************
    // Powered Exponential Correlation Function with 1<power<2
    // 1D POW_EXP_CORR_FUNC r(x1,x2) is once differential everywhere (and
    // twice+ differentiable where x1!=x2)
    // *******************************************************************
    double neg_theta_pow=-powExpCorrFuncPow*correlations(Ider,0); //save
    //matrix dereference for speed
    double pow_m_1=powExpCorrFuncPow-1.0; //for speed
    double delta_x;
    for(int j=0; j<nptsxr; ++j)
      for(int i=0; i<numPointsKeep; ++i) {
	delta_x=xr(Ider,j)-XRreorder(Ider,i);
	dr(i,j)=r(i,j)*dsign(delta_x)*neg_theta_pow*
	  std::pow(std::fabs(delta_x),pow_m_1);
      }
  } else if((corrFunc==MATERN_CORR_FUNC)&&(maternCorrFuncNu==1.5)) {
    // *******************************************************************
    // Matern 3/2 Correlation Function
    // 1D MATERN_CORR_FUNC 1.5 is once differentiable everywhere (and
    // twice+ differentiable where x1!=x2, while not twice differentiable
    // at x1==x2 the limit of the 2nd derivative is defined and is the
    // same from both sides see Lockwood and Anitescu)
    // *******************************************************************
    double theta=correlations(Ider,0); //save matrix dereference for speed
    for(j=0; j<nptsxr; ++j)
      for(i=0; i<numPointsKeep; ++i)
	dr(i,j)=r(i,j)*
	  matern_1pt5_d1_mult_r(theta,xr(Ider,j)-XRreorder(Ider,i));
  } else if((corrFunc==MATERN_CORR_FUNC)&&(maternCorrFuncNu==2.5)) {
    // *******************************************************************
    // Matern 5/2 Correlation Function
    // 1D MATERN_CORR_FUNC 2.5 is twice differentiable everywhere (and
    // twice+ differentiable where x1!=x2)
    // *******************************************************************
    double theta=correlations(Ider,0); //save matrix dereference for speed
    for(j=0; j<nptsxr; ++j)
      for(i=0; i<numPointsKeep; ++i)
	dr(i,j)=r(i,j)*
	  matern_2pt5_d1_mult_r(theta,xr(Ider,j)-XRreorder(Ider,i));
  } else{
    std::cerr << "unknown corrFunc in MtxDbl& KrigingModel::eval_kriging_dcorrelation_matrix_dxI(MtxDbl& dr, const MtxDbl& r, const MtxDbl& xr, int Ider) const\n";
    assert(false);
  }
  return dr;
}
///Ider is the variable/dimension not the point
MtxDbl& KrigingModel::eval_gek_dcorrelation_matrix_dxI(MtxDbl& dr, const MtxDbl& r, const MtxDbl& xr, int Ider) const
{
  if(buildDerOrder!=1) {
    std::cerr << "You should only call eval_gek_dcorrelation_matrix_dxI when you want to evaluate Gradient Enhanced Kriging's first derivative\n";
    assert(buildDerOrder==1);
  }
  int nptsxr=xr.getNCols();
#ifdef __KRIG_ERR_CHECK__
  assert((r.getNCols()==nptsxr)&&(r.getNRows()==numRowsR)&&
	 (xr.getNRows()==numVarsr)&&(0<=Ider)&&(Ider<numVarsr));
#endif
  dr.newSize(numRowsR,nptsxr);
  int neqn_per_pt=1+numVarsr;
  int i; //row index of r & dr
  int j; //column index of r & dr, also the point index of xr
  int k; //dimension index
  int ipt; //point index of reordered XR

  if(corrFunc==GAUSSIAN_CORR_FUNC) {
    // *******************************************************************
    // Gaussian Correlation Function
    // GAUSSIAN_CORR_FUNC is infinitely differentiable
    // *******************************************************************
    double two_theta=2.0*correlations(Ider,0); //save matrix dereference for speed
    double neg_two_theta_dx;
    if(numVarsr==1)
      for(j=0; j<nptsxr; ++j) {
	for(ipt=0, i=0; ipt<numWholePointsKeep; ++ipt, i+=2) {
	  neg_two_theta_dx=two_theta*(XRreorder(Ider,ipt)-xr(Ider,j));
	  dr(i  ,j)=r(i,j)*neg_two_theta_dx;
	  dr(i+1,j)=r(i,j)*two_theta + r(i+1,j)*neg_two_theta_dx;
	}
	// since there is only one dimension if there is a partial point
	// it will be a function value only, and actually recalculating it
	// will likely be faster on average then checking if there's a
	// partial point and calculating it if needed
	ipt=numPointsKeep-1;
	i=ipt*2;
	dr(i  ,j)=r(i,j)*two_theta*(XRreorder(Ider,ipt)-xr(Ider,j));
      }
    else{
      for(j=0; j<nptsxr; ++j) {
	for(ipt=0, i=0; ipt<numWholePointsKeep; ++ipt, i+=neqn_per_pt) {
	  neg_two_theta_dx=two_theta*(XRreorder(Ider,ipt)-xr(Ider,j));
	  dr(i,j)=r(i,j)*neg_two_theta_dx;
	  for(k=0; k<numVarsr; ++k)
	    dr(i+1+k,j)=r(i+1+k,j)*neg_two_theta_dx;
	  dr(i+1+Ider,j)+=r(i,j)*two_theta;
	}
	if(numPointsKeep>numWholePointsKeep) {
	  //ipt and i should be what we need them to be
#ifdef __KRIG_ERR_CHECK__
	  assert((ipt==numWholePointsKeep)&&
		 (ipt==numPointsKeep-1)&&
		 (i==neqn_per_pt*numWholePointsKeep));
#endif
	  neg_two_theta_dx=two_theta*(XRreorder(Ider,ipt)-xr(Ider,j));
	  dr(i,j)=r(i,j)*neg_two_theta_dx;
	  for(k=0; k<numExtraDerKeep; ++k)
	    dr(i+1+k,j)=r(i+1+k,j)*neg_two_theta_dx;
	  if(Ider<numExtraDerKeep)
	    dr(i+1+Ider,j)+=r(i,j)*two_theta;
	}
      }
    }
  } else if(corrFunc==EXP_CORR_FUNC) {
    std::cerr << "The exponential correlation function is not a valid correlation function for gradient enhanced Kriging\n";
      assert(false);
  } else if(corrFunc==POW_EXP_CORR_FUNC) {
    std::cerr << "The powered exponential (with power < 2) correlation function is not a valid correlation function for gradient enhanced Kriging\n";
      assert(false);
  } else if((corrFunc==MATERN_CORR_FUNC)&&(maternCorrFuncNu==1.5)) {
    // *******************************************************************
    // Matern 3/2 Correlation Function
    // 1D MATERN_CORR_FUNC 1.5 is once differentiable everywhere (and
    // twice+ differentiable where x1!=x2, while not twice differentiable
    // at x1==x2 the limit of the 2nd derivative is defined and is the
    // same from both sides see Lockwood and Anitescu)
    // *******************************************************************
    double theta=correlations(Ider,0); //save matrix dereference for speed
    double neg_theta_squared=-theta*theta;
    double deltax;
    double matern_coef;
    if(numVarsr==1)
      for(j=0; j<nptsxr; ++j) {
	for(ipt=0, i=0; ipt<numWholePointsKeep; ++ipt, i+=2) {
	  deltax=(xr(Ider,j)-XRreorder(Ider,ipt));
	  matern_coef=1.0+theta*std::fabs(deltax);
	  dr(i  ,j)=r(i,j)*neg_theta_squared*deltax/matern_coef;
	  dr(i+1,j)=r(i,j)*neg_theta_squared*(1.0-2.0/matern_coef);
	}
	// since there is only one dimension if there is a partial point
	// it will be a function value only, and actually recalculating it
	// will likely be faster on average then checking if there's a
	// partial point and calculating it if needed
	ipt=numPointsKeep-1;
	i=ipt*2;
	deltax=(xr(Ider,j)-XRreorder(Ider,ipt));
	dr(i,j)=r(i,j)*neg_theta_squared*deltax/(1.0+theta*std::fabs(deltax));
      }
    else{
      double matern_d1_mult_r;
      for(j=0; j<nptsxr; ++j) {
	for(ipt=0, i=0; ipt<numWholePointsKeep; ++ipt, i+=neqn_per_pt) {
	  deltax=(xr(Ider,j)-XRreorder(Ider,ipt));
	  matern_coef=1.0+theta*std::fabs(deltax);
	  matern_d1_mult_r=neg_theta_squared*deltax/matern_coef;
	  dr(i  ,j)=r(i,j)*matern_d1_mult_r;
	  for(k=0; k<numVarsr; ++k)
	    dr(i+1+k,j)=r(i+1+k,j)*matern_d1_mult_r;
	  dr(i+1+Ider,j)=r(i,j)*neg_theta_squared*(1.0-2.0/matern_coef);
	}
	if(numPointsKeep>numWholePointsKeep) {
	  //ipt and i should be what we need them to be
#ifdef __KRIG_ERR_CHECK__
	  assert((ipt==numWholePointsKeep)&&
		 (ipt==numPointsKeep-1)&&
		 (i==neqn_per_pt*numWholePointsKeep));
#endif
	  deltax=(xr(Ider,j)-XRreorder(Ider,ipt));
	  matern_coef=1.0+theta*std::fabs(deltax);
	  matern_d1_mult_r=neg_theta_squared*deltax/matern_coef;
	  dr(i  ,j)=r(i,j)*matern_d1_mult_r;
	  for(k=0; k<numExtraDerKeep; ++k)
	    dr(i+1+k,j)=r(i+1+k,j)*matern_d1_mult_r;
	  if(Ider<numExtraDerKeep)
	    dr(i+1+Ider,j)=r(i,j)*neg_theta_squared*(1.0-2.0/matern_coef);
	}
      }
    }
  } else if((corrFunc==MATERN_CORR_FUNC)&&(maternCorrFuncNu==2.5)) {
    // *******************************************************************
    // Matern 5/2 Correlation Function
    // 1D MATERN_CORR_FUNC 2.5 is twice differentiable everywhere (and
    // twice+ differentiable where x1!=x2)
    // *******************************************************************
    double theta=correlations(Ider,0); //save matrix dereference for speed
    double theta_squared=theta*theta;
    double theta_abs_dx;
    double deltax;
    if(numVarsr==1) {
      double r_theta_squared_div_3_matern_coef;
      for(j=0; j<nptsxr; ++j) {
	for(ipt=0, i=0; ipt<numWholePointsKeep; ++ipt, i+=2) {
	  deltax=(xr(Ider,j)-XRreorder(Ider,ipt));
	  theta_abs_dx=theta*std::fabs(deltax);
	  r_theta_squared_div_3_matern_coef=r(i,j)*theta_squared/
	    (3.0*(1.0+theta_abs_dx)+theta_abs_dx*theta_abs_dx);
	  dr(i  ,j)=-r_theta_squared_div_3_matern_coef*
	    deltax*(1.0+theta_abs_dx);
	  dr(i+1,j)=r_theta_squared_div_3_matern_coef*
	    (1.0+theta_abs_dx-theta_abs_dx*theta_abs_dx);
	}
	// since there is only one dimension if there is a partial point
	// it will be a function value only, and actually recalculating it
	// will likely be faster on average then checking if there's a
	// partial point and calculating it if needed
	ipt=numPointsKeep-1;
	i=ipt*2;
	deltax=(xr(Ider,j)-XRreorder(Ider,ipt));
	theta_abs_dx=theta*std::fabs(deltax);
	dr(i,j)=-r(i,j)*theta_squared/
	  (3.0*(1.0+theta_abs_dx)+theta_abs_dx*theta_abs_dx)*
	  deltax*(1.0+theta_abs_dx);
      }
    } else{
      double theta_squared_div_3_matern_coef;
      double matern_d1_mult_r;
      for(j=0; j<nptsxr; ++j) {
	for(ipt=0, i=0; ipt<numWholePointsKeep; ++ipt, i+=neqn_per_pt) {
	  deltax=xr(Ider,j)-XRreorder(Ider,ipt);
	  theta_abs_dx=theta*std::fabs(deltax);
	  theta_squared_div_3_matern_coef=theta_squared/
	    (3.0*(1.0+theta_abs_dx)+theta_abs_dx*theta_abs_dx);
	  matern_d1_mult_r=-theta_squared_div_3_matern_coef*
	    deltax*(1.0+theta_abs_dx);
	  dr(i  ,j)=r(i,j)*matern_d1_mult_r;
	  for(k=0; k<numVarsr; ++k)
	    dr(i+1+k,j)=r(i+1+k,j)*matern_d1_mult_r;
	  dr(i+1+Ider,j)=r(i,j)*theta_squared_div_3_matern_coef*
	    (1.0+theta_abs_dx-theta_abs_dx*theta_abs_dx);
	}
	if(numPointsKeep>numWholePointsKeep) {
	  //ipt and i should be what we need them to be
#ifdef __KRIG_ERR_CHECK__
	  assert((ipt==numWholePointsKeep)&&
		 (ipt==numPointsKeep-1)&&
		 (i==neqn_per_pt*numWholePointsKeep));
#endif
	  deltax=xr(Ider,j)-XRreorder(Ider,ipt);
	  theta_abs_dx=theta*std::fabs(deltax);
	  theta_squared_div_3_matern_coef=theta_squared/
	    (3.0*(1.0+theta_abs_dx)+theta_abs_dx*theta_abs_dx);
	  matern_d1_mult_r=-theta_squared_div_3_matern_coef*
	    deltax*(1.0+theta_abs_dx);
	  dr(i  ,j)=r(i,j)*matern_d1_mult_r;
	  for(k=0; k<numExtraDerKeep; ++k)
	    dr(i+1+k,j)=r(i+1+k,j)*matern_d1_mult_r;
	  if(Ider<numExtraDerKeep)
	    dr(i+1+Ider,j)=r(i,j)*theta_squared_div_3_matern_coef*
	      (1.0+theta_abs_dx-theta_abs_dx*theta_abs_dx);
	}
      }
    }
  } else{
    std::cerr << "unknown corrFunc in MtxDbl& KrigingModel::eval_gek_dcorrelation_matrix_dxI(MtxDbl& dr, const MtxDbl& r, const MtxDbl& xr, int Ider) const\n";
    assert(false);
  }

  return dr;
}



MtxDbl& KrigingModel::eval_kriging_d2correlation_matrix_dxIdxJ(MtxDbl& d2r, const MtxDbl& drI, const MtxDbl& r, const MtxDbl& xr, int Ider, int Jder) const
{
  if(buildDerOrder!=0) {
    std::cerr << "You should only call eval_kriging_correlation_matrix when you want to evaluate regular Kriging (not GEK)\n";
    assert(buildDerOrder==0);
  }

  int nptsxr=xr.getNCols(); //points at which we are evalutating the model
  d2r.newSize(numPointsKeep,nptsxr);

#ifdef __KRIG_ERR_CHECK__
  assert((r.getNCols()==nptsxr)&&(r.getNRows()==numPointsKeep)&&
	 (xr.getNRows()==numVarsr)&&(0<=Jder)&&(Jder<numVarsr));
#endif

  int i; //row index of r, d1r, & d2r; also the point index of reordered XR
  int j; //column index of r, d1r, & d2r; also the point index of xr

  if(corrFunc==GAUSSIAN_CORR_FUNC) {
    // *********************************************************************
    // The GAUSSIAN CORRELATION FUNCTION
    // is infinitely differentiable, i.e. is C^infinity continuous
    // *********************************************************************
    double neg_two_theta_J=-2.0*correlations(Jder,0);
    if(Ider==Jder) {
      // taking the 2nd derivative of the 1D correlation function
      for(j=0; j<nptsxr; ++j)
	for(i=0; i<numPointsKeep; ++i)
	  d2r(i,j)=neg_two_theta_J*
	    ((xr(Jder,j)-XRreorder(Jder,i))*drI(i,j)+r(i,j));
    } else {
      // taking the product of the 1st derivative of 2 independent 1D
      // correlation functions
      for(j=0; j<nptsxr; ++j)
	for(i=0; i<numPointsKeep; ++i)
	  d2r(i,j)=neg_two_theta_J*
	    (xr(Jder,j)-XRreorder(Jder,i))*drI(i,j);
    }
  } else if(corrFunc==EXP_CORR_FUNC) {
    // *********************************************************************
    // The EXPONENTIAL CORRELATION FUNCTION
    // the first derivative WRT theta(J) is
    //     drJ=-theta(J)*sign(xr(J)-XR(J))*r
    // if away from xr(J)==XR(J) then d(sign(xr(J)-XR(J))/dxr(J)=0
    // it at xr(J)==XR(J) then derivative of sign (a.k.a step function) is
    // two times the delta function (or a rectangle with area 2, whose base
    // width is a point, i.e. zero, meaning the delta function is infinite).
    // The following is correct as long as xr(J)=/=XR(J), i.e. as long as
    // the evaluation point doesn't share a coordinate with any build point.
    // *********************************************************************
    double neg_theta_J=-correlations(Jder,0);
    for(j=0; j<nptsxr; ++j)
      for(i=0; i<numPointsKeep; ++i)
	d2r(i,j)=neg_theta_J*dsign(xr(Jder,j)-XRreorder(Jder,i))*drI(i,j);
  } else if(corrFunc==POW_EXP_CORR_FUNC) {
    // *********************************************************************
    // The POWERED EXPONENTIAL CORRELATION FUNCTION with 1<power<2
    //
    // the 1st derivative with respect to xr of the 1D correlation function
    // is defined
    //
    // The 2nd derivative with respect to xr of the 1D correlation function
    // *is undefined at xr==XR,
    // *approaches negative infinity as xr approaches XR from below, and
    // *approaches positive infinity as xr approaches XR from above
    // when xr==XR we use the average of the second derivative from above and
    // the second derivative from below, that average is exactly zero
    // *********************************************************************
    double neg_thetaJ_pow=-correlations(Jder,0)*powExpCorrFuncPow;
    double pow_minus_1=powExpCorrFuncPow-1.0;
    double abs_dx;
    double deltax;
    if(Ider==Jder) {
      // taking the 2nd derivative of the 1D correlation function
      double pow_minus_2=powExpCorrFuncPow-2.0;
      for(j=0; j<nptsxr; ++j)
	for(i=0; i<numPointsKeep; ++i) {
	  deltax=xr(Jder,j)-XRreorder(Jder,i);
	  if(deltax==0) {
	    d2r(i,j)=0.0;
	    std::cerr << "the 2nd derivative of the powered exponential correlation function (with 1<power<2) is undefined when a coordinate of the evaluation point equals the coordinate of a build point, using the zero as the average of + infinity (from above) and - infinity (from below)\n";
	  }
	  else{
	    abs_dx=std::fabs(deltax);
	    d2r(i,j)=neg_thetaJ_pow*dsign(deltax)*
	      (pow_minus_1*std::pow(abs_dx,pow_minus_2)*r(i,j)+
	       std::pow(abs_dx,pow_minus_1)*drI(i,j));
	  }
	}
    } else {
      //we are taking the product of the first derivatives of 2 independent
      //1D correlation functions so we don't have to worry about the 2nd
      //derivative of a 1D correlation function being undefined
      for(j=0; j<nptsxr; ++j)
	for(i=0; i<numPointsKeep; ++i) {
	  deltax=xr(Jder,j)-XRreorder(Jder,i);
	  d2r(i,j)=neg_thetaJ_pow*dsign(deltax)*
	    std::pow(std::fabs(deltax),pow_minus_1)*drI(i,j);
	}
    }
  } else if((corrFunc==MATERN_CORR_FUNC)&&(maternCorrFuncNu==1.5)) {
    // *********************************************************************
    // The MATERN 3/2 CORRELATION FUNCTION
    //
    // the 1st derivative with respect to xr of the 1D correlation function
    // is defined
    //
    // The 2nd derivative with respect to xr of the 1D correlation function
    // *is undefined at xr==XR,
    // *is -theta^2*(1-theta*|xr-XR|)*exp(-theta*|xr-XR|) at xr=/=XR
    // *approaches -theta^2 from above and below
    // when xr==XR we use the limit, -theta^2
    // this follows the approach of
    //   Lockwood, Brian A. and Anitescu, Mihai, "Gradient-Enhanced
    //      Universal Kriging for Uncertainty Proagation"
    //      Preprint ANL/MCS-P1808-1110
    // *********************************************************************
    double thetaJ=correlations(Jder,0);
    double neg_thetaJ_squared=-thetaJ*thetaJ;
    double deltax;
    if(Ider==Jder) {
      // taking the 2nd derivative of the 1D correlation function
      for(j=0; j<nptsxr; ++j)
	for(i=0; i<numPointsKeep; ++i) {
	  deltax=xr(Jder,j)-XRreorder(Jder,i);
	  d2r(i,j)=neg_thetaJ_squared*
	    (2.0/(1.0+thetaJ*std::fabs(deltax))-1.0)*r(i,j);
	}
    }else{
      // taking the product of the 1st derivative of 2 independent 1D
      // correlation functions
      for(j=0; j<nptsxr; ++j)
	for(i=0; i<numPointsKeep; ++i) {
	  deltax=xr(Jder,j)-XRreorder(Jder,i);
	  d2r(i,j)=neg_thetaJ_squared*deltax/(1.0+thetaJ*std::fabs(deltax))*
	    drI(i,j);
	}
    }
  } else if((corrFunc==MATERN_CORR_FUNC)&&(maternCorrFuncNu==2.5)) {
    // *********************************************************************
    // The MATERN 5/2 CORRELATION FUNCTION
    //
    // the 1st and 2nd derivatives with respect to xr of the 1D correlation
    // function are defined, no special treatment is required
    // *********************************************************************
    double thetaJ=correlations(Jder,0);
    double neg_thetaJ_squared=-thetaJ*thetaJ;
    double deltax;
    double thetaJ_abs_dx;
    if(Ider==Jder) {
      // taking the 2nd derivative of the 1D correlation function
      for(j=0; j<nptsxr; ++j)
	for(i=0; i<numPointsKeep; ++i) {
	  deltax=xr(Jder,j)-XRreorder(Jder,i);
	  thetaJ_abs_dx=thetaJ*std::fabs(deltax);
	  d2r(i,j)=neg_thetaJ_squared*
	    (1.0+thetaJ_abs_dx-thetaJ_abs_dx*thetaJ_abs_dx)/
	    (3.0*(1.0+thetaJ_abs_dx)+thetaJ_abs_dx*thetaJ_abs_dx)*
	    r(i,j);
	}
    }
    else {
      // taking the product of the 1st derivative of 2 independent 1D
      // correlation functions
      for(j=0; j<nptsxr; ++j)
	for(i=0; i<numPointsKeep; ++i) {
	  deltax=xr(Jder,j)-XRreorder(Jder,i);
	  thetaJ_abs_dx=thetaJ*std::fabs(deltax);
	  d2r(i,j)=neg_thetaJ_squared*deltax*(1.0+thetaJ_abs_dx)/
	    (3.0*(1.0+thetaJ_abs_dx)+thetaJ_abs_dx*thetaJ_abs_dx)*
	    drI(i,j);
	}
    }
  } else{
    std::cerr << "unknown corrFunc in MtxDbl& KrigingModel::eval_kriging_d2correlation_matrix_dxIdxJ(MtxDbl& d2r, const MtxDbl& drI, const MtxDbl& r, const MtxDbl& xr, int Ider, int Jder) const\n";
    assert(false);
  }
  return d2r;
}
MtxDbl& KrigingModel::eval_gek_d2correlation_matrix_dxIdxJ(MtxDbl& d2r, const MtxDbl& drI, const MtxDbl& r, const MtxDbl& xr, int Ider, int Jder) const
{
  if(buildDerOrder!=1) {
    std::cerr << "You should only call eval_gek_dcorrelation_matrix_dxI when you want to evaluate Gradient Enhanced Kriging's second derivative\n";
    assert(buildDerOrder==1);
  }
  int nptsxr=xr.getNCols(); //points at which we are evalutating the model
  d2r.newSize(numRowsR,nptsxr);

#ifdef __KRIG_ERR_CHECK__
  assert((r.getNCols()==nptsxr)&&(r.getNRows()==numPointsKeep)&&
	 (xr.getNRows()==numVarsr)&&(0<=Jder)&&(Jder<numVarsr));
#endif

  int i; //row index of r, d1r, & d2r
  int j; //column index of r, d1r, & d2r; also the point index of xr
  int k; //dimension index
  int ipt; //point index of reordered XR
  int neqn_per_pt=numVarsr+1;
  double deltax;

  if(corrFunc==GAUSSIAN_CORR_FUNC) {
    // *********************************************************************
    // The GAUSSIAN CORRELATION FUNCTION
    // is infinitely differentiable, i.e. is C^infinity continuous
    // the reuse lower order derivates formulas are derived by taking
    // derivatives in the reverse order of occurance and not expanding
    // derivatives or r of d1r (here called drI)
    // *********************************************************************
    double neg_two_thetaJ=-2.0*correlations(Jder,0);
    if(numVarsr==1) {
      // if there is only one input variable we are taking the 2nd derivative
      // of the 1D correlation function
      // AND WE KNOW THAT Ider=Jder=k so we don't have to "if" to add the
      // the extra terms
      for(j=0; j<nptsxr; ++j) {
	for(ipt=0, i=0; ipt<numWholePointsKeep; ++ipt, i+=2) {
	  deltax=xr(Jder,j)-XRreorder(Jder,ipt);
	  d2r(i  ,j)=neg_two_thetaJ*(deltax*drI(i  ,j)+r(i  ,j));
	  d2r(i+1,j)=neg_two_thetaJ*(deltax*drI(i+1,j)+r(i+1,j)-drI(i,j));
	}
	// since there is only one dimension if there is a partial point
	// it will be a function value only, and actually recalculating it
	// will likely be faster on average then checking if there's a
	// partial point and calculating it if needed
	ipt=numPointsKeep-1;
	i=ipt*2;
	d2r(i,j)=neg_two_thetaJ*
	  ((xr(Jder,j)-XRreorder(Jder,ipt))*drI(i,j)+r(i,j));
      }
    } else if(Ider==Jder) {
      // taking the 2nd derivative of the 1D correlation function
      // the extra term is -2*theta(J)*r (remember r is for GEK so
      // the k loop part of it contains derivatives of the Kriging
      // correlation function with respect to XR)
      //      std::cout << "size(r)=[" << r.getNRows() << "," << r.getNCols() << "]\n"
      //		<< "size(drI)=[" << drI.getNRows() << "," << drI.getNCols() << "\n"
      //	<< "size(d2r)=[" << d2r.getNRows() << "," << d2r.getNCols()
      //	<< std::endl;
      for(j=0; j<nptsxr; ++j) {
	for(ipt=0, i=0; ipt<numWholePointsKeep; ++ipt, i+=neqn_per_pt) {
	  deltax=xr(Jder,j)-XRreorder(Jder,ipt);
	  d2r(i,j)=neg_two_thetaJ*(deltax*drI(i,j)+r(i,j));
	  for(k=0; k<numVarsr; ++k) {
	    //std::cout << "i=" << i << " j=" << j << " k=" << k << std::endl;
	    d2r(i+1+k,j)=neg_two_thetaJ*(deltax*drI(i+1+k,j)+r(i+1+k,j));
	  }
	  d2r(i+1+Jder,j)-=neg_two_thetaJ*drI(i,j); //minus a negative is
	  //a positive is correct, this extra term is for Jder=k
	}
	if(numPointsKeep>numWholePointsKeep) {
	  //ipt and i should already have the values we need them to
#ifdef __KRIG_ERR_CHECK__
	  assert((ipt==numWholePointsKeep)&&
		 (ipt==numPointsKeep-1)&&
		 (i==neqn_per_pt*numWholePointsKeep));
#endif
	  deltax=xr(Jder,j)-XRreorder(Jder,ipt);
	  d2r(i,j)=neg_two_thetaJ*(deltax*drI(i,j)+r(i,j));
	  for(k=0; k<numExtraDerKeep; ++k)
	    d2r(i+1+k,j)=neg_two_thetaJ*(deltax*drI(i+1+k,j)+r(i+1+k,j));
	  if(Jder<numExtraDerKeep)
	    d2r(i+1+Jder,j)-=neg_two_thetaJ*drI(i,j); //minus a negative is
  	    //a positive is correct, this extra term is for Jder=k
	}
      }
    } else {
      // taking the product of the 1st derivative of 2 independent 1D
      // correlation functions, (actually because this is for GEK, the
      // k loop is 2nd derivative of the Kriging r, in dimensions
      // independent of the one we're now taking the 1st derivative of)
      for(j=0; j<nptsxr; ++j) {
	for(ipt=0, i=0; ipt<numWholePointsKeep; ++ipt, i+=neqn_per_pt) {
	  deltax=xr(Jder,j)-XRreorder(Jder,ipt);
	  d2r(i,j)=neg_two_thetaJ*deltax*drI(i,j);
	  for(k=0; k<numVarsr; ++k)
	    d2r(i+1+k,j)=neg_two_thetaJ*deltax*drI(i+1+k,j);
	  d2r(i+1+Jder,j)-=neg_two_thetaJ*drI(i,j); //actually one element
	  //of the k loop is the dimension we're taking a derivative with
	  //respect to, it gets an extra term added to it. minus a negative
	  //is a positive is correct, this extra term is for Jder=k
	}
	if(numPointsKeep>numWholePointsKeep) {
	  //ipt and i should already have the values we need them to
#ifdef __KRIG_ERR_CHECK__
	  assert((ipt==numWholePointsKeep)&&
		 (ipt==numPointsKeep-1)&&
		 (i==neqn_per_pt*numWholePointsKeep));
#endif
	  deltax=xr(Jder,j)-XRreorder(Jder,ipt);
	  d2r(i,j)=neg_two_thetaJ*deltax*drI(i,j);
	  for(k=0; k<numExtraDerKeep; ++k)
	    d2r(i+1+k,j)=neg_two_thetaJ*deltax*drI(i+1+k,j);
	  if(Jder<numExtraDerKeep)
	    d2r(i+1+Jder,j)-=neg_two_thetaJ*drI(i,j); //minus a negative is
 	    //a positive is correct, this extra term is for Jder=k
	}
      }
    }
  } else if(corrFunc==EXP_CORR_FUNC) {
    std::cerr << "The exponential correlation function is not a valid correlation function for gradient enhanced Kriging\n";
      assert(false);
  } else if(corrFunc==POW_EXP_CORR_FUNC) {
    std::cerr << "The powered exponential (with power < 2) correlation function is not a valid correlation function for gradient enhanced Kriging\n";
      assert(false);
  } else if((corrFunc==MATERN_CORR_FUNC)&&(maternCorrFuncNu==1.5)) {
    // *********************************************************************
    // The MATERN 3/2 CORRELATION FUNCTION
    //
    // the 1st derivative with respect to xr of the 1D correlation function
    // is defined
    //
    // The 2nd derivative with respect to xr of the 1D correlation function
    // *is undefined at xr==XR,
    // *is -theta^2*(1-theta*|xr-XR|)*exp(-theta*|xr-XR|) at xr=/=XR
    // *approaches -theta^2 from above and below
    // when xr==XR we use the limit, -theta^2
    // this follows the approach of
    //   Lockwood, Brian A. and Anitescu, Mihai, "Gradient-Enhanced
    //      Universal Kriging for Uncertainty Proagation"
    //      Preprint ANL/MCS-P1808-1110
    // *********************************************************************
    double thetaJ=correlations(Jder,0);
    double thetaJ_squared=thetaJ*thetaJ;
    double thetaJ_abs_dx;
    if(numVarsr==1) {
      // if there is only one input variable we are taking the 2nd derivative
      // of the 1D GEK correlation function (which contains first derivatives
      // of the Kriging r with respect to XR) AND WE KNOW THAT Ider=Jder=k so
      // we don't have to "if" to known when to give the 2nd derivative of
      // GEK r = 3rd derivative of Kriging r, special treatment
      double thetaJ_cubed=thetaJ_squared*thetaJ;
      double r_div_matern_coef;
      for(j=0; j<nptsxr; ++j) {
	for(ipt=0, i=0; ipt<numWholePointsKeep; ++ipt, i+=2) {
	  deltax=xr(Jder,j)-XRreorder(Jder,ipt);
	  thetaJ_abs_dx=thetaJ*std::fabs(deltax);
	  r_div_matern_coef=r(i,j)/(1.0+thetaJ_abs_dx);
	  d2r(i  ,j)=r_div_matern_coef*thetaJ_squared*(thetaJ_abs_dx-1.0);
	  d2r(i+1,j)=r_div_matern_coef*thetaJ_cubed*dsign(deltax)*
	    (thetaJ_abs_dx-2.0);
	}
	// since there is only one dimension if there is a partial point
	// it will be a function value only, and actually recalculating it
	// will likely be faster on average then checking if there's a
	// partial point and calculating it if needed
	ipt=numPointsKeep-1;
	i=ipt*2;
	thetaJ_abs_dx=thetaJ*std::fabs(xr(Jder,j)-XRreorder(Jder,ipt));
	d2r(i,j)=r(i,j)/(1.0+thetaJ_abs_dx)*thetaJ_squared*(thetaJ_abs_dx-1.0);
      }
    } else if(Ider==Jder) {
      // taking the 2nd derivative of the 1D correlation function of the GEK
      // (not Kriging) r, which itself contains derivative of the Kriging r
      // with respect to XR, but this 2nd derivative is indepedent of those
      // first derivatives in all but one dimension
      double thetaJ_cubed=thetaJ_squared*thetaJ;
      double matern_coef;
      double d2_mult_r;
      for(j=0; j<nptsxr; ++j) {
	for(ipt=0, i=0; ipt<numWholePointsKeep; ++ipt, i+=neqn_per_pt) {
	  deltax=xr(Jder,j)-XRreorder(Jder,ipt);
	  thetaJ_abs_dx=thetaJ*std::fabs(deltax);
	  matern_coef=(1.0+thetaJ_abs_dx);
	  d2_mult_r=thetaJ_squared*(thetaJ_abs_dx-1.0)/matern_coef;
	  d2r(i,j)=r(i,j)*d2_mult_r;
	  for(k=0; k<numVarsr; ++k) //this k loop assumes that the current
	    //dimension is independent of the one that the XR derivative was
	    //taken with respect to, it's correct for all but one k
	    d2r(i+1+k,j)=r(i+1+k,j)*d2_mult_r;
	  //rather than having an if inside the loop which is slow, we're just
	  //going to reassign the d2r for k==Jder like this
	  d2r(i+1+Jder,j)=r(i,j)* //indexes of r(i,j) are correct
	    (thetaJ_cubed*dsign(deltax)*(thetaJ_abs_dx-2.0)/matern_coef);
	}

	if(numPointsKeep>numWholePointsKeep) {
	  //ipt and i should already have the values we need them to
#ifdef __KRIG_ERR_CHECK__
	  assert((ipt==numWholePointsKeep)&&
		 (ipt==numPointsKeep-1)&&
		 (i==neqn_per_pt*numWholePointsKeep));
#endif
	  deltax=xr(Jder,j)-XRreorder(Jder,ipt);
	  thetaJ_abs_dx=thetaJ*std::fabs(deltax);
	  matern_coef=(1.0+thetaJ_abs_dx);
	  d2_mult_r=thetaJ_squared*(thetaJ_abs_dx-1.0)/matern_coef;
	  d2r(i,j)=r(i,j)*d2_mult_r;
	  for(k=0; k<numExtraDerKeep; ++k) //this k loop assumes that the
	    //current dimension is independent of the one that the XR
	    //derivative was taken with respect to
	    d2r(i+1+k,j)=r(i+1+k,j)*d2_mult_r;
	  if(Jder<numExtraDerKeep) //if the dimension we're now taking a
	    //derivative with respect to wasn't clipped from the partial point
	    //we need to correct/reassign it for k==Jder
	    d2r(i+1+Jder,j)=r(i,j)* //indexes of r(i,j) are correct
	      (thetaJ_cubed*dsign(deltax)*(thetaJ_abs_dx-2.0)/matern_coef);
	}
      }
    } else {
      // taking the product of the 1st derivative (for GEK) of 2 independent
      // 1D correlation functions (they're independent because Jder!=Ider).
      // But since the GEK r contains first derivatives of the Kriging r,
      // there is one dimension, k==Jder, that needs special treatment
      double matern_coef;
      double d1_mult_r;
      for(j=0; j<nptsxr; ++j) {
	for(ipt=0, i=0; ipt<numWholePointsKeep; ++ipt, i+=neqn_per_pt) {
	  deltax=xr(Jder,j)-XRreorder(Jder,ipt);
	  thetaJ_abs_dx=thetaJ*std::fabs(deltax);
	  matern_coef=1.0+thetaJ_abs_dx;
	  d1_mult_r=-thetaJ_squared*deltax/matern_coef;
	  d2r(i,j)=drI(i,j)*d1_mult_r;
	  for(k=0; k<numVarsr; ++k)  //this k loop assumes that the
	    //current dimension is independent of the one that the XR
	    //derivative was taken with respect to
	    d2r(i+1+k,j)=drI(i+1+k,j)*d1_mult_r;
	  //rather than having an if inside the loop which is slow, we're just
	  //going to reassign the d2r for k==Jder like this
	  d2r(i+1+Jder,j)=drI(i,j)* //indexes of drI(i,j) are correct
	    thetaJ_squared*(1.0-thetaJ_abs_dx)/matern_coef; //sign is
	    //opposite the numVarsr==1 d2r(i,j) because one of the 2
	    //derivatives is taken with respect to XR instead of xr
	}
	if(numPointsKeep>numWholePointsKeep) {
	  //ipt and i should already have the values we need them to
#ifdef __KRIG_ERR_CHECK__
	  assert((ipt==numWholePointsKeep)&&
		 (ipt==numPointsKeep-1)&&
		 (i==neqn_per_pt*numWholePointsKeep));
#endif
	  deltax=xr(Jder,j)-XRreorder(Jder,ipt);
	  thetaJ_abs_dx=thetaJ*std::fabs(deltax);
	  matern_coef=1.0+thetaJ_abs_dx;
	  d1_mult_r=-thetaJ_squared*deltax/matern_coef;
	  d2r(i,j)=drI(i,j)*d1_mult_r;
	  for(k=0; k<numExtraDerKeep; ++k) //this k loop assumes that the
	    //current dimension is independent of the one that the XR
	    //derivative was taken with respect to
	    d2r(i+1+k,j)=drI(i+1+k,j)*d1_mult_r;
	  if(Jder<numExtraDerKeep)  //if the dimension we're now taking a
	    //derivative with respect to wasn't clipped from the partial point
	    //we need to correct/reassign it for k==Jder
	    d2r(i+1+Jder,j)=drI(i,j)* //indexes of drI(i,j) are correct
	      thetaJ_squared*(1.0-thetaJ_abs_dx)/matern_coef;
	}
      }
    }
  } else if((corrFunc==MATERN_CORR_FUNC)&&(maternCorrFuncNu==2.5)) {
    // *********************************************************************
    // The MATERN 5/2 CORRELATION FUNCTION
    //
    // the 1st and 2nd derivatives with respect to xr of the 1D correlation
    // function are defined,
    // 3rd derivative of Kriging r technically not defined at xr==XR (it is
    // defined everywhere else) but the limit from both sides is defined and
    // goes to zero at xr==XR (which means the limit from both sides agree)
    // so we'll use the else where defined 3rd derivative even at xr==XR
    // *********************************************************************
    double thetaJ=correlations(Jder,0);
    double thetaJ_squared=thetaJ*thetaJ;
    double thetaJ_abs_dx;
    if(numVarsr==1) {
      // if there is only one input variable we are taking the 2nd derivative
      // of the 1D GEK correlation function (which contains first derivatives
      // of the Kriging r with respect to XR) AND WE KNOW THAT Ider=Jder=k so
      // we don't have to "if" to known when to give the 2nd derivative of
      // GEK r = 3rd derivative of Kriging r, special treatment
      double r_thetaJ_squared_div_3_matern_coef;
      for(j=0; j<nptsxr; ++j) {
	for(ipt=0, i=0; ipt<numWholePointsKeep; ++ipt, i+=2) {
	  deltax=xr(Jder,j)-XRreorder(Jder,ipt);
	  thetaJ_abs_dx=thetaJ*std::fabs(deltax);
	  r_thetaJ_squared_div_3_matern_coef=r(i,j)*thetaJ_squared/
	    (3.0*(1.0+thetaJ_abs_dx)+thetaJ_abs_dx*thetaJ_abs_dx);
	  d2r(i  ,j)=r_thetaJ_squared_div_3_matern_coef*
	    -(1.0+thetaJ_abs_dx-thetaJ_abs_dx*thetaJ_abs_dx);
	  d2r(i+1,j)=r_thetaJ_squared_div_3_matern_coef*
	    -thetaJ_squared*deltax*(3.0-thetaJ_abs_dx);
	}
	// since there is only one dimension, if there is a partial point
	// it will be a function value only, and actually recalculating it
	// will likely be faster on average then checking if there's a
	// partial point and calculating it if needed
	ipt=numPointsKeep-1;
	i=ipt*2;
	thetaJ_abs_dx=thetaJ*std::fabs(xr(Jder,j)-XRreorder(Jder,ipt));
	d2r(i,j)=r(i,j)*thetaJ_squared/
	  (3.0*(1.0+thetaJ_abs_dx)+thetaJ_abs_dx*thetaJ_abs_dx)*
	  -(1.0+thetaJ_abs_dx-thetaJ_abs_dx*thetaJ_abs_dx);
      }
    } else if(Ider==Jder) {
      // taking the 2nd derivative of the 1D correlation function of the GEK
      // (not Kriging) r, which itself contains derivative of the Kriging r
      // with respect to XR, but this 2nd derivative is indepedent of those
      // first derivatives in all but one dimension
      double neg_thetaJ_squared_div_3_matern_coef;
      double d2_mult_r;
      for(j=0; j<nptsxr; ++j) {
	for(ipt=0, i=0; ipt<numWholePointsKeep; ++ipt, i+=neqn_per_pt) {
	  deltax=xr(Jder,j)-XRreorder(Jder,ipt);
	  thetaJ_abs_dx=thetaJ*std::fabs(deltax);
	  neg_thetaJ_squared_div_3_matern_coef=-thetaJ_squared/
	    (3.0*(1.0+thetaJ_abs_dx)+thetaJ_abs_dx*thetaJ_abs_dx);
	  d2_mult_r=neg_thetaJ_squared_div_3_matern_coef*
	    (1.0+thetaJ_abs_dx-thetaJ_abs_dx*thetaJ_abs_dx);
	  d2r(i,j)=r(i,j)*d2_mult_r;
	  for(k=0; k<numVarsr; ++k) //this k loop assumes that the current
	    //dimension is independent of the one that the XR derivative was
	    //taken with respect to, it's correct for all but one k
	    d2r(i+1+k,j)=r(i+1+k,j)*d2_mult_r;
	  //rather than having an if inside the loop which is slow, we're just
	  //going to reassign the d2r for k==Jder like this
	  d2r(i+1+Jder,j)=r(i,j)* //indexes of r(i,j) are correct
	    neg_thetaJ_squared_div_3_matern_coef*
	    thetaJ_squared*deltax*(3.0-thetaJ_abs_dx);
	}

	if(numPointsKeep>numWholePointsKeep) {
	  //ipt and i should already have the values we need them to
#ifdef __KRIG_ERR_CHECK__
	  assert((ipt==numWholePointsKeep)&&
		 (ipt==numPointsKeep-1)&&
		 (i==neqn_per_pt*numWholePointsKeep));
#endif
	  deltax=xr(Jder,j)-XRreorder(Jder,ipt);
	  thetaJ_abs_dx=thetaJ*std::fabs(deltax);
	  neg_thetaJ_squared_div_3_matern_coef=-thetaJ_squared/
	    (3.0*(1.0+thetaJ_abs_dx)+thetaJ_abs_dx*thetaJ_abs_dx);
	  d2_mult_r=neg_thetaJ_squared_div_3_matern_coef*
	    (1.0+thetaJ_abs_dx-thetaJ_abs_dx*thetaJ_abs_dx);
	  d2r(i,j)=r(i,j)*d2_mult_r;
	  for(k=0; k<numExtraDerKeep; ++k) //this k loop assumes that the
	    //current dimension is independent of the one that the XR
	    //derivative was taken with respect to
	    d2r(i+1+k,j)=r(i+1+k,j)*d2_mult_r;
	  if(Jder<numExtraDerKeep) //if the dimension we're now taking a
	    //derivative with respect to wasn't clipped from the partial point
	    //we need to correct/reassign it for k==Jder
	    d2r(i+1+Jder,j)=r(i,j)* //indexes of r(i,j) are correct
	      neg_thetaJ_squared_div_3_matern_coef*
	      thetaJ_squared*deltax*(3.0-thetaJ_abs_dx);
	}
      }
    } else {
      // taking the product of the 1st derivative (for GEK) of 2 independent
      // 1D correlation functions (they're independent because Jder!=Ider).
      // But since the GEK r contains first derivatives of the Kriging r,
      // there is one dimension, k==Jder, that needs special treatment
      double thetaJ_squared_div_3_matern_coef;
      double d1_mult_r;
      for(j=0; j<nptsxr; ++j) {
	for(ipt=0, i=0; ipt<numWholePointsKeep; ++ipt, i+=neqn_per_pt) {
	  deltax=xr(Jder,j)-XRreorder(Jder,ipt);
	  thetaJ_abs_dx=thetaJ*std::fabs(deltax);
	  thetaJ_squared_div_3_matern_coef=thetaJ_squared/
	    (3.0*(1.0+thetaJ_abs_dx)+thetaJ_abs_dx*thetaJ_abs_dx);
	  d1_mult_r=-thetaJ_squared_div_3_matern_coef*
	    deltax*(1.0+thetaJ_abs_dx);
	  d2r(i,j)=drI(i,j)*d1_mult_r;
	  for(k=0; k<numVarsr; ++k)  //this k loop assumes that the
	    //current dimension is independent of the one that the XR
	    //derivative was taken with respect to
	    d2r(i+1+k,j)=drI(i+1+k,j)*d1_mult_r;
	  //rather than having an if inside the loop which is slow, we're just
	  //going to reassign the d2r for k==Jder like this
	  d2r(i+1+Jder,j)=drI(i,j)* //indexes of drI(i,j) are correct
	    thetaJ_squared_div_3_matern_coef*
	    (1.0+thetaJ_abs_dx-thetaJ_abs_dx*thetaJ_abs_dx); //sign is
	    //opposite the numVarsr==1 d2r(i,j) because one of the 2
	    //derivatives is taken with respect to XR instead of xr
	}
	if(numPointsKeep>numWholePointsKeep) {
	  //ipt and i should already have the values we need them to
#ifdef __KRIG_ERR_CHECK__
	  assert((ipt==numWholePointsKeep)&&
		 (ipt==numPointsKeep-1)&&
		 (i==neqn_per_pt*numWholePointsKeep));
#endif
	  deltax=xr(Jder,j)-XRreorder(Jder,ipt);
	  thetaJ_abs_dx=thetaJ*std::fabs(deltax);
	  thetaJ_squared_div_3_matern_coef=thetaJ_squared/
	    (3.0*(1.0+thetaJ_abs_dx)+thetaJ_abs_dx*thetaJ_abs_dx);
	  d1_mult_r=-thetaJ_squared_div_3_matern_coef*
	    deltax*(1.0+thetaJ_abs_dx);
	  d2r(i,j)=drI(i,j)*d1_mult_r;
	  for(k=0; k<numExtraDerKeep; ++k) //this k loop assumes that the
	    //current dimension is independent of the one that the XR
	    //derivative was taken with respect to
	    d2r(i+1+k,j)=drI(i+1+k,j)*d1_mult_r;
	  if(Jder<numExtraDerKeep)  //if the dimension we're now taking a
	    //derivative with respect to wasn't clipped from the partial point
	    //we need to correct/reassign it for k==Jder
	    d2r(i+1+Jder,j)=drI(i,j)* //indexes of drI(i,j) are correct
	      thetaJ_squared_div_3_matern_coef*
	      (1.0+thetaJ_abs_dx-thetaJ_abs_dx*thetaJ_abs_dx);
	}
      }
    }
  } else{
    std::cerr << "unknown corrFunc in MtxDbl& KrigingModel::eval_gek_d2correlation_matrix_dxIdxJ(MtxDbl& d2r, const MtxDbl& drI, const MtxDbl& r, const MtxDbl& xr, int Ider, int Jder) const\n";
    assert(false);
  }
  return d2r;
}



/** this function is typically used during emulator construction, the below
    the diagonal portion of R = exp(Z^T*theta), where R is symmetric with 1's
    on the diagonal, theta is the vector of correlations and the Z matrix is
    defined as Z(k,ij)=-(XR(k,i)-XR(k,j))^2 where ij counts downward within
    columns of R starting from the element below the diagonal and continues
    from one column to the next, Z^T*theta is matrix vector multiplication to
    be performed efficiently by BLAS, V=Z^T*theta is a vector with
    nchoosek(numPoints,2) elements.  We need to copy exp(V(ij)) to R(i,j)
    and R(j,i) to produce R. The Z matrix is produced by
    KrigingModel::gen_Z_matrix()     KRD wrote this */
void KrigingModel::correlation_matrix(const MtxDbl& theta)
{
  int ncolsZ=Z.getNCols();
  //printf("nrowsZ=%d; numPoints=%d; ''half'' numPoints^2=%d; numVarsr=%d; theta.getNRows()=%d\n",
  //	 ncolsZ,numPoints,nchoosek(numPoints,2),numVarsr,theta.getNRows());
  //fflush(stdout);
#ifdef __KRIG_ERR_CHECK__
  assert((ncolsZ==nchoosek(numPoints,2))&&
	 (numVarsr==Z.getNRows())&&
	 (numVarsr==theta.getNRows())&&
	 (1==theta.getNCols()));
#endif

  Ztran_theta.newSize(ncolsZ,1); //Z transpose because subsequent access of a
  //column vector should be marginally faster than a row vector
  matrix_mult(Ztran_theta,Z,theta,0.0,1.0,'T','N');

  if(buildDerOrder==0)
    numRowsR=numPoints;
  else if(buildDerOrder==1)
    numRowsR=numPoints*nDer;
  else{
    std::cerr << "buildDerOrder=" << buildDerOrder << " in void KrigingModel::correlation_matrix(const MtxDbl& theta).  It must either be 0 for Kriging or 1 for Gradient Enhanced Kriging.  Higher order build derivatives, (e.g. Hessian Enhanced Kriging) have not been implemented." << std::endl;
    assert(false);
  }
  R.newSize(numRowsR,numRowsR);

  //Do the regular (Der0) Kriging Portion of the Correlation matrix first
  double Rij_temp;
  int ij=0;
  if((corrFunc==GAUSSIAN_CORR_FUNC)||
     (corrFunc==EXP_CORR_FUNC)||
     (corrFunc==POW_EXP_CORR_FUNC)) {
    for(int j=0; j<numPoints-1; ++j) {
      R(j,j)=1.0;
      for(int i=j+1; i<numPoints; ++i, ++ij) {
	Rij_temp=std::exp(Ztran_theta(ij,0));
	R(i,j)=Rij_temp;
	R(j,i)=Rij_temp;
      }
    }
  } else if((corrFunc==MATERN_CORR_FUNC)&&(maternCorrFuncNu==1.5)){
    //for matern Z(k,ij)=-|XR(k,i)-XR(k,j)| we want to feed
    //theta(k,0)*|XR(k,i)-XR(k,j)| to matern_1pt5_coef so we need to
    //negate the already negative quantity
    if(numVarsr==1)
      for(int j=0; j<numPoints-1; ++j) {
	R(j,j)=1.0;
	for(int i=j+1; i<numPoints; ++i, ++ij) {
	  Rij_temp=std::exp(Ztran_theta(ij,0))*
	    matern_1pt5_coef(-Ztran_theta(ij,0));
	  R(i,j)=Rij_temp;
	  R(j,i)=Rij_temp;
	}
      }
    else
      for(int j=0; j<numPoints-1; ++j) {
	R(j,j)=1.0;
	for(int i=j+1; i<numPoints; ++i, ++ij) {
	  Rij_temp=std::exp(Ztran_theta(ij,0))*
	    matern_1pt5_coef(-Z(0,ij)*theta(0,0));
	  for(int k=1; k<numVarsr; ++k)
	    Rij_temp*=matern_1pt5_coef(-Z(k,ij)*theta(k,0));
	  R(i,j)=Rij_temp;
	  R(j,i)=Rij_temp;
	}
      }
  } else if((corrFunc==MATERN_CORR_FUNC)&&(maternCorrFuncNu==2.5)){
    //for matern Z(k,ij)=-|XR(k,i)-XR(k,j)| we want to feed
    //theta(k,0)*|XR(k,i)-XR(k,j)| to matern_2pt5_coef so we need to
    //negate the already negative quantity
    if(numVarsr==1)
      for(int j=0; j<numPoints-1; ++j) {
	R(j,j)=1.0;
	for(int i=j+1; i<numPoints; ++i, ++ij) {
	  Rij_temp=std::exp(Ztran_theta(ij,0))*
	    matern_2pt5_coef(-Ztran_theta(ij,0));
	  R(i,j)=Rij_temp;
	  R(j,i)=Rij_temp;
	}
      }
    else
      for(int j=0; j<numPoints-1; ++j) {
	R(j,j)=1.0;
	for(int i=j+1; i<numPoints; ++i, ++ij) {
	  Rij_temp=std::exp(Ztran_theta(ij,0))*
	    matern_2pt5_coef(-Z(0,ij)*theta(0,0));
	  for(int k=1; k<numVarsr; ++k)
	    Rij_temp*=matern_2pt5_coef(-Z(k,ij)*theta(k,0));
	  R(i,j)=Rij_temp;
	  R(j,i)=Rij_temp;
	}
      }
  }else{
    std::cerr << "unknown corrFunc in void KrigingModel::correlation_matrix(const MtxDbl& theta)\n";
    assert(false);
  }
  R(numPoints-1,numPoints-1)=1.0;

  /*
  FILE *fp=fopen("km_Rmat_check.txt","w");
  for(int i=0; i<numPoints; ++i) {
      fprintf(fp,"%-12.6g", R(i,0));
      for(int j=1; j<numPoints; ++j)
	fprintf(fp," %-12.6g", R(i,j));
      fprintf(fp,"\n");
  }
  fclose(fp);
  */

  if(buildDerOrder>0) {
    //Gaussian Matern1.5 and Matern2.5 are valid correlation functions for
    //Gradient Enhanced Kriging
    double temp_double;
    int zij, j;


    if(corrFunc==GAUSSIAN_CORR_FUNC) {
      //now handle the first order derivative submatrices, indiviually the first
      //order derivative SUBmatrices are anti-symmetric but the whole matrix is
      //symmetric
      int Ii, Ij, Jj, Ji; //first letter identifies index OF derivative submatrix
      //second letter identifies index INTO derivative SUBmatrix
      for(int Ider=0; Ider<numVarsr; ++Ider) {
	zij=0;
	double two_theta_Ider=2.0*theta(Ider,0);
	for(j=0; j<numPoints-1; ++j) {//j<numPoints-1 avoids an i loop of length 0
	  //diagonal (_j,_j) of off diagonal (I_, _) submatrix
	  Ij=(Ider+1)*numPoints+j;
	  R(Ij, j)=0.0;
	  R( j,Ij)=0.0;
	  //Ij=(Ider+1)*numPoints+j;
	  for(int i=j+1; i<numPoints; ++i, ++zij) {
	    //off diagonal (_i,_j) of off-diagonal (I_, _) submatrix
	    Ii=(Ider+1)*numPoints+i;
	    temp_double=-two_theta_Ider*deltaXR(zij,Ider)*R( i, j);
	    //here  temp_double=
	    //                  R(Ii, j) = dR(i,j)/dXR1(Ider,i)
	    //                  R( j,Ii) = dR(i,j)/dXR2(Ider,j)
	    //and  -temp_double=
	    //                  R(Ij, i) = dR(i,j)/dXR1(Ider,j)
	    //                  R( i,Ij) = dR(i,j)/dXR2(Ider,i)
	    //where XR1 is the first argument of the correlation function
	    //and XR2 is the second argument of the correlation function
	    R(Ii, j)= temp_double;
	    R( j,Ii)= temp_double; //whole R matrix is symmetric
	    R(Ij, i)=-temp_double;
	    R( i,Ij)=-temp_double; //off-diagonal 1st order (actually all odd
	    //order) derivative SUBmatrices are anti-symmetric
	  }
	}
	//diagonal (_j,_j) of off diagonal (I_, _) submatrix
	j=numPoints-1; //avoids an i loop of length 0
	Ij=(Ider+1)*numPoints+j;
	R(Ij,j)=0.0;
	R(j,Ij)=0.0;
      }

      //note that all 2nd order (actually all even order) derivative SUBmatrices
      //are symmetric because the hadamard product of 2 (actually any even
      //number of) anti-symmetric matrices is a symmetric matrix
      double two_theta_Jder;
      for(int Jder=0; Jder<numVarsr; ++Jder) {
	//do the on diagonal (J_,J_) submatrix
	two_theta_Jder=2.0*theta(Jder,0);
	zij=0;
	for(j=0; j<numPoints-1; ++j) { //j<numPoints-1 avoids an i loop of length 0
	  //diagonal (_j,_j) of on diagonal (J_,J_) submatrix
	  Jj=(Jder+1)*numPoints+j;
	  R(Jj,Jj)=two_theta_Jder; //R(Jj,Jj)=2*theta(Jder,0)*R(j,j); R(j,j)=1;
	  for(int i=j+1; i<numPoints; ++i) {
	    //off diagonal (_i,_j) of on-diagonal (J_,J_) submatrix
	    Ji=(Jder+1)*numPoints+i;
	    temp_double=two_theta_Jder*deltaXR(zij,Jder)*R(Ji, j)+
	    two_theta_Jder*R( i, j);
	    R(Ji,Jj)=temp_double;
	    R(Jj,Ji)=temp_double;
	    ++zij;
	  }
	}
	//diagonal (_j,_j) of on diagonal (J_,J_) submatrix
	j=numPoints-1; //avoids an i loop of length 0
	Jj=(Jder+1)*numPoints+j;
	R(Jj,Jj)=two_theta_Jder; //R(j,j)=1 R(Jj,Jj)=2*theta(Jder,0)*R(j,j)


	//do the off diagonal (I_,J_) submatrices
	for(int Ider=Jder+1; Ider<numVarsr; ++Ider) {
	  //off diagonal (I_,J_) submatrix
	  zij=0;
	  for(j=0; j<numPoints-1; ++j) {//j<numPoints-1 avoids an i loop of length 0
	    //diagonal (_j,_j) of off-diagonal (I_,J_) submatrix
	    Jj=(Jder+1)*numPoints+j;
	    Ij=(Ider+1)*numPoints+j;
	    R(Ij,Jj)=0.0;
	    R(Jj,Ij)=0.0;


	    for(int i=j+1; i<numPoints; ++i) {
	      //off diagonal (_i,_j) of off-diagonal (I_,J_) submatrix
	      Ii=(Ider+1)*numPoints+i;
	      Ji=(Jder+1)*numPoints+i;
	      temp_double=two_theta_Jder*deltaXR(zij,Jder)*R(Ii, j);
	      R(Ii,Jj)= temp_double;
	      R(Ij,Ji)= temp_double;
	      R(Ji,Ij)= temp_double;
	      R(Jj,Ii)= temp_double;
	      ++zij;
	    }
	  }
	  //diagonal (_j,_j) of off-diagonal (I_,J_) submatrix
	  j=numPoints-1; //avoids an i loop of length 0
	  Ij=(Ider+1)*numPoints+j;
	  Jj=(Jder+1)*numPoints+j;
	  R(Ij,Jj)=0.0;
	  R(Jj,Ij)=0.0;
	}
      }
    } else if((corrFunc==MATERN_CORR_FUNC)&&(maternCorrFuncNu==1.5)) {
      //The second derivative of the Matern1.5 correlation function
      //is not strictly defined at XR(i,Ider)==XR(j,Jder) but the limit
      //of the second derivative from both sides is defined and is the same
      //this follows
      //Lockwood, Brian A. and Anitescu, Mihai, "Gradient-Enhanced
      //    Universal Kriging for Uncertainty Proagation"
      //    Preprint ANL/MCS-P1808-1110
      //
      //d2r_dXIdXJ with Ider==Jder
      // = theta^2*exp(-theta*|XI-XJ|)-theta^3*|XI-XJ|*exp(-theta*|XI-XJ|)
      // = -theta^2*(1-2/matern_1pt5_coef)*r(XI,XJ)
      // = -matern_1pt5_d2_mult_r(theta,+/-(XI-XJ))*r(XI,XJ) (note the
      //    negative sign, it should be here, but does not appear when
      //    evalutation 2nd derivative of GP, because there it is second
      //    derivative with respect to the SAME argument, here it is the
      //    second derivative with respect to different arguments)

      //now handle the first order derivative submatrices, indiviually the first
      //order derivative SUBmatrices are anti-symmetric but the whole matrix is
      //symmetric
      int Ii, Ij, Jj, Ji; //first letter identifies index OF derivative
      //submatrix second letter identifies index INTO derivative SUBmatrix
      for(int Ider=0; Ider<numVarsr; ++Ider) {
	zij=0;
	double theta_Ider=theta(Ider,0);
	for(j=0; j<numPoints-1; ++j) {//j<numPoints-1 avoids an i loop of
	  //length 0

	  //diagonal (_j,_j) of off diagonal (I_, _) submatrix
	  Ij=(Ider+1)*numPoints+j;
	  R(Ij, j)=0.0;
	  R( j,Ij)=0.0;
	  //Ij=(Ider+1)*numPoints+j;
	  for(int i=j+1; i<numPoints; ++i) {
	    //off diagonal (_i,_j) of off-diagonal (I_, _) submatrix
	    Ii=(Ider+1)*numPoints+i;
	    temp_double=
	      matern_1pt5_d1_mult_r(theta_Ider,deltaXR(zij,Ider))*R( i, j);
	    R(Ii, j)= temp_double;
	    R( j,Ii)= temp_double; //whole R matrix is symmetric
	    R(Ij, i)=-temp_double;
	    R( i,Ij)=-temp_double; //off-diagonal 1st order (actually all odd
	    //order) derivative SUBmatrices are anti-symmetric
	    ++zij;
	  }
	}
	//diagonal (_j,_j) of off diagonal (I_, _) submatrix
	j=numPoints-1; //avoids an i loop of length 0
	Ij=(Ider+1)*numPoints+j;
	R(Ij, j)=0.0;
	R( j,Ij)=0.0;
      }

      //note that all 2nd order (actually all even order) derivative SUBmatrices
      //are symmetric because the hadamard product of 2 (actually any even
      //number of) anti-symmetric matrices is a symmetric matrix
      double theta_Jder;
      double theta_Jder_squared;
      for(int Jder=0; Jder<numVarsr; ++Jder) {
	//do the on diagonal (J_,J_) submatrix
	theta_Jder=theta(Jder,0);
	theta_Jder_squared=theta_Jder*theta_Jder;
	zij=0;
	for(j=0; j<numPoints-1; ++j) { //j<numPoints-1 avoids an i loop of length 0
	  //diagonal (_j,_j) of on diagonal (J_,J_) submatrix
	  Jj=(Jder+1)*numPoints+j;
	  R(Jj,Jj)=theta_Jder_squared;
	  for(int i=j+1; i<numPoints; ++i) {
	    //off diagonal (_i,_j) of on-diagonal (J_,J_) submatrix
	    Ji=(Jder+1)*numPoints+i;
	    temp_double=//neg sign because d^2/dXR1dXR2 instead of d^2/dXR1^2
	      -matern_1pt5_d2_mult_r(theta_Jder,deltaXR(zij,Jder))*R( i, j);
	    R(Ji,Jj)=temp_double;
	    R(Jj,Ji)=temp_double;
	    ++zij;
	  }
	}
	//diagonal (_j,_j) of on diagonal (J_,J_) submatrix
	j=numPoints-1; //avoids an i loop of length 0
	Jj=(Jder+1)*numPoints+j;
	R(Jj,Jj)=theta_Jder_squared;


	//do the off diagonal (I_,J_) submatrices
	for(int Ider=Jder+1; Ider<numVarsr; ++Ider) {
	  //off diagonal (I_,J_) submatrix
	  zij=0;
	  for(j=0; j<numPoints-1; ++j) {//j<numPoints-1 avoids an i loop of length 0
	    //diagonal (_j,_j) of off-diagonal (I_,J_) submatrix
	    Jj=(Jder+1)*numPoints+j;
	    Ij=(Ider+1)*numPoints+j;
	    R(Ij,Jj)=0.0;
	    R(Jj,Ij)=0.0;

	    for(int i=j+1; i<numPoints; ++i) {
	      //off diagonal (_i,_j) of off-diagonal (I_,J_) submatrix
	      Ii=(Ider+1)*numPoints+i;
	      Ji=(Jder+1)*numPoints+i;
	      temp_double=
		matern_1pt5_d1_mult_r(theta_Jder,-deltaXR(zij,Jder))*R(Ii, j);
	      R(Ii,Jj)= temp_double;
	      R(Ij,Ji)= temp_double;
	      R(Ji,Ij)= temp_double;
	      R(Jj,Ii)= temp_double;
	      ++zij;
	    }
	  }
	  //diagonal (_j,_j) of off-diagonal (I_,J_) submatrix
	  j=numPoints-1; //avoids an i loop of length 0
	  Ij=(Ider+1)*numPoints+j;
	  Jj=(Jder+1)*numPoints+j;
	  R(Ij,Jj)=0.0;
	  R(Jj,Ij)=0.0;
	}
      }

    } else if((corrFunc==MATERN_CORR_FUNC)&&(maternCorrFuncNu==2.5)) {
      //now handle the first order derivative submatrices, indiviually the first
      //order derivative SUBmatrices are anti-symmetric but the whole matrix is
      //symmetric
      int Ii, Ij, Jj, Ji; //first letter identifies index OF derivative
      //submatrix second letter identifies index INTO derivative SUBmatrix
      for(int Ider=0; Ider<numVarsr; ++Ider) {
	zij=0;
	double theta_Ider=theta(Ider,0);
	for(j=0; j<numPoints-1; ++j) {//j<numPoints-1 avoids an i loop of
	  //length 0

	  //diagonal (_j,_j) of off diagonal (I_, _) submatrix
	  Ij=(Ider+1)*numPoints+j;
	  R(Ij, j)=0.0;
	  R( j,Ij)=0.0;
	  //Ij=(Ider+1)*numPoints+j;
	  for(int i=j+1; i<numPoints; ++i) {
	    //off diagonal (_i,_j) of off-diagonal (I_, _) submatrix
	    Ii=(Ider+1)*numPoints+i;
	    temp_double=
	      matern_2pt5_d1_mult_r(theta_Ider,deltaXR(zij,Ider))*R( i, j);
	    R(Ii, j)= temp_double;
	    R( j,Ii)= temp_double; //whole R matrix is symmetric
	    R(Ij, i)=-temp_double;
	    R( i,Ij)=-temp_double; //off-diagonal 1st order (actually all odd
	    //order) derivative SUBmatrices are anti-symmetric
	    ++zij;
	  }
	}
	//diagonal (_j,_j) of off diagonal (I_, _) submatrix
	j=numPoints-1; //avoids an i loop of length 0
	Ij=(Ider+1)*numPoints+j;
	R(Ij, j)=0.0;
	R( j,Ij)=0.0;
      }

      //note that all 2nd order (actually all even order) derivative SUBmatrices
      //are symmetric because the hadamard product of 2 (actually any even
      //number of) anti-symmetric matrices is a symmetric matrix
      double theta_Jder;
      double theta_Jder_squared_div_3;
      for(int Jder=0; Jder<numVarsr; ++Jder) {
	//do the on diagonal (J_,J_) submatrix
	theta_Jder=theta(Jder,0);
	theta_Jder_squared_div_3=theta_Jder*theta_Jder/3.0;
	zij=0;
	for(j=0; j<numPoints-1; ++j) { //j<numPoints-1 avoids an i loop of length 0
	  //diagonal (_j,_j) of on diagonal (J_,J_) submatrix
	  Jj=(Jder+1)*numPoints+j;
	  R(Jj,Jj)=theta_Jder_squared_div_3;
	  for(int i=j+1; i<numPoints; ++i) {
	    //off diagonal (_i,_j) of on-diagonal (J_,J_) submatrix
	    Ji=(Jder+1)*numPoints+i;
	    temp_double=//neg sign because d^2/dXR1dXR2 instead of d^2/dXR1^2
	      -matern_2pt5_d2_mult_r(theta_Jder,deltaXR(zij,Jder))*R( i, j);
	    R(Ji,Jj)=temp_double;
	    R(Jj,Ji)=temp_double;
	    ++zij;
	  }
	}
	//diagonal (_j,_j) of on diagonal (J_,J_) submatrix
	j=numPoints-1; //avoids an i loop of length 0
	Jj=(Jder+1)*numPoints+j;
	R(Jj,Jj)=theta_Jder_squared_div_3;


      	//do the off diagonal (I_,J_) submatrices
	for(int Ider=Jder+1; Ider<numVarsr; ++Ider) {
	  //off diagonal (I_,J_) submatrix
	  zij=0;
	  for(j=0; j<numPoints-1; ++j) {//j<numPoints-1 avoids an i loop of length 0
	    //diagonal (_j,_j) of off-diagonal (I_,J_) submatrix
	    Jj=(Jder+1)*numPoints+j;
	    Ij=(Ider+1)*numPoints+j;
	    R(Ij,Jj)=0.0;
	    R(Jj,Ij)=0.0;

	    for(int i=j+1; i<numPoints; ++i) {
	      //off diagonal (_i,_j) of off-diagonal (I_,J_) submatrix
	      Ii=(Ider+1)*numPoints+i;
	      Ji=(Jder+1)*numPoints+i;
	      temp_double=
		matern_2pt5_d1_mult_r(theta_Jder,-deltaXR(zij,Jder))*R(Ii, j);
	      R(Ii,Jj)= temp_double;
	      R(Ij,Ji)= temp_double;
	      R(Ji,Ij)= temp_double;
	      R(Jj,Ii)= temp_double;
	      ++zij;
	    }
	  }
	  //diagonal (_j,_j) of off-diagonal (I_,J_) submatrix
	  j=numPoints-1; //avoids an i loop of length 0
	  Ij=(Ider+1)*numPoints+j;
	  Jj=(Jder+1)*numPoints+j;
	  R(Ij,Jj)=0.0;
	  R(Jj,Ij)=0.0;
	}
      }
    } else{
      std::cerr << "Unknown or Invalid Correlation function for Gradient Enhanced Kriging in void KrigingModel::correlation_matrix(const MtxDbl& theta)\n";
      assert(false);
    }
    /*
    printf("theta=[%14.8g",theta(0,0));
    for(int k=1; k<numVarsr; ++k)
      printf(" %14.8g",theta(k,0));
    printf("]^T\n");
    printf("M3/2 GEK R=\n");
    for(int i=0; i<numEqnAvail; ++i) {
      for(int j=0; j<numEqnAvail; ++j)
	printf("%14.8g ",R(i,j));
      printf("\n");
    }
    printf("\n\n");
    */

  }

  return;
}

/** the Z matrix is defined as Z(k,ij)=-(XR(i,k)-XR(j,k))^2 where
    ij=i+j*XR.getNRows(), it enables the efficient repeated calculation
    of the R matrix during model construction:
    R=reshape(exp(Z*theta),XR.getNRows(),XR.getNRows()) where theta is
    the vector of correlations and * is matrix vector multiplication,
    note that the Z matrix is independent of the correlation vector so
    it can be formed once and later during the search for a good
    correlation vector, the matrix vector product Z*theta can be
    performed efficiently (for each member of the set of candidate
    theta vectors) by calling BLAS. Z and XR are member variables so
    they don't need to be passed in, KRD wrote this,  */
MtxDbl& KrigingModel::gen_Z_matrix()
{
#ifdef __KRIG_ERR_CHECK__
  assert((XR.getNRows()==numVarsr)&&(XR.getNCols()==numPoints));
#endif
  int ncolsZ=nchoosek(numPoints,2);
  Z.newSize(numVarsr,ncolsZ);

  if(buildDerOrder>0) {
    //deltaXR is only needed for GEK
    deltaXR.newSize(ncolsZ,numVarsr); //this ordering (transpose of Z)
    //is useful for constructing the GEK R matrix
  }

  int ij=0;
  if(corrFunc==GAUSSIAN_CORR_FUNC)  {
    // ****************************************************************
    // The Gausssian Correlation Function
    // can be used for Gradient Enhanced Kriging (GEK)
    // (or more generally, for derivative enhanced Kriging) because
    // it is C infinity continuous.
    // It uses Z(k,ij) = -(XR(k,i)-XR(k,j))^2
    // if GEK is used we will also compute and store
    // deltaXR(ij,k) = XR(k,i)-XR(k,j), the order of indexes is correct
    // this transposed ordering is useful for efficient computation of
    // the GEK R matrix (using the SUBmatrix construction process)
    // ****************************************************************
    double dXR; //a temporary variable to make squaring easier
    if(buildDerOrder>0)
      for(int j=0; j<numPoints-1; ++j)
	for(int i=j+1; i<numPoints; ++i, ++ij)
	  for(int k=0; k<numVarsr; k++) {
	    dXR=XR(k,i)-XR(k,j);
	    deltaXR(ij,k)=dXR;
	    Z(k,ij)=-dXR*dXR;
	  }
    else
      for(int j=0; j<numPoints-1; ++j)
	for(int i=j+1; i<numPoints; ++i, ++ij)
	  for(int k=0; k<numVarsr; k++) {
	    dXR=XR(k,i)-XR(k,j);
	    Z(k,ij)=-dXR*dXR;
	  }
  } else if((corrFunc==EXP_CORR_FUNC)||(corrFunc==MATERN_CORR_FUNC)) {
    // ****************************************************************
    // The Exponential and Matern 3/2 and 5/2 Correlation Functions
    // all use Z(k,ij) = -|XR(k,i)-XR(k,j)|
    // the Exponential Correlation Function
    //     can NOT be used for Gradient Enhanced Kriging (GEK)
    // the Matern 3/2 and 5/2 Correlation Functions
    //     CAN be used for Gradient Enhanced Kriging (GEK)
    // if GEK is used we will also compute and store
    //     deltaXR(ij,k) = XR(k,i)-XR(k,j), the order of indexes is
    //     correct this transposed ordering is useful for efficient
    //     computation of the GEK R matrix (using the SUBmatrix
    //     construction process)
    // ****************************************************************
    if(buildDerOrder>0) {
      if(corrFunc==EXP_CORR_FUNC) {
	std::cerr << "the exponential correlation function is not a valid choice for Gradient Enhanced Kriging\n";
	assert(!((corrFunc==EXP_CORR_FUNC)&&(buildDerOrder>0)));
      }
      for(int j=0; j<numPoints-1; ++j)
	for(int i=j+1; i<numPoints; ++i, ++ij)
	  for(int k=0; k<numVarsr; k++) {
	    deltaXR(ij,k)=XR(k,i)-XR(k,j);
	    Z(k,ij)=-std::fabs(deltaXR(ij,k));
	  }
    } else
      for(int j=0; j<numPoints-1; ++j)
	for(int i=j+1; i<numPoints; ++i, ++ij)
	  for(int k=0; k<numVarsr; k++)
	    Z(k,ij)=-std::fabs(XR(k,i)-XR(k,j));
  } else if(corrFunc==POW_EXP_CORR_FUNC) {
    // ****************************************************************
    // The Powered Exponential Correlation Function
    // uses Z(k,ij) = -(|XR(k,i)-XR(k,j)|^powExpCorrFuncPow)
    // where 1.0<powExpCorrFuncPow<2.0
    // It can NOT be used for Gradient Enhanced Kriging (GEK)
    // ****************************************************************
    if(buildDerOrder>0) {
      std::cerr << "the powered exponential correlation function is not a valid choice for Gradient Enhanced Kriging\n";
      assert(!((corrFunc==POW_EXP_CORR_FUNC)&&(buildDerOrder>0)));
    }
    for(int j=0; j<numPoints-1; ++j)
      for(int i=j+1; i<numPoints; ++i, ++ij)
	for(int k=0; k<numVarsr; k++)
	  Z(k,ij)=-std::pow(std::fabs(XR(k,i)-XR(k,j)),powExpCorrFuncPow);
  } else{
    std::cerr << "unknown Correlation Function in MtxDbl& KrigingModel::gen_Z_matrix()\n";
    assert(false);
  }
  return Z;
}

void KrigingModel::reorderCopyRtoRChol() {
  numRowsR=numEqnAvail;
  RChol.newSize(numRowsR,numRowsR);

  if(buildDerOrder==0) {
    //Kriging
    for(int jpt=0; jpt<numPoints; ++jpt) {
      int jsrc=iPtsKeep(jpt,0);
      for(int ipt=0; ipt<numPoints; ++ipt)
	RChol(ipt,jpt)=R(iPtsKeep(ipt,0),jsrc);
    }
  } else if(buildDerOrder==1) {
    //Gradient Enhanced Kriging, R is blocked into (1+numVarsr) by (1+numVarsr)
    //submatrices.  Each submatrix has numPoints by numPoints elements, i.e.
    //the same size as the Kriging R matrix, in fact the upper-left-most
    //submatrix is the Kriging R matrix, we need to reorder this so that
    //"Whole Points" (a function value immediately followed by its gradient)
    //are listed in the order given in iPtsKeep

    int i, j;
    for(int jpt=0, j=0; jpt<numPoints; ++jpt)
      for(int jder=-1; jder<numVarsr; ++jder, ++j) {
	int jsrc=iPtsKeep(jpt,0)+(jder+1)*numPoints;
	for(int ipt=0, i=0; ipt<numPoints; ++ipt)
	  for(int ider=-1; ider<numVarsr; ++ider, ++i)
	    RChol(i,j)=R(iPtsKeep(ipt,0)+(ider+1)*numPoints,jsrc);
      }
  } else {
    std::cerr << "buildDerOrder=" << buildDerOrder
	      << " in void KrigingModel::reorderCopyRtoRChol(); "
	      << "for Kriging buildDerOrder must be 0; "
	      << "for Gradient Enhanced Kriging buildDerOrder must be 1; "
	      << "Higher order derivative enhanced Kriging "
	      << "(e.g Hessian Enhanced Kriging) has not been implemented"
	      << std::endl;
    assert(false);
  }
  return;
}

void KrigingModel::nuggetSelectingCholR(){
  if(buildDerOrder==0)
    numExtraDerKeep=0;
  else if(buildDerOrder==1)
    numExtraDerKeep=numVarsr; //the last point will have all of the gradient
  else{
    std::cerr << "buildDerOrder=" << buildDerOrder
	      << " in void KrigingModel::nuggetSelectingCholR(); "
	      << "for Kriging buildDerOrder must be 0; "
	      << "for Gradient Enhanced Kriging buildDerOrder must be 1; "
	      << "Higher order derivative enhanced Kriging "
	      << "(e.g Hessian Enhanced Kriging) has not been implemented"
	      << std::endl;
    assert(false);
  }
  numWholePointsKeep=numPointsKeep=numPoints;

  double min_allowed_rcond=1.0/maxCondNum;
  int ld_RChol=RChol.getNRowsAct();
  rcondDblWork.newSize(3*ld_RChol,1);
  rcondIntWork.newSize(ld_RChol,1);
  scaleRChol.newSize(numEqnAvail,1); //scaling/equilibrating is only
  //necessary if GEK is used (because Kriging already has all ones on
  //the diagonal of R; GEK doesn't) but the generic Cholesky
  //factorization won't know in advance whether it's needed or not
  //you can calculate rcond essentially "for free" if you do it at the
  //same time as the Cholesky factorization
  int chol_info;

  //point order is the default point order
  for(int ipt=0; ipt<numPointsKeep; ++ipt)
    iPtsKeep(ipt,0)=ipt;
  if(ifAssumeRcondZero==true)
    rcondR=0.0;
  else {
    //but if GEK is used I still need to reorder from derivative submatrix
    //blocks to whole point order
    reorderCopyRtoRChol();

    //See the end of the KrigingModel constructor for why Y and Gtran are
    //what we already need them to be.
    //the maximumAllowedPolyOrder given the number of Points is already
    //selected, and Gtran is already what we need it to be

    nug=0.0;
    Chol_fact_workspace(RChol,scaleRChol,rcondDblWork,rcondIntWork,
			chol_info,rcondR);
  }

  //this rcondR is for the equilibrated R/RChol (so pretend it has all
  //ones on the diagonal)
  if(rcondR<=min_allowed_rcond) {
    double dbl_num_eqn=static_cast<double>(numEqnAvail);
    double sqrt_num_eqn=std::sqrt(dbl_num_eqn);
    min_allowed_rcond*=sqrt_num_eqn; //one norm is within a factor of N^0.5
    //of 2 norm
    rcondR/=sqrt_num_eqn; //one norm is within a factor of N^0.5 of 2 norm
    double min_eig_worst=(rcondR*dbl_num_eqn)/(1.0+(dbl_num_eqn-1.0)*rcondR);
    double max_eig_worst=dbl_num_eqn-(dbl_num_eqn-1.0)*min_eig_worst;
    nug=(min_allowed_rcond*max_eig_worst-min_eig_worst)/
      (1.0-min_allowed_rcond);
    //this nugget will make the worst case scenario meet (with an ==)
    //the maxCondNum constraint, I (KRD) don't expect this to
    //ever == fail because I don't expect rcond to be *N^-0.5 without
    //nugget and be *N^0.5 with nugget while the maximum eigen value
    //of R (without nugget) is N-(N-1)*min_eigval (that comes from
    //assumming all eigenvalues except the largest are the smallest
    //possible for the given rcond) note that rcond is the LAPACK
    //ESTIMATE of the 1 norm condition number so there are no 100%
    //guarantees.
    apply_nugget_build(); //multiply the diagonal elements by (1.0+nug)
    reorderCopyRtoRChol();

    Chol_fact_workspace(RChol,scaleRChol,rcondDblWork,rcondIntWork,
			chol_info,rcondR);
  }
  return;
}


/* use Pivoted Cholesky to efficiently select an optimal subset
   of available build points from which to construct the Gaussian
   Process.  Here "optimal" means that, given the current set of
   assumed correlation parameters, this subset maximizes the
   amount of unique information content in R, note that this is
   equivalent to a "best spaced" (for the chosen correlation
   function and its parameters) set of points and the output at
   those points does is not considered.  Thus if you have 2 points
   that are very close together but on opposite sides of a
   discontinutity it is highly likely that at least one of them
   will get discarded */
void KrigingModel::equationSelectingCholR(){
  if(!((buildDerOrder==0)||(buildDerOrder==1))) {
    std::cerr << "buildDerOrder=" << buildDerOrder
	      << " in void KrigingModel::equationSelectingCholR().  "
	      << "For Kriging buildDerOrder must equal 0.  "
	      << "For Gradient Enhanced Kriging (GEK) buildDerOrder "
	      << "must equal 1.  Higher order derivative enhanced "
	      << "Kriging (e.g. Hessian Enhanced Kriging) has not "
	      << "been implemented." << std::endl;
    assert(false);
  }

  //polyOrder=polyOrderRequested;
  nTrend=numTrend(polyOrderRequested,0);
  Rinv_Gtran.newSize(numEqnAvail,nTrend);


  //printf("Entered equationSelectingCholR()\n");
  double min_allowed_rcond=1.0/maxCondNum;
  //printf("min_allowed_rcond=%g\n",min_allowed_rcond);
  //exit(0);
  //double min_allowed_pivot_est_rcond=256.0/maxCondNum;

  int ld_RChol=RChol.getNRowsAct();
  //printf("ld_RChol=%d\n",ld_RChol);
  int chol_info;
  RChol.newSize(numPoints,numPoints);
  scaleRChol.newSize(numEqnAvail,3); //maximum space needed
  rcondDblWork.newSize(3*ld_RChol,1);
  rcondIntWork.newSize(ld_RChol,1);
  ld_RChol=RChol.getNRowsAct();

  iPtsKeep.newSize(numPoints,1);
  //assign the default order to points
  for(int ipt=0; ipt<numPoints; ++ipt)
    iPtsKeep(ipt,0)=ipt;

  if(buildDerOrder==0) {
    //We're using regular Kriging not GEK and for large matrices
    //the pivoted Cholesky algorithm is nowhere close to as fast
    //as the highly optimized level 3 LAPACK Cholesky so to make
    //this run faster on average we're going to attempt to use
    //the LAPACK cholesky take a look at the rcondR and then only
    //do Pivoted Cholesky if we actually need to
    RChol.copy(R);

    //no scaling is necessary since have all ones on the diagonal
    //of the Kriging R but the generic equilibrated Cholesky
    //factoriztion function doesn't know that in advance
    Chol_fact_workspace(RChol,scaleRChol,rcondDblWork,rcondIntWork,
			chol_info,rcondR);
    if(min_allowed_rcond<rcondR) {
      numRowsR=numWholePointsKeep=numPointsKeep=numPoints;

      Y.copy(Yall);

      nTrend=numTrend(polyOrderRequested,0);
      Gtran.newSize(numPoints,nTrend);
      for(int itrend=0; itrend<nTrend; ++itrend)
	for(int ipt=0; ipt<numPoints; ++ipt)
	  Gtran(ipt,itrend)=Gall(itrend,ipt);

      return;
    }
  }

  // *******************************************************************
  // in this section I need to get an optimally reordered Cholesky
  // factorization of the Kriging or GEK R matrix and I need to compute
  // the one norm for all sizes of that optimally reordered R matrix
  // *******************************************************************
  // We got here in one of two ways
  //
  // 1) We're doing Kriging and we actually need to do Pivoted Cholesky
  //
  // 2) We're doing Gradient Enhanced Kriging, and reordering "whole
  //    points" (function value immediately followed by its derivatives)
  //    according to the order of the Pivoted Cholesky on the Kriging R
  //    works a lot better and is a lot faster than doing Pivoted
  //    Cholesky on the GEK R.  In fact because the GEK R is so much
  //    larger than the Kriging R, the cost of doing pivoted Cholesky on
  //    the Kriging R will be insignificant compared to the cost of
  //    doing LAPACK Cholesky on the GEK R so I'm just going to go ahead
  //    and reorder the GEK R whether I need to or not (which I don't
  //    know at this point anyway) rather than risk having to do the
  //    LAPACK Cholesky twice

  //if the user specifies an anchor point it must be the first point to
  //prevent it from being pivoted away
  if(ifHaveAnchorPoint&&(iAnchorPoint!=0)) {
    iPtsKeep(iAnchorPoint,0)=0;
    iPtsKeep(0,0)=iAnchorPoint;
  }
  else iAnchorPoint=0;

  for(int jpt=0; jpt<numPoints; ++jpt) {
    int jsrc=iPtsKeep(jpt,0);
    for(int ipt=0; ipt<numPoints; ++ipt)
      RChol(ipt,jpt)=R(iPtsKeep(ipt,0),jsrc);
  }

  int info=0;
  char uplo='B'; //'B' means we have both halves of R in RChol so the
  //fortran doesn't have to copy one half to the other, having both
  //halves makes the memory access faster (can always go down columns)
  numPointsKeep=numPoints;
  PIVOTCHOL_F77(&uplo, &numPoints, RChol.ptr(0,0), &ld_RChol,
    		iPtsKeep.ptr(0,0), &numPointsKeep, &min_allowed_rcond,
		&info);

  //for(int ipt=0; ipt<numPoints; ++ipt)
  //printf("F77 iPtsKeep(%d)=%d\n",ipt,iPtsKeep(ipt,0));
  //printf("\n");

  //printf("*********************************\n");

  if(ifHaveAnchorPoint&&(iAnchorPoint!=0)) {
    iPtsKeep(0,0)=iAnchorPoint;
    for(int ipt=1; ipt<numPoints; ++ipt) {
      iPtsKeep(ipt,0)-=1; //Fortran indices start at 1 not zero so
      //we have to convert to C++ indices which start from 0
      if(iPtsKeep(ipt,0)==iAnchorPoint)
	iPtsKeep(ipt,0)=0;
    }
  }
  else {
    for(int ipt=0; ipt<numPoints; ++ipt) {
      iPtsKeep(ipt,0)-=1; //Fortran indices start at 1 not zero so
      //we have to convert to C++ indices which start from 0
      //printf("iPtsKeep(%2d,0)=%d\n",ipt,iPtsKeep(ipt,0));
    }
    //printf("\n");
  }

  //if I feed LAPACK a one norm of R and a Cholesky factorization of
  //R it will give me back an rcond for O(N^2) ops which is practically
  //free compared to the O(N^3) ops that the Cholesky factorization
  //costs, the wonderful thing is if I just drop equations off the end
  //of a pivoted Cholesky factorization I can get the rcondR for any
  //number of rows/columns of the pivoted R matrix, but to make this
  //efficient I need to get the one norms of R cheaply (easily doable,
  //that's what happens in the next if Kriging els if GEK statement)
  //and I'll use bisection to find the last equation I can retain
  int iprev_lapack_rcondR;
  int icurr_lapack_rcondR=numEqnAvail-1;
  int num_eqn_keep=numEqnAvail;
  oneNormR.newSize(numEqnAvail,1);
  sumAbsColR.newSize(numEqnAvail,1);
  if(buildDerOrder==0) {
    //Kriging we need to compute the one-norms for the reordered Kriging
    //R matrix
    //the one norm is the largest of the sums of the absolute value of
    //any of the columns, of course how many rows there are affects what
    //the sums of absolute value of the columns are and how many columns
    //there are affects which is the largest but we can build this up in
    //such a way as to reuse the information from smaller numbers of
    //rows/columns for larger numbers of rows/columns

    int jsrc=iPtsKeep(0,0);
    for(int ipt=0; ipt<numPoints; ++ipt)
      sumAbsColR(ipt,0)=std::fabs(R(iPtsKeep(ipt,0),jsrc));
    oneNormR(0,0)=sumAbsColR(0,0); //this is the one norm for the 1 by
    //1 reordered R matrix

    double tempdouble;
    for(int jpt=1; jpt<numPoints; ++jpt) {
      jsrc=iPtsKeep(jpt,0);
      for(int ipt=0; ipt<numPoints; ++ipt)
	sumAbsColR(ipt,0)+=std::fabs(R(iPtsKeep(ipt,0),jsrc));
      tempdouble=sumAbsColR(0,0);
      for(int ipt=1; ipt<=jpt; ++ipt)
	if(tempdouble<sumAbsColR(ipt,0))
	  tempdouble=sumAbsColR(ipt,0);
      oneNormR(jpt,0)=tempdouble; //this is the one norm for the
      //jpt by jpt reordered R matrix
    }
    uplo='L'; //get it into the same state as GEK
    iprev_lapack_rcondR=0; //a 1 by 1 matrix has a condition number of 1

  } else if(buildDerOrder==1){
    //Gradient Enhanced Kriging
    //it works better (and is a lot faster) if we reorder whole points
    //according to the Pivoted Cholesky ON THE KRIGING R order in iPtsKeep
    //so we'll calculate the one norm for all sizes of the reordered
    //GEK R and then Cholesky factorize the GEK R
    reorderCopyRtoRChol();
    /*
    printf("R=\n");
    for(int i=0; i<numEqnAvail; ++i) {
      for(int j=0; j<numEqnAvail; ++j)
	printf("%14.8g ",RChol(i,j));
      printf("\n");
    }
    printf("\n");
    */
    scaleRChol.newSize(numEqnAvail,2);
    for(int i=0; i<numEqnAvail; ++i) {
      scaleRChol(i,1)=std::sqrt(RChol(i,i));
      scaleRChol(i,0)=1.0/scaleRChol(i,1);
    }

    //equilibrate RChol
    for(int j=0; j<numEqnAvail; ++j) {
      for(int i=0; i<numEqnAvail; ++i)
	RChol(i,j)*=scaleRChol(i,0)*scaleRChol(j,0);
      RChol(j,j)=1.0; //there is zero correlation between an individual
      //point's function value and its derivatives so we know how to fix
      //round of error so just do it
    }
    /*
    printf("RE=\n");
    for(int i=0; i<numEqnAvail; ++i) {
      for(int j=0; j<numEqnAvail; ++j)
	printf("%14.8g ",RChol(i,j));
      printf("\n");
    }
    printf("\n");
    */
    //the one norm number is the largest of the sums of the absolute
    //value of any of the columns of the matrix, of course how many rows
    //there are affects what the sums of absolute value of the columns
    //are and how many columns there are affects which is the largest
    //but we can build this up in such a way as to reuse the information
    //from smaller numbers of rows/columns for larger numbers of rows/
    //columns

    //right now RChol holds the reordered R matrix
    for(int i=0; i<numEqnAvail; ++i)
      sumAbsColR(i,0)=std::fabs(RChol(i,0));
    oneNormR(0,0)=sumAbsColR(0,0); //this is the one norm for the 1 by
    //1 reordered R matrix

    double tempdouble;
    for(int j=1; j<numEqnAvail; ++j) {
      for(int i=0; i<numEqnAvail; ++i)
	sumAbsColR(i,0)+=std::fabs(RChol(i,j));
      tempdouble=sumAbsColR(0,0);
      for(int i=1; i<=j; ++i)
	if(tempdouble<sumAbsColR(i,0))
	  tempdouble=sumAbsColR(i,0);
      oneNormR(j,0)=tempdouble;  //this is the one norm for the
      //j by j reordered R matrix
    }

    //do the (highly optimized) LAPACK Cholesky Decomposition of all
    //the equations (but sorted into the point order determined by
    //the pivoting cholesky above)
    uplo='L';
    DPOTRF_F77(&uplo,&numEqnAvail,RChol.ptr(0,0),&ld_RChol,&info);

    //Kriging already has the rcondR so to get GEK into an equivalent
    //state we will feed LAPACK the one norm of the full GEK R (after
    //the reordering and equilibration) and it will give me back GEK's
    //rcondR
    DPOCON_F77(&uplo,&numEqnAvail,RChol.ptr(0,0),&ld_RChol,
	       oneNormR.ptr(icurr_lapack_rcondR,0),
	       &rcondR,rcondDblWork.ptr(0,0),rcondIntWork.ptr(0,0),
	       &info);

    //printf("rcond(RE)=%g icurr_lapack_rcondR=%d\n",rcondR,icurr_lapack_rcondR);

    //the first derivatives of the correlation at a point are uncorrelated
    //with the correlation function at the same point, i.e. the (1+numVarsr)
    //by (1+numVarsr) correlation matrix has a condition number of 1
    iprev_lapack_rcondR=numVarsr; //no 1+ because C++ indexes start at zero
  }

  // *****************************************************************
  // in this section we will efficiently determine the maximum number
  // of equations that we can retain by doing a bisection search using
  // O(log2(N)) calls of LAPACK's rcond estimate (each of which cost
  // only O(n^2) ops where n is the number of eqns in the current
  // subset
  // *****************************************************************

  lapackRcondR.newSize(numEqnAvail,1);
  lapackRcondR(iprev_lapack_rcondR,0)=1.0; //since the condition number
  //is one at iprev_lapack_rcondR we know we can keep at least this many
  //equations

  lapackRcondR(icurr_lapack_rcondR,0)=rcondR; //the maximum number
  //of equations we can keep is icurr_lapack_rcondR=numEqnAvail-1
  //and we know the rcondR for that many equations

  //note num_eqn_keep is now numEqnAvail
  int inext_lapack_rcondR=icurr_lapack_rcondR; //the last available eqn
  if((rcondR<=min_allowed_rcond)&&
     (inext_lapack_rcondR-iprev_lapack_rcondR==1)) {
    //at this point the previous lapack rcondR==1.0
    rcondR=1.0;
    inext_lapack_rcondR=iprev_lapack_rcondR;
    //printf("if1\n");
  }

  //do the bisection search if necessary, at most ceil(log2()) more
  //calls to the LAPACK rcond function
  int rcond_iter=0;
  int max_rcond_iter=
    std::ceil(std::log(static_cast<double>
		       (inext_lapack_rcondR-iprev_lapack_rcondR))
	      /std::log(2.0));
  while((lapackRcondR(inext_lapack_rcondR,0)<=min_allowed_rcond)&&
        (inext_lapack_rcondR>iprev_lapack_rcondR)) {
    //printf("inWhile\n");
    ++rcond_iter;
    icurr_lapack_rcondR=(iprev_lapack_rcondR+inext_lapack_rcondR)/2;
    num_eqn_keep=icurr_lapack_rcondR+1;

    //the LAPACK rcond function
    DPOCON_F77(&uplo,&num_eqn_keep,RChol.ptr(0,0),&ld_RChol,
	       oneNormR.ptr(icurr_lapack_rcondR,0),
	       &rcondR,rcondDblWork.ptr(0,0),rcondIntWork.ptr(0,0),
	       &info);
    lapackRcondR(icurr_lapack_rcondR,0)=rcondR;
    //printf("rcond_iter=%d icurr_lapack_rcondR=%d rcondR=%g\n",
    //rcond_iter,icurr_lapack_rcondR,rcondR);

    if(rcondR<min_allowed_rcond)
      inext_lapack_rcondR=icurr_lapack_rcondR;
    else if(min_allowed_rcond<rcondR)
      iprev_lapack_rcondR=icurr_lapack_rcondR;
    else if(min_allowed_rcond==rcondR) {
      //num_eqn_keep=icurr_lapack_rcondR+1;
      break;
    }
    if((inext_lapack_rcondR-iprev_lapack_rcondR==1)||
       (max_rcond_iter<rcond_iter)) {
      num_eqn_keep=iprev_lapack_rcondR+1;
      rcondR=lapackRcondR(iprev_lapack_rcondR,0);
      break;
    }
  }
  //printf(" pivoted_rcondR=%g numRowsR=%d\n",rcondR,num_eqn_keep);

  numRowsR=num_eqn_keep; //this is the maximum number of equations that
  //we can keep

  // ***************************************************************
  // in this section we downsize the arrays being retained and keep
  // only the optimal subset, in or working copies
  // ***************************************************************

  RChol.resize(num_eqn_keep,num_eqn_keep); //resize() instead of newSize()
  //because we want to keep the current contents in the same 2D
  //order
  /*
  if(num_eqn_keep>=10) {
    printf("RChol(1:10,1:10)=\n");
    for(int i=0; i<10; ++i) {
      for(int j=0; j<10; ++j)
	printf("%12.6g ",RChol(i,j));
      printf("\n");
    }
    printf("\n\n");
  }
  */

  //polyOrder=polyOrderRequested; //redundant but for clarity

  //the following while loop was commented out when adaptive selection of
  //the trend basis functions via pivote cholesky factorization of G*R^-1*G^T
  //was implemented
  //while((numRowsR<=numTrend(polyOrder,0))&&(polyOrder>0))
  //--polyOrder;

  //nTrend=numTrend(polyOrder,0); //commented out because we now select a subset of Poly based using a Pivoted Cholesky factorization of G*R^-1*G^T which happens when trendSelectingPivotedCholesy is called by materObjectiveAndConstraints() (we no longer select SOLELY on polynomial order and number of points

  nTrend=numTrend(polyOrderRequested,0);

  //printf("num_eqn_keep=%d numRowsR=%d polyOrder=%d nTrend=%d rcondR=%g lapackRcondR(num_eqn_keep-1,0)=%g\n",num_eqn_keep,numRowsR,polyOrder,nTrend,rcondR,lapackRcondR(num_eqn_keep-1,0));

  //we need to downsize Gtran now but we only need to downsize Poly at
  //the end of create()
  Gtran.newSize(num_eqn_keep,nTrend); //newSize() because we don't care
 //about the current contents of Gtran

  Y.newSize(num_eqn_keep,1); //newSize() because we don't care about
  //the current contents of Y

  if(buildDerOrder==0) {
    // keep only the useful parts for Kriging
    numWholePointsKeep=numPointsKeep=num_eqn_keep;
    numExtraDerKeep=0;

    /*
    if(numPointsKeep>10) {

      MtxDbl RCholDEBUG(numPointsKeep,numPointsKeep);
      for(int jpt=0; jpt<numPointsKeep; ++jpt) {
	int jsrc=iPtsKeep(jpt,0);
	for(int ipt=0; ipt<numPointsKeep; ++ipt)
	  RCholDEBUG(ipt,jpt)=R(iPtsKeep(ipt,0),jsrc);
      }
      double rcondRDEBUG;
      int chol_info_debug;

      printf("Rreorder(1:10,1:10)=\n");
      for(int ipt=0; ipt<10; ++ipt) {
	for(int jpt=0; jpt<10; ++jpt)
	  printf("%12.6g ",RCholDEBUG(ipt,jpt));
	printf("\n");
      }
      printf("\n\n");

      Chol_fact_workspace(RCholDEBUG,scaleRChol,rcondDblWork,rcondIntWork,
			  chol_info_debug,rcondRDEBUG);

      printf("RChol(1:10,1:10)=\n");
      for(int ipt=0; ipt<10; ++ipt) {
	for(int jpt=0; jpt<10; ++jpt)
	  printf("%12.6g ",RChol(ipt,jpt));
	printf("\n");
      }
      printf("\n\n");

      printf("RCholDEBUG(1:10,1:10)=\n");
      for(int ipt=0; ipt<10; ++ipt) {
	for(int jpt=0; jpt<10; ++jpt)
	  printf("%12.6g ",RCholDEBUG(ipt,jpt));
	printf("\n");
      }
      printf("\n\n");

      printf("[RChol-RCholDEBUG](1:10,1:10)=\n");
      for(int ipt=0; ipt<10; ++ipt) {
	for(int jpt=0; jpt<10; ++jpt)
	  printf("%12.6g ",RChol(ipt,jpt)-RCholDEBUG(ipt,jpt));
	printf("\n");
      }
      printf("\n\n");

      printf("rcondR=%g rcondRDEBUG=%g numPointsKeep=%d\nErrorRChol=\n",
	     rcondR,rcondRDEBUG,numPointsKeep);

    }
    */

    //keep the useful part of Y
    for(int ipt=0; ipt<numPointsKeep; ++ipt)
      Y(ipt,0)=Yall(iPtsKeep(ipt,0),0);

    //keep the useful part of G (actually G^T)
    for(int itrend=0; itrend<nTrend; ++itrend)
      for(int ipt=0; ipt<numPointsKeep; ++ipt)
	Gtran(ipt,itrend)=Gall(itrend,iPtsKeep(ipt,0));

    /*
    for(int ipt=0; ipt<numPointsKeep; ++ipt) {
      printf("Gtran(%3d,:)=[%12.6g",ipt,Gtran(ipt,0));
      for(int itrend=1; itrend<nTrend; ++itrend)
	printf(", %12.6g",Gtran(ipt,itrend));
      printf("] XR(:,%3d)=[%12.6g",ipt,XR(0,iPtsKeep(ipt,0)));
      for(int k=1; k<numVarsr; ++k)
	printf(", %12.6g",XR(k,iPtsKeep(ipt,0)));
      printf("]^T Y(%3d,0)=%12.6g\n",ipt,Y(ipt,0));
    }
    printf("\n");
    */


  } else if(buildDerOrder==1) {
    // keep on the useful parts for Gradient Ehanced Kriging

    //integer division automatically rounds down
    numWholePointsKeep=num_eqn_keep/(1+numVarsr);

    //we also need to round up
    numPointsKeep=
      static_cast<int>(std::ceil(static_cast<double>(num_eqn_keep)/
				 static_cast<double>(1+numVarsr)));

    if(numPointsKeep==numWholePointsKeep) {
      //perhaps a better name would be numLastDerKeep... this is the number
      //of derivatives retained for the last point.
      numExtraDerKeep==numVarsr;
    } else
      numExtraDerKeep=num_eqn_keep-(1+numWholePointsKeep*(1+numVarsr));

    //we need to undo the equilibration of RChol, recall that scaleRChol
    //is already in the pivoted Cholesky order
    for(int j=0; j<num_eqn_keep; ++j)
      for(int i=j; i<num_eqn_keep; ++i)
	RChol(i,j)*=scaleRChol(i,1); //note that this assumes that the
    //nkm::SurfMat class uses the lower triangular part of of RChol
    //otherwise (if this function uses the lower triangular part but
    //surfmat uses the upper triangular part) you'd need
    //RChol(j,i)=RChol(i,j)*scaleRChol(i,j) but you could do a
    //RChol(j,i)=RChol(i,j)*=ScaleRChol(i,1); just to be safe

    //keep the useful part of G (actually G^T)
    for(int itrend=0; itrend<nTrend; ++itrend) {
      int i, ipt;
      for(ipt=0, i=0; ipt<numWholePointsKeep; ++ipt) {
	int isrc=iPtsKeep(ipt,0)*(1+numVarsr);
	for(int k=0; k<1+numVarsr; ++k, ++isrc, ++i)
	  Gtran(i,itrend)=Gall(itrend,isrc);
      }
      if(numPointsKeep>numWholePointsKeep) {
#ifdef __KRIG_ERR_CHECK__
	assert((ipt==numWholePointsKeep)&&
	       (i==(numWholePointsKeep*(1+numVarsr)))&&
	       (numWholePointsKeep+1==numPointsKeep));
#endif
	int isrc=iPtsKeep(ipt,0)*(1+numVarsr);
	Gtran(i,itrend)=Gall(itrend,isrc);
	++i;
	++isrc;
	for(int k=0; k<numExtraDerKeep; ++k, ++isrc, ++i)
	  Gtran(i,itrend)=Gall(itrend,isrc);
      }
    }

    //keep the useful part of Y, also we need to undo the
    //equilibration of RChol
    { //{ to impose scope
      int i, ipt;
      for(ipt=0, i=0; ipt<numWholePointsKeep; ++ipt) {
	int isrc=iPtsKeep(ipt,0)*(1+numVarsr);
	for(int k=0; k<1+numVarsr; ++k, ++isrc, ++i)
	  Y(i,0)=Yall(isrc,0);
      }
      if(numPointsKeep>numWholePointsKeep) {
#ifdef __KRIG_ERR_CHECK__
	assert((ipt==numWholePointsKeep)&&
	       (i==(numWholePointsKeep*(1+numVarsr)))&&
	       (numWholePointsKeep+1==numPointsKeep));
#endif
	int isrc=iPtsKeep(ipt,0)*(1+numVarsr);
	Y(i,0)=Yall(isrc,0);
	++i;
	++isrc;
	for(int k=0; k<numExtraDerKeep; ++k, ++isrc, ++i)
	  Y(i,0)=Yall(isrc,0);
      }
    } //{} to impose scope
  } //else if(buildDerOrder==1) {

  iPtsKeep.resize(numPointsKeep,1);

  return;
}

/** G_Rinv_Gtran must be filled with G*R^-1*G^T prior to calling
    void KrigingModel::trendSelectingPivotedCholesky() */
void KrigingModel::trendSelectingPivotedCholesky(){

  //nTrend=numTrend(polyOrderRequested,0);
  iTrendKeep.newSize(nTrend,1);
  double min_allowed_rcond=1.0/maxCondNum;
  int num_trend_want=numTrend(polyOrderRequested,0);

  int max_trend_to_keep=numRowsR-1-2*numVarsr;
  if(numRowsR/2<max_trend_to_keep)
    max_trend_to_keep=numRowsR/2;
  if(num_trend_want<max_trend_to_keep)
    max_trend_to_keep=num_trend_want;
  if(max_trend_to_keep<1)
    max_trend_to_keep=1;

  //std::cout << "numRowrR=" << numRowsR << " nTrend=" << nTrend << " max_trend_to_keep=" << max_trend_to_keep << std::endl;

  if(nTrend<=max_trend_to_keep) {
    //do a LAPACK cholesky first
    G_Rinv_Gtran_Chol.copy(G_Rinv_Gtran);
    int chol_info;
    Chol_fact_workspace(G_Rinv_Gtran_Chol,G_Rinv_Gtran_Chol_Scale,
			G_Rinv_Gtran_Chol_DblWork,G_Rinv_Gtran_Chol_IntWork,
			chol_info,rcond_G_Rinv_Gtran);
    if(min_allowed_rcond<rcond_G_Rinv_Gtran) {
      //the LAPACK Cholesky was not ill-conditioned with the desired
      //full set of trend basis functions so we'll use them all
      for(int itrend=0; itrend<nTrend; ++itrend)
	iTrendKeep(itrend,0)=itrend;
      return;
    }
  }

  // *****************************************************************
  // if we got here we need to do a Pivoted Cholesky factorization of
  // G*R^-1*G^T to select an optimal subset of trend basis functions
  // *****************************************************************

  // we have a slight problem, basically the one trend function that
  // is guaranteed to be selected is the one that has the largest
  // diagonal element in G*R^-1*G^T, or if they are all the same it
  // will be the first one.  we want to guarantee that the constant
  // basis function (polynomial of order 0) is retained so we need to
  // equilibrate G_Rinv_Gtran so that it has all ones on the diagonal
  // and therefore that the constant trend basis function will be
  // retained

  // I realize this has R in the name but I don't want to allocate another
  // variable
  scaleRChol.newSize(nTrend,3);
  for(int itrend=0; itrend<nTrend; ++itrend) {
    scaleRChol(itrend,1)=std::sqrt(G_Rinv_Gtran(itrend,itrend));
    scaleRChol(itrend,0)=1.0/scaleRChol(itrend,1);
  }


  for(int jtrend=0; jtrend<nTrend; ++jtrend) {
    double tempdouble=scaleRChol(jtrend,0);
    for(int itrend=0; itrend<nTrend; ++itrend)
      G_Rinv_Gtran(itrend,jtrend)*=scaleRChol(itrend,0)*tempdouble;
    G_Rinv_Gtran(jtrend,jtrend)=1.0; //numerical round off error could
    //kill our guarantee of keeping the constant trend basis function
    //so we'll just assign the right correct answer which is 1.0
  }

  G_Rinv_Gtran_Chol.copy(G_Rinv_Gtran);

  int ld_G_Rinv_Gtran_Chol=G_Rinv_Gtran_Chol.getNRowsAct();
  int info=0;
  char uplo='B'; //'B' means we have both halves of G*R^-1*G^T in
  //G_Rinv_Gtran_Chol so the FORTRAN doesn't have to copy one half to the
  //other, having both halves makes the memory access faster (can always
  //go down columns)
  int num_trend_keep=-max_trend_to_keep; //RANK<0 on input tells PIVOTCHOL_F77
  //that we only want to keep the first abs(RANK) entries so it can stop early
  PIVOTCHOL_F77(&uplo, &nTrend, G_Rinv_Gtran_Chol.ptr(0,0),
		&ld_G_Rinv_Gtran_Chol, iTrendKeep.ptr(0,0), &num_trend_keep,
		&min_allowed_rcond, &info);

  nTrend=num_trend_keep; //this is the maximum number of trend functions
  //we could keep, we might not be able to keep this many

  //FORTRAN indices start at 1 not zero so we have to convert to C++ indices
  //which start at zero
  for(int itrend=0; itrend<nTrend; ++itrend)
    iTrendKeep(itrend,0)-=1;

  // ************************************************************
  // in this section I need to calculate the one norm for subsets
  // of the first jtrend (reordered) rows/columns of the
  // equilibrated G*R^-1*G^T , for jtrend=1 to nTrend
  // ************************************************************

  //I realize that this says R but I don't want to allocate another variable
  oneNormR.newSize(nTrend,1);
  //I realize that this says R but I don't want to allocate another variable
  sumAbsColR.newSize(nTrend,1);
  int jsrc=iTrendKeep(0,0);
  for(int itrend=0; itrend<nTrend; ++itrend)
    sumAbsColR(itrend,0)=std::fabs(G_Rinv_Gtran(iTrendKeep(itrend,0),jsrc));
  oneNormR(0,0)=sumAbsColR(0,0); //this is the one norm for the 1 by
  //1 reordered G_Rinv_Gtran matrix

  double tempdouble;
  for(int jtrend=1; jtrend<nTrend; ++jtrend) {
    jsrc=iTrendKeep(jtrend,0);
    for(int itrend=0; itrend<nTrend; ++itrend)
      sumAbsColR(itrend,0)+=std::fabs(G_Rinv_Gtran(iTrendKeep(itrend,0),jsrc));
    tempdouble=sumAbsColR(0,0);
    for(int itrend=1; itrend<=jtrend; ++itrend)
      if(tempdouble<sumAbsColR(itrend,0))
	tempdouble=sumAbsColR(itrend,0);
    oneNormR(jtrend,0)=tempdouble; //this is the one norm for the
    //jtrend by jtrend reordered G_Rinv_Gtran matrix
  }

  ld_G_Rinv_Gtran_Chol=G_Rinv_Gtran_Chol.getNRowsAct(); //this probably hasn't changed
  //but better safe than sorry
  rcondDblWork.newSize(3*ld_G_Rinv_Gtran_Chol,1);
  rcondIntWork.newSize(ld_G_Rinv_Gtran_Chol,1);
  int icurr_lapack_rcond=nTrend-1;
  uplo='L';
  DPOCON_F77(&uplo,&nTrend,G_Rinv_Gtran_Chol.ptr(0,0),&ld_G_Rinv_Gtran_Chol,
	     oneNormR.ptr(icurr_lapack_rcond,0),&rcond_G_Rinv_Gtran,
	     rcondDblWork.ptr(0,0),rcondIntWork.ptr(0,0),&info);

  //need to do a bisection search for the last trend function that we can keep
  lapackRcondR(icurr_lapack_rcond,0)=rcond_G_Rinv_Gtran; //the maximum number
  //of trend basis functions we can keep is icurr_lapack_rcond=nTrend-1
  //and we know the rcond for that many equations
  int iprev_lapack_rcond=0;
  lapackRcondR(iprev_lapack_rcond,0)=1.0; //the constant trend funcation by
  //itself is guaranteed to have condition # = 1.0

  //note num_trend_keep is now nTrend
  int inext_lapack_rcond=icurr_lapack_rcond; //the last available basis
  //function
  if((rcond_G_Rinv_Gtran<=min_allowed_rcond)&&
     (inext_lapack_rcond-iprev_lapack_rcond==1)) {
    //at this point the previous lapack rcondR==1.0
    rcond_G_Rinv_Gtran=1.0;
    inext_lapack_rcond=iprev_lapack_rcond;
    //printf("if1\n");
  }

  //do the bisection search if necessary, at most ceil(log2()) more
  //calls to the LAPACK rcond function
  int rcond_iter=0;
  int max_rcond_iter=
    std::ceil(std::log(static_cast<double>
		       (inext_lapack_rcond-iprev_lapack_rcond))
	      /std::log(2.0));
  while((lapackRcondR(inext_lapack_rcond,0)<=min_allowed_rcond)&&
        (inext_lapack_rcond>iprev_lapack_rcond)) {
    //printf("inWhile\n");
    ++rcond_iter;
    icurr_lapack_rcond=(iprev_lapack_rcond+inext_lapack_rcond)/2;
    num_trend_keep=icurr_lapack_rcond+1;

    //the LAPACK rcond function
    DPOCON_F77(&uplo,&num_trend_keep,G_Rinv_Gtran_Chol.ptr(0,0),
	       &ld_G_Rinv_Gtran_Chol,oneNormR.ptr(icurr_lapack_rcond,0),
	       &rcond_G_Rinv_Gtran,rcondDblWork.ptr(0,0),
	       rcondIntWork.ptr(0,0),&info);
    lapackRcondR(icurr_lapack_rcond,0)=rcond_G_Rinv_Gtran;
    //printf("rcond_iter=%d icurr_lapack_rcondR=%d rcondR=%g\n",
    //rcond_iter,icurr_lapack_rcondR,rcondR);

    if(rcond_G_Rinv_Gtran<min_allowed_rcond)
      inext_lapack_rcond=icurr_lapack_rcond;
    else if(min_allowed_rcond<rcond_G_Rinv_Gtran)
      iprev_lapack_rcond=icurr_lapack_rcond;
    else if(min_allowed_rcond==rcond_G_Rinv_Gtran) {
      //num_trend_keep=icurr_lapack_rcond+1;
      break;
    }
    if((inext_lapack_rcond-iprev_lapack_rcond==1)||
       (max_rcond_iter<rcond_iter)) {
      num_trend_keep=iprev_lapack_rcond+1;
      rcond_G_Rinv_Gtran=lapackRcondR(iprev_lapack_rcond,0);
      break;
    }
  }
  //printf(" pivoted_rcondR=%g numRowsR=%d\n",rcondR,num_eqn_keep);

  nTrend=num_trend_keep; //this is the maximum number of trend basis functions
  //that we can keep

  iTrendKeep.resize(nTrend,1);
  //we don't want the basis functions in their optimal order we want them
  //in their logical order, minus the ones we couldn't keep
  iTrendKeep.sortRows();

  //were going to copy the portion of the equilibrated G*R^-1*G^T matrix
  //that were keeping, in its logical order, into G_Rinv_Gtran_Chol and then
  //do a LAPACK cholesky on it, and then undo the equilibration
  G_Rinv_Gtran_Chol.newSize(nTrend,nTrend);
  for(int jtrend=0; jtrend<nTrend; ++jtrend) {
    int jsrc=iTrendKeep(jtrend,0);
    scaleRChol(jtrend,2)=scaleRChol(jsrc,1);
    for(int itrend=0; itrend<nTrend; ++itrend)
      G_Rinv_Gtran_Chol(itrend,jtrend)=G_Rinv_Gtran(iTrendKeep(itrend,0),jsrc);
  }

  //do the LAPACK cholesky factorization
  int chol_info;
  Chol_fact_workspace(G_Rinv_Gtran_Chol,G_Rinv_Gtran_Chol_Scale,
		      G_Rinv_Gtran_Chol_DblWork,G_Rinv_Gtran_Chol_IntWork,
		      chol_info,rcond_G_Rinv_Gtran);

  //we still need to undo the equilbration because we copied the
  //equilibrated G*R^-1*G^T matrix into G_Rinv_Gtran_Chol before doing
  //the cholesky factorization
  for(int jtrend=0; jtrend<nTrend; ++jtrend)
    for(int itrend=jtrend; itrend<nTrend; ++itrend)
      G_Rinv_Gtran_Chol(itrend,jtrend)*=scaleRChol(itrend,2);



  for(int itrend=1; itrend<nTrend; ++itrend) {
    int isrc=iTrendKeep(itrend,0);
    if(itrend<isrc)
      for(int i=0; i<numRowsR; ++i) {
	Gtran(i,itrend)=Gtran(i,isrc);
	Rinv_Gtran(i,itrend)=Rinv_Gtran(i,isrc);
      }
  }
  Gtran.resize(numRowsR,nTrend);
  Rinv_Gtran.resize(numRowsR,nTrend);

  //and were done with the pivoted cholesky, but at the end of the optimization
  //we will still need to discard the subset of Poly that was not selected by
  //trendSelectingPivotedCholesky()
  return;
}





/** this function calculates the objective function (negative "per equation"
    log likelihood) and/or the constraint (reciprocal condition number)
    functions using a precompute and store (store across sequential calls
    to this function) strategy to reduce the computational cost, make sure
    only to COPY OUT results from member variables so the state is not
    changed
*/
void KrigingModel::masterObjectiveAndConstraints(const MtxDbl& theta, int obj_der_mode, int con_der_mode)
{
  // if(obj_der_mode=1) (1=2^0=> 0th derivative) calculate objective function
  // if(con_der_mode=1) (1=2^0=> 0th derivative) calculate the constraint
  //functions
  // ERROR if(con_der_mode>=2) (2=2^1 = 1st derivative) this function does not
  //                           support analytical derivatives of the objective
  //                           function
  // ERROR if(con_der_mode>=2) (2=2^1 = 1st derivative) this function does not
  //                           support analytical derivatives of the constraint
  //                           function

  //printf("maxConDerMode=%d con_der_mode=%d maxObjDerMode=%d obj_der_mode=%d\n",
  //maxConDerMode,con_der_mode,maxObjDerMode,obj_der_mode);

  //might want to replace this with a thrown exception
  assert((maxObjDerMode<=1)&&(maxConDerMode<=1)&&
	 (0<=obj_der_mode)&&(obj_der_mode<=maxObjDerMode)&&
	 (0<=con_der_mode)&&(con_der_mode<=maxConDerMode)&&
	 ((1<=obj_der_mode)||(1<=con_der_mode)));

  //if theta was the same as the last time we called this function than we can reuse some of the things we calculated last time

  if(prevTheta.getNElems()!=numTheta) {
    //different number of elements means we can't reuse
    prevTheta.newSize(numTheta,1);
    prevObjDerMode=prevConDerMode=0;
  }
  else
    for(int k=0; k<numTheta; ++k)
      if(prevTheta(k,0)!=theta(k,0)) {
	//some parameter changed so we can't reuse
	prevObjDerMode=prevConDerMode=0;
	break;
      }

  if((obj_der_mode<=prevObjDerMode)&&
     (con_der_mode<=prevConDerMode)) {
    //we've already calculated everything you just asked for so reuse it
    return;
  }

  //record the current theta as the previous theta so next time we can tell
  //if we should reuse the stuff we calculate this time
  if((prevObjDerMode==0)&&(prevConDerMode==0))
    for(int k=0; k<numTheta; ++k)
      prevTheta(k,0)=theta(k,0);

  if(prevObjDerMode==0) {
    //fill R with the build data "correlation matrix" (R is a member variable)
    //for Kriging R is actually a correlation matrix (it is real, symmetric,
    //positive definite and has all ones on the diagonal) for GEK it is
    //real symmetric, and positive definite but does not have all ones on the
    //diagonal, but the GEK R can be equilibrated/scaled to an honest to
    //goodness correlation matrix.
    correlation_matrix(theta);

    //we need to perform a LU decomposition of R and calculate the
    //determinant of R, I have replaced LU with Cholesky because it's
    //better/faster, see
    //http://en.wikipedia.org/wiki/Determinant#Determinant_from_LU_decomposition
    //for how to efficiently compute the determinant from an LU factorization

    int chol_info=0;
    if(ifPrescribedNug==true) {
      //the user prescribed a nugget for us to use, e.g. for when there is
      //measurement error of known magnitude
      apply_nugget_build(); //modify R by a nugget in place
      reorderCopyRtoRChol();

      Chol_fact_workspace(RChol,scaleRChol,rcondDblWork,rcondIntWork,
			  chol_info,rcondR);
      //Pivoted Cholesky o G*R^-1*G^T does not require pivoted Cholesky of R
      //so the size of Gtran could have changed
      nTrend=numTrend(polyOrderRequested,0);
      if(Gtran.getNCols() < nTrend) {
	Gtran.newSize(numEqnAvail,nTrend);
	for(int itrend=0; itrend<nTrend; ++itrend)
	  for(int i=0; i<numEqnAvail; ++i)
	    Gtran(i,itrend)=Gall(itrend,i);
      }

    } else if(ifChooseNug==true) {
      //the user wants us to select a small nugget to fix ill-conditioning of R
      nuggetSelectingCholR();
      //Pivoted Cholesky o G*R^-1*G^T does not require pivoted Cholesky of R
      //so the size of Gtran could have changed
      nTrend=numTrend(polyOrderRequested,0);
      if(Gtran.getNCols() < nTrend) {
	Gtran.newSize(numEqnAvail,nTrend);
	for(int itrend=0; itrend<nTrend; ++itrend)
	  for(int i=0; i<numEqnAvail; ++i)
	    Gtran(i,itrend)=Gall(itrend,i);
      }
    }else {
      //the user wants us to fix ill-conditioning of R by using Pivoted Cholesky
      //to select an optimal subset of points from which to build the Kriging
      //(or Gradient Enhanced Kriging) model
      equationSelectingCholR();
    }
    double min_allowed_rcond=1.0/maxCondNum;
    //nTrend=numTrend(polyOrder,0);
    nTrend=numTrend(polyOrderRequested,0);

    if((rcondR<=min_allowed_rcond)) { //||(numRowsR<=nTrend)) {
      printf("singular correlation matrix rcondR=%g numRowsR=%d numTrend=%d numEqnAvail=%d\n",
	     rcondR,numRowsR,nTrend,numEqnAvail);
      MtxDbl corr_len_temp(numVarsr,1);
      get_corr_len_from_theta(corr_len_temp, theta);
      printf("corr_len=[%g",corr_len_temp(0,0));
      for(int kk=1; kk<numVarsr; ++kk)
	printf(",%g",corr_len_temp(kk,0));
      printf("]^T\n");

      obj=HUGE_VAL; //the objective would actually be infinite, but it might
      //say nan if we let it continue and we don't want to trust the optimizer
      //to handle nan's correctly

      con.newSize(numConFunc,1);
      con(0,0)=1.0-rcondR*maxCondNum;
      //there should only be 1 constraint but just in case we'll fill the rest
      //as being violated
      for(int i=1; i<numConFunc; ++i)
	con(i,0)=1.0; //say the constraints are violated,

      //no point in wasting computation on something useless by continuing so
      //return early
      return;
    }

    double log_determinant_R = 0.0; //need to do this to avoid underflow error for large numbers of points, log(0)=-inf
    for (int i = 0; i < numRowsR; ++i)
      log_determinant_R += std::log(RChol(i,i));
    log_determinant_R *= 2.0; //only multiply by 2 for Cholesky factorization
    //of R because det(L)=det(U) and det(R)=det(L)*det(U)=det(L)^2
    //so log(det(R))=2*log(det(L))

    //if a future developer wants to switch back from cholesky to LU (and I
    //strongly recommend against that) you'll need to do a
    //determinant_R=std::fabs(determinant_R); //for LU factorization
    //because "The determinant of a positive definite matrix is always
    //positive" http://mathworld.wolfram.com/PositiveDefiniteMatrix.html and
    //det(R)=det(pivot Mtx)*det(L)*det(U); det(L)=1, det(U) is what we'd
    //calculated above for LU and det(pivot Mtx)= +/- 1, which is why you'd
    //need to do the fabs(det(U)) if you used LU decomp instead of Cholesky

    //Do the generalized (by R^-1) least squares using min # of ops
    //printf("numPoints=%d numPointsKeep=%d numRowsR=%d nTrend=%d\n",
    //   numPoints,numPointsKeep,numRowsR,nTrend);


    Rinv_Gtran.newSize(numRowsR,nTrend); //precompute and store
    solve_after_Chol_fact(Rinv_Gtran,RChol,Gtran);

    G_Rinv_Gtran.newSize(nTrend,nTrend);
    matrix_mult(G_Rinv_Gtran,Gtran,Rinv_Gtran,0.0,1.0,'T','N');

    trendSelectingPivotedCholesky();
    //Chol_fact_workspace(G_Rinv_Gtran_Chol,G_Rinv_Gtran_Chol_Scale,G_Rinv_Gtran_Chol_DblWork,G_Rinv_Gtran_Chol_IntWork,chol_info,rcond_G_Rinv_Gtran);
    if((rcond_G_Rinv_Gtran<min_allowed_rcond)||(numRowsR<=nTrend)) {
      //we could instead use pivoted cholesky to adaptively selected an optimal
      //subset of trend basis functions (i.e. it could be lower in some
      //dimensions than in others or have quadratic but not linear in certain
      //dimensions etc) then we wouldn't have to worry about this

      std::cerr << "R is not singular but G*R^-1*G^T is numerically "
		<< "singular.  This is probably\ndue to you not having "
		<< "enough UNIQUE values in one of your input dimensions\n"
		<< "to support the utilized trend function even though "
		<< "the total number of\npoints would normally be "
		<< "sufficient for the selected trend." << std::endl;
      obj=HUGE_VAL; //the objective would actually be infinite, but it might
      //say nan if we let it continue and we don't want to trust the optimizer
      //to handle nan's correctly

      con.newSize(numConFunc,1);

      //there should only be 1 constraint but just in case we'll fill them all
      //as being violated
      for(int i=0; i<numConFunc; ++i)
	con(i,0)=1.0; //say the constraints are violated,

      //no point in wasting computation on something useless by continuing so
      //return early
      return;
    }

#ifdef __KRIG_ERR_CHECK__
    assert(chol_info==0);  //for debug, do something else for production
#endif

    double log_determinant_G_Rinv_Gtran=0.0;
    for (int itrend = 0; itrend < nTrend; ++itrend)
      log_determinant_G_Rinv_Gtran +=
	std::log(G_Rinv_Gtran_Chol(itrend,itrend));
    log_determinant_G_Rinv_Gtran *= 2.0; //only for Cholesky factorization

    G_Rinv_Y.newSize(nTrend,1);
    matrix_mult(G_Rinv_Y, Rinv_Gtran, Y, 0.0, 1.0, 'T', 'N');
    betaHat.newSize(nTrend,1);

    solve_after_Chol_fact(betaHat,G_Rinv_Gtran_Chol,G_Rinv_Y); //O(nTrend^2) ops
    eps.copy(Y); //this will be eps=epsilon=Y-G(XR)^T*betaHat
    matrix_mult(eps, Gtran, betaHat, 1.0, -1.0, 'N', 'N'); //eps=Y-G(XR)^T*betaHat
    rhs.newSize(numRowsR,1);
    solve_after_Chol_fact(rhs,RChol,eps);


    //it's actually the log likelihood, which we want to maximize
    //likelihood = -0.5*(numPoints*(std::log(4.0*std::acos(0.0))+std::log(estVarianceMLE)+1)
    //		       +std::log(determinant_R)); //from Koehler and Owen

#ifdef __NKM_UNBIASED_LIKE__
    //derived following: C. E. Rasmussen & C. K. I. Williams, Gaussian Processes for Machine Learning, the MIT Press, 2006, ISBN 026218253X. c 2006 Massachusetts Institute of Technology. www.GaussianProcess.org/gpml...  we assume a "vague prior" (i.e. that we don't know anything) for betaHat, then like "Koehler and Owen" we replace the covariance matrix K with (unadjusted variance)*R (where R is the correlation matrix) and find unadjusted variance and betaHat through maximum likelihood.

    //the unbiased estimate of unadjusted variance
    estVarianceMLE = dot_product(eps,rhs)/(numRowsR-nTrend);

    //the "per equationt" unbiased log(likelihood)
    likelihood = -0.5*(std::log(estVarianceMLE)+(log_determinant_R+log_determinant_G_Rinv_Gtran)/(numRowsR-nTrend));
#else
    //derived the "Koehler and Owen" way (assumes we know the trend function, and is therefore biased

    //the estimate of unadjusted variance
    estVarianceMLE = dot_product(eps,rhs)/numRowsR; //the "Koehler and Owen" way

    //the "per equation" log(likelihood)
    likelihood = -0.5*(std::log(estVarianceMLE)+log_determinant_R/numRowsR);
#endif

    //if(likelihood>=DBL_MAX)
    //printf("[estVarianceMLE=%g determinant_R=%g]",estVarianceMLE,determinant_R);

    //the objective function being MINIMIZED is the negative of the log
    //likelihood (on a per equation basis so numbers will be comparable
    //regardless of how many equations there are)
    obj=-likelihood;
    //printf("[obj=%g]",obj);

    prevObjDerMode=1; //increase prevObjDerMode to the current value
    if((obj_der_mode==1)&&(con_der_mode<=prevConDerMode)) {
      //we have everything we need so exit early
      return;
    }
  }

  if((prevConDerMode==0)&&(1<=con_der_mode)) {
    //calculate the constraint on reciprocal condition number that ensures
    //that the correlation matrix is well conditioned.
    con.newSize(numConFunc,1);

    if(!(1<=prevObjDerMode))
      std::cerr << "We need to have already calculated rcondR (during the "
		<< "calculation of the\nobjective function) in order to "
		<< "calculate the constraint (on rcondR)\nfunction (where "
		<< "rcondR is the reciprocal of the condition number of R,\n"
		<< "and R is the ''correlation matrix'')." << std::endl;
    else if(!(numConFunc==1))
      std::cerr << "The calling function is asking us for more than one "
		<< "constraint function\nbut we only have one constraint "
		<< "function; only rcondR (the reciprocal of\nthe "
		<< "condition number of the ''correlation matrix'', R) is "
		<< "constrained." << std::endl;
    assert((1<=prevObjDerMode)&&(numConFunc==1));

    //the matrix is considered "ill-conditioned" if the following constraint
    //equation is greater than zero
    con(0,0)=1.0-rcondR*maxCondNum;

    prevConDerMode=1; //increase prevConDerMode to current value
    if((con_der_mode==1)&&(obj_der_mode<=prevObjDerMode)) {
      //we have everything we need so exit early
      return;
    }
  }

  return;
}


void KrigingModel::getRandGuess(MtxDbl& guess) const
{
  int mymod = 1048576; //2^20 instead of 10^6 to be kind to the computer
  guess.newSize(numVarsr,1);
  for(int k=0; k<numVarsr; k++) {
    guess(k,0) = (std::rand() % mymod)*(maxNatLogCorrLen-minNatLogCorrLen)/mymod+
      minNatLogCorrLen; //this returns a random nat_log_corr_len which is the space we need to search in
  }
  return;
}

// BMA TODO: These need to be moved to optimizer and then any defauls
// overridden here

void KrigingModel::set_conmin_parameters(OptimizationProblem& opt) const
{
  //set conmin specific parameters for this problem
  if((maxObjDerMode==1)&&(maxConDerMode==1)) {
    //use numerical  gradients of objective and constraints
    opt.conminData.nfdg = 0;
  } else {
    std::cerr << "This Kriging/Gradient-Enhanced-Kriging model does not "
	      << "support analytical\nderivatives of the objective "
	      << "(negative per equation log likelihood) or\nconstraint "
	      << "(reciprocal condition number) functions." << std::endl;
    assert(false);
  }

  opt.conminData.iprint = 0; //ammount of to screen output from Conmin
  opt.conminData.itmax  = maxTrials; //maximum # of Conmin iterations
  opt.conminData.fdch   = 1.0e-2; //Relative finite difference step size.
  opt.conminData.fdchm  = 1.0e-2; //Absolute finite difference step size.
  opt.conminData.ct     = -0.1; // Constraint thickness parameter, The absolute value of CT decreases in magnitude during optimization.
  opt.conminData.ctmin  = 0.004; //Minimum absolute value of CT used during optimization.
  opt.conminData.ctl    = -0.01; //Constraint thickness parameter for linear and side constraints.
  opt.conminData.ctlmin = 0.001; //Minimum value of CTL used during optimization.
  opt.conminData.delfun = 0.001; //Relative convergence criterion threshold, Threshold for the minimum relative change in the objective function
  opt.conminData.dabfun = 0.001; //Absolute convergence criterion threshold. Threshold for the minimum relative change in the objective function
  opt.conminData.nside  = 1; //side constraints parameter
  opt.conminData.itrm   = 3; //diminishing return criterion iteration number
  opt.conminData.icndir = numTheta+1; //conjugate direction restart parameter
}

void KrigingModel::set_direct_parameters(OptimizationProblem& opt) const
{
  opt.directData.minBoxSize = -1.0;
  opt.directData.volBoxSize = -1.0;
  //opt.directData.minBoxSize = 1.0e-15;
  //opt.directData.volBoxSize = 1.0e-15;
  //opt.directData.minBoxSize = 1.0e-3;
  //opt.directData.volBoxSize = 1.0e-5;
  opt.directData.solutionTarget = -DBL_MAX;
  opt.directData.convergenceTol = 1.0e-4;
  opt.directData.maxFunctionEvals = maxTrials;
  opt.directData.maxIterations = 1000;
  opt.directData.verboseOutput = false;
  opt.directData.constraintsPresent = true;
}

} // end namespace nkm
