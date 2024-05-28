#include "MainEffectsExcelOutput.h"

    std::ostringstream ss;

    MainEffectsExcelOutput::MainEffectsExcelOutput(){;}

    MainEffectsExcelOutput::~MainEffectsExcelOutput(){;}

    std::string MainEffectsExcelOutput::outputColumnHeaders
                    (int numInputs, int numOutputs) {

        std::ostringstream ss;

         /* input variables */
         for (int i=0; i<numInputs; i++) {
         	if (ss.str()!="") ss << ",";
         	ss << "in(" << i << ")";
         }

         /* output variables */
         for (int i=0; i<numOutputs; i++) {
         	if (ss.str()!="") ss << ",";
         	ss << "out(" << i << ")";
         }

         /* number of observations */
         ss << ",nObservations";

         /* sum of all observations */
         ss << ",sumOfAllObservations";

         /* average of all observations */
         ss << ",avgOfAllObservation";

         /* sum of squares of all observations */
         ss << ",sumOfSquaresOfAllObservations";

         /* degrees of freedom of all observations */
         ss << ",degreesOfFreedomOfAllObservations";

         /* variance of all Observations */
         ss << ",varianceOfAllObservations";

         /* sum */
         ss << ",sum";

         /* average */
         ss << ",average";

         /* sum of squares */
         ss << ",sumOfSquares";

         /* variance */
         ss << ",variance";

         /* sum of squares between groups */
         ss << ",sumOfSquaresBetweenGroups";

         /* degrees of freedom between groups */
         ss << ",degreesOfFreedomBetweenGroups";

         /* variance between groups */
         ss << ",varianceBetweenGroups";

         /* sum of squares between groups */
         ss << ",sumOfSquaresWithinGroups";

         /* degrees of freedom between groups */
         ss << ",degreesOfFreedomWithinGroups";

         /* variance between groups */
         ss << ",varianceWithinGroups";

         /* F */
         ss << ",F";

         ss << "\n";

         return(ss.str());
    }


    std::string MainEffectsExcelOutput::outputMainEffects
         (int inputVarIndex, int numInputs,
          int outputVarIndex, int numOutputs,
          DDaceMainEffects::Factor factor) {

        std::ostringstream ss;
        int numberOfLevels =
            factor.getNumberOfLevels();

        for (int i=0; i<numberOfLevels; i++) {
        	ss << outputMainEffects (inputVarIndex, numInputs,
                   outputVarIndex, numOutputs, factor, i);
        }

        return(ss.str());
     }

    std::string MainEffectsExcelOutput::outputMainEffects
         (int inputVarIndex, int numInputs,
          int outputVarIndex, int numOutputs,
          DDaceMainEffects::Factor factor,
          int indexOfInputValue) {

        std::ostringstream ss;

         /* put an F under the input variable */
         for (int i=0; i<numInputs; i++) {
         	if (ss.str()!="") ss << ",";
         	if (i==inputVarIndex) ss << "F";
         }

         /* put an R under the output variable */
         for (int i=0; i<numOutputs; i++) {
         	if (ss.str()!="") ss << ",";
         	if (i==outputVarIndex) ss << "R";
         }

         /* number of observations */
         ss << ",";
         if (indexOfInputValue==0) {
             ss << factor.getNumberOfObservations();
         }


         /* sum of all observations */
         ss << ",";
         if (indexOfInputValue==0) {
             DDaceMainEffects::Response response = factor.getResponse();
             std::vector<double> responses = response.responses_;
             ss << response.getSumPop();
         }

         /* average of all Observation */
         ss << ",";
         if (indexOfInputValue==0) {
          	 DDaceMainEffects::Response response = factor.getResponse();
             ss << response.getAveragePop();
         }

         /* sum of squares for Observation */
         ss << ",";
         if (indexOfInputValue==0) {
          	 DDaceMainEffects::Response response = factor.getResponse();
             ss << response.getSumOfSquaresPop();
         }

         /* degrees of freedom for all Observation*/
         ss << ",";
         if (indexOfInputValue==0) {
             int dTotal = factor.getNumberOfObservations() - 1;
             ss << dTotal;
         }

         /* variance for all Observations */
         ss << ",";
         if (indexOfInputValue==0) {
         	DDaceMainEffects::Response response = factor.getResponse();
             ss << response.getVariancePop();
         }

         /* sum */
         ss << "," << factor.getLevelSum(indexOfInputValue);

         /* average */
         ss << "," << factor.getLevelAverage(indexOfInputValue);

         /* sum of squares */
         ss << "," << factor.getLevelSumOfSquares(indexOfInputValue);

         /* variance */
         ss << "," << factor.getLevelVariance(indexOfInputValue);

         /* sum of squares between groups */
         ss << ",";
         if (indexOfInputValue==0) {
             ss << factor.sumOfSquaresBetweenGroups();
         }

         /* degrees of freedom between groups */
         ss << ",";
         if (indexOfInputValue==0) {
             ss << factor.doFBetween();
         }

         /* variance between groups */
         ss << ",";
         if (indexOfInputValue==0) {
             ss << factor.varianceBetweenGroups();
         }

         /* sum of squares within groups */
         ss << ",";
         if (indexOfInputValue==0) {
             ss << factor.sumOfSquaresWithinGroups();
         }

         /* degrees of freedom within groups */
         ss << ",";
         if (indexOfInputValue==0) {
             ss << factor.doFWithin();
         }

         /* variance within groups */
         ss << ",";
         if (indexOfInputValue==0) {
             ss << factor.varianceWithinGroups();
         }

         /* F */
         ss << ",";
         if (indexOfInputValue==0) {
             ss << factor.Fdata();
         }


         ss << "\n";

         return(ss.str());
     }

    std::string MainEffectsExcelOutput::computeExcelOutput
        (std::vector<std::vector<double> > vectorInputData,
         std::vector<std::vector<double> > vectorOutputData){

        std::ostringstream ss;

    	/* error check */
    	if (vectorInputData.size() == 0) return("");
    	if (vectorOutputData.size() == 0) return("");


    	MainEffectsConverter converter;

         /* Replace every INPUT data value with a counting number */
         VectorCountingNumbersAndCount vectorCountingNumbersAndCount =
             converter.convertAllDoublesToCountingNumbers(vectorInputData);
         std::vector<std::vector<int> > vectorInputIndicies =
                vectorCountingNumbersAndCount.vectorCountingNumbers;
         int numberOfCountingNumbers = vectorCountingNumbersAndCount.count;

        /* How many columns are in the input table? */
        int numInputs = vectorInputData[0].size();

        /* How many columns are in the output table? */
        int numOutputs = vectorOutputData[0].size();

        /* output the column headers */
        ss << outputColumnHeaders (numInputs, numOutputs);

        /* pair input column 1 with output column 1 */
        /* pair input column 1 with output column 2 */
        /* pair input column 1 with output column 3 */
        /* etc.                                     */
        /* pair input column 2 with output column 1 */
        /* pair input column 2 with output column 2 */
        /* pair input column 2 with output column 3 */
        /* etc.                                     */
        for (int indexInput=0; indexInput<numInputs; indexInput++) {
        for (int indexOutput=0; indexOutput<numOutputs; indexOutput++) {

             /* slice out the selected input var & selected output var */
             DDaceMainEffects::Factor factor =
                 converter.sliceOutOneInputVarAndOneOutputVar
                      (vectorInputIndicies,    //data from all input vars
                       vectorOutputData, //data from all output vars
                       indexInput,             //slice out this input var
                       indexOutput,            //slice out this output var
                       numberOfCountingNumbers);  //# of different input values



             ss << outputMainEffects (indexInput, numInputs, indexOutput,
                    numOutputs, factor);
	     std::cout << ss.str() << std::endl;


        }//for indexOutput
        }//for indexInput


    	return(ss.str());
    }
