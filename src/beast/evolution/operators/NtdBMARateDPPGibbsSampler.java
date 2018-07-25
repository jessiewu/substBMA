package beast.evolution.operators;

import beast.core.Description;
import beast.core.Input;
import beast.core.Loggable;
import beast.core.Operator;
import beast.core.parameter.*;
import beast.evolution.likelihood.DPTreeLikelihood;
import beast.evolution.likelihood.TempTreeLikelihood;
import beast.math.distributions.CompoundDirichletProcess;
import beast.math.distributions.ParametricDistribution;
import beast.util.Randomizer;

import java.io.PrintStream;
import java.util.List;

/**
 * @author Chieh-Hsi Wu
 */
@Description("Gibbs sampler with DPP for NtdBMA and rate having the same partition scheme.")
public class NtdBMARateDPPGibbsSampler  extends Operator implements Loggable {
    public Input<DPPointer> parameterPointersInput = new Input<DPPointer>(
            "parameterPointers",
            "array which points a set of unique parameter values",
            Input.Validate.REQUIRED
    );
    public Input<ParameterList> parameterListInput = new Input<ParameterList>(
            "parameterList",
            "points at which the density is calculated",
            Input.Validate.REQUIRED
    );

    public Input<DPPointer> modelPointersInput = new Input<DPPointer>(
            "modelPointers",
            "array which points a set of unique parameter values",
            Input.Validate.REQUIRED
    );
    public Input<ParameterList> modelListInput = new Input<ParameterList>(
            "modelList",
            "points at which the density is calculated",
            Input.Validate.REQUIRED
    );

    public Input<DPPointer> freqPointersInput = new Input<DPPointer>(
            "freqPointers",
            "array which points a set of unique parameter values",
            Input.Validate.REQUIRED
    );
    public Input<ParameterList> freqsListInput = new Input<ParameterList>(
            "freqsList",
            "points at which the density is calculated",
            Input.Validate.REQUIRED
    );

    public Input<DPPointer> ratesPointersInput = new Input<DPPointer>(
            "ratesPointers",
            "array which points a set of unique parameter values",
            Input.Validate.REQUIRED
    );
    public Input<ParameterList> ratesListInput = new Input<ParameterList>(
            "ratesList",
            "points at which the density is calculated",
            Input.Validate.REQUIRED
    );


    public Input<CompoundDirichletProcess> dpInput = new Input<CompoundDirichletProcess>(
            "dirichletProcess",
            "An object of Dirichlet Process",
            Input.Validate.REQUIRED
    );

    public Input<Integer> sampleSizeInput = new Input<Integer>(
            "sampleSize",
            "The number of prelimiary proposals",
            Input.Validate.REQUIRED
    );


    public Input<TempTreeLikelihood> tempLikelihoodInput = new Input<TempTreeLikelihood>(
            "tempLikelihood",
            "The temporary likelihood given the data at site i",
            Input.Validate.REQUIRED
    );

    public Input<DPValuable> dpValuableInput = new Input<DPValuable>(
            "dpVal",
            "reports the counts in each cluster",
            Input.Validate.REQUIRED
    );

    public Input<DPTreeLikelihood> dpTreeLikelihoodInput = new Input<DPTreeLikelihood>(
            "dpTreeLik",
            "Tree likelihood that handle DPP",
            Input.Validate.REQUIRED
    );

    int tempCount = 1;
    private CompoundDirichletProcess dp;
    private int sampleSize;
    private ParametricDistribution paramBaseDistr;
    private ParametricDistribution modelBaseDistr;
    private ParametricDistribution freqsBaseDistr;
    private ParametricDistribution ratesBaseDistr;
    private DPValuable dpVal;
    private DPTreeLikelihood dpTreeLikelihood;
    private int addClusterCount =0;
    public void initAndValidate(){
        dp = dpInput.get();
        List<ParametricDistribution> distrs = dp.getBaseDistributions();
        paramBaseDistr = distrs.get(0);
        modelBaseDistr = distrs.get(1);
        freqsBaseDistr = distrs.get(2);
        ratesBaseDistr = distrs.get(3);

        sampleSize = sampleSizeInput.get();
        dpVal = dpValuableInput.get();
        dpTreeLikelihood = dpTreeLikelihoodInput.get();

    }


    public double proposal(){
        //Get the pointer and the list of unique values
        //DPPointer paramPointers = parameterPointersInput.get(this);
        //DPPointer freqsPointers = freqPointersInput.get(this);
        //DPPointer modelPointers = modelPointersInput.get(this);
        //DPPointer ratesPointers = ratesPointersInput.get(this);

        DPPointer paramPointers = parameterPointersInput.get();
        DPPointer freqsPointers = freqPointersInput.get();
        DPPointer modelPointers = modelPointersInput.get();
        DPPointer ratesPointers = ratesPointersInput.get();

        ParameterList paramList = parameterListInput.get();
        ParameterList modelList = modelListInput.get();
        ParameterList freqsList = freqsListInput.get();
        ParameterList ratesList = ratesListInput.get();


        //Randomly pick an index to update, gets it's current value and its position in the parameter list
        int dimPointer = paramPointers.getDimension();
        int index = Randomizer.nextInt(dimPointer);



        int listIndex = paramPointers.indexInList(index,paramList);
        RealParameter currParamVal = paramList.getParameter(listIndex);
        RealParameter currModelVal = modelList.getParameter(listIndex);
        RealParameter currFreqsVal = freqsList.getParameter(listIndex);
        RealParameter currRatesVal = ratesList.getParameter(listIndex);




        TempTreeLikelihood tempLik = tempLikelihoodInput.get();

        //Get the dimension of the parameter
        int paramDimValue = paramList.getParameterDimension();
        int modelDimValue = modelList.getParameterDimension();
        int freqsDimValue = freqsList.getParameterDimension();
        int ratesDimValue = ratesList.getParameterDimension();

        RealParameter curr = freqsPointers.getParameter(index);

        //Count the number of items in each cluster but excluding the one about to be updated
        int[] clusterCounts = dpVal.getClusterCounts();
        clusterCounts[listIndex] =  clusterCounts[listIndex]-1;

        QuietRealParameter[] existingParamVals = new QuietRealParameter[clusterCounts.length];
        QuietRealParameter[] existingModelVals = new QuietRealParameter[clusterCounts.length];
        QuietRealParameter[] existingFreqsVals = new QuietRealParameter[clusterCounts.length];
        QuietRealParameter[] existingRatesVals = new QuietRealParameter[clusterCounts.length];
        int[] existingCluster = new int[clusterCounts.length];

        int counter = 0;
        int zeroCount = -1;
        //System.out.println("clusterCounts.length: "+dpVal.getDimension());
        for(int i = 0; i < clusterCounts.length;i++){
            if(clusterCounts[i]>0){
                clusterCounts[counter] = clusterCounts[i];
                existingParamVals[counter] = paramList.getParameter(i);
                existingModelVals[counter] = modelList.getParameter(i);
                existingFreqsVals[counter] = freqsList.getParameter(i);
                existingRatesVals[counter] = ratesList.getParameter(i);
                existingCluster[counter] = i;
                counter++;

            }else{
                zeroCount = i;
            }
        }

        try{

            //Generate a sample of proposals
//            QuietRealParameter[] paramPreProposals = getSamples(paramBaseDistr,currParamVal);
//            QuietRealParameter[] modelPreProposals = getSamples(modelBaseDistr,currModelVal);
//            QuietRealParameter[] freqsPreProposals = getSamples(freqsBaseDistr,currFreqsVal);
//            QuietRealParameter[] ratesPreProposals = getSamples(ratesBaseDistr,currRatesVal);
            QuietRealParameter[] paramPreProposals = QuietRealParameter.getSamples(paramBaseDistr,
                    sampleSize, currParamVal.getUpper(), currParamVal.getLower());
            QuietRealParameter[] modelPreProposals = QuietRealParameter.getSamples(modelBaseDistr,
                    sampleSize, currModelVal.getUpper(), currModelVal.getLower());
            QuietRealParameter[] freqsPreProposals = QuietRealParameter.getSamples(freqsBaseDistr,
                    sampleSize, currFreqsVal.getUpper(), currFreqsVal.getLower());
            QuietRealParameter[] ratesPreProposals = QuietRealParameter.getSamples(ratesBaseDistr,
                    sampleSize, currRatesVal.getUpper(), currRatesVal.getLower());

            //System.err.println("zero count: "+zeroCount);
            //If the a singleton has been picked
            if(zeroCount > -1){
                paramPreProposals[0] = paramList.getParameter(zeroCount);
                modelPreProposals[0] = modelList.getParameter(zeroCount);
                freqsPreProposals[0] = freqsList.getParameter(zeroCount);
                ratesPreProposals[0] = ratesList.getParameter(zeroCount);

            }
            //int dimList = paramList.getDimension();
            int i;
            double concVal =dp.getConcParameter();

            double[] logFullCond = new double[counter+sampleSize];

            for(i = 0; i < counter; i++){

                logFullCond[i] = Math.log(clusterCounts[i]/(dimPointer - 1 + concVal));
                /*double temp1 = tempLik.calculateLogP(
                                existingParamVals[i],
                                existingModelVals[i],
                                existingFreqsVals[i],
                                existingRatesVals[i],
                                index
                        );
                double temp2 = dpTreeLikelihood.getSiteLogLikelihood(existingCluster[i],index);
                if(temp1 != temp2){
                    System.out.println("existingParamVals: "+existingParamVals[i]);
                    System.out.println("existingModelVals: "+existingModelVals[i]);
                    System.out.println("existingFreqsVals: "+existingFreqsVals[i]);
                    throw new RuntimeException(temp1+" "+temp2);
                } */
                //System.err.println("clusterCounts: "+clusterCounts[i]);
                //logFullCond[i] = logFullCond[i]+ dpTreeLikelihood.getSiteLogLikelihood(existingCluster[i],index);
                logFullCond[i] = logFullCond[i]+getSiteLogLikelihood(
                        existingParamVals[i],
                        existingModelVals[i],
                        existingFreqsVals[i],
                        existingRatesVals[i],
                        existingCluster[i],
                        index,
                        tempLik
                );
            }

            if(zeroCount > -1){
                logFullCond[i] = Math.log(concVal/sampleSize/(dimPointer - 1 + concVal));
                /*double temp1 = tempLik.calculateLogP(
                                paramPreProposals[0],
                                modelPreProposals[0],
                                freqsPreProposals[0],
                                ratesPreProposals[0],
                                index
                        );
                double temp2 = dpTreeLikelihood.getSiteLogLikelihood(zeroCount,index);
                if(temp1 != temp2){
                    throw new RuntimeException(temp1+" "+temp2);
                }*/
                //logFullCond[i] = logFullCond[i]+ dpTreeLikelihood.getSiteLogLikelihood(zeroCount,index);
                logFullCond[i] = logFullCond[i]+getSiteLogLikelihood(
                        paramPreProposals[0],
                        modelPreProposals[0],
                        freqsPreProposals[0],
                        ratesPreProposals[0],
                        zeroCount,
                        index,
                        tempLik
                );
                i++;
            }
            for(; i < logFullCond.length; i++){

                logFullCond[i] = Math.log(concVal/sampleSize/(dimPointer - 1 + concVal));
                /*double temp1 = tempLik.calculateLogP(
                                paramPreProposals[i-counter],
                                modelPreProposals[i-counter],
                                freqsPreProposals[i-counter],
                                index
                        );*/

                logFullCond[i] = logFullCond[i]+tempLik.calculateLogP(
                        paramPreProposals[i-counter],
                        modelPreProposals[i-counter],
                        freqsPreProposals[i-counter],
                        ratesPreProposals[i-counter],
                        index
                );
                if(Double.isNaN(logFullCond[i])){

                    System.err.println(logFullCond[i]+" "+logFullCond[i]);
                    System.err.println(paramPreProposals[i-counter]);
                    System.err.println(modelPreProposals[i-counter]);
                    System.err.println(freqsPreProposals[i-counter]);
                    dpTreeLikelihoodInput.get().getTree().log(-1,System.err);
                    return Double.NEGATIVE_INFINITY;

                }

            }




            //System.err.println("smallestVal2: "+smallestVal);

            double[] fullConditional = new double[logFullCond.length];

            for(i = 0; i < fullConditional.length;i++){
                fullConditional[i] = Math.exp(logFullCond[i]);

                //System.out.println("fullConditional[i]: "+fullConditional[i]+" logFullCond[i]: "+logFullCond[i]);
                if(Double.isNaN(fullConditional[i])){
                    return Double.NEGATIVE_INFINITY;
                }else if(fullConditional[i] == Double.NEGATIVE_INFINITY){

                    fullConditional[i] = 0.0;
                    
                    //return Double.NEGATIVE_INFINITY;
                    //throw new RuntimeException("Huh?");
                    //System.out.println("reject directly, existing parameter vals: "+existingParamVals.length);
                }else if(fullConditional[i] == Double.POSITIVE_INFINITY){
                    //System.out.println("accept directly, existing parameter vals: "+existingParamVals.length);
                    fullConditional[i] = 1.0;

                }
            }

            double total = Randomizer.getTotal(fullConditional);
            if(total < 1e-100){
                /*System.out.println("flag1");
                for(int j = 0; j < fullConditional.length; j++){
                    //fullConditional[j] = Math.exp(logFullCond[j] - smallestVal);
                    System.out.println(fullConditional[j]+" "+logFullCond[j]);
                }*/
                double maxVal = logFullCond[0];
                for(int j = 1; j < logFullCond.length;j++){
                    if(maxVal < logFullCond[j]) {
                        maxVal = logFullCond[j];
                    }
                }
                for(int j = 0; j < fullConditional.length; j++){
                    fullConditional[j] = Math.exp(logFullCond[j] - maxVal);
                    //System.out.println(logFullCond[j] - smallestVal);
                }

            }

            /*for(int j = 0; j < fullConditional.length;j++){
                System.out.println(j+" "+fullConditional[j]);
            }*/


            int proposedIndex = Randomizer.randomChoicePDF(fullConditional);
            //System.out.println("proposedIndex: "+proposedIndex);
            QuietRealParameter paramProposal;
            QuietRealParameter modelProposal;
            QuietRealParameter freqsProposal;
            QuietRealParameter ratesProposal;

            if(proposedIndex < counter){
                paramProposal = existingParamVals[proposedIndex];
                modelProposal = existingModelVals[proposedIndex];
                freqsProposal = existingFreqsVals[proposedIndex];
                ratesProposal = existingRatesVals[proposedIndex];
            }else{
                paramProposal = paramPreProposals[proposedIndex-counter];
                modelProposal = modelPreProposals[proposedIndex-counter];
                freqsProposal = freqsPreProposals[proposedIndex-counter];
                ratesProposal = ratesPreProposals[proposedIndex-counter];
            }

            if(curr != freqsProposal){//Takes a different value from the current
                //System.out.println("before sampler: "+ratesList.getID()+" "+ ratesList.getDimension());
                if(proposedIndex >= counter && zeroCount > -1){//Singleton takes new value

                    paramList = parameterListInput.get(this);
                    modelList = modelListInput.get(this);
                    freqsList = freqsListInput.get(this);
                    ratesList = ratesListInput.get(this);

                    int paramListIndex = paramPointers.indexInList(index,paramList);
                    //System.err.println("paramListIndex: "+paramListIndex);
                    for(i = 0; i < paramDimValue;i++){
                        paramList.setValue(paramListIndex,i,paramProposal.getValue(i));
                    }

                    modelList.setValue(paramListIndex,0,modelProposal.getValue());

                    for(i = 0; i < freqsDimValue; i++){
                        freqsList.setValue(paramListIndex,i,freqsProposal.getValue(i));
                    }
                    ratesList.setValue(paramListIndex,0,ratesProposal.getValue());
                    zeroCount = -1;

                }else{
                    //Singleton takes existing value or
                    //non-singleton takes new or existing value
                    paramPointers = parameterPointersInput.get(this);
                    freqsPointers = freqPointersInput.get(this);
                    modelPointers = modelPointersInput.get(this);
                    ratesPointers = ratesPointersInput.get(this);
                    paramPointers.point(index, paramProposal);
                    modelPointers.point(index, modelProposal);
                    freqsPointers.point(index, freqsProposal);
                    ratesPointers.point(index, ratesProposal);

                    //Non singleton takes new value
                    if(proposedIndex >= counter){
                        //addClusterCount++;
                        paramList = parameterListInput.get(this);
                        modelList = modelListInput.get(this);
                        freqsList = freqsListInput.get(this);
                        ratesList = ratesListInput.get(this);
                        paramList.addParameter(paramProposal);
                        modelList.addParameter(modelProposal);
                        freqsList.addParameter(freqsProposal);
                        ratesList.addParameter(ratesProposal);
                        //System.out.println("add "+ (++tempCount));
                        //System.out.println("sampler: "+ratesList.getID()+" "+ ratesList.getDimension());
                       //System.out.println("add "+ ratesList.getDimension());
                    }

                }

                //If any cluster has no member then it is removed.
                if(zeroCount > -1){

                    paramList = parameterListInput.get(this);
                    modelList = modelListInput.get(this);
                    freqsList = freqsListInput.get(this);
                    ratesList = ratesListInput.get(this);
                    paramList.removeParameter(zeroCount);
                    modelList.removeParameter(zeroCount);
                    freqsList.removeParameter(zeroCount);
                    ratesList.removeParameter(zeroCount);
                    //System.out.println("remove "+ (--tempCount));
                    //System.out.println("sampler: "+ratesList.getID()+" "+ ratesList.getDimension());

                }

            } /*else{
                System.out.println("No change");
            }*/


        }catch(Exception e){
                throw new RuntimeException(e);
            }
        return Double.POSITIVE_INFINITY;
    }

    public double getSiteLogLikelihood(
            RealParameter param,
            RealParameter model,
            RealParameter freqs,
            RealParameter rates,
            int clusterIndex,
            int siteIndex,
            TempTreeLikelihood tempLik){


            //System.out.println("flag 2");
            double siteLogLik =  dpTreeLikelihood.getSiteLogLikelihood(clusterIndex,siteIndex);
            if(Double.isNaN(siteLogLik)){
                //System.out.println("flag 2");
                siteLogLik =  tempLik.calculateLogP(
                        param,
                        model,
                        freqs,
                        rates,
                        siteIndex
                );//todo need a new temp likelihood

            }

        return siteLogLik;
    }

//    public QuietRealParameter[] getSamples(ParametricDistribution distr, RealParameter example) throws Exception{
//        QuietRealParameter[] samples = new QuietRealParameter[sampleSize];
//        Double[][] sampleVals = distr.sample(sampleSize);
//        for(int i = 0; i < samples.length;i++){
//            samples[i] = new QuietRealParameter(sampleVals[i]);
//            samples[i].setUpper(example.getUpper());
//            samples[i].setLower(example.getLower());
//        }
//        return samples;
//    }

    @Override
	public void init(PrintStream out) {
        out.print("count(Add cluster)\t");
    }

    @Override
	public void log(int nSample, PrintStream out) {
    	out.print((double)addClusterCount/(double)nSample+ "\t");
	}


    @Override
	public void close(PrintStream out) {
	}


}
