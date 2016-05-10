package beast.evolution.operators;

import beast.core.parameter.*;
import beast.core.*;
import beast.evolution.likelihood.*;
import beast.math.distributions.ParametricDistribution;
import beast.math.distributions.CompoundDirichletProcess;
import beast.util.Randomizer;
import beast.evolution.sitemodel.DPNtdRateSepSiteModel;

import java.io.PrintStream;
import java.util.List;

/**
 * @author Chieh-Hsi Wu
 */
@Description("Gibbs sampler with DPP for NtdBMA.")
public class NtdBMADPPGibbsSampler extends Operator implements Loggable {
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

    public Input<DPPointer> freqsPointersInput = new Input<DPPointer>(
            "freqPointers",
            "array which points a set of unique parameter values",
            Input.Validate.REQUIRED
    );
    public Input<ParameterList> freqsListInput = new Input<ParameterList>(
            "freqsList",
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

    public Input<Boolean> testCorrectInput = new Input<Boolean>(
            "testCorrect",
            "Whether to check whether the likelihood calculation is consistent.",
            false
    );


    private CompoundDirichletProcess dp;
    private int sampleSize;
    private ParametricDistribution paramBaseDistr;
    private ParametricDistribution modelBaseDistr;
    private ParametricDistribution freqsBaseDistr;
    private DPValuable dpVal;
    private DPTreeLikelihood dpTreeLikelihood;
    private int addClusterCount =0;
    private boolean testCorrect;
    public void initAndValidate(){
        testCorrect = testCorrectInput.get();
        dp = dpInput.get();
        List<ParametricDistribution> distrs = dp.getBaseDistributions();
        paramBaseDistr = distrs.get(0);
        modelBaseDistr = distrs.get(1);
        freqsBaseDistr = distrs.get(2);
        
        sampleSize = sampleSizeInput.get();
        dpVal = dpValuableInput.get();
        dpTreeLikelihood = dpTreeLikelihoodInput.get();

    }


    public double proposal(){
        //System.out.println("NtdBMADPPGibbsSampler proposing! "+getID());
        //Get the pointer and the list of unique values
        //DPPointer paramPointers = parameterPointersInput.get(this);
        //DPPointer freqPointers = freqPointersInput.get(this);
        //DPPointer modelPointers = modelPointersInput.get(this);
        DPPointer paramPointers = parameterPointersInput.get();
        DPPointer freqPointers = freqsPointersInput.get();
        DPPointer modelPointers = modelPointersInput.get();

        ParameterList paramList = parameterListInput.get();
        ParameterList freqsList = freqsListInput.get();
        ParameterList modelList = modelListInput.get();


        //Randomly pick an index to update, gets it's current value and its position in the parameter list
        int dimPointer = paramPointers.getDimension();
        int index = Randomizer.nextInt(dimPointer);



        int listIndex = paramPointers.indexInList(index,paramList);
        RealParameter currParamVal = paramList.getParameter(listIndex);
        RealParameter currModelVal = modelList.getParameter(listIndex);
        RealParameter currFreqsVal = freqsList.getParameter(listIndex);




        TempTreeLikelihood tempLik = tempLikelihoodInput.get();

        //Get the dimension of the parameter
        int paramDimValue = paramList.getParameterDimension();
        int modelDimValue = modelList.getParameterDimension();
        int freqDimValue = freqsList.getParameterDimension();

        RealParameter curr = freqPointers.getParameter(index);

        //Count the number of items in each cluster but excluding the one about to be updated
        int[] clusterCounts = dpVal.getClusterCounts();
        clusterCounts[listIndex] =  clusterCounts[listIndex]-1;

        QuietRealParameter[] existingParamVals = new QuietRealParameter[clusterCounts.length];
        QuietRealParameter[] existingModelVals = new QuietRealParameter[clusterCounts.length];
        QuietRealParameter[] existingFreqsVals = new QuietRealParameter[clusterCounts.length];
        int[] existingCluster = new int[clusterCounts.length];

        int counter = 0;
        int zeroCount = -1;

        for(int i = 0; i < clusterCounts.length;i++){
            if(clusterCounts[i]>0){
                clusterCounts[counter] = clusterCounts[i];
                existingParamVals[counter] = paramList.getParameter(i);
                existingModelVals[counter] = modelList.getParameter(i);
                existingFreqsVals[counter] = freqsList.getParameter(i);
                existingCluster[counter] = i;
                counter++;

            }else{
                zeroCount = i;
            }
        }

        try{

            //Generate a sample of proposals
            QuietRealParameter[] paramPreProposals = getSamples(paramBaseDistr,currParamVal);
            QuietRealParameter[] modelPreProposals = getSamples(modelBaseDistr,currModelVal);
            QuietRealParameter[] freqsPreProposals = getSamples(freqsBaseDistr,currFreqsVal);

            //System.err.println("zero count: "+zeroCount);
            //If the a singleton has been picked
            if(zeroCount > -1){
                paramPreProposals[0] = paramList.getParameter(zeroCount);
                modelPreProposals[0] = modelList.getParameter(zeroCount);
                freqsPreProposals[0] = freqsList.getParameter(zeroCount);

            }
            //int dimList = paramList.getDimension();
            int i;
            double concVal =dp.getConcParameter();

            double[] logFullCond = new double[counter+sampleSize];

            for(i = 0; i < counter; i++){

                logFullCond[i] = Math.log(clusterCounts[i]/(dimPointer - 1 + concVal));

                //logFullCond[i] = logFullCond[i]+ dpTreeLikelihood.getSiteLogLikelihood(existingCluster[i],index);

                if(testCorrect){
                    double temp1 = ((SepTempTreeLikelihood)tempLik).calculateLogP(
                                    existingParamVals[i],
                            existingModelVals[i],
                            existingFreqsVals[i],
                                    index
                            );
                    double temp2 = getSiteLogLikelihood(
                            existingParamVals[i],
                            existingModelVals[i],
                            existingFreqsVals[i],
                            existingCluster[i],
                            index,
                            tempLik
                    );

                    if(temp1 != temp2){
                        System.out.println(existingParamVals[i]);
                        System.out.println(existingModelVals[i]);
                        System.out.println(existingFreqsVals[i]);
                        //((DPSepTreeLikelihood)dpTreeLikelihood).printThings();
                        throw new RuntimeException(temp1+" "+temp2);
                    }
                }
                //logFullCond[i] = logFullCond[i]+temp2;
                logFullCond[i] = logFullCond[i]+getSiteLogLikelihood(
                        existingParamVals[i],
                        existingModelVals[i],
                        existingFreqsVals[i],
                        existingCluster[i],
                        index,
                        tempLik
                );
            }

            if(zeroCount > -1){
                logFullCond[i] = Math.log(concVal/sampleSize/(dimPointer - 1 + concVal));

                if(testCorrect){
                    double temp1 = tempLik.calculateLogP(
                                    paramPreProposals[0],
                                    modelPreProposals[0],
                                    freqsPreProposals[0],
                                    index
                            );
                    double temp2 = getSiteLogLikelihood(
                            paramPreProposals[0],
                            modelPreProposals[0],
                            freqsPreProposals[0],
                            zeroCount,
                            index,
                            tempLik
                    );
                    if(temp1 != temp2){
                        throw new RuntimeException(temp1+" "+temp2);
                    }
                }
                //logFullCond[i] = logFullCond[i]+temp2;
                logFullCond[i] = logFullCond[i]+getSiteLogLikelihood(
                        paramPreProposals[0],
                        modelPreProposals[0],
                        freqsPreProposals[0],
                        zeroCount,
                        index,
                        tempLik
                );
                //logFullCond[i] = logFullCond[i]+ dpTreeLikelihood.getSiteLogLikelihood(zeroCount,index);
                i++;
            }
            for(; i < logFullCond.length; i++){

                logFullCond[i] = Math.log(concVal/sampleSize/(dimPointer - 1 + concVal));

                logFullCond[i] = logFullCond[i]+tempLik.calculateLogP(
                                paramPreProposals[i-counter],
                                modelPreProposals[i-counter],
                                freqsPreProposals[i-counter],
                                index
                );
                if(Double.isNaN(logFullCond[i])){

                    System.err.println(logFullCond[i]+" "+logFullCond[i]);
                    System.err.println(paramPreProposals[i-counter]);
                    System.err.println(modelPreProposals[i-counter]);
                    System.err.println(freqsPreProposals[i-counter]);
                    dpTreeLikelihoodInput.get().getTree().log(-1,System.err);
                    //return Double.NEGATIVE_INFINITY;
                    logFullCond[i] = Double.NEGATIVE_INFINITY;

                }

            }

            /*double smallestVal = logFullCond[0];
            for(i = 1; i < logFullCond.length;i++){
                if(smallestVal > logFullCond[i])
                    smallestVal = logFullCond[i];
            }*/


            double[] fullConditional = new double[logFullCond.length];
            //double sum = 0.0;
            for(i = 0; i < fullConditional.length;i++){
                fullConditional[i] = Math.exp(logFullCond[i]);
                //-smallestVal);
                //System.err.println("fullConditional[i]: "+fullConditional[i]);
                if(Double.isNaN(fullConditional[i])||fullConditional[i]==Double.POSITIVE_INFINITY){
                    //System.err.println(smallestVal);
                    for(int j = 0;j <logFullCond.length;j++){
                        System.err.println("logFullCond[i]: "+logFullCond[j]);

                    }
                    for(int j = 0; j < existingParamVals.length;j++){
                        System.err.println(existingParamVals[j]+" "+existingModelVals[j]+" "+existingFreqsVals[j]);

                    }
                    for(int j = 0; j < paramPreProposals.length;j++){
                        System.err.println(paramPreProposals[j]+" "+modelPreProposals[j]+" "+freqsPreProposals[j]);

                    }
                    //return Double.NEGATIVE_INFINITY;
                }

                //sum+=fullConditional[i];
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




            int proposedIndex = Randomizer.randomChoicePDF(fullConditional);
            //System.err.println("proposedIndex: "+proposedIndex);
            QuietRealParameter paramProposal;
            QuietRealParameter modelProposal;
            QuietRealParameter freqsProposal;

            if(proposedIndex < counter){
                paramProposal = existingParamVals[proposedIndex];
                modelProposal = existingModelVals[proposedIndex];
                freqsProposal = existingFreqsVals[proposedIndex];
            }else{
                paramProposal = paramPreProposals[proposedIndex-counter];
                modelProposal = modelPreProposals[proposedIndex-counter];
                freqsProposal = freqsPreProposals[proposedIndex-counter];

            }

            if(curr != freqsProposal){//Takes a different value from the current

                if(proposedIndex >= counter && zeroCount > -1){//Singleton takes new value
                    //System.out.println(getID()+"value changed");

                    paramList = parameterListInput.get(this);
                    modelList = modelListInput.get(this);
                    freqsList = freqsListInput.get(this);

                    int paramListIndex = paramPointers.indexInList(index,paramList);
                    //System.err.println("paramListIndex: "+paramListIndex);
                    for(i = 0; i < paramDimValue;i++){
                        paramList.setValue(paramListIndex,i,paramProposal.getValue(i));
                    }

                    modelList.setValue(paramListIndex,0,modelProposal.getValue());

                    for(i = 0; i < freqDimValue; i++){
                        freqsList.setValue(paramListIndex,i,freqsProposal.getValue(i));
                    }

                    zeroCount = -1;

                }else{
                    //Singleton takes existing value or
                    //non-singleton takes new or existing value
                    paramPointers = parameterPointersInput.get(this);
                    freqPointers = freqsPointersInput.get(this);
                    modelPointers = modelPointersInput.get(this);
                    paramPointers.point(index, paramProposal);
                    modelPointers.point(index, modelProposal);
                    freqPointers.point(index, freqsProposal);
                    //System.out.println(getID()+"pointer changed");

                    //Non singleton takes new value
                    if(proposedIndex >= counter){
                        //addClusterCount++;
                        paramList = parameterListInput.get(this);
                        modelList = modelListInput.get(this);
                        freqsList = freqsListInput.get(this);
                        paramList.addParameter(paramProposal);
                        modelList.addParameter(modelProposal);
                        freqsList.addParameter(freqsProposal);

                        /*for(i = 0; i < logFullCond.length;i++){
                            System.out.println(i+" "+ logFullCond[i]+" "+proposedIndex+" U: "+u);
                        }
                        

                        System.out.println(getID()+"added"+" "+sum);*/
                    }

                }

                //If any cluster has no member then it is removed.
                if(zeroCount > -1){
                    paramList = parameterListInput.get(this);
                    modelList = modelListInput.get(this);
                    freqsList = freqsListInput.get(this);
                    paramList.removeParameter(zeroCount);
                    modelList.removeParameter(zeroCount);
                    freqsList.removeParameter(zeroCount);
                    //System.out.println(getID()+"removed U "+u);
                    /*for(i = 1; i < logFullCond.length;i++){
                            System.out.println(i+" "+ fullConditional[i]);
                        }*/

                }

            }/*else{
                System.out.println(getID()+"No change");
            }*/

            
        }catch(Exception e){

                throw new RuntimeException(e);
            }
        return Double.POSITIVE_INFINITY;
    }


    public QuietRealParameter[] getSamples(ParametricDistribution distr, RealParameter example) throws Exception{
        QuietRealParameter[] samples = new QuietRealParameter[sampleSize];
        Double[][] sampleVals = distr.sample(sampleSize);
        for(int i = 0; i < samples.length;i++){
            samples[i] = new QuietRealParameter(sampleVals[i]);
            samples[i].setUpper(example.getUpper());
            samples[i].setLower(example.getLower());
        }
        return samples;
    }

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

    public double getSiteLogLikelihood(
            QuietRealParameter param,
            QuietRealParameter model,
            QuietRealParameter freqs,
            int clusterIndex,
            int siteIndex,
            TempTreeLikelihood tempLik){

        double siteLogLik = Double.NaN;
        if(dpTreeLikelihood instanceof DPSepTreeLikelihood){
            //System.out.println("flag 1");
            siteLogLik = ((DPSepTreeLikelihood)dpTreeLikelihood).getSiteLogLikelihood(
                    DPNtdRateSepSiteModel.NTDBMA,
                    freqs.getIDNumber(),
                    siteIndex
            );

            if(Double.isNaN(siteLogLik)){
                //System.out.println("flag 2");
                siteLogLik =  tempLik.calculateLogP(
                        param,
                        model,
                        freqs,
                        siteIndex
                );//todo need a new temp likelihood

            }
        }else if(dpTreeLikelihood instanceof SlowDPSepTreeLikelihood){
            siteLogLik = ((SlowDPSepTreeLikelihood)dpTreeLikelihood).getSiteLogLikelihood(
                    DPNtdRateSepSiteModel.NTDBMA,
                    freqs.getIDNumber(),
                    siteIndex
            );

            if(Double.isNaN(siteLogLik)){
                siteLogLik =  tempLik.calculateLogP(
                        param,
                        model,
                        freqs,
                        siteIndex
                );//todo need a new temp likelihood

            }
        }else{
            //System.out.println("flag 2");
            siteLogLik =  dpTreeLikelihood.getSiteLogLikelihood(clusterIndex,siteIndex);
            if(Double.isNaN(siteLogLik)){
                //System.out.println("flag 2");
                siteLogLik =  tempLik.calculateLogP(
                        param,
                        model,
                        freqs,
                        siteIndex
                );//todo need a new temp likelihood

            }
        }

        return siteLogLik;
    }

    public static void main(String[] args){
        System.out.println(Math.exp(Double.NEGATIVE_INFINITY));
    }
}