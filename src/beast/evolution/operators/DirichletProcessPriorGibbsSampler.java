package beast.evolution.operators;

import beast.core.Operator;
import beast.core.Input;
import beast.core.Distribution;
import beast.core.Description;
import beast.core.parameter.*;
import beast.math.distributions.DirichletProcess;
import beast.math.distributions.ParametricDistribution;
import beast.util.Randomizer;
import org.apache.commons.math.MathException;


/**
 * @author Chieh-Hsi Wu
 */
@Description("Gibbs sampler with a Dirichlet process prior.")
public class DirichletProcessPriorGibbsSampler extends Operator {
    public Input<DPPointer> pointersInput = new Input<DPPointer>(
            "pointers",
            "array which points a set of unique parameter values",
            Input.Validate.REQUIRED
    );
    public Input<ParameterList> xListInput = new Input<ParameterList>(
            "xList",
            "points at which the density is calculated",
            Input.Validate.REQUIRED
    );

    public Input<DirichletProcess> dpInput = new Input<DirichletProcess>(
            "dirichletProcess",
            "An object of Dirichlet Process",
            Input.Validate.REQUIRED
    );

    public Input<Integer> sampleSizeInput = new Input<Integer>(
            "sampleSize",
            "The number of prelimiary proposals",
            Input.Validate.REQUIRED
    );

    public Input<Distribution> likelihoodInput = new Input<Distribution>(
            "likelihood",
            "The likelihood given the data",
            Input.Validate.REQUIRED
    );

    public Input<DPValuable> dpValuableInput = new Input<DPValuable>(
            "dpVal",
            "reports the counts in each cluster",
            Input.Validate.REQUIRED
    );

    private DirichletProcess dp;
    private int sampleSize;
    private ParametricDistribution baseDistr;
    private DPValuable dpVal;
    public void initAndValidate(){
        dp = dpInput.get();
        sampleSize = sampleSizeInput.get();
        baseDistr = dp.getBaseDistribution();
        dpVal = dpValuableInput.get();

    }
    public double proposal(){
        //System.err.println("DPPGBS");

        //Get the pointer and the list of unique values
        DPPointer pointers = pointersInput.get(this);
        ParameterList paramList = xListInput.get();

        //Randomly pick an index to update, gets it's current value and its position in the parameter list
        int dimPointer = pointers.getDimension();
        int index = Randomizer.nextInt(dimPointer);
        //System.err.println("index: "+index+" "+pointers.getParameterValue(0)+" "+pointers.getParameterValue(1));
        RealParameter currVal = paramList.getParameter(pointers.indexInList(index,paramList));
        int listIndex = paramList.indexOf(currVal);


        Distribution lik = likelihoodInput.get();

        //Get the dimension of the parameter
        int dimValue = paramList.getParameterDimension();

        //Count the number of items in each cluster but excluding the one about to be updated
        int[] clusterCounts = dpVal.getClusterCounts();
        clusterCounts[listIndex] =  clusterCounts[listIndex]-1;

        QuietRealParameter[] existingVals = new QuietRealParameter[clusterCounts.length];
        int counter = 0;
        int zeroCount = -1;
        for(int i = 0; i < clusterCounts.length;i++){
            if(clusterCounts[i]>0){
                clusterCounts[counter] = clusterCounts[i];
                existingVals[counter++] = paramList.getParameter(i);
                //System.err.println(i+", cluster.counts: "+clusterCounts[i]);


            }else{
                zeroCount = i;

            }
        }
        //System.err.println("paramList size1: "+paramList.getDimension());
        try{

            //Generate a sample of proposals
            QuietRealParameter[] preliminaryProposals = sampleFromBaseDistribution(paramList);
            //System.err.println("zero count: "+zeroCount);
            //If the a singleton has been picked
            if(zeroCount > -1){
                preliminaryProposals[0] = paramList.getParameter(zeroCount);

            }
            //int dimList = paramList.getDimension();
            int i;
            double concVal =dp.getConcParameter();

            double[] logFullCond = new double[counter+sampleSize];
            //System.err.println("pointers: "+pointers);
            //System.err.println("counter: "+counter);
            //System.err.println();
            for(i = 0; i < counter; i++){
                //n_{-index,i}/(n - 1 + alpha)
                logFullCond[i] = Math.log(clusterCounts[i]/(dimPointer - 1 + concVal));
                pointers.point(index, existingVals[i]);
                //System.err.println("clusterCounts[i]: "+clusterCounts[i]);
                //System.err.println("lgc and lik: "+logFullCond[i]+" "+lik.calculateLogP());
                logFullCond[i] = logFullCond[i]+lik.calculateLogP();
                //System.err.println(logFullCond[i]);
                //System.err.println("fullConditional[i]: "+fullConditional[i]+", val: "+paramList.getParameter(i));
                //System.err.println("logFullCond[i]: "+logFullCond[i]+", val: "+existingVals[i]);
            }
            for(; i < logFullCond.length; i++){
                //alpha/m/(n - 1 + alpha)
                logFullCond[i] = Math.log(concVal/sampleSize/(dimPointer - 1 + concVal));
                //System.err.println("lgc and lik: "+logFullCond[i]+" "+lik.calculateLogP());
                pointers.point(index, preliminaryProposals[i-counter]);
                logFullCond[i] = logFullCond[i]+lik.calculateLogP();
                //System.err.println("lgc and lik: "+logFullCond[i]+" "+lik.calculateLogP());
                //System.err.println("logFullCond[i]: "+logFullCond[i]+", val: "+preliminaryProposals[i-counter]);
            }

            double smallestVal = logFullCond[0];
            for(i = 1; i < logFullCond.length;i++){
               // if(smallestVal == Double.NEGATIVE_INFINITY ||
               //         logFullCond[i] < smallestVal & logFullCond[i]!=Double.NEGATIVE_INFINITY){
                if(smallestVal > logFullCond[i])
                    smallestVal = logFullCond[i];
                //}
            }
            //System.err.println("smallestVal2: "+smallestVal);

            double[] fullConditional = new double[logFullCond.length];
            for(i = 0; i < fullConditional.length;i++){
                //if(logFullCond[i] == Double.NEGATIVE_INFINITY){
                //    fullConditional[i] = 0.0;
                //}else{
                //if(this.counter>3240)
                //    System.err.println("logfc: "+logFullCond[i]+" "+i);
                    fullConditional[i] = Math.exp(logFullCond[i]-smallestVal);
                //}
                //System.err.println("fc: "+fullConditional[i]);
            }


            int proposedIndex = Randomizer.randomChoicePDF(fullConditional);
            //System.err.println("proposedIndex: "+proposedIndex);
            if(proposedIndex < counter){
                //take up an existing value
                pointers.point(index, existingVals[proposedIndex]);
                //if(this.counter>6000)
                    //System.err.println("proposedIndex: "+proposedIndex+"; existingVal: "+existingVals[proposedIndex]);
            }else{
                RealParameter proposal = preliminaryProposals[proposedIndex-counter];
                if(zeroCount > -1){
                    int paramListIndex = pointers.indexInList(index,paramList);

                    for(i = 0; i < dimValue;i++){
                        paramList.setValue(paramListIndex,i,proposal.getValue(i));


                    }
                    zeroCount = -1;

                }else{
                //take up a new value
                    pointers.point(index, preliminaryProposals[proposedIndex-counter]);
                    paramList = xListInput.get(this);
                    paramList.addParameter(preliminaryProposals[proposedIndex-counter]);
                    //System.err.println("add cluster");
                }
                //System.err.println("paramList size2: "+paramList.getDimension());
                //if(this.counter>6000)
                    //System.err.println("proposedIndex: "+proposedIndex+"; newVal: "+preliminaryProposals[proposedIndex-counter]);
            }

            /*for(i = 0; i < clusterCounts.length;i++){
                System.err.println(paramList.getParameter(i).getValue()+", cluster counts: "+clusterCounts[i]);
            }*/

            //If any cluster has no member then it is removed.
            if(zeroCount > -1){
                paramList = xListInput.get(this);
                paramList.removeParameter(zeroCount);
                //if(this.counter>3240)
                  //System.err.println("Removing cluster");

            }

            //System.err.println(counter);

        }catch(Exception e){
                throw new RuntimeException(e);
            }
        return Double.POSITIVE_INFINITY;
    }

    public QuietRealParameter[] sampleFromBaseDistribution(ParameterList paramList) throws MathException {
//        QuietRealParameter[] preliminaryProposals = getSamples(baseDistr);
//        for(RealParameter parameter: preliminaryProposals){
//            parameter.setUpper(paramList.getUpper());
//            parameter.setLower(paramList.getLower());
//        }
        QuietRealParameter[] preliminaryProposals = QuietRealParameter.getSamples(baseDistr,
                sampleSize, paramList.getUpper(), paramList.getLower());

        /*RealParameter[] preliminaryProposals = new RealParameter[sampleSize];

        for(int i = 0; i < sampleSize; i++){

            Double[] sample = new Double[dimValue];
            for(int j = 0;j < dimValue;j++){

                sample[j] = baseDistr.inverseCumulativeProbability(Randomizer.nextDouble());
                double s= Randomizer.nextDouble();
                if(s<0.4){
                    sample[j] = Randomizer.nextDouble()*0.5;
                }else{
                    sample[j] = Randomizer.nextDouble()*0.5+0.5;

                }
                //System.err.print(sample[j]+" ");
            }
            preliminaryProposals[i] = new RealParameter(sample);
            preliminaryProposals[i].setUpper(paramList.getUpper());
            preliminaryProposals[i].setLower(paramList.getLower());

        }*/
        //System.err.println();
        return preliminaryProposals;

    }

//    public QuietRealParameter[] getSamples(ParametricDistribution distr) throws Exception{
//        QuietRealParameter[] samples = new QuietRealParameter[sampleSize];
//        Double[][] sampleVals = distr.sample(sampleSize);
//        for(int i = 0; i < samples.length;i++){
//            samples[i] = new QuietRealParameter(sampleVals[i]);
//        }
//        return samples;
//    }

}
