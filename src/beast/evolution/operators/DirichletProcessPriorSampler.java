package beast.evolution.operators;

import beast.core.Operator;
import beast.core.Input;
import beast.core.Description;
import beast.core.parameter.*;
import beast.math.distributions.DirichletProcess;
import beast.math.distributions.ParametricDistribution;
import beast.util.Randomizer;

/**
 * @author Chieh-Hsi Wu
 *
 */
@Description("Performs gibbs sampling when a dirichlet process prior is used.")
public class DirichletProcessPriorSampler extends Operator {
    //assignment
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

    /*public Input<Distribution> likelihoodInput = new Input<Distribution>(
            "likelihood",
            "The likelihood given the data",
            Input.Validate.REQUIRED
    );*/

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
        //Get the pointer and the list of unique values
        DPPointer pointers = pointersInput.get(this);
        ParameterList paramList = xListInput.get(this);

        //Randomly pick an index to update, gets it's current value and its position in the parameter list
        int dimPointer = pointers.getDimension();
        int index = Randomizer.nextInt(dimPointer);
        RealParameter currVal = paramList.getParameter(pointers.indexInList(index,paramList));
        int listIndex = paramList.indexOf(currVal);


        //Distribution lik = likelihoodInput.get();

        //Get the dimension of the parameter
        int dimValue = paramList.getParameterDimension();

        //Count the number of items in each cluster but excluding the one about to be updated
        int[] clusterCounts = dpVal.getClusterCounts();
        clusterCounts[listIndex] =  clusterCounts[listIndex]-1;

        try{

            //Generate a sample of proposals
            QuietRealParameter[] preliminaryProposals = sampleFromBaseDistribution(dimValue);
            int dimList = paramList.getDimension();
            int i;
            double concVal =dp.getConcParameter();
            
            double[] fullConditional = new double[dimList+sampleSize];

            for(i = 0; i < dimList; i++){
                //n_{-index,i}/(n - 1 + alpha)
                fullConditional[i] = clusterCounts[i]/(dimPointer - 1 + concVal);
                //System.err.println("fullConditional[i]: "+fullConditional[i]+", val: "+paramList.getParameter(i));

            }
            for(; i < fullConditional.length; i++){
                //alpha/m/(n - 1 + alpha)
                fullConditional[i] = concVal/sampleSize/(dimPointer - 1 + concVal);
                //System.err.println("fullConditional[i]: "+fullConditional[i]+", val: "+i);
            }

            int proposedIndex = Randomizer.randomChoicePDF(fullConditional);
            //System.err.println("proposedIndex: "+proposedIndex);
            if(proposedIndex < dimList){
                //take up an existing value
                pointers.point(index, paramList.getParameter(proposedIndex));
            }else{
                //take up a new value
                pointers.point(index, preliminaryProposals[proposedIndex-dimList]);
                paramList.addParameter(preliminaryProposals[proposedIndex-dimList]);
            }

            /*for(i = 0; i < clusterCounts.length;i++){
                System.err.println(paramList.getParameter(i).getValue()+", cluster counts: "+clusterCounts[i]);
            }*/

            //If any cluster has no member then it is removed.
            for(i = 0; i < clusterCounts.length;i++){
                if(clusterCounts[i] ==0){
                    paramList.removeParameter(i);
                    break;
                }

            }

        }catch(Exception e){
                throw new RuntimeException(e);
            }
        return Double.POSITIVE_INFINITY;
    }

    public QuietRealParameter[] sampleFromBaseDistribution(int dimValue) throws Exception{
        QuietRealParameter[] preliminaryProposals = new QuietRealParameter[sampleSize];

        for(int i = 0; i < sampleSize; i++){

            Double[] sample = new Double[dimValue];
            for(int j = 0;j < dimValue;j++){

                sample[j] = baseDistr.inverseCumulativeProbability(Randomizer.nextDouble());
                //System.err.print(sample[j]+" ");
            }
            preliminaryProposals[i] = new QuietRealParameter(sample);

        }
        //System.err.println();
        return preliminaryProposals;

    }
}
