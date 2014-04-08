package beast.evolution.operators;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.DPPointer;
import beast.core.parameter.DPValuable;
import beast.math.util.MathUtil;
import beast.util.Randomizer;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Chieh-Hsi Wu
 */
@Description("This operator random selections a x sites from two clusters and switches the assignment.")
public class MultiPointerSwapOperator extends Operator {
    public Input<Double> deltaInput = new Input<Double>(
            "delta",
            "Magnitude of change for two randomly picked categories.",
            Input.Validate.REQUIRED
    );

    public Input<DPValuable> dpValInput = new Input<DPValuable>(
            "dpVal",
            "The Dirichlet process valuable object that records information on the clusters.",
            Input.Validate.REQUIRED
    );


    public Input<List<DPPointer>> pointersInput = new Input<List<DPPointer>>(
            "pointers",
            "The Dirichlet process valuable object that records information on the clusters.",
            new ArrayList<DPPointer>(),
            Input.Validate.REQUIRED
    );

    public final Input<Boolean> input_autoOptimize = new Input<Boolean>(
            "autoOptimize",
            "if true, window size will be adjusted during the MCMC run to improve mixing.",
            true
    );


    private double delta;
    private DPValuable dpVal;
    private List<DPPointer> pointers;
    private boolean autoOptimize;

    public void initAndValidate(){
        delta = deltaInput.get();
        dpVal = dpValInput.get();
        pointers = pointersInput.get();
        autoOptimize = input_autoOptimize.get();

    }

    public double proposal(){
        int categoryCount = dpVal.getCategoryCount();
        if(categoryCount <=1){
            return Double.NEGATIVE_INFINITY;
        }
        //Randomly pick two categories
        int icat = Randomizer.nextInt(categoryCount);
        int jcat = 0;
        while(icat == jcat){
            jcat = Randomizer.nextInt(categoryCount);
        }


        int[] isites = dpVal.getClusterSites(icat);
        int[] jsites = dpVal.getClusterSites(jcat);

        //Pick the number of sites to be moved
        int size = Randomizer.nextInt((int)Math.round(delta));



        //The number of the sites to be moved is greater or equal to number of sites in the clusters.
        //As a the (log) prior probability on this clustering structure is 0 (-infinity).
        while(size > isites.length + jsites.length){
            size = Randomizer.nextInt((int)Math.round(delta));
        }

        int[] sites = new int[isites.length + jsites.length];
        System.arraycopy(isites, 0, sites, 0, isites.length);
        System.arraycopy(jsites, 0, sites, isites.length, jsites.length);



        int[] siteIndicesToBeMoved = new int[size];
        int[] sitesToBeMoved = MathUtil.sample(size, sites, false, siteIndicesToBeMoved);
        int[] desPointerIndex = new int[siteIndicesToBeMoved.length];

        int siteFromITOJCount = 0;
        int siteFromJTOICount = 0;
        for(int i = 0; i < desPointerIndex.length; i++){

            if(siteIndicesToBeMoved[i] < isites.length){
                desPointerIndex[i] = jsites[0];
                siteFromITOJCount++;
            }else{
                desPointerIndex[i] = isites[0];
                siteFromJTOICount++;
            }
        }
        if(siteFromJTOICount >= jsites.length || siteFromITOJCount >= isites.length ){
            return Double.NEGATIVE_INFINITY;
        }
        //System.out.println();




        for(DPPointer pointer: pointers){
            pointer.multiPointerChanges(sitesToBeMoved, desPointerIndex);
        }

       return 0.0;

    }

    @Override
    public double getCoercableParameterValue() {
        return delta;
    }

    @Override
    public void setCoercableParameterValue(final double fValue) {
        delta = fValue;
    }

    /**
     * called after every invocation of this operator to see whether
     * a parameter can be optimised for better acceptance hence faster
     * mixing
     *
     * @param logAlpha difference in posterior between previous state & proposed state + hasting ratio
     */
    @Override
    public void optimize(final double logAlpha) {
        // must be overridden by operator implementation to have an effect
        if (autoOptimize) {
            double fDelta = calcDelta(logAlpha);
            fDelta += Math.log(delta);
            delta = Math.exp(fDelta);
        }

    }

    @Override
    public final String getPerformanceSuggestion() {
        final double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        final double targetProb = getTargetAcceptanceProbability();

        double ratio = prob / targetProb;
        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        // new scale factor
        final double newDelta = delta * ratio;

        final DecimalFormat formatter = new DecimalFormat("#.###");
        if (prob < 0.10) {
            return "Try setting delta to about " + formatter.format(newDelta);
        } else if (prob > 0.40) {
            return "Try setting delta to about " + formatter.format(newDelta);
        } else return "";
    }
}
