package beast.evolution.operators;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.DPPointer;
import beast.core.parameter.DPValuable;
import beast.core.parameter.RealParameter;
import beast.math.util.MathUtil;
import beast.util.Randomizer;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Chieh-Hsi Wu
 */
@Description("This class re-assigns a random number of sites in a randomly picked category to another randomly picked category.")
public class CategoryDeltaExchangeOperator extends Operator {
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


        //Pick the number of sites to be moved
        int size = Randomizer.nextInt((int)Math.round(delta));

        //The number of the sites to be moved is greater or equal to number of sites in the clusters.
        //As a the (log) prior probability on this clustering structure is 0 (-infinity).
        if(size >= dpVal.getClusterSize(icat)){
            return Double.NEGATIVE_INFINITY;
        }


        int[] sites = dpVal.getClusterSites(icat);
        /*System.out.print("cati: ");
        for(int site:sites){
            System.out.print(site+" ");
        }
        System.out.println();
        System.out.print("catj: ");
        int[] sitesj = dpVal.getClusterSites(jcat);
        for(int sitej:sitesj){
            System.out.print(sitej+" ");
        }
        System.out.println();*/

        int[] sitesToBeMoved = MathUtil.sample(size, sites, false);
        /*System.out.print("sitesToBeMoved: ");

        for(int siteToBeMoved:sitesToBeMoved){
            System.out.print(siteToBeMoved+" ");
        }
        System.out.println();  */
        int sitej = dpVal.getOneClusterSite(jcat);
        for(DPPointer pointer: pointers){
            pointer.multiPointerChanges(sitesToBeMoved, sitej);
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
