package beast.math.distributions;

import beast.core.Input;
import beast.core.Valuable;
import beast.util.Randomizer;
import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ContinuousDistribution;
import org.apache.commons.math.distribution.Distribution;
import org.apache.commons.math.distribution.IntegerDistribution;

import javax.management.RuntimeErrorException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * @aurthor Chieh-Hsi Wu
 */
public class ConditionalParametricDistribution extends ParametricDistribution  {
    public Input<List<ParametricDistribution>> distrListInput = new Input<List<ParametricDistribution>>(
            "distr",
            "Conditional distribution of a parameter",
            new ArrayList<ParametricDistribution>()

    );

    public Input<String> integerConditionValuesInput = new Input<String>(
            "integerConditionValues",
            "A list of values that provides the conditions.",
            Input.Validate.REQUIRED


    );

    private ParametricDistribution[] distr;
    private HashMap<Integer,Integer> valueMap;

    public void initAndValidate(){

        List<ParametricDistribution> distrList = distrListInput.get();

        int[] intVals = convertStringToInteger(integerConditionValuesInput.get());
        if(intVals.length != distrList.size()){
            throw new RuntimeException("The number of Parametric distribution objects, "+distrList.size()+" don't match the number of integer values "+intVals.length+".");
        }
        distr = new ParametricDistribution[intVals.length];
        valueMap = new HashMap<Integer, Integer>();
        for(int i = 0; i < intVals.length; i++){
            valueMap.put(intVals[i],i);
            distr[i] = distrList.get(i);

        }

    }

    public double calcLogP(final Valuable x) throws Exception {
        throw new RuntimeException("Not applicable.");
    }



    /**
     * Calculate log probability of a valuable x for this distribution.
     * If x is multidimensional, the components of x are assumed to be independent,
     * so the sum of log probabilities of all elements of x is returned as the prior.
     */
    public double calcLogP(final Valuable x, int conditionVal) throws Exception {

        return distr[valueMap.get(conditionVal)].calcLogP(x);
    }


    public Double[][] sample(final int size) throws Exception {
        throw new RuntimeException("Not applicable.");

    }

    /*
     * This implemenatation is only suitable for univariate distributions.
     * Must be overwritten for multivariate ones.
     */
    public Double[][] sample(final int size, int conditionVal) throws Exception {

        return distr[valueMap.get(conditionVal)].sample(size);

    }


    public double inverseCumulativeProbability(final double p) throws MathException {
        throw new RuntimeException("Not applicable.");
    }



    /**
     * For this distribution, X, this method returns x such that P(X &lt; x) = p.
     *
     * @param p the cumulative probability.
     * @return x.
     * @throws org.apache.commons.math.MathException if the inverse cumulative probability can not be
     *                       computed due to convergence or other numerical errors.
     */
    //@Override
    public double inverseCumulativeProbability(final double p, int conditionVal) throws MathException {
        return distr[valueMap.get(conditionVal)].inverseCumulativeProbability(p);
    }



    /**
     * Return the probability density for a particular point.
     * NB this does not take offset in account
     *
     * @param x The point at which the density should be computed.
     * @return The pdf at point x.
     */
    //@Override
    public double density(final double x) {
        throw new RuntimeException("Not applicable.");
    }

    /**
     * Return the probability density for a particular point.
     * NB this does not take offset in account
     *
     * @param x The point at which the density should be computed.
     * @return The pdf at point x.
     */
    //@Override
    public double density(final double x, int conditionVal) {
        return distr[valueMap.get(conditionVal)].density(x);
    }


    //@Override
    /** NB logDensity does not take offset in account **/
    public double logDensity(final double x) {
        throw new RuntimeException("Not applicable.");

    }



    //@Override
    /** NB logDensity does not take offset in account **/
    public double logDensity(final double x, int conditionVal) {
        return distr[valueMap.get(conditionVal)].logDensity(x);
    }

    public double cumulativeProbability(final double x) throws MathException {
        throw new RuntimeException("Not applicable.");
    }

    public double cumulativeProbability(final double x, int conditionVal) throws MathException {
        return distr[valueMap.get(conditionVal)].cumulativeProbability(x);
    }


    public double cumulativeProbability(final double x0, final double x1) throws MathException {
        throw new RuntimeException("Not applicable.");
    }

    public double cumulativeProbability(final double x0, final double x1, int conditionVal) throws MathException {
        return distr[valueMap.get(conditionVal)].cumulativeProbability(x0, x1);
    }








    private int[] convertStringToInteger(String str){
        try{
            String[] strArray = str.trim().split("\\s+");
            int[] intVals = new int[strArray.length];
            for(int i = 0; i < intVals.length; i++){
                intVals[i] = Integer.parseInt(strArray[i]);
            }

            return intVals;
        }catch (NumberFormatException e){
            throw new NumberFormatException("Number format exception:\n"+
                "Can only handle integer condition values at the moment."
            );

        }

    }



    public Distribution getDistribution(){
        throw new RuntimeException("Not applicable");
    }
}
