package beast.evolution.operators;

import beast.core.parameter.IntegerParameter;
import beast.core.parameter.ParameterList;
import beast.core.Input;
import beast.core.Operator;
import beast.core.Description;
import beast.util.Randomizer;

/**
 * @author Chieh-Hsi Wu
 */
@Description("A generic operator for use with sum-constrained (possibly weighted) vector parameters in a list.")
public class PLDeltaExchangeOperator extends Operator {
    public Input<ParameterList> parameterListInput = new Input<ParameterList>("parameterList", "if specified, this parameter is scaled", Input.Validate.REQUIRED);

    public Input<Double> input_delta = new Input<Double>("delta", "Magnitude of change for two randomly picked values.", 1.0);
    public Input<Boolean> input_autoOptimize =
            new Input<Boolean>("autoOptimize", "if true, window size will be adjusted during the MCMC run to improve mixing.", true);
    public Input<Boolean> input_isIntegerOperator = new Input<Boolean>("integer", "if true, changes are all integers.", false);
    public Input<IntegerParameter> input_parameterWeights = new Input<IntegerParameter>("parameter", "weights on a vector parameter",
            Input.Validate.OPTIONAL);

    private boolean autoOptimize;
    private double delta;
    private boolean isIntegerOperator;
    private int[] parameterWeights;

   
    public void initAndValidate() {

        autoOptimize = input_autoOptimize.get();
        delta = input_delta.get();
        isIntegerOperator = input_isIntegerOperator.get();

        parameterWeights = new int[parameterListInput.get().getParameterDimension()];

        if(input_parameterWeights.get() != null){
            for (int i = 0; i < parameterWeights.length; i++) {
                parameterWeights[i] = input_parameterWeights.get().getValue(i);
            }
        }else{
            for (int i = 0; i < parameterWeights.length; i++) {
                parameterWeights[i] = 1;
            }
        }

        if (isIntegerOperator && delta != Math.round(delta)) {
            throw new IllegalArgumentException("Can't be an integer operator if delta is not integer");
        }
    }

    public final double proposal() {
        //System.err.println("PLDE");
        double logq = 0.0;
        // get two dimensions
        ParameterList parameterList = parameterListInput.get(this);
        int iParam = Randomizer.nextInt(parameterList.getDimension());
        int dimVal = parameterList.getParameterDimension();
        final int iValue1 = Randomizer.nextInt(dimVal);
        int iValue2 = iValue1;
        while (iValue1 == iValue2) {
            iValue2 = Randomizer.nextInt(dimVal);
        }

        double scalar1 = parameterList.getValue(iParam,iValue1);
        double scalar2 = parameterList.getValue(iParam,iValue2);

        if (isIntegerOperator) {
            int d = Randomizer.nextInt((int) Math.round(delta)) + 1;

            if (parameterWeights[iValue1] != parameterWeights[iValue2]) throw new RuntimeException();
            scalar1 = Math.round(scalar1 - d);
            scalar2 = Math.round(scalar2 + d);
        } else {

            // exchange a random delta
            final double d = Randomizer.nextDouble() * delta;
            scalar1 -= d;
            if (parameterWeights[iValue1] != parameterWeights[iValue2]) {
                scalar2 += d * (double) parameterWeights[iValue1] / (double) parameterWeights[iValue2];
            } else {
                scalar2 += d;
            }

        }


        if (scalar1 < parameterList.getLower() ||
                scalar1 > parameterList.getUpper() ||
                scalar2 < parameterList.getLower() ||
                scalar2 > parameterList.getUpper()) {
            logq = Double.NEGATIVE_INFINITY;
        }else{

            parameterList.setValue(iParam,iValue1, scalar1);
            parameterList.setValue(iParam,iValue2, scalar2);
        }
        //System.err.println("apply deltaEx");
        // symmetrical move so return a zero hasting ratio
        return logq;
    }

    /**
     * called after every invocation of this operator to see whether
     * a parameter can be optimised for better acceptance hence faster
     * mixing
     *
     * @param logAlpha difference in posterior between previous state & proposed state + hasting ratio
     */
    public void optimize(double logAlpha) {
        // must be overridden by operator implementation to have an effect
        if(autoOptimize){
            double fDelta = calcDelta(logAlpha);
            fDelta += Math.log(delta);
            delta = Math.exp(fDelta);
        }

    }
}
