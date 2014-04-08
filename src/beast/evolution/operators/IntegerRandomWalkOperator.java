package beast.evolution.operators;

import beast.core.Operator;
import beast.core.Input;
import beast.core.Description;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.util.Randomizer;

/**
 * @author Chieh-Hsi Wu
 */
@Description("A random walk operator that selects a random dimension of the real parameter and perturbs the value a " +
        "random amount within +/- integer windowSize.")
public class IntegerRandomWalkOperator extends Operator {
    public Input<Integer> windowSizeInput =
            new Input<Integer>("windowSize", "the size of the window both up and down when using uniform interval OR standard deviation when using Gaussian", Input.Validate.REQUIRED);

    public Input<RealParameter> parameterInput =
            new Input<RealParameter>("parameter", "the parameter to operate a random walk on.");
    public Input<IntegerParameter> intparameterInput =
            new Input<IntegerParameter>("intparameter", "the parameter to operate a random walk on.", Input.Validate.XOR,parameterInput);
    public final Input<Boolean> input_autoOptimize =
                new Input<Boolean>("autoOptimize", "if true, window size will be adjusted during the MCMC run to improve mixing.", true);

    protected int windowSize;
    double logq;
    double range;
    boolean autoOptimize;
    public void initAndValidate() throws Exception {
        autoOptimize = input_autoOptimize.get();
        windowSize = windowSizeInput.get();
        if(intparameterInput.get() != null){

            range = intparameterInput.get().getUpper() - intparameterInput.get().getLower();
            System.out.println(intparameterInput.get()+" "+ intparameterInput.get().getUpper());
            System.out.println(intparameterInput.get()+" "+ intparameterInput.get().getLower());
            if(windowSize < 1 || windowSize > (intparameterInput.get().getUpper() - intparameterInput.get().getLower())){
                throw new RuntimeException("Window size must be an integer great than "+0+" and less than "+range+".");
            }

        }else{
            range = parameterInput.get().getUpper() - parameterInput.get().getLower();
            System.out.println(parameterInput.get()+" "+ parameterInput.get().getUpper());
            System.out.println(parameterInput.get()+" "+ parameterInput.get().getLower());
            if(windowSize < 1 || windowSize > (parameterInput.get().getUpper() - parameterInput.get().getLower())){
                throw new RuntimeException("Window size must be an integer great than "+0+" and less than "+range+".");
            }
        }

    }

    /**
     * override this for proposals,
     * returns log of hastingRatio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
    @Override

    public double proposal() {
        logq = 0.0;

        if(intparameterInput.get() != null){
            IntegerParameter parameter = intparameterInput.get(this);


            // a random dimension to perturb
            int index = Randomizer.nextInt(parameter.getDimension()); // use getSize(), which = getDimension()

            int newValue = (int)calculateNewValue(parameter.getValue(index),parameter.getUpper(),parameter.getLower());
            parameter.setValue(index,newValue);
        }else{
            RealParameter parameter = parameterInput.get(this);


            // a random dimension to perturb
            int index = Randomizer.nextInt(parameter.getDimension()); // use getSize(), which = getDimension()

            int newValue = (int)calculateNewValue(parameter.getValue(index),parameter.getUpper(),parameter.getLower());
            parameter.setValue(index,(double)newValue);

        }

        return logq;
    }

    protected double calculateNewValue(double oldValue, double upper, double lower) {


        if (upper == lower) return upper;


        double newValue;
        int roll = Randomizer.nextInt(2 * windowSize); // windowSize="1"; roll = {0, 1}
        if (roll >= windowSize) { // roll = 1
            //roll - window is the positive step size
            int step = 1 + (roll - windowSize);
            newValue = oldValue + step;

            if (newValue > upper){
                newValue = 2 * upper - newValue; //reflect down
            }

        } else {  // roll = 0
            newValue = oldValue - 1 - roll;

            if (newValue < lower){
                newValue = 2 * lower - newValue; //reflect up
            }

        }


        //New and seemingly correct (according to the running MCMC with uniform prior)
        //calculation of the hastings ratio --CHW
        int newToOldCount = 0;
        int oldToNewCount = 0;
        if(newValue != oldValue){
            oldToNewCount = oldToNewCount +1;
            newToOldCount = newToOldCount +1;
        }
        double temp = oldValue + windowSize;
        if(temp > upper){
            if((2*upper - temp) <= newValue && newValue != upper){
                oldToNewCount = oldToNewCount+1;
            }
        }

        temp = oldValue - windowSize;
        if(temp < lower){
            if((2*lower - temp) >= newValue && newValue != lower){
                oldToNewCount = oldToNewCount+1;
            }
        }

        temp = newValue + windowSize;
        if( temp > upper){
            if((2*upper - temp) <= oldValue && oldValue != upper){
                newToOldCount = newToOldCount+1;
            }
        }

        temp = newValue - windowSize;
        if( temp < lower){
            if((2*lower - temp) >= oldValue && oldValue != lower){
                newToOldCount = newToOldCount+1;
            }
        }


        logq = Math.log(newToOldCount)- Math.log(oldToNewCount);

        return newValue;
    }


    /**
     * called after every invocation of this operator to see whether
     * a parameter can be optimised for better acceptance hence faster
     * mixing
     *
     * @param logAlpha difference in posterior between previous state & proposed state + hasting ratio
     */
    public void optimize(double logAlpha) {
        if(autoOptimize){
            // must be overridden by operator implementation to have an effect
            double fDelta = calcDelta(logAlpha);

            fDelta += Math.log(windowSize);
            windowSize = (int)Math.exp(fDelta);
            if(windowSize < 1){
                windowSize = 1;
            }else if(windowSize > range){
                windowSize = (int) range;
            }
        }
    }

}
