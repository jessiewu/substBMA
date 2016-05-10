package beast.evolution.operators;

import beast.core.Input;
import beast.core.Description;
import beast.core.parameter.RealParameter;
import beast.util.Randomizer;

/**
 * @author Chieh-Hsi Wu
 */
@Description("This random walk operator performs random-walk with different window sizes on different indices of a parameter of dimension > 1.")
public class ExtendedRealRandomWalkOperator extends RealRandomWalkOperator{
    public Input<Boolean> autoOptimizeInput = new Input<Boolean>("autoOptimize", "Boolean to indicate whether to auto-optimize the operator.", true);

    public Input<String> windowSizesInput =
            new Input<String>("windowSizes", "the size of the window both up and down when using uniform interval OR standard deviation when using Gaussian", Input.Validate.REQUIRED);
    /*public Input<RealParameter> cpInput =
                new Input<RealParameter>("parameter", "the parameter to operate a random walk on.", Input.Validate.REQUIRED);
    */
    protected double[] windowSizes;
    protected int lastChangedValueIndex;

    public ExtendedRealRandomWalkOperator(){
        windowSizeInput.setRule(Input.Validate.OPTIONAL);
        parameterInput.setRule(Input.Validate.OPTIONAL);
    }

    public void initAndValidate() {
        windowSizes = new double[parameterInput.get().getDimension()];
        String[] windowSizesStr =  windowSizesInput.get().split("\\s+");
        setupWindowSizes(windowSizesStr);
    }

    protected void setupWindowSizes(String[] windowSizesStr){
        if(windowSizesStr.length == 1){
            double windowSize = Double.parseDouble(windowSizesStr[0]);
            for(int i = 0; i < windowSizes.length; i++){
                windowSizes[i] = windowSize;
            }
        }else if(windowSizesStr.length == windowSizes.length){
            for(int i = 0; i < windowSizes.length; i++){
                windowSizes[i] = Double.parseDouble(windowSizesStr[i]);
            }
        }else{
            throw new RuntimeException("The number of window sizes should be either 1 or the same as the dimension of the parameter");
        }

    }

    /**
     * override this for proposals,
     * returns log of hastingRatio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
    @Override
    public double proposal() {
        RealParameter param = parameterInput.get(this);
        int i = Randomizer.nextInt(param.getDimension());
        double value = param.getValue(i);
        double newValue = getProposedVal(
                value, //original value
                i,     //index of the value
                useGaussian, //whether to use gaussian moves
                windowSizes);   //window sizes

        if (newValue < param.getLower() || newValue > param.getUpper()) {
        	return Double.NEGATIVE_INFINITY;
        }
        if (newValue == value) {
        	// this saves calculating the posterior
        	return Double.NEGATIVE_INFINITY;
        }

        param.setValue(i, newValue);
        lastChangedValueIndex = i;


        return 0.0;
    }

    public double getProposedVal(
            double value,
            int iValue,
            boolean useGaussian,
            double[] windowSizes){

        double newValue = value;
        if (useGaussian) {
        	newValue += Randomizer.nextGaussian()* windowSizes[iValue];
        } else {
        	newValue += Randomizer.nextDouble() * 2 * windowSizes[iValue] - windowSizes[iValue];
        }


        return newValue;

    }

    public void optimize(double logAlpha) {
        if(autoOptimizeInput.get()){
        // must be overridden by operator implementation to have an effect
        double fDelta = calcDelta(logAlpha);

        fDelta += Math.log(windowSizes[lastChangedValueIndex]);

        //double fScaleFactor = m_pScaleFactor.get();
        //fDelta += Math.log(1.0 / windowSize - 1.0);
        //windowSize = 1.0 / (Math.exp(fDelta) + 1.0);
        windowSizes[lastChangedValueIndex] = Math.exp(fDelta);
        //System.out.println(windowSize);
        }
    }

}
