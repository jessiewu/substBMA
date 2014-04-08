package beast.evolution.operators;

import beast.core.Input;
import beast.core.Description;
import beast.core.parameter.ParameterList;
import beast.util.Randomizer;

/**
 * @author Chieh-Hsi Wu
 */
@Description("This random walk operator performs is adjusted so that it can perform" +
        " extended random-walk on a list of parameters.")
public class PLExtendedRealRandomWalkOperator extends ExtendedRealRandomWalkOperator{
     public Input<ParameterList> parameterListInput =
                new Input<ParameterList>("parameters", "the parameter to operate a random walk on.", Input.Validate.REQUIRED);
     
    public PLExtendedRealRandomWalkOperator(){
        super();
        parameterInput.setRule(Input.Validate.OPTIONAL);
    }

    private ParameterList parameterList;

    int lastChangedParameterIndex;

    public void initAndValidate() {
        parameterList = parameterListInput.get();
        //windowSize = windowSizeInput.get();
        String[] windowSizesStr =  windowSizesInput.get().split("\\s+");
        windowSizes = new double[parameterList.getParameterDimension()];
        setupWindowSizes(windowSizesStr);
    }

  

    /**
     * override this for proposals,
     * returns log of hastingRatio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
    @Override
    public double proposal() {
        //System.err.println("PLERRWO");
        parameterList = parameterListInput.get(this);
        int iParam = Randomizer.nextInt(parameterList.getDimension());
        int iValue = Randomizer.nextInt(parameterList.getParameterDimension());
        double value = parameterList.getValue(iParam,iValue);
        /*double newValue = getProposedVal(
                value, //original value
                iValue,     //index of the value
                m_bUseGaussian, //whether to use gaussian moves
                windowSizes);   //window sizes
                */

        double newValue = value;
        double step= Randomizer.nextDouble() * 2 * windowSizes[iValue] - windowSizes[iValue];
        newValue += step;

        /*for(int i = 0;i < parameterList.getDimension();i++){
            System.err.print(parameterList.getParameter(i).getValue()+" ");
        }*/
        /*System.err.println("old: "+value+", new: "+newValue+", step: "+step);
        System.err.println("iParam: "+iParam+
                ", lower: "+parameterList.getLower(iParam)+
                ", upper: "+parameterList.getUpper(iParam));*/
        if (newValue < parameterList.getParameterLower(iParam) || newValue > parameterList.getParameterUpper(iParam)) {
            //System.err.println("flag1");
        	return Double.NEGATIVE_INFINITY;
        }
        if (newValue == value) {
        	// this saves calculating the posterior
        	return Double.NEGATIVE_INFINITY;
        }

        parameterList.setValue(iParam, iValue, newValue);


        lastChangedParameterIndex = iParam;
        lastChangedValueIndex = iValue;


        return 0.0;
    }

}
