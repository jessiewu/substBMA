package beast.evolution.operators;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.ParameterList;
import beast.util.Randomizer;

/**
 * @author Chieh-Hsi Wu
 */
@Description("A random walk operator that selects a random dimension of a random RealParameter from ParameterList " +
        "and perturbs the value a random amount within +/- integer windowSize.")
public class PLIntegerRandomWalkOperator extends IntegerRandomWalkOperator{
    public Input<ParameterList> parameterListInput =
            new Input<ParameterList>(
                    "parameters",
                    "the parameter to operate a random walk on.",
                    Input.Validate.REQUIRED
            );
    public PLIntegerRandomWalkOperator(){
        parameterInput.setRule(Input.Validate.OPTIONAL);
        intparameterInput.setRule(Input.Validate.OPTIONAL);
    }
    public void initAndValidate() throws Exception {
        //super.initAndValidate();
        windowSize = windowSizeInput.get();
        range = parameterListInput.get().getUpper() - parameterListInput.get().getLower();
        if(windowSize < 1 || windowSize > range){
            throw new RuntimeException("Window size must be an integer great than "+0+" and less than "+range+".");
        }
    }

    @Override
    public double proposal() {
        logq =0.0;

        ParameterList paramList = parameterListInput.get(this);


        // a random parameter to perturb
        int iParam = Randomizer.nextInt(paramList.getDimension());

        //a random dimesion to perturb
        int iValue = Randomizer.nextInt(paramList.getParameterDimension());
        double newValue = calculateNewValue(
                paramList.getValue(iParam,iValue),
                paramList.getUpper(),
                paramList.getLower()
        );
        paramList.setValue(iParam,iValue,newValue);

        return logq;
    }



}
