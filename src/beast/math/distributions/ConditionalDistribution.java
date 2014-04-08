package beast.math.distributions;

import beast.core.Distribution;
import beast.core.Input;
import beast.core.Plugin;
import beast.core.State;
import beast.core.parameter.RealParameter;

import java.util.List;
import java.util.Random;

/**
 * @author Chieh-Hsi Wu
 */
public class ConditionalDistribution extends Distribution {

    public Input<RealParameter> conditionParameterInput = new Input<RealParameter>(
            "conditionParameter",
            "The parameter that provides the condition.",
            Input.Validate.REQUIRED
    );

    public Input<RealParameter> parameterInput = new Input<RealParameter>(
            "parameter",
            "value at which the density is calculated.",
            Input.Validate.REQUIRED
    );

    public Input<ConditionalParametricDistribution> conditionalParametricDistrInput = new Input<ConditionalParametricDistribution>(
            "distr",
            "Conditional distribution used to calculate prior",
            Input.Validate.REQUIRED

    );

    private RealParameter conditionParameter;
    private RealParameter parameter;
    protected ConditionalParametricDistribution conditionalParametricDistr;

    public void initAndValidate(){
        conditionParameter = conditionParameterInput.get();
        parameter = parameterInput.get();
        conditionalParametricDistr = conditionalParametricDistrInput.get();
        calculateLogP();
    }

    public double calculateLogP(){
        try{
            //System.out.println("conditional: "+conditionParameter);
            logP =  conditionalParametricDistr.calcLogP(parameter,(int)(double)conditionParameter.getValue());
            return logP;
        }catch(Exception e){
            throw new RuntimeException(e);
        }
    }

    /** return name of the parameter this prior is applied to **/
    public String getParameterName() {
    	if (parameter instanceof Plugin) {
    		return ((Plugin) parameter).getID();
    	}
    	return parameter + "";
    }

    @Override
    public void sample(State state, Random random) {
    }

    @Override
    public List<String> getArguments() {
        return null;
    }

    @Override
    public List<String> getConditions() {
        return null;
    }
}

