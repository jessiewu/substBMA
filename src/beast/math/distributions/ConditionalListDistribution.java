package beast.math.distributions;

import beast.core.Input;
import beast.core.parameter.ParameterList;
import beast.core.parameter.RealParameter;

/**
 * Created by IntelliJ IDEA.
 * User: Jessie Wu
 * Date: 10/09/13
 * Time: 3:25 PM
 * To change this template use File | Settings | File Templates.
 */
public class ConditionalListDistribution extends ConditionalDistribution{
    public Input<ParameterList> conditionParameterListInput = new Input<ParameterList>(
            "conditionXList",
            "The parameter that provides the condition.",
            Input.Validate.REQUIRED
    );

    public Input<ParameterList> parameterListInput = new Input<ParameterList>(
            "xList",
            "A list value at which the density is calculated.",
            Input.Validate.REQUIRED
    );

    public ConditionalListDistribution(){
        parameterInput.setRule(Input.Validate.OPTIONAL);
        conditionParameterInput.setRule(Input.Validate.OPTIONAL);

    }

    private ParameterList conditionParameterList;
    private ParameterList parameterList;

    public void initAndValidate(){
        conditionParameterList = conditionParameterListInput.get();
        parameterList = parameterListInput.get();
        conditionalParametricDistr = conditionalParametricDistrInput.get();
        calculateLogP();

    }

    public double calculateLogP(){
        try{
            logP = 0.0;
            double dim = parameterList.getDimension();
            for(int i = 0; i < dim; i++){
                logP += conditionalParametricDistr.calcLogP(
                        parameterList.getParameter(i),
                        (int)(double)conditionParameterList.getParameter(i).getValue()
                );
            }

            //System.out.println("conditional: "+conditionParameter);
            return logP;
        }catch(Exception e){
            throw new RuntimeException(e);
        }
    }

}
