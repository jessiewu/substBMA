package beast.evolution.operators;

import beast.core.Input;
import beast.core.Description;
import beast.core.parameter.ParameterList;
import beast.util.Randomizer;

/**
 * @author Chieh-Hsi Wu
 */
@Description("Performs a random walk on a network which a single connected component for items in a list.")
public class PLNetworkIntRandomWalkOperator extends NetworkIntRandomWalkOperator{
    public Input<ParameterList> parameterListInput =
                    new Input<ParameterList>("parameters", "the parameter to operate a random walk on within a network.", Input.Validate.REQUIRED);

    public PLNetworkIntRandomWalkOperator(){
        parameterInput.setRule(Input.Validate.OPTIONAL);
    }

    public double proposal(){
        ParameterList paramList = parameterListInput.get(this);
        int iParam = Randomizer.nextInt(paramList.getDimension());
        int currVertex = (int)paramList.getValue(iParam,0) - offset ;
        int nextVertex = neighbours[currVertex][Randomizer.nextInt(neighbours[currVertex].length)];
        paramList.setValue(iParam,0,nextVertex+ offset);
        return logHastingsRatios[currVertex][nextVertex];

    }
}
