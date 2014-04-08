package beast.evolution.operators;

import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.RealParameter;
import beast.math.distributions.ParametricDistribution;
import beast.util.Randomizer;

import java.util.Random;

/**
 * Created by IntelliJ IDEA.
 * User: Jessie Wu
 * Date: 26/01/14
 * Time: 2:05 PM
 * To change this template use File | Settings | File Templates.
 */
public class DistributionSampler extends Operator {
    public Input<RealParameter> parameterInput = new Input<RealParameter>(
            "parameter",
            "the parameter to operate a random walk on.",
            Input.Validate.REQUIRED
    );

    public Input<ParametricDistribution> distrInput = new Input<ParametricDistribution>(
            "distr",
            "The proposal distribution of the parameter.",
            Input.Validate.REQUIRED
    );


    private RealParameter parameter;
    private ParametricDistribution distr;
    public void initAndValidate(){
        parameter =  parameterInput.get();
        distr = distrInput.get();

    }

    public double proposal(){
        try{
            int dim = Randomizer.nextInt(parameter.getDimension());
            double oldValue = parameter.getValue(dim);
            double newValue = distr.inverseCumulativeProbability(Randomizer.nextDouble());
            parameter.setValue(dim,newValue);
            return Math.log(distr.density(oldValue)/distr.density(newValue));
        }catch(Exception e){
            throw new RuntimeException(e);
        }


    }
}
