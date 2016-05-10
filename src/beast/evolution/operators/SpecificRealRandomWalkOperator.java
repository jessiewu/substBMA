package beast.evolution.operators;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.util.Randomizer;

import java.security.PublicKey;

/**
 * @author Chieh-Hsi Wu
 */
@Description("Only apply random walk on the value of a specific index.")
public class SpecificRealRandomWalkOperator extends RealRandomWalkOperator {
    public Input<Integer> indexInput = new Input<Integer>("index", "The index of value in the real parameter to be operated on.", Input.Validate.REQUIRED);
    private int index;
    public void initAndValidate(){
        super.initAndValidate();
        index = indexInput.get();

    }
    public double proposal() {

        RealParameter param = parameterInput.get(this);


        double value = param.getValue(index);
        double newValue = value;
        if (useGaussian) {
            newValue += Randomizer.nextGaussian() * windowSize;
        } else {
            newValue += Randomizer.nextDouble() * 2 * windowSize - windowSize;
        }

        if (newValue < param.getLower() || newValue > param.getUpper()) {
            return Double.NEGATIVE_INFINITY;
        }
        if (newValue == value) {
            // this saves calculating the posterior
            return Double.NEGATIVE_INFINITY;
        }

        param.setValue(index, newValue);

        return 0.0;
    }
}
