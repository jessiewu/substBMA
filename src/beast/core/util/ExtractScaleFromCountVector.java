package beast.core.util;

import beast.core.Input;
import beast.core.Valuable;
import beast.core.parameter.RealParameter;

/**
 * Created by IntelliJ IDEA.
 * User: cwu080
 * Date: 24/04/13
 * Time: 5:16 PM
 * To change this template use File | Settings | File Templates.
 */
public class ExtractScaleFromCountVector implements Valuable {
    public Input<RealParameter> parameterInput = new Input<RealParameter>("parameter", "The vector as a parameter.", Input.Validate.REQUIRED);
    private RealParameter parameter;
    public int getDimension(){
        return 1;
    }

    public double getArrayValue(){
        double sum = 0.0;
        for(int i = 0; i < parameter.getDimension(); i++){
            sum += parameter.getValue(i);
        }
        return sum;

    }

    public double getArrayValue(int i){
        double sum = 0.0;
        for(int j = 0; j < parameter.getDimension(); j++){
            sum += parameter.getValue(j);
        }
        return sum;

    }


}
