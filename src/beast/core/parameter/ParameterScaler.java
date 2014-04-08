package beast.core.parameter;

import beast.core.Input;
import beast.evolution.tree.Scaler;

/**
 * Created by IntelliJ IDEA.
 * User: cwu080
 * Date: 16/11/13
 * Time: 5:20 PM
 * To change this template use File | Settings | File Templates.
 */
public class ParameterScaler extends Scaler {
    public Input<RealParameter> parameterInput = new Input<RealParameter>(
            "parameter",
            "The parameter to be object as a scaler.",
            Input.Validate.REQUIRED
    );

    private RealParameter parameter;

    public void initAndValidate(){
        parameter = parameterInput.get();

    }

    public double getScaleFactor(){
        return parameter.getValue();

    }




}
