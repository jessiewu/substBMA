package beast.core.parameter;

import beast.core.Description;
import beast.core.Input;

/**
 * @author Chieh-Hsi Wu
 */
@Description("Wrap around a parameter but only look at one value.")
public class RealParameterWrapper extends QuietRealParameter{
    public Input<RealParameter> parameterInput = new Input<RealParameter>(
            "parameter",
            "Parameter to be wrapped around",
            Input.Validate.REQUIRED
    );
    public Input<Integer> indexInput = new Input<Integer>(
            "index",
            "The index the wrapper should refer to", 
            Input.Validate.REQUIRED
    );
    private RealParameter parameter;
    private int index;

    public void initAndValidate(){
        parameter = parameterInput.get();
        index = indexInput.get();
    }
    public RealParameterWrapper(RealParameter parameter, int index){
        this.parameter = parameter;
        this.index = index;

    }

    @Override
    public Double getValue(){
        return parameter.getValue(index);

    }

    public boolean somethingIsDirty(){
        return parameter.isDirty(index);
    }

    @Override
    public void setValue(Double val){
        parameter.setValue(index,val);
    }

    @Override
    public Double getValue(int index){
        throw new RuntimeException("This is a one dimensional parameter");
    }

    @Override
    public void setValue(int index, Double value){
        throw new RuntimeException("This is a one dimensional parameter");
    }


}
