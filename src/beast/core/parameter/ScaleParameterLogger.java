package beast.core.parameter;

import beast.core.BEASTObject;
import beast.core.Function;
import beast.core.Input;
import beast.core.Loggable;
import beast.evolution.tree.Scaler;

import java.io.PrintStream;

/**
 * Created by IntelliJ IDEA.
 * User: cwu080
 * Date: 20/12/13
 * Time: 4:10 PM
 * To change this template use File | Settings | File Templates.
 */
public class ScaleParameterLogger extends BEASTObject implements Loggable{
    //ParameterList
    public Input<Function> parameterInput = new Input<>(
            "parameter",
            "A parameter to be scaled",
            Input.Validate.REQUIRED
    );

    //assignment
    public Input<Scaler> scalerInput = new Input<Scaler>(
            "scaler",
            "A scaler to scale a parameter",
            Input.Validate.REQUIRED
    );

    public Input<Boolean> divideInput = new Input<Boolean>(
            "divide",
            "Whether to divide the valuable by the scaler",
            false
    );

    private Function parameter;
    private Scaler scaler;
    private boolean divide;

    public void initAndValidate(){
        parameter = parameterInput.get();
        scaler = scalerInput.get();
        divide = divideInput.get();
    }

    public void init(PrintStream out){
        //System.out.println(getID());
        int dim = parameter.getDimension();
        for(int i = 0; i < dim; i++){
            out.print(getID()+(i+1)+"\t");
        }
    }


    public void log(int nSample, PrintStream out){
        int dim = parameter.getDimension();
        if(divide){
            for(int i = 0; i < dim; i++){
                out.print(parameter.getArrayValue(i)/scaler.getScaleFactor()+"\t");
            }
        }else{
            for(int i = 0; i < dim; i++){
                out.print(parameter.getArrayValue(i)*scaler.getScaleFactor()+"\t");

            }
        }


    }

    public void close(PrintStream out){}
}
