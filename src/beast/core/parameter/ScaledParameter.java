package beast.core.parameter;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Loggable;
import beast.core.Valuable;

import java.io.PrintStream;

/**
 * Created by IntelliJ IDEA.
 * User: Jessie Wu
 * Date: 29/01/14
 * Time: 3:31 PM
 * To change this template use File | Settings | File Templates.
 */
public class ScaledParameter extends CalculationNode implements Valuable , Loggable {
    public Input<RealParameter> parameterInput = new Input<RealParameter>(
            "parameter",
            "The parameter that is to be scaled.",
            Input.Validate.REQUIRED
    );

    public Input<RealParameter> scalerInput = new Input<RealParameter>(
            "scaler",
            "The parameter used to scale another parameter.",
            Input.Validate.REQUIRED
    );

    private RealParameter parameter;
    private RealParameter scaler;
    public void initAndValidate(){
        parameter = parameterInput.get();
        scaler = scalerInput.get();

    }

    public int getDimension(){
        return parameter.getDimension();

    }

    public double getArrayValue(){
        return parameter.getArrayValue()*scaler.getArrayValue();
    }

    public double getArrayValue(int i){
        return parameter.getArrayValue(i)*scaler.getArrayValue();
    }

    public void init(PrintStream out){
        //System.out.println(getID());
        int dim = parameter.getDimension();
        for(int i = 0; i < dim; i++)
        out.print(getID()+(i+1)+"\t");
    }

    public void log(int nSample, PrintStream out){

        int dim = parameter.getDimension();
        for(int i = 0; i < dim;i++){
            out.print(getArrayValue(i)+"\t");
        }

    }

    public void close(PrintStream out){}

}
