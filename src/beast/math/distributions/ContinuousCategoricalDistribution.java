package beast.math.distributions;

import beast.core.Valuable;
import beast.core.Input;
import beast.core.Description;
import beast.core.parameter.DPPointer;
import beast.core.parameter.ParameterList;

/**
 * @author Chieh-Hsi Wu
 *
 */
@Description("For testing only. Don't laugh.")
public class ContinuousCategoricalDistribution extends ParametricDistribution{

    public void initAndValidate(){}
    public ParametricDistribution getDistribution(){
        return null;

    }

    public double calcLogP(Valuable val){
        //System.err.println("logDensity");
        if(!(val instanceof DPPointer)){
            throw new RuntimeException("This distriubtion only applies to DPValuable");
        }
        DPPointer pointers = ((DPPointer)val);
        double logL = 0.0;
        //if(pointers.getLastDirty() == -1){

            for(int i = 0; i < pointers.getDimension();i++){
                if(pointers.getParameterValue(i)<0.5){
                    logL+=Math.log(0.2);
                }else{
                    logL+=Math.log(0.8);
                }
                //System.err.print(pointers.getParameterValue(i)+" ");
            }
        //System.err.println();
        //System.err.println("logL: "+logL);
            return logL;
        //}
    



    }
}
