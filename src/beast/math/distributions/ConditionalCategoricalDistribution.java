package beast.math.distributions;

import beast.core.Input;
import beast.core.Plugin;
import beast.core.parameter.RealParameter;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Chieh-Hsi Wu
 */
public class ConditionalCategoricalDistribution extends Plugin {
    public Input<List<RealParameter>> conditionalDensityInput =
            new Input<List<RealParameter>>(
                    "condProb",
                    "Probabilities of each integer value",
                    new ArrayList<RealParameter>(),
                    Input.Validate.REQUIRED
            );

    public Input<Integer> offsetInput =
            new Input<Integer>(
                    "offset",
                    "Offset of the integer value",
                    0
            );

    private double[][] conditionalDensity;
    private double[][] logCondDensity;
    private int offset;

    public void initAndValidate(){
        List<RealParameter> conditionalDistrList = conditionalDensityInput.get();
        conditionalDensity = new double[conditionalDistrList.size()][];
        logCondDensity = new double[conditionalDensity.length][];
        for(int i = 0; i < conditionalDensity.length;i++){
            conditionalDensity[i] = new double[conditionalDistrList.get(i).getDimension()];
            logCondDensity[i] = new double[conditionalDensity.length];
            if(conditionalDensity.length != conditionalDensity[i].length){
                throw new RuntimeException("State space dont' match!");
            }
            for(int j = 0; j < conditionalDensity[i].length; j++){
                conditionalDensity[i][j] = conditionalDistrList.get(i).getValue(j);
                logCondDensity[i][j] = Math.log(conditionalDensity[i][j]);
            }
        }
        offset = offsetInput.get();
    }

    public double conditionalDensity(int given, int x){
        return conditionalDensity[given - offset][x - offset];
    }

    public double[] conditionalDensities(int given){
        return conditionalDensity[given - offset].clone();
    }

    public double logConditionalDensity(int given, int x){
        return logCondDensity[given-offset][x-offset];
    }

    public int getOffset(){
        return offset;
    }




}
