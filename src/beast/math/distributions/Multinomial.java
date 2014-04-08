package beast.math.distributions;

import beast.core.Input;
import beast.core.Valuable;
import beast.core.Description;
import beast.core.parameter.RealParameter;
import org.apache.commons.math.distribution.Distribution;

/**
 * @author Chieh-Hsi Wu
 **/
@Description("This class implements the multinomial distribution.")
public class Multinomial extends ParametricDistribution{
    public Input<RealParameter> m_probs = new Input<RealParameter>("probs","coefficients of the Multinomial distribution and should add up to one", Input.Validate.REQUIRED);
    private Double[] probs;

    public void initAndValidate(){
        refresh();


    }

    public void refresh(){
        probs = m_probs.get().getValues();
        //check that probs sum up to 1
        double sumP = 0;
        for(Double prob: probs){
            sumP +=prob;
        }
        if((sumP - 1.0) > 1e-10){
            throw new RuntimeException("Probabilities don't sum to one");
        }

    }


    public boolean requiresRecalculation(){
        refresh();
        return super.requiresRecalculation();
    }

    public void restore(){
        refresh();
    }


    @Override
	public Distribution getDistribution() {
		return null;
	}

    @Override
	public double calcLogP(Valuable xVec) throws Exception {
        

        int xVecDim  = xVec.getDimension();
        double xVecSum = 0;
        double sumLogXElFac =0;
        double sumXLogP = 0;

        for(int i = 0; i < xVecDim; i++){
            double xEl = xVec.getArrayValue(i);
            xVecSum += xEl;
            sumLogXElFac += org.apache.commons.math.special.Gamma.logGamma(xEl+1);
            sumXLogP += xEl*Math.log(probs[i]);

        }
        double logNFac = org.apache.commons.math.special.Gamma.logGamma(xVecSum+1);

        double fLogP = logNFac - sumLogXElFac + sumXLogP;
		return fLogP;
	}

    

}
