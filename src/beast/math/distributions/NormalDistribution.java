package beast.math.distributions;

import beast.core.Input;
import beast.core.parameter.RealParameter;

/**
 * Created by IntelliJ IDEA.
 * User: cwu080
 * Date: 9/07/13
 * Time: 4:21 PM
 * To change this template use File | Settings | File Templates.
 */
public class NormalDistribution extends Normal {
    public Input<RealParameter> precisionInput = new Input<RealParameter>("precision", "variance of the normal distribution, defaults to 1");

    void refresh() {
        double fMean;
        double fSigma;
        if (meanInput.get() == null) {
            fMean = 0;
        } else {
            fMean = meanInput.get().getValue();
        }
        if (sigmaInput.get() != null) {
            fSigma = sigmaInput.get().getValue();
        }else if(precisionInput.get() != null){
            fSigma = 1.0/precisionInput.get().getValue();

        } else {
            fSigma = 1;
        }

        dist.setMean(fMean);
        dist.setStandardDeviation(fSigma);
    }
}
