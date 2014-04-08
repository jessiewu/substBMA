package beast.math.distributions;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import org.apache.commons.math.distribution.ContinuousDistribution;
import org.apache.commons.math.distribution.GammaDistributionImpl;

/**
 * @author Chieh-Hsi Wu
 */
public class GammaDistribution extends ParametricDistribution {
    public Input<RealParameter> shapeInput = new Input<RealParameter>(
            "shape",
            "shape parameter of the gamma distribution",
            Input.Validate.REQUIRED);
    public Input<RealParameter> ratesInput = new Input<RealParameter>(
            "rate",
            "rate parameter of the gamma distribution",
            Input.Validate.REQUIRED);

    static org.apache.commons.math.distribution.GammaDistribution m_dist = new GammaDistributionImpl(1, 1);

    @Override
    public void initAndValidate() {
        refresh();
    }

    /**
     * make sure internal state is up to date *
     */
    void refresh() {
        double shape;
        double rate;
        shape = shapeInput.get().getValue();
        rate = ratesInput.get().getValue();

        m_dist.setAlpha(shape);
        m_dist.setBeta(1.0/rate);
    }

    @Override
    public ContinuousDistribution getDistribution() {
        refresh();
        return m_dist;
    }

    @Override
    public double getMean() {
    	return m_offset.get() + m_dist.getAlpha() / m_dist.getBeta();
    }
}
