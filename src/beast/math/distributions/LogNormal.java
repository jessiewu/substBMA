package beast.math.distributions;

import beast.core.Description;
import beast.core.Input;
import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ContinuousDistribution;
import org.apache.commons.math.distribution.Distribution;
import org.apache.commons.math.distribution.NormalDistributionImpl;

/**
 * @author Chieh-Hsi Wu
 */
@Description("Lognormal distribution which takes the parameter precision instead of standard deviation in log space.")
public class LogNormal extends LogNormalDistributionModel {
    public Input<Boolean> sIsPrecInput = new Input<Boolean>("sIsPrec", "Whether the M parameter is in real space, or in log-transformed space. Default false = log-transformed.", false);

    boolean m_bMeanInRealSpace ;
    boolean sIsPrec;
    LogNormalImpl m_dist = new LogNormalImpl(0, 1);

    public void initAndValidate() {
        sIsPrec = sIsPrecInput.get();
        super.initAndValidate();
    	m_bMeanInRealSpace = hasMeanInRealSpaceInput.get();

        if (MParameterInput.get() != null) {
            if (MParameterInput.get().getLower() == null) {
                MParameterInput.get().setLower(Double.NEGATIVE_INFINITY);
            }
            if (MParameterInput.get().getUpper() == null) {
                MParameterInput.get().setUpper(Double.POSITIVE_INFINITY);
            }
        }

        if (SParameterInput.get() != null) {
            if (SParameterInput.get().getLower() == null) {
                SParameterInput.get().setLower(0.0);
            }
            if (SParameterInput.get().getUpper() == null) {
                SParameterInput.get().setUpper(Double.POSITIVE_INFINITY);
            }
        }
        refresh();
    }

    public boolean requiresRecalculation(){
        if(MParameterInput.get().somethingIsDirty() || SParameterInput.get().somethingIsDirty()){
            refresh();
        }
        return super.requiresRecalculation();

    }

    public void store(){
        m_dist.store();
        super.store();
    }

    public void restore(){
        m_dist.restore();
        super.restore();
    }

    @Override
    public Distribution getDistribution() {
        refresh();
        return m_dist;
    }

    public void setMean(double mean){
        MParameterInput.get().setValue(0,mean);
        refresh();
    }

	/** make sure internal state is up to date **/
	void refresh() {
		double fMean;
		double fSigma;
		if (SParameterInput.get() == null) {
			fSigma = 1;
		} else {
			fSigma = SParameterInput.get().getValue();
		}
		if (MParameterInput.get() == null) {
			fMean = 0;
		} else {
			fMean = MParameterInput.get().getValue();
		}

        if (sIsPrec) {
			fSigma = Math.sqrt(1/fSigma);
		}

		if (m_bMeanInRealSpace) {
			fMean = Math.log(fMean) - (0.5 * fSigma * fSigma);
		}

        //System.err.println(fMean+" "+ fSigma);
		m_dist.setMeanAndStdDev(fMean, fSigma);
	}



	class LogNormalImpl implements ContinuousDistribution {
		double m_fMean;
		double m_fStdDev;
        double storedMean;
		double storedStdDev;
	    NormalDistributionImpl m_normal = new NormalDistributionImpl(0, 1);
	    LogNormalImpl(double fMean, double fStdDev) {
	    	setMeanAndStdDev(fMean, fStdDev);
	    }
		void setMeanAndStdDev(double fMean, double fStdDev) {
			m_fMean = fMean;
			m_fStdDev = fStdDev;
			m_normal.setMean(fMean);
			m_normal.setStandardDeviation(fStdDev);
		}

		@Override
		public double cumulativeProbability(double x) throws MathException {
			return m_normal.cumulativeProbability(Math.log(x));
		}

		@Override
		public double cumulativeProbability(double x0, double x1) throws MathException {
			return cumulativeProbability(x1) - cumulativeProbability(x0);
		}

		@Override
		public double inverseCumulativeProbability(double p) throws MathException {
            return Math.exp(m_normal.inverseCumulativeProbability(p));
		}

		public double density(double fX) {

	        return m_normal.density(Math.log(fX)) / fX;
		}

		public double logDensity(double fX) {
            return m_normal.logDensity(Math.log(fX)) - Math.log(fX);
		}

        public void store(){
            storedMean = m_fMean;
            storedStdDev = m_fStdDev;
        }

        public void restore(){
            m_fMean = storedMean;
            m_fStdDev = storedStdDev;
            m_normal.setMean(storedMean);
			m_normal.setStandardDeviation(storedStdDev);
        }
	} // class LogNormalImpl
}
