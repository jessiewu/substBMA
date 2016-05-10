package beast.math.distributions;

import beast.core.Input;
import beast.core.Description;
import org.apache.commons.math.distribution.Distribution;
import org.apache.commons.math.distribution.ContinuousDistribution;
import org.apache.commons.math.MathException;


/**
 * @author Chieh-Hsi Wu
 */
@Description("Uniform distribution.  The one in BEAST 2 is not completed.")
public class UniformDistribution extends ParametricDistribution{
	public Input<Double> lowerInput = new Input<Double>("lower","lower bound on the interval, defaul 0", 0.0);
	public Input<Double> upperInput = new Input<Double>("upper","lower bound on the interval, defaul 1", 1.0);

    UniformImpl dist = new UniformImpl();

	double m_fLower, m_fUpper;

	@Override
	public void initAndValidate() {
		refresh();
	}

    public void refresh(){
        m_fLower = lowerInput.get();
		m_fUpper = upperInput.get();
		if (m_fLower >= m_fUpper) {
			throw new RuntimeException("Upper value should be higher than lower value");
		}
        dist.setBounds(m_fLower,m_fUpper);

    }

    class UniformImpl implements ContinuousDistribution {
        private double lower;
        private double upper;

        public void setBounds(double lower, double upper){
            this.lower = lower;
            this.upper = upper;
        }

		@Override
		public double cumulativeProbability(double x) throws MathException {
            if(x< lower|| x> upper){
                throw new RuntimeException("Value x ("+x+") out of bounds ("+lower+","+upper+").");
            }
			return (x-lower)/(upper-lower);
		}

		@Override
		public double cumulativeProbability(double x0, double x1) throws MathException {
            if(x0< lower|| x0> upper){
                throw new RuntimeException("Value x ("+x0+") out of bounds ("+lower+","+upper+").");
            }
            if(x1< lower|| x1> upper){
                throw new RuntimeException("Value x ("+x1+") out of bounds ("+lower+","+upper+").");
            }
			return (x1-x0)/(upper-lower);
		}

		@Override
		public double inverseCumulativeProbability(double p) throws MathException {
            if(p< 0.0|| p> 1.0){
                throw new RuntimeException("Where have you ever seen a probility greater than 1 or less than 0? Definitely not in my world.");
            }
			return (upper-lower)*p+lower;
		}


		public double density(double x) {
            if(x > upper || x < lower){
                return 0;
            }
			return 1.0/(upper-lower);
		}

		public double logDensity(double x) {
            if(x > upper || x < lower){
                return Double.NEGATIVE_INFINITY;
            }
			return -Math.log(upper-lower);
		}
	} // class OneOnXImpl


	@Override
	public Distribution getDistribution() {

		return dist;
	}

   


}
