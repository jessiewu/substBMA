package beast.math.distributions;

import beast.core.Description;
import beast.core.Input;
import org.apache.commons.math.distribution.ContinuousDistribution;
import org.apache.commons.math.distribution.Distribution;
import org.apache.commons.math.MathException;

/**
 * @author Chieh-Hsi Wu
 */
@Description("Uniform distribution for integers.")
public class IntegerUniformDistribution extends ParametricDistribution{
    public Input<Double> lowerInput = new Input<Double>("lower","lower bound on the interval, defaul 0", Input.Validate.REQUIRED);
	public Input<Double> upperInput = new Input<Double>("upper","lower bound on the interval, defaul 1", Input.Validate.REQUIRED);


    IntegerUniformImpl dist = new IntegerUniformImpl();

	double lower, upper;

	@Override
	public void initAndValidate() throws Exception {
		refresh();
	}

    public void refresh(){
        lower = lowerInput.get();
		upper = upperInput.get();
		if (lower >= upper) {
			throw new RuntimeException("Upper value should be higher than lower value");
		}
        dist.setBounds(lower,upper);

    }

    class IntegerUniformImpl implements ContinuousDistribution {
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
			return (x1-x0+1)/(upper-lower+1);
		}


		public double inverseCumulativeProbability(double p) throws MathException {
            if(p< 0.0|| p> 1.0){
                throw new RuntimeException("Where have you ever seen a probility greater than 1 or less than 0? Definitely not in my world.");
            }
			return Math.floor((upper-lower+1)*p)+lower;
		}



		public double density(double x) {


			return 1.0/(upper-lower+1);
		}

		public double logDensity(double x) {
            if(x > upper || x < lower)
                return Double.NEGATIVE_INFINITY;
			return -Math.log(upper-lower+1);
		}
    }

	@Override
	public Distribution getDistribution() {
		return dist;
	}
}
