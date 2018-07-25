package beast.core.parameter;

import beast.core.Description;
import beast.core.util.Log;
import beast.math.distributions.ParametricDistribution;
import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ContinuousDistribution;
import org.apache.commons.math.distribution.IntegerDistribution;

/**
 * @author Chieh-Hsi Wu
 */
@Description("A real parameter class that allows values to be set quietly. It's convenience is beyond imagination.")
public class QuietRealParameter extends RealParameter{
    private int idNumber = -1;
    public QuietRealParameter(){

    }

    public void initAndValidate() {
        super.initAndValidate();
    }

    public QuietRealParameter(Double value) {
        super(new Double[]{value});

    }



    public QuietRealParameter(Double[] values) {
        super(values);

    }


    public void setValueQuietly(int dim, Double value){
        values[dim] = value;
        m_bIsDirty[dim] = true;
        m_nLastDirty = dim;
    }

    public void setIDNumber(int idNumber){
        this.idNumber = idNumber;
    }

    public int getIDNumber(){
        return idNumber;
    }

    @Override
    public QuietRealParameter copy() {
        try {
            @SuppressWarnings("unchecked") final
            QuietRealParameter copy = (QuietRealParameter) this.clone();
            copy.values = values.clone();//new Boolean[values.length];
            copy.m_bIsDirty = new boolean[values.length];
            return copy;
        } catch (Exception e) {
            e.printStackTrace();
        }
        return null;
    }

    /**
     * create a {@link QuietRealParameter} sampled from a parametric distribution.
     * If <code>distr.getDistribution()</code> not {@link IntegerDistribution},
     * DO NOT apply offset to values in QuietRealParameter.
     * @see QuietRealParameter#getQuietRealParameter(Double[], double, double, double)
     * @param distr given a parametric distribution
     * @param upper upper bound of {@link QuietRealParameter}
     * @param lower lower bound of {@link QuietRealParameter}
     * @return {@link QuietRealParameter}
     * @throws MathException
     */
    public static QuietRealParameter getSample(ParametricDistribution distr,
                                               double upper, double lower) throws MathException {

        Double[][] sampleVals = distr.sample(1);

        double offset = distr.getOffset();
        // if .getDistribution() not IntegerDistribution, DO NOT apply offset to values in QuietRealParameter
        if ( !(distr.getDistribution() instanceof IntegerDistribution) )
            offset = 0;

        return getQuietRealParameter(sampleVals[0], offset, upper, lower);
    }


    /**
     * create a {@link QuietRealParameter} for sampler
     * {@link QuietRealParameter#getSample(ParametricDistribution, double, double)} and
     * {@link QuietRealParameter#getSamples(ParametricDistribution, int, double, double)}.
     * @param sampleVal the sampled values
     * @param offset the offset defined in a distribution. If NOT {@link IntegerDistribution}, set it to 0.
     *               IntegerDistribution uses offset to adjust the model index
     *               to match the array index of categorical probabilities,
     *               which behaves differently to {@link ContinuousDistribution}
     *               in ParametricDistribution density calculations.
     * @param upper upper bound of {@link QuietRealParameter}
     * @param lower lower bound of {@link QuietRealParameter}
     * @return
     */
    public static QuietRealParameter getQuietRealParameter(Double[] sampleVal, double offset,
                                                            double upper, double lower) {
        for (int i = 0; i < sampleVal.length; i++) {
            // handle the offset issue https://github.com/jessiewu/substBMA/issues/3
            // after CompEvol/beast2#658, ParametricDistribution#sample(int) does not handle offset,
            // it makes a problem for CategoricalDistribution$CategoricalImpl,
            // which uses offseted value as the array index of prob[].
            // make sure this is only used by CategoricalDistribution$CategoricalImpl,
            // otherwise, especially, ContinuousDistribution will count it twice.
            // TODO a better solution for CategoricalDistribution$CategoricalImpl?
            sampleVal[i] += offset;
            if (sampleVal[i] < lower || sampleVal[i] > upper)
                Log.err("Sampled value[" + i + "] = " + sampleVal[i] + " is out of bounds [" +
                        lower + ", " + upper + "] !");
        }

        QuietRealParameter sampleParameter = new QuietRealParameter(sampleVal);
        sampleParameter.setUpper(upper);
        sampleParameter.setLower(lower);

        return sampleParameter;
    }

    /**
     * create an array of {@link QuietRealParameter}s sampled from a parametric distribution given a size.
     * If <code>distr.getDistribution()</code> not {@link IntegerDistribution},
     * DO NOT apply offset to values in QuietRealParameter.
     * @see QuietRealParameter#getQuietRealParameter(Double[], double, double, double)
     * @param distr given a parametric distribution
     * @param sampleSize sample size for {@link ParametricDistribution#sample(int)}
     * @param upper upper bound of all {@link QuietRealParameter}s
     * @param lower lower bound of all {@link QuietRealParameter}s
     * @return an array of {@link QuietRealParameter}s
     * @throws MathException
     */
    public static QuietRealParameter[] getSamples(ParametricDistribution distr, int sampleSize,
                                                  double upper, double lower) throws MathException {
        Double[][] sampleVals = distr.sample(sampleSize);

        QuietRealParameter[] samples = new QuietRealParameter[sampleSize];
        for(int i = 0; i < samples.length;i++){
            double offset = distr.getOffset();
            // if .getDistribution() not IntegerDistribution, DO NOT apply offset to values in QuietRealParameter
            if ( !(distr.getDistribution() instanceof IntegerDistribution) )
                offset = 0;
            samples[i] = getQuietRealParameter(sampleVals[i], offset, upper, lower);
        }
        return samples;
    }

}