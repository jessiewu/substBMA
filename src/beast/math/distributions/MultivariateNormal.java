package beast.math.distributions;


import beast.core.Function;
import beast.core.Input;
import beast.core.Description;
import beast.core.parameter.RealParameter;
import beast.math.matrixAlgebra1.*;
import beast.math.matrixAlgebra1.SymmetricMatrix;
import beast.util.Randomizer;
import org.apache.commons.math.distribution.ContinuousDistribution;

/**
 * @author Marc Suchard
 *
 */

@Description("This class implements the multivariate nromal distribution.")
public class MultivariateNormal extends ParametricDistribution implements MultivariateDistribution{
    public Input<RealParameter> meanVec = new Input<RealParameter>("mean","Mean vector of the multivariate normal distribution", Input.Validate.REQUIRED);
    public Input<RealParameter> precisionMatrixInput = new Input<RealParameter>("precision", "Precision matrix of the multivariate normal distribution", Input.Validate.REQUIRED);
    private double[] mean;
    private double[][] precision;
    private double[][] variance = null;
    private double[][] cholesky = null;
    private double logDet;

    private RealParameter precisionMatrix;

    @Override
	public void initAndValidate() {
        System.err.println("Using multivariate normal");
        precisionMatrix = precisionMatrixInput.get();

        refresh();
	}

    public boolean requiresRecalculation(){
        refresh();
        return true;
    }

    public void refresh(){
        RealParameter meanVec = this.meanVec.get();
        mean = new double[meanVec.getDimension()];
        for(int i = 0; i < mean.length;i++){
            mean[i] = meanVec.getValue(i);
            
        }

        int dim = (int)Math.sqrt(precisionMatrix.getDimension());
        int k = 0;
        precision = new double[dim][dim];
        for(int i = 0; i < precision.length;i++){
            for(int j = 0; j < precision[i].length;j++){
                precision[i][j] = precisionMatrix.getValue(k++);
            }

        }

        logDet = Math.log(calculatePrecisionMatrixDeterminate(precision));
        variance = new SymmetricMatrix(precision).inverse().toComponents();
        cholesky = getCholeskyDecomposition(variance);

        storedMean = new double[mean.length];
        storedPrecision = new double[precision.length][precision.length];
        storedVariance = new double[variance.length][variance.length];
        storedCholesky = new double[cholesky.length][cholesky.length];




    }

    private double[] storedMean;
    private double[][] storedPrecision;
    private double storedLogDet;
    private double[][] storedVariance;
    private double[][] storedCholesky;

    public void store(){
        System.arraycopy(mean,0,storedMean,0,mean.length);

        for(int i = 0; i < storedPrecision.length; i++){
            System.arraycopy(precision[i],0,storedPrecision[i],0,precision[i].length);
        }

        storedLogDet = logDet;

        for(int i = 0; i < storedVariance.length; i++){
            System.arraycopy(variance[i],0,storedVariance[i],0,variance[i].length);
        }

        for(int i = 0; i < storedCholesky.length; i++){
            System.arraycopy(cholesky[i],0,storedCholesky[i],0,cholesky[i].length);
        }

    }

    public void restore(){
        double[] tmp1 = mean;
        mean = storedMean;
        storedMean = tmp1;

        double[][] tmp2 = precision;
        precision = storedPrecision;
        storedPrecision = tmp2;

        double[][] tmp3 = storedVariance;
        variance = storedVariance;
        storedVariance = tmp3;

        double[][] tmp4 = storedCholesky;
        cholesky = storedCholesky;
        storedCholesky = tmp4;

        logDet = storedLogDet;


    }

    public double calcLogP(Function x) {
        if(requiresRecalculation()){
            refresh();
        }
        double[] xVals = new double[x.getDimension()];
        for(int i = 0;i < xVals.length;i++){
            xVals[i] = x.getArrayValue(i);
        }
		return logPdf(xVals);
    }

	@Override
	public ContinuousDistribution getDistribution() {
		return null;
	}


    public static final String TYPE = "MultivariateNormal";


    public String getType() {
        return TYPE;
    }


    public double[][] getScaleMatrix() {
        return precision;
    }

    public double[] getMeanVector() {
        return mean;
    }

    public Double[][] sample(int size){
        Double[][] samples = new Double[size][];
        try{
            for(int i =0; i < samples.length;i++){
                double[] sample = nextMultivariateNormal();
                Double[] sampleVals = new Double[sample.length];
                for(int j = 0; j < sampleVals.length;j++){
                    sampleVals[j] = sample[j];
                }
                samples[i] = sampleVals;
            }
        }catch(Exception e){
            throw new RuntimeException(e);

        }
        return samples;

    }

    public Double[][] sample(int size, double scale){
        Double[][] samples = new Double[size][];
        try{
            for(int i =0; i < samples.length;i++){
                double[] sample = nextScaledMultivariateNormal(mean,scale);
                Double[] sampleVals = new Double[sample.length];
                for(int j = 0; j < sampleVals.length;j++){
                    sampleVals[j] = sample[j];
                }
                samples[i] = sampleVals;
            }
        }catch(Exception e){
            throw new RuntimeException(e);

        }
        return samples;

    }

    public double[] sample(){
        return nextMultivariateNormal();
    }

    public double[] nextMultivariateNormal() {
        return nextMultivariateNormalCholesky(mean, cholesky, 1.0);
    }

    public double[] nextMultivariateNormal(double[] x) {
        return nextMultivariateNormalCholesky(x, cholesky, 1.0);
    }

    // Scale lives in variance-space
    public double[] nextScaledMultivariateNormal(double[] mean, double scale) {
        return nextMultivariateNormalCholesky(mean, cholesky, Math.sqrt(scale));
    }

    // Scale lives in variance-space
    public void nextScaledMultivariateNormal(double[] mean, double scale, double[] result) {
        nextMultivariateNormalCholesky(mean, cholesky, Math.sqrt(scale), result);
    }


    public static double calculatePrecisionMatrixDeterminate(double[][] precision) {
        try {
            return new Matrix(precision).determinant();
        } catch (IllegalDimension e) {
            throw new RuntimeException(e.getMessage());
        }
    }


    public double logPdf(double[] x, double scale){
        return  logPdf(x, mean, precision, logDet, scale);

    }
    public double logPdf(double[] x) {
        return logPdf(x, mean, precision, logDet, 1.0);
    }

    public static double logPdf(double[] x, double[] mean, double[][] precision,
                                double logDet, double scale) {

        if (logDet == Double.NEGATIVE_INFINITY)
            return logDet;

        final int dim = x.length;
        final double[] delta = new double[dim];
        final double[] tmp = new double[dim];

        for (int i = 0; i < dim; i++) {
            delta[i] = x[i] - mean[i];
        }

        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                tmp[i] += delta[j] * precision[j][i];
            }
        }

        double SSE = 0;

        for (int i = 0; i < dim; i++)
            SSE += tmp[i] * delta[i];

        return dim * logNormalize + 0.5 * (logDet - dim * Math.log(scale) - SSE / scale);   // There was an error here.
        // Variance = (scale * Precision^{-1})
    }

    /* Equal precision, independent dimensions */
    public static double logPdf(double[] x, double[] mean, double precision, double scale) {

        final int dim = x.length;

        double SSE = 0;
        for (int i = 0; i < dim; i++) {
            double delta = x[i] - mean[i];
            SSE += delta * delta;
        }

        return dim * logNormalize + 0.5 * (dim * (Math.log(precision) - Math.log(scale)) - SSE * precision / scale);
    }


    /* Exactly the same as the one above.
     * It's just that in BEAST2 things tend to be objevts rather than primatives. */
    public static double logPdf(Double[] x, Double[] mean, double[][] precision, double logDet, double scale) {

        if (logDet == Double.NEGATIVE_INFINITY)
            return logDet;

        final int dim = x.length;
        final double[] delta = new double[dim];
        final double[] tmp = new double[dim];

        for (int i = 0; i < dim; i++) {
            delta[i] = x[i] - mean[i];
        }

        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                tmp[i] += delta[j] * precision[j][i];
            }
        }

        double SSE = 0;

        for (int i = 0; i < dim; i++)
            SSE += tmp[i] * delta[i];

        return dim * logNormalize + 0.5 * (logDet - dim * Math.log(scale) - SSE / scale);   // There was an error here.
    }

    private static double[][] getInverse(double[][] x) {
        return new SymmetricMatrix(x).inverse().toComponents();
    }

    private static double[][] getCholeskyDecomposition(double[][] variance) {
        double[][] cholesky;
        try {
            cholesky = (new CholeskyDecomposition(variance)).getL();
        } catch (IllegalDimension illegalDimension) {
            throw new RuntimeException("Attempted Cholesky decomposition on non-square matrix");
        }
        return cholesky;
    }



    public static double[] nextMultivariateNormalPrecision(double[] mean, double[][] precision) {
        return nextMultivariateNormalVariance(mean, getInverse(precision));
    }

    public static double[] nextMultivariateNormalVariance(double[] mean, double[][] variance) {
        return nextMultivariateNormalVariance(mean, variance, 1.0);
    }

    public static double[] nextMultivariateNormalVariance(double[] mean, double[][] variance, double scale) {
        return nextMultivariateNormalCholesky(mean, getCholeskyDecomposition(variance), Math.sqrt(scale));
    }

    public static double[] nextMultivariateNormalCholesky(double[] mean, double[][] cholesky) {
        return nextMultivariateNormalCholesky(mean, cholesky, 1.0);
    }

    public static double[] nextMultivariateNormalCholesky(double[] mean, double[][] cholesky, double sqrtScale) {

        double[] result = new double[mean.length];
        nextMultivariateNormalCholesky(mean, cholesky, sqrtScale, result);
        return result;
    }


    public static void nextMultivariateNormalCholesky(double[] mean, double[][] cholesky, double sqrtScale, double[] result) {

        final int dim = mean.length;

        System.arraycopy(mean, 0, result, 0, dim);

        double[] epsilon = new double[dim];
        for (int i = 0; i < dim; i++)
            epsilon[i] = Randomizer.nextGaussian() * sqrtScale;

        for (int i = 0; i < dim; i++) {
            for (int j = 0; j <= i; j++) {
                result[i] += cholesky[i][j] * epsilon[j];
                // caution: decomposition returns lower triangular
            }
        }
    }

    // TODO should be a junit test
    public static void main(String[] args) {
        testPdf();
        testRandomDraws();
    }

    public static void testPdf() {
        double[] start = {1, 2};
        double[] stop = {0, 0};
        double[][] precision = {{2, 0.5}, {0.5, 1}};
        double scale = 0.2;
        System.err.println("logPDF = " + logPdf(start, stop, precision, Math.log(calculatePrecisionMatrixDeterminate(precision)), scale));
        System.err.println("Should = -19.94863\n");

        System.err.println("logPDF = " + logPdf(start, stop, 2, 0.2));
        System.err.println("Should = -24.53529\n");
    }

    public static void testRandomDraws() {

        double[] start = {1, 2};
        double[][] precision = {{2, 0.5}, {0.5, 1}};
        int length = 10000000;


        System.err.println("Random draws (via precision) ...");
        double[] mean = new double[2];
        double[] SS = new double[2];
        double[] var = new double[2];
        double ZZ = 0;
        for (int i = 0; i < length; i++) {
            double[] draw = nextMultivariateNormalPrecision(start, precision);
            for (int j = 0; j < 2; j++) {
                mean[j] += draw[j];
                SS[j] += draw[j] * draw[j];
            }
            ZZ += draw[0] * draw[1];
        }

        for (int j = 0; j < 2; j++) {
            mean[j] /= length;
            SS[j] /= length;
            var[j] = SS[j] - mean[j] * mean[j];
        }
        ZZ /= length;
        ZZ -= mean[0] * mean[1];

        System.err.println("Mean: " + new Vector(mean));
        System.err.println("TRUE: [ 1 2 ]\n");
        System.err.println("MVar: " + new Vector(var));
        System.err.println("TRUE: [ 0.571 1.14 ]\n");
        System.err.println("Covv: " + ZZ);
        System.err.println("TRUE: -0.286");
    }

    public static final double logNormalize = -0.5 * Math.log(2.0 * Math.PI);


}
