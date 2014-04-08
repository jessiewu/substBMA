package beast.math.distributions;


import beast.core.Description;
import beast.core.Input;
import beast.core.Valuable;
import beast.core.parameter.RealParameter;
import beast.math.GammaFunction;
import beast.math.matrixAlgebra1.CholeskyDecomposition;
import beast.math.matrixAlgebra1.IllegalDimension;
import beast.math.matrixAlgebra1.Matrix;
import beast.util.Randomizer;
import org.apache.commons.math.distribution.ContinuousDistribution;

/**
 * @author Marc Suchard (Migrated by CHW)
 */
@Description("This class implements the Wishart distribution.")
public class Wishart extends ParametricDistribution implements MultivariateDistribution {
    public Input<RealParameter> scaleMatrixInput = new Input<RealParameter>("scaleMatrix","Scale matrix of the Wishart distribution", Input.Validate.REQUIRED);
    public Input<Integer> dfInput = new Input<Integer>("df","Degrees of freedom of the Wishart distribution", Input.Validate.REQUIRED);

    public static final String TYPE = "Wishart";

    private int df;
    private int dim;
    private double[][] scaleMatrix;
    private double[] Sinv;
    private Matrix SinvMat;
    private double logNormalizationConstant;

    public void initAndValidate(){
        df = dfInput.get();
        RealParameter scaleMatrix1d = scaleMatrixInput.get();
        dim = (int)Math.sqrt(scaleMatrix1d.getDimension());
        scaleMatrix = new double[dim][dim];
        int k = 0;
        for(int i = 0; i < dim;i++){
            for(int j = 0; j < dim; j++){
                scaleMatrix[i][j] =  scaleMatrix1d.getValue(k++);
                //System.out.print(scaleMatrix[i][j]+" ");
            }
            //System.out.println();
        }
        refresh();

    }

    public Wishart(){}


    /**
     * A Wishart distribution class for \nu degrees of freedom and scale matrix S
     * Expectation = \nu * S
     *
     * @param df          degrees of freedom
     * @param scaleMatrix scaleMatrix
     */

    public Wishart(int df, double[][] scaleMatrix) {
        this.df = df;
        this.scaleMatrix = scaleMatrix;
        this.dim = scaleMatrix.length;

        refresh();

    }

    private void refresh(){
        SinvMat = new Matrix(scaleMatrix).inverse();
        double[][] tmp = SinvMat.toComponents();
        Sinv = new double[dim * dim];
        for (int i = 0; i < dim; i++) {
            System.arraycopy(tmp[i], 0, Sinv, i * dim, dim);
        }

        computeNormalizationConstant();

    }

    public Wishart(int dim) { // returns a non-informative (unormalizable) density
        this.df = 0;
        this.scaleMatrix = null;
        this.dim = dim;
        logNormalizationConstant = 1.0;
    }


    private void computeNormalizationConstant() {
        logNormalizationConstant = computeNormalizationConstant(new Matrix(scaleMatrix), df, dim);
    }

    public static double computeNormalizationConstant(Matrix Sinv, int df, int dim) {

        double logNormalizationConstant = 0;
        try {
            logNormalizationConstant = -df / 2.0 * Math.log(Sinv.determinant());
        } catch (IllegalDimension illegalDimension) {
            illegalDimension.printStackTrace();
        }
        //System.out.println("flag1: "+logNormalizationConstant);
        logNormalizationConstant -= df * dim / 2.0 * Math.log(2);
        //System.out.println("flag2: "+logNormalizationConstant);
        logNormalizationConstant -= dim * (dim - 1) / 4.0 * Math.log(Math.PI);
        //System.out.println("flag3: "+logNormalizationConstant);
        for (int i = 1; i <= dim; i++) {
            logNormalizationConstant -= GammaFunction.lnGamma((df + 1 - i) / 2.0);
            //System.out.println((df + 1 - i) / 2.0+" "+logNormalizationConstant);
        }
        //System.out.println("flag4: "+logNormalizationConstant);
        return logNormalizationConstant;
    }


    public String getType() {
        return TYPE;
    }

    public double[][] getScaleMatrix() {
        return scaleMatrix;
    }

    public double[] getMeanVector() {
        return null;
    }

    public void testMe() {

        int length = 100000;

        double save1 = 0;
        double save2 = 0;
        double save3 = 0;
        double save4 = 0;

        for (int i = 0; i < length; i++) {

            double[][] draw = nextWishart();
            save1 += draw[0][0];
            save2 += draw[0][1];
            save3 += draw[1][0];
            save4 += draw[1][1];

        }

        save1 /= length;
        save2 /= length;
        save3 /= length;
        save4 /= length;

        System.err.println("S1: " + save1);
        System.err.println("S2: " + save2);
        System.err.println("S3: " + save3);
        System.err.println("S4: " + save4);


    }

    public int df() {
        return df;
    }

    public double[][] scaleMatrix() {
        return scaleMatrix;
    }


    public double[] sample(){
        return null;
    }
    public Double[][] sample(int size){ throw new RuntimeException("");}

    public double[][] nextWishart() {
        return nextWishart(df, scaleMatrix);
    }

    /**
     * Generate a random draw from a Wishart distribution
     * Follows Odell and Feiveson (1996) JASA 61, 199-203
     * <p/>
     * Returns a random variable with expectation = df * scaleMatrix
     *
     * @param df          degrees of freedom
     * @param scaleMatrix scaleMatrix
     * @return a random draw
     */
    public static double[][] nextWishart(double df, double[][] scaleMatrix) {

        int dim = scaleMatrix.length;
        double[][] draw = new double[dim][dim];

        double[][] z = new double[dim][dim];

        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < i; j++) {
                z[i][j] = Randomizer.nextGaussian();
            }
        }

        for (int i = 0; i < dim; i++)
            z[i][i] = Math.sqrt(Randomizer.nextGamma((df - i) * 0.5, 0.5));   // sqrt of chisq with df-i dfs

        double[][] cholesky = new double[dim][dim];
        for (int i = 0; i < dim; i++) {
            for (int j = i; j < dim; j++)
                cholesky[i][j] = cholesky[j][i] = scaleMatrix[i][j];
        }

        try {
            cholesky = (new CholeskyDecomposition(cholesky)).getL();
            // caution: this returns the lower triangular form
        } catch (IllegalDimension illegalDimension) {
            throw new RuntimeException("Numerical exception in WishartDistribution");
        }

        double[][] result = new double[dim][dim];

        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {     // lower triangular
                for (int k = 0; k < dim; k++)     // can also be shortened
                    result[i][j] += cholesky[i][k] * z[k][j];
            }
        }

        for (int i = 0; i < dim; i++) {           // lower triangular, so more efficiency is possible
            for (int j = 0; j < dim; j++) {
                for (int k = 0; k < dim; k++)
                    draw[i][j] += result[i][k] * result[j][k];   // transpose of 2nd element
            }
        }

        return draw;
    }

    public double calcLogP(Valuable x) throws Exception {
        int vDim = x.getDimension();
        if(vDim == 4){

            double[] val = new double[x.getDimension()];
            for(int i = 0; i < val.length; i++){
                val[i] = x.getArrayValue(i);
            }
            return logPdf2D(val, Sinv, df, dim, logNormalizationConstant);

        }else{

            int matDim = (int)Math.sqrt(vDim);
            double[][] val = new double[matDim][matDim];
            int k = 0;
            for(int i = 0; i < matDim; i++){
                for(int j = 0; j < matDim; j++){
                    val[i][j] = x.getArrayValue(k++);
                    //System.out.print(val[i][j]+" ");
                }
                //System.out.println();

            }
            Matrix W = new Matrix(val);
            double logP = logPdf(W, SinvMat, df, dim, logNormalizationConstant);
            /*System.out.println("logP: "+logP);
            double[][] a = SinvMat.toComponents();
            for(int i = 0; i < SinvMat.rows();i++){
                for(int j = 0;j < SinvMat.columns();j++){
                    System.out.print(a[i][j]+" ");
                }
                System.out.println();
            }
            System.out.println("df: "+df+", dim: "+dim+", logNormalizationConstant: "+logNormalizationConstant); */


            return logP;
        }

    }

    public double logPdf(double[] x) {
        if (x.length == 4) { // bivariate
            return logPdf2D(x, Sinv, df, dim, logNormalizationConstant);
        } else {
            return logPdfSlow(x);
        }
    }



    public double logPdfSlow(double[] x) {
        Matrix W = new Matrix(x, dim, dim);
        return logPdf(W, SinvMat, df, dim, logNormalizationConstant);
    }

    public double logPdfSlow(double[][] x) {
        Matrix W = new Matrix(x);
        return logPdf(W, SinvMat, df, dim, logNormalizationConstant);
    }

    public static double logPdf2D(double[] W, double[] Sinv, int df, int dim, double logNormalizationConstant) {

        final double det = W[0] * W[3] - W[1] * W[2];
        if (det <= 0) {
            return Double.NEGATIVE_INFINITY;
        }

        double logDensity = Math.log(det);
        logDensity *= 0.5 * (df - dim - 1);

        // logDensity -= 0.5 * tr(Sinv %*% W)
        final double trace = Sinv[0] * W[0] + Sinv[1] * W[2] + Sinv[2] * W[1] + Sinv[3] * W[3];
        logDensity -= 0.5 * trace;

        logDensity += logNormalizationConstant;

        return logDensity;
    }
    public ContinuousDistribution getDistribution() {
            return null;
        }



    public static double logPdf(Matrix W, Matrix Sinv, int df, int dim, double logNormalizationConstant) {

        double logDensity = 0;

        try {
            if (!W.isPD())
                return Double.NEGATIVE_INFINITY;

            final double det = W.determinant();

            if (det <= 0) {
                return Double.NEGATIVE_INFINITY;
            }

            logDensity = Math.log(det);

            logDensity *= 0.5;
            logDensity *= df - dim - 1;

            // need only diagonal, no? seems a waste to compute
            // the whole matrix
            if (Sinv != null) {
                Matrix product = Sinv.product(W);

                for (int i = 0; i < dim; i++)
                    logDensity -= 0.5 * product.component(i, i);
            }

        } catch (IllegalDimension illegalDimension) {
            illegalDimension.printStackTrace();
        }

        logDensity += logNormalizationConstant;
        return logDensity;
    }

    public static void testBivariateMethod() {
        System.out.println("Testing new computations ...");
        Wishart wd = new Wishart(5, new double[][]{{2.0, -0.5}, {-0.5, 2.0}});
        double[] W = new double[]{4.0, 1.0, 1.0, 3.0};
        System.out.println("Fast logPdf = " + wd.logPdf(W));
        System.out.println("Slow logPdf = " + wd.logPdfSlow(W));
    }

    public static void main(String[] argv) {
        try{
            /*Wishart wd = new Wishart(2, new double[][]{{500.0}});*/
            Wishart wd = new Wishart();
            wd.initByName(
                    "scaleMatrix", new RealParameter(new Double[]{500.0}),
                    "df", 2);

            // The above is just an approximation
            Gamma gd = new Gamma();
            gd.initByName(
                    "alpha", new RealParameter(new Double[]{1.0 / 1000.0}),
                    "beta", new RealParameter(new Double[]{1000.0})
            );

            RealParameter x = new RealParameter(new Double[]{1.0});
            System.out.println("Wishart, df=2, scale = 500, PDF(1.0): " + wd.calcLogP(x));
            System.out.println("Gamma, shape = 1/1000, scale = 1000, PDF(1.0): " + gd.logDensity(x.getValue()));

            //wd = new Wishart(4, new double[][]{{5.0}});
            wd.initByName(
                    "scaleMatrix", new RealParameter(new Double[]{5.0}),
                    "df", 4);
            gd.initByName(
                    "alpha",new RealParameter(new Double[]{2.0}),
                    "beta", new RealParameter(new Double[]{10.0})
            );
            x = new RealParameter(new Double[]{1.0});
            System.out.println("Wishart, df=4, scale = 5, PDF(1.0): " + wd.calcLogP(x));
            System.out.println("Gamma, shape = 1/1000, scale = 10, PDF(1.0): " + gd.logDensity(x.getValue()));
            // These tests show the correspondence between a 1D Wishart and a Gamma

            wd = new Wishart(1);
            x = new RealParameter(new Double[]{0.1});
            System.out.println("Wishart, uninformative, PDF(0.1): " + wd.calcLogP(x));
            x = new RealParameter(new Double[]{1.0});
            System.out.println("Wishart, uninformative, PDF(1.0): " + wd.calcLogP(x));
            x = new RealParameter(new Double[]{10.0});
            System.out.println("Wishart, uninformative, PDF(10.0): " + wd.calcLogP(x));
            // These tests show the correspondence between a 1D Wishart and a Gamma
            testBivariateMethod();
        }catch(Exception e){
            throw new RuntimeException(e);
        }

    }
}
