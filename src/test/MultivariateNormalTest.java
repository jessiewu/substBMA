package test;

import beast.core.parameter.RealParameter;
import beast.math.distributions.MultivariateNormal;
import beast.math.distributions.Normal;
import junit.framework.TestCase;

/**
 * @author Chieh-Hsi Wu
 */
public class MultivariateNormalTest extends TestCase {
    interface Instance {

        String getMean();
        RealParameter getPrecision() throws Exception;
        String[] getX();
        double[] getFLogX();
    }

    Instance test0 = new Instance() {
        public String getMean(){

            return "1.2 2.1 1.4 1.6";
        }

        public RealParameter getPrecision()throws Exception{



            RealParameter precision = new RealParameter();
            precision.initByName(
                    "value",
                    "1.0000000 0.9529819 0.3554853 0.1695793 " +
                            "0.9529819 2.0000000 0.3658894 0.6902582 " +
                            "0.3554853 0.3658894 1.5000000 0.8099044 " +
                            "0.1695793 0.6902582 0.8099044 1.3000000"
            );
            //System.out.println(precision==null);
            return precision;
        }

        public String[] getX(){
            return new String[]{
                    "1.6 2.3 1.6 0.9",
                    "0.16926 4.51521 0.18029 0.90557",
                    "2.88483 1.31626 -0.24917 2.6704"
            };
        }

        public double[] getFLogX(){
            return new double[]{
                    -4.04250951908,
                    -8.15220554717,
                    -5.0541914572
            };
        }

    };


    Instance test1 = new Instance() {
        public String getMean(){

            return "7 9 5 10 2 9";
        }

        public RealParameter getPrecision()throws Exception{
            RealParameter precision = new RealParameter();
            precision.initByName(
                    "value",
                    "8 6 4 1 6 3 6 12 6 1 6 3 4 6 10 6 7 7 1 1 6 21 9 10 6 6 7 9 22 3 3 3 7 10 3 12"
            );

            return precision;
        }

        public String[] getX(){
            return new String[]{
                    "7.841304 8.933441 5.031969 10.48876 1.645417 8.288646",
                    "6.22691 8.99749 5.23868 9.79864 2.35423 8.77663",
                    "7.47852 8.62351 4.60867 9.6396 2.33618 9.53719"
            };
        }

        public double[] getFLogX(){
            return new double[]{-1.31033090602,-1.68784377673,-1.33586909163};
        }

    };

    Instance test2 = new Instance() {
        public String getMean(){

            return "0.88773 0.20646 0.88280 0.57498 0.64055 0.09056";
        }

        public RealParameter getPrecision()throws Exception{
            RealParameter precision = new RealParameter();
            precision.initByName(
                    "value",
                            "2.87755 0.86139 0.04848 0.34234 0.73282 0.64613 " +
                            "0.86139 1.32000 0.60056 0.40075 0.62026 0.78024 " +
                            "0.04848 0.60056 1.92977 0.69593 0.16473 0.43684 " +
                            "0.34234 0.40075 0.69593 1.97503 0.03434 0.07275 " +
                            "0.73282 0.62026 0.16473 0.03434 1.89613 0.53328 " +
                            "0.64613 0.78024 0.43684 0.07275 0.53328 1.21555"
            );


            return precision;
        }

        public String[] getX(){
            return new String[]{
                    "0.2261449 1.752875 0.8292972 0.4661388 0.8885903 -0.5040412",
                    "0.6499052 0.7368242 -0.4418067 0.0197093 0.9622208 -0.4388211",
                    "0.8624 -0.46159 0.97389 -0.1198 1.12903 0.52234"

            };
        }

        public double[] getFLogX(){
            return new double[]{
                    -5.53715146247,
                    -6.94744854554,
                    -5.3236834678
            };
        }

    };

   

    Instance[] all = new Instance[]{test0,test1,test2};

    public void testMultivariateNormal() throws Exception{
        for(Instance test: all){
            RealParameter mean = new RealParameter();

            mean.initByName(
                    "value",test.getMean(),
                    "lower",Double.NEGATIVE_INFINITY,
                    "upper",Double.POSITIVE_INFINITY);

            RealParameter precision = test.getPrecision();

            MultivariateNormal multiNorm = new MultivariateNormal();
            multiNorm.initByName(
                    "mean",mean,
                    "precision", precision);

            String[] xStrs = test.getX();
            double[] expectedFLogX = test.getFLogX();
            for(int i = 0; i < xStrs.length;i++){
                RealParameter x  =  new RealParameter();
                x.initByName(
                        "value",xStrs[i],
                        "lower",Double.NEGATIVE_INFINITY,
                        "upper",Double.POSITIVE_INFINITY);
                double fLogX1 = multiNorm.calcLogP(x);
                assertEquals(fLogX1, expectedFLogX[i], 1e-10);


                double fLogX2 =  MultivariateNormal.logPdf(
                        x.getValues(),
                        mean.getValues(),
                        multiNorm.getScaleMatrix(),
                        Math.log(MultivariateNormal.calculatePrecisionMatrixDeterminate(multiNorm.getScaleMatrix())),
                        1.0
                );

                assertEquals(fLogX2, expectedFLogX[i], 1e-10);
            }




        }
    }

}
