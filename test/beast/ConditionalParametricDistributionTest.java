package beast;

import beast.core.parameter.RealParameter;
import beast.math.distributions.*;
import junit.framework.TestCase;

/**
 * Created by IntelliJ IDEA.
 * User: cwu080
 * Date: 9/09/13
 * Time: 11:23 AM
 * To change this template use File | Settings | File Templates.
 */
public class ConditionalParametricDistributionTest extends TestCase {
    interface Instance {

        String getIntegerConditionValues();
        ParametricDistribution[] getParametricDistributions();
        double[] getXValues();
        int[] getConditionValues();
        double[] getLogProbs();
    }

    Instance test0 = new Instance() {
        public String getIntegerConditionValues(){

            return "3 2 5 7";
        }

        public ParametricDistribution[] getParametricDistributions(){
            try{
                ParametricDistribution[] distr = new ParametricDistribution[4];
                distr[0] = new Normal();
                distr[0].initByName(
                        "mean", new RealParameter(new Double[]{0.5}),
                        "sigma", new RealParameter(new Double[]{2.0})
                );

                distr[1] = new Normal();
                distr[1].initByName(
                        "mean", new RealParameter(new Double[]{-1.0}),
                        "sigma", new RealParameter(new Double[]{0.5})
                );

                distr[2] = new Gamma();
                distr[2].initByName(
                        "alpha", new RealParameter(new Double[]{2.0}),
                        "beta", new RealParameter(new Double[]{0.5})
                );


                distr[3] = new Gamma();
                distr[3].initByName(
                        "alpha", new RealParameter(new Double[]{0.1}),
                        "beta", new RealParameter(new Double[]{10.0})
                );

                return distr;
            }catch(Exception e ){
                throw new RuntimeException(e);
            }
        }

        public double[] getXValues(){
            return new double[]{6.382002, -1.821,-0.8656,0.0001,-0.3005,0.3214,0.69969,-1.055, -0.7303,2.97367};

        }

        public int[] getConditionValues(){
            return new int[]{7,3,2,7,3,2,5,3,3,7};

        }

        public double[] getLogProbs(){
            return new double[]{
                    -4.789305018147259,-2.285465838764618,-0.261918072644727,5.806325173744954,
                    -1.692185745014618,-3.717987272644728,-0.370203538051884,-1.914338838764618,
                    -1.801290475014618,-3.761155353225738
            };
        }

    };

    Instance[] tests = new Instance[]{test0};

    public void testConditionalParametricDistribution(){
        try{
            for(Instance test:tests){
                String conditionalVal = test.getIntegerConditionValues();
                ParametricDistribution[] distrs = test.getParametricDistributions();
                double[] x = test.getXValues();
                int[] y = test.getConditionValues();
                double[] logP = test.getLogProbs();

                ConditionalParametricDistribution cpd = new ConditionalParametricDistribution();
                cpd.initByName(
                        "integerConditionValues", conditionalVal,
                        "distr", distrs[0],
                        "distr", distrs[1],
                        "distr", distrs[2],
                        "distr", distrs[3]

                );

                for(int i = 0; i < x.length; i++){
                    assertEquals(cpd.calcLogP(new RealParameter(new Double[]{x[i]}), y[i]),logP[i], 1e-10);
                }


            }
        }catch(Exception e){
            throw new RuntimeException(e);
        }

    }
}
