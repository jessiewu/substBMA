package test;

import org.apache.commons.math.distribution.Distribution;
import junit.framework.TestCase;
import beast.math.distributions.IntegerUniformDistribution;
import org.apache.commons.math.distribution.ContinuousDistribution;

/**
 * @author Chieh-Hsi Wu
 */
public class IntegerUniformDistributionTest extends TestCase {
    interface Instance {

        double getUpper();
        double getLower();
        double getExpectedLogDensity();

    }

    Instance test0 = new Instance(){

        public double getUpper(){
            return 10.0;
        }

        public double getLower(){
            return 1.0;
        }

        public double getExpectedLogDensity(){
            return Math.log(1.0/10.0);
        }

    };


    Instance test1 = new Instance(){

        public double getUpper(){
            return 53.0;
        }

        public double getLower(){
            return 28.0;
        }

        public double getExpectedLogDensity(){
            return Math.log(1.0/26.0);
        }

    };



    Instance test2 = new Instance(){

        public double getUpper(){
            return 103.0;
        }

        public double getLower(){
            return 74.0;
        }

        public double getExpectedLogDensity(){
            return Math.log(1.0/30.0);
        }

    };

    Instance[] all = new Instance[]{test0,test1,test2};
    public void testDiscrete() throws Exception{
        for(Instance test: all){
            IntegerUniformDistribution iu = new IntegerUniformDistribution();

            iu.initByName(
                    "lower",test.getLower(),
                    "upper",test.getUpper()
            );

            Distribution iuDist = iu.getDistribution();


            for(int i = (int)test.getLower();i <= (int)test.getUpper();i++){
                ContinuousDistribution c = ((ContinuousDistribution) iuDist);
                assertEquals(
                        ((ContinuousDistribution) iuDist).logDensity((double) i),
                        test.getExpectedLogDensity(),
                        1e-10
                );

            }



        }
    }
    


}
