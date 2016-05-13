package beast;

import beast.core.parameter.RealParameter;
import beast.math.distributions.DirichletDistribution;
import junit.framework.TestCase;

/**
 * @author Chieh-Hsi Wu
 */
public class DirichletDistributionTest extends TestCase {
    interface Instance {

        Double[] getAlpha();
        Double getScale();
        Double[] getX();
        double getFLogX();
    }

    Instance test0 = new Instance() {
        @Override
        public Double[] getAlpha() {
            return new Double[]{0.375, 0.25, 0.125, 0.25};

        }

        @Override
        public Double getScale() {
            return 8.0;
        }

        @Override
        public Double[] getX() {
            return new Double[]{0.35, 0.2, 0.15, 0.3};  //To change body of implemented methods use File | Settings | File Templates.
        }

        @Override
        public double getFLogX() {
            return 2.91895921474808;  //To change body of implemented methods use File | Settings | File Templates.
        }
    };


    Instance[] all = new Instance[]{test0};
    public void testDirichlet(){
        for(Instance test:all){
            double logf1 = DirichletDistribution.logPDF(test.getX(),test.getAlpha(),test.getScale());
            assertEquals(logf1,test.getFLogX(),1e-10);
        }

    }
}
