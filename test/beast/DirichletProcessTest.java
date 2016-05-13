package beast;

import beast.math.distributions.*;
import beast.core.parameter.*;
import junit.framework.TestCase;

/**
 * @author Chieh-Hsi Wu
 */
public class DirichletProcessTest extends TestCase {

    interface Instance {
        public void setup();
        public ParametricDistribution getBaseDistribution();
        public ParameterList getParameterList();
        public RealParameter getAlpha() ;
        public DPValuable getDPValuable();
        public double getExpectedLogP();
    }

    Instance test0 = new Instance(){
        Normal normal;
        ParameterList paramList;
        RealParameter alpha;
        DPValuable dpValuable;
        DPPointer dpPointer;
        public void setup(){
            normal = new Normal();
            paramList = new ParameterList();
            alpha = new RealParameter();
            dpValuable = new DPValuable();
            try{
                normal.initByName(
                        "mean", "0",
                        "sigma", "3"
                );

                QuietRealParameter param1 = new QuietRealParameter(new Double[]{-0.3});
                QuietRealParameter param2 = new QuietRealParameter(new Double[]{-0.1});
                QuietRealParameter param3 = new QuietRealParameter(new Double[]{0.5});
                paramList.initByName(
                        "parameter", param1,
                        "parameter", param2,
                        "parameter", param3
                );

                IntegerParameter assignment = new IntegerParameter();
                assignment.initByName("value","0 1 0 2 0 0 1 0");

                dpPointer = new DPPointer();
                dpPointer.initByName(
                        "uniqueParameter",param1,
                        "uniqueParameter",param2,
                        "uniqueParameter",param3,
                        "initialAssignment",assignment
                );

                alpha.initByName("value","0.5");
                dpValuable = new DPValuable();
                dpValuable.initByName(
                        "paramList",paramList,
                        "pointers", dpPointer
                );
            }catch(Exception e){
                throw new RuntimeException(e);
            }
        }

        public ParametricDistribution getBaseDistribution(){
            return normal;
        }

        public ParameterList getParameterList() {
            return paramList;
        }

        public RealParameter getAlpha () {
            return alpha;
        }

        public DPValuable getDPValuable(){
            dpPointer.setEverythingDirty(true);
            return dpValuable;

        }

        public double getExpectedLogP(){
            return -13.9503869358;
        }


    };

    Instance test1 = new Instance(){
        Normal normal;
        ParameterList paramList;
        RealParameter alpha;
        DPValuable dpValuable;
        DPPointer dpPointer;
        public void setup(){
            normal = new Normal();
            paramList = new ParameterList();
            alpha = new RealParameter();
            dpValuable = new DPValuable();
            try{
                normal.initByName(
                        "mean", "0",
                        "sigma", "3"
                );

                QuietRealParameter param1 = new QuietRealParameter(new Double[]{-0.3,0.2,-0.1});
                QuietRealParameter param2 = new QuietRealParameter(new Double[]{-0.1,0.03,0.15});
                QuietRealParameter param3 = new QuietRealParameter(new Double[]{0.5,0.4,-0.01});

                paramList.initByName(
                        "parameter", param1,
                        "parameter", param2,
                        "parameter", param3
                );
                alpha.initByName("value","0.5");
                IntegerParameter assignment = new IntegerParameter();
                assignment.initByName("value","0 1 0 2 0 0 1 0");
                dpPointer = new DPPointer();
                dpPointer.initByName(
                        "uniqueParameter",param1,
                        "uniqueParameter",param2,
                        "uniqueParameter",param3,
                        "initialAssignment",assignment
                );
                dpValuable = new DPValuable();
                dpValuable.initByName(
                        "paramList",paramList,
                        "pointers", dpPointer
                );
            }catch(Exception e){
                throw new RuntimeException(e);
            }
        }

        public ParametricDistribution getBaseDistribution(){
            return normal;
        }

        public ParameterList getParameterList() {
            return paramList;
        }

        public RealParameter getAlpha () {
            return alpha;
        }

        public DPValuable getDPValuable(){
            dpPointer.setEverythingDirty(true);
            return dpValuable;

        }

        public double getExpectedLogP(){
            return -26.0686640892;
        }


    };

    Instance test2 = new Instance(){
        MultivariateNormal mvnorm;
        ParameterList paramList;
        RealParameter alpha;
        DPValuable dpValuable;
        DPPointer dpPointer;
        public void setup(){
            mvnorm = new MultivariateNormal();
            paramList = new ParameterList();
            alpha = new RealParameter();
            try{
                RealParameter mean = new RealParameter();
                mean.initByName("value", "0.04 -0.01 0.1");


                RealParameter precisionMatrix = new RealParameter();
                precisionMatrix.initByName(
                        "value",
                        "1.0359869138495 -0.1690294438386 -0.0436205016358 -0.1690294438386  0.8696837513631 -0.0981461286805 -0.0436205016358 -0.0981461286805  1.2649945474373");
                mvnorm.initByName(
                        "mean", mean,
                        "precision", precisionMatrix
                );

                QuietRealParameter param1 = new QuietRealParameter(new Double[]{-0.3,0.2,-0.1});
                QuietRealParameter param2 = new QuietRealParameter(new Double[]{-0.1,0.03,0.15});
                QuietRealParameter param3 = new QuietRealParameter(new Double[]{0.5,0.4,-0.01});
                paramList.initByName(
                        "parameter", param1,
                        "parameter", param2,
                        "parameter", param3
                );
                alpha.initByName("value","0.5");
                IntegerParameter assignment = new IntegerParameter();
                assignment.initByName("value","0 1 0 2 0 0 1 0");

                dpPointer = new DPPointer();
                dpPointer.initByName(
                        "uniqueParameter",param1,
                        "uniqueParameter",param2,
                        "uniqueParameter",param3,
                        "initialAssignment",assignment
                );
                dpValuable = new DPValuable();
                dpValuable.initByName(
                        "paramList",paramList,
                        "pointers", dpPointer
                );
            }catch(Exception e){
                throw new RuntimeException(e);
            }
        }

        public ParametricDistribution getBaseDistribution(){
            return mvnorm;
        }

        public ParameterList getParameterList() {
            return paramList;
        }

        public RealParameter getAlpha () {
            return alpha;
        }

        public DPValuable getDPValuable(){
            dpPointer.setEverythingDirty(true);
            return dpValuable;

        }

        public double getExpectedLogP(){
            return -16.3149436859;
        }


    };

    Instance[] all = new Instance[]{test0,test1,test2};

    public void testDirichletProcessPrior(){
        try{
            for(Instance test:all){
                test.setup();
                DPValuable dpValuable = test.getDPValuable();
                RealParameter alpha = test.getAlpha();
                ParameterList paramList = test.getParameterList();
                ParametricDistribution distr = test.getBaseDistribution();
                DirichletProcess dp = new DirichletProcess();

                dp.initByName(
                        "dpVal", dpValuable,
                        "alpha", alpha,
                        "baseDistr", distr
                );

                ParameterListPrior dpp = new ParameterListPrior();
                dpp.initByName(
                        "distr", dp,
                        "xList", paramList,
                        "applyToList", true
                );

                //System.err.println(dpp.calculateLogP());
                assertEquals(dpp.calculateLogP(), test.getExpectedLogP(), 5e-10);

            }

        }catch(Exception e){
            throw new RuntimeException(e);

        }
    }
}
