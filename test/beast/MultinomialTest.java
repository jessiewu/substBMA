package beast;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.math.distributions.Multinomial;
import junit.framework.TestCase;

/**
 * @author Chieh-Hsi Wu
 */
public class MultinomialTest extends TestCase{

    interface Instance {

        String getProbs();
        String getX();
        double getFLogX();

        double[] x = new double[]{3,4,5,6};
    }

    Instance test0 = new Instance() {
        public String getProbs (){

            return "0.25 0.25 0.25 0.25";
        }

        public String getX(){
            return "3 4 5 6";
        }

        public double getFLogX(){
            return -4.89440954649;
        }

    };


    Instance test1 = new Instance() {
        public String getProbs (){

            return "0.15 0.23 0.14 0.07 0.25 0.16";
        }

        public String getX(){
            return "3 2 1 6 2 1";
        }

        public double getFLogX(){
            return -13.0155888173;
        }

    };
    
    Instance test2 = new Instance() {
        public String getProbs (){

            return "0.11 0.22 0.12 0.24 0.3 0.01";
        }

        public String getX(){
            return "1 2 4 6 7 0";
        }

        public double getFLogX(){
            return -7.3470894106;
        }

    };

    Instance test3 = new Instance() {
        public String getProbs (){

            return "0.45 0.72 0.64 0.32 0.82 0.16 0.16 0.15 0.08 0.49";
        }

        public String getX(){
            return "1 8 4 0 1 0 2 0 1 3";
        }

        public double getFLogX(){
            return -15.3488523568;
        }

    };

    Instance[] all = new Instance[]{test0,test1,test2};

    public void testMultinomial() throws Exception{
        for(Instance test: all){
            RealParameter probs = new RealParameter();
            IntegerParameter x  = new IntegerParameter();
            probs.initByName(
                    "value",test.getProbs(),
                    "lower","0",
                    "upper","1");
            x.initByName("value", test.getX(),
                    "lower", "0",
                    "upper", "100000000");
            Multinomial multiNom = new Multinomial();
            multiNom.initByName("probs", probs);
            double fLogX = multiNom.calcLogP(x);
            assertEquals(fLogX, test.getFLogX(), 1e-10);

        }
    }







}
