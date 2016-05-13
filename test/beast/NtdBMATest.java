package beast;

import beast.core.parameter.QuietRealParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.substitutionmodel.NtdBMA;
import junit.framework.TestCase;

/**
 * @author Chieh-Hsi Wu
 *
 * JUnit test for NtdBMA model
 */
public class NtdBMATest extends TestCase {
    interface Instance {
        String getPi();

        double getLogKappa();

        double getLogTN();

        double getLogAC();

        double getLogAT();
        double getLogGC();
        double getLogGT();
        QuietRealParameter getModelChoose() throws Exception;

        double getDistance();

        double[] getExpectedResult();
    }

    
    //A HKY model
    Instance test0 = new Instance() {
        public String getPi() {
            return "0.25 0.25 0.25 0.25";
        }

        public double getLogKappa() {
            return Math.log(2);
        }

        public double getLogTN(){
            return Math.log(1.2);
        }

        public double getLogAC(){
            return Math.log(0.5);
        }

        public double getLogAT(){
            return Math.log(0.5);
        }

        public double getLogGC(){
            return Math.log(0.5);
        }
        public double getLogGT(){
            return Math.log(0.5);
        }

        public QuietRealParameter getModelChoose() throws Exception{
            QuietRealParameter modelChoose = new QuietRealParameter();
            modelChoose.initByName(
                    "value", "3",
                    "lower", "0",
                    "upper", "5");
            return modelChoose;
        }

        public double getDistance() {
            return 0.1;
        }

        public double[] getExpectedResult() {
            return new double[]{
                    0.906563342722, 0.023790645491, 0.045855366296, 0.023790645491,
                    0.023790645491, 0.906563342722, 0.023790645491, 0.045855366296,
                    0.045855366296, 0.023790645491, 0.906563342722, 0.023790645491,
                    0.023790645491, 0.045855366296, 0.023790645491, 0.906563342722
            };
        }
    };

    //A TN93 model
    Instance test1 = new Instance() {
        public String getPi() {
            return "0.1 0.2 0.3 0.4";
        }

        public double getLogKappa() {
            return Math.log(3);
        }

        public double getLogTN(){
            return Math.log(1.5);
        }

        public double getLogAC(){
            return Math.log(0.5);
        }

        public double getLogAT(){
            return Math.log(0.5);
        }

        public double getLogGC(){
            return Math.log(0.5);
        }
        public double getLogGT(){
            return Math.log(0.5);
        }

        public QuietRealParameter getModelChoose() throws Exception{
            QuietRealParameter modelChoose = new QuietRealParameter();
            modelChoose.initByName(
                    "value", "4",
                    "upper", "0",
                    "lower", "5");
            return modelChoose;
        }


        public double getDistance() {
            return 0.1;
        }

        public double[] getExpectedResult() {
            return new double[]{
                0.8978002167574,0.0139801109390,0.0602594504256,0.027960221878,
                0.0069900554695,0.8565503184732,0.0209701664085,0.115489459649,
                0.0200864834752,0.0139801109390,0.9379731837078,0.027960221878,
                0.0069900554695,0.0577447298244,0.0209701664085,0.914295048298
            };
        }
    };

    //GTR example
    Instance test2 = new Instance() {
        public String getPi() {
            return "0.20 0.30 0.25 0.25";
        }


        public double getLogKappa() {
            return Math.log(3);
        }

        public double getLogTN(){
            return Math.log(1.5);
        }

        public double getLogAC(){
            return Math.log(1.2);
        }

        public double getLogAT(){
            return Math.log(0.6);
        }

        public double getLogGC(){
            return Math.log(0.5);
        }
        public double getLogGT(){
            return Math.log(0.8);
        }

        public QuietRealParameter getModelChoose() throws Exception{
            QuietRealParameter modelChoose = new QuietRealParameter();
            modelChoose.initByName(
                    "value", "5",
                    "lower", "0",
                    "upper", "5");
            return modelChoose;
        }


        public double getDistance() {
            return 0.1;
        }

        public double[] getExpectedResult() {
            return new double[]{
                    0.91402704834815,0.0244356629701,0.05034219302632,0.0111950956554,
                    0.01629044198007,0.9014141683003,0.00940551704632,0.0728898726733,
                    0.04027375442105,0.0112866204556,0.93135093054649,0.0170886945769,
                    0.00895607652434,0.0874678472080,0.01708869457688,0.8864873816908
            };
        }
    };

    Instance test3 = new Instance() {
        public String getPi() {
            return "0.25 0.25 0.25 0.25";
        }


        public double getLogKappa() {
            return Math.log(2);
        }

        public double getLogTN(){
            return Math.log(1.2);
        }

        public double getLogAC(){
            return Math.log(0.5);
        }

        public double getLogAT(){
            return Math.log(0.5);
        }

        public double getLogGC(){
            return Math.log(0.5);
        }
        public double getLogGT(){
            return Math.log(0.5);
        }

        public QuietRealParameter getModelChoose() throws Exception{
            QuietRealParameter modelChoose = new QuietRealParameter();
            modelChoose.initByName(
                    "value", "3",
                    "lower", "0",
                    "upper", "5");
            return modelChoose;
        }


        public double getDistance() {
            return 1.8;
        }

        public double[] getExpectedResult() {
            return new double[]{
                    0.324927478425,0.208675277945,0.257721965686,0.208675277945,
                    0.208675277945,0.324927478425,0.208675277945,0.257721965686,
                    0.257721965686,0.208675277945,0.324927478425,0.208675277945,
                    0.208675277945,0.257721965686,0.208675277945,0.324927478425
            };
        }
    };

    Instance test4 = new Instance() {

        public String getPi() {
            return "0.1 0.2 0.3 0.4";
        }


        public double getLogKappa() {
            return Math.log(3);
        }

        public double getLogTN(){
            return Math.log(1.5);
        }

        public double getLogAC(){
            return Math.log(0.5);
        }

        public double getLogAT(){
            return Math.log(0.5);
        }

        public double getLogGC(){
            return Math.log(0.5);
        }
        public double getLogGT(){
            return Math.log(0.5);
        }

        public QuietRealParameter getModelChoose() throws Exception{
            QuietRealParameter modelChoose = new QuietRealParameter();
            modelChoose.initByName(
                    "value", "4",
                    "lower", "0",
                    "upper", "5");
            return modelChoose;
        }



        public double getDistance() {
            return 2.5;
        }

        public double[] getExpectedResult() {
            return new double[]{
                    0.1532752904967,0.167321310649,0.344760777556,0.334642621298,
                    0.0836606553245,0.224212046022,0.250981965973,0.441145332680,
                    0.1149202591854,0.167321310649,0.383115808868,0.334642621298,
                    0.0836606553245,0.220572666340,0.250981965973,0.444784712362
            };
        }
    };

    //GTR example
    Instance test5 = new Instance() {
        public String getPi() {
            return "0.20 0.30 0.25 0.25";
        }


        public double getLogKappa() {
            return Math.log(3);
        }

        public double getLogTN(){
            return Math.log(1.5);
        }

        public double getLogAC(){
            return Math.log(1.2);
        }

        public double getLogAT(){
            return Math.log(0.6);
        }

        public double getLogGC(){
            return Math.log(0.5);
        }
        public double getLogGT(){
            return Math.log(0.8);
        }

        public QuietRealParameter getModelChoose() throws Exception{
            QuietRealParameter modelChoose = new QuietRealParameter();
            modelChoose.initByName(
                    "value", "5",
                    "upper", "0",
                    "lower", "5");
            return modelChoose;
        }
       

        public double getDistance() {
            return 2.5;
        }

        public double[] getExpectedResult() {
            return new double[]{
                    0.267057968497,0.239880956276,0.297062866519,0.195998208708,
                    0.159920637517,0.359793005963,0.186316991846,0.293969364674,
                    0.237650293216,0.223580390216,0.347708298048,0.191061018520,
                    0.156798566967,0.352763237608,0.191061018520,0.299377176905
            };
        }
    };

    Instance[] all = {test0,test1,test2,test3,test4,test5};

    public void testNtdBMA() throws Exception{
        for (Instance test : all) {

            QuietRealParameter logKappa = new QuietRealParameter();
            QuietRealParameter logTN = new QuietRealParameter();
            QuietRealParameter logAC = new QuietRealParameter();
            QuietRealParameter logAT = new QuietRealParameter();
            QuietRealParameter logGC = new QuietRealParameter();
            QuietRealParameter logGT = new QuietRealParameter();

            logKappa.init(Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, ""+test.getLogKappa(),1);
            logTN.init(Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, ""+test.getLogTN(),1);
            logAC.init(Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, ""+test.getLogAC(),1);
            logAT.init(Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, ""+test.getLogAT(),1);
            logGC.init(Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, ""+test.getLogGC(),1);
            logGT.init(Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, ""+test.getLogGT(),1);
            QuietRealParameter modelChoose = test.getModelChoose();
           
            QuietRealParameter f = new QuietRealParameter();
            f.init(0.0,1.0,test.getPi(),4);

            
            NtdBMA ntdBMA = new NtdBMA();
            ntdBMA.initByName(
                    "logKappa",logKappa,
                    "logTN", logTN,
                    "logAC",logAC,
                    "logAT",logAT,
                    "logGC",logGC,
                    "modelChoose",modelChoose,
                    "frequenciesParameter",f
            );

            double distance = test.getDistance();

            double[] mat = new double[4 * 4];
            ntdBMA.getTransitionProbabilities(null,distance, 0,1.0, mat);
            final double[] result = test.getExpectedResult();

            for (int k = 0; k < mat.length; ++k) {
                assertEquals(mat[k], result[k], 5e-10);
            }
        }
    }
}
