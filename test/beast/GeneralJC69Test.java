package beast;

import beast.core.parameter.RealParameter;
import beast.evolution.sitemodel.GeneralJC69;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.GeneralSubstitutionModel;
import junit.framework.TestCase;

/**
 * @author Chieh-Hsi Wu
 */
public class GeneralJC69Test extends TestCase {
    public void testGeneralJC69(){
        try{
            for(int stateCount = 2; stateCount < 10; stateCount++){
                Frequencies freqs = new Frequencies();
                Double[] freqVals = new Double[stateCount];
                for(int i = 0; i < freqVals.length;i++){
                    freqVals[i] = 1.0/stateCount;
                }
                RealParameter freqsParam = new RealParameter(freqVals);
                freqs.initByName("frequencies", freqsParam);

                Double[] rateVals = new Double[stateCount*(stateCount - 1)];
                for(int i = 0; i < rateVals.length; i++){
                    rateVals[i] = 1.0;
                }
                RealParameter rates = new RealParameter(rateVals);

                GeneralSubstitutionModel gsm = new GeneralSubstitutionModel();
                gsm.initByName(
                        "rates", rates,
                        "frequencies", freqs
                );
                double[] gsmProbs = new double[stateCount*stateCount];
                gsm.getTransitionProbabilities(null, 0.5, 0, 1.0, gsmProbs);

                GeneralJC69 jc69 = new GeneralJC69();
                jc69.initByName(
                        "stateCount", stateCount
                );
                double[] jc69Probs = new double[stateCount*stateCount];
                jc69.getTransitionProbabilities(null, 0.5, 0, 1.0, jc69Probs);
                for(int i = 0; i < jc69Probs.length; i++){
                    assertEquals(gsmProbs[i],jc69Probs[i],1e-10);

                }
            }

        }catch(Exception e){
            throw new RuntimeException(e);

        }

    }
}
