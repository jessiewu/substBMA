package beast.evolution.sitemodel;

import beast.core.Description;
import beast.core.Input;
import beast.evolution.substitutionmodel.EigenDecomposition;
import beast.evolution.substitutionmodel.GeneralSubstitutionModel;
import beast.evolution.tree.Node;

/**
 * @author Chieh-Hsi Wu
 */
@Description("A more general version of JC69.")
public class GeneralJC69 extends GeneralSubstitutionModel {

    public Input<Integer> stateCountInput = new Input<Integer>(
            "stateCount",
            "The number of states there is in the data to be modeled",
            Input.Validate.REQUIRED
    );

    private int stateCount;
    private double[] frequencies;
    public GeneralJC69(){
        frequenciesInput.setRule(Input.Validate.OPTIONAL);
        m_rates.setRule(Input.Validate.OPTIONAL);

    }

    public void initAndValidate(){
        stateCount = stateCountInput.get();
        double freqVal = 1.0/stateCount;
        frequencies = new double[stateCount];
        for(int i = 0; i < frequencies.length;i++){
            frequencies[i] = freqVal;
        }
    }

    public double[] getFrequencies() {
        return frequencies;
    }



    public void getTransitionProbabilities(Node node, double fStartTime, double fEndTime, double fRate, double[] matrix){

        double distance = (fStartTime-fEndTime)*fRate;
        double pii = 1.0/stateCount + (stateCount - 1.0)/stateCount*Math.exp(-stateCount/(stateCount - 1.0)*distance);
        double pij = 1.0/stateCount - 1.0/stateCount*Math.exp(-stateCount/(stateCount - 1.0)*distance);
        //System.out.println(matrix.length);
        int k = 0;
        for(int i = 0; i < stateCount; i++){
            for(int j = 0; j < stateCount; j++){
                if(i == j){
                    matrix[k++] = pii;
                }else{
                    matrix[k++] = pij;
                }
            }
        }
        //System.out.println("distance: "+distance+" "+matrix[1]);
    }
    public EigenDecomposition getEigenDecomposition(Node node){
        return null;

    }
}
