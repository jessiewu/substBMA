package beast.evolution.likelihood;

import beast.core.Description;

/**
 * @author Chieh-Hsi Wu
 */
@Description("This class is for testing purpose.")
public class DummyDPTreeLikelihood extends DPTreeLikelihood{
    public void initAndValidate() {
        super.initAndValidate();
    }

    @Override
    public double calculateLogP(){
        double sum =0.0;
        for(NewWVTreeLikelihood treeLik:treeLiks){
            sum+=treeLik.weightSum();
        }
        if(sum != alignment.getSiteCount()){
          throw new RuntimeException("weights do not add up correctly");
        }
        /*else{
            System.out.println("good");
        }*/
        logP = 0.0;
        return logP;
    }
    public double getSiteLogLikelihood(int iCluster, int iSite){
        return 0.0;
    }

    /*public void addTreeLikelihood(){
        //super.addTreeLikelihood();
    }
    public void removeTreeLikelihood(int index){

    }

    protected void updateWeights(){

    } */

    /*public boolean requiresRecalculation(){
        return true;
    }*/
}
