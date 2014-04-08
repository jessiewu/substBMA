package beast.evolution.likelihood;

import beast.core.Description;

/**
 * @author Chieh-Hsi Wu
 */
@Description("This class is for testing purpose.")
public class DummyDPSepTreeLikelihood extends DPSepTreeLikelihood{
    

    @Override
    public double calculateLogP(){
        return 0.0;
    }
    public double getSiteLogLikelihood(int iCluster, int iSite){
        return 0.0;
    }

    public double getSiteLogLikelihood(int changeType, int iCluster, int iSite){
        return 0.0;
    }

    
}
