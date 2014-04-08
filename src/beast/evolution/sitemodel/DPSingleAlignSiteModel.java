package beast.evolution.sitemodel;

import beast.core.Description;
import beast.core.parameter.ChangeType;

import java.util.ArrayList;

/**
 * @author Chieh-Hsi Wu
 */
@Description("Parent class of all DP site models that can only handle a single alignment.")
public abstract class DPSingleAlignSiteModel extends DPSiteModel {

    public abstract int getRemovedIndex();
    public abstract int getLastDirtySite();
    protected abstract void addSiteModel();
    public ChangeType getChangeType(){
        return changeType;
    }

    public abstract int getPrevCluster(int index);

    public abstract int getCurrCluster(int index);



    public abstract int getDirtySiteModelIndex();

    public QuietSiteModel getSiteModelOfSiteIndex(int siteIndex){
        return siteModels.get(getCurrCluster(siteIndex));

    }




}
