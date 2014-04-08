package beast.evolution.sitemodel;

import beast.core.Description;
import beast.core.Input;
import beast.core.MCMCNodeFactory;
import beast.core.PluginList;
import beast.core.parameter.*;
import beast.evolution.substitutionmodel.SubstitutionModel;

import java.util.ArrayList;

/**
 * @author Chieh-Hsi Wu
 */
@Description("This class facilitates partition selection by DPP on site rates.")

public class DPRateSiteModel extends DPSingleAlignSiteModel implements PluginList {

    protected DPPointer ratePointers;
    protected ParameterList rateList;



    //ParameterList
    public Input<ParameterList> rateListInput = new Input<ParameterList>(
            "ratesList",
            "A list of unique site rate values",
            Input.Validate.REQUIRED
    );


    //assignment
    public Input<DPPointer> ratePointersInput = new Input<DPPointer>(
            "ratesPointers",
            "array which points a set of unique parameter values",
            Input.Validate.REQUIRED
    );


    public Input<SubstitutionModel.Base> substModelInput = new Input<SubstitutionModel.Base>(
            "substModel",
            "substitution model along branches in the beast.tree",
            Input.Validate.REQUIRED
    );


    private SubstitutionModel.Base substModel;
    private int[] pointerIndices;
    public void initAndValidate() throws Exception{
        super.initAndValidate();
        substModel = substModelInput.get();
        ratePointers = ratePointersInput.get();
        rateList = rateListInput.get();

        int dim = rateList.getDimension();

        for(int i = 0;i < dim; i++){
            QuietRealParameter muParameter = rateList.getParameter(ratePointers.indexInList(i,rateList));
            QuietSiteModel siteModel = new QuietSiteModel(substModel, muParameter);
            siteModels.add(siteModel);
        }
        pointerIndices = new int[ratePointers.getDimension()];

    }



    protected void addSiteModel(){
        try{
            QuietRealParameter muParameter = rateList.getParameter(rateList.getDimension()-1);
            QuietSiteModel siteModel = new QuietSiteModel(substModel, muParameter);
            siteModels.add(siteModel);
        }catch(Exception e){
            throw new RuntimeException(e);
        }

    }

    public int getPrevCluster(int index){
        return ratePointers.storedIndexInList(index,rateList);
    }

    public int getCurrCluster(int index){
        return ratePointers.indexInList(index,rateList);
    }

    public int getPrevCategoryIDNumber(int siteIndex){
        return ratePointers.getStoredParameterIDNumber(siteIndex);
    }

    public int getCurrCategoryIDNumber(int siteIndex){
        return ratePointers.getParameterIDNumber(siteIndex);
    }



    public int getLastDirtySite(){
        return ratePointers.getLastDirty();
    }

    public int getRemovedIndex(){
        return rateList.getRemovedIndex();
    }

    public int getLastAddedIndex(){
        return rateList.getLastAddedIndex();
    }

    public int getDirtySiteModelIndex(){
        return rateList.getDirtyIndex();
    }

    public int[] getPointerIndices(){
        return pointerIndices;
    }

    public int[] getSwappedSites(){
        return ratePointers.getSwappedSites();
    }

    public int[] getLastDirtySites(){
        return ratePointers.getLastDirtySites();
    }

    public void resetAllPointerIndices(){
        for(int i = 0; i < pointerIndices.length;i++){
            pointerIndices[i] = ratePointers.indexInList(i,rateList);
        }
    }

    int prevCluster = -1;
    public void setupPointerIndices(){
        //System.err.println("changeType: "+changeType);
        if(changeType == ChangeType.ADDED || changeType == ChangeType.POINTER_CHANGED){
            int changedIndex = ratePointers.getLastDirty();
            prevCluster = pointerIndices[changedIndex];
            pointerIndices[changedIndex] = ratePointers.indexInList(changedIndex,rateList);
        }else if (changeType == ChangeType.REMOVED){
            resetAllPointerIndices();
        }else if (changeType == ChangeType.ALL){
            resetAllPointerIndices();
        }
    }

    public boolean requiresRecalculation(){

        boolean recalculate = false;
        //System.err.println("dirty0");
        if(substModel.isDirtyCalculation()){
            //System.err.println("substModel is dirty ");
            this.changeType = ChangeType.ALL;
            recalculate = true;
            for(SiteModel siteModel:siteModels){
                MCMCNodeFactory.checkDirtiness(siteModel);
            }
        }else if(rateList.somethingIsDirty()){
            changeType = rateList.getChangeType();
            //System.err.println("changeType: "+changeType);

            if(changeType == ChangeType.ADDED){
                addSiteModel();
                setupPointerIndices();
            }else if(changeType == ChangeType.REMOVED){
                removeSiteModel(rateList.getRemovedIndex());
                setupPointerIndices();
            }else if(changeType == ChangeType.SPLIT){
                addSiteModel();
            }else if(changeType == ChangeType.MERGE){                
                removeSiteModel(rateList.getRemovedIndex());
            }else if(changeType == ChangeType.VALUE_CHANGED){
                for(SiteModel siteModel:siteModels){
                    MCMCNodeFactory.checkDirtiness(siteModel);
                }

            }else {
                changeType = ChangeType.ALL;
                for(SiteModel siteModel:siteModels){
                    MCMCNodeFactory.checkDirtiness(siteModel);
                }
                setupPointerIndices();
            }
            recalculate = true;
            //System.err.println("dirty1");

        }else if (ratePointers.somethingIsDirty()){
            //When there is a change in cluster number or a simple pointer change,
            //so the logic has to be written in a way so that the option that there
            //is a change in the number of clusters is eliminated.

            //System.err.println("pointer changed ");
            recalculate = true;
            changeType = ChangeType.POINTER_CHANGED;
            setupPointerIndices();

        }

        return recalculate;
    }




    protected void store(){
        prevCluster = -1;
        super.store();
        //System.out.println("storing");
    }

    public void restore(){
        prevCluster = -1;
        super.restore();
        //System.out.println("restoring");
    }




}