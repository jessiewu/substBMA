package beast.evolution.sitemodel;

import beast.core.Description;
import beast.core.Input;
import beast.core.MCMCNodeFactory;
import beast.core.parameter.ChangeType;
import beast.core.parameter.QuietRealParameter;
import beast.evolution.substitutionmodel.DPNtdBMA;

/**
 * @author Chieh-Hsi Wu
 */
@Description("This site model class facilitates partition selection by DPP on nucleotide substitution models.")

public class DPNtdSiteModel extends DPSingleAlignSiteModel {

    private DPNtdBMA dpNtdBMA;
    private QuietRealParameter muParameter;


    public Input<DPNtdBMA> dpNtdBMAInput = new Input<DPNtdBMA>(
            "ntdBMAList",
            "array which points a set of unique parameter values",
            Input.Validate.REQUIRED
    );


    public Input<QuietRealParameter> muParameterInput = new Input<QuietRealParameter>(
            "mu",
            "subsitution rate for all partitions in an alignment",
            Input.Validate.REQUIRED
    );


    int counter = 0;
    public void initAndValidate() throws Exception{
        super.initAndValidate();
        dpNtdBMA = dpNtdBMAInput.get();
        muParameter = muParameterInput.get();
        int ntdBMACount = dpNtdBMA.getDimension();



        for(int i = 0;i < ntdBMACount; i++){
            QuietSiteModel siteModel = new QuietSiteModel(dpNtdBMA.getModel(i), muParameter);
            siteModels.add(siteModel);
            siteModel.setID(""+counter++);
        }

    }

    public int[] getSubstModelPointerIndices(){
        return dpNtdBMA.getPointerIndices();
    }

    public int getRemovedIndex(){
        return dpNtdBMA.getRemovedIndex();
    }

    public int getLastAddedIndex(){
        return dpNtdBMA.getLastAddedIndex();
    }

    public int getDirtySiteModelIndex(){
        return dpNtdBMA.getDirtyModelIndex();
    }





    public int getLastDirtySite(){
        return dpNtdBMA.getLastDirtySite();
    }

    public int[] getLastDirtySites(){
        return dpNtdBMA.getLastDirtySites();
    }


    public int getSiteModelCount(){
        return siteModels.size();
    }



    protected void addSiteModel(){
        try{
            QuietSiteModel siteModel = new QuietSiteModel(dpNtdBMA.getModel(dpNtdBMA.getDimension()-1), muParameter);
            siteModels.add(siteModel);
            siteModel.setID(""+counter++);
        }catch(Exception e){
            throw new RuntimeException(e);
        }

    }

    public int[] getSwappedSites(){
        return dpNtdBMA.getSwappedSites();

    }




    


    public int getPrevCluster(int index){
        return dpNtdBMA.getPrevCluster(index);

    }

    public int getCurrCluster(int index){
        return dpNtdBMA.getCurrCluster(index);

    }

    public int getPrevCategoryIDNumber(int siteIndex){
        return dpNtdBMA.getPrevModelIDNumber(siteIndex);
    }

    public int getCurrCategoryIDNumber(int siteIndex){
            return dpNtdBMA.getCurrModelIDNumber(siteIndex);

        }


    public boolean requiresRecalculation(){

        boolean recalculate = false;
        if(dpNtdBMA.isDirtyCalculation()){
            changeType = dpNtdBMA.getChangeType();
            //System.err.println("dpNtd: "+changeType);
            if(changeType == ChangeType.ADDED || changeType == ChangeType.SPLIT){
                addSiteModel();

            }else if(changeType == ChangeType.REMOVED || changeType == ChangeType.MERGE){
                removeSiteModel(dpNtdBMA.getRemovedIndex());

            }else if(changeType == ChangeType.VALUE_CHANGED){

                for(SiteModel siteModel:siteModels){
                    MCMCNodeFactory.checkDirtiness(siteModel);
                }


            }else if(changeType == ChangeType.POINTER_CHANGED){
                this.changeType = ChangeType.POINTER_CHANGED;

            }else{
                this.changeType = ChangeType.ALL;
                for(SiteModel siteModel:siteModels){
                    MCMCNodeFactory.checkDirtiness(siteModel);

                }
            }
            recalculate = true;

        }else if(muParameterInput.get().somethingIsDirty()){
            this.changeType = ChangeType.ALL;
            for(SiteModel siteModel:siteModels){
                MCMCNodeFactory.checkDirtiness(siteModel);
            }
            recalculate = true;

        }

        return recalculate;
    }






}
