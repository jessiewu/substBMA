package beast.evolution.sitemodel;

import beast.core.Description;
import beast.core.MCMCNodeFactory;
import beast.core.parameter.*;
import beast.core.Input;
import beast.evolution.substitutionmodel.DPNtdBMA;

import java.util.ArrayList;

/**
 * @author Chieh-Hsi Wu
 */
@Description("Manages a list of simple gamma site models. Can create and remove site models on the fly. Used when a DPP is applied to the partitioning scheme of the substitution model and rate together.")
public class DPNtdRateSiteModel extends DPRateSiteModel{
    protected DPNtdBMA dpNtdBMA;
    public Input<DPNtdBMA> dpNtdBMAInput = new Input<DPNtdBMA>(
            "ntdBMAList",
            "array which points a set of unique parameter values",
            Input.Validate.REQUIRED
    );

    public DPNtdRateSiteModel(){
        substModelInput.setRule(Input.Validate.OPTIONAL);
    }

    public void initAndValidate() throws Exception{
        //super.initAndValidate();
        dpNtdBMA = dpNtdBMAInput.get();

        int ntdBMACount = dpNtdBMA.getDimension();

        ratePointers = ratePointersInput.get();
        rateList = rateListInput.get();

        if(rateList.getDimension() != dpNtdBMA.getDimension()){
            throw new RuntimeException("The number of clustres for rates and substution models must be the same.");
        }

        //Setting up site models
        siteModels = new ArrayList<QuietSiteModel>();

        for(int i = 0;i < ntdBMACount; i++){
            QuietRealParameter muParameter = rateList.getParameter(ratePointers.indexInList(i,rateList));
            QuietSiteModel siteModel = new QuietSiteModel(dpNtdBMA.getModel(i),muParameter);
            siteModels.add(siteModel);
        }
    }



    protected void addSiteModel(){
        try{

            QuietRealParameter muParameter = rateList.getParameter(rateList.getDimension()-1);

            QuietSiteModel siteModel = new QuietSiteModel(dpNtdBMA.getModel(dpNtdBMA.getDimension()-1),muParameter);
            siteModels.add(siteModel);
        }catch(Exception e){
            throw new RuntimeException(e);
        }

    }

    public int[] getPointerIndices(){
        return dpNtdBMA.getPointerIndices();
    }

    public int[] getSwappedSites(){
        return ratePointers.getSwappedSites();
    }




    public boolean requiresRecalculation(){

        boolean recalculate = false;
        //System.err.println("dirty0");
        ChangeType substModelChangeType = dpNtdBMA.getChangeType();
        if(rateList.somethingIsDirty() && dpNtdBMA.isDirtyCalculation()){
            changeType = rateList.getChangeType();
            if(changeType != substModelChangeType){
                throw new RuntimeException("Can only handle same type of changes to subst and rate at one time.");
            }


            if(changeType == ChangeType.ADDED||changeType == ChangeType.SPLIT){
                addSiteModel();

                //setupPointerIndices();
            }else if(changeType == ChangeType.REMOVED || changeType == ChangeType.MERGE){
                removeSiteModel(rateList.getRemovedIndex());
                //setupPointerIndices();
            }else if(changeType == ChangeType.VALUE_CHANGED){

                for(SiteModel siteModel:siteModels){
                    MCMCNodeFactory.checkDirtiness(siteModel);
                }

            }else {
                this.changeType = ChangeType.ALL;
                for(SiteModel siteModel:siteModels){
                    MCMCNodeFactory.checkDirtiness(siteModel);
                }
                //setupPointerIndices();
            }
            recalculate = true;
            //System.err.println("dirty1");

        }else if (ratePointers.somethingIsDirty()){

            changeType = ratePointers.getChangeType();
            if(changeType != substModelChangeType){
                System.out.println(changeType+" "+substModelChangeType);
                throw new RuntimeException("Can only handle same type of changes to subst and rate at one time.");
            }
            recalculate = true;

            //this.changeType = ChangeType.POINTER_CHANGED;


            //setupPointerIndices();

        }else if(rateList.somethingIsDirty()){
            this.changeType = rateList.getChangeType();
            for(SiteModel siteModel:siteModels){
                MCMCNodeFactory.checkDirtiness(siteModel);
            }
            recalculate = true;
        }else if(dpNtdBMA.isDirtyCalculation() && dpNtdBMA.getChangeType() == ChangeType.VALUE_CHANGED){
            this.changeType = ChangeType.VALUE_CHANGED;
            for(SiteModel siteModel:siteModels){
                MCMCNodeFactory.checkDirtiness(siteModel);
            }
            recalculate = true;
        }

        return recalculate;
    }








}
