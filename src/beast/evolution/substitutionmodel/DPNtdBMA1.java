package beast.evolution.substitutionmodel;

import beast.core.*;
import beast.core.parameter.*;

import java.util.ArrayList;

/**
 * @author Chieh-Hsi Wu
 */
@Description("A list of SwitchingNtdBMA classes for DPP.")
public class DPNtdBMA1 extends DPNtdBMA implements PluginList, Recycle {




    public Input<ParameterList> pointersRefListInput = new Input<ParameterList>(
                "pointersRefList",
                "A list of unique model values",
                Input.Validate.REQUIRED
        );









    private ParameterList pointerRefList;

    public void initAndValidate(){
        paramList = paramListInput.get();
        modelList = modelListInput.get();
        freqsList = freqsListInput.get();
        pointerRefList = pointersRefListInput.get();
        pointers = pointersInput.get();

        int dimParamList = paramList.getDimension();
        for(int i = 0; i < dimParamList;i++){
            ntdBMAs.add(
                    createSwitchingNtdBMA(
                            paramList.getParameter(i),
                            modelList.getParameter(i),
                            freqsList.getParameter(i)
                    )
            );

        }
        //System.err.println("model counts: "+ntdBMAs.size());
        pointerIndices = new int[pointers.getDimension()];
        resetAllPointerIndices();
    }


    public int getLastAddedIndex(){
        //System.out.println("getLastAddedIndex() "+pointerRefList.getLastAddedIndex());
        return pointerRefList.getLastAddedIndex();
    }



    public int[] getPointerIndices(){
        return pointerIndices;
    }







    int prevCluster = -1;
    public void setupPointerIndices(){
        //System.out.println("changeType: "+changeType);
        if(changeType == ChangeType.ADDED || changeType == ChangeType.POINTER_CHANGED){
            int changedIndex = pointers.getLastDirty();
            prevCluster = pointerIndices[changedIndex];
            pointerIndices[changedIndex] = pointers.indexInList(changedIndex,pointerRefList);
        }else if (changeType == changeType.POINTERS_SWAPPED){
            int[] changedSites = pointers.getSwappedSites();
            pointerIndices[changedSites[0]] = pointers.indexInList(changedSites[0],pointerRefList);
            pointerIndices[changedSites[1]] = pointers.indexInList(changedSites[1],pointerRefList);
        }else if (changeType == ChangeType.REMOVED){
            resetAllPointerIndices();
        }else if (changeType == ChangeType.ALL){
            resetAllPointerIndices();
        }
    }



    public int getPrevCluster(int index){
        return pointers.storedIndexInList(index,pointerRefList);
    }

    public int getCurrCluster(int index){
        return pointers.indexInList(index,pointerRefList);
    }


 

    public void resetAllPointerIndices(){
        for(int i = 0; i < pointerIndices.length;i++){
            pointerIndices[i] = pointers.indexInList(i,pointerRefList);
        }
    }







    public boolean requiresRecalculation(){

        boolean recalculate = false;
        //System.err.println("dirty0");
        if(paramList.somethingIsDirty()||modelList.somethingIsDirty()||freqsList.somethingIsDirty()){
            //System.out.println("dirty0: "+freqsList.getChangeType());
            //ChangeType changeType;
            if(paramList.somethingIsDirty()){
                changeType = paramList.getChangeType();

            }else if(modelList.somethingIsDirty()){
                changeType = modelList.getChangeType();

            }else{
                changeType = freqsList.getChangeType();
            }

            //System.out.println("dirty0: "+changeType);
            if(changeType == ChangeType.ADDED){
                //System.out.println(getID()+"model added");
                addModel();
                //setupPointerIndices();

            }else if(changeType == ChangeType.REMOVED){
                //System.out.println(getID()+"model removed, "+paramList.getRemovedIndex());
                removeModel(freqsList.getRemovedIndex());
                //setupPointerIndices();

            }else if(changeType == ChangeType.VALUE_CHANGED){
                //System.out.println(getID()+": model changed");
                for(SwitchingNtdBMA ntdBMA:ntdBMAs){

                    MCMCNodeFactory.checkDirtiness(ntdBMA);
                    //System.out.println(ntdBMA.getIDNumber()+" is dirty? "+ntdBMA.isDirtyCalculation());
                }
            }else if(changeType == ChangeType.SPLIT){
                addModel();
            }else if(changeType == ChangeType.SPLIT_AND_VALUE_CHANGE){

                addModel();
                MCMCNodeFactory.checkDirtiness(ntdBMAs.get(freqsList.getDirtyIndex()));
            }else if(changeType == ChangeType.MERGE){
                removeModel(freqsList.getRemovedIndex());
            }else if(changeType == ChangeType.MERGE_AND_VALUE_CHANGE){

                removeModel(freqsList.getRemovedIndex());
                MCMCNodeFactory.checkDirtiness(ntdBMAs.get(freqsList.getDirtyIndex()));

            }else{
                this.changeType = ChangeType.ALL;
                setupPointerIndices();
                for(SwitchingNtdBMA ntdBMA:ntdBMAs){
                    MCMCNodeFactory.checkDirtiness(ntdBMA);
                    //System.out.println(ntdBMA.getIDNumber()+" is dirty? "+ntdBMA.isDirtyCalculation());
                }
            }
            recalculate = true;
            //System.err.println("dirty1");

        }else if (pointers.somethingIsDirty()){
            recalculate = true;
            //System.out.println(getID()+": pointer changed");
            changeType = pointers.getChangeType();

            //setupPointerIndices();
        }

        //System.out.println("recalculate: "+recalculate+" "+changeType);

        return recalculate;
    }

    protected void accept() {
        for(SwitchingNtdBMA ntdBMA:ntdBMAs){
            ntdBMA.makeAccept();
        }
    	super.accept();
    }




}
