package beast.evolution.sitemodel;

import beast.core.Description;
import beast.core.Input;
import beast.core.MCMCNodeFactory;
import beast.core.parameter.*;
import beast.evolution.alignment.Alignment;
import beast.evolution.substitutionmodel.SubstitutionModel;

import javax.swing.plaf.basic.BasicInternalFrameTitlePane;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Set;

/**
 * @author Chieh-Hsi Wu
 */
@Description("Site model that can handle multiple alignments as well as DPM on rates.")
public class DPMultiAlignRateSiteModel extends DPMultiAlignSiteModel {

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

    public Input<List<SubstitutionModel.Base>> substModelsInput = new Input<List<SubstitutionModel.Base>>(
            "substModel",
            "A list of substitution models.",
            new ArrayList<SubstitutionModel.Base>(),
            Input.Validate.REQUIRED
    );



    public Input<DPValuable> dpValRateInput = new Input<DPValuable>(
            "dpVal",
            "The DPValuable object for rates.",
            Input.Validate.REQUIRED
    );

    protected DPPointer ratePointers;
    protected ParameterList rateList;
    private List<SubstitutionModel.Base> substModels;

    private HashMap<Integer,QuietSiteModel>[]  siteModelMatrix;
    private HashMap<Integer,QuietSiteModel>[] storedSiteModelMatrix;
    private HashMap<Integer, Integer>[] siteModelCount;
    private HashMap<Integer, Integer>[] storedSiteModelCount;
    private DPValuable dpValRate;
    private int[] siteModelMap;
    private int[] storedSiteModelMap;




    public void initAndValidate() throws Exception{
        super.initAndValidate();
        alignments = alignmentsInput.get();
        substModels = substModelsInput.get();
        if(alignments.size() != substModels.size()){
            throw new RuntimeException("The number of assignments ("+alignments.size()+") " +
                    "does not match up with the the number of models ("+substModels.size()+").");
        }

        ratePointers = ratePointersInput.get();
        rateList = rateListInput.get();
         dpValRate =  dpValRateInput.get();

        int rateCount = rateList.getDimension();

        int alignmentCount = alignments.size();


        siteModelMatrix = (HashMap<Integer, QuietSiteModel>[]) new HashMap[alignmentCount];
        storedSiteModelMatrix = (HashMap<Integer, QuietSiteModel>[]) new HashMap[alignmentCount];

        siteModelCount = (HashMap<Integer,Integer>[]) new HashMap[alignmentCount];
        storedSiteModelCount = (HashMap<Integer,Integer>[]) new HashMap[alignmentCount];

        int totalSiteCount = ratePointers.getDimension();
        siteModelMap = new int[totalSiteCount];
        storedSiteModelMap = new int[totalSiteCount];

        alignmentRef = new int[totalSiteCount];
        int k = 0;
        for(int i = 0; i < alignmentCount; i++){
            siteModelMatrix[i] =  new HashMap<Integer,QuietSiteModel>();
            siteModelCount[i] = new HashMap<Integer, Integer>();

            Alignment alignment = alignments.get(i);
            int siteCount = alignment.getSiteCount();
            ArrayList<Integer> rateIndices = new ArrayList<Integer>();

            for(int j = 0; j < siteCount; j++){
                alignmentRef[k] = i;
                //System.out.println("site: "+k+" "+alignmentRef[k]);
                siteModelMap[k] = ratePointers.getParameter(k).getIDNumber();
                int rateIndex = rateList.indexOf(ratePointers.getParameter(k));
                if(!rateIndices.contains(rateIndex)){
                    rateIndices.add(rateIndex);
                }
                if(siteModelCount[i].containsKey(siteModelMap[k])){
                    siteModelCount[i].put(siteModelMap[k], siteModelCount[i].get(siteModelMap[k]) + 1);

                }else{
                    siteModelCount[i].put(siteModelMap[k],1);
                }
                //System.out.println(i+" "+siteModelCount[i].get(siteModelMap[k]));
                k++;

            }


            for(int j = 0; j < rateCount; j++){
                if(rateIndices.contains(j)){
                    QuietRealParameter muParameter = rateList.getParameter(j);
                    QuietSiteModel siteModel = new QuietSiteModel(substModels.get(i),muParameter);
                    siteModelMatrix[i].put(muParameter.getIDNumber(),siteModel);
                    siteModels.add(siteModel);


                }
            }

        }




    }

    public int getAlignmentCount(){
        return siteModelMap.length;
    }

    public int getAlignmentCategoryCount(int index){
        return siteModelCount[index].size();
    }

    public int getLastAddedIndex(){
        throw new RuntimeException("Not applicable!");
    }
    protected void addSiteModel(){
        /*for(int i = 0; i < siteModelCount.length; i++){
                System.out.print("siteModelCount, "+siteModelCount[i].size()+": ");
                Set<Integer> keys = siteModelCount[i].keySet();
                for(int key:keys){
                    System.out.print(key+" ");
                }
                System.out.println();

            }


        System.out.println(siteModels.size());*/
        try{
            //Get the new mutation parameter
            QuietRealParameter muParameter = rateList.getParameter(rateList.getDimension()-1);
            //Get the index of the alignment given the index of the pointer last updated
            int dirtySite = ratePointers.getLastDirty();
            int alignmentIndex = alignmentRef[dirtySite];


            QuietSiteModel siteModel = new QuietSiteModel(substModels.get(alignmentIndex),muParameter);
            int rateID = muParameter.getIDNumber();
            siteModelMatrix[alignmentIndex].put(rateID,siteModel);
            siteModelCount[alignmentIndex].put(rateID,1);
            int oldRateID = siteModelMap[ratePointers.getLastDirty()];
            int oldClusterCount = siteModelCount[alignmentIndex].get(oldRateID)- 1;
            siteModelCount[alignmentIndex].put(oldRateID, oldClusterCount);
            siteModels.add(siteModel);
            //System.out.println(siteModels.size());
            siteModelMap[ratePointers.getLastDirty()] = rateID;
            //System.out.println("siteModelCount: "+siteModelCount[alignmentIndex].size());
            if(siteModelCount[alignmentIndex].get(oldRateID) == 0){
                siteModelCount[alignmentIndex].remove(oldRateID);
                siteModels.remove(siteModelMatrix[alignmentIndex].remove(oldRateID));
            }
            /*System.out.println(siteModels.size());
            for(int i = 0; i < siteModelCount.length; i++){
                System.out.print("siteModelCount, "+siteModelCount[i].size()+": ");
                Set<Integer> keys = siteModelCount[i].keySet();
                for(int key:keys){
                    System.out.print(key+" "+siteModelCount[i].get(key));
                }
                System.out.println();

            } */
        }catch(Exception e){
            throw new RuntimeException(e);
        }

    }

    public int[] getSwappedSites(){
        return ratePointers.getSwappedSites();//todo Check this
    }

    protected void splitSiteModel(){


        try{
            //Get the new mutation parameter
            QuietRealParameter muParameter = rateList.getParameter(rateList.getDimension()-1);
            //Get the indices of the alignment given the index of the pointer last updated
            int[] siteIndices = dpValRate.getClusterSites(rateList.getDimension()-1);
            /*for(int siteIndex: siteIndices){
                System.out.print(siteIndex+" ");
            }
            System.out.println();*/
            dirtySites = siteIndices;

            ArrayList<Integer> alignmentIndices = new ArrayList<Integer>();
            for(int i = 0; i < siteIndices.length; i++){
                //System.out.println("alignmentRef[siteIndices[i]]: "+alignmentRef[siteIndices[i]]);
                if(!alignmentIndices.contains(alignmentRef[siteIndices[i]])){
                    alignmentIndices.add(alignmentRef[siteIndices[i]]);
                }
            }
            //System.out.println("size: "+alignmentIndices.size());

            int rateID = muParameter.getIDNumber();
            for(int i = 0; i < alignmentIndices.size(); i++){
                QuietSiteModel siteModel = new QuietSiteModel(substModels.get(alignmentIndices.get(i)),muParameter);
                siteModels.add(siteModel);
                siteModelMatrix[alignmentIndices.get(i)].put(rateID,siteModel);
                siteModelCount[alignmentIndices.get(i)].put(rateID,0);
                //System.out.println("put rateID: "+rateID);

            }

            int oldRateID = siteModelMap[siteIndices[0]];
            int[] addCount = new int[siteModelCount.length];
            for(int i = 0; i < siteIndices.length; i++){
                addCount[alignmentRef[siteIndices[i]]]++;
                siteModelMap[siteIndices[i]] = rateID;
            }
            for(int i = 0; i < siteModelCount.length;i++){
                if(addCount[i] > 0){
                    siteModelCount[i].put(rateID, addCount[i]);

                    //System.out.println("oldRateID: "+oldRateID+" "+siteModelCount[i].get(oldRateID));
                    int updatedOldCount = siteModelCount[i].get(oldRateID) - addCount[i];
                    siteModelCount[i].put(oldRateID, updatedOldCount);
                    //System.out.println("oldRateID: "+oldRateID+" "+updatedOldCount);
                    if(updatedOldCount == 0){
                        siteModelCount[i].remove(oldRateID);
                        siteModels.remove(siteModelMatrix[i].remove(oldRateID));
                    }
                }

            }

            /*for(int i = 0; i < siteModelCount.length; i++){
                System.out.print("siteModelCount, "+siteModelCount[i].size()+": ");
                Set<Integer> keys = siteModelCount[i].keySet();
                for(int key:keys){
                    System.out.print(key+" "+siteModelCount[i].get(key)+" ");
                }
                System.out.println();

            }*/

        }catch(Exception e){
            throw new RuntimeException(e);
        }

        /*for(int i = 0; i < siteModelCount.length; i++){
                System.out.print("siteModelCount, "+siteModelCount[i].size()+": ");
                Set<Integer> keys = siteModelCount[i].keySet();
                for(int key:keys){
                    System.out.print(key+" "+siteModelCount[i].get(key)+" ");
                }
                System.out.println();

            }
        System.out.println("siteModels2: "+siteModels.size());*/
    }



    int[] dirtySites;

    public int[] getLastDirtySites(){
        return dirtySites;
    }
    public void mergeSiteModel(){


        //Get the indices of the alignment given the index of the pointer last updated
        int[] siteIndices = dpValRate.getStoredClusterSites(rateList.getRemovedIndex());

        //System.out.println();
        dirtySites = siteIndices;
        //Get the id number of the remaining category
        int mergedRateID = ratePointers.getParameterIDNumber(siteIndices[0]);
        int removedRateID = rateList.getRemovedIDNumber();
        //System.out.println("mergedRateID: "+mergedRateID+" removedRateID: "+removedRateID);
        int[] siteCounts = new int[alignments.size()];
        ArrayList<Integer> alignmentIndices = new ArrayList<Integer>();
        for(int i = 0; i < siteIndices.length; i++){
            if(!alignmentIndices.contains(alignmentRef[siteIndices[i]])){
                alignmentIndices.add(alignmentRef[siteIndices[i]]);
            }
            siteModelMap[siteIndices[i]] = mergedRateID;
            siteCounts[alignmentRef[siteIndices[i]]]++;
            //siteModelMap[siteIndices[i]] = ratePointers.getParameterIDNumber(siteIndices[i]); // sanity checking
        }
        //System.out.println("alignmentIndices: "+alignmentIndices.size());

        QuietRealParameter mergedParameter = ratePointers.getParameter(siteIndices[0]);
        try{
            for(int i = 0; i < alignmentIndices.size();i++){

                int alignmentIndex = alignmentIndices.get(i);
                //System.out.println("alignmentIndex: "+alignmentIndex+" mergedRateID: "+mergedRateID);
                if(siteModelCount[alignmentIndex].containsKey(mergedRateID)){
                    //System.out.println("merged0: "+siteModelCount[alignmentIndex].get(mergedRateID));
                    siteModelCount[alignmentIndex].put(
                            mergedRateID,
                            siteModelCount[alignmentIndex].get(mergedRateID) + siteCounts[alignmentIndex]
                    );
                    //System.out.println("merged1: "+siteModelCount[alignmentIndex].get(mergedRateID));
                }else{
                    siteModelCount[alignmentIndex].put(mergedRateID, siteCounts[alignmentIndex]);
                    QuietSiteModel siteModel = new QuietSiteModel(substModels.get(alignmentIndex),mergedParameter);
                    siteModelMatrix[alignmentIndex].put(mergedRateID,siteModel);
                    siteModels.add(siteModel);
                    //System.out.println("merged2: "+siteModelCount[alignmentIndex].get(mergedRateID));

                }
            }
        }catch(Exception e){
            throw new RuntimeException(e);
        }


        //System.out.println(siteModels.size());
        //System.out.println("removedRateID: "+removedRateID);
        int alignmentIndex;
        for(int i = 0; i < alignmentIndices.size(); i++){
            //System.out.println("i: "+i +" as: "+alignmentIndices.size());
            alignmentIndex = alignmentIndices.get(i);
            //System.out.println(alignmentIndex+" "+removedRateID);
            int updatedOldCount = siteModelCount[alignmentIndex].get(removedRateID) - siteCounts[alignmentIndex];
            if(updatedOldCount == 0){
                //System.out.println("HI!");
                siteModelCount[alignmentIndex].remove(removedRateID);
                SiteModel removedSiteModel = siteModelMatrix[alignmentIndex].remove(removedRateID);
                siteModels.remove(removedSiteModel);
                //System.out.println(siteModels.size()+" "+(removedSiteModel == null));
            }else{
               siteModelCount[alignmentIndex].put(removedRateID, updatedOldCount);
            }
        }

        /*for(int i = 0; i < siteModelCount.length; i++){
                System.out.print("siteModelCount, "+siteModelCount[i].size()+": ");
                Set<Integer> keys = siteModelCount[i].keySet();
                for(int key:keys){
                    System.out.print(key+" "+siteModelCount[i].get(key)+" ");
                }
                System.out.println();

            }
        System.out.println("siteModels0: "+siteModels.size());
        for(int i = 0; i < siteModelMatrix.length; i++){
                System.out.print("siteModelMatrix, "+siteModelMatrix[i].size()+": ");
                Set<Integer> keys = siteModelMatrix[i].keySet();
                for(int key:keys){
                    System.out.print(key+" ");
                }
                System.out.println();

            }
        System.out.println("siteModels1: "+siteModels.size());
        System.out.println("End merge");*/
    }

    protected void removeSiteModel(int index){

        /*System.out.println(siteModels.size());
            for(int i = 0; i < siteModelCount.length; i++){
                System.out.print("siteModelCount, "+siteModelCount[i].size()+": ");
                Set<Integer> keys = siteModelCount[i].keySet();
                for(int key:keys){
                    System.out.print(key+" "+siteModelCount[i].get(key));
                }
                System.out.println();

            }*/


        int siteIndex = ratePointers.getLastDirty();

        //System.out.println("removed: "+index+" "+siteModelMap[siteIndex]);
        siteModelMap[siteIndex] = ratePointers.getParameterIDNumber(siteIndex); // update sideModelMap

        int oldRateID = rateList.getRemovedIDNumber();
        siteModelCount[alignmentRef[siteIndex]].remove(oldRateID);
        SiteModel removedSiteModel = siteModelMatrix[alignmentRef[siteIndex]].remove(oldRateID);
        siteModels.remove(removedSiteModel);

        Integer currCount = siteModelCount[alignmentRef[siteIndex]].get(siteModelMap[siteIndex]);
        //System.out.println("currCount: "+currCount+" "+alignmentRef[siteIndex]+" old: "+oldRateID);
        if(currCount == null){
            siteModelCount[alignmentRef[siteIndex]].put(siteModelMap[siteIndex],1);
            try{
                QuietSiteModel siteModel = new QuietSiteModel(substModels.get(alignmentRef[siteIndex]),rateList.getParameter(ratePointers.indexInList(index,rateList)));
                siteModels.add(siteModel);
                siteModelMatrix[alignmentRef[siteIndex]].put(siteModelMap[siteIndex],siteModel);
            }catch(Exception e){
                throw new RuntimeException(e);
            }


        }else{
            siteModelCount[alignmentRef[siteIndex]].put(siteModelMap[siteIndex],currCount+1);
        }




        /*for(int i = 0; i < siteModelMap.length; i++){
            System.out.print(siteModelMap[i]+" ");
        }
        System.out.println();

        for(int i = 0; i < siteModelCount.length; i++)  {
            Set<Integer> keys = siteModelCount[i].keySet();
            System.out.print("key: ");
            for(Integer key: keys){
                System.out.print(key);
            }
            System.out.println();
        } */
    }

    public void switchSiteModel(){
        int index = ratePointers.getLastDirty();  //index of the site changed
        int oldRateID = siteModelMap[index];
        /*System.out.print(index+" ");
        for(int i = 0; i < siteModelMap.length; i++){
            System.out.print(siteModelMap[i]+" ");
        }
        System.out.println();
        System.out.println("size: "+siteModelCount[alignmentRef[index]].size());
        Set<Integer> keys = siteModelCount[alignmentRef[index]].keySet();
        System.out.print("key: ");
        for(Integer key: keys){
            System.out.print(key+" ");
        }
        System.out.println();
        System.out.println(siteModelMap[index]+" "+storedSiteModelMap[index]);*/
        int oldClusterCount = siteModelCount[alignmentRef[index]].get(oldRateID) - 1;


        //remove clusters with no members
        if(oldClusterCount == 0){
            //System.out.println("Removing: "+oldRateID);
            siteModelCount[alignmentRef[index]].remove(oldRateID);
            SiteModel removedSiteModel = siteModelMatrix[alignmentRef[index]].remove(oldRateID);
            siteModels.remove(removedSiteModel);
        }else{
            siteModelCount[alignmentRef[index]].put(oldRateID,oldClusterCount);
        }

        int newRateID = ratePointers.getParameterIDNumber(index);
        siteModelMap[index] = newRateID;
        // If we don't have siteModel with this rate value for this alignment.
        //System.out.println("alignmentRef[index]]: "+alignmentRef[index]+" newRateID: "+newRateID);
        if(!siteModelCount[alignmentRef[index]].containsKey(newRateID)){


            //Create new site count of this rate for the alignment
            siteModelCount[alignmentRef[index]].put(newRateID,1);

            //Create new site model of this rate for the alignment
            try{
                //Get the new rate parameter
                QuietRealParameter muParameter = ratePointers.getParameter(index);
                QuietSiteModel siteModel = new QuietSiteModel(substModels.get(alignmentRef[index]),muParameter);
                siteModelMatrix[alignmentRef[index]].put(newRateID,siteModel);
                siteModels.add(siteModel);

            }catch(Exception e){
                throw new RuntimeException(e);

            }
        }else{//If the appropriate site model already exists

            //Only the site counts of this rate category for this alignment needs to be updated.
            siteModelCount[alignmentRef[index]].put(newRateID,siteModelCount[alignmentRef[index]].get(newRateID)+1);

        }



    }

    public void store(){

        for(int i = 0; i < siteModelMatrix.length;i++){

            storedSiteModelMatrix[i] = (HashMap<Integer,QuietSiteModel>)siteModelMatrix[i].clone();
        }
        for(int i = 0; i < siteModelCount.length; i++){
            storedSiteModelCount[i] = (HashMap<Integer,Integer>)siteModelCount[i].clone();
        }
        System.arraycopy(siteModelMap,0,storedSiteModelMap,0,siteModelMap.length);

        super.store();
    }

    public int getSiteModelWeight(int alignmentIndex, int category){
        return siteModelCount[alignmentIndex].get(category) == null ? 0:siteModelCount[alignmentIndex].get(category);
    }

    public void restore(){
        HashMap<Integer,QuietSiteModel>[] temp1 = siteModelMatrix;
        siteModelMatrix = storedSiteModelMatrix;
        storedSiteModelMatrix = temp1;

        HashMap<Integer, Integer>[] temp2 = siteModelCount;
        siteModelCount = storedSiteModelCount;
        storedSiteModelCount = temp2;

        int[] temp3 = siteModelMap;
        siteModelMap = storedSiteModelMap;
        storedSiteModelMap = temp3;

        super.restore();

    }



    public int getLastDirtySite(){
        return ratePointers.getLastDirty();
    }


    public int getCurrCluster(int index){
        return ratePointers.indexInList(index,rateList);
    }

    public int getPrevCluster(int index){
        return ratePointers.storedIndexInList(index,rateList);
    }

    public int getCurrCategoryIDNumber(int siteIndex){
        return ratePointers.getParameterIDNumber(siteIndex);
    }

    public int getPrevCategoryIDNumber(int siteIndex){
        return ratePointers.getStoredParameterIDNumber(siteIndex);
    }

    public int getCategoryIDNumberFromIndex(int categoryIndex){
        return rateList.getParameter(categoryIndex).getIDNumber();
    }

    public int getAlignmentIndex(int index){
        return alignmentRef[index];
    }

    public QuietSiteModel getSiteModel(int alignmentIndex, int siteModelID){
        return siteModelMatrix[alignmentIndex].get(siteModelID);
    }



    public boolean requiresRecalculation(){

        boolean recalculate = false;
        boolean substitutionModel = false;
        for(SubstitutionModel.Base substModel:substModels){
            if(substModel.isDirtyCalculation()){
                substitutionModel = true;
            }
        }

        if(substitutionModel){

            this.changeType = ChangeType.ALL;
            recalculate = true;

        }else if(rateList.somethingIsDirty()){

            changeType = rateList.getChangeType();
            //System.out.println("changeType: "+changeType);
            if(changeType == ChangeType.ADDED){

                addSiteModel();
                //setupPointerIndices();

            }else if(changeType == ChangeType.REMOVED){

                removeSiteModel(rateList.getRemovedIndex());

                //setupPointerIndices();
            }else if(changeType == ChangeType.SPLIT){

                splitSiteModel();

            }else if(changeType == ChangeType.MERGE){

                mergeSiteModel();

            }else if(changeType == ChangeType.VALUE_CHANGED){
                for(SiteModel siteModel:siteModels){
                    MCMCNodeFactory.checkDirtiness(siteModel);
                }

            }else {
                changeType = ChangeType.ALL;
                for(SiteModel siteModel:siteModels){
                    MCMCNodeFactory.checkDirtiness(siteModel);
                }
                //setupPointerIndices();
            }
            recalculate = true;
            //System.err.println("dirty1");

        }else if (ratePointers.somethingIsDirty()){
            //When there is a change in cluster number or a simple pointer change,
            //so the logic has to be written in a way so that the option that there
            //is a change in the number of clusters is eliminated.

            //System.out.println("pointer changed ");
            recalculate = true;
            changeType = ChangeType.POINTER_CHANGED;
            switchSiteModel();
            //setupPointerIndices();

        }
        /*int counter=0;
        for(int i = 0; i < siteModelCount.length; i++)  {
            Set<Integer> keys = siteModelCount[i].keySet();

            for(Integer key: keys){
                System.out.println("key: "+key+" "+siteModelCount[i].get(key));
                counter+=siteModelCount[i].get(key);
            }

        }
        if(counter != ratePointers.getDimension()){
            throw new RuntimeException("Uh-oh");
        }



        for(int i = 0; i < siteModelMap.length; i++){
            System.out.print(siteModelMap[i]+" ");
        }
        System.out.println(); */

        return recalculate;
    }



    public static void main(String[] args){
        HashMap<Integer, String> map1 = new HashMap<Integer, String>();
        map1.put(0,"WTF");
        map1.put(1,"bloody hell!");
        HashMap<Integer,String> map2 = (HashMap<Integer,String>) map1.clone();
        map1.remove(0);
        System.out.println(map2.get(0));


        HashMap<Integer, Integer> map3 = new HashMap<Integer, Integer>();
        map3.put(0,0);
        map3.put(1,1);
        HashMap<Integer,Integer> map4 = (HashMap<Integer,Integer>) map3.clone();
        map3.put(0,2);
        System.out.println(map3.get(0)+" "+map4.get(0));


    }

    public int getDimension(){
        ArrayList<Integer> categoryCount = new ArrayList<Integer>();
        for(int i = 0; i < siteModelMatrix.length; i++){
            Set<Integer> keys = siteModelMatrix[i].keySet();
            for(int key:keys){
                //System.out.print(key +" ");
                if(!categoryCount.contains(key)){
                    categoryCount.add(key);
                }
            }
            //System.out.println();
        }
        if(rateList.getDimension() != categoryCount.size())
            throw new RuntimeException(categoryCount.size()+" "+rateList.getDimension());
        return categoryCount.size();
        //return rateList.getDimension();
    }




}
