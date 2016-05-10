package beast.evolution.sitemodel;

import beast.core.Description;
import beast.core.MCMCNodeFactory;
import beast.core.parameter.*;
import beast.core.Input;
import beast.evolution.substitutionmodel.DPNtdBMA;
import beast.evolution.substitutionmodel.SwitchingNtdBMA;
import beast.evolution.substitutionmodel.SubstitutionModel;

import java.util.Arrays;

/**
 * @author Chieh-Hsi Wu
 */
@Description("Manages a list of simple gamma site models. Can create and remove site models on the fly. Used when independent DPPs are applied to the partitioning schemes of the substitution model and rate.")
public class DPNtdRateSepSiteModel extends DPSingleAlignSiteModel {

    protected int lastDirtySite = -1;
    public static final int NTDBMA = 0;
    public static final int RATES = 1;
    protected DPNtdBMA dpNtdBMA;
    protected DPPointer ratesPointers;
    protected ParameterList ratesList;
    public Input<DPNtdBMA> dpNtdBMAInput = new Input<DPNtdBMA>(
            "ntdBMAList",
            "array which points a set of unique parameter values",
            Input.Validate.REQUIRED
    );
    //ParameterList
    public Input<ParameterList> ratesListInput = new Input<ParameterList>(
            "ratesList",
            "A list of unique site rate values",
            Input.Validate.REQUIRED
    );


    //assignment
    public Input<DPPointer> ratesPointersInput = new Input<DPPointer>(
            "ratesPointers",
            "array which points a set of unique parameter values",
            Input.Validate.REQUIRED
    );

    public Input<Integer> substClusterLimitInput = new Input<Integer>(
            "substClusterLimit",
            "Initial limit placed on the dimension of siteModel matrix for substitution models",
            50
    );
    public Input<Integer> ratesClusterLimitInput = new Input<Integer>(
            "ratesClusterLimit",
            "Initial limit placed on the dimension of siteModel matrix for rates",
            50
    );

    public Input<DPValuable> dpValSubstInput = new Input<DPValuable>(
        "dpValSubst",
        "Dirichlet process valuable object that records the info of clustering"
    );

    public Input<DPValuable> dpValRateInput = new Input<DPValuable>(
        "dpValRate",
        "Dirichlet process valuable object that records the info of clustering"
    );

    protected int eltCount;
    protected int substClusterLimit;
    protected int ratesClusterLimit;
    protected QuietSiteModel[][] siteModelsMatrix;
    protected QuietSiteModel[][] storedSiteModelsMatrix;
    protected int[][] siteModelWeights;
    protected int[][] storedSiteModelWeights;
    //private int[] substWeights;
    //private int[] ratesWeights;
    protected int[][] clusterMap;
    protected int[][] storedClusterMap;
    protected DPValuable dpValSubst;
    protected DPValuable dpValRate;



    public void initAndValidate() {
        super.initAndValidate();
        if(siteModels == null){
            throw new RuntimeException("Site model list is not instantiated.");
        }
        dpNtdBMA = dpNtdBMAInput.get();
        ratesPointers = ratesPointersInput.get();
        ratesList = ratesListInput.get();
        eltCount = ratesPointers.getDimension();
        if(dpNtdBMA.getCount() != eltCount){
            throw new RuntimeException("The number of model elements is different to the number of rate elements.");
        }
        dpValSubst = dpValSubstInput.get();
        dpValRate = dpValRateInput.get();

        substClusterLimit = substClusterLimitInput.get();
        ratesClusterLimit = ratesClusterLimitInput.get();
        setup();
    }

    /*
     * Setup the site model and weight matrices.
     */
    public void setup(){
        siteModelsMatrix = new QuietSiteModel[substClusterLimit][ratesClusterLimit];
        storedSiteModelsMatrix = new QuietSiteModel[substClusterLimit][ratesClusterLimit];
        siteModelWeights = new int[substClusterLimit][ratesClusterLimit];
        storedSiteModelWeights = new int[substClusterLimit][ratesClusterLimit];
        //substWeights = new int[substClusterLimit];
        //ratesWeights = new int[ratesClusterLimit];
        clusterMap = new int[2][ratesPointers.getDimension()];
        storedClusterMap = new int[2][ratesPointers.getDimension()];

        int[] substModelPointerIndicies = dpNtdBMA.getPointerIndices();
        int substModelIndex = -1;
        int rateIndex = -1;

        for(int i = 0; i < eltCount; i++){
            substModelIndex = substModelPointerIndicies[i];
            //System.out.println("eltCount: "+substModelIndex);
            rateIndex = ratesPointers.indexInList(i,ratesList);
            if(siteModelsMatrix[substModelIndex][rateIndex] == null){

                try{

                    QuietSiteModel siteModel = new QuietSiteModel(dpNtdBMA.getModel(substModelIndex),ratesList.getParameter(rateIndex));
                    siteModelsMatrix[substModelIndex][rateIndex] = siteModel;
                    storedSiteModelsMatrix[substModelIndex][rateIndex] = siteModel;
                    siteModels.add(siteModel);

                }catch(Exception e){
                    throw new RuntimeException(e);
                }

            }
            siteModelWeights[substModelIndex][rateIndex]++;
            storedSiteModelWeights[substModelIndex][rateIndex]++;
            //substWeights[substModelIndex]++;
            //ratesWeights[rateIndex]++;
            clusterMap[NTDBMA][i] = substModelIndex;
            clusterMap[RATES][i] = rateIndex;
            storedClusterMap[NTDBMA][i] = substModelIndex;
            storedClusterMap[RATES][i] = rateIndex;
        }

    }

    public int getLastAddedIndex(){
        throw new RuntimeException("Not applicable!");
    }


    public int getCombinationWeight(int substID, int rateID){
        /*System.out.println("getCombinationWeight");
        for(int i = 0; i < siteModelWeights.length;i++){
            for(int j = 0; j < siteModelWeights[i].length;j++){
                System.out.print(siteModelWeights[i][j]+" ");
            }
            System.out.println();
        }*/
        return siteModelWeights[substID][rateID];
    }
    public int getSubstClusterLimit(){
        return substClusterLimit;
    }

    public int getRatesClusterLimit(){
        return ratesClusterLimit;
    }

    public SubstitutionModel getModel(int siteIndex){
        //System.out.println(siteModelsMatrix[clusterMap[NTDBMA][siteIndex]][clusterMap[RATES][siteIndex]]==null);
        return siteModelsMatrix[clusterMap[NTDBMA][siteIndex]][clusterMap[RATES][siteIndex]].getSubstitutionModel();
    }

    public QuietRealParameter getRate(int unitIndex){
        return siteModelsMatrix[clusterMap[NTDBMA][unitIndex]][clusterMap[RATES][unitIndex]].getRateParameter();
    }

    int[] clusterSites;

    protected void handleSplit(int changedInput) throws Exception{

        int ntdBMAIDNum;
        int muIDNum;
        if(changedInput == NTDBMA){
            //Retreive the new NtdBMA model
            ntdBMAIDNum = dpNtdBMA.getLastModel().getIDNumber();
            //Retrieve all the (dirty) sites that uses the newly added model
            clusterSites = dpValSubst.getClusterSites(dpNtdBMA.getDimension()-1);

            for(int dirtySite: clusterSites){
                //Get the rate of the each of the dirty sites
                muIDNum = ratesPointers.getParameter(dirtySite).getIDNumber();
                //Update the substModel and rate combination
                updateMap(dirtySite, ntdBMAIDNum, muIDNum);
            }

        }else if(changedInput == RATES){

            muIDNum = ratesList.getLastParameter().getIDNumber();
            clusterSites = dpValRate.getClusterSites(ratesList.getDimension()-1);
            for(int dirtySite:clusterSites){
                ntdBMAIDNum = dpNtdBMA.getModel(dpNtdBMA.getCurrCluster(dirtySite)).getIDNumber();
                updateMap(dirtySite, ntdBMAIDNum, muIDNum);
            }
        }else{
            throw new RuntimeException("Can only add clusters for either ntdBMA model or rates, changedInput: "+changedInput);
        }

    }


    /*
     * Creates a new site model when there is a new combination of substitution model and mutation rate
     */
    public void addSiteModel(int changedInput, int lastDirtySite) throws Exception{
        int ntdBMACluster;
        int rateCluster;
        if(changedInput == NTDBMA){
            ntdBMACluster = dpNtdBMA.getLastModel().getIDNumber();
            rateCluster = ratesPointers.getParameter(lastDirtySite).getIDNumber();

        }else if(changedInput == RATES){
            ntdBMACluster = dpNtdBMA.getModel(dpNtdBMA.getCurrCluster(lastDirtySite)).getIDNumber();
            rateCluster = ratesList.getLastParameter().getIDNumber();
        }else{
            throw new RuntimeException("Can only add clusters for either ntdBMA model or rates, changedInput: "+changedInput);
        }

        updateMap(lastDirtySite, ntdBMACluster, rateCluster);

    }

    void handleMerge(int changedInput) throws Exception{

        int ntdBMAIDNum;
        int muIDNum;
        if(changedInput == NTDBMA){
            int removedSubstModelIndex = dpNtdBMA.getRemovedIndex();
            //System.out.println(removedSubstModelIndex);
            clusterSites = dpValSubst.getStoredClusterSites(removedSubstModelIndex);
            //This the substModel cluster that takes up the members of the removed substModel cluster
            ntdBMAIDNum = dpNtdBMA.getDirtyModelIDNumber();
            //System.out.println("ntdBMAIDNum: "+dpNtdBMA.getDimension());
            for(int dirtySite: clusterSites){
                muIDNum = ratesPointers.getParameterIDNumber(dirtySite);
                updateMap(dirtySite, ntdBMAIDNum, muIDNum);
            }

        }else if(changedInput == RATES){
            //The rate cluster that takes up the member of the removed rate cluster
            muIDNum = ratesList.getDirtyParameterIDNumber();

            int removedRateIndex = ratesList.getRemovedIndex();
            clusterSites = dpValRate.getStoredClusterSites(removedRateIndex);

            for(int dirtySite:clusterSites){
                ntdBMAIDNum = dpNtdBMA.getModel(dpNtdBMA.getCurrCluster(dirtySite)).getIDNumber();
                //System.out.println(ntdBMAIDNum +" " +muIDNum);
                updateMap(dirtySite, ntdBMAIDNum, muIDNum);

            }
        }else{
            throw new RuntimeException("Can only add clusters for either ntdBMA model or rates, changedInput: "+changedInput);
        }

    }



    public void removeSiteModel(int changedInput, int lastDirtySite) throws Exception {
        //int[] siteClusterMap = new int[2];
        int ntdBMACluster;
        int rateCluster;
        if(changedInput == NTDBMA){

            ntdBMACluster = dpNtdBMA.getModel(dpNtdBMA.getCurrCluster(lastDirtySite)).getIDNumber();
            rateCluster = ratesPointers.getParameter(lastDirtySite).getIDNumber();

        }else if(changedInput == RATES){

            ntdBMACluster = dpNtdBMA.getModel(dpNtdBMA.getCurrCluster(lastDirtySite)).getIDNumber();
            rateCluster = ratesPointers.getParameterIDNumber(lastDirtySite);

        }else{
            throw new RuntimeException("Can only remove clusters for either ntdBMA model or rates.");
        }


        //System.out.println("changedInput: "+changedInput+" "+ntdBMAId+" "+rateId);
        //siteModels.remove(siteModels.indexOf(siteModelsMatrix[ntdBMAId][rateId]));
        //siteModelsMatrix[ntdBMAId][rateId] = null;

        updateMap(lastDirtySite, ntdBMACluster, rateCluster);

    }



    void handleMultiPointerChanges(int changedInput) throws Exception{

        if(changedInput == NTDBMA){
            clusterSites = dpValSubst.getLastDirtySites();
        }else if(changedInput == RATES){
            clusterSites = dpValRate.getLastDirtySites();
        }else{
            throw new RuntimeException("Can only add clusters for either ntdBMA model or rates, changedInput: "+changedInput);
        }

        for(int dirtySite: clusterSites){
            SwitchingNtdBMA ntdBMA = dpNtdBMA.getModel(dpNtdBMA.getCurrCluster(dirtySite));
            QuietRealParameter muParameter = ratesPointers.getParameter(dirtySite);

            int ntdBMACluster = ntdBMA.getIDNumber();
            int rateCluster = muParameter.getIDNumber();

            //update mapping and weights
            updateWeights(dirtySite, ntdBMACluster, rateCluster);



            //updateModelMatrix(dirtySite);



        }
        for(int dirtySite: clusterSites){
            updateModelMatrix(dirtySite);

        }

    }

    public int getDirtySiteModelIndex(){
        throw new RuntimeException("Getting a single dirty site model index is not always applicable in this case!");
    }

    public int getDirtySubstModelIDNumber() {
        return dpNtdBMA.getDirtyModelIDNumber();
    }

    public int getDirtySiteModelIDNumber(){
        return ratesList.getDirtyParameterIDNumber();
    }

    public void checkSiteModelsDirtiness(int changedInput){
        int fixedIndex = -1;

        if(changedInput == NTDBMA){

            fixedIndex = dpNtdBMA.getModel(dpNtdBMA.getDirtyModelIndex()).getIDNumber(); //todo use getDirtyModelIDNumber
            for(int i = 0; i < siteModelsMatrix[fixedIndex].length;i++){
                if(siteModelsMatrix[fixedIndex][i] != null){
                    MCMCNodeFactory.checkDirtiness(siteModelsMatrix[fixedIndex][i]);
                    //siteModelsMatrix[fixedIndex][i].checkDirtiness();
                    //System.out.println(i+" "+fixedIndex+" "+siteModelsMatrix[fixedIndex][i].isDirtyCalculation()+" "+siteModelsMatrix[fixedIndex][i].getRateParameter().getIDNumber());
                }
            }
            //System.out.println(getID()+" fixedIndex: "+fixedIndex);
        }else if(changedInput == RATES){

            fixedIndex = ratesList.getParameter(ratesList.getDirtyIndex()).getIDNumber(); //todo use getDirtyParameterIDNumber
            //System.out.println("fixedIndex: "+fixedIndex);
            for(int i = 0; i < siteModelsMatrix.length;i++){
                if(siteModelsMatrix[i][fixedIndex] != null){
                    //System.out.println(i);
                    MCMCNodeFactory.checkDirtiness(siteModelsMatrix[i][fixedIndex]);
                    //System.out.println(i+" "+fixedIndex+" "+siteModelsMatrix[i][fixedIndex].isDirtyCalculation()+" "+siteModelsMatrix[i][fixedIndex].getRateParameter().getIDNumber());
                    //System.out.println(siteModelsMatrix[i][fixedIndex]+" "+siteModelsMatrix[i][fixedIndex].isDirtyCalculation());
                }
            }

        }else{
            throw new RuntimeException("Can only remove clusters for either ntdBMA model or rates.");
        }
        /*for(int i = 0; i < siteModels.size();i++){
            siteModels.get(i).checkDirtiness();
            if(siteModels.get(i).isDirtyCalculation()){
                SwitchingNtdBMA temp =  ((SwitchingNtdBMA)siteModels.get(i).getSubstitutionModel());
                System.out.println(getID()+" " +temp.getModelChooseIDNumber());

            }
        }*/

    }

    /*
     * Handles changes in pointers that does not involve changes in the number of clusters
     */
    protected void handlePointerChange(int lastDirtySite) throws Exception{

        //Previous cluster ids of the last dirty site
        int prevNtdBMAIdNum = clusterMap[NTDBMA][lastDirtySite];
        int prevRateIdNum = clusterMap[RATES][lastDirtySite];

        //Current cluster ids of the last dirty site
        SwitchingNtdBMA ntdBMA = dpNtdBMA.getModel(dpNtdBMA.getCurrCluster(lastDirtySite)); //todo
        QuietRealParameter muParameter = ratesPointers.getParameter(lastDirtySite);
        //int[] siteClusterMap = new int[2];
        int ntdBMACluster = ntdBMA.getIDNumber();
        int rateCluster = muParameter.getIDNumber();



        //update mapping and weights
        updateMap(lastDirtySite, ntdBMACluster, rateCluster);



        //If the propose combination is new then create a new site model
        if(siteModelsMatrix[ntdBMACluster][rateCluster] == null){
            QuietSiteModel siteModel = new QuietSiteModel(ntdBMA,muParameter);
            siteModelsMatrix[ntdBMACluster][rateCluster] = siteModel;
            siteModels.add(siteModel);
        }

        //If the previous combination no longer has any weight then remove
        //System.out.println(prevNtdBMAIdNum+" "+prevRateIdNum);

    }

    public void handlePointerSwap(int siteIndex1, int siteIndex2) throws Exception{



        SwitchingNtdBMA ntdBMA1 = dpNtdBMA.getModel(dpNtdBMA.getCurrCluster(siteIndex1));
        QuietRealParameter muParameter1 = ratesPointers.getParameter(siteIndex1);

        int ntdBMACluster1 = ntdBMA1.getIDNumber();
        int rateCluster1 = muParameter1.getIDNumber();

        //update mapping and weights
        updateWeights(siteIndex1, ntdBMACluster1, rateCluster1);

        SwitchingNtdBMA ntdBMA2 = dpNtdBMA.getModel(dpNtdBMA.getCurrCluster(siteIndex2));
        QuietRealParameter muParameter2 = ratesPointers.getParameter(siteIndex2);

        int ntdBMACluster2 = ntdBMA2.getIDNumber();
        int rateCluster2 = muParameter2.getIDNumber();

        //update mapping and weights
        updateWeights(siteIndex2, ntdBMACluster2, rateCluster2);

        updateModelMatrix(siteIndex1);
        updateModelMatrix(siteIndex2);



        //If the propose combination is new then create a new site model
        /*if(siteModelsMatrix[ntdBMACluster][rateCluster] == null){
            SiteModel siteModel = new SiteModel();
            siteModel.initByName(
                    "substModel", ntdBMA,
                    "mutationRate", muParameter
            );
            siteModelsMatrix[ntdBMACluster][rateCluster] = siteModel;
            siteModels.add(siteModel);
        }*/



    }



    boolean mappingChanged = false;

    protected void updateWeights(int siteIndex, int ntdBMACluster, int rateCluster){
        mappingChanged = true;
        //Update the weight matrix
        moveWeight(
                clusterMap[NTDBMA][siteIndex],  //current substModel
                clusterMap[RATES][siteIndex],   //current rate
                ntdBMACluster,                  //potentially new substModel
                rateCluster,                    //potentially new rate
                1
        );
        int prevNtdBMAId = clusterMap[NTDBMA][siteIndex];
        int prevRateId = clusterMap[RATES][siteIndex];

        //Update mapping
        clusterMap[NTDBMA][siteIndex] = ntdBMACluster;
        clusterMap[RATES][siteIndex] = rateCluster;
    }

    protected void updateModelMatrix(int siteIndex) throws Exception{
        //Remove site model if it has zero pattern weight
        if(siteModelWeights[storedClusterMap[NTDBMA][siteIndex]][storedClusterMap[RATES][siteIndex]] == 0
                && siteModelsMatrix[storedClusterMap[NTDBMA][siteIndex]][storedClusterMap[RATES][siteIndex]] != null){
            siteModels.remove(siteModels.indexOf(siteModelsMatrix[storedClusterMap[NTDBMA][siteIndex]][storedClusterMap[RATES][siteIndex]]));
            siteModelsMatrix[storedClusterMap[NTDBMA][siteIndex]][storedClusterMap[RATES][siteIndex]] = null;
        }

        //Add site model if this is a new substModel and rate combination
        if(siteModelsMatrix[clusterMap[NTDBMA][siteIndex]][clusterMap[RATES][siteIndex]] == null){
            //System.out.println("add: "+clusterMap[NTDBMA][siteIndex]+" "+clusterMap[RATES][siteIndex]);
            QuietSiteModel siteModel = new QuietSiteModel(dpNtdBMA.getModel(dpNtdBMA.getCurrCluster(siteIndex)),ratesPointers.getParameter(siteIndex));
            siteModelsMatrix[clusterMap[NTDBMA][siteIndex]][clusterMap[RATES][siteIndex]] = siteModel;
            siteModels.add(siteModel);
        }

    }

    /*
     * Update cluster mapping
     */
    public void updateMap(int siteIndex, int ntdBMACluster, int rateCluster) throws Exception{
        /*for(int i = 0; i < clusterMap.length;i++){
            System.out.print("NTDBMA: ");
            for(int j = 0; j < clusterMap[i].length;j++){
                System.out.print(clusterMap[i][j]+" ");
            }
            System.out.println();
        }

        for(int i = 0; i < storedClusterMap.length;i++){
            System.out.print("stored NTDBMA: ");
            for(int j = 0; j < storedClusterMap[i].length;j++){
                System.out.print(storedClusterMap[i][j]+" ");
            }
            System.out.println();
        }*/

        mappingChanged = true;

        //Update the weight matrix
        moveWeight(
                clusterMap[NTDBMA][siteIndex],  //current substModel
                clusterMap[RATES][siteIndex],   //current rate
                ntdBMACluster,                  //potentially new substModel
                rateCluster,                    //potentially new rate
                1
        );

        int prevNtdBMAId = clusterMap[NTDBMA][siteIndex];
        int prevRateId = clusterMap[RATES][siteIndex];

        //Update mapping
        clusterMap[NTDBMA][siteIndex] = ntdBMACluster;
        clusterMap[RATES][siteIndex] = rateCluster;

        //Stored map
        //storedClusterMap[NTDBMA][siteIndex] = prevNtdBMAId;
        //storedClusterMap[RATES][siteIndex] = prevRateId;

        //System.out.println("flag1: "+siteModelWeights[prevNtdBMAId][prevRateId]);


        //Remove site model if it has zero pattern weight
        if(siteModelWeights[prevNtdBMAId][prevRateId] == 0){
            //System.out.println("remove");
            //System.out.println("siteModels: "+siteModels.size());
            //System.out.println(siteModelsMatrix[prevNtdBMAId][prevRateId]);
            siteModels.remove(siteModels.indexOf(siteModelsMatrix[prevNtdBMAId][prevRateId]));
            siteModelsMatrix[prevNtdBMAId][prevRateId] = null;
        }

        //Add site model if this is a new substModel and rate combination
        if(siteModelsMatrix[clusterMap[NTDBMA][siteIndex]][clusterMap[RATES][siteIndex]] == null){
            //System.out.println("add: "+clusterMap[NTDBMA][siteIndex]+" "+clusterMap[RATES][siteIndex]);
            QuietSiteModel siteModel = new QuietSiteModel(dpNtdBMA.getModel(dpNtdBMA.getCurrCluster(siteIndex)),ratesPointers.getParameter(siteIndex));
            siteModelsMatrix[clusterMap[NTDBMA][siteIndex]][clusterMap[RATES][siteIndex]] = siteModel;
            siteModels.add(siteModel);
        }

        /*for(int i = 0; i < clusterMap.length;i++){
            System.out.print("RATE: ");
            for(int j = 0; j < clusterMap[i].length;j++){
                System.out.print(clusterMap[i][j]+" ");
            }
            System.out.println();
        }
        for(int i = 0; i < storedClusterMap.length;i++){
            System.out.print("stored RATE: ");
            for(int j = 0; j < storedClusterMap[i].length;j++){
                System.out.print(storedClusterMap[i][j]+" ");
            }
            System.out.println();
        }*/




    }





    /*
     * Moving weight from one cluster combination to another
     */
    public void moveWeight(
            int ntdBMAIndexFrom,
            int rateIndexFrom,
            int ntdBMAIndexTo,
            int rateIndexTo,
            int weight){
        /*for(int i = 0; i < 3;i++){
            for(int j = 0; j < 3;j++){
                System.out.print(siteModelWeights[i][j]+" ");
            }
            System.out.println();
        }

        for(int i = 0; i < 3;i++){
            for(int j = 0; j < 3;j++){
                System.out.print(storedSiteModelWeights[i][j]+" ");
            }
            System.out.println();
        }


        if(!storeWeightsMatrix[ntdBMAIndexFrom][rateIndexFrom]){
            storedSiteModelWeights[ntdBMAIndexFrom][rateIndexFrom] = siteModelWeights[ntdBMAIndexFrom][rateIndexFrom];
            storedSiteModelWeights[ntdBMAIndexTo][rateIndexTo] = siteModelWeights[ntdBMAIndexTo][rateIndexTo];
            storeWeightsMatrix[ntdBMAIndexFrom][rateIndexFrom] = true;
            storeWeightsMatrix[ntdBMAIndexTo][rateIndexTo] = true;
            storeWeights = true;
        }*/

        //Change site model weights
        siteModelWeights[ntdBMAIndexFrom][rateIndexFrom] -= weight;
        siteModelWeights[ntdBMAIndexTo][rateIndexTo] += weight;


        /*if(siteModelWeights[ntdBMAIndexFrom][rateIndexFrom] < 0){
            throw new RuntimeException("NEGATIVE WEIGHT");

        }*/

        //Change marginal substitution model weights
        //substWeights[ntdBMAIndexFrom] -= weight;
        //substWeights[ntdBMAIndexTo] += weight;

        //Change marginal rate weights
        //ratesWeights[rateIndexFrom] -= weight;
        //ratesWeights[rateIndexTo] += weight;

    }

    public QuietSiteModel getSiteModelOfSiteIndex(int siteIndex){
        return siteModelsMatrix[clusterMap[NTDBMA][siteIndex]][clusterMap[RATES][siteIndex]];
    }

    public int getPrevCluster(int siteIndex){
        throw new RuntimeException("Retrieving the previous cluster index of a given site is not suitable in this case.");
    }

    public int getCurrCluster(int siteIndex){
        throw new RuntimeException("Retrieving the current cluster index of a given site is not suitable in this case.");
    }
    public int getPrevCategoryIDNumber(int siteIndex){
        throw new RuntimeException("Retrieving the previous cluster index of a given site is not suitable in this case.");
    }

    public int getCurrCategoryIDNumber(int siteIndex){
        throw new RuntimeException("Retrieving the current cluster index of a given site is not suitable in this case.");
    }

    public int getRemovedIndex(){
         throw new RuntimeException("Retrieving the removed cluster index of a given site is not suitable in this case.");
    }

    public int[] getCurrClusters(int siteIndex){
        return new int[]{clusterMap[NTDBMA][siteIndex],clusterMap[RATES][siteIndex]};
    }

    public int[] getPrevClusters(int siteIndex){
        int[] prevClusters = new int[2];
        //System.out.println("prev: "+dpNtdBMA.getPrevCluster(siteIndex));
        //prevClusters[NTDBMA] = dpNtdBMA.getPrevModelIDNumber(siteIndex);
        //prevClusters[RATES] = ratesPointers.getStoredParameterIDNumber(siteIndex);
        //System.out.println("prev: "+dpNtdBMA.getPrevCluster(siteIndex));

        return new int[]{storedClusterMap[NTDBMA][siteIndex],storedClusterMap[RATES][siteIndex]};
    }

    public void addSiteModel(){}

    public int getLastDirtySite(){
        return lastDirtySite;
    }

    public QuietSiteModel getSiteModel(int ntdBMAID, int rateID){
        return siteModelsMatrix[ntdBMAID][rateID];
    }

    public QuietSiteModel getStoredSiteModel(int ntdBMAID, int rateID){
        return storedSiteModelsMatrix[ntdBMAID][rateID];
    }

    public int getChangedInput(){
        return changedInput;
    }


    public int[] getLastDirtySites(){
        return Arrays.copyOf(clusterSites, clusterSites.length);
    }

    public int[] getSwappedSites(){
        return getLastDirtySites();
    }


    protected int changedInput = -1;
    public boolean requiresRecalculation(){

        boolean recalculate = false;
        mappingChanged = false;
        try{

            //ChangeType substModelChangeType = dpNtdBMA.getChangeType();
            if(ratesList.somethingIsDirty() ||ratesPointers.somethingIsDirty()|| dpNtdBMA.isDirtyCalculation()){
                //ChangeType inputChangeType = null;


                //Check whether subst model or rates have changed and
                //retrieve change type.
                if(dpNtdBMA.isDirtyCalculation()){
                    changeType = dpNtdBMA.getChangeType();
                    changedInput = NTDBMA;
                    lastDirtySite = dpNtdBMA.getLastDirtySite();
                }else{
                    if(ratesList.somethingIsDirty()){
                        changeType = ratesList.getChangeType();
                    }else{
                        changeType = ratesPointers.getChangeType();
                    }
                    lastDirtySite = ratesPointers.getLastDirty();
                    changedInput = RATES;
                }


                //System.out.println("changeType: "+this.changeType);
                if(changeType == ChangeType.ADDED){
                    //System.out.println("ADD");
                    addSiteModel(changedInput, lastDirtySite);

                }else if(changeType == ChangeType.SPLIT){

                    handleSplit(changedInput);

                }else if(changeType == ChangeType.REMOVED){
                    //System.out.println("REMOVED");
                    removeSiteModel(changedInput, lastDirtySite);

                }else if(changeType == ChangeType.MERGE){

                    handleMerge(changedInput);

                }else if(changeType == ChangeType.POINTER_CHANGED){
                    //System.out.println("POINTER_CHANGED");

                    handlePointerChange(lastDirtySite);
                }else if(changeType == ChangeType.MULTIPLE_POINTERS_CHANGED){

                    handleMultiPointerChanges(changedInput);


                }else if(changeType == ChangeType.POINTERS_SWAPPED){
                    if(changedInput == NTDBMA){
                        clusterSites = dpNtdBMA.getSwappedSites();

                    }else{
                        clusterSites = ratesPointers.getSwappedSites();
                    }

                    if(clusterSites[0] != clusterSites[1]){
                        handlePointerSwap(clusterSites[0], clusterSites[1]);
                        //handlePointerChange(clusterSites[0]);
                        //handlePointerChange(clusterSites[1]);

                    }

                }else if(changeType == ChangeType.VALUE_CHANGED){
                    //System.out.println("VALUE_CHANGED");
                    //When there's a value change, the cluster assignment doesn't change.
                    checkSiteModelsDirtiness(changedInput);

                }else {
                    this.changeType = ChangeType.ALL;
                    for(int i = 0; i < siteModelsMatrix.length; i++){
                        for(int j = 0; j < siteModelsMatrix[i].length; j++){
                            if(siteModelsMatrix[i][j] != null){
                                MCMCNodeFactory.checkDirtiness(siteModelsMatrix[i][j]);

                            }

                        }
                    }
                    //setupPointerIndices();
                }
                recalculate = true;
                //System.err.println("dirty1");

            }



        }catch(Exception e){
            throw new RuntimeException(e);
        }
        //System.out.println("siteModel changeType: "+changeType);
        return recalculate;
    }



    public void store(){
        /*System.out.println("WEIGHTS");
        for(int i = 0; i < siteModelWeights.length;i++){
            for(int j = 0; j < siteModelWeights[i].length;j++){
                System.out.print(siteModelWeights[i][j]+" ");
            }
            System.out.println();
        }*/
        if(!(this instanceof EfficientDPNtdRateSepSiteModel)){
            changedInput = -1;
            for(int i = 0; i < storedClusterMap.length;i++){
                System.arraycopy(clusterMap[i],0,storedClusterMap[i],0,clusterMap[i].length);
            }
            for(int i = 0; i < storedSiteModelWeights.length;i++){
                System.arraycopy(siteModelWeights[i],0,storedSiteModelWeights[i],0,siteModelWeights[i].length);
            }

            for(int i = 0; i < storedSiteModelsMatrix.length;i++){
                System.arraycopy(siteModelsMatrix[i],0,storedSiteModelsMatrix[i],0,siteModelsMatrix[i].length);
            }
        }
        //for(int i = 0; i < clusterMap.length; i++){
        //    System.arraycopy(clusterMap[i],0,storedClusterMap[i],0,clusterMap[i].length);
        //}
        super.store();
    }

    public void restore(){
        if(!(this instanceof EfficientDPNtdRateSepSiteModel)){
            if(mappingChanged){
                int[][] tmp = clusterMap;
                clusterMap = storedClusterMap;
                storedClusterMap = tmp;
            }
            int[][] tmp1 = siteModelWeights;
            siteModelWeights = storedSiteModelWeights;
            storedSiteModelWeights = tmp1;

            QuietSiteModel[][] tmp2 = siteModelsMatrix;
            siteModelsMatrix = storedSiteModelsMatrix;
            storedSiteModelsMatrix = tmp2;
        }
        
        super.restore();
    }





    



}
