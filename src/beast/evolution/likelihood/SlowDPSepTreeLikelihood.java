package beast.evolution.likelihood;

import beast.core.Description;
import beast.core.MCMCNodeFactory;
import beast.evolution.sitemodel.DPNtdRateSepSiteModel;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.SwitchingNtdBMA;
import beast.core.parameter.ChangeType;
import beast.core.Input;

import java.util.ArrayList;

/**
 * @author Chieh-Hsi Wu
 */
@Description("Does a lot of unnecessary calculation (calculates likelihoods at sites that are not required). Used for testing.")
public class SlowDPSepTreeLikelihood extends SlowDPTreeLikelihood{

    private DPNtdRateSepSiteModel dpSiteModel;
    private WVTreeLikelihood[][] treeLiksMatrix;
    private WVTreeLikelihood[][] storedTreeLiksMatrix;
    private ChangeType changeType = ChangeType.ALL;



    public SlowDPSepTreeLikelihood(){
        /*
         * We need DPNtdRateSepSiteModel specifically,
         * therefore an Input<DPSiteModel> is not useful.
         */
        dpValInput.setRule(Input.Validate.OPTIONAL);
    }

    public void initAndValidate() throws Exception{


        alignment = m_data.get();
        int patternCount = alignment.getPatternCount();
        if(!(m_pSiteModel.get() instanceof DPNtdRateSepSiteModel)){
            throw new RuntimeException("DPNtdRateSepSiteModel object is required for site model");
        }
        dpSiteModel = (DPNtdRateSepSiteModel) m_pSiteModel.get();
        int siteModelCount = dpSiteModel.getSiteModelCount();

        /*
         * Set up a 3 dimensional array to store the weights of each ntdBMA/rate combination (dimensions: ntdBMA, rates, weights).
         */
        int[][][]clusterWeights = new int[dpSiteModel.getSubstClusterLimit()][dpSiteModel.getRatesClusterLimit()][patternCount];
        int[] clusterIds;
        int siteCount = alignment.getSiteCount();
        for(int i = 0; i < siteCount; i++){
            clusterIds = dpSiteModel.getCurrClusters(i);
            clusterWeights[clusterIds[DPNtdRateSepSiteModel.NTDBMA]][clusterIds[DPNtdRateSepSiteModel.RATES]][alignment.getPatternIndex(i)]++;
        }

        /*
         *  Traverse through the list site models and create tree likelihoods.
         *  The tree likelihoods are then added to a list and the matrix.
         *  The order of likelihoods in the list corresponds to that of site models in the DPSiteModel list.
         */
        //treeLiksMatrix = new WVTreeLikelihood[dpSiteModel.getSubstClusterLimit()][dpSiteModel.getRatesClusterLimit()];
        treeLiksMatrix = new WVTreeLikelihood[dpSiteModel.getSubstClusterLimit()][dpSiteModel.getRatesClusterLimit()];
        storedTreeLiksMatrix = new WVTreeLikelihood[dpSiteModel.getSubstClusterLimit()][dpSiteModel.getRatesClusterLimit()];
        for(int i = 0; i < siteModelCount;i++){

            //Get the ids and hence the positions of siteModel/treeLikelihoods in a list
            int ntdBMAId = ((SwitchingNtdBMA)dpSiteModel.getSiteModel(i).getSubstitutionModel()).getIDNumber();
            int ratesId = dpSiteModel.getSiteModel(i).getRateParameter().getIDNumber();

            //Create the tree likelihood
            //WVTreeLikelihood treeLik = new WVTreeLikelihood(clusterWeights[ntdBMAId][ratesId]);
            WVTreeLikelihood treeLik = new WVTreeLikelihood(clusterWeights[ntdBMAId][ratesId]);
            treeLik.initByName(
                    "data", alignment,
                    "tree", m_tree.get(),
                    "siteModel", dpSiteModel.getSiteModel(i),
                    "branchRateModel", m_pBranchRateModel.get(),
                    "useAmbiguities",useAmbiguitiesInput.get()
            );

            //Add to list and matrix for the convenience of processesing
            treeLiks.add(treeLik);
            treeLiksMatrix[ntdBMAId][ratesId] = treeLik;
            storedTreeLiksMatrix[ntdBMAId][ratesId] = treeLik;
        }
    }

    /*public double calculateLogP() throws Exception{
        int sum = getSumWeight();
        if(sum != alignment.getSiteCount()){
            throw new RuntimeException("Weights ("+sum+") and site count("+alignment.getSiteCount()+").");
        }
        return super.calculateLogP();
    }*/

    private void update(){
        int dirtySite = dpSiteModel.getLastDirtySite();
        int[] currClusterIds = dpSiteModel.getCurrClusters(dirtySite);
        int[] prevClusterIds = dpSiteModel.getPrevClusters(dirtySite);

        //Create a new tree likelihood (with zero weights) when there is a new site model.
        if(treeLiksMatrix[currClusterIds[DPNtdRateSepSiteModel.NTDBMA]][currClusterIds[DPNtdRateSepSiteModel.RATES]] == null){
            //System.out.println("Add clusters at: "+currClusterIds[DPNtdRateSepSiteModel.NTDBMA]+" "+currClusterIds[DPNtdRateSepSiteModel.RATES]);
            addTreeLikelihood(currClusterIds[DPNtdRateSepSiteModel.NTDBMA],currClusterIds[DPNtdRateSepSiteModel.RATES]);
        }
        //System.out.println(prevClusterIds[DPNtdRateSepSiteModel.NTDBMA]+" "+prevClusterIds[DPNtdRateSepSiteModel.RATES]);
        //Move weight
        moveWeight(
            prevClusterIds[DPNtdRateSepSiteModel.NTDBMA],
            prevClusterIds[DPNtdRateSepSiteModel.RATES],
            currClusterIds[DPNtdRateSepSiteModel.NTDBMA],
            currClusterIds[DPNtdRateSepSiteModel.RATES],
            dirtySite,
            1
        );

        //Remove likelihoods that have zero weights
        if(dpSiteModel.getCombinationWeight(prevClusterIds[DPNtdRateSepSiteModel.NTDBMA], prevClusterIds[DPNtdRateSepSiteModel.RATES]) == 0){
            treeLiks.remove(treeLiksMatrix[prevClusterIds[DPNtdRateSepSiteModel.NTDBMA]][prevClusterIds[DPNtdRateSepSiteModel.RATES]]);
            treeLiksMatrix[prevClusterIds[DPNtdRateSepSiteModel.NTDBMA]][prevClusterIds[DPNtdRateSepSiteModel.RATES]] = null;
        }
    }

    private void updates(){
        int[] dirtySites =  dpSiteModel.getLastDirtySites();
        for(int dirtySite:dirtySites){
            update(dirtySite);
        }

        for(int dirtySite:dirtySites){
        //int dirtySite = dpSiteModel.getLastDirtySite();
            int[] prevClusterIds = dpSiteModel.getPrevClusters(dirtySite);

        //Remove likelihoods that have zero weights
            if(dpSiteModel.getCombinationWeight(prevClusterIds[DPNtdRateSepSiteModel.NTDBMA], prevClusterIds[DPNtdRateSepSiteModel.RATES]) == 0){
                treeLiks.remove(treeLiksMatrix[prevClusterIds[DPNtdRateSepSiteModel.NTDBMA]][prevClusterIds[DPNtdRateSepSiteModel.RATES]]);
                treeLiksMatrix[prevClusterIds[DPNtdRateSepSiteModel.NTDBMA]][prevClusterIds[DPNtdRateSepSiteModel.RATES]] = null;
            }
        }
    }

    private void update(int dirtySite){
        //int dirtySite = dpSiteModel.getLastDirtySite();
        int[] currClusterIds = dpSiteModel.getCurrClusters(dirtySite);
        int[] prevClusterIds = dpSiteModel.getPrevClusters(dirtySite);

        //Create a new tree likelihood (with zero weights) when there is a new site model.
        if(treeLiksMatrix[currClusterIds[DPNtdRateSepSiteModel.NTDBMA]][currClusterIds[DPNtdRateSepSiteModel.RATES]] == null){
            //System.out.println("Add cluster at: "+currClusterIds[DPNtdRateSepSiteModel.NTDBMA]+" "+currClusterIds[DPNtdRateSepSiteModel.RATES]);
            addTreeLikelihood(currClusterIds[DPNtdRateSepSiteModel.NTDBMA],currClusterIds[DPNtdRateSepSiteModel.RATES]);
        }
        //System.out.println(prevClusterIds[DPNtdRateSepSiteModel.NTDBMA]+" "+prevClusterIds[DPNtdRateSepSiteModel.RATES]);
        //Move weight


        moveWeight(
            prevClusterIds[DPNtdRateSepSiteModel.NTDBMA],
            prevClusterIds[DPNtdRateSepSiteModel.RATES],
            currClusterIds[DPNtdRateSepSiteModel.NTDBMA],
            currClusterIds[DPNtdRateSepSiteModel.RATES],
            dirtySite,
            1
        );

        /*    for(int i = 0; i < treeLiksMatrix.length;i++){
                for(int j = 0; j < treeLiksMatrix[i].length; j++){
                    System.out.print(i +" "+j+": ");
                    if(treeLiksMatrix[i][j] == null){
                        System.out.println("null");

                    }else{
                    int[] weights = treeLiksMatrix[i][j].getPatternWeights();
                    for(int k = 0; k < weights.length;k++){
                        System.out.print(weights[k]+" ");
                    }
                    System.out.println();
                    }
                }
            }

         */


    }


    /*
     * Moving weight from one cluster combination to another
     */
    public void moveWeight(
            int ntdBMAIndexFrom,
            int rateIndexFrom,
            int ntdBMAIndexTo,
            int rateIndexTo,
            int dirtySite,
            int weight){
        /*for(int i = 0; i < treeLiksMatrix.length;i++){
                for(int j = 0; j < treeLiksMatrix[i].length; j++){
                    System.out.print(i +" "+j+": ");
                    if(treeLiksMatrix[i][j] == null){
                        System.out.println("null");

                    }else{
                    int[] weights = treeLiksMatrix[i][j].getPatternWeights();
                    for(int k = 0; k < weights.length;k++){
                        System.out.print(weights[k]+" ");
                    }
                    System.out.println();
                    }
                }
            }
        System.out.println(ntdBMAIndexFrom+" "+rateIndexFrom+" "+(treeLiksMatrix[ntdBMAIndexFrom][rateIndexFrom] ==null));*/
        //try{
        treeLiksMatrix[ntdBMAIndexFrom][rateIndexFrom].removeWeight(alignment.getPatternIndex(dirtySite),weight);
        treeLiksMatrix[ntdBMAIndexTo][rateIndexTo].addWeight(alignment.getPatternIndex(dirtySite),weight);
        //}catch(Exception e){



        //}



    }


    public void addTreeLikelihood(int substModelID, int rateID){
        //SiteModel siteModel = dpSiteModel.getLastAdded();
        SiteModel siteModel = dpSiteModel.getSiteModel(substModelID, rateID);


        int[] patternWeights = new int[alignment.getPatternCount()];

        //WVTreeLikelihood treeLik = new WVTreeLikelihood(patternWeights);
        WVTreeLikelihood treeLik = new WVTreeLikelihood(patternWeights);
        try{
            treeLik.initByName(
                    "data", alignment,
                    "tree", m_tree.get(),
                    "siteModel", siteModel,
                    "branchRateModel", m_pBranchRateModel.get(),
                    "useAmbiguities",useAmbiguitiesInput.get()
            );
            treeLik.calculateLogP();
            treeLik.store();
            treeLiks.add(treeLik);
            treeLiksMatrix[substModelID][rateID] =treeLik;
        }catch(Exception e){
            throw new RuntimeException(e);
        }

    }

    public double getSiteLogLikelihood(int iCluster, int iSite){
        throw new RuntimeException("Retrieving the site given a site index and cluster index is not applicable in this case.");
    }

    /*
     * @param inputType either ntdBMA or rates
     * @clusterID id of the cluster
     * @siteIndex the index of a site in an alignment
     */
    public double getSiteLogLikelihood(int inputType, int clusterID, int siteIndex){
        //System.out.println("pattern: "+alignment.getPatternIndex(siteIndex));
        int[] currClusters = dpSiteModel.getCurrClusters(siteIndex);
        //System.out.println("clusterID: "+clusterID+" "+prevClusters[DPNtdRateSepSiteModel.RATES]+" "+alignment.getPatternIndex(siteIndex));
        if(inputType == DPNtdRateSepSiteModel.NTDBMA){
            if(treeLiksMatrix[clusterID][currClusters[DPNtdRateSepSiteModel.RATES]] != null){

                //WVTreeLikelihood tmpTL = treeLiksMatrix[clusterID][prevClusters[DPNtdRateSepSiteModel.RATES]];
                WVTreeLikelihood tmpTL = treeLiksMatrix[clusterID][currClusters[DPNtdRateSepSiteModel.RATES]];
                //System.out.println("clusterID: "+clusterID+" "+prevClusters[DPNtdRateSepSiteModel.RATES]+" "+alignment.getPatternIndex(siteIndex));
                //tmpTL.printThings();
                return tmpTL.getPatternLogLikelihood(alignment.getPatternIndex(siteIndex));

            }
        }else{

            if(treeLiksMatrix[currClusters[DPNtdRateSepSiteModel.NTDBMA]][clusterID] != null){
                //WVTreeLikelihood tmpTL = treeLiksMatrix[prevClusters[DPNtdRateSepSiteModel.NTDBMA]][clusterID];
                //System.out.println("hi!!");
                WVTreeLikelihood tmpTL = treeLiksMatrix[currClusters[DPNtdRateSepSiteModel.NTDBMA]][clusterID];
                return tmpTL.getPatternLogLikelihood(alignment.getPatternIndex(siteIndex));
            }
        }

        return Double.NaN;

    }




    public void printThings(){
        //for(WVTreeLikelihood treeLik:treeLiks){
        for(WVTreeLikelihood treeLik:treeLiks){
            treeLik.printThings();
            for(int i = 0; i < alignment.getPatternCount();i++){
                System.out.print(treeLik.getPatternLogLikelihood(alignment.getPatternIndex(i))+" ");
            }
            System.out.println();
        }
        System.out.print(0+" "+0+" ");
        treeLiksMatrix[0][0].printThings();
        System.out.print(0+" "+1+" ");
        treeLiksMatrix[0][1].printThings();
        System.out.print(0+" "+2+" ");
        treeLiksMatrix[0][2].printThings();
        System.out.print(0+" "+3+" ");
        treeLiksMatrix[0][3].printThings();
    }


    public int getSumWeight(){
        int sum = 0;
        //for(WVTreeLikelihood treeLik: treeLiks){
        for(WVTreeLikelihood treeLik: treeLiks){
            sum +=treeLik.weightSum();
        }
        return sum;
    }


    boolean storeTreeLikelihoods = false;
    @Override
    protected boolean requiresRecalculation() {
        boolean recalculate = false;
        if(dpSiteModel.isDirtyCalculation()){

            changeType = dpSiteModel.getChangeType();
            //System.out.println("treeLik requires recal!!"+changeType);
            if(changeType == ChangeType.ADDED || changeType == ChangeType.REMOVED || changeType == ChangeType.POINTER_CHANGED){
                //System.out.println("changeType: "+changeType);
                update();
            }else if(changeType == ChangeType.SPLIT || changeType == ChangeType.MERGE){
                //storeTreeLikelihoods();

                updates();
            }else if(changeType == ChangeType.VALUE_CHANGED){

                changeType = ChangeType.VALUE_CHANGED;

            }else{
                changeType = ChangeType.ALL;
            }
            recalculate = true;
        }else if(m_tree.get().somethingIsDirty()){
            recalculate = true;

        }else if(m_pBranchRateModel.get().isDirtyCalculation()){
            recalculate = true;
        }
        
        if(recalculate){
            int counter = 0;
            for(TreeLikelihood treeLik:treeLiks){
                MCMCNodeFactory.checkDirtiness(treeLik);
            }
        }
        return recalculate;
    }
    public void store(){
        storeTreeLikelihoods = true;
        for(int i = 0; i < treeLiksMatrix.length;i++){
            System.arraycopy(treeLiksMatrix[i],0,storedTreeLiksMatrix[i],0, treeLiksMatrix[i].length);
        }
        super.store();

    }

    /*public void storeTreeLikelihoods(){
        storeTreeLikelihoods = true;
        for(int i = 0; i < treeLiksMatrix.length;i++){
            System.arraycopy(treeLiksMatrix[i],0,storedTreeLiksMatrix[i],0, treeLiksMatrix[i].length);
        }
        //super.store();
    }*/

    public void restore(){
        //if(storeTreeLikelihoods){
        WVTreeLikelihood[][] tmp = treeLiksMatrix;
        treeLiksMatrix = storedTreeLiksMatrix;
        storedTreeLiksMatrix = tmp;
            storeTreeLikelihoods = false;
        //}
        super.restore();
    }



}
