package beast.evolution.likelihood;

import beast.app.BeastMCMC;
import beast.core.Description;
import beast.core.MCMCNodeFactory;
import beast.evolution.sitemodel.DPNtdRateSepSiteModel;
import beast.evolution.sitemodel.QuietSiteModel;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.SwitchingNtdBMA;
import beast.core.parameter.ChangeType;
import beast.core.Input;
import beast.evolution.tree.Tree;

import javax.sound.midi.SysexMessage;

/**
 * @author Chieh-Hsi Wu
 */
@Description("Tree likelihood that supports independent Dirichlet process priors on the partitionings of substitution model and rate.")
public class DPSepTreeLikelihood extends DPTreeLikelihood{

    protected DPNtdRateSepSiteModel dpSiteModel;
    protected NewWVTreeLikelihood[][] treeLiksMatrix;
    protected NewWVTreeLikelihood[][] storedTreeLiksMatrix;
    protected ChangeType changeType = ChangeType.ALL;



    public DPSepTreeLikelihood(){
        /*
         * We need DPNtdRateSepSiteModel specifically,
         * therefore an Input<DPSiteModel> is not useful.
         */

        dpValInput.setRule(Input.Validate.OPTIONAL);
    }

    public void initAndValidate() {
        useThreads = useThreadsInput.get() && (BeastMCMC.m_nThreads > 1);
        useThreadsEvenly = useThreadsEvenlyInput.get() && (BeastMCMC.m_nThreads > 1);

        alignment = dataInput.get();
        int patternCount = alignment.getPatternCount();
        if(!(siteModelInput.get() instanceof DPNtdRateSepSiteModel)){
            throw new RuntimeException("DPNtdRateSepSiteModel required for site model.");
        }
        dpSiteModel = (DPNtdRateSepSiteModel)siteModelInput.get();
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
        treeLiksMatrix = new NewWVTreeLikelihood[dpSiteModel.getSubstClusterLimit()][dpSiteModel.getRatesClusterLimit()];
        storedTreeLiksMatrix = new NewWVTreeLikelihood[dpSiteModel.getSubstClusterLimit()][dpSiteModel.getRatesClusterLimit()];
        for(int i = 0; i < siteModelCount;i++){

            //Get the ids and hence the positions of siteModel/treeLikelihoods in a list
            int ntdBMAId = ((SwitchingNtdBMA)dpSiteModel.getSiteModel(i).getSubstitutionModel()).getIDNumber();
            int ratesId = dpSiteModel.getSiteModel(i).getRateParameter().getIDNumber();

            //Create the tree likelihood
            //WVTreeLikelihood treeLik = new WVTreeLikelihood(clusterWeights[ntdBMAId][ratesId]);

            NewWVTreeLikelihood treeLik = new NewWVTreeLikelihood(
                    clusterWeights[ntdBMAId][ratesId],
                    alignment,
                    (Tree) treeInput.get(),
                    useAmbiguitiesInput.get(),
                    dpSiteModel.getSiteModel(i),
                    branchRateModelInput.get());

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

    protected void update(){
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

    public void checkLikelihoodDirtiness(int changedInput){
int fixedIndex = -1;

        if(changedInput == DPNtdRateSepSiteModel.NTDBMA){

            fixedIndex = dpSiteModel.getDirtySubstModelIDNumber(); //todo use getDirtyModelIDNumber
            for(int i = 0; i < treeLiksMatrix[fixedIndex].length;i++){
                if(treeLiksMatrix[fixedIndex][i] != null){
                    MCMCNodeFactory.checkDirtiness(treeLiksMatrix[fixedIndex][i]);
                    //siteModelsMatrix[fixedIndex][i].checkDirtiness();
                    //System.out.println(i+" "+fixedIndex+" "+siteModelsMatrix[fixedIndex][i].isDirtyCalculation()+" "+siteModelsMatrix[fixedIndex][i].getRateParameter().getIDNumber());
                }
            }
            //System.out.println(getID()+" fixedIndex: "+fixedIndex);
        }else if(changedInput == DPNtdRateSepSiteModel.RATES){

            fixedIndex = dpSiteModel.getDirtySiteModelIDNumber(); //todo use getDirtyParameterIDNumber
            //System.out.println("fixedIndex: "+fixedIndex);
            for(int i = 0; i < treeLiksMatrix.length;i++){
                if(treeLiksMatrix[i][fixedIndex] != null){
                    //System.out.println(i);
                    MCMCNodeFactory.checkDirtiness(treeLiksMatrix[i][fixedIndex]);
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

    protected void updates(){
        int[] dirtySites =  dpSiteModel.getLastDirtySites();
        for(int dirtySite:dirtySites){
            update(dirtySite);            
        }

        for(int dirtySite:dirtySites){
        //int dirtySite = dpSiteModel.getLastDirtySite();
            int[] prevClusterIds = dpSiteModel.getPrevClusters(dirtySite);

        //Remove likelihoods that have zero weights
            if(dpSiteModel.getCombinationWeight(prevClusterIds[DPNtdRateSepSiteModel.NTDBMA], prevClusterIds[DPNtdRateSepSiteModel.RATES]) == 0){
                //System.out.println("Remove?");
                treeLiks.remove(treeLiksMatrix[prevClusterIds[DPNtdRateSepSiteModel.NTDBMA]][prevClusterIds[DPNtdRateSepSiteModel.RATES]]);
                treeLiksMatrix[prevClusterIds[DPNtdRateSepSiteModel.NTDBMA]][prevClusterIds[DPNtdRateSepSiteModel.RATES]] = null;
            }
        }
    }









    protected void update(int dirtySite){
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
    protected void moveWeight(
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
        //NewWVTreeLikelihood treeLik = new NewWVTreeLikelihood(patternWeights);
        NewWVTreeLikelihood treeLik = new NewWVTreeLikelihood(
                    patternWeights,
                    alignment,
                    (Tree) treeInput.get(),
                    useAmbiguitiesInput.get(),
                    siteModel,
                    branchRateModelInput.get());
        try{



            treeLik.calculateLogP();
            treeLik.store();
            treeLiks.add(treeLik);
            treeLiksMatrix[substModelID][rateID] = treeLik;
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
                NewWVTreeLikelihood tmpTL = treeLiksMatrix[clusterID][currClusters[DPNtdRateSepSiteModel.RATES]];
                //System.out.println("clusterID: "+clusterID+" "+prevClusters[DPNtdRateSepSiteModel.RATES]+" "+alignment.getPatternIndex(siteIndex));
                //tmpTL.printThings();
                return tmpTL.getPatternLogLikelihood(alignment.getPatternIndex(siteIndex));

            }
        }else{

            if(treeLiksMatrix[currClusters[DPNtdRateSepSiteModel.NTDBMA]][clusterID] != null){
                //WVTreeLikelihood tmpTL = treeLiksMatrix[prevClusters[DPNtdRateSepSiteModel.NTDBMA]][clusterID];
                //System.out.println("hi!!");
                NewWVTreeLikelihood tmpTL = treeLiksMatrix[currClusters[DPNtdRateSepSiteModel.NTDBMA]][clusterID];
                return tmpTL.getPatternLogLikelihood(alignment.getPatternIndex(siteIndex));
            }
        }

        return Double.NaN;

    }




    public void printThings(){
        //for(WVTreeLikelihood treeLik:treeLiks){
        for(NewWVTreeLikelihood treeLik:treeLiks){
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
        for(NewWVTreeLikelihood treeLik: treeLiks){
            sum +=treeLik.weightSum();
        }
        return sum;
    }


    //boolean storeTreeLikelihoods = false;

    @Override
    protected boolean requiresRecalculation() {

         /*QuietSiteModel q = dpSiteModel.getSiteModel(0,2);
        if(q != null && treeLiksMatrix[0][2] !=null){

              if(treeLiksMatrix[0][2].getSiteModel() != q){
                  System.out.println("whoa!");
                  treeLiksMatrix[0][2].printThings();
                  throw new RuntimeException("whoa!");
              }


        }  */
        /*for(int i = 0; i < treeLiksMatrix.length; i++){
                for(int j = 0; j < treeLiksMatrix[i].length; j++){
                    if(treeLiksMatrix[i][j] != null)
                        System.out.println("treeLik: "+i+" "+j+" "+treeLiksMatrix[i][j].m_siteModel.isDirtyCalculation()+" "
                                +(dpSiteModel.storedSiteModelsMatrix[i][j] == treeLiksMatrix[i][j].m_siteModel));
                }
            } */
        boolean recalculate = false;
        if(dpSiteModel.isDirtyCalculation()){

            changeType = dpSiteModel.getChangeType();
            //System.out.println("treeLik requires recal!!"+changeType);
            if(changeType == ChangeType.ADDED || changeType == ChangeType.REMOVED || changeType == ChangeType.POINTER_CHANGED){
                //System.out.println("changeType: "+changeType);
                update();
            }else if(changeType == ChangeType.SPLIT || changeType == ChangeType.MERGE || changeType == ChangeType.MULTIPLE_POINTERS_CHANGED){
                //storeTreeLikelihoods();
                //System.out.println(changeType);
                updates();
            }else if(changeType == ChangeType.SPLIT_AND_VALUE_CHANGE || changeType == ChangeType.MERGE_AND_VALUE_CHANGE){
                updates();
                checkLikelihoodDirtiness(dpSiteModel.getChangedInput());


            }else if(changeType == ChangeType.POINTERS_SWAPPED){
                int[] dirtySites = dpSiteModel.getLastDirtySites();
                if(dirtySites[0] !=  dirtySites[1]){
                    updates();
                }
            }else if(changeType == ChangeType.VALUE_CHANGED){

                changeType = ChangeType.VALUE_CHANGED;

            }else{
                changeType = ChangeType.ALL;
            }
            recalculate = true;
        }else if(treeInput.get().somethingIsDirty()){
            recalculate = true;

        }else if(branchRateModelInput.get().isDirtyCalculation()){
            recalculate = true;
        }
        
        if(recalculate){
            //SiteModel siteModel = (SiteModel)treeLiks.get(0).m_siteModel;
            for(TreeLikelihood treeLik:treeLiks){
                MCMCNodeFactory.checkDirtiness(treeLik);
                //System.out.println(treeLik.m_siteModel.isDirtyCalculation());

            }


        }

        /*q = dpSiteModel.getSiteModel(0,2);
        if(q != null && treeLiksMatrix[0][2] !=null){

              if(treeLiksMatrix[0][2].getSiteModel() != q){
                  System.out.println("whoa!");
                  treeLiksMatrix[0][2].printThings();
                  throw new RuntimeException("whoa!");
              }


        }    */

        return recalculate;
    }


    public void store(){
        //storeTreeLikelihoods = true;
        storeTreeLiksMatrix();

        super.store();

    }

    protected void storeTreeLiksMatrix(){
        for(int i = 0; i < treeLiksMatrix.length;i++){
            System.arraycopy(treeLiksMatrix[i],0,storedTreeLiksMatrix[i],0, treeLiksMatrix[i].length);
        }

    }

    /*public void storeTreeLikelihoods(){
        storeTreeLikelihoods = true;
        for(int i = 0; i < treeLiksMatrix.length;i++){
            System.arraycopy(treeLiksMatrix[i],0,storedTreeLiksMatrix[i],0, treeLiksMatrix[i].length);
        }
        //super.store();
    }*/

    public void restore(){
        restoreTreeLiksMatrix();

        super.restore();
    }

    protected void restoreTreeLiksMatrix(){
        //if(storeTreeLikelihoods){
        NewWVTreeLikelihood[][] tmp = treeLiksMatrix;
        treeLiksMatrix = storedTreeLiksMatrix;
        storedTreeLiksMatrix = tmp;
        //storeTreeLikelihoods = false;
        //}

    }

    /*public double calculateLogP()throws Exception{
        QuietSiteModel q = dpSiteModel.getSiteModel(0,2);
        if(q != null && treeLiksMatrix[0][2] !=null){

              if(treeLiksMatrix[0][2].getSiteModel() != q){
                  System.out.println("whoa!");
                  treeLiksMatrix[0][2].printThings();
                  throw new RuntimeException("whoa!");
              }


        }
        return super.calculateLogP();
    }  */



}
