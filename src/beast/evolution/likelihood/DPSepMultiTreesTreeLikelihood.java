package beast.evolution.likelihood;

import beast.app.BeastMCMC;
import beast.core.Description;
import beast.core.Input;
import beast.evolution.alignment.Alignment;
import beast.evolution.sitemodel.DPNtdRateSepSiteModel;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.SwitchingNtdBMA;
import beast.evolution.tree.Tree;
import org.omg.CORBA.PUBLIC_MEMBER;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Chieh-Hsi Wu
 *
 */
@Description("Tree likelihood that supports independent Dirichlet process priors on the partitionings of substitution model and rate with multiple alignments/trees.")
public class DPSepMultiTreesTreeLikelihood extends DPSepTreeLikelihood {
    private NewWVTreeLikelihood[][][] treeLiksMatrix;
    private NewWVTreeLikelihood[][][] storedTreeLiksMatrix;
    private int[][][] treeLikWeightMatrix;
    private int[][][] storedTreeLikWeightMatrix;

    public Input<List<Alignment>> alignmentsInput = new Input<List<Alignment>>(
            "dataAlignment",
            "A list of nucleotide alignments used as a the data.",
            new ArrayList<Alignment>(),
            Input.Validate.REQUIRED
    );

    public Input<List<Tree>> treesInput = new Input<List<Tree>>(
            "trees",
            "A list of trees which describes the phylogenies under which the data has been generated.",
            new ArrayList<Tree>(),
            Input.Validate.REQUIRED
    );



    public DPSepMultiTreesTreeLikelihood(){
        /*
         * We need DPNtdRateSepSiteModel specifically,
         * therefore an Input<DPSiteModel> is not useful.
         */

        dpValInput.setRule(Input.Validate.OPTIONAL);
        dataInput.setRule(Input.Validate.REQUIRED);
        treeInput.setRule(Input.Validate.REQUIRED);
    }

    private List<Alignment> alignments;
    private List<Tree> trees;
    private int[] alignmentStartingIndex;
    private int[] alignmentIndexBySite;
    public void initAndValidate() {
        useThreads = useThreadsInput.get() && (BeastMCMC.m_nThreads > 1);
        useThreadsEvenly = useThreadsEvenlyInput.get() && (BeastMCMC.m_nThreads > 1);

        alignments = alignmentsInput.get();
        alignmentStartingIndex = new int[alignments.size()];
        alignmentStartingIndex[0] = 0;
        for(int i = 1; i < alignmentStartingIndex.length; i++){
            alignmentStartingIndex[i] = alignments.get(i - 1).getSiteCount()+alignmentStartingIndex[0];
        }
        int totalSiteCount = 0;
        for(int i = 0; i < alignments.size(); i++){
            totalSiteCount += alignments.get(i).getSiteCount();
        }

        alignmentIndexBySite = new int[totalSiteCount];

        int m = 0;
        for(int i = 1; i < alignmentStartingIndex.length;i++){
            while(m < alignmentStartingIndex[i]){
                alignmentIndexBySite[m++] = i - 1;
            }
        }
        for(; m < alignmentIndexBySite.length; m++){
            alignmentIndexBySite[m]= alignmentStartingIndex.length -1;
        }

        //int patternCount = alignment.getPatternCount();
        if(!(siteModelInput.get() instanceof DPNtdRateSepSiteModel)){
            throw new RuntimeException("DPNtdRateSepSiteModel required for site model.");
        }
        dpSiteModel = (DPNtdRateSepSiteModel)siteModelInput.get();
        int siteModelCount = dpSiteModel.getSiteModelCount();

        /*
         * Set up a 4 dimensional array to store the weights of each ntdBMA/rate combination (dimensions: ntdBMA, rates, weights).
         */
        int[][][][]clusterWeights = new int[alignments.size()][dpSiteModel.getSubstClusterLimit()][dpSiteModel.getRatesClusterLimit()][];

        int patternCount;
        for(int i = 0; i < clusterWeights.length; i++){
            patternCount = alignments.get(i).getPatternCount();
            for(int j = 0; j < clusterWeights.length; j++){
                for(int k = 0; k < clusterWeights.length; k++){
                    clusterWeights[i][j][k] = new int[patternCount];
                }
            }
        }



        treeLikWeightMatrix = new int[alignments.size()][dpSiteModel.getSubstClusterLimit()][dpSiteModel.getRatesClusterLimit()];
        storedTreeLikWeightMatrix = new int[alignments.size()][dpSiteModel.getSubstClusterLimit()][dpSiteModel.getRatesClusterLimit()];
        int[] clusterIds;
        int siteCount;
        for(int i = 0; i < alignments.size(); i++){
            siteCount = alignments.get(i).getSiteCount();
            for(int j = 0; j < siteCount; j++){
                clusterIds = dpSiteModel.getCurrClusters(i);
                clusterWeights[i][clusterIds[DPNtdRateSepSiteModel.NTDBMA]][clusterIds[DPNtdRateSepSiteModel.RATES]][alignmentStartingIndex[i]+j]++;
                treeLikWeightMatrix[i][clusterIds[DPNtdRateSepSiteModel.NTDBMA]][clusterIds[DPNtdRateSepSiteModel.RATES]]++;
            }

        }


        /*
         *  Traverse through the list site models and create tree likelihoods.
         *  The tree likelihoods are then added to a list and the matrix.
         *  The order of likelihoods in the list corresponds to that of site models in the DPSiteModel list.
         */
        treeLiksMatrix = new NewWVTreeLikelihood[alignments.size()][dpSiteModel.getSubstClusterLimit()][dpSiteModel.getRatesClusterLimit()];
        storedTreeLiksMatrix = new NewWVTreeLikelihood[alignments.size()][dpSiteModel.getSubstClusterLimit()][dpSiteModel.getRatesClusterLimit()];


        for(int i = 0; i < clusterWeights.length; i++){
            for(int j = 0; j < clusterWeights[i].length; j++){
                for(int k = 0; k < clusterWeights[i][j].length; k++){
                    if(treeLikWeightMatrix[i][j][k] > 0){
                        NewWVTreeLikelihood treeLik = new NewWVTreeLikelihood(
                                clusterWeights[i][j][k],
                                alignments.get(i),
                                trees.get(i),
                                useAmbiguitiesInput.get(),
                                dpSiteModel.getSiteModel(i,j),
                                branchRateModelInput.get()

                        );

                        treeLiks.add(treeLik);
                        treeLiksMatrix[i][j][k] = treeLik;
                        storedTreeLiksMatrix[i][j][k] = treeLik;
                    }
                }
            }

        }

    }



    protected void update(){
        int dirtySite = dpSiteModel.getLastDirtySite();
        int[] currClusterIds = dpSiteModel.getCurrClusters(dirtySite);
        int[] prevClusterIds = dpSiteModel.getPrevClusters(dirtySite);


        //Create a new tree likelihood (with zero weights) when there is a new site model.
        if(treeLiksMatrix[alignmentIndexBySite[dirtySite]][currClusterIds[DPNtdRateSepSiteModel.NTDBMA]][currClusterIds[DPNtdRateSepSiteModel.RATES]] == null){
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
        if(treeLikWeightMatrix[alignmentIndexBySite[dirtySite]][prevClusterIds[DPNtdRateSepSiteModel.NTDBMA]][prevClusterIds[DPNtdRateSepSiteModel.RATES]] == 0){
            treeLiks.remove(treeLiksMatrix[alignmentIndexBySite[dirtySite]][prevClusterIds[DPNtdRateSepSiteModel.NTDBMA]][prevClusterIds[DPNtdRateSepSiteModel.RATES]]);
            treeLiksMatrix[alignmentIndexBySite[dirtySite]][prevClusterIds[DPNtdRateSepSiteModel.NTDBMA]][prevClusterIds[DPNtdRateSepSiteModel.RATES]] = null;
        }
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

            if(treeLikWeightMatrix[alignmentIndexBySite[dirtySite]][prevClusterIds[DPNtdRateSepSiteModel.NTDBMA]][prevClusterIds[DPNtdRateSepSiteModel.RATES]] == 0){
                //System.out.println("Remove?");
                treeLiks.remove(treeLiksMatrix[prevClusterIds[DPNtdRateSepSiteModel.NTDBMA]][prevClusterIds[DPNtdRateSepSiteModel.RATES]]);
                treeLiksMatrix[alignmentIndexBySite[dirtySite]][prevClusterIds[DPNtdRateSepSiteModel.NTDBMA]][prevClusterIds[DPNtdRateSepSiteModel.RATES]] = null;
            }
        }
    }

    protected void update(int dirtySite){
        //int dirtySite = dpSiteModel.getLastDirtySite();
        int[] currClusterIds = dpSiteModel.getCurrClusters(dirtySite);
        int[] prevClusterIds = dpSiteModel.getPrevClusters(dirtySite);

        //Create a new tree likelihood (with zero weights) when there is a new site model.
        if(treeLiksMatrix[alignmentIndexBySite[dirtySite]][currClusterIds[DPNtdRateSepSiteModel.NTDBMA]][currClusterIds[DPNtdRateSepSiteModel.RATES]] == null){
            //System.out.println("Add cluster at: "+currClusterIds[DPNtdRateSepSiteModel.NTDBMA]+" "+currClusterIds[DPNtdRateSepSiteModel.RATES]);
            addTreeLikelihood(
                    alignmentIndexBySite[dirtySite],
                    currClusterIds[DPNtdRateSepSiteModel.NTDBMA],
                    currClusterIds[DPNtdRateSepSiteModel.RATES]
            );
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

    public void addTreeLikelihood(int alignmentIndex, int substModelID, int rateID){
        //SiteModel siteModel = dpSiteModel.getLastAdded();
        SiteModel siteModel = dpSiteModel.getSiteModel(substModelID, rateID);


        int[] patternWeights = new int[alignments.get(alignmentIndex).getPatternCount()];

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
            treeLiksMatrix[alignmentIndex][substModelID][rateID] = treeLik;
        }catch(Exception e){
            throw new RuntimeException(e);
        }

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
        int siteIndexInAlignment = dirtySite - alignmentStartingIndex[alignmentIndexBySite[dirtySite]];
        int patternIndex = alignments.get(alignmentIndexBySite[dirtySite]).getPatternIndex(siteIndexInAlignment);
        treeLiksMatrix[alignmentIndexBySite[dirtySite]][ntdBMAIndexFrom][rateIndexFrom].removeWeight(patternIndex,weight);
        treeLiksMatrix[alignmentIndexBySite[dirtySite]][ntdBMAIndexTo][rateIndexTo].addWeight(patternIndex,weight);

        //Update the weights for each treelikelihood.
        treeLikWeightMatrix[alignmentIndexBySite[dirtySite]][ntdBMAIndexFrom][rateIndexFrom]-=weight;
        treeLikWeightMatrix[alignmentIndexBySite[dirtySite]][ntdBMAIndexTo][rateIndexTo]+=weight;

        assert treeLikWeightMatrix[alignmentIndexBySite[dirtySite]][ntdBMAIndexFrom][rateIndexFrom] > -1;
        assert treeLikWeightMatrix[alignmentIndexBySite[dirtySite]][ntdBMAIndexTo][rateIndexTo] > -1;
        //}catch(Exception e){
        //}
    }

    /*
     * @param inputType either ntdBMA or rates
     * @clusterID id of the cluster
     * @siteIndex the index of a site in an alignment
     */
    @Override
    public double getSiteLogLikelihood(int inputType, int clusterID, int siteIndex){
        //System.out.println("pattern: "+alignment.getPatternIndex(siteIndex));
        int[] currClusters = dpSiteModel.getCurrClusters(siteIndex);
        //System.out.println("clusterID: "+clusterID+" "+prevClusters[DPNtdRateSepSiteModel.RATES]+" "+alignment.getPatternIndex(siteIndex));
        if(inputType == DPNtdRateSepSiteModel.NTDBMA){
            if(treeLiksMatrix[alignmentIndexBySite[siteIndex]][clusterID][currClusters[DPNtdRateSepSiteModel.RATES]] != null){

                //WVTreeLikelihood tmpTL = treeLiksMatrix[clusterID][prevClusters[DPNtdRateSepSiteModel.RATES]];
                NewWVTreeLikelihood tmpTL = treeLiksMatrix[alignmentIndexBySite[siteIndex]][clusterID][currClusters[DPNtdRateSepSiteModel.RATES]];
                //System.out.println("clusterID: "+clusterID+" "+prevClusters[DPNtdRateSepSiteModel.RATES]+" "+alignment.getPatternIndex(siteIndex));
                //tmpTL.printThings();
                return tmpTL.getPatternLogLikelihood(alignment.getPatternIndex(siteIndex));

            }
        }else{

            if(treeLiksMatrix[alignmentIndexBySite[siteIndex]][currClusters[DPNtdRateSepSiteModel.NTDBMA]][clusterID] != null){
                //WVTreeLikelihood tmpTL = treeLiksMatrix[prevClusters[DPNtdRateSepSiteModel.NTDBMA]][clusterID];
                //System.out.println("hi!!");
                NewWVTreeLikelihood tmpTL = treeLiksMatrix[alignmentIndexBySite[siteIndex]][currClusters[DPNtdRateSepSiteModel.NTDBMA]][clusterID];
                return tmpTL.getPatternLogLikelihood(alignment.getPatternIndex(siteIndex));
            }
        }

        return Double.NaN;

    }


    public void store(){
        for(int i = 0; i < treeLikWeightMatrix.length;i++){
            for(int j = 0; j < treeLikWeightMatrix[i].length; j++){
                System.arraycopy(treeLikWeightMatrix[i][j],0,storedTreeLikWeightMatrix[i][j],0, treeLikWeightMatrix[i][j].length);
            }
        }

        super.store();
    }

    protected void storeTreeLiksMatrix(){
        for(int i = 0; i < treeLiksMatrix.length;i++){
            for(int j = 0; j < treeLiksMatrix[i].length; j++){
                System.arraycopy(treeLiksMatrix[i][j],0,storedTreeLiksMatrix[i][j],0, treeLiksMatrix[i][j].length);
            }
        }
    }

    public void restore(){
        int[][][] tmp  = treeLikWeightMatrix;
        treeLikWeightMatrix = storedTreeLikWeightMatrix;
        storedTreeLikWeightMatrix = tmp;
        super.restore();
    }

    protected void restoreTreeLiksMatrix(){
        //if(storeTreeLikelihoods){
        NewWVTreeLikelihood[][][] tmp = treeLiksMatrix;
        treeLiksMatrix = storedTreeLiksMatrix;
        storedTreeLiksMatrix = tmp;
        //storeTreeLikelihoods = false;
        //}
    }
}
