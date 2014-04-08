package beast.evolution.likelihood;

import beast.core.Description;
import beast.core.Input;
import beast.core.MCMCNodeFactory;
import beast.core.parameter.ChangeType;
import beast.evolution.alignment.Alignment;
import beast.evolution.sitemodel.DPMultiAlignSiteModel;
import beast.evolution.sitemodel.DPSiteModel;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.tree.Tree;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * @author Chieh-Hsi Wu
 */
@Description("DPTreeLikelihood that allows to take up multiple alignments/trees.")
public class DPMultiTreesTreeLikelihood extends DPTreeLikelihood{
    public Input<List<Tree>> treesInput = new Input<List<Tree>>(
            "trees",
            "A list of unlinked trees.",
            new ArrayList<Tree>(),
            Input.Validate.REQUIRED
    );

    public Input<List<Alignment>> alignmentsInput = new Input<List<Alignment>>(
            "dataAlignment",
            "sequence data for the beast.tree",
            new ArrayList<Alignment>(),
            Input.Validate.REQUIRED
    );

    public DPMultiTreesTreeLikelihood(){
        m_tree.setRule(Input.Validate.OPTIONAL);
        m_data.setRule(Input.Validate.OPTIONAL);
    }

    private HashMap<Integer,NewWVTreeLikelihood>[] treeLikelihoodMap;
    private HashMap<Integer,NewWVTreeLikelihood>[] storedTreeLikelihoodMap;
    private List<Alignment> alignments;
    private List<Tree> trees;
    private int[] siteIndexWithinAlignment;
    private int[] alignmentIndexBySite;
    private HashMap<Integer, Integer>[] likelihoodWeight;     //todo store restore
    private HashMap<Integer, Integer>[] storedLikelihoodWeight;
    public void initAndValidate(){
        //Get the required DPSiteModel
        if(!(m_pSiteModel.get() instanceof DPSiteModel)){
            throw new RuntimeException("DPSiteModel required for site model.");
        }
        dpSiteModel = (DPSiteModel) m_pSiteModel.get();


        trees = treesInput.get();
        alignments = alignmentsInput.get();
        dpVal = dpValInput.get();
        int alignmentCount = alignments.size();



        int[] prevAlignEndIndex = new int[alignments.size()];

        prevAlignEndIndex[0] = 0;

        for(int i = 1; i < prevAlignEndIndex.length; i++){
            prevAlignEndIndex[i] = prevAlignEndIndex[i - 1] + alignments.get(i - 1).getSiteCount();
        }


        int siteCount = prevAlignEndIndex[prevAlignEndIndex.length - 1] +
                alignments.get(alignments.size() - 1).getSiteCount();

        alignmentIndexBySite = new int[siteCount];

        int k = 0;
        for(int i = 1; i< prevAlignEndIndex.length;i++){
            while(k < prevAlignEndIndex[i]){
                alignmentIndexBySite[k++] = i - 1;
            }
        }
        while(k < siteCount){
            alignmentIndexBySite[k++] = prevAlignEndIndex.length - 1;
        }

        siteIndexWithinAlignment = new int[siteCount];
        for(int i = 0; i < siteIndexWithinAlignment.length; i++){
            siteIndexWithinAlignment[i] = i - prevAlignEndIndex[alignmentIndexBySite[i]];

        }



        int categoryCount = dpVal.getCategoryCount();
        int[][][] clusterPatternWeights = new int[alignmentCount][categoryCount][];
        for(int i = 0; i < categoryCount; i++){


            //Get all the sites in that category
            int[] sites = dpVal.getClusterSites(i);

            for(int j = 0; j < sites.length; j++){

                int alignmentIndex = alignmentIndexBySite[sites[j]];

                //Only initializing the array when the alignment has site(s) in that category.

                if(clusterPatternWeights[alignmentIndex][i] == null){
                    int patternCount = alignments.get(alignmentIndex).getPatternCount();
                    clusterPatternWeights[alignmentIndex][i] = new int[patternCount];
                }

                //Already initialized
                int site = siteIndexWithinAlignment[sites[j]];
                int patternIndex = alignments.get(alignmentIndex).getPatternIndex(site);
                clusterPatternWeights[alignmentIndex][i][patternIndex]++;
            }


        }


        treeLikelihoodMap =  (HashMap<Integer,NewWVTreeLikelihood>[])new HashMap[alignmentCount];
        storedTreeLikelihoodMap =  (HashMap<Integer,NewWVTreeLikelihood>[])new HashMap[alignmentCount];

        try{
            for(int i = 0; i < treeLikelihoodMap.length; i++){
                Alignment alignment = alignments.get(i);
                treeLikelihoodMap[i] = new HashMap<Integer,NewWVTreeLikelihood>();
                //System.out.println(i+" "+(treeLikelihoodMap[i] == null));         DP
                for(int j = 0; j < clusterPatternWeights[i].length; j++){
                    if(clusterPatternWeights[i][j] != null){
                        NewWVTreeLikelihood treeLik = new NewWVTreeLikelihood(
                                clusterPatternWeights[i][j],
                                alignment,
                                trees.get(i),
                                useAmbiguitiesInput.get(),
                                dpSiteModel.getSiteModel(i),
                                m_pBranchRateModel.get()
                        );

                        treeLik.calculateLogP();
                        treeLik.store();
                        treeLiks.add(treeLik);
                        //System.out.println(i + " " + (treeLikelihoodMap[i] == null));
                        treeLikelihoodMap[i].put(dpVal.getCategoryIDNumber(j), treeLik);

                        int likWeight = 0;
                        for(int w = 0; w < clusterPatternWeights[i][j].length; w++){
                            likWeight+=clusterPatternWeights[i][j][w];
                        }
                        likelihoodWeight[i].put(dpVal.getCategoryIDNumber(j),likWeight);

                    }
                }

            }
        }catch(Exception e){
            throw new RuntimeException(e);
        }

    }



    public int getAlignmentIndexBySite(int siteIndex){
        return alignmentIndexBySite[siteIndex];
    }

    private void update(){
        int dirtySite = dpSiteModel.getLastDirtySite();
        update(dirtySite);
        int prevCategory = dpSiteModel.getPrevCategoryIDNumber(dirtySite);
        checkAndRemove(alignmentIndexBySite[dirtySite],prevCategory);
    }

    private void update(int[] dirtySites){
        for(int dirtySite:dirtySites){
            update(dirtySite);
        }

        int prevCategory = dpSiteModel.getPrevCategoryIDNumber(dirtySites[0]);

        //Remove likelihoods that have zero weights
        for(int i = 0; i < alignments.size();i++){
            checkAndRemove(i,prevCategory);
        }


    }

    private void update(int dirtySite){

        int alignmentIndex = alignmentIndexBySite[dirtySite];
        int currCatIDNum = dpSiteModel.getCurrCategoryIDNumber(dirtySite);
        int prevCatIDNum = dpSiteModel.getPrevCategoryIDNumber(dirtySite);
        //System.out.println(getClass()+": "+dirtySite+" "+currCategory+" alignmentIndex: "+alignmentIndex);
        //System.out.println(treeLiks.size());
        if(!treeLikelihoodMap[alignmentIndex].containsKey(currCatIDNum)){
            addTreeLikelihood(alignmentIndex, currCatIDNum,dirtySite);
            //System.out.println("Add? "+treeLikelihoodMap[alignmentIndex].get(currCategory).m_data.get()+" "+alignments.get(alignmentIndex).getPatternIndex(siteIndexWithinAlignment[dirtySite]));
        }
        //System.out.println(treeLiks.size());

        moveWeight(
            alignmentIndex,
            currCatIDNum,
            prevCatIDNum,
            siteIndexWithinAlignment[dirtySite],
            1
        );



    }

    public void moveWeight(
            int alignmentIndex,
            int currCatIDNum,
            int prevCatIDNum,
            int dirtySite,
            int weight){
        //System.out.println(alignmentIndex+" "+prevCategory+" "+ currCategory+" "+weight);
        //System.out.println("prevCat: "+prevCategory+" "+(treeLikelihoodMap[alignmentIndex].get(prevCategory)==null));

        treeLikelihoodMap[alignmentIndex].get(prevCatIDNum).removeWeight(
                alignments.get(alignmentIndex).getPatternIndex(dirtySite),
                weight
        );

        treeLikelihoodMap[alignmentIndex].get(currCatIDNum).addWeight(
                alignments.get(alignmentIndex).getPatternIndex(dirtySite),
                weight
        );

        if(likelihoodWeight[alignmentIndex].containsKey(currCatIDNum)){
            int currCatOldWeight = likelihoodWeight[alignmentIndex].get(currCatIDNum);
            likelihoodWeight[alignmentIndex].put(currCatIDNum,currCatOldWeight+weight);

            int prevCatOldWeight = likelihoodWeight[alignmentIndex].get(prevCatIDNum);
            likelihoodWeight[alignmentIndex].put(currCatIDNum,prevCatOldWeight-weight);


            //This is here for sanity checks
            if(likelihoodWeight[alignmentIndex].get(prevCatIDNum) == 0){
                likelihoodWeight[alignmentIndex].remove(prevCatIDNum);
            }



        }else{
            likelihoodWeight[alignmentIndex].put(currCatIDNum, 1);

        }

    }



    public void addTreeLikelihood(int alignmentIndex, int categoryID, int dirtySite){
        //System.out.println("alignmentIndex"+alignmentIndex);
        //Get the site required siteModel
        SiteModel siteModel = dpSiteModel.getSiteModelOfSiteIndex(dirtySite);
        int[] patternWeights = new int[alignments.get(alignmentIndex).getPatternCount()];

        try{
            NewWVTreeLikelihood treeLik = new NewWVTreeLikelihood(patternWeights,
                    alignments.get(alignmentIndex),
                    trees.get(alignmentIndex),
                    useAmbiguitiesInput.get(),
                    siteModel,
                    m_pBranchRateModel.get());

            treeLik.calculateLogP();
            treeLik.store();
            treeLiks.add(treeLik);
            treeLikelihoodMap[alignmentIndex].put(categoryID, treeLik);
        }catch(Exception e){
            throw new RuntimeException(e);
        }

    }

    private void checkAndRemove(int alignmentIndex, int prevCatIDNum){

        if(likelihoodWeight[alignmentIndex].get(prevCatIDNum) == 0){
            treeLiks.remove(treeLikelihoodMap[alignmentIndex].get(prevCatIDNum));
            treeLikelihoodMap[alignmentIndex].remove(prevCatIDNum);
        }
    }

    private void handleSplit(){
        int[] dirtySites = dpVal.getClusterSites(dpVal.getDimension() - 1);
        update(dirtySites);
    }

    private void handleMerge(){
        int[] dirtySites = dpVal.getMergedSites();
        update(dirtySites);
    }

    private void handleMultiPointerChanges(){
        int[] dirtySites = dpVal.getLastDirtySites();
        update(dirtySites);
    }

    public void store(){
        for(int i = 0; i < treeLikelihoodMap.length;i++){
            storedTreeLikelihoodMap[i] = (HashMap<Integer,NewWVTreeLikelihood>)treeLikelihoodMap[i].clone();
        }
        for(int i = 0; i < likelihoodWeight.length; i++){
            storedLikelihoodWeight[i] = (HashMap<Integer,Integer>)likelihoodWeight[i].clone();
        }
        super.store();
    }

    public void restore(){
        HashMap<Integer,NewWVTreeLikelihood>[] tmp = treeLikelihoodMap;
        treeLikelihoodMap = storedTreeLikelihoodMap;
        storedTreeLikelihoodMap = tmp;

        HashMap<Integer,Integer>[] tmp1 = likelihoodWeight;
        likelihoodWeight = storedLikelihoodWeight;
        storedLikelihoodWeight = tmp1;
        super.restore();

    }

    @Override
    protected boolean requiresRecalculation() {

        boolean recalculate = false;

        if(dpSiteModel.isDirtyCalculation()){

            changeType = dpSiteModel.getChangeType();
            //System.out.println("treeLik requires recal!!"+changeType);
            if(changeType == ChangeType.ADDED || changeType == ChangeType.REMOVED || changeType == ChangeType.POINTER_CHANGED){
                //System.out.println("changeType: "+changeType);

                update(); //Update the single dirty site.

            }else if(changeType == ChangeType.SPLIT){
                //storeTreeLikelihoods();
                handleSplit();

            }else if(changeType == ChangeType.MERGE){

                handleMerge();
            }else if(changeType == ChangeType.MULTIPLE_POINTERS_CHANGED){
                handleMultiPointerChanges();
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
            for(TreeLikelihood treeLik:treeLiks){
                MCMCNodeFactory.checkDirtiness(treeLik);
            }
        }

        return recalculate;
    }
}
