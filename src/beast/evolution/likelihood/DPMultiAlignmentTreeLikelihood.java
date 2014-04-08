package beast.evolution.likelihood;

import beast.core.Description;
import beast.core.Input;
import beast.core.MCMCNodeFactory;
import beast.core.parameter.ChangeType;
import beast.evolution.alignment.Alignment;
import beast.evolution.sitemodel.DPMultiAlignSiteModel;
import beast.evolution.sitemodel.SiteModel;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Set;

/**
 * @author Chieh-Hsi
 */
@Description("A class that is able compute the likelihood for multiple alignments when the across-site rates are also partitioned.")
public class DPMultiAlignmentTreeLikelihood extends DPTreeLikelihood {
    private HashMap<Integer,NewWVTreeLikelihood>[] treeLikelihoodMap;
    private HashMap<Integer,NewWVTreeLikelihood>[] storedTreeLikelihoodMap;
    private List<Alignment> alignments;


    public Input<List<Alignment>> alignmentsInput = new Input<List<Alignment>>(
            "dataAlignment",
            "A list of alignments that cannot be combined into one alignment",
            new ArrayList<Alignment>(),
            Input.Validate.REQUIRED
    );

    public Input<DPMultiAlignSiteModel> dpMultiAlignSiteModelInput = new Input<DPMultiAlignSiteModel>(
            "dpMultiAlignSiteModel",
            "A site model that can handle multiple alignments that cannot be combined.",
            Input.Validate.REQUIRED
    );



    private DPMultiAlignSiteModel dpSiteModel;
    private int[] siteIndexWithinAlignment;

    public DPMultiAlignmentTreeLikelihood(){
        m_data.setRule(Input.Validate.OPTIONAL);
        m_pSiteModel.setRule(Input.Validate.OPTIONAL);

    }

    public void initAndValidate () throws Exception{
        dpSiteModel = dpMultiAlignSiteModelInput.get();
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

        siteIndexWithinAlignment = new int[siteCount];
        for(int i = 0; i < siteIndexWithinAlignment.length; i++){
            siteIndexWithinAlignment[i] = i - prevAlignEndIndex[dpSiteModel.getAlignmentIndex(i)];

        }



        int categoryCount = dpVal.getCategoryCount();
        int[][][] clusterPatternWeights = new int[alignmentCount][categoryCount][];
        for(int i = 0; i < categoryCount; i++){


            //Get all the sites in that category
            int[] sites = dpVal.getClusterSites(i);

            for(int j = 0; j < sites.length; j++){

                int alignmentIndex = dpSiteModel.getAlignmentIndex(sites[j]);

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
        for(int i = 0; i < treeLikelihoodMap.length; i++){
            Alignment alignment = alignments.get(i);
            treeLikelihoodMap[i] = new HashMap<Integer,NewWVTreeLikelihood>();
            //System.out.println(i+" "+(treeLikelihoodMap[i] == null));         DP
            for(int j = 0; j < clusterPatternWeights[i].length; j++){
                if(clusterPatternWeights[i][j] != null){
                    NewWVTreeLikelihood treeLik = new NewWVTreeLikelihood(
                            clusterPatternWeights[i][j],
                            alignment,
                            m_tree.get(),
                            useAmbiguitiesInput.get(),
                            dpSiteModel.getSiteModel(i, dpVal.getCategoryIDNumber(j)),
                            m_pBranchRateModel.get()
                    );

                    treeLik.calculateLogP();
                    treeLik.store();
                    treeLiks.add(treeLik);
                    //System.out.println(i + " " + (treeLikelihoodMap[i] == null));
                    treeLikelihoodMap[i].put(dpVal.getCategoryIDNumber(j), treeLik);

                }
            }

        }

    }

    private void update(){
        int dirtySite = dpSiteModel.getLastDirtySite();
        update(dirtySite);
        int prevCategory = dpSiteModel.getPrevCategoryIDNumber(dirtySite);
        checkAndRemove(dpSiteModel.getAlignmentIndex(dirtySite),prevCategory);
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

        int alignmentIndex = dpSiteModel.getAlignmentIndex(dirtySite);
        int currCategory = dpSiteModel.getCurrCategoryIDNumber(dirtySite);
        int prevCategory = dpSiteModel.getPrevCategoryIDNumber(dirtySite);
        //System.out.println(getClass()+": "+dirtySite+" "+currCategory+" alignmentIndex: "+alignmentIndex);
        //System.out.println(treeLiks.size());
        if(!treeLikelihoodMap[alignmentIndex].containsKey(currCategory)){
            addTreeLikelihood(alignmentIndex, currCategory);
            //System.out.println("Add? "+treeLikelihoodMap[alignmentIndex].get(currCategory).m_data.get()+" "+alignments.get(alignmentIndex).getPatternIndex(siteIndexWithinAlignment[dirtySite]));
        }
        //System.out.println(treeLiks.size());

        moveWeight(
            alignmentIndex,
            currCategory,
            prevCategory,
            siteIndexWithinAlignment[dirtySite],
            1
        );



    }

    public void moveWeight(
            int alignmentIndex,
            int currCategory,
            int prevCategory,
            int dirtySite,
            int weight){
        //System.out.println(alignmentIndex+" "+prevCategory+" "+ currCategory+" "+weight);
        //System.out.println("prevCat: "+prevCategory+" "+(treeLikelihoodMap[alignmentIndex].get(prevCategory)==null));

        treeLikelihoodMap[alignmentIndex].get(prevCategory).removeWeight(
                alignments.get(alignmentIndex).getPatternIndex(dirtySite),
                weight
        );

        treeLikelihoodMap[alignmentIndex].get(currCategory).addWeight(
                alignments.get(alignmentIndex).getPatternIndex(dirtySite),
                weight
        );

    }

    public void addTreeLikelihood(int alignmentIndex, int categoryID){
        //System.out.println("alignmentIndex"+alignmentIndex);
        //Get the site required siteModel
        SiteModel siteModel = dpSiteModel.getSiteModel(alignmentIndex,categoryID);
        int[] patternWeights = new int[alignments.get(alignmentIndex).getPatternCount()];

        try{
            NewWVTreeLikelihood treeLik = new NewWVTreeLikelihood(patternWeights,
                    alignments.get(alignmentIndex),
                    m_tree.get(),
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

    private void checkAndRemove(int alignmentIndex, int prevCategory){

        if(dpSiteModel.getSiteModelWeight(alignmentIndex, prevCategory) == 0){
            treeLiks.remove(treeLikelihoodMap[alignmentIndex].get(prevCategory));
            treeLikelihoodMap[alignmentIndex].remove(prevCategory);
        }
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

            }else if(changeType == ChangeType.SPLIT || changeType == ChangeType.MERGE){
                //storeTreeLikelihoods();
                update(dpSiteModel.getLastDirtySites());


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

    public void store(){

        for(int i = 0; i < treeLikelihoodMap.length;i++){
            storedTreeLikelihoodMap[i] = (HashMap<Integer,NewWVTreeLikelihood>)treeLikelihoodMap[i].clone();
        }
        super.store();
    }

    public void restore(){
        //if(storeTreeLikelihoods){
        HashMap<Integer,NewWVTreeLikelihood>[] tmp = treeLikelihoodMap;
        treeLikelihoodMap = storedTreeLikelihoodMap;
        storedTreeLikelihoodMap = tmp;
        //storeTreeLikelihoods = false;
        //}
        super.restore();
    }


    public double getSiteLogLikelihood(int categoryID, int siteIndex){

        //System.out.println("pattern: "+alignment.getPatternIndex(siteIndex));
        int currCategoryIDNum = dpSiteModel.getCategoryIDNumberFromIndex(categoryID);
        int alignmentIndex = dpSiteModel.getAlignmentIndex(siteIndex);
        int sitePatternIndex = alignments.get(alignmentIndex).getPatternIndex(siteIndexWithinAlignment[siteIndex]);
        /*System.out.println("flag1: "+(treeLikelihoodMap[alignmentIndex]==null));

        System.out.println("currCategoryIDNum: "+currCategoryIDNum);  */
        //System.out.println(treeLikelihoodMap[alignmentIndex].containsKey(currCategoryIDNum));
        /*for(int i = 0 ; i < treeLikelihoodMap.length; i++ ){
            Set<Integer> keys = treeLikelihoodMap[i].keySet();
            for(int key:keys){
                System.out.print(key+" ");
            }
            System.out.println();
        }  */

        //dpSiteModel.getCategoryIDNumberFromIndex(categoryID);
        //System.out.println("siteIndex: "+siteIndex+" alignmentIndex: "+alignmentIndex+ " currCategoryIDNum: "+currCategoryIDNum+", "+
        //dpSiteModel.getCategoryIDNumberFromIndex(categoryID));
        if(treeLikelihoodMap[alignmentIndex].get(currCategoryIDNum) == null){
            return Double.NaN;
        }
        return treeLikelihoodMap[alignmentIndex].get(currCategoryIDNum).getPatternLogLikelihood(sitePatternIndex);

    }

    public int getDimension(){
        ArrayList<Integer> categoryCount = new ArrayList<Integer>();
        for(int i = 0; i < treeLikelihoodMap.length; i++){
            Set<Integer> keys = treeLikelihoodMap[i].keySet();
            for(int key:keys){
                if(!categoryCount.contains(key)){
                    categoryCount.add(key);
                }
            }
        }
        return categoryCount.size();
        //return rateList.getDimension();
    }

}
