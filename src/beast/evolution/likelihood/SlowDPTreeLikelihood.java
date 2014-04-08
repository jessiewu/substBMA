package beast.evolution.likelihood;

import beast.core.*;
import beast.core.parameter.ChangeType;
import beast.core.parameter.DPValuable;
import beast.evolution.sitemodel.DPNtdSiteModel;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.sitemodel.DPSiteModel;
import beast.evolution.alignment.Alignment;
import beast.evolution.substitutionmodel.SwitchingNtdBMA;
import beast.evolution.tree.Tree;
import beast.evolution.branchratemodel.BranchRateModel;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

/**
 * @author Chieh-Hsi Wu
 */

@Description("Does a lot of unnecessary calculation (calculates likelihoods at sites that are not required). Used for testing.")
public class SlowDPTreeLikelihood extends DPTreeLikelihood implements PluginList {

    private DPSiteModel dpSiteModel;
    protected ArrayList<WVTreeLikelihood> treeLiks = new ArrayList<WVTreeLikelihood>();
    protected ArrayList<WVTreeLikelihood> storedTreeLiks = new ArrayList<WVTreeLikelihood>();
    protected ChangeType changeType = ChangeType.ALL;

    /*public Input<DPSiteModel> dpSiteModelInput = new Input<DPSiteModel>(
            "siteModelList",
            "array which points a set of unique parameter values",
            Input.Validate.REQUIRED
    );

    public Input<Alignment> alignmentInput = new Input<Alignment>(
            "data",
            "sequence data for the beast.tree",
            Input.Validate.REQUIRED
    );

    public Input<Tree> treeInput = new Input<Tree>(
            "tree",
            "phylogenetic beast.tree with sequence data in the leafs",
            Input.Validate.REQUIRED
    );

    public Input<BranchRateModel.Base> branchRateModelInput = new Input<BranchRateModel.Base>(
            "branchRateModel",
            "A model describing the rates on the branches of the beast.tree."
    );
    public Input<Boolean> useAmbiguitiesInput = new Input<Boolean>(
            "useAmbiguities",
            "flag to indicate leafs that sites containing ambigue states should be handled instead of ignored (the default)",
            false
    );

    public Input<DPValuable> dpValInput = new Input<DPValuable>(
            "dpVal",
            "The object that stores the information on clustering.",
            Input.Validate.REQUIRED
    );*/

    /** calculation engine **/

    //private ArrayList<int[]> clusterWeights;

    protected Alignment alignment;
    private DPValuable dpVal;

    public void initAndValidate() throws Exception{
        if(!(m_pSiteModel.get() instanceof DPSiteModel)){
            throw new RuntimeException("DPSiteModel object required for site model.");
        }
        dpSiteModel = (DPSiteModel)m_pSiteModel.get();


        alignment = m_data.get();
        int patternCount = alignment.getPatternCount();



        int[][] clusterWeights = new int[dpSiteModel.getDimension()][patternCount];

        int siteModelCount = dpSiteModel.getSiteModelCount();

        int siteCount = alignment.getSiteCount();
        for(int i = 0; i < siteCount; i++){
            //System.err.println("substModelIndices[i]: "+dpNtdSiteModel.getSubstCurrCluster(i));
            //System.out.println(dpSiteModel.getCurrCluster(i)+" "+alignment.getPatternIndex(i));
            clusterWeights[dpSiteModel.getCurrCluster(i)][alignment.getPatternIndex(i)]++;
        }

        for(int i = 0; i < siteModelCount;i++){
            //WVTreeLikelihood treeLik = new WVTreeLikelihood(clusterWeights[i]);
            /*System.out.print("cluster weights:");
            for(int j = 0; j<clusterWeights[i].length;j++){
                System.out.print(clusterWeights[i][j]+" ");

            }
            System.out.println(); */
            WVTreeLikelihood treeLik = new WVTreeLikelihood(clusterWeights[i]);
            treeLik.initByName(
                    "data", alignment,
                    "tree", m_tree.get(),
                    "siteModel", dpSiteModel.getSiteModel(i),
                    "branchRateModel", m_pBranchRateModel.get(),
                    "useAmbiguities",useAmbiguitiesInput.get()
            );
            treeLiks.add(treeLik);

            
        }
        dpVal = dpValInput.get();

    }

    public int getDimension(){
        return treeLiks.size();
    }

    @Override
    public double calculateLogP() throws Exception{
        logP = 0.0;
        //System.out.println("hello: "+treeLiks.size());
        //double sum = 0.0;
        //for(WVTreeLikelihood treeLik : treeLiks) {
        for(WVTreeLikelihood treeLik : treeLiks) {
            //sum+=treeLik.weightSum();
            //System.err.println("hello?"+treeLik.isDirtyCalculation());
            //System.out.println("siteModel id: "+ ((SwitchingNtdBMA)treeLik.m_substitutionModel).getIDNumber());
        	if (treeLik.isDirtyCalculation()) {
                //System.err.println("hello?");
                /*logP += treeLik.calculateLogP();*/
                double tmp = treeLik.calculateLogP();
        		logP += tmp;
                //System.out.println("calcLogP: "+tmp);
        	} else {
        		logP += treeLik.getCurrentLogP();
                //System.out.println("currLogP: "+treeLik.calculateLogP());
        	}
            if (Double.isInfinite(logP) || Double.isNaN(logP)) {
            	return logP;
            }
        }

        /*if(sum != alignment.getSiteCount()){
            throw new RuntimeException("WAY WRONG");
        }*/
        //System.out.println("logP: "+logP);
        return logP;
    }


    public double getSiteLogLikelihood(int iCluster, int iSite){
        return treeLiks.get(iCluster).getPatternLogLikelihood(alignment.getPatternIndex(iSite));
    }


    public void store(){
        storedTreeLiks = new ArrayList<WVTreeLikelihood>();
        for(WVTreeLikelihood treeLik : treeLiks) {
            //System.out.println("siteModel id: "+ ((SwitchingNtdBMA)treeLik.m_substitutionModel).getIDNumber());
            storedTreeLiks.add(treeLik);
            treeLik.store();
        }

        super.store();
        //System.out.println("storedLogP: "+logP);
    }

    public void restore(){
        treeLiks = storedTreeLiks;
        storedTreeLiks = null;
        for(WVTreeLikelihood treeLik : treeLiks) {
            treeLik.restore();
        }
        super.restore();
        //System.out.println("restoredLogP: "+logP);

    }

    @Override
    public List<String> getConditions() {
        List<String> conditions = new ArrayList<String>();
        for(TreeLikelihood treeLik : treeLiks) {
            conditions.addAll(treeLik.getConditions());
        }
        return conditions;
    }

    public void addTreeLikelihood(){
        //SiteModel siteModel = dpSiteModel.getLastAdded();
        SiteModel siteModel = dpSiteModel.getSiteModel(dpSiteModel.getLastAddedIndex());

        int[] patternWeights = new int[alignment.getPatternCount()];
        patternWeights[alignment.getPatternIndex(dpSiteModel.getLastDirtySite())] +=1;

        //OldWVAlignment wvalign = new OldWVAlignment(alignment, patternWeights);
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
        }catch(Exception e){
            throw new RuntimeException(e);
        }

    }

    public void splitTreeLikelihood(){
        //alignment.printPattern();
        //SiteModel siteModel = dpSiteModel.getLastAdded();
        SiteModel siteModel = dpSiteModel.getSiteModel(dpSiteModel.getLastAddedIndex());

        int[] clusterSites = dpVal.getClusterSites(dpSiteModel.getLastAddedIndex());
        int[] patternWeights = new int[alignment.getPatternCount()];

        int prevCluster = dpSiteModel.getDirtySiteModelIndex();
        //System.out.println("prevCluster: "+prevCluster);
        WVTreeLikelihood prevTreeLikelihood = treeLiks.get(prevCluster);
        for(int i = 0; i < clusterSites.length;i++){

            int patternIndex = alignment.getPatternIndex(clusterSites[i]);
            //System.out.println("clusterSites: "+clusterSites[i]+" PatIndex: "+patternIndex);
            patternWeights[patternIndex]++;
            //System.out.println(i+" splitting: "+patternWeights[patternIndex]);
            prevTreeLikelihood.removeWeight(patternIndex,1);
        }
        //OldWVAlignment wvalign = new OldWVAlignment(alignment, patternWeights);
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
        }catch(Exception e){
            throw new RuntimeException(e);
        }



    }

    public void mergeTreeLikelihoods(int removedIndex){
        WVTreeLikelihood removedTreeLikelihood =  treeLiks.remove(removedIndex);
        int [] patternWeights = removedTreeLikelihood.getPatternWeights();
        WVTreeLikelihood mergedTreeLikelihood =  treeLiks.get(dpSiteModel.getDirtySiteModelIndex());
        for(int i = 0;i < patternWeights.length;i++){
            mergedTreeLikelihood.addWeight(i,patternWeights[i]);
        }
    }

    public void removeTreeLikelihood(int removedIndex){
        treeLiks.remove(removedIndex);
    }

    protected void updateWeights(){
        int dirtySite = dpSiteModel.getLastDirtySite();
        //System.out.println("changeType: "+changeType);
        if(changeType == ChangeType.ADDED){
            int prevCluster = dpSiteModel.getPrevCluster(dirtySite);
            //System.out.println(prevCluster);
            treeLiks.get(prevCluster).removeWeight(alignment.getPatternIndex(dirtySite),1);
        }else if(changeType==ChangeType.POINTER_CHANGED){
            int prevCluster = dpSiteModel.getPrevCluster(dirtySite);
            treeLiks.get(prevCluster).removeWeight(alignment.getPatternIndex(dirtySite),1);

            int currCluster = dpSiteModel.getCurrCluster(dirtySite);
            treeLiks.get(currCluster).addWeight(alignment.getPatternIndex(dirtySite),1);
            //System.out.println(prevCluster+" "+currCluster);
        }else if(changeType == ChangeType.REMOVED){
            int currCluster = dpSiteModel.getCurrCluster(dirtySite);
            treeLiks.get(currCluster).addWeight(alignment.getPatternIndex(dirtySite),1);
        }
    }


    @Override
    protected boolean requiresRecalculation() {
        boolean recalculate = false;
        if(dpSiteModel.isDirtyCalculation()){

            ChangeType changeType = dpSiteModel.getChangeType();
            //System.out.println("treeLik requires recal!!"+changeType);
            if(changeType == ChangeType.ADDED){
                //System.out.println("added!!");
                addTreeLikelihood();
                this.changeType = ChangeType.ADDED;
                updateWeights();

            }else if(changeType == ChangeType.REMOVED){
                //System.out.println("removed!!");
                removeTreeLikelihood(dpSiteModel.getRemovedIndex());
                this.changeType = ChangeType.REMOVED;
                updateWeights();
            }else if(changeType == ChangeType.SPLIT){
                splitTreeLikelihood();
                this.changeType = ChangeType.SPLIT;
            }else if(changeType == ChangeType.MERGE){
                mergeTreeLikelihoods(dpSiteModel.getRemovedIndex());
                this.changeType = ChangeType.MERGE;
            }else if (changeType == ChangeType.POINTER_CHANGED){
                this.changeType = ChangeType.POINTER_CHANGED;
                updateWeights();
            }else if(changeType == ChangeType.VALUE_CHANGED){
                this.changeType = ChangeType.VALUE_CHANGED;

            }else{
                this.changeType = ChangeType.ALL;
            }


            
            recalculate = true;
        }else if(m_tree.get().somethingIsDirty()){
            recalculate = true;

        }else if(m_pBranchRateModel.get().isDirtyCalculation()){
            recalculate = true;
        }
        if(recalculate){
            for(WVTreeLikelihood treeLik:treeLiks){
                MCMCNodeFactory.checkDirtiness(treeLik);
            }
        }
        return recalculate;
    }





    @Override
    public void sample(State state, Random random) {
        throw new UnsupportedOperationException("Can't sample a fixed alignment!");
    }

    @Override
    public List<String> getArguments() {
        List<String> arguments = new ArrayList<String>();
        for(TreeLikelihood treeLik : treeLiks) {
            arguments.addAll(treeLik.getArguments());
        }
        return arguments;
    }

    public int[][] getClusterWeights(){
        int[][] clusterWeights = new int[treeLiks.size()][];
        for(int i = 0; i < clusterWeights.length;i++){
            clusterWeights[i] = treeLiks.get(i).getPatternWeights();
        }
        return clusterWeights;
    }

    public Tree getTree(){
        return m_tree.get();
    }

    public int[] getClusterWeights(int clusterIndex){
        return Arrays.copyOf(treeLiks.get(clusterIndex).getPatternWeights(), alignment.getPatternCount());

    }



    protected void accept() {
        //for(WVTreeLikelihood treeLik:treeLiks){
        for(WVTreeLikelihood treeLik:treeLiks){
            treeLik.makeAccept();
        }
    	super.accept();
    }
}
