package beast.evolution.likelihood;

import beast.app.BeastMCMC;
import beast.core.*;
import beast.core.parameter.ChangeType;
import beast.core.parameter.DPValuable;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.sitemodel.DPSiteModel;
import beast.evolution.alignment.Alignment;
import beast.evolution.tree.Tree;

import javax.sound.midi.SysexMessage;
import java.lang.Runnable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.Executors;
import java.util.concurrent.RejectedExecutionException;

/**
 * @author Chieh-Hsi Wu
 */
@Description("Tree likelihood that supports Dirichlet process mixture model on the partitionings of substitution model and rate together.")
public class DPTreeLikelihood extends GenericTreeLikelihood implements PluginList {

    protected DPSiteModel dpSiteModel;
    //protected ArrayList<WVTreeLikelihood> treeLiks = new ArrayList<WVTreeLikelihood>();
    protected ArrayList<NewWVTreeLikelihood> treeLiks = new ArrayList<NewWVTreeLikelihood>();
    protected ArrayList<NewWVTreeLikelihood> storedTreeLiks = new ArrayList<NewWVTreeLikelihood>();
    protected ChangeType changeType = ChangeType.ALL;

    public Input<Boolean> useAmbiguitiesInput = new Input<Boolean>(
            "useAmbiguities",
            "flag to indicate leafs that sites containing ambigue states should be handled instead of ignored (the default)",
            false
    );

    public Input<DPValuable> dpValInput = new Input<DPValuable>(
            "dpVal",
            "The object that stores the information on clustering.",
            Input.Validate.REQUIRED
    );

    public Input<Boolean> useThreadsInput = new Input<Boolean>("useThreads", "calculated the distributions in parallel using threads (default false)", false);
    public Input<Boolean> useThreadsEvenlyInput = new Input<Boolean>("useThreadsEvenly", "calculated the distributions in parallel using threads (default false)", false);
    boolean useThreads;
    boolean useThreadsEvenly;
    /** calculation engine **/

    //private ArrayList<int[]> clusterWeights;

    protected Alignment alignment;
    protected DPValuable dpVal;

    public void initAndValidate() {
        useThreads = useThreadsInput.get() && (BeastMCMC.m_nThreads > 1);
        useThreadsEvenly = useThreadsEvenlyInput.get() && (BeastMCMC.m_nThreads > 1);
        dpVal = dpValInput.get();
        if(!(siteModelInput.get() instanceof DPSiteModel)){
            throw new RuntimeException("DPSiteModel required for site model.");
        }
        dpSiteModel = (DPSiteModel) siteModelInput.get();


        alignment = dataInput.get();
        int patternCount = alignment.getPatternCount();



        int[][] clusterWeights = new int[dpSiteModel.getSiteModelCount()][patternCount];

        int siteModelCount = dpSiteModel.getSiteModelCount();

        int siteCount = alignment.getSiteCount();
        for(int i = 0; i < siteCount; i++){
            //System.err.println("substModelIndices[i]: "+dpNtdSiteModel.getSubstCurrCluster(i));
            //System.out.println(dpSiteModel.getCurrCluster(i)+" "+alignment.getPatternIndex(i));
            //clusterWeights[dpSiteModel.getCurrCluster(i)][alignment.getPatternIndex(i)]++;
            clusterWeights[dpVal.getCurrCategory(i)][alignment.getPatternIndex(i)]++;
        }

        for(int i = 0; i < siteModelCount;i++){
            //System.out.println("Hi!");
            //WVTreeLikelihood treeLik = new WVTreeLikelihood(clusterWeights[i]);
            /*System.out.print("cluster weights:");
            for(int j = 0; j<clusterWeights[i].length;j++){
                System.out.print(clusterWeights[i][j]+" ");

            }
            System.out.println(); */
            NewWVTreeLikelihood treeLik = new NewWVTreeLikelihood(
                    clusterWeights[i],
                    alignment,
                    (Tree) treeInput.get(),
                    useAmbiguitiesInput.get(),
                    dpSiteModel.getSiteModel(i),
                    branchRateModelInput.get());
            /*treeLik.initByName(
                    "data", alignment,
                    "tree", treeInput.get(),
                    "siteModel", dpSiteModel.getSiteModel(i),
                    "branchRateModel", branchRateModelInput.get(),
                    "useAmbiguities",useAmbiguitiesInput.get()
            );*/
            treeLiks.add(treeLik);

            
        }


    }

    public int getDimension(){
        return treeLiks.size();
    }

    /*@Override
    public double calculateLogP() throws Exception{
        logP = 0.0;
        //System.out.println("hello: "+treeLiks.size());
        //double sum = 0.0;
        //for(WVTreeLikelihood treeLik : treeLiks) {
        for(NewWVTreeLikelihood treeLik : treeLiks) {
            //sum+=treeLik.weightSum();
            //System.err.println("hello?"+treeLik.isDirtyCalculation());
            //System.out.println("siteModel id: "+ ((SwitchingNtdBMA)treeLik.m_substitutionModel).getIDNumber());
        	if (treeLik.isDirtyCalculation()) {
                //System.err.println("hello?");
                //logP += treeLik.calculateLogP();
                double tmp = treeLik.calculateLogP();
        		logP += tmp;
                //System.out.println("calcLogP: "+tmp);
                //System.out.println("calcLogP: "+tmp+treeLik.m_data.get()+" "+((SiteModel)treeLik.m_siteModel).getRateParameter());
        	} else {
        		logP += treeLik.getCurrentLogP();
                //System.out.println("currLogP: "+treeLik.calculateLogP());
        	}
            if (Double.isInfinite(logP) || Double.isNaN(logP)) {
            	return logP;
            }
        }

        //if(sum != alignment.getSiteCount()){
        //    throw new RuntimeException("WAY WRONG");
        //}
        //System.out.println("logP: "+logP);
        return logP;
    }  */


    @Override
    public double calculateLogP() {

        logP = 0.0;

        int nrOfDirtyDistrs = 0;
        /*for (Distribution dists : treeLiks) {
                if (dists.isDirtyCalculation()) {
                    nrOfDirtyDistrs++;
                }
            }
        System.out.println("nrOfDirtyDistrs:"+nrOfDirtyDistrs); */
        if(useThreadsEvenly){
            logP = calculateLogPUsingThreadsEvenly();

        }if (useThreads) {

            logP = calculateLogPUsingThreads();
        }else{
            //int counter = 0;
            double sum = 0.0;
            for(NewWVTreeLikelihood treeLik : treeLiks) {

               //System.out.println("start: ");

            	if (treeLik.isDirtyCalculation()) {
                    double tmp = treeLik.calculateLogP();
            		logP += tmp;
                    //System.out.println("calculateLogP: "+tmp);
                    //counter++;
            	} else {
            		logP += treeLik.getCurrentLogP();
                    //System.out.println("CurrentLogP: "+treeLik.getCurrentLogP());
            	}

                if (Double.isInfinite(logP) || Double.isNaN(logP)) {
                	return logP;
                }
                sum+=treeLik.weightSum();

            }

            if(sum != alignment.getSiteCount()){
                System.out.println("Sum: "+sum+" "+alignment.getSiteCount());
                throw new RuntimeException("Sum wrong");
            }

            //System.out.println(counter++);

        }

    //System.out.println("End Compute Likelihood");
        /*if(changeType == ChangeType.SPLIT){
            System.out.println("p: "+logP+" "+storedLogP);
        }*/

        //System.out.println(treeLiks.size()+" p: "+logP+" "+storedLogP);

        //System.out.println(getID()+": "+logP);
        //treeLiks.get(0).getTree().log(-1, System.out);
        return logP;

    }

    CountDownLatch m_nCountDown;
    private double calculateLogPUsingThreads() {
        try {
            //System.out.println("What?");
            int nrOfDirtyDistrs = 0;
            for (Distribution dists : treeLiks) {
                if (dists.isDirtyCalculation()) {
                    nrOfDirtyDistrs++;
                }
            }
            //System.out.println("nrOfDirtyDistrs: "+nrOfDirtyDistrs);
            m_nCountDown = new CountDownLatch(nrOfDirtyDistrs);
            // kick off the threads
            for (Distribution dists : treeLiks) {
                if (dists.isDirtyCalculation()) {
                    //System.out.println("isDirty");
                    CoreRunnable coreRunnable = new CoreRunnable(dists);
                    BeastMCMC.g_exec.execute(coreRunnable);
                }
            }
            m_nCountDown.await();
            logP = 0;
            for (Distribution distr : treeLiks) {
                logP += distr.getCurrentLogP();
            }
            return logP;
        } catch (RejectedExecutionException e) {
            useThreads = false;
            //System.err.println("Stop using threads: " + e.getMessage());
            // refresh thread pool
            BeastMCMC.g_exec = Executors.newFixedThreadPool(BeastMCMC.m_nThreads);
            return calculateLogP();
        } catch (InterruptedException e) {
            e.printStackTrace();
            throw new RuntimeException("Interrupted.");
        }
    }


    private double calculateLogPUsingThreadsEvenly() {

        try {

            //int nrOfDirtyDistrs = 0;
            ArrayList<Distribution> tempDists = new ArrayList<>();
            for (Distribution dists : treeLiks) {
                if (dists.isDirtyCalculation()) {
                    tempDists.add(dists);
                    //nrOfDirtyDistrs++;
                }
            }


            if(tempDists.size() < BeastMCMC.m_nThreads){
                m_nCountDown = new CountDownLatch(tempDists.size());
                for(Distribution distr:tempDists){
                    CoreRunnable coreRunnable = new CoreRunnable(distr);
                    BeastMCMC.g_exec.execute(coreRunnable);
                }
            }else{
            //m_nCountDown = new CountDownLatch(nrOfDirtyDistrs);
                m_nCountDown = new CountDownLatch(BeastMCMC.m_nThreads);
                // kick off the threads
                int batchSize = tempDists.size()/ BeastMCMC.m_nThreads;

                int start, end;
                for (int i = 0; i < BeastMCMC.m_nThreads; i++) {
                    //if (dists.isDirtyCalculation()) {    //todo
                    start = batchSize*i;
                    end = batchSize*(i + 1)-1;
                    if(i == BeastMCMC.m_nThreads -1){
                        end += tempDists.size()% BeastMCMC.m_nThreads;
                    }

                    EvenlyCoreRunnable coreRunnable = new EvenlyCoreRunnable(start, end, tempDists);

                    BeastMCMC.g_exec.execute(coreRunnable);
                    //}

                }
            }

            m_nCountDown.await();

            logP = 0;
            for (Distribution distr : treeLiks) {
                logP += distr.getCurrentLogP();
            }
            return logP;
        } catch (RejectedExecutionException e) {
            useThreadsEvenly = false;
            System.err.println("Stop using threads: " + e.getMessage());
            // refresh thread pool
            BeastMCMC.g_exec = Executors.newFixedThreadPool(BeastMCMC.m_nThreads);
            return calculateLogP();
        } catch (InterruptedException e) {
            e.printStackTrace();
            throw new RuntimeException("Interrupted.");
        }
    }

    class EvenlyCoreRunnable implements Runnable {
        int start, end;
        ArrayList<Distribution> dists;

        EvenlyCoreRunnable(int start, int end, ArrayList<Distribution> dists) {
            this.start = start;
            this.end = end;
            this.dists = dists;
        }

        public void run() {
            Distribution distr = null;
            try {
                System.out.println("start: "+start+" end: "+end);
                for(int i = start; i < end; i++){
                    distr = dists.get(i);
                    if (distr.isDirtyCalculation()) {

                        //logP += distr.calculateLogP();
                        distr.calculateLogP();
                    } else {
                        //logP += distr.getCurrentLogP();

                        distr.getCurrentLogP();
                    }

                }

            } catch (Exception e) {
                System.err.println("Something went wrong in a calculation of " + distr.getID());
                e.printStackTrace();
                System.exit(0);
            }
            m_nCountDown.countDown();
        }

    }

    class CoreRunnable implements Runnable {
        Distribution distr;

        CoreRunnable(Distribution core) {
            distr = core;
        }

        public void run() {

            try {
                //long startTime = System.currentTimeMillis();
                if (distr.isDirtyCalculation()) {
                    //logP += distr.calculateLogP();
                    distr.calculateLogP();
                } else {
                    //logP += distr.getCurrentLogP();
                    distr.getCurrentLogP();
                }
                //long endTime = System.currentTimeMillis();

                //System.out.println("That took " + (endTime - startTime) + " milliseconds");
            } catch (Exception e) {
                System.err.println("Something went wrong in a calculation of " + distr.getID());
                e.printStackTrace();
                System.exit(0);
            }
            m_nCountDown.countDown();
        }

    } // CoreRunnable

    public double getSiteLogLikelihood(int iCluster, int iSite){
        return treeLiks.get(iCluster).getPatternLogLikelihood(alignment.getPatternIndex(iSite));
    }


    public void store(){
        storedTreeLiks = new ArrayList<NewWVTreeLikelihood>();
        for(NewWVTreeLikelihood treeLik : treeLiks) {
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
        for(NewWVTreeLikelihood treeLik : treeLiks) {
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
        SiteModel siteModel = dpSiteModel.getSiteModel(dpSiteModel.getLastAddedIndex());

        int[] patternWeights = new int[alignment.getPatternCount()];
        patternWeights[alignment.getPatternIndex(dpSiteModel.getLastDirtySite())] +=1;

        //OldWVAlignment wvalign = new OldWVAlignment(alignment, patternWeights);
        //WVTreeLikelihood treeLik = new WVTreeLikelihood(patternWeights);
        NewWVTreeLikelihood treeLik = new NewWVTreeLikelihood(
                patternWeights,
                alignment,
                (Tree) treeInput.get(),
                useAmbiguitiesInput.get(),
                siteModel,
                branchRateModelInput.get());
        try{
            /*treeLik.initByName(
                    "data", alignment,
                    "tree", treeInput.get(),
                    "siteModel", siteModel,
                    "branchRateModel", branchRateModelInput.get(),
                    "useAmbiguities",useAmbiguitiesInput.get()
            );*/

            treeLik.calculateLogP();
            treeLik.store();
            treeLiks.add(dpSiteModel.getLastAddedIndex(),treeLik);
        }catch(Exception e){
            throw new RuntimeException(e);
        }

    }

    public void splitTreeLikelihood(){

        /*for(int i = 0; i < treeLiks.size(); i++){
            int[] weights = treeLiks.get(i).getPatternWeights();
            for(int j = 0; j < weights.length; j++){
                System.out.print(weights[j]+" ");
            }
            System.out.println();
        }*/

        //alignment.printPattern();
        SiteModel siteModel = dpSiteModel.getSiteModel(dpSiteModel.getLastAddedIndex());

        int[] clusterSites = dpVal.getClusterSites(dpSiteModel.getLastAddedIndex());

        //System.out.println("Last: "+dpSiteModel.getLastAddedIndex());

        int[] patternWeights = new int[alignment.getPatternCount()];

        int prevCluster = dpSiteModel.getDirtySiteModelIndex();
        //System.out.println("prevCluster: "+prevCluster);
        NewWVTreeLikelihood prevTreeLikelihood = treeLiks.get(prevCluster);
        for(int i = 0; i < clusterSites.length;i++){

            int patternIndex = alignment.getPatternIndex(clusterSites[i]);
            //System.out.println("clusterSites: "+clusterSites[i]+" PatIndex: "+patternIndex);
            patternWeights[patternIndex]++;
            //System.out.println(i+" splitting: "+patternWeights[patternIndex]);
            prevTreeLikelihood.removeWeight(patternIndex,1);
        }
        //OldWVAlignment wvalign = new OldWVAlignment(alignment, patternWeights);
        //WVTreeLikelihood treeLik = new WVTreeLikelihood(patternWeights);
        NewWVTreeLikelihood treeLik = new NewWVTreeLikelihood(
                patternWeights,
                alignment,
                (Tree) treeInput.get(),
                useAmbiguitiesInput.get(),
                siteModel,
                branchRateModelInput.get());
        try{
            /*treeLik.initByName(
                    "data", alignment,
                    "tree", treeInput.get(),
                    "siteModel", siteModel,
                    "branchRateModel", branchRateModelInput.get(),
                    "useAmbiguities",useAmbiguitiesInput.get()
            );*/
            treeLik.calculateLogP();
            treeLik.store();
            treeLiks.add(dpSiteModel.getLastAddedIndex(),treeLik);
        }catch(Exception e){
            throw new RuntimeException(e);
        }



    }


    public void mergeTreeLikelihoods(int removedIndex){

        NewWVTreeLikelihood removedTreeLikelihood =  treeLiks.remove(removedIndex);
        NewWVTreeLikelihood mergedTreeLikelihood =  treeLiks.get(dpSiteModel.getDirtySiteModelIndex());
        int [] patternWeights = removedTreeLikelihood.getPatternWeights();
        for(int i = 0;i < patternWeights.length;i++){
            mergedTreeLikelihood.addWeight(i,patternWeights[i]);
        }
    }


    /*public void mergeTreeLikelihoods(int removedIndex, boolean valueChanged){

        NewWVTreeLikelihood removedTreeLikelihood =  treeLiks.remove(removedIndex);
        NewWVTreeLikelihood mergedTreeLikelihood =  treeLiks.get(dpSiteModel.getDirtySiteModelIndex());
        if(valueChanged){
            int[] initPatternWeights = mergedTreeLikelihood.getPatternWeights();
            for(int i = 0;i < initPatternWeights.length;i++){
                if(initPatternWeights[i] > 0){
                    mergedTreeLikelihood.addZeroWeight(i);
                }
            }
        }

        int [] patternWeights = removedTreeLikelihood.getPatternWeights();
        for(int i = 0;i < patternWeights.length;i++){
            mergedTreeLikelihood.addWeight(i,patternWeights[i]);
        }
    }*/




    /*public void mergeTreeLikelihoods(int removedIndex, boolean valueChanged){
        NewWVTreeLikelihood mergedTreeLikelihood =  treeLiks.get(dpSiteModel.getDirtySiteModelIndex());
        if(valueChanged){
            int[] initPatternWeights = mergedTreeLikelihood.getPatternWeights();
            for(int i = 0;i < initPatternWeights.length;i++){
                if(initPatternWeights[i] > 0){
                    mergedTreeLikelihood.addZeroWeight(i);
                }
            }
        }

        NewWVTreeLikelihood removedTreeLikelihood =  treeLiks.remove(removedIndex);
        int [] patternWeights = removedTreeLikelihood.getPatternWeights();
        for(int i = 0;i < patternWeights.length;i++){
            mergedTreeLikelihood.addWeight(i,patternWeights[i]);
        }
    } */

    public void removeTreeLikelihood(int removedIndex){
        treeLiks.remove(removedIndex);
    }

    private void handlePointerChange(int dirtySite){
        int prevCluster = dpVal.getPrevCategory(dirtySite);
        int currCluster = dpVal.getCurrCategory(dirtySite);
        //System.out.println("dirtySite: "+dirtySite+" prevCluster: "+prevCluster+" currCluster"+currCluster);

        //System.out.println(dirtySite+" "+prevCluster+" "+currCluster);
        //System.out.println(treeLiks.get(2).getPatternWeights()[alignment.getPatternIndex(307)]);

        treeLiks.get(prevCluster).removeWeight(alignment.getPatternIndex(dirtySite),1);


        treeLiks.get(currCluster).addWeight(alignment.getPatternIndex(dirtySite),1);


    }

    private void handlePointersChange(){
        int[] dirtySites = dpSiteModel.getLastDirtySites();
        for(int dirtySite: dirtySites){
            //System.out.println(dirtySite);
            handlePointerChange(dirtySite);
        }


    }

    protected void updateWeights(){
        int dirtySite = dpSiteModel.getLastDirtySite();
        //System.out.println("changeType: "+changeType);
        if(changeType == ChangeType.ADDED){
            int prevCluster = dpSiteModel.getPrevCluster(dirtySite);
            int currCluster = dpSiteModel.getCurrCluster(dirtySite);
            //System.out.println(currCluster+" "+prevCluster);
            treeLiks.get(prevCluster).removeWeight(alignment.getPatternIndex(dirtySite),1);
        }else if(changeType==ChangeType.POINTER_CHANGED){
            int prevCluster = dpSiteModel.getPrevCluster(dirtySite);
            treeLiks.get(prevCluster).removeWeight(alignment.getPatternIndex(dirtySite),1);

            int currCluster = dpSiteModel.getCurrCluster(dirtySite);
            treeLiks.get(currCluster).addWeight(alignment.getPatternIndex(dirtySite),1);
        }else if(changeType == ChangeType.POINTERS_SWAPPED){
            int[] swappedSites = dpSiteModel.getSwappedSites();

            for(int i = 0; i < swappedSites.length; i++){
                dirtySite = swappedSites[i];
                int prevCluster = dpSiteModel.getPrevCluster(dirtySite);
                treeLiks.get(prevCluster).removeWeight(alignment.getPatternIndex(dirtySite),1);

                int currCluster = dpSiteModel.getCurrCluster(dirtySite);
                treeLiks.get(currCluster).addWeight(alignment.getPatternIndex(dirtySite),1);
            }

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

            }else if(changeType == ChangeType.SPLIT_AND_VALUE_CHANGE){

                splitTreeLikelihood();
                MCMCNodeFactory.checkDirtiness(treeLiks.get(dpSiteModel.getDirtySiteModelIndex()));
                //System.out.println("changed:"+dpSiteModel.getDirtySiteModelIndex());
                this.changeType = ChangeType.SPLIT_AND_VALUE_CHANGE;

            }else if(changeType == ChangeType.MERGE){

                mergeTreeLikelihoods(dpSiteModel.getRemovedIndex());
                this.changeType = ChangeType.MERGE;

            }else if(changeType == ChangeType.MERGE_AND_VALUE_CHANGE){

                mergeTreeLikelihoods(dpSiteModel.getRemovedIndex());
                MCMCNodeFactory.checkDirtiness(treeLiks.get(dpSiteModel.getDirtySiteModelIndex()));
                this.changeType = ChangeType.MERGE_AND_VALUE_CHANGE;

            }else if (changeType == ChangeType.POINTER_CHANGED||  changeType == ChangeType.POINTERS_SWAPPED){
                this.changeType = changeType;
                updateWeights();
            }else if(changeType == ChangeType.MULTIPLE_POINTERS_CHANGED){
                handlePointersChange();
                this.changeType = changeType;
            }else if(changeType == ChangeType.VALUE_CHANGED){
                this.changeType = ChangeType.VALUE_CHANGED;

            }else{
                this.changeType = ChangeType.ALL;
            }


            
            recalculate = true;
        }else if(treeInput.get().somethingIsDirty()){
            recalculate = true;

        }else if(branchRateModelInput.get().isDirtyCalculation()){
            recalculate = true;
        }
        if(recalculate){

            for(NewWVTreeLikelihood treeLik:treeLiks){

                MCMCNodeFactory.checkDirtiness(treeLik);
            }
        }
        //System.out.println("changeType: "+changeType);
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
        return (Tree) treeInput.get();
    }

    public int[] getClusterWeights(int clusterIndex){
        return Arrays.copyOf(treeLiks.get(clusterIndex).getPatternWeights(), alignment.getPatternCount());

    }



    protected void accept() {
        //for(WVTreeLikelihood treeLik:treeLiks){
        for(NewWVTreeLikelihood treeLik:treeLiks){
            treeLik.makeAccept();
        }
    	super.accept();
    }
}
