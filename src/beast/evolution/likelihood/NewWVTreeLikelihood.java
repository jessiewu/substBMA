package beast.evolution.likelihood;

import beast.core.Description;
import beast.core.parameter.QuietRealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.AscertainedAlignment;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.branchratemodel.StrictClockModel;
//import beast.evolution.sitemodel.QuietGammaSiteBMA;
import beast.evolution.sitemodel.QuietGammaSiteBMA;
import beast.evolution.sitemodel.QuietSiteModel;
import beast.evolution.sitemodel.SiteModelInterface;
import beast.evolution.substitutionmodel.NtdBMA;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.substitutionmodel.SwitchingNtdBMA;
import beast.evolution.sitemodel.SiteModel;

import java.util.Arrays;

/**
 * @author Chieh-Hsi Wu
 */
@Description("Tree likelihood supports alignment pattern weights to change during the MCMC.")
public class NewWVTreeLikelihood extends QuietTreeLikelihood {
    WVLikelihoodCore m_likelihoodCore;
   
    protected int[] patternWeights;
    protected int[] storedPatternWeights;
    protected boolean weightsChanged;
    //protected int weightSum = 0;
    //protected int storedWeightSum = 0;
    //protected boolean nonZeroPatternIncreased;
    protected double[] storedPatternLogLikelihoods;

    public NewWVTreeLikelihood(){

    }

    public NewWVTreeLikelihood(int[] patternWeights){
        this.patternWeights = patternWeights;
        storedPatternWeights = new int[patternWeights.length];

    }

    public NewWVTreeLikelihood(int[] patternWeights,
                               Alignment data,
                               Tree tree,
                               boolean useAmbiguities,
                               SiteModel siteModel,
                               BranchRateModel.Base branchRateModel){
        this.patternWeights = patternWeights;
        storedPatternWeights = new int[patternWeights.length];
        this.data = data;
        this.tree = tree;
        this.useAmbiguities = useAmbiguities;
        m_siteModel = siteModel;
        this.branchRateModel = branchRateModel;
        this.substitutionModel = (SubstitutionModel.Base)m_siteModel.getSubstitutionModel();
        setup();
    }

    public String getSiteModelID(){
        return m_siteModel.getID();
    }

    public void initAndValidate() throws Exception{

        super.initAndValidate();
        //setup();

    }





    @Override
    protected void setup() {
        int nodeCount = tree.getNodeCount();
        //m_siteModel = m_pSiteModel.get();
        m_siteModel.setDataType(data.getDataType());
        //m_substitutionModel = m_siteModel.m_pSubstModel.get();
        //m_substitutionModel = ((QuietSiteModel)m_siteModel).getSubstitutionModel();

        /*if (m_pBranchRateModel.get() != null) {
        	m_branchRateModel = m_pBranchRateModel.get();
        } else {
            m_branchRateModel = new StrictClockModel();
        }*/
    	m_branchLengths = new double[nodeCount];
    	storedBranchLengths = new double[nodeCount];

        int nStateCount = data.getMaxStateCount();
        //System.out.println("nStateCount: "+nStateCount);
        int nPatterns = data.getPatternCount();
        //System.out.println("patternWeights.length:"+patternWeights.length);

        boolean[] unmasked = new boolean[patternWeights.length];
        for(int i = 0; i < unmasked.length;i++){
                ///System.out.println(patternWeights[i]);
                unmasked[i] = patternWeights[i] > 0;
            }

        if (nStateCount == 4) {

            m_likelihoodCore = new WVLikelihoodCore4(unmasked);
        } else {
            m_likelihoodCore = new WVLikelihoodCore(nStateCount,unmasked);
        }

        //System.out.println("m_likelihoodCore: "+(m_likelihoodCore==null));
        //System.err.println("TreeLikelihood uses " + m_likelihoodCore.getClass().getName());

        proportionInvariant = m_siteModel.getProportionInvariant();
        m_siteModel.setPropInvariantIsCategory(false);
        if(m_siteModel instanceof QuietGammaSiteBMA){
            calcConstantPatternIndices(nPatterns, nStateCount);

        }else if (proportionInvariant > 0) {
        	calcConstantPatternIndices(nPatterns, nStateCount);
        }
        addedPatternIds = new int[nPatterns];
        initCore();

        patternLogLikelihoods = new double[nPatterns];
        storedPatternLogLikelihoods = new double[nPatterns];
        m_fRootPartials = new double[nPatterns * nStateCount];
        matrixSize = (nStateCount +1)* (nStateCount+1);
        probabilities = new double[(nStateCount +1)* (nStateCount+1)];
        //System.out.println("m_fProbabilities: "+m_fProbabilities.length+" "+m_data.get());
        Arrays.fill(probabilities, 1.0);

        if (data instanceof AscertainedAlignment) {
            useAscertainedSitePatterns = true;
        }
    }



    void initCore() {
        final int nodeCount = tree.getNodeCount();
        m_likelihoodCore.initialize(
                nodeCount,
                data.getPatternCount(),
                m_siteModel.getCategoryCount(),
                true, useAmbiguities
        );

        final int extNodeCount = nodeCount / 2 + 1;
        final int intNodeCount = nodeCount / 2;
        //System.out.println("m_likelihoodCore: "+(m_likelihoodCore == null));
        if (useAmbiguities) {
            setPartials(tree.getRoot(), data.getPatternCount());
        } else {
            setStates(tree.getRoot(), data.getPatternCount());
        }
        hasDirt = Tree.IS_FILTHY;
        for (int i = 0; i < intNodeCount; i++) {
            m_likelihoodCore.createNodePartials(extNodeCount + i);
        }
    }

    /**
     * set leaf partials in likelihood core *
     */
    void setPartials(Node node, int patternCount) {
        if (node.isLeaf()) {
            //Alignment data = m_data.get();
            int nStates = data.getDataType().getStateCount();
            double[] partials = new double[patternCount * nStates];

            int k = 0;
            int iTaxon = data.getTaxonIndex(node.getID());
            for (int iPattern = 0; iPattern < patternCount; iPattern++) {
                int nState = data.getPattern(iTaxon, iPattern);
                boolean[] stateSet = data.getStateSet(nState);
                for (int iState = 0; iState < nStates; iState++) {
                    partials[k++] = (stateSet[iState] ? 1.0 : 0.0);
                }
            }
            //System.out.println(m_likelihoodCore == null);
            m_likelihoodCore.setNodePartials(node.getNr(), partials);

        } else {
            //System.out.println(m_likelihoodCore == null);
            setPartials(node.getLeft(), patternCount);
            setPartials(node.getRight(), patternCount);
        }
    }

    /** set leaf states in likelihood core **/
    /**
     * set leaf states in likelihood core *
     */
    void setStates(Node node, int patternCount) {
        if (node.isLeaf()) {
            int i;
            int[] states = new int[patternCount];
            int iTaxon = data.getTaxonIndex(node.getID());
            for (i = 0; i < patternCount; i++) {
                states[i] = data.getPattern(iTaxon, i);
            }
            m_likelihoodCore.setNodeStates(node.getNr(), states);

        } else {
            setStates(node.getLeft(), patternCount);
            setStates(node.getRight(), patternCount);
        }
    }


    boolean valueChanged = false;
    @Override
    public double calculateLogP() {

        //((NtdBMA)m_substitutionModel).printDetails();
        boolean[] trueUnmasked = m_likelihoodCore.getUnmasked();
        //System.out.println("m_nHasDirt: "+m_nHasDirt);
        //System.out.println("valueChanged: "+valueChanged);
        if(hasDirt > -1){
            //System.out.println("addedPatternIdCount: "+addedPatternIdCount);
            if(addedPatternIdCount > 0){
                boolean[] tempUnmasked = new boolean[patternWeights.length];
                //System.out.println(addedPatternIdCount);
                for(int i = 0; i < addedPatternIdCount; i++){
                    if(!valueChanged){
                        tempUnmasked[addedPatternIds[i]] = true;
                    }
                    trueUnmasked[addedPatternIds[i]] = true;
                    //System.out.println("addedPatternId: "+addedPatternIds[i]);
                }
                //m_likelihoodCore.setUnmasked(trueUnmasked);
                if(valueChanged){
                    m_likelihoodCore.setUnmasked(trueUnmasked);

                }else{
                    m_likelihoodCore.setUnmasked(tempUnmasked);
                }
            }
            /*for(int i = 0; i < trueUnmasked.length;i++){
                trueUnmasked[i] = patternWeights[i] > 0;
            } */

            /*boolean[] temp = m_likelihoodCore.getUnmasked();
            for(int i = 0; i < temp.length; i++){
                System.out.print(temp[i]+" ");
            }
            System.out.println(); */
            for(int i = 0; i < patternWeights.length; i++){
                if(patternWeights[i] == 0 && trueUnmasked[i]){
                    System.out.println(addedPatternIdCount);
                    for(int j = 0; j < addedPatternIdCount;j++){
                        System.out.println(addedPatternIds[j]+" "+patternWeights[addedPatternIds[j]]);

                    }
                    throw new RuntimeException("what? "+i);

                }
            }
            /*for(int i = 0; i < patternWeights.length;i++){
                System.out.println(i+" "+m_likelihoodCore.getUnmasked(i));
            } */
       	    traverse(tree.getRoot());
        }

        /*boolean[] U = m_likelihoodCore.getUnmasked();
        System.out.print("Unmasked: ");
        for(boolean u:U){
            System.out.print(u+" ");
        }
        System.out.println(); */

        calcLogP();


        m_nScale++;
        if (logP > 0 || (m_likelihoodCore.getUseScaling() && m_nScale > X)) {
            //System.out.println("Switch off scaling");
            m_likelihoodCore.setUseScaling(1.0);
            m_likelihoodCore.unstore();
            hasDirt = Tree.IS_FILTHY;
            X *= 2;
           	traverse(tree.getRoot());
            calcLogP();
        } else if (logP == Double.NEGATIVE_INFINITY && m_fScale < 10) { // && !m_likelihoodCore.getUseScaling()) {
        	m_nScale = 0;
        	m_fScale *= 1.01;
            //System.out.println("Turning on scaling to prevent numeric instability " + m_fScale);
            m_likelihoodCore.setUseScaling(m_fScale);
            m_likelihoodCore.unstore();
            hasDirt = Tree.IS_FILTHY;
           	traverse(tree.getRoot());
            calcLogP();
        }

        if(addedPatternIdCount > 0 && !valueChanged){
            m_likelihoodCore.setUnmasked(trueUnmasked);
        }


        return logP;
    }

    protected void calcLogP() {
        /*if(m_siteModel instanceof  QuietSiteModel){
            printThings();
        System.out.println("mu: "+((QuietSiteModel)m_siteModel).getRateParameter()+" "+((QuietSiteModel)m_siteModel).getRateParameter().getIDNumber());
        }*/
        logP = 0.0;
        if (useAscertainedSitePatterns) {
            double ascertainmentCorrection = ((AscertainedAlignment)data).getAscertainmentCorrection(patternLogLikelihoods);
            for (int i = 0; i < data.getPatternCount(); i++) {
            	logP += (patternLogLikelihoods[i] - ascertainmentCorrection) * patternWeights[i];

            }
        } else {
            //System.err.println("Likelihood: ");
	        for (int i = 0; i < data.getPatternCount(); i++) {
                //System.out.println("logKappa: "+(((SwitchingNtdBMA)m_substitutionModel).getLogKappa()).getValue());
                //System.out.println((((SwitchingNtdBMA)m_substitutionModel).getIDNumber()));
                /*if(patternWeights[i]>0 && m_fPatternLogLikelihoods[i]==0){
                    ((SwitchingNtdBMA)m_substitutionModel).printDetails();
                    int pattern[] = m_data.get().getPattern(i);
                    for(int j = 0;j < pattern.length;j++){
                        System.out.print(pattern[j]+" ");
                    }
                    System.out.println();
                    throw new RuntimeException("WRONG!");
                }*/
                //System.out.println(m_data.get().getID()+" pattern "+i+" Likelihood: "+m_fPatternLogLikelihoods[i]+" "+patternWeights[i]);
	            logP += patternLogLikelihoods[i] * patternWeights[i];
                //System.out.println(m_fPatternLogLikelihoods[i] +" "+ patternWeights[i]);
                //printThings();
                //System.out.println("m_fPatternLogLikelihoods "+i+": "+m_fPatternLogLikelihoods[i]+" "+ patternWeights[i]);
	        }

            //printThings();

        }
    }

    /* Assumes there IS a branch rate model as opposed to traverse() */
    int traverse(Node node) {

        int update = (node.isDirty()| hasDirt);
        //System.out.println("update node: "+update);
        int iNode = node.getNr();

        double branchRate = branchRateModel.getRateForBranch(node);
        double branchTime = node.getLength() * branchRate;
        m_branchLengths[iNode] = branchTime;

        // First update the transition probability matrix(ices) for this branch
        if (!node.isRoot() && (update != Tree.IS_CLEAN || branchTime != storedBranchLengths[iNode])) {
            Node parent = node.getParent();
            m_likelihoodCore.setNodeMatrixForUpdate(iNode);
            for (int i = 0; i < m_siteModel.getCategoryCount(); i++) {
                double jointBranchRate = m_siteModel.getRateForCategory(i, node) * branchRate;
                //System.out.println(getID()+" mu: "+m_siteModel.getRateForCategory(i, node)+" "+branchRate);
                //System.out.println(m_data.get()+ " update node: "+jointBranchRate);
                substitutionModel.getTransitionProbabilities(node, parent.getHeight(), node.getHeight(), jointBranchRate, probabilities);
                m_likelihoodCore.setNodeMatrix(iNode, i, probabilities);
            }
            update |= Tree.IS_DIRTY;
        }

        // If the node is internal, update the partial likelihoods.
        if (!node.isLeaf()) {

            // Traverse down the two child nodes
            Node child1 = node.getLeft(); //Two children
            int update1 = traverse(child1);

            Node child2 = node.getRight();
            int update2 = traverse(child2);

            // If either child node was updated then update this node too
            if (update1 != Tree.IS_CLEAN || update2 != Tree.IS_CLEAN) {

                int childNum1 = child1.getNr();
                int childNum2 = child2.getNr();
                if(addedPatternIdCount == 0 || valueChanged){

                    m_likelihoodCore.setNodePartialsForUpdate(iNode);
                }
                update |= (update1|update2);
                if (update >= Tree.IS_FILTHY) {
                    m_likelihoodCore.setNodeStatesForUpdate(iNode);
                }

                if (m_siteModel.integrateAcrossCategories()) {

                    m_likelihoodCore.calculatePartials(childNum1, childNum2, iNode);
                } else {
                    throw new RuntimeException("Error TreeLikelihood 201: Site categories not supported");
                    //m_pLikelihoodCore->calculatePartials(childNum1, childNum2, nodeNum, siteCategories);
                }

                if (node.isRoot()) {
                    // No parent this is the root of the beast.tree -
                    // calculate the pattern likelihoods
                    double[] frequencies = //m_pFreqs.get().
                            substitutionModel.getFrequencies();

                    double[] proportions = m_siteModel.getCategoryProportions(node);
                    m_likelihoodCore.integratePartials(node.getNr(), proportions, m_fRootPartials);
                    /*System.out.println(getID()+":");
                    for(int i = 0;i < m_fRootPartials.length;i++){
                        System.out.println("m_fRootPartials: "+m_fRootPartials[i]);
                    }
                    System.out.println();*/

                    //System.out.println(m_iConstantPattern == null);
                    if (constantPattern != null) { // && !SiteModel.g_bUseOriginal) {
                    	proportionInvariant = m_siteModel.getProportionInvariant();
                        //System.out.println("m_fProportionInvariant:"+m_fProportionInvariant);

                    	// some portion of sites is invariant, so adjust root partials for this
                    	for (int i : constantPattern) {
                            //System.out.println(i+" "+m_fRootPartials[i]+" "+m_fProportionInvariant);
                			m_fRootPartials[i] += proportionInvariant;
                    	}
                    }

                    m_likelihoodCore.calculateLogLikelihoods(m_fRootPartials, frequencies, patternLogLikelihoods);
                }

            }
        }
        return update;
    }

    public double getPatternLogLikelihood(int iPat){
        //System.out.println(iPat+" "+m_fPatternLogLikelihoods[iPat]);
        if(patternWeights[iPat] == 0){
            //throw new RuntimeException("Zero weight.");
            return Double.NaN;
        }

        return patternLogLikelihoods[iPat];
    }
    /**
     * check state for changed variables and update temp results if necessary *
     */
    @Override
    protected boolean requiresRecalculation() {
        //printThings();

        hasDirt = Tree.IS_CLEAN;

        if(weightsChanged){

            if(addedPatternIdCount > 0){
                hasDirt = Tree.IS_FILTHY;
            }else{
                hasDirt = -1;
            }

            if (m_siteModel.isDirtyCalculation()) {
                 hasDirt = Tree.IS_DIRTY;
                valueChanged = true;
            }
            //System.out.println("flag2");
            return true;
        }

        if (branchRateModel != null && branchRateModel.isDirtyCalculation()) {
            hasDirt = Tree.IS_FILTHY;
            //System.out.println("flag3");
            return true;
        }

        if (data.isDirtyCalculation()) {
            hasDirt = Tree.IS_FILTHY;
            //System.out.println("flag4");
            return true;
        }

        if (m_siteModel.isDirtyCalculation()) {
            valueChanged = true;
            hasDirt = Tree.IS_DIRTY;
            //System.out.println("flag5");
            return true;
        }

        return tree.somethingIsDirty();
    }

    int[] addedPatternIds;
    int addedPatternIdCount = 0;
    public void addWeight(int patId, int dweight){
        if(patternWeights[patId] == 0 && dweight > 0){
            //System.out.println("patId: "+patId);
            addedPatternIds[addedPatternIdCount++] = patId;
        }
        patternWeights[patId] += dweight;
        //System.err.println("patId: "+patId+", patternWeights[patId]: "+patternWeights[patId]);
        weightsChanged = true;
    }



    public void removeWeight(int patId, int dweight){
        //System.out.println("patId: "+patId+", patternWeights[patId]: "+patternWeights[patId]);
        patternWeights[patId] -= dweight;

        if(patternWeights[patId] == 0){
            m_likelihoodCore.setUnmasked(patId,false);
        }else if(patternWeights[patId] < 0){
            /*for(int i = 0; i < patternWeights.length;i++){
                System.out.print(patternWeights[i]+" ");
            }
            System.out.println();  */
            throw new RuntimeException("NEGATIVE WEIGHTS: "+patternWeights[patId]+" "+dweight);
        }

        weightsChanged = true;
    }

    public void setWeight(int patId, int weight){
        patternWeights[patId] = weight;
        weightsChanged = true;
    }



    @Override
    public void store(){
        m_likelihoodCore.store();
        weightsChanged = false;
        valueChanged = false;
        addedPatternIds = new int[patternWeights.length];
        addedPatternIdCount = 0;
        System.arraycopy(patternWeights,0,storedPatternWeights,0,patternWeights.length);
        System.arraycopy(patternLogLikelihoods,0,storedPatternLogLikelihoods,0,patternLogLikelihoods.length);
        super.store();
    }

    @Override
    public void restore(){
        if(addedPatternIdCount == 0 || valueChanged){
            //System.out.println("restore");
            m_likelihoodCore.restore();
        }else{
            //System.out.println("Restore?");
            m_likelihoodCore.restoreUnmaskedOnly();
        }
        int[] temp1 = patternWeights;
        patternWeights = storedPatternWeights;
        storedPatternWeights = temp1;

        double[] temp2 = patternLogLikelihoods;
        patternLogLikelihoods = storedPatternLogLikelihoods;
        storedPatternLogLikelihoods = temp2;

        weightsChanged = false;
        valueChanged = false;
        addedPatternIdCount = 0;
        super.restore();
    }

    public int[] getPatternWeights(){
        return patternWeights;
    }

    public void setPatternWeights(int[] newPatternWeights){
        System.arraycopy(newPatternWeights, 0, patternWeights, 0, newPatternWeights.length);
        boolean[] unmasked = new boolean[patternWeights.length];
        for(int i = 0;i < patternWeights.length;i++){
            unmasked[i] = patternWeights[i] > 0;
        }
        m_likelihoodCore.setUnmasked(unmasked);
    }

    public int weightSum(){
        int sum = 0;
        for(int weight:patternWeights){
            sum += weight;
        }
        return sum;
    }

    protected void makeAccept() {
    	super.accept();
    }

    public void printThings(){
        System.out.println("modelID: "+((SwitchingNtdBMA)substitutionModel).getIDNumber());
        ((SwitchingNtdBMA)substitutionModel).printDetails();
        System.out.println("modelID: " + ((QuietSiteModel) m_siteModel).getRateParameter().getIDNumber()+" "+((QuietSiteModel) m_siteModel).isDirtyCalculation());
        ((QuietSiteModel) m_siteModel).printDetails();
    }

    public int getWeight(int patternIndex){
        return patternWeights[patternIndex];

    }

    public SiteModel.Base getSiteModel(){
        return m_siteModel;
    }
}
