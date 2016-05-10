package beast.evolution.likelihood;

import beast.core.Description;
import beast.core.parameter.QuietRealParameter;
import beast.evolution.sitemodel.QuietSiteModel;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.tree.Tree;
import beast.evolution.alignment.AscertainedAlignment;
import beast.evolution.branchratemodel.StrictClockModel;
import beast.evolution.substitutionmodel.SwitchingNtdBMA;

import java.util.Arrays;

/**
 * @author Chieh-Hsi Wu
 */
@Description("This class is used to tree likelihood with variable pattern weights.")
public class WVTreeLikelihood extends TreeLikelihood{

    protected int[] patternWeights;
    protected int[] storedPatternWeights;
    protected boolean weightsChanged;
    protected double[] storedPatternLogLikelihoods;

    public WVTreeLikelihood(){

    }

    public WVTreeLikelihood(int[] patternWeights){
        this.patternWeights = patternWeights;
        storedPatternWeights = new int[patternWeights.length];
        //System.err.println("Number of pattern: "+patternWeights.length);

    }

@Override
    public void initAndValidate() {
    	// sanity check: alignment should have same #taxa as tree
    	if (dataInput.get().getNrTaxa() != treeInput.get().getLeafNodeCount()) {
    		throw new RuntimeException("The number of nodes in the tree does not match the number of sequences");
    	}
    	//m_beagle = null;
    	/*m_beagle = new BeagleTreeLikelihood();
    	m_beagle.initByName(
                "data", m_data.get(),
                "tree", m_tree.get(),
                "siteModel", m_pSiteModel.get(),
                "branchRateModel", m_pBranchRateModel.get(),
                "useAmbiguities", m_useAmbiguities.get());
    	if (m_beagle.beagle != null) {
    		//a Beagle instance was found, so we use it
    		return;
    	}
		// No Beagle instance was found, so we use the good old java likelihood core
    	m_beagle = null;*/

        int nodeCount = treeInput.get().getNodeCount();
        if (!(siteModelInput.get() instanceof SiteModel.Base)) {
        	throw new RuntimeException ("siteModel input should be of type SiteModel.Base");
        }
        m_siteModel = (SiteModel.Base) siteModelInput.get();
        m_siteModel.setDataType(dataInput.get().getDataType());
        substitutionModel = m_siteModel.substModelInput.get();

        if (branchRateModelInput.get() != null) {
        	branchRateModel = branchRateModelInput.get();
        } else {
            branchRateModel = new StrictClockModel();
        }
    	m_branchLengths = new double[nodeCount];
    	storedBranchLengths = new double[nodeCount];

        int nStateCount = dataInput.get().getMaxStateCount();
        int nPatterns = dataInput.get().getPatternCount();
        if (nStateCount == 4) {
            likelihoodCore = new BeerLikelihoodCore4();
        } else {
            likelihoodCore = new BeerLikelihoodCore(nStateCount);
        }
        //System.err.println("TreeLikelihood uses " + m_likelihoodCore.getClass().getName());

        proportionInvariant = m_siteModel.getProportionInvariant();
        m_siteModel.setPropInvariantIsCategory(false);
        if (proportionInvariant > 0) {
        	calcConstantPatternIndices(nPatterns, nStateCount);
        }

        initCore();

        patternLogLikelihoods = new double[nPatterns];
        storedPatternLogLikelihoods = new double[nPatterns];
        m_fRootPartials = new double[nPatterns * nStateCount];
        matrixSize = (nStateCount +1)* (nStateCount+1);
        probabilities = new double[(nStateCount +1)* (nStateCount+1)];
        Arrays.fill(probabilities, 1.0);

        if (dataInput.get() instanceof AscertainedAlignment) {
            useAscertainedSitePatterns = true;
        }

        if(patternWeights.length != dataInput.get().getPatternCount()){
            throw new RuntimeException("Pattern counts differ!");
        }
    }

    @Override
    public double calculateLogP() {
    	/*if (m_beagle != null) {
    		logP =  m_beagle.calculateLogP();
    		return logP;
    	}*/
        Tree tree = (Tree) treeInput.get();
        //System.err.println("m_nHasDirt: "+m_nHasDirt);
        if(hasDirt > -1){
       	    traverse(tree.getRoot());
        }

        calcLogP();

        m_nScale++;
        if (logP > 0 || (likelihoodCore.getUseScaling() && m_nScale > X)) {
            //System.err.println("Switch off scaling");
            likelihoodCore.setUseScaling(1.0);
            likelihoodCore.unstore();
            hasDirt = Tree.IS_FILTHY;
            X *= 2;
           	traverse(tree.getRoot());
            calcLogP();
            return logP;
        } else if (logP == Double.NEGATIVE_INFINITY && m_fScale < 10) { // && !m_likelihoodCore.getUseScaling()) {
        	m_nScale = 0;
        	m_fScale *= 1.01;
            System.err.println("Turning on scaling to prevent numeric instability " + m_fScale);
            likelihoodCore.setUseScaling(m_fScale);
            likelihoodCore.unstore();
            hasDirt = Tree.IS_FILTHY;
           	traverse(tree.getRoot());
            calcLogP();
            //System.err.println(getID()+": "+logP);
            return logP;
        }

        //System.err.println(getID()+": "+logP);
        return logP;
    }

    protected void calcLogP() {
        logP = 0.0;
        if (useAscertainedSitePatterns) {
            double ascertainmentCorrection = ((AscertainedAlignment)dataInput.get()).getAscertainmentCorrection(patternLogLikelihoods);
            for (int i = 0; i < dataInput.get().getPatternCount(); i++) {
            	logP += (patternLogLikelihoods[i] - ascertainmentCorrection) * patternWeights[i];
            }
        } else {

	        for (int i = 0; i < dataInput.get().getPatternCount(); i++) {
	            logP += patternLogLikelihoods[i] * patternWeights[i];
	        }
        }
    }


    public double getPatternLogLikelihood(int iPat){
        return patternLogLikelihoods[iPat];
    }
    /**
     * check state for changed variables and update temp results if necessary *
     */
    @Override
    protected boolean requiresRecalculation() {
        //System.err.println("flag1");
    	/*if (m_beagle != null) {
    		return m_beagle.requiresRecalculation();
    	} */
        hasDirt = Tree.IS_CLEAN;

        if(weightsChanged){
            hasDirt = -1;
            //System.err.println("flag2");
            return true;            
        }

        if (branchRateModel != null && branchRateModel.isDirtyCalculation()) {
            hasDirt = Tree.IS_FILTHY;
            //System.err.println("flag3");
            return true;
        }

        if (dataInput.get().isDirtyCalculation()) {
            hasDirt = Tree.IS_FILTHY;
            //System.err.println("flag4");
            return true;
        }

        if (m_siteModel.isDirtyCalculation()) {
            hasDirt = Tree.IS_DIRTY;
            //System.err.println("flag5");
            return true;
        }

        return treeInput.get().somethingIsDirty();
    }

    public void addWeight(int patId, int dweight){
        patternWeights[patId] += dweight;
        //System.out.println("patId: "+patId+", patternWeights[patId]: "+patternWeights[patId]);
        weightsChanged = true;
    }

    public void removeWeight(int patId, int dweight){
        patternWeights[patId] -= dweight;
        //System.out.println("patId: "+patId+", patternWeights[patId]: "+patternWeights[patId]);
        weightsChanged = true;
    }

    public void setWeight(int patId, int weight){
        patternWeights[patId] = weight;
        weightsChanged = true;
    }

    @Override
    public void store(){
        weightsChanged = false;
        System.arraycopy(patternWeights,0,storedPatternWeights,0,patternWeights.length);
        System.arraycopy(patternLogLikelihoods,0,storedPatternLogLikelihoods,0,patternLogLikelihoods.length);
        super.store();
    }

    @Override
    public void restore(){
        int[] temp1 = patternWeights;
        patternWeights = storedPatternWeights;
        storedPatternWeights = temp1;

        double[] temp2 = patternLogLikelihoods;
        patternLogLikelihoods = storedPatternLogLikelihoods;
        storedPatternLogLikelihoods = temp2;
        
        weightsChanged = false;

        super.restore();
    }

    public int[] getPatternWeights(){
        return patternWeights;
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
        System.out.println("modelID: "+((QuietRealParameter)(((QuietSiteModel)m_siteModel)).getRateParameter()).getIDNumber());
    }

}
