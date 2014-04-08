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
    public void initAndValidate() throws Exception {
    	// sanity check: alignment should have same #taxa as tree
    	if (m_data.get().getNrTaxa() != m_tree.get().getLeafNodeCount()) {
    		throw new Exception("The number of nodes in the tree does not match the number of sequences");
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

        int nodeCount = m_tree.get().getNodeCount();
        if (!(m_pSiteModel.get() instanceof SiteModel.Base)) {
        	throw new Exception ("siteModel input should be of type SiteModel.Base");
        }
        m_siteModel = (SiteModel.Base) m_pSiteModel.get();
        m_siteModel.setDataType(m_data.get().getDataType());
        m_substitutionModel = m_siteModel.m_pSubstModel.get();

        if (m_pBranchRateModel.get() != null) {
        	m_branchRateModel = m_pBranchRateModel.get();
        } else {
            m_branchRateModel = new StrictClockModel();
        }
    	m_branchLengths = new double[nodeCount];
    	m_StoredBranchLengths = new double[nodeCount];

        int nStateCount = m_data.get().getMaxStateCount();
        int nPatterns = m_data.get().getPatternCount();
        if (nStateCount == 4) {
            m_likelihoodCore = new BeerLikelihoodCore4();
        } else {
            m_likelihoodCore = new BeerLikelihoodCore(nStateCount);
        }
        //System.err.println("TreeLikelihood uses " + m_likelihoodCore.getClass().getName());

        m_fProportionInvariant = m_siteModel.getProportianInvariant();
        m_siteModel.setPropInvariantIsCategory(false);
        if (m_fProportionInvariant > 0) {
        	calcConstantPatternIndices(nPatterns, nStateCount);
        }

        initCore();

        m_fPatternLogLikelihoods = new double[nPatterns];
        storedPatternLogLikelihoods = new double[nPatterns];
        m_fRootPartials = new double[nPatterns * nStateCount];
        m_nMatrixSize = (nStateCount +1)* (nStateCount+1);
        m_fProbabilities = new double[(nStateCount +1)* (nStateCount+1)];
        Arrays.fill(m_fProbabilities, 1.0);

        if (m_data.get() instanceof AscertainedAlignment) {
            m_bAscertainedSitePatterns = true;
        }

        if(patternWeights.length != m_data.get().getPatternCount()){
            throw new RuntimeException("Pattern counts differ!");
        }
    }

    @Override
    public double calculateLogP() throws Exception {
    	/*if (m_beagle != null) {
    		logP =  m_beagle.calculateLogP();
    		return logP;
    	}*/
        Tree tree = m_tree.get();
        //System.err.println("m_nHasDirt: "+m_nHasDirt);
        if(m_nHasDirt > -1){
       	    traverse(tree.getRoot());
        }

        calcLogP();

        m_nScale++;
        if (logP > 0 || (m_likelihoodCore.getUseScaling() && m_nScale > X)) {
            //System.err.println("Switch off scaling");
            m_likelihoodCore.setUseScaling(1.0);
            m_likelihoodCore.unstore();
            m_nHasDirt = Tree.IS_FILTHY;
            X *= 2;
           	traverse(tree.getRoot());
            calcLogP();
            return logP;
        } else if (logP == Double.NEGATIVE_INFINITY && m_fScale < 10) { // && !m_likelihoodCore.getUseScaling()) {
        	m_nScale = 0;
        	m_fScale *= 1.01;
            System.err.println("Turning on scaling to prevent numeric instability " + m_fScale);
            m_likelihoodCore.setUseScaling(m_fScale);
            m_likelihoodCore.unstore();
            m_nHasDirt = Tree.IS_FILTHY;
           	traverse(tree.getRoot());
            calcLogP();
            //System.err.println(getID()+": "+logP);
            return logP;
        }

        //System.err.println(getID()+": "+logP);
        return logP;
    }

    protected void calcLogP() throws Exception {
        logP = 0.0;
        if (m_bAscertainedSitePatterns) {
            double ascertainmentCorrection = ((AscertainedAlignment)m_data.get()).getAscertainmentCorrection(m_fPatternLogLikelihoods);
            for (int i = 0; i < m_data.get().getPatternCount(); i++) {
            	logP += (m_fPatternLogLikelihoods[i] - ascertainmentCorrection) * patternWeights[i];
            }
        } else {

	        for (int i = 0; i < m_data.get().getPatternCount(); i++) {
	            logP += m_fPatternLogLikelihoods[i] * patternWeights[i];
	        }
        }
    }


    public double getPatternLogLikelihood(int iPat){
        return m_fPatternLogLikelihoods[iPat];
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
        m_nHasDirt = Tree.IS_CLEAN;

        if(weightsChanged){
            m_nHasDirt = -1;
            //System.err.println("flag2");
            return true;            
        }

        if (m_branchRateModel != null && m_branchRateModel.isDirtyCalculation()) {
            m_nHasDirt = Tree.IS_FILTHY;
            //System.err.println("flag3");
            return true;
        }

        if (m_data.get().isDirtyCalculation()) {
            m_nHasDirt = Tree.IS_FILTHY;
            //System.err.println("flag4");
            return true;
        }

        if (m_siteModel.isDirtyCalculation()) {
            m_nHasDirt = Tree.IS_DIRTY;
            //System.err.println("flag5");
            return true;
        }

        return m_tree.get().somethingIsDirty();
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
        System.arraycopy(m_fPatternLogLikelihoods,0,storedPatternLogLikelihoods,0,m_fPatternLogLikelihoods.length);
        super.store();
    }

    @Override
    public void restore(){
        int[] temp1 = patternWeights;
        patternWeights = storedPatternWeights;
        storedPatternWeights = temp1;

        double[] temp2 = m_fPatternLogLikelihoods;
        m_fPatternLogLikelihoods = storedPatternLogLikelihoods;
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
        System.out.println("modelID: "+((SwitchingNtdBMA)m_substitutionModel).getIDNumber());
        System.out.println("modelID: "+((QuietRealParameter)(((QuietSiteModel)m_siteModel)).getRateParameter()).getIDNumber());
    }

}
