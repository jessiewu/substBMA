package beast.evolution.likelihood;

import beast.core.Description;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.AscertainedAlignment;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.branchratemodel.StrictClockModel;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * @author Chieh-Hsi Wu
 */
@Description("The tree likelihood used and created on the fly during the MCMC when using the DPM model.")
public class QuietTreeLikelihood extends TreeLikelihood{
    protected Alignment data;
    protected Tree tree;
    protected boolean useAmbiguities;
    private Scaling scale;



    public QuietTreeLikelihood(){}

    public QuietTreeLikelihood(Alignment data,
                               Tree tree,
                               boolean useAmbiguities,
                               SiteModel siteModel,
                               BranchRateModel.Base branchRateModel){
        this.data = data;
        this.tree = tree;
        this.useAmbiguities = useAmbiguities;
        m_siteModel = siteModel;
        m_branchRateModel = branchRateModel;
        m_substitutionModel = (SubstitutionModel.Base)m_siteModel.getSubstitutionModel();
        setup();

    }

    @Override
    public void initAndValidate() throws Exception {

        data = m_data.get();
        tree = m_tree.get();
        if (!(m_pSiteModel.get() instanceof SiteModel.Base)) {
        	throw new Exception ("siteModel input should be of type SiteModel.Base");
        }
        m_siteModel = (SiteModel.Base) m_pSiteModel.get();
        if (m_pBranchRateModel.get() != null) {
            m_branchRateModel = m_pBranchRateModel.get();
        } else {
            m_branchRateModel = new StrictClockModel();
        }
        useAmbiguities = m_useAmbiguities.get();
        scale = scaling.get();
        m_substitutionModel = (SubstitutionModel.Base) m_siteModel.getSubstitutionModel();
        System.out.println("ID: "+m_siteModel.getID());
        // sanity check: alignment should have same #taxa as tree
        if (data.getNrTaxa() != tree.getLeafNodeCount()) {
            throw new Exception("The number of nodes in the tree does not match the number of sequences");
        }

        /*m_beagle = null;
        m_beagle = new BeagleTreeLikelihood();
        m_beagle.initByName("data", data, "tree", tree, "siteModel", m_siteModel, "branchRateModel", m_pBranchRateModel.get(), "useAmbiguities", m_useAmbiguities.get(), "scaling", scale.toString());
        if (m_beagle.beagle != null) {
            //a Beagle instance was found, so we use it
            return;
        }
        // No Beagle instance was found, so we use the good old java likelihood core
        m_beagle = null;*/
        setup();


    }

    protected void setup(){



        int nodeCount = tree.getNodeCount();

        m_siteModel.setDataType(data.getDataType());



        m_branchLengths = new double[nodeCount];
        m_StoredBranchLengths = new double[nodeCount];

        int nStateCount = data.getMaxStateCount();
        int nPatterns = data.getPatternCount();
        if (nStateCount == 4) {
            m_likelihoodCore = new CHWLikelihoodCore4();
        } else {
            m_likelihoodCore = new CHWLikelihoodCore(nStateCount);
        }
        System.out.println("TreeLikelihood uses " + m_likelihoodCore.getClass().getName());

        m_fProportionInvariant = m_siteModel.getProportianInvariant();
        m_siteModel.setPropInvariantIsCategory(false);
        if (m_fProportionInvariant > 0) {
            calcConstantPatternIndices(nPatterns, nStateCount);
        }

        initCore();

        m_fPatternLogLikelihoods = new double[nPatterns];
        m_fRootPartials = new double[nPatterns * nStateCount];
        m_nMatrixSize = (nStateCount + 1) * (nStateCount + 1);
        m_fProbabilities = new double[(nStateCount + 1) * (nStateCount + 1)];
        Arrays.fill(m_fProbabilities, 1.0);

        if (data instanceof AscertainedAlignment) {
            m_bAscertainedSitePatterns = true;
        }

    }

    /**
     * Determine indices of m_fRootProbabilities that need to be updates
     * // due to sites being invariant. If none of the sites are invariant,
     * // the 'site invariant' category does not contribute anything to the
     * // root probability. If the site IS invariant for a certain character,
     * // taking ambiguities in account, there is a contribution of 1 from
     * // the 'site invariant' category.
     */
    void calcConstantPatternIndices(final int nPatterns, final int nStateCount) {
        m_iConstantPattern = new ArrayList<Integer>();
        for (int i = 0; i < nPatterns; i++) {
            final int[] pattern = data.getPattern(i);
            final boolean[] bIsInvariant = new boolean[nStateCount];
            Arrays.fill(bIsInvariant, true);
            for (final int state : pattern) {
                final boolean[] bStateSet = data.getStateSet(state);
                if (useAmbiguities || !data.getDataType().isAmbiguousState(state)) {
                    for (int k = 0; k < nStateCount; k++) {
                        bIsInvariant[k] &= bStateSet[k];
                    }
                }
            }
            for (int k = 0; k < nStateCount; k++) {
                if (bIsInvariant[k]) {
                    m_iConstantPattern.add(i * nStateCount + k);
                }
            }
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

        if (useAmbiguities) {
            setPartials(tree.getRoot(), data.getPatternCount());
        } else {
            setStates(tree.getRoot(), data.getPatternCount());
        }
        m_nHasDirt = Tree.IS_FILTHY;
        for (int i = 0; i < intNodeCount; i++) {
            m_likelihoodCore.createNodePartials(extNodeCount + i);
        }
    }

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

    @Override
    public double calculateLogP() throws Exception {
        if (m_beagle != null) {
            logP = m_beagle.calculateLogP();
            return logP;
        }

        traverse(tree.getRoot());
        calcLogP();

        m_nScale++;
        if (logP > 0 || (m_likelihoodCore.getUseScaling() && m_nScale > X)) {
//            System.err.println("Switch off scaling");
//            m_likelihoodCore.setUseScaling(1.0);
//            m_likelihoodCore.unstore();
//            m_nHasDirt = Tree.IS_FILTHY;
//            X *= 2;
//            traverse(tree.getRoot());
//            calcLogP();
//            return logP;
        } else if (logP == Double.NEGATIVE_INFINITY && m_fScale < 10) { // && !m_likelihoodCore.getUseScaling()) {
            m_nScale = 0;
            m_fScale *= 1.01;
            System.err.println("Turning on scaling to prevent numeric instability " + m_fScale);
            m_likelihoodCore.setUseScaling(m_fScale);
            m_likelihoodCore.unstore();
            m_nHasDirt = Tree.IS_FILTHY;
            traverse(tree.getRoot());
            calcLogP();
            return logP;
        }
        return logP;
    }

    void calcLogP() throws Exception {
        logP = 0.0;
        if (m_bAscertainedSitePatterns) {
            final double ascertainmentCorrection = ((AscertainedAlignment) data).getAscertainmentCorrection(m_fPatternLogLikelihoods);
            for (int i = 0; i < data.getPatternCount(); i++) {
                logP += (m_fPatternLogLikelihoods[i] - ascertainmentCorrection) * data.getPatternWeight(i);
            }
        } else {
            for (int i = 0; i < data.getPatternCount(); i++) {
                logP += m_fPatternLogLikelihoods[i] * data.getPatternWeight(i);
            }
        }
    }


    /** CalculationNode methods **/

    /**
     * check state for changed variables and update temp results if necessary *
     */
    @Override
    protected boolean requiresRecalculation() {
        if (m_beagle != null) {
            return m_beagle.requiresRecalculation();
        }
        m_nHasDirt = Tree.IS_CLEAN;

        if (data.isDirtyCalculation()) {
            m_nHasDirt = Tree.IS_FILTHY;
            return true;
        }
        if (m_siteModel.isDirtyCalculation()) {
            m_nHasDirt = Tree.IS_DIRTY;
            return true;
        }
        if (m_branchRateModel != null && m_branchRateModel.isDirtyCalculation()) {
            m_nHasDirt = Tree.IS_DIRTY;
            return true;
        }
        return tree.somethingIsDirty();
    }

    /**
     * @return a list of unique ids for the state nodes that form the argument
     */
    public List<String> getArguments() {
        return Collections.singletonList(data.getID());
    }

    public Tree getTree(){
        return tree;
    }
}
