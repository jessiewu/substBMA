package beast.evolution.likelihood;

import beast.core.Description;
import beast.core.Input;
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

    public Input<GenericTreeLikelihood> trueLikelihoodInput = new Input<>("trueLikelihood",
            "If present, use this as source of branch rate model.");


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
        this.branchRateModel = branchRateModel;
        substitutionModel = m_siteModel.getSubstitutionModel();
        setup();

    }

    @Override
    public void initAndValidate() {

        data = dataInput.get();
        tree = (Tree) treeInput.get();
        if (!(siteModelInput.get() instanceof SiteModel.Base)) {
        	throw new RuntimeException ("siteModel input should be of type SiteModel.Base");
        }
        m_siteModel = (SiteModel.Base) siteModelInput.get();
        if (branchRateModelInput.get() != null) {
            branchRateModel = branchRateModelInput.get();
        } else {
            if (trueLikelihoodInput.get() != null) {
                branchRateModel = trueLikelihoodInput.get().branchRateModelInput.get();
            } else {
                branchRateModel = new StrictClockModel();
            }
        }
        useAmbiguities = m_useAmbiguities.get();
        scale = scaling.get();
        substitutionModel = m_siteModel.getSubstitutionModel();
        System.out.println("ID: "+m_siteModel.getID());
        // sanity check: alignment should have same #taxa as tree
        if (data.getNrTaxa() != tree.getLeafNodeCount()) {
            throw new RuntimeException("The number of nodes in the tree does not match the number of sequences");
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
        storedBranchLengths = new double[nodeCount];

        int nStateCount = data.getMaxStateCount();
        int nPatterns = data.getPatternCount();
        if (nStateCount == 4) {
            likelihoodCore = new CHWLikelihoodCore4();
        } else {
            likelihoodCore = new CHWLikelihoodCore(nStateCount);
        }
        System.out.println("TreeLikelihood uses " + likelihoodCore.getClass().getName());

        proportionInvariant = m_siteModel.getProportionInvariant();
        m_siteModel.setPropInvariantIsCategory(false);
        if (proportionInvariant > 0) {
            calcConstantPatternIndices(nPatterns, nStateCount);
        }

        initCore();

        patternLogLikelihoods = new double[nPatterns];
        m_fRootPartials = new double[nPatterns * nStateCount];
        matrixSize = (nStateCount + 1) * (nStateCount + 1);
        probabilities = new double[(nStateCount + 1) * (nStateCount + 1)];
        Arrays.fill(probabilities, 1.0);

        if (data instanceof AscertainedAlignment) {
            useAscertainedSitePatterns = true;
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
        constantPattern = new ArrayList<>();
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
                    constantPattern.add(i * nStateCount + k);
                }
            }
        }
    }

    protected void initCore() {
        final int nodeCount = tree.getNodeCount();
        likelihoodCore.initialize(
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
        hasDirt = Tree.IS_FILTHY;
        for (int i = 0; i < intNodeCount; i++) {
            likelihoodCore.createNodePartials(extNodeCount + i);
        }
    }

    /**
     * set leaf states in likelihood core *
     */
    protected void setStates(Node node, int patternCount) {
        if (node.isLeaf()) {
            int i;
            int[] states = new int[patternCount];
            int iTaxon = data.getTaxonIndex(node.getID());
            for (i = 0; i < patternCount; i++) {
                states[i] = data.getPattern(iTaxon, i);
            }
            likelihoodCore.setNodeStates(node.getNr(), states);

        } else {
            setStates(node.getLeft(), patternCount);
            setStates(node.getRight(), patternCount);
        }
    }


    /**
     * set leaf partials in likelihood core *
     */
    protected void setPartials(Node node, int patternCount) {
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
            likelihoodCore.setNodePartials(node.getNr(), partials);

        } else {
            //System.out.println(m_likelihoodCore == null);
            setPartials(node.getLeft(), patternCount);
            setPartials(node.getRight(), patternCount);
        }
    }

    @Override
    public double calculateLogP() {
        if (beagle != null) {
            logP = beagle.calculateLogP();
            return logP;
        }

        traverse(tree.getRoot());
        calcLogP();

        m_nScale++;
        if (logP > 0 || (likelihoodCore.getUseScaling() && m_nScale > X)) {
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
            likelihoodCore.setUseScaling(m_fScale);
            likelihoodCore.unstore();
            hasDirt = Tree.IS_FILTHY;
            traverse(tree.getRoot());
            calcLogP();
            return logP;
        }
        return logP;
    }

    void calcLogP() {
        logP = 0.0;
        if (useAscertainedSitePatterns) {
            final double ascertainmentCorrection = ((AscertainedAlignment) data).getAscertainmentCorrection(patternLogLikelihoods);
            for (int i = 0; i < data.getPatternCount(); i++) {
                logP += (patternLogLikelihoods[i] - ascertainmentCorrection) * data.getPatternWeight(i);
            }
        } else {
            for (int i = 0; i < data.getPatternCount(); i++) {
                logP += patternLogLikelihoods[i] * data.getPatternWeight(i);
            }
        }
    }


    /** CalculationNode methods **/

    /**
     * check state for changed variables and update temp results if necessary *
     */
    @Override
    protected boolean requiresRecalculation() {
        if (beagle != null) {
            return beagle.requiresRecalculation();
        }
        hasDirt = Tree.IS_CLEAN;

        if (data.isDirtyCalculation()) {
            hasDirt = Tree.IS_FILTHY;
            return true;
        }
        if (m_siteModel.isDirtyCalculation()) {
            hasDirt = Tree.IS_DIRTY;
            return true;
        }
        if (branchRateModel != null && branchRateModel.isDirtyCalculation()) {
            hasDirt = Tree.IS_DIRTY;
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
