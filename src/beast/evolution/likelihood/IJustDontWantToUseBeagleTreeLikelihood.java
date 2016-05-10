package beast.evolution.likelihood;

import beast.evolution.alignment.AscertainedAlignment;
import beast.evolution.branchratemodel.StrictClockModel;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.SubstitutionModel;

import java.util.Arrays;

/**
 * Created by IntelliJ IDEA.
 * User: Jessie Wu
 * Date: 25/01/13
 * Time: 3:56 PM
 * To change this template use File | Settings | File Templates.
 */
public class IJustDontWantToUseBeagleTreeLikelihood extends TreeLikelihood {
    @Override
    public void initAndValidate() {
        // sanity check: alignment should have same #taxa as tree
        if (dataInput.get().getNrTaxa() != treeInput.get().getLeafNodeCount()) {
            throw new RuntimeException("The number of nodes in the tree does not match the number of sequences");
        }

        // No Beagle instance was found, so we use the good old java likelihood core
        beagle = null;

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

        if (dataInput.get() instanceof AscertainedAlignment) {
            useAscertainedSitePatterns = true;
        }
    }
}
