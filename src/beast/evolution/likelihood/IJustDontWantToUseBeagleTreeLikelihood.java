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
    public void initAndValidate() throws Exception {
        // sanity check: alignment should have same #taxa as tree
        if (m_data.get().getNrTaxa() != m_tree.get().getLeafNodeCount()) {
            throw new Exception("The number of nodes in the tree does not match the number of sequences");
        }

        // No Beagle instance was found, so we use the good old java likelihood core
        m_beagle = null;

        int nodeCount = m_tree.get().getNodeCount();
        if (!(m_pSiteModel.get() instanceof SiteModel.Base)) {
        	throw new Exception ("siteModel input should be of type SiteModel.Base");
        }
        m_siteModel = (SiteModel.Base) m_pSiteModel.get();
        m_siteModel.setDataType(m_data.get().getDataType());
        m_substitutionModel = (SubstitutionModel.Base) m_siteModel.m_pSubstModel.get();

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

        if (m_data.get() instanceof AscertainedAlignment) {
            m_bAscertainedSitePatterns = true;
        }
    }
}
