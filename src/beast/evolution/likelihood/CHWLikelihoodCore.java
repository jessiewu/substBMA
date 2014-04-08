package beast.evolution.likelihood;

/**
 * Created by IntelliJ IDEA.
 * User: cwu080
 * Date: 17/09/12
 * Time: 3:33 PM
 * To change this template use File | Settings | File Templates.
 */
public class CHWLikelihoodCore extends BeerLikelihoodCore {
    private double m_fScalingThreshold = 1.0E-100;
    public CHWLikelihoodCore(int nStateCount) {
        super(nStateCount);
    }

    protected void scalePartials(int iNodeIndex) {
//        int v = 0;
//    	double [] fPartials = m_fPartials[m_iCurrentPartialsIndices[iNodeIndex]][iNodeIndex];
//        for (int i = 0; i < m_nPatternCount; i++) {
//            for (int k = 0; k < m_nMatrixCount; k++) {
//                for (int j = 0; j < m_nStateCount; j++) {
//                	fPartials[v] *= SCALE;
//                	v++;
//                }
//            }
//        }
        int u = 0;

        for (int i = 0; i < m_nPatterns; i++) {

            double scaleFactor = 0.0;
            int v = u;
            for (int k = 0; k < m_nMatrices; k++) {
                for (int j = 0; j < m_nStates; j++) {
                    if (m_fPartials[m_iCurrentPartials[iNodeIndex]][iNodeIndex][v] > scaleFactor) {
                        scaleFactor = m_fPartials[m_iCurrentPartials[iNodeIndex]][iNodeIndex][v];
                    }
                    v++;
                }
                v += (m_nPatterns - 1) * m_nStates;
            }

            if (scaleFactor < m_fScalingThreshold && scaleFactor > 0.0) {

                v = u;
                for (int k = 0; k < m_nMatrices; k++) {
                    for (int j = 0; j < m_nStates; j++) {
                        m_fPartials[m_iCurrentPartials[iNodeIndex]][iNodeIndex][v] /= scaleFactor;
                        v++;
                    }
                    v += (m_nPatterns - 1) * m_nStates;
                }
                m_fScalingFactors[m_iCurrentPartials[iNodeIndex]][iNodeIndex][i] = Math.log(scaleFactor);

            } else {
                m_fScalingFactors[m_iCurrentPartials[iNodeIndex]][iNodeIndex][i] = 0.0;
            }
            u += m_nStates;


        }
    }
}
