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

        for (int i = 0; i < nrOfPatterns; i++) {

            double scaleFactor = 0.0;
            int v = u;
            for (int k = 0; k < nrOfMatrices; k++) {
                for (int j = 0; j < nrOfStates; j++) {
                    if (partials[currentPartialsIndex[iNodeIndex]][iNodeIndex][v] > scaleFactor) {
                        scaleFactor = partials[currentPartialsIndex[iNodeIndex]][iNodeIndex][v];
                    }
                    v++;
                }
                v += (nrOfPatterns - 1) * nrOfStates;
            }

            if (scaleFactor < m_fScalingThreshold && scaleFactor > 0.0) {

                v = u;
                for (int k = 0; k < nrOfMatrices; k++) {
                    for (int j = 0; j < nrOfStates; j++) {
                        partials[currentPartialsIndex[iNodeIndex]][iNodeIndex][v] /= scaleFactor;
                        v++;
                    }
                    v += (nrOfPatterns - 1) * nrOfStates;
                }
                scalingFactors[currentPartialsIndex[iNodeIndex]][iNodeIndex][i] = Math.log(scaleFactor);

            } else {
                scalingFactors[currentPartialsIndex[iNodeIndex]][iNodeIndex][i] = 0.0;
            }
            u += nrOfStates;


        }
    }
}
