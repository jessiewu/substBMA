package beast.evolution.operators;

import beast.core.parameter.BooleanParameter;
import beast.core.parameter.ParameterList;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;

/**
 * Created by IntelliJ IDEA.
 * User: Jessie Wu
 * Date: 25/01/14
 * Time: 3:41 PM
 * To change this template use File | Settings | File Templates.
 */
public class PLScaleOperator2 extends PLScaleOperator{
    @Override
    public double proposal() {
    	try {

            double hastingsRatio;
            final double scale = getScaler();

            if (m_bIsTreeScaler) {
                Tree tree = treeInput.get(this);
                // scale the beast.tree
                final int nInternalNodes = tree.scale(scale);
                return Math.log(scale) * (nInternalNodes - 2);
            }

            // not a tree scaler, so scale a parameter
            final boolean bScaleAll = scaleAllInput.get();
            final int nDegreesOfFreedom = degreesOfFreedomInput.get();
            final boolean bScaleAllIndependently = scaleAllIndependentlyInput.get();

            final ParameterList param = parameterListInput.get(this);

            int pIndex = 0;
            if(param.getDimension()>1)
                pIndex = Randomizer.nextInt(param.getDimension());

            //assert param.getLower() != null  && param.getUpper() != null;

            final int dim = param.getDimension();

            if (bScaleAllIndependently) {
                // update all dimensions independently.
                hastingsRatio = 0;
                for (int i = 0; i < dim; i++) {

                    final double scaleOne = getScaler();
                    final double newValue = scaleOne * param.getValue(pIndex,i);

                    hastingsRatio -= Math.log(scaleOne);

                    if ( outsideBounds(newValue, param) ) {
                        return Double.NEGATIVE_INFINITY;
                    }

                    param.setValue(pIndex,i, newValue);
                }
            } else if (bScaleAll) {
                // update all dimensions
                // hasting ratio is dim-2 times of 1dim case. would be nice to have a reference here
                // for the proof. It is supposed to be somewhere in an Alexei/Nicholes article.
                final int df = (nDegreesOfFreedom > 0) ? -nDegreesOfFreedom  : dim - 2;
                hastingsRatio = df * Math.log(scale);

                // all Values assumed independent!
                for (int i = 0; i < dim; i++) {
                    final double newValue = param.getValue(pIndex,i) * scale;

                    if ( outsideBounds(newValue, param) ) {
                        return Double.NEGATIVE_INFINITY;
                    }
                    param.setValue(pIndex,i, newValue);
                }
            } else {
                hastingsRatio = -0.5*Math.log(scale);

                // which position to scale
                int index = -1;
                final BooleanParameter indicators = indicatorInput.get();
                if ( indicators != null ) {
                    final int nDim = indicators.getDimension();
                    Boolean [] indicator = indicators.getValues();
                    final boolean impliedOne = nDim == (dim - 1);

                    // available bit locations. there can be hundreds of them. scan list only once.
                    int[] loc = new int[nDim + 1];
                    int nLoc = 0;

                    if (impliedOne) {
                        loc[nLoc] = 0;
                        ++nLoc;
                    }
                    for (int i = 0; i < nDim; i++) {
                        if ( indicator[i] ) {
                            loc[nLoc] = i + (impliedOne ? 1 : 0);
                            ++nLoc;
                        }
                    }

                    if (nLoc > 0) {
                        final int rand = Randomizer.nextInt(nLoc);
                        index = loc[rand];
                    } else {
                        return Double.NEGATIVE_INFINITY; // no active indicators
                    }

                } else {
                    // any is good
                    index = Randomizer.nextInt(dim);
                }

                final double oldValue = param.getValue(pIndex,index);

                if (oldValue == 0) {
                    // Error: parameter has value 0 and cannot be scaled
                    return Double.NEGATIVE_INFINITY;
                }

                final double newValue = scale * oldValue;

                if ( outsideBounds(newValue, param) ) {
                    // reject out of bounds scales
                    return Double.NEGATIVE_INFINITY;
                }

                param.setValue(pIndex,index, newValue);
                // provides a hook for subclasses
                //cleanupOperation(newValue, oldValue);
            }

            return hastingsRatio;

        } catch (Exception e) {
            // whatever went wrong, we want to abort this operation...
            return Double.NEGATIVE_INFINITY;
        }
    }
}
