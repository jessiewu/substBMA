package beast.evolution.operators;

import beast.core.Description;
import beast.core.parameter.ParameterList;
import beast.core.parameter.RealParameter;
import beast.core.parameter.BooleanParameter;
import beast.core.Input;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;

/**
 * @author Chieh-Hsi Wu
 */
@Description("ScaleOperator performed on a list of parameters.")
public class PLScaleOperator extends ScaleOperator{
    public Input<ParameterList> parameterListInput =
            new Input<ParameterList>(
                    "parameters",
                    "the parameter to operate a random walk on.",
                    Input.Validate.REQUIRED
            );
    public PLScaleOperator(){
        m_pParameter.setRule(Input.Validate.OPTIONAL);
        m_pTree.setRule(Input.Validate.OPTIONAL);
    }

    public void initAndValidate() throws Exception {
        super.initAndValidate();
    }



    protected boolean outsideBounds(double value, ParameterList param) {
        final Double l = param.getLower();
        final Double h = param.getUpper();

        return ( value < l || value > h );
        //return (l != null && value < l || h != null && value > h);
    }
/** override this for proposals,
	 * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should not be accepted **/
    @Override
    public double proposal() {
    	try {

            double hastingsRatio;
            final double scale = getScaler();

            if (m_bIsTreeScaler) {
                Tree tree = m_pTree.get(this);
                // scale the beast.tree
                final int nInternalNodes = tree.scale(scale);
                return Math.log(scale) * (nInternalNodes - 2);
            }

            // not a tree scaler, so scale a parameter
            final boolean bScaleAll = m_pScaleAll.get();
            final int nDegreesOfFreedom = m_pDegreesOfFreedom.get();
            final boolean bScaleAllIndependently = m_pScaleAllIndependently.get();

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
                hastingsRatio = -Math.log(scale);

                // which position to scale
                int index = -1;
                final BooleanParameter indicators = m_indicator.get();
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
