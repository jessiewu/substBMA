package beast.evolution.operators;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.core.StateNode;
import beast.core.parameter.RealParameter;
import beast.util.Randomizer;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Chieh-Hsi Wu
 */
@Description("Same as an updown operator but can scale parameters in exponential space.")
public class ExtendedUpDownOperator extends Operator {
    public Input<Double> m_scaleFactor = new Input<Double>("scaleFactor",
			"magnitude factor used for scaling", Input.Validate.REQUIRED);
	public Input<List<StateNode>> m_up = new Input<List<StateNode>>("up",
			"zero or more items to scale upwards", new ArrayList<StateNode>());
	public Input<List<StateNode>> m_down = new Input<List<StateNode>>("down",
			"zero or more items to scale downwards", new ArrayList<StateNode>());
    public Input<Boolean> m_bOptimise = new Input<Boolean>("optimise","flag to indicate that the scale factor is automatically changed in order to acheive a good acceptance rate (default true)", true);

    public Input<String> isExpInput = new Input<String>("isExp","Whether to scale the exponential of the parameter value", Input.Validate.REQUIRED);
	double m_fScaleFactor;
    boolean[] isExpUp;
    boolean[] isExpDown;

	@Override
	public void initAndValidate() throws Exception {
		m_fScaleFactor = m_scaleFactor.get();
		// sanity checks
		if (m_up.get().size() + m_down.get().size() == 0) {
			throw new Exception("At least one up or down item must be specified");
		}
		if (m_up.get().size() == 0 || m_down.get().size() == 0) {
			System.err.println("WARNING: no " + (m_up.get().size() == 0 ? "up" : "down") + " item specified in UpDownOperator");
		}
        String isExpStr = isExpInput.get();
        String[] isExpStrElt = isExpStr.split("\\s+");
        isExpUp = new boolean[m_up.get().size()];
        isExpDown = new boolean[m_down.get().size()];
        if(isExpUp.length + isExpDown.length == isExpStrElt.length){

            int k = 0;

            for(int i = 0; i < isExpUp.length; i++){
                isExpUp[i] = Boolean.parseBoolean(isExpStrElt[k++]);
            }

            for(int i = 0; i < isExpDown.length; i++){
                isExpDown[i] = Boolean.parseBoolean(isExpStrElt[k++]);
            }


        }else if(isExpStrElt.length == 1){
            boolean isExpTemp = Boolean.parseBoolean(isExpStrElt[0]);
            for(int i = 0; i < isExpUp.length; i++){
                isExpUp[i] = isExpTemp;
            }

            for(int i = 0; i < isExpDown.length; i++){
                isExpDown[i] = isExpTemp;
            }

        }else{
            throw new RuntimeException(
                    "The dimension of the isExp specification ("+isExpStrElt.length+
                    ") should equal to either 1 or the dimension of the number of up ("+isExpUp.length+") + down("+isExpDown.length+")."
            );

        }


	}

	/**
	 * override this for proposals,
	 *
	 * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal
	 *         should not be accepted
	 **/
	@Override
	public final double proposal() {
        double r = Randomizer.nextDouble();
		final double scale = (m_fScaleFactor + (r * ((1.0 / m_fScaleFactor) - m_fScaleFactor)));
        final double scaleLog = (m_fScaleFactor + r*-2*m_fScaleFactor);
		int goingUp = 0, goingDown = 0;

		try {
            List<StateNode> upList = m_up.get();
			for (int i = 0; i < upList.size(); i++) {
				StateNode up = upList.get(i).getCurrentEditable(this);
                if(isExpUp[i]){
                    if(up instanceof RealParameter){
                        RealParameter upParameter = (RealParameter)up;
                        int upParameterDim = upParameter.getDimension();
                        for(int j = 0; j < upParameterDim; j++){
                            double upVal = upParameter.getValue(j);

                            upParameter.setValue(Math.log(Math.exp(upVal)*scale));

                        }

                    }else{
                        throw new RuntimeException("If isExp is true, then the state node must be a real parameter.");
                    }

                }else{


				    goingUp += up.scale(scale);
                }
			}

            List<StateNode> downList = m_down.get();

			for (int i = 0; i < downList.size(); i++) {
				StateNode down = downList.get(i).getCurrentEditable(this);
                if(isExpDown[i]){
                    if(down instanceof RealParameter){
                        RealParameter downParameter = (RealParameter)down;
                        int downParameterDim = downParameter.getDimension();
                        for(int j = 0; j < downParameterDim; j++){
                            double downVal = downParameter.getValue(j);

                            downParameter.setValue(Math.log(Math.exp(downVal)*scale));

                        }

                    }else{
                        throw new RuntimeException("If isExp is true, then the state node must be a real parameter.");
                    }
                }else{
				    goingDown += down.scale(1.0 / scale);
                }
			}
		} catch (Exception e) {
			// scale resulted in invalid StateNode, abort proposal
			return Double.NEGATIVE_INFINITY;
		}
        //System.out.println(goingUp - goingDown );
		return (goingUp - goingDown - 2) * Math.log(scale);
        //return -2 * Math.log(scale);
	}

	/**
	 * automatic parameter tuning *
	 */
	@Override
	public void optimize(double logAlpha) {
		if (m_bOptimise.get()) {
	        double fDelta = calcDelta(logAlpha);
	        fDelta += Math.log(1.0 / m_fScaleFactor - 1.0);
	        m_fScaleFactor = 1.0 / (Math.exp(fDelta) + 1.0);
		}
	}

    @Override
    public double getCoercableParameterValue() {
        return m_fScaleFactor;
    }

    @Override
    public String getPerformanceSuggestion() {
        double prob = m_nNrAccepted/(m_nNrAccepted+m_nNrRejected+0.0);
        double targetProb = getTargetAcceptanceProbability();

        double ratio = prob / targetProb;
        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        // new scale factor
        double sf = Math.pow(m_fScaleFactor, ratio);

        DecimalFormat formatter = new DecimalFormat("#.###");
        if (prob < 0.10) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else if (prob > 0.40) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else return "";
    }

}
