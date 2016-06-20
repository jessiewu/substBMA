package beast.math.distributions;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.parameter.ParameterList;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * @author Chieh-Hsi Wu
 */
@Description("This class is a wrapper that provides prior to parameter list.")
public class ParameterListsPrior extends Distribution {

    final public Input<ParametricDistribution> distInput = new Input<>(
            "distr",
            "distribution used to calculate prior, e.g. normal, beta, gamma.",
            Input.Validate.REQUIRED);

    public Input<List<ParameterList>> xListsInput = new Input<>(
            "xList",
            "points at which the density is calculated",
            new ArrayList<ParameterList>(),
            Input.Validate.REQUIRED
    );

    public Input<Boolean> applyToListInput = new Input<>(
            "applyToList",
            "Whether the prior is applied to the entire list",
            Input.Validate.REQUIRED
    );

    protected ParametricDistribution dist;

    public ParameterListsPrior() { }

    boolean applyToList;
    public void initAndValidate(){

        applyToList = applyToListInput.get();

        dist = distInput.get();
        calculateLogP();
    }

    @Override
	public double calculateLogP() {
        List<ParameterList> parameterLists = xListsInput.get();
        logP = ((CompoundDirichletProcess)dist).calcLogP(parameterLists);

		return logP;
	}

    @Override
    public List<String> getArguments() { return null; }

    @Override
    public List<String> getConditions() { return null; }

    @Override
    public void sample(State state, Random random) {

    }

    public String getParameterName() {
    	String name = "";
    	for (ParameterList list : xListsInput.get()) {
    		name += list.getID() + "/";
    	}
    	if (name.length() > 1) {
    		name = name.substring(0, name.length() - 1);
    	}
    	return name;
    }

}