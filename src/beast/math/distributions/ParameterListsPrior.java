package beast.math.distributions;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.ParameterList;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Chieh-Hsi Wu
 */
@Description("This class is a wrapper that provides prior to parameter list.")
public class ParameterListsPrior extends Prior{

    public Input<List<ParameterList>> xListsInput = new Input<List<ParameterList>>(
            "xList",
            "points at which the density is calculated",
            new ArrayList<ParameterList>(),
            Input.Validate.REQUIRED
    );

    public Input<Boolean> applyToListInput = new Input<Boolean>(
            "applyToList",
            "Whether the prior is applied to the entire list",
            Input.Validate.REQUIRED
    );

    public ParameterListsPrior(){
        m_x.setRule(Input.Validate.OPTIONAL);
    }

    boolean applyToList;
    public void initAndValidate(){

        applyToList = applyToListInput.get();
        super.initAndValidate();

    }

    @Override
	public double calculateLogP() throws Exception {
        List<ParameterList> parameterLists = xListsInput.get();
        logP = ((CompoundDirichletProcess)m_dist).calcLogP(parameterLists);

		return logP;
	}

    @Override
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