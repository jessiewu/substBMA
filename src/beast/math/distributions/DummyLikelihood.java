package beast.math.distributions;

import beast.core.*;


import java.util.Random;
import java.util.List;
import java.util.ArrayList;

/**
 * @author Chieh-Hsi Wu
 */
@Description("A class that is well uselful for testing things.")
public class DummyLikelihood extends Distribution {
    public Input<List<Plugin>> calcNodes = new  Input<List<Plugin>>("plugin","Some sort of input", new ArrayList<Plugin>());

	@Override
	public void initAndValidate() {}

	@Override
	public double calculateLogP() throws Exception {
		return 0.0;
	}

	@Override
	public boolean requiresRecalculation() {
		return true;
	}

	@Override public void sample(State state, Random random) {}
	@Override public List<String> getArguments() {return null;}
	@Override public List<String> getConditions() {return null;}

}
