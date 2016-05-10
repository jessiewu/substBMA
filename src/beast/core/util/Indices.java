package beast.core.util;

import beast.core.*;
import beast.evolution.substitutionmodel.DPNtdBMA;

import java.io.PrintStream;

/**
 * @author Chieh-Hsi Wu
 */
@Description("Reports the index of the parameter to which each of the pointer index points.")
public class Indices extends CalculationNode implements Loggable {

    public Input<DPNtdBMA> dpNtdBMAInput = new Input<DPNtdBMA>("dpNtdBMA", "a list of unique ntdBMA objects in a DP", Input.Validate.REQUIRED);
    DPNtdBMA dpNtdBMA;
    @Override
	public void initAndValidate() {
		dpNtdBMA = dpNtdBMAInput.get();
	}

    @Override
	public void init(PrintStream out) {
        int[] indices = dpNtdBMA.getPointerIndices();
        for(int i = 0; i < indices.length;i++){
            out.print("Indices("+dpNtdBMA.getID() + ")."+i+"\t");
        }

    }

    @Override
	public void log(int nSample, PrintStream out) {
        int[] indices = dpNtdBMA.getPointerIndices();
        for(int index:indices){
    	    out.print(index + "\t");
        }
	}

    @Override
	public void close(PrintStream out) {
		// nothing to do
	}
}
