package beast.core.util;

import beast.core.Description;
import beast.core.parameter.DPPointer;
import beast.core.parameter.ParameterList;
import beast.core.Input;

import java.io.PrintStream;

/**
 * @author Chieh-Hsi Wu
 */
@Description("This class prints the ID numbers of a list of parameters to which a DPPointer refers.")
public class PrintPointerIDNumber extends PrintListIDNumber{

    public Input<DPPointer> pointersInput = new Input<DPPointer>("pointers", "list of items to be counted", Input.Validate.REQUIRED);

    public DPPointer pointers = null;
    public void initAndValidate(){
         pointers = pointersInput.get();

    }
    @Override
	public void init(PrintStream out) throws Exception {
        for(int i = 0; i < pointers.getDimension(); i++){
            out.print(pointers.getID() + ".idNum." + i + "\t");
        }
    }

    @Override
	public void log(int nSample, PrintStream out) {
        int dim = pointers.getDimension();
        for(int i = 0; i < dim; i++){
    	    out.print(pointers.getParameter(i).getIDNumber() + "\t");
        }
	}

}
