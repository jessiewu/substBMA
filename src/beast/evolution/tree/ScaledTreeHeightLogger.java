package beast.evolution.tree;

import beast.core.Description;
import beast.core.Input;
import beast.core.Loggable;
import beast.evolution.tree.Scaler;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeHeightLogger;

import java.io.PrintStream;

/**
 * @author Chieh-Hsi Wu
 */
@Description("Scales the tree with the provided scaler and then log it.")
public class ScaledTreeHeightLogger extends TreeHeightLogger {
    public Input<Scaler> m_scaler = new Input<Scaler>("scaler","To scale the branch lengths by some number.", Input.Validate.REQUIRED);
    Scaler scaler;

	public void initAndValidate() {
        super.initAndValidate();
        scaler = m_scaler.get();
	}



    public void log (int nSample, PrintStream out){
        final Tree tree = treeInput.get();
		out.print(tree.getRoot().getHeight()*scaler.getScaleFactor() + "\t");
    }



}
