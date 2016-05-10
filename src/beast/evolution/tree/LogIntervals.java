package beast.evolution.tree;

import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.Loggable;
import beast.core.parameter.RealParameter;

import java.io.PrintStream;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

/**
 * Created by IntelliJ IDEA.
 * User: cwu080
 * Date: 25/03/14
 * Time: 12:13 PM
 * To change this template use File | Settings | File Templates.
 */
public class LogIntervals extends BEASTObject implements Loggable{
    public Input<Tree> treeInput = new Input<Tree>(
            "tree",
            "the tree with intervals to be logged"
    );

    public Input<RealParameter> muInput = new Input<RealParameter>(
            "scaler",
            "scaling the tree heights later."
    );


    Tree tree;
    RealParameter mu;
    public void initAndValidate(){
        tree = treeInput.get();
        mu=muInput.get();

    }

    public void init(PrintStream out) {
        int count = tree.getInternalNodeCount();
        for(int i = 0; i < count; i++){
            out.print("treeInterval." + (i + 1) + "\t");
        }
    }

    /**
     * log this sample for current state to PrintStream,
     * e.g. value of a parameter, list of parameters or Newick tree
     *
     * @param nSample chain sample number
     * @param out     log stream
     */
    public void log(int nSample, PrintStream out){
        int count = tree.getInternalNodeCount();
        List<Node> internalNodes = tree.getInternalNodes();
        SortedSet<Double> set = new TreeSet();
        set.add(0.0);
        for(int i = 0; i < count; i++){
            set.add(internalNodes.get(i).getHeight());

        }
        Double[] heights = set.toArray(new Double[set.size()]);
        for(int i = 1; i < heights.length; i++){
            out.print((heights[i] - heights[i-1])*mu.getArrayValue()+"\t");
        }
    }

    /**
     * close log. An end of log message can be left (as in End; for Nexus trees)
     *
     * @param out log stream
     */
    public void close(PrintStream out){}
}
