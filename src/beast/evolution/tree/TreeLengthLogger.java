package beast.evolution.tree;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;
import beast.core.Loggable;

import java.io.PrintStream;

/**
 * @author Chieh-Hsi Wu
 */
@Description("Reports the length of the tree.")
public class TreeLengthLogger extends CalculationNode implements Loggable{
    public Input<Tree> treeInput = new Input<Tree>("tree", "tree to report length for.", Input.Validate.REQUIRED);
    public Input<Scaler> scalerInput = new Input<Scaler>("scaler","Provide scaleFactor to scale the tree length.");
    private Tree tree;
    private Scaler scaler;

    public void initAndValidate(){
        tree = treeInput.get();
        scaler = scalerInput.get();
    }

    public void init(PrintStream out){
        if(getID() == null){
            out.print("treeLength("+tree.getID()+")+\t");
        }else{
            out.print(getID()+"\t");
        }
    }

    public void log(int nSample, PrintStream out){
        double treeLength = calculateTreeLength();

            out.print(treeLength+"\t");

    }

    private double calculateTreeLength(){
        int nodeCount = tree.getNodeCount();
        double treeLength = 0.0;
        for(int i = 0; i < nodeCount; i++){
            treeLength += tree.getNode(i).getLength();
        }
        if(scaler != null){
            //make sure that the scaler is updated!
            requiresRecalculation();
            treeLength = treeLength*scaler.getScaleFactor();
        }
        return treeLength;
    }

    public void close(PrintStream out){}





}
