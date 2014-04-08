package beast.evolution.tree;

import beast.core.*;

/**
 * @author Chieh-Hsi Wu
 */
@Description("Extract the tree height and return it as a valuable.")
public class TreeHeightValuable extends CalculationNode implements Valuable {
    public Input<Tree> treeInput = new Input<Tree>("tree", "The tree, where the root height of which is of interest.", Input.Validate.REQUIRED);


    private Tree tree;
    public void initAndValidate(){
        tree = treeInput.get();
    }


    public int getDimension(){
        return 1;
    }

    public double getArrayValue(){
        return tree.getRoot().getHeight();
    }

    public double getArrayValue(int dim){
        if(dim < 1){
            return getArrayValue();

        }else {
            throw new RuntimeException("This is a unidimensional valuable, so please use getArrayValue() instead to retrieve the root height.");
        }
    }

}
