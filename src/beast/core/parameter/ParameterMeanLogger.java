package beast.core.parameter;

import beast.core.*;
import beast.evolution.tree.TreeHeightValuable;

import java.io.PrintStream;
import java.util.logging.StreamHandler;

/**
 * @author Chieh-Hsi Wu
 */
@Description("Use this class to print the mean value of a ParameterList.")
public class ParameterMeanLogger extends BEASTObject implements Loggable {
    //ParameterList
    public Input<ParameterList> paramListInput = new Input<ParameterList>(
            "paramList",
            "A list of unique parameter values",
            Input.Validate.REQUIRED
    );

    //assignment
    public Input<DPPointer> pointersInput = new Input<DPPointer>(
            "pointers",
            "An array which points a set of unique parameter values",
            Input.Validate.REQUIRED
    );

    public Input<DPValuable> dpValInput = new Input<DPValuable>(
            "dpVal",
            "Returns the number of sites in each cluster",
            Input.Validate.REQUIRED
    );

    public Input<Function> treeHeightInput = new Input<Function>(
            "tree",
            "Returns the number of sites in each cluster",
            Input.Validate.REQUIRED
    );

    public Input<Boolean> divideInput = new Input<Boolean>(
            "divide",
            "Whether to divide the valuable by the scaler",
            false
    );

    private DPValuable dpVal;
    private DPPointer pointers;
    private ParameterList paramList;
    private Function treeHeight;
    private boolean divide;

    public void initAndValidate(){
        dpVal = dpValInput.get();
        pointers = pointersInput.get();
        paramList = paramListInput.get();
        treeHeight = treeHeightInput.get();
        divide = divideInput.get();
    }

    public void init(PrintStream out){
        //System.out.println(getID());
        int dim = treeHeight.getDimension();
        for(int i = 0; i < dim; i++)
        out.print(getID()+(i+1)+"\t");
    }


    public void log(int nSample, PrintStream out){

        double mean = 0.0;
        int[] clusterCounts = dpVal.getClusterCounts();
        for(int i = 0; i < clusterCounts.length;i++){
            mean += clusterCounts[i]*paramList.getValue(i,0);
        }
        mean /= pointers.getDimension();
        int dim = treeHeight.getDimension();
        if(divide){
            for(int i = 0; i < dim; i++){
                out.print(treeHeight.getArrayValue(i)/mean+"\t");
            }


        }else{
            for(int i = 0; i < dim; i++){
                out.print(mean*treeHeight.getArrayValue(i)+"\t");
            }
        }

    }

    public void close(PrintStream out){}
}
