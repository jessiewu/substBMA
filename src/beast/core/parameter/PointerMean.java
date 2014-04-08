package beast.core.parameter;

import beast.core.Description;
import beast.core.Input;
import beast.core.Loggable;
import beast.core.Valuable;
import beast.evolution.tree.Scaler;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Chieh-Hsi Wu
 */
@Description("This object returns the mean values of a vector parameters to which the DPPointer refers.")
public class PointerMean extends Scaler implements Loggable, Valuable{
    //ParameterList
    public Input<ParameterList> paramListInput = new Input<ParameterList>(
            "paramList",
            "A list of unique parameter values"
    );

    //assignment
    public Input<DPPointer> pointersInput = new Input<DPPointer>(
            "pointers",
            "An array which points a set of unique parameter values"
    );

    public Input<DPValuable> dpValInput = new Input<DPValuable>(
            "dpVal",
            "Returns the number of sites in each cluster"
    );


    public Input<List<RealParameter>> parametersInput = new Input<List<RealParameter>>(
            "parameter",
            "A list of unique parameter values",
            new ArrayList<RealParameter>()
    );

    private DPValuable dpVal;
    private DPPointer pointers;
    private ParameterList paramList;
    private List<RealParameter> parameters;

    public void initAndValidate(){
        dpVal = dpValInput.get();
        pointers = pointersInput.get();
        paramList = paramListInput.get();
        parameters = parametersInput.get();
    }



    private boolean updateMean = false;
    public boolean requiresRecalculation(){
        //System.out.println("hey");
        updateMean = true;
        return true;
    }


    private double mean;
    public double getMean(){
        int[] clusterCounts;
        //if(paramList.somethingIsDirty() || pointers.somethingIsDirty()  ){
            mean = 0.0;
            if(dpVal != null){
                clusterCounts = dpVal.getClusterCounts();
                for(int i = 0; i < clusterCounts.length;i++){
                    mean += clusterCounts[i]*paramList.getValue(i,0);
                    //System.out.println(clusterCounts[i]+" "+paramList.getValue(i,0));
                }
                mean /= pointers.getDimension();


            }else if(parameters != null){

                for(RealParameter parameter: parameters){
                    mean += parameter.getValue();
                }
                mean /= parameters.size();
            }
            updateMean = false;
        //}

        return mean;
    }

    public int getDimension(){
        return 1;
    }

    public double getArrayValue(){
        return getMean();
    }

    public double getArrayValue(int index){
        return getArrayValue();
    }




    public double getScaleFactor(){
        return getMean();
    }



    public void restore(){
        updateMean = true;
    }


    @Override
	public void init(PrintStream out) throws Exception {
        out.print(getID() +"\t");

    }

    @Override
	public void log(int nSample, PrintStream out) {
        out.print(getMean() + "\t");
	}

    @Override
	public void close(PrintStream out) {
		// nothing to do
	}
}
