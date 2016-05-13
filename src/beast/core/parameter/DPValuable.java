package beast.core.parameter;

import beast.core.*;

import java.io.PrintStream;
import java.util.Arrays;

/**
 * @author Chieh-Hsi Wu
 */

@Description("Summarizes a set of valuables from a Dirichlet prior process")
public class DPValuable extends CalculationNode implements Function, Loggable{

    //ParameterList
    public Input<ParameterList> paramListInput = new Input<ParameterList>(
            "paramList",
            "A list of unique parameter values",
            Input.Validate.REQUIRED
    );


    //assignment
    public Input<DPPointer> pointersInput = new Input<DPPointer>(
            "pointers",
            "array which points a set of unique parameter values",
            Input.Validate.REQUIRED
    );

    private ParameterList paramList;
    private DPPointer pointers;
    private boolean pointersChanged;
    private int[] clusterCounts;
    private int[] storedClusterCounts;
    private int[][] storedClusterSites;


    public void initAndValidate(){
        paramList = paramListInput.get();
        pointers = pointersInput.get();
        pointersChanged = true;
        getClusterCounts();
    }

    public DPValuable(){
    }

    public DPValuable(ParameterList paramList, DPPointer pointers){
        this.paramList = paramList;
        this.pointers = pointers;
        pointersChanged = true;
    }

    public String printParameterList(){
        return paramList.toString();
    }


    /** CalculationNode methods **/
    // smarter vesions to be computed
	@Override
	public void store() {

        storedClusterCounts = new int[clusterCounts.length];
        System.arraycopy(clusterCounts,0,storedClusterCounts,0,clusterCounts.length);
        storedClusterSites = new int[clusterCounts.length][];
        for(int i = 0;i < clusterCounts.length;i++){
            storedClusterSites[i] = new int[clusterCounts[i]];
            System.arraycopy(clusterSites[i],0,storedClusterSites[i],0,clusterCounts[i]);
        }
		super.store();
	}
	@Override
	public void restore() {
        int[] tmp1 = clusterCounts;
        clusterCounts = storedClusterCounts;
        storedClusterCounts = tmp1;

        int[][] tmp2 = clusterSites;
        clusterSites = storedClusterSites;
        storedClusterSites = tmp2;


        super.restore();
	}

    @Override
	public boolean requiresRecalculation() {
        //System.out.println(pointers.getID()+" Hello?? "+ pointers.somethingIsDirty());
        pointersChanged =pointers.somethingIsDirty();
        if(pointers.somethingIsDirty()){
            update();
        }
		return pointers.somethingIsDirty();
	}

    /** Valuable implementation follows **/
	@Override
	public int getDimension() {
		return paramList.getDimension();
	}

    @Override
	public double getArrayValue() {
		return clusterCounts[0];
	}

    public double getArrayValue(int dim){
        return clusterCounts[dim];
    }

    public int getSiteID(int index){
        return pointers.getParameterIDNumber(index);
    }

    public int getOneClusterSite(int clusterIndex){
        return clusterSites[clusterIndex][0];
    }

    public int[] getLastDirtySites(){
        return pointers.getLastDirtySites();
    }

    /*public int getCategoryID(int index){
        return paramList.getParameter(index).getIDNumber();

    }*/

    private int[][] clusterSites;
    public int[] getClusterSites(int index){
        /*if(pointersChanged){
            update();
        }*/

        
        int[] clusterSites = new int[clusterCounts[index]];
        System.arraycopy(this.clusterSites[index],0,clusterSites,0,clusterSites.length);
        return clusterSites;
    }

    public int[] getMergedSites(){
        int index = paramList.getDirtyIndex();
        int[] storedClusterSites = new int[storedClusterCounts[index]];
        System.arraycopy(this.storedClusterSites[index],0,storedClusterSites,0,storedClusterSites.length);
        return storedClusterSites;

    }

    public int[] getStoredClusterSites(int index){
        int[] storedClusterSites = new int[storedClusterCounts[index]];
        System.arraycopy(this.storedClusterSites[index],0,storedClusterSites,0,storedClusterSites.length);
        return storedClusterSites;
    }

    public int[] getClusterCounts(){
        if(pointersChanged){
            update();
        }

        return Arrays.copyOf(clusterCounts, clusterCounts.length);
    }

    public int getClusterSize(int clusterIndex){
        if(pointersChanged){
            update();
        }

        return clusterCounts[clusterIndex];
    }

    public int getCategoryCount(){
        return paramList.getDimension();
    }

    public int getCategoryIDNumber(int index){
        return paramList.getParameter(index).getIDNumber();
    }

    public int getRemovedCategoryIDNumber(){
        return paramList.getRemovedIDNumber();
    }



    public void update(){
        //System.out.println("dpval update");
        clusterCounts = new int[paramList.getDimension()];

        clusterSites = new int[paramList.getDimension()][pointers.getDimension()];
        int index;
        for(int i = 0; i < pointers.getDimension();i++ ){
            index = pointers.indexInList(i,paramList);
            if(index == -1){
                System.out.println(i+" "+pointers.getParameter(i)+" "+paramList);
            }
            clusterSites[index][clusterCounts[index]] = i;
            clusterCounts[index]++;
        }

        /*for(int i = 0; i < clusterSites.length;i++){
            for(int j = 0; j < clusterSites[i].length;j++){
                System.out.print(clusterSites[i][j]+" ");
            }
            System.out.println();
        }*/
        pointersChanged = false;
    }

    public int getPrevCategory(int pointerIndex){
        return pointers.storedIndexInList(pointerIndex,paramList);
    }

    public int getCurrCategory(int pointerIndex){
        return pointers.indexInList(pointerIndex,paramList);
    }

    public int getPrevCategoryIDNumber(int pointerIndex){
        return pointers.getStoredParameterIDNumber(pointerIndex);
    }

    public int getCurrCategoryIDNumber(int pointerIndex){
        return pointers.getParameterIDNumber(pointerIndex);
    }

    public int getPointerDimension(){
        return pointers.getDimension();
    }

     public void init(PrintStream out) {
        out.print(getID()+"\t");

    }

    public void log(int nSample, PrintStream out){
        //
        /*System.out.println(clusterSites.length);
        String s ="{";
        for(int i =0; i < clusterSites.length;i++){
            s+=clusterSites[i].length+": ";
            for(int j = 0; j < clusterCounts[i]; j++){
                s+=clusterSites[i][j];
                if(j == clusterSites[i].length - 1){
                    s+=";";
                } else{
                    s+=",";
                }
            }

        }
        System.out.println(s);*/
        out.print(paramList.getDimension() + "\t");
    }

    public void close(PrintStream out){
    }


}
