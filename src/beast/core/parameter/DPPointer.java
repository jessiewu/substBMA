package beast.core.parameter;

import beast.core.Input;
import beast.core.StateNode;
import beast.core.Description;

import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;
import java.io.PrintStream;

import org.w3c.dom.Node;

/**
 * @author Chieh-Hsi Wu
 */
@Description("Array that points to some set of parameters")
public class DPPointer extends StateNode {
    public Input<List<QuietRealParameter>> uniqueParametersInput = new Input<List<QuietRealParameter>>(
            "uniqueParameter",
            "refrence to a parameter",
            new ArrayList<QuietRealParameter>(),
            Input.Validate.REQUIRED
    );

    public Input<IntegerParameter> assignmentInput = new Input<IntegerParameter>(
            "initialAssignment",
            "a parameter which specifies the assignment of elements to clusters",
            Input.Validate.REQUIRED
    );





    protected QuietRealParameter[] parameters;
    private QuietRealParameter[] storedParameters;
    protected int lastDirty = -1;
    private int[] lastSwappedSites;
    private int[] lastDirtySites;
    private int storedLastDirty = -1;
    protected ChangeType changeType;

    public DPPointer(){
        
    }
    public DPPointer(int dim){
        parameters = new QuietRealParameter[dim];
        storedParameters = new QuietRealParameter[dim];
        lastSwappedSites = new int[2];

    }

    public void initAndValidate(){
        IntegerParameter initialAssignment = assignmentInput.get();
        List<QuietRealParameter> uniqueParameters = uniqueParametersInput.get();
        parameters = new QuietRealParameter[initialAssignment.getDimension()];
        storedParameters = new QuietRealParameter[initialAssignment.getDimension()];
        for(int i = 0; i < parameters.length;i++){
            parameters[i] = uniqueParameters.get(initialAssignment.getValue(i));
        }
        System.arraycopy(parameters,0,storedParameters,0,parameters.length);
        lastSwappedSites = new int[2];
    }

    public void point(int dim, QuietRealParameter parameter){
        //System.out.println("pointer: "+dim);
        startEditing(null);
        lastDirty = dim;
        parameters[dim] = parameter;
        changeType = ChangeType.POINTER_CHANGED;

    }

    public void pointQuitely(int dim, QuietRealParameter parameter){
        parameters[dim] = parameter;

    }

    public void swapPointers(int siteIndex1, int siteIndex2){
        startEditing(null);
        lastSwappedSites[0] = siteIndex1;
        lastSwappedSites[1] = siteIndex2;
        QuietRealParameter temp = parameters[siteIndex1];
        parameters[siteIndex1] = parameters[siteIndex2];
        parameters[siteIndex2] = temp;
        changeType = ChangeType.POINTERS_SWAPPED;

    }

    public void multiPointerChanges(int dim, QuietRealParameter parameter){
        //System.out.println("pointer: "+dim);
        startEditing(null);
        lastDirty = dim;
        parameters[dim] = parameter;
        changeType = ChangeType.MULTIPLE_POINTERS_CHANGED;

    }

    public void multiPointerChanges(int[] fromSites, QuietRealParameter parameter){
        startEditing(null);

        for(int fromSite: fromSites){

            //if(fromSite == toSite) throw new RuntimeException("Same site! fromSite: "+fromSite+" "+toSite);
            parameters[fromSite] = parameter;
        }
        //System.out.println();
        lastDirtySites = new int[fromSites.length];
        System.arraycopy(fromSites,0,lastDirtySites,0,fromSites.length);
        changeType = ChangeType.MULTIPLE_POINTERS_CHANGED;
    }

    public void multiPointerChanges(int[] fromSites, int toSite){
        startEditing(null);

        for(int fromSite: fromSites){

            //if(fromSite == toSite) throw new RuntimeException("Same site! fromSite: "+fromSite+" "+toSite);
            parameters[fromSite] = parameters[toSite];
        }
        //System.out.println();
        lastDirtySites = new int[fromSites.length];
        System.arraycopy(fromSites,0,lastDirtySites,0,fromSites.length);
        changeType = ChangeType.MULTIPLE_POINTERS_CHANGED;
    }

    public void multiPointerChanges(int[] fromSites, int[] toSites){
        startEditing(null);

        for(int i = 0; i < fromSites.length; i++){
            parameters[fromSites[i]] = storedParameters[toSites[i]];
        }
        //System.out.println();
        lastDirtySites = new int[fromSites.length];
        System.arraycopy(fromSites,0,lastDirtySites,0,fromSites.length);
        changeType = ChangeType.MULTIPLE_POINTERS_CHANGED;
    }


    public void multiPointerChanges(int[][] fromSites, QuietRealParameter[] ref){
        startEditing(null);
        //System.out.println(fromSites.length);


        for(int i = 0; i < fromSites.length; i++){
            for(int j  = 0; j < fromSites[i].length; j++){
                parameters[fromSites[i][j]] = ref[i];
                //System.out.print(fromSites[i][j]+" ");
            }

            //System.out.println();
            //System.out.println(fromSites[i].length);
        }
        //System.out.println();

        int totalLength = 0;
        for(int[] array:fromSites){
            totalLength += array.length;
        }
        //System.out.println("totalLength: "+totalLength);
        lastDirtySites = new int[totalLength];
        int k = 0;
        for(int i = 0; i < fromSites.length; i++){
            System.arraycopy(fromSites[i],0,lastDirtySites,k,fromSites[i].length);
            k += fromSites[i].length;
            //System.out.println("k: "+k);


        }
        /*System.out.println("----------------");
        for(int i =0; i < lastDirtySites.length; i++){
            System.out.print(lastDirtySites[i]+" ");
        }
        System.out.println();
        System.out.println("----------------");   */
        changeType = ChangeType.MULTIPLE_POINTERS_CHANGED;
    }

    public void multiPointerChanges(int[][] fromSites, int toSite[]){
        startEditing(null);
        QuietRealParameter[] toSitePointers = new QuietRealParameter[toSite.length];
        for(int i = 0; i < toSitePointers.length; i++){
            toSitePointers[i] = parameters[toSite[i]];
        }

        for(int i = 0; i < fromSites.length; i++){
            for(int j  = 0; j < fromSites[i].length; j++){
                parameters[fromSites[i][j]] = toSitePointers[i];
            }
        }
        //System.out.println();

        int totalLength = 0;
        for(int[] array:fromSites){
            totalLength += array.length;
        }
        lastDirtySites = new int[totalLength];
        int k = 0;
        for(int[] fromSitesArray:fromSites){
            System.arraycopy(fromSitesArray,0,lastDirtySites,k,fromSitesArray.length);
            k = fromSitesArray.length;

        }

        changeType = ChangeType.MULTIPLE_POINTERS_CHANGED;
    }

    public ChangeType getChangeType(){
        return changeType;
    }

    public int[] getSwappedSites(){
        return Arrays.copyOf(lastSwappedSites, lastSwappedSites.length);
    }

    public int[] getLastDirtySites(){
        return Arrays.copyOf(lastDirtySites, lastDirtySites.length);
    }



    public int getDimension(){
        return parameters.length;
    }

    public DPPointer copy(){
        DPPointer copy = new DPPointer(getDimension());
        for(int i = 0; i < copy.getDimension(); i++){
            copy.pointQuitely(i, parameters[i]);
        }
        return copy;
    }

    public void setEverythingDirty(boolean isDirty){
        setSomethingIsDirty(isDirty);
        //System.err.println("isDirty: "+isDirty);
    }

    public boolean isParameterDirty(int index){
        return parameters[index].somethingIsDirty();

    }

    public void store(){
        //System.out.println("store?");
        storedLastDirty = lastDirty;
        System.arraycopy(parameters,0,storedParameters,0,parameters.length);
        /*for(int i = 0; i < storedParameters.length; i++){
            System.out.println(getID()+" stored param "+i+": "+storedParameters[i]);
        }*/
    }

    @Override
	public void restore() {
        setEverythingDirty(false);
        //System.err.println("restore?");
        lastDirty = storedLastDirty;
        QuietRealParameter[] temp = parameters;
        parameters = storedParameters;
        storedParameters = temp;
	}

    public int storedIndexInList(int index, ParameterList paramList){
        /*for(int i = 0; i < storedParameters.length; i++){
            System.out.println(getID()+" param "+i+": "+storedParameters[i]);
        }
        for(int i = 0; i < paramList.getDimension(); i++){
            System.out.println("list:"+(paramList.getParameter(i)));
        }
        System.out.println("target: "+storedParameters[index]);*/
        //System.out.println("Hi!");
        //System.out.println();
        return paramList.storedIndexOf(storedParameters[index]);

    }



    public int scale(double fScale){
        throw new RuntimeException("Scaling simply does not make sense in this case");

    }

    public int getLastDirty(){
        //System.out.println("getLastDirty: "+lastDirty);
        return lastDirty;
    }

    public double getArrayValue(){
        return -1;

    }

    public double getArrayValue(int dim){
        return -1;
    }

        /** other := this
     *  Assign all values of this to other **/
    public void assignTo(StateNode other){
        //todo

    }

    /** this := other
     * Assign all values of other to this
     *
     **/
    public void assignFrom(StateNode other){
        //todo

    }

    /** As assignFrom, but only those parts are assigned that
     * are variable, for instance for parameters bounds and dimension
     * do not need to be copied.
     */
    public void assignFromFragile(StateNode other){
        //todo
    }

    /** StateNode implementation **/
    @Override
    public void fromXML(Node node) {
        //todo
    }

    public QuietRealParameter getParameter(int dim){
        return parameters[dim];

    }

    public double getParameterValue(int dim){
        return parameters[dim].getValue();

    }

    public int indexInList(int dim, ParameterList paramList){
        /*System.err.println("parameter: "+parameters[dim]);
        System.err.println("parameterList: "+paramList.getDimension());
        for(int i = 0; i < paramList.getDimension();i++){
            System.err.print(paramList.getParameter(i)+", ");
            System.err.println(paramList.getParameter(i) == parameters[dim]);
        }
        System.err.println(paramList.indexOf(parameters[dim]));*/
        //System.out.println("dim: "+dim );
        return paramList.indexOf(parameters[dim]);
    }

    public int getParameterIDNumber(int siteIndex){
        return parameters[siteIndex].getIDNumber();

    }

    public int getStoredParameterIDNumber(int siteIndex){
        return storedParameters[siteIndex].getIDNumber();

    }


    public boolean sameParameter(int dim, RealParameter parameter){
        return parameters[dim].getValue() == parameter.getValue();
    }

        /** Loggable implementation **/
    @Override
    public void log(int nSample, PrintStream out) {
        DPPointer dPPointer = (DPPointer) getCurrent();
        int dim = dPPointer.getDimension();
        for(int i = 0; i < dim; i++){
            parameters[i].log(nSample,out);
        }
        //todo
    }

    public void close(PrintStream out){

    }

    @Override
    public void init(PrintStream out) {
        int dimParam = getDimension();
        for (int iParam = 0; iParam < dimParam; iParam++) {
            int dimValue = getParameter(iParam).getDimension();
            for(int iValue = 0; iValue < dimValue; iValue++){
                out.print(getID()+"."+iParam +"."+ iValue + "\t");
            }
        }

    }

    public boolean pointEqual(int index1, int index2){
        return parameters[index1] == parameters[index2];
    }

    public boolean pointEqual(int index1, RealParameter param){
        return parameters[index1] == param;
    }

    public String toString(){
        final StringBuffer buf = new StringBuffer();
        int dimParam = getDimension();
        buf.append(getID()+"["+dimParam+"]:\n");
        for(int iParam = 0; iParam < dimParam; iParam++){
            buf.append(getParameter(iParam).getIDNumber()).append("\t");
        }
        return buf.toString();
    }
}
