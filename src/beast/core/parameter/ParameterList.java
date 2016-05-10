package beast.core.parameter;

import beast.core.*;

import java.util.ArrayList;
import java.util.List;
import java.io.PrintStream;

import org.omg.CORBA.PUBLIC_MEMBER;
import org.w3c.dom.Node;

/**
 * @author Chieh-Hsi Wu
 */
@Description("This class stores a list of parameters and the size of the list can change.")
public class ParameterList extends StateNode implements PluginList, Recycle {
    public Input<List<QuietRealParameter>> parametersInput =
                new Input<List<QuietRealParameter>>(
                        "parameter",
                        "refrence to a parameter",
                        new ArrayList<QuietRealParameter>()
                );

    private ArrayList<QuietRealParameter> parameterList;
    private ArrayList<QuietRealParameter> storedParameterList;
    private int changedIndex = -1;
    private int removedIndex = -1;

    ChangeType changeType = ChangeType.ALL;
    
    public ParameterList(){
        parameterList = new ArrayList<QuietRealParameter>();
        storedParameterList = new ArrayList<QuietRealParameter>();
    }



    public void initAndValidate(){
        //System.err.println("initialize parameter list");
        List<QuietRealParameter> parameterList  = parametersInput.get();
        for(QuietRealParameter parameter: parameterList){
            parameter.setIDNumber(createID());
            this.parameterList.add(parameter);
        }
    }

    private int lastAddedIndex;
    public void addParameter(QuietRealParameter parameter){

        startEditing(null);
        parameter.setIDNumber(createID());
        lastAddedIndex = parameterList.size();
        parameterList.add(parameter);
        changeType = ChangeType.ADDED;

        //System.out.println("add parameter: "+getID()+" "+getDimension());//+" " +toString());
    }



    public void addParameter(int pIndex, QuietRealParameter parameter){
        startEditing(null);
        parameter.setIDNumber(createID());
        lastAddedIndex = pIndex;
        parameterList.add(pIndex,parameter);
        changeType = ChangeType.ADDED;

    }

    public int getLastAddedIndex(){
        //System.out.println(getID()+": "+lastAddedIndex+" "+parameterList.get(lastAddedIndex).getIDNumber());
        return lastAddedIndex;
    }

    ArrayList<Integer> idPool = new ArrayList<Integer>();
    private int newIDCount = 0;
    private int createID(){
        //Recycle id
        if(idPool.size() > 0){
            return idPool.remove(0);
        }
        
        //If no ids available to be recycled
        int newID = newIDCount;
        newIDCount++;
        return newID;
    }

    //Stores the ids of removed objects
    private void storeID(int spareID){
        idPool.add(spareID);
    }

    //Get how many ids have been created so far
    public int getObjectCreatedCount(){
        return newIDCount;
    }




    public void addParameterQuietly(QuietRealParameter parameter){
        parameterList.add(parameter);
        changeType = ChangeType.ADDED;
        //System.err.println(getID()+":added");

    }

    public void removeParameter(QuietRealParameter p){
        startEditing(null);
        //System.out.println("size: "+parameterList.size());
        storeID(p.getIDNumber());
        parameterList.remove(p);
        changeType = ChangeType.REMOVED;
        removedIndex = storedParameterList.indexOf(p);
        //System.out.println(getID()+": removed "+getDimension()+" "+toString());
    }


    public void removeParameter(int pIndex){
        startEditing(null);
        //System.out.println("size: "+parameterList.size());
        storeID(parameterList.get(pIndex).getIDNumber());
        parameterList.remove(pIndex);
        changeType = ChangeType.REMOVED;
        removedIndex = pIndex;
        //System.out.println(getID()+": removed "+getDimension()+" "+toString());
    }

    public void removeParameterQuietly(int pIndex){
        storeID(parameterList.get(pIndex).getIDNumber());
        parameterList.remove(pIndex);
        removedIndex = pIndex;
    }

    public void setValue(int pIndex, int dim, double value) {
        startEditing(null);

        parameterList.get(pIndex).setValueQuietly(dim,value);
        parameterList.get(pIndex).setEverythingDirty(true);
        changedIndex = pIndex;
        changeType = ChangeType.VALUE_CHANGED;
        //System.out.println(getID()+": changed "+changedIndex);
    }

    public void setValues(int pIndex, Double[] values) {
        startEditing(null);
        for(int i = 0; i < values.length; i++){
            parameterList.get(pIndex).setValueQuietly(i,values[i]);
        }
        parameterList.get(pIndex).setEverythingDirty(true);
        changedIndex = pIndex;
        changeType = ChangeType.VALUE_CHANGED;
        //System.out.println(getID()+": changed "+changedIndex);
    }

    public void splitParameter(int pIndex, QuietRealParameter parameter){
        startEditing(null);
        parameter.setIDNumber(createID());
        parameterList.add(parameter);
        changedIndex = pIndex;
        lastAddedIndex = parameterList.size() - 1;
        changeType = ChangeType.SPLIT;

    }

    public void splitParameter(int pIndex, int newPIndex, double newValue)throws Exception{
        startEditing(null);
        QuietRealParameter newParameter = new QuietRealParameter(new Double[]{newValue});
        newParameter.setBounds(getLower(), getUpper());

        newParameter.setIDNumber(createID());
        parameterList.add(newPIndex,newParameter);
        changedIndex = pIndex;
        lastAddedIndex = newPIndex;
        //System.out.println("newPIndexID:  "+newParameter.getIDNumber());
        changeType = ChangeType.SPLIT;

    }

    public void splitParameter(int pIndex, int newPIndex, Double[] newValue) throws Exception{
        startEditing(null);
        QuietRealParameter newParameter = new QuietRealParameter(newValue);
        newParameter.setBounds(getLower(), getUpper());

        newParameter.setIDNumber(createID());
        parameterList.add(newPIndex, newParameter);
        changedIndex = pIndex;
        lastAddedIndex = newPIndex;
        changeType = ChangeType.SPLIT;

    }

    public void splitParameter(int pIndex, double value1, int newPIndex, double value2) throws Exception{
        startEditing(null);
        QuietRealParameter newParameter = new QuietRealParameter(new Double[]{value2});
        newParameter.setBounds(getLower(), getUpper());

        newParameter.setIDNumber(createID());
        setValue(pIndex,0, value1);
        parameterList.add(newPIndex, newParameter);
        //System.out.println("newPIndexID:  "+newParameter.getIDNumber()+" "+newPIndex);
        changedIndex = pIndex;
        lastAddedIndex = newPIndex;
        changeType = ChangeType.SPLIT_AND_VALUE_CHANGE;

    }



    public void splitParameter(int pIndex, int dim, double value1, int newPIndex, Double[] values2) throws Exception{
        startEditing(null);
        QuietRealParameter newParameter = new QuietRealParameter(values2);
        newParameter.setBounds(getLower(), getUpper());
        newParameter.setIDNumber(createID());
        setValue(pIndex,dim, value1);
        parameterList.add(newPIndex, newParameter);
        changedIndex = pIndex;
        lastAddedIndex = newPIndex;
        changeType = ChangeType.SPLIT_AND_VALUE_CHANGE;

    }

    public void splitParameter(int pIndex,  Double[] values1, int newPIndex, Double[] values2) throws Exception{
        startEditing(null);
        QuietRealParameter newParameter = new QuietRealParameter(values2);
        newParameter.setBounds(getLower(), getUpper());

        newParameter.setIDNumber(createID());
        setValues(pIndex,values1);
        parameterList.add(newPIndex, newParameter);

        changedIndex = pIndex;
        lastAddedIndex = newPIndex;
        changeType = ChangeType.SPLIT_AND_VALUE_CHANGE;
    }

    /*
     * @param pIndex1 The index of the parameter to be removed
     * @param pIndex The index of the parameter that is retained.
     */
    public void mergeParameter(int pIndex1, int pIndex2){
        startEditing(null);
        storeID(parameterList.get(pIndex1).getIDNumber());
        parameterList.remove(pIndex1);
        changedIndex = pIndex2 < pIndex1? pIndex2:(pIndex2-1);
        removedIndex = pIndex1;
        changeType = ChangeType.MERGE;
    }

    public void mergeParameter(int pIndex1, int pIndex2, double newValue){
        startEditing(null);
        storeID(parameterList.get(pIndex1).getIDNumber());
        setValue(pIndex2,0,newValue);
        parameterList.remove(pIndex1);
        changedIndex = pIndex2 < pIndex1? pIndex2:(pIndex2-1);
        removedIndex = pIndex1;
        changeType = ChangeType.MERGE_AND_VALUE_CHANGE;
    }

    public void mergeParameter(int pIndex1, int pIndex2, int dim, double newValue){
        startEditing(null);
        storeID(parameterList.get(pIndex1).getIDNumber());
        setValue(pIndex2,dim,newValue);
        parameterList.remove(pIndex1);
        changedIndex = pIndex2 < pIndex1? pIndex2:(pIndex2-1);
        removedIndex = pIndex1;
        changeType = ChangeType.MERGE_AND_VALUE_CHANGE;
    }

    public void mergeParameter(int pIndex1, int pIndex2, Double[] newValues){
        startEditing(null);
        storeID(parameterList.get(pIndex1).getIDNumber());
        setValues(pIndex2, newValues);
        parameterList.remove(pIndex1);
        changedIndex = pIndex2 < pIndex1? pIndex2:(pIndex2-1);
        removedIndex = pIndex1;
        changeType = ChangeType.MERGE_AND_VALUE_CHANGE;
    }


    public int getDirtyIndex(){
        //System.out.println(getID()+changedIndex);
        return changedIndex;
    }

    public int getDirtyParameterIDNumber(){
        return parameterList.get(getDirtyIndex()).getIDNumber();
    }

    public int getParameterIDNumber(int index){

        return parameterList.get(index).getIDNumber();
    }


    public int getRemovedIndex(){
        return removedIndex;
    }

    public int getRemovedIDNumber(){
        return idPool.get(idPool.size() - 1);
    }

    public double getValue(int pIndex){
        //System.err.println(pIndex);
        return parameterList.get(pIndex).getValue();

    }

    public double getValue(int pIndex, int dim) {
        //System.out.println(pIndex+" "+dim);
        return parameterList.get(pIndex).getValue(dim);

    }

    public Double[] getValues(int pIndex){
        Double[] values = new Double[parameterList.get(pIndex).getDimension()];
        for(int i = 0;i < values.length;i++){
            values[i] = parameterList.get(pIndex).getValue(i);
        }

        return values;
    }

    public double getUpper(){
        return getParameter(0).getUpper();
    }

    public double getLower(){
        return getParameter(0).getLower();
    }

    /*
     * Get the last parameter added.
     */
    public QuietRealParameter getLastParameter(){
        return parameterList.get(parameterList.size() - 1);
    }


    public double getParameterUpper(int iParam){
        return getParameter(iParam).getUpper();
    }

    public double getParameterLower(int iParam){
        return getParameter(iParam).getLower();
    }

    

    public QuietRealParameter getParameter(int pIndex){
        return parameterList.get(pIndex);
    }

    public int getParameterDimension(){
        return getParameter(0).getDimension();
    }

    int storedNewIDCount;
    ArrayList<Integer> storedIDPool;
    protected void store(){
        storedNewIDCount = newIDCount;
        storedParameterList = new ArrayList<QuietRealParameter>();

        for(QuietRealParameter parameter:parameterList){
            parameter.store();
            storedParameterList.add(parameter);
        }
        storedIDPool = new ArrayList<Integer>();
        for(int id:idPool){
            storedIDPool.add(id);

        }
        //System.out.println("storing "+storedParameterList.size());

    }



    public void restore(){
        setEverythingDirty(false);
        //System.err.println("restore, storedListSize: "+storedParameterList.size());
        hasStartedEditing = false;
        newIDCount = storedNewIDCount;
        ArrayList<QuietRealParameter> tempList = storedParameterList;
        storedParameterList = parameterList;
        parameterList = tempList;
        for(RealParameter parameter:parameterList){
            parameter.restore();
        }
        changeType = ChangeType.ALL;
        idPool = storedIDPool;
        storedIDPool = null;


    }

    public ChangeType getChangeType(){
        //System.out.println(getID() +" "+  changeType);
        return changeType;
    }


    public ParameterList copy(){
        ParameterList copy = new ParameterList();
        for(QuietRealParameter parameter: parameterList){
            copy.addParameterQuietly(parameter.copy());
        }
        return copy;
    }

    public void setEverythingDirty(boolean isDirty){
        //System.err.println("Parameter list: "+getID()+" "+isDirty);
        setSomethingIsDirty(isDirty);
        //System.err.println("list size: "+parameterList.size());
        for(QuietRealParameter parameter : parameterList){
            parameter.setEverythingDirty(isDirty);
        }
        changeType = ChangeType.ALL;

    }

    public int indexOf(RealParameter param){
        return parameterList.indexOf(param);
    }

    public int storedIndexOf(RealParameter param){
        if(somethingIsDirty())
            return storedParameterList.indexOf(param);
        else
            return parameterList.indexOf(param);
    }



    public double getArrayValue(){
        throw new RuntimeException("As suggested by the name, this class represents a list, so it's not appropriate to return a single value");

    }

    public double getArrayValue(int dim){
        throw new RuntimeException("As suggested by the name, this class represents a list, so it's not appropriate to return a single value");

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

    /** Loggable implementation **/
    @Override
    public void log(int nSample, PrintStream out) {
        ParameterList paramList = (ParameterList) getCurrent();
        int dim = paramList.getDimension();
        if(dim == 0){
            out.print("-\t");
            return;
        }

        for(int i = 0; i < dim; i++){
            paramList.getParameter(i).log(nSample,out);
        }
        //todo
    }

    public void close(PrintStream out){

    }

    @Override
    public void init(PrintStream out) {
        int dimParam = getDimension();
        if(dimParam == 0){
            out.print(getID()+"\t");
            return;
        }

        for (int iParam = 0; iParam < dimParam; iParam++) {
            int dimValue = getParameter(iParam).getDimension();
            for(int iValue = 0; iValue < dimValue; iValue++){
                out.print(getID()+"."+iParam +"."+ iValue + "\t");
            }
        }

    }

    /** Note that changing toString means fromXML needs to be changed as well,
     * since it parses the output of toString back into a parameter.
     */
    public String toString() {
        final StringBuffer buf = new StringBuffer();
        int dimParam = getDimension();
        buf.append(getID()+"["+dimParam+"]:\n");
        for(int iParam = 0; iParam < dimParam; iParam++){
            buf.append(getParameter(iParam).getIDNumber()+" "+getParameter(iParam).toString()).append("\n");
        }
        return buf.toString();
    }    

    public int getDimension(){
        //System.out.println(getID()+ ": "+parameterList.size());
        return parameterList.size();
    }

    public int scale(double fScale) {
        if(parameterList.size() > 0){
            for(RealParameter parameter:parameterList){
                parameter.scale(fScale);
            }
            changeType = ChangeType.ALL;
            return getDimension()*getParameterDimension();
        }else{
            return 0;
        }
    }


}
