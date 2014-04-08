package beast.evolution.substitutionmodel;

import beast.core.*;
import beast.core.parameter.*;
import beast.evolution.sitemodel.SiteModel;

import java.util.ArrayList;

/**
 * @author Chieh-Hsi Wu
 */
@Description("A list of SwitchingNtdBMA classes for DPP.")
public class DPNtdBMA extends CalculationNode implements PluginList, Recycle {
    ArrayList<SwitchingNtdBMA> ntdBMAs = new ArrayList<SwitchingNtdBMA>();
    ArrayList<SwitchingNtdBMA> storedNtdBMAs = new ArrayList<SwitchingNtdBMA>();


    public Input<ParameterList> paramListInput = new Input<ParameterList>(
            "paramList",
            "A list of unique parameter values",
            Input.Validate.REQUIRED
    );

    public Input<ParameterList> freqsListInput = new Input<ParameterList>(
            "freqsList",
            "A list of unique parameter values",
            Input.Validate.REQUIRED
    );

    public Input<ParameterList> modelListInput = new Input<ParameterList>(
            "modelList",
            "A list of unique model values",
            Input.Validate.REQUIRED
    );




    public Input<DPPointer> pointersInput = new Input<DPPointer>(
            "pointers",
            "array which points a set of unique model",
            Input.Validate.REQUIRED
    );


    protected ParameterList paramList;
    protected ParameterList modelList;
    protected ParameterList freqsList;
    protected DPPointer pointers;
    //private DPPointer modelPointers;
    //private DPPointer freqPointers;
    ArrayList<SwitchingNtdBMA> ntdBMAPool = new ArrayList<SwitchingNtdBMA>();

    protected int[] pointerIndices;
    protected int removedIndex = -1;

    protected ChangeType changeType;
    public void initAndValidate(){
        paramList = paramListInput.get();
        modelList = modelListInput.get();
        freqsList = freqsListInput.get();
        pointers = pointersInput.get();

        int dimParamList = paramList.getDimension();
        for(int i = 0; i < dimParamList;i++){
            ntdBMAs.add(
                    createSwitchingNtdBMA(
                            paramList.getParameter(i),
                            modelList.getParameter(i),
                            freqsList.getParameter(i)
                    )
            );

        }
        //System.err.println("model counts: "+ntdBMAs.size());
        pointerIndices = new int[pointers.getDimension()];
        resetAllPointerIndices();
    }

    public int getDimension(){
        return ntdBMAs.size();
    }

    public int getCount(){
        return pointers.getDimension();
    }

    public SwitchingNtdBMA getModel(int index){
        return ntdBMAs.get(index);
    }

    public int getLastAddedIndex(){
        return freqsList.getLastAddedIndex();
    }



    public int[] getPointerIndices(){
        return pointerIndices;
    }

    public int getDirtyModelIndex(){
        if(paramList.somethingIsDirty()){
            return paramList.getDirtyIndex();
        }else if(modelList.somethingIsDirty()){
            return modelList.getDirtyIndex();
        }else {
            return freqsList.getDirtyIndex();

        }
    }

    public int getDirtyModelIDNumber(){
        return ntdBMAs.get(getDirtyModelIndex()).getIDNumber();
    }

    public int getModelBySiteIndexIDNumber(int index){

        return pointers.getParameterIDNumber(index);
    }


    public int getModelByListIndexIDNumber(int index){

        return ntdBMAs.get(index).getIDNumber();
    }



    protected void addModel(){
        /*QuietRealParameter parameter = paramList.getParameter(paramList.getDimension()-1);
        QuietRealParameter model = modelList.getParameter(modelList.getDimension()-1);
        QuietRealParameter freqs = freqsList.getParameter(freqsList.getDimension()-1); */
        QuietRealParameter parameter = paramList.getParameter(paramList.getLastAddedIndex());
        QuietRealParameter model = modelList.getParameter(modelList.getLastAddedIndex());
        QuietRealParameter freqs = freqsList.getParameter(freqsList.getLastAddedIndex());

        /*if(paramList.getDimension()-1 != paramList.getLastAddedIndex()){

            throw new RuntimeException(paramList.getChangeType()+" "+(paramList.getDimension()-1)+" "+paramList.getLastAddedIndex());
        }

        if(modelList.getDimension()-1 != modelList.getLastAddedIndex()){
            throw new RuntimeException((modelList.getDimension()-1)+" "+modelList.getLastAddedIndex());
        }

        if(freqsList.getDimension()-1 != freqsList.getLastAddedIndex()){
            throw new RuntimeException((freqsList.getDimension()-1)+" "+freqsList.getLastAddedIndex());
        } */


        ntdBMAs.add(paramList.getLastAddedIndex(),createSwitchingNtdBMA(parameter,model,freqs));


    }

    protected void removeModel(int index){
        //System.out.println("removeNtdModel: "+index+" list size: "+ntdBMAs.size());
        pushPool(ntdBMAs.remove(index));
        removedIndex = index;
        //System.out.println("removeNtdModel: "+index+" list size: "+ntdBMAs.size());
    }

    private void pushPool(SwitchingNtdBMA ntdBMA){
        ntdBMAPool.add(ntdBMA);
    }

    private SwitchingNtdBMA popPool(){
        if(ntdBMAPool.size()>0){
            return ntdBMAPool.remove(0);
        }else{
            return null;
        }
    }

    public int getLastDirtySite(){
        return pointers.getLastDirty();
    }

    public int[] getLastDirtySites(){
        return pointers.getLastDirtySites();
    }

    /*
     * Returns the last added model
     */
    public SwitchingNtdBMA getLastModel(){
        return ntdBMAs.get(ntdBMAs.size() - 1);
    }


    public int getRemovedIndex(){
        return removedIndex;
    }

    public int getRemovedModelIDNumber(){
        //System.out.println(storedNtdBMAs.get(removedIndex).getIDNumber());
        return storedNtdBMAs.get(removedIndex).getIDNumber();
    }

    public ChangeType getChangeType(){
        return changeType;
    }

    int prevCluster = -1;
    public void setupPointerIndices(){
        //System.out.println("changeType: "+changeType);
        if(changeType == ChangeType.ADDED || changeType == ChangeType.POINTER_CHANGED){
            int changedIndex = pointers.getLastDirty();
            prevCluster = pointerIndices[changedIndex];
            pointerIndices[changedIndex] = pointers.indexInList(changedIndex,freqsList);
        }else if (changeType == changeType.POINTERS_SWAPPED){
            int[] changedSites = pointers.getSwappedSites();
            pointerIndices[changedSites[0]] = pointers.indexInList(changedSites[0],freqsList);
            pointerIndices[changedSites[1]] = pointers.indexInList(changedSites[1],freqsList);
        }else if (changeType == ChangeType.REMOVED){
            resetAllPointerIndices();
        }else if (changeType == ChangeType.ALL){
            resetAllPointerIndices();
        }
    }

    public int[] getSwappedSites(){
        return pointers.getSwappedSites();
    }

    public int getPrevCluster(int index){
        return pointers.storedIndexInList(index,freqsList);
    }

    public int getCurrCluster(int index){
        return pointers.indexInList(index,freqsList);
    }

    public int getCurrModelIDNumber(int index){
        return pointers.getParameterIDNumber(index);
    }

    public int getPrevModelIDNumber(int index){
        return pointers.getStoredParameterIDNumber(index);
    }

 

    public void resetAllPointerIndices(){
        for(int i = 0; i < pointerIndices.length;i++){
            pointerIndices[i] = pointers.indexInList(i,freqsList);
        }
    }

    int newModelCreatedCount = 0;
    int storedNewModelCreatedCount;
    public SwitchingNtdBMA createSwitchingNtdBMA(
            QuietRealParameter modelParameters,
            QuietRealParameter modelCode,
            QuietRealParameter frequencies){
        QuietRealParameter logKappa = new RealParameterWrapper(modelParameters, 0);
        QuietRealParameter logTN = new RealParameterWrapper(modelParameters,1);
        QuietRealParameter logAC = new RealParameterWrapper(modelParameters,2);
        QuietRealParameter logAT = new RealParameterWrapper(modelParameters,3);
        QuietRealParameter logGC = new RealParameterWrapper(modelParameters,4);
        //RealParameter logGT = new RealParameterWrapper(modelParameters,5);
       
        //System.err.println();
        SwitchingNtdBMA recycled = popPool();
        if(recycled == null){
            SwitchingNtdBMA  ntdBMA = new SwitchingNtdBMA(
                logKappa,
                logTN,
                logAC,
                logAT,
                logGC,
                //logGT,
                modelCode,
                frequencies
            );
            ntdBMA.setIDNumber(newModelCreatedCount);
            newModelCreatedCount++;
            return ntdBMA;
        }else{
            recycled.initialize(
                    logKappa,
                    logTN,
                    logAC,
                    logAT,
                    logGC,
                    //logGT,
                    modelCode,
                    frequencies
            );
            return recycled;
        }

    }

    public int getObjectCreatedCount(){
        return newModelCreatedCount;
    }


    ArrayList storedPool;
    public void store(){
        storedNewModelCreatedCount = newModelCreatedCount;
        storedNtdBMAs = new ArrayList<SwitchingNtdBMA>();

        for(SwitchingNtdBMA ntdBMA:ntdBMAs){
            //System.out.println("stored ntdBMA: "+ntdBMA.getIDNumber());
            storedNtdBMAs.add(ntdBMA);
            ntdBMA.store();
        }

        storedPool = new ArrayList<SwitchingNtdBMA>();
        for(SwitchingNtdBMA ntdBMA:ntdBMAPool){
            //System.out.println("pool ntdBMA: "+ntdBMA.getIDNumber());
            storedPool.add(ntdBMA);
        }
        prevCluster = -1;
        super.store();
    }

    public void restore(){
        newModelCreatedCount = storedNewModelCreatedCount;
        ntdBMAs = storedNtdBMAs;
        storedNtdBMAs = null;
        for(SwitchingNtdBMA ntdBMA:ntdBMAs){
            ntdBMA.restore();
        }
        ntdBMAPool = storedPool;
        storedPool = null;
        prevCluster = -1;
        super.restore();
        //System.out.println("ntd restoring");
    }



    public boolean requiresRecalculation(){

        boolean recalculate = false;
        //System.err.println("dirty0");
        if(paramList.somethingIsDirty()||modelList.somethingIsDirty()||freqsList.somethingIsDirty()){
            //System.out.println("dirty0: "+freqsList.getChangeType());
            //ChangeType changeType;
            if(paramList.somethingIsDirty()){
                changeType = paramList.getChangeType();

            }else if(modelList.somethingIsDirty()){
                changeType = modelList.getChangeType();

            }else{
                changeType = freqsList.getChangeType();
            }

            //System.out.println("dirty0: "+changeType);
            if(changeType == ChangeType.ADDED){
                //System.out.println(getID()+"model added");
                addModel();
                //setupPointerIndices();

            }else if(changeType == ChangeType.REMOVED){
                //System.out.println(getID()+"model removed, "+paramList.getRemovedIndex());
                removeModel(freqsList.getRemovedIndex());
                //setupPointerIndices();

            }else if(changeType == ChangeType.VALUE_CHANGED){
                //System.out.println(getID()+": model changed");
                for(SwitchingNtdBMA ntdBMA:ntdBMAs){

                    MCMCNodeFactory.checkDirtiness(ntdBMA);
                    //System.out.println(ntdBMA.getIDNumber()+" is dirty? "+ntdBMA.isDirtyCalculation());
                }
            }else if(changeType == ChangeType.SPLIT){
                addModel();
            }else if(changeType == ChangeType.SPLIT_AND_VALUE_CHANGE){

                addModel();
                MCMCNodeFactory.checkDirtiness(ntdBMAs.get(freqsList.getDirtyIndex()));
            }else if(changeType == ChangeType.MERGE){
                removeModel(freqsList.getRemovedIndex());
            }else if(changeType == ChangeType.MERGE_AND_VALUE_CHANGE){

                removeModel(freqsList.getRemovedIndex());
                MCMCNodeFactory.checkDirtiness(ntdBMAs.get(freqsList.getDirtyIndex()));

            }else{
                this.changeType = ChangeType.ALL;
                setupPointerIndices();
                for(SwitchingNtdBMA ntdBMA:ntdBMAs){
                    MCMCNodeFactory.checkDirtiness(ntdBMA);
                    //System.out.println(ntdBMA.getIDNumber()+" is dirty? "+ntdBMA.isDirtyCalculation());
                }
            }
            recalculate = true;
            //System.err.println("dirty1");

        }else if (pointers.somethingIsDirty()){
            recalculate = true;
            //System.out.println(getID()+": pointer changed");
            changeType = pointers.getChangeType();

            //setupPointerIndices();
        }

        //System.out.println("recalculate: "+recalculate+" "+changeType);

        return recalculate;
    }

    protected void accept() {
        for(SwitchingNtdBMA ntdBMA:ntdBMAs){
            ntdBMA.makeAccept();
        }
    	super.accept();
    }




}
