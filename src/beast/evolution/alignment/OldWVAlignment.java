package beast.evolution.alignment;

import beast.core.Description;
import beast.core.Input;
import beast.evolution.tree.Tree;

/**
 * @author Chieh-Hsi Wu
 */
@Description("This class should allow weight patterns to vary, but I think it is incomplete and it is not used.")
public class OldWVAlignment extends Alignment{
    public Input<Alignment> alignmentInput = new Input<Alignment>("data", "sequence data for the beast.tree", Input.Validate.REQUIRED);
    public Input<Tree> m_tree = new Input<Tree>("tree", "phylogenetic beast.tree with sequence data in the leafs", Input.Validate.REQUIRED);
        
    public boolean weightsChanged = false;
    public int[] storedWeight;
    public void initAndValidate(){
        this.alignment = alignmentInput.get();

    }

    private Alignment alignment;
    

    public OldWVAlignment(Alignment alignment){
        this(alignment,alignment.getWeights());

    }


    public OldWVAlignment(Alignment alignment, int[] weights){
        this.alignment = alignment;
        patternWeight = weights;
        storedWeight = new int[weights.length];
    }

     public int getNrTaxa() {
        return alignment.getNrTaxa();
    }

    public int getTaxonIndex(String sID) {
        return alignment.getTaxonIndex(sID);
    }

    public int getPatternCount() {
        return alignment.getPatternCount();
    }

    public int[] getPattern(int id) {
        return alignment.getPattern(id);
    }

    public int getPattern(int iTaxon, int id) {
        return alignment.getPattern(iTaxon,id);
    }

    public int getPatternWeight(int id) {
        return alignment.getPatternWeight(id);
    }

    public int getMaxStateCount() {
        return alignment.getMaxStateCount();
    }

    public int getPatternIndex(int iSite) {
        return alignment.getPatternIndex(iSite);
    }

    public int getSiteCount() {
        return alignment.getSiteCount();
    }

    public int [] getWeights() {
        return patternWeight; //TODO TO BE ONVERRIDEN!
    }

    public boolean[] getStateSet(int iState) {
    	return alignment.getStateSet(iState);
    }

    public void setWeight(int patId, int weight){
        patternWeight[patId] = weight;
        weightsChanged = true;
    }

    public void addWeight(int patId, int dweight){
        patternWeight[patId] += dweight;
        //System.out.println("patId: "+patId+", patternWeights[patId]: "+patternWeights[patId]);
        weightsChanged = true;
    }

    public void removeWeight(int patId, int dweight){
        patternWeight[patId] -= dweight;
        //System.out.println("patId: "+patId+", patternWeights[patId]: "+patternWeights[patId]);
        weightsChanged = true;
    }

    @Override
    public void store(){
        weightsChanged = false;
        System.arraycopy(patternWeight,0,storedWeight,0,patternWeight.length);
        super.store();
    }

    @Override
    public void restore(){
        weightsChanged = false;
        
        int[] temp = patternWeight;
        patternWeight = storedWeight;
        storedWeight = temp;

        super.restore();
    }

    @Override
    public boolean requiresRecalculation(){
        return weightsChanged; 
    }
   
}
