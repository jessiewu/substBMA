package beast.evolution.alignment;

import beast.core.Description;
import beast.core.Input;

/**
 * @author Chieh-Hsi Wu
 */

@Description("The class serves as a wrapper and only returns a subset of columns defined by the user.")

public class AlignmentSubset extends Alignment{

    public Input<Alignment> alignmentInput = new Input<Alignment>("data", "full sequence data for the beast.tree", Input.Validate.REQUIRED);
    public Input<Integer> siteIndexInput = new Input<Integer>("siteIndex", "site index of the sequnence data interested", Input.Validate.REQUIRED);
    public Alignment alignment;
    public int siteIndex;

    public AlignmentSubset(){

    }
    public AlignmentSubset(Alignment alignment, int siteIndex){
        this.alignment = alignment;
        this.siteIndex = siteIndex;
        //System.err.println("siteIndex: "+siteIndex);

    }

    public void initAndValidate(){
        alignment = alignmentInput.get();
        siteIndex = siteIndexInput.get();
    }

    public int getNrTaxa() {
        return alignment.getNrTaxa();
    }

    public int getPatternCount() {
        return 1; //todo assuming a subset only has a single column
    }

    public int getTaxonIndex(String sID) {
        return alignment.getTaxonIndex(sID);
    }

    public int[] getPattern(int id) {
        return alignment.getPattern(alignment.getPatternIndex(siteIndex)); //todo assuming a subset only has a single column
    }

    public int getPattern(int iTaxon, int id) {
        //System.err.println("siteIndex: "+siteIndex);
        alignment.getPatternIndex(siteIndex);
        int[] patterns = alignment.getPattern(alignment.getPatternIndex(siteIndex));  //todo assuming a subset only has a single column
        return patterns[iTaxon];
    }

    public int getPatternWeight(int id) {
        return 1; //todo assuming a subset only has a single column
    }

    public int getMaxStateCount() {
        return alignment.getMaxStateCount();
    }

    public int getPatternIndex(int iSite) {
        throw new RuntimeException("Does not make much sense to retrive the pattern index as the colume of interest has already been defined");

    }

    public int getSiteCount() {
        return 1; //todo assuming a subset only has a single column
    }

    public int [] getWeights() {
        return new int[]{1};
    }

    public boolean[] getStateSet(int iState) {
    	return alignment.getStateSet(iState);
    }

    boolean isAmbiguousState(int state) {
    	return (state >=0 && state < getMaxStateCount());
    }
}
