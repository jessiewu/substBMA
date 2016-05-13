package beast.evolution.sitemodel;

import beast.core.Description;
import beast.core.Input;
import beast.evolution.alignment.Alignment;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Chieh-Hsi Wu
 */
@Description("Site model that can handle multiple alignments as well as DPM.")
public abstract class DPMultiAlignSiteModel extends DPSiteModel {

    protected int[] alignmentRef;
    protected List<Alignment> alignments;


    public Input<List<Alignment>> alignmentsInput = new Input<List<Alignment>>(
            "alignment",
            "A list of alignments.",
            new ArrayList<Alignment>(),
            Input.Validate.REQUIRED
    );

    public void initAndValidate() {
        super.initAndValidate();
    }





    public abstract int getSiteModelWeight(int alignmentIndex, int category);


    public int getDirtySiteModelIndex(){
        throw new RuntimeException("Does not make sense in this case");
    }

    public int getRemovedIndex(){
         throw new RuntimeException("Retrieving the removed cluster index of a given site is not suitable in this case.");
    }

    public abstract int[] getLastDirtySites();

    public abstract int getAlignmentCount();
    public abstract int getAlignmentCategoryCount(int index);

    public abstract int getCurrCategoryIDNumber(int siteIndex);
    public abstract int getPrevCategoryIDNumber(int siteIndex);

    public abstract int getAlignmentIndex(int index);
    public abstract int getCategoryIDNumberFromIndex(int categoryIndex);

    public abstract QuietSiteModel getSiteModel(int alignmentIndex, int siteModelID);
    public QuietSiteModel getSiteModelOfSiteIndex(int siteIndex){
        return getSiteModel(getAlignmentIndex(siteIndex),getCurrCategoryIDNumber(siteIndex));
    }



}
