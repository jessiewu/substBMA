package beast.evolution.operators;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.DPPointer;
import beast.core.parameter.DPValuable;
import beast.evolution.alignment.Alignment;
import beast.util.Randomizer;

/**
 * @author Chieh-Hsi Wu
 */
@Description("This operator swaps sites from different substitution model categories.")
public class NtdBMAPointersSwapOperator extends Operator {
    //assignment
    public Input<DPPointer> paramPointersInput = new Input<DPPointer>(
            "paramPointers",
            "array which points a set of unique parameter values",
            Input.Validate.REQUIRED
    );

    public Input<DPPointer> modelPointersInput = new Input<DPPointer>(
            "modelPointers",
            "array which points a set of unique model",
            Input.Validate.REQUIRED
    );

    public Input<DPPointer> freqPointersInput = new Input<DPPointer>(
            "freqPointers",
            "array which points a set of unique model",
            Input.Validate.REQUIRED
    );

    public Input<Boolean> checkPatternInput = new Input<Boolean>(
            "checkPattern",
            "check whether the sites swaps have identical pattern(s).",
            true
    );
    public Input<DPValuable> dpValInput = new Input<DPValuable>("dpVal", "A valuable that stores information about the category assignments", Input.Validate.REQUIRED);
    public Input<Alignment> alignmentInput = new Input<Alignment>("alignment", "A valuable that stores information about the category assignments", Input.Validate.REQUIRED);

    protected DPValuable dpVal;
    protected Alignment alignment;
    private boolean checkPattern;

    public void initAndValidate(){
        dpVal = dpValInput.get();
        alignment = alignmentInput.get();
        checkPattern = checkPatternInput.get();
    }

    public double proposal(){
        int categoryCount = dpVal.getCategoryCount();
        if(categoryCount == 1){
            return Double.NEGATIVE_INFINITY;
        }

        int categoryIndex1 = (int)(Randomizer.nextDouble()*categoryCount);
        int categoryIndex2 = categoryIndex1;
        while(categoryIndex1 == categoryIndex2){
            int something = (int)(Randomizer.nextDouble()*categoryCount);
            //System.out.println("Hi! "+categoryCount+" "+something);
            categoryIndex2 = something;
        }

        //System.out.println(categoryIndex1+" "+categoryIndex2);

        int[] sitesInCategoryIndex1 = dpVal.getClusterSites(categoryIndex1);
        int[] sitesInCategoryIndex2 = dpVal.getClusterSites(categoryIndex2);
        int site1 = sitesInCategoryIndex1[Randomizer.nextInt(sitesInCategoryIndex1.length)];
        int site2 = sitesInCategoryIndex2[Randomizer.nextInt(sitesInCategoryIndex2.length)];

        if(checkPattern && alignment.getPatternIndex(site1) == alignment.getPatternIndex(site2)){
            return Double.NEGATIVE_INFINITY;
        }


        DPPointer paramPointers = paramPointersInput.get(this);
        paramPointers.swapPointers(site1, site2);
        DPPointer modelPointers = modelPointersInput.get(this);
        modelPointers.swapPointers(site1, site2);
        DPPointer freqPointers = freqPointersInput.get(this);
        freqPointers.swapPointers(site1, site2);

        return 0.0;
    }

}
