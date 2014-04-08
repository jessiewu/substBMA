package beast.evolution.alignment;

import beast.core.Input;
import beast.core.Loggable;
import beast.core.Plugin;
import beast.core.Valuable;

import java.io.PrintStream;

/**
 * @author Chieh-Hsi Wu
 */
public class MultinomialTestStatistic extends Plugin implements Loggable {
    public Input<Alignment> alignmentInput = new Input<Alignment>(
            "alignment",
            "The alignment which its multinomial statistic is computed.",
            Input.Validate.REQUIRED
    );

    private Alignment alignment;
    private double multinomialTestStatistic;
    public void initAndValidate(){
        alignment = alignmentInput.get();
        computeMultinomalTestStatistic();
    }

    public void computeMultinomalTestStatistic(){
        int[] patternWeights = new int[alignment.getPatternCount()];
        multinomialTestStatistic = 0.0;
        for(int i = 0; i < patternWeights.length; i++){
            patternWeights[i]=alignment.getPatternWeight(i);
            multinomialTestStatistic += patternWeights[i]*Math.log(patternWeights[i]);
            //System.out.println(patternWeights[i]+" "+Math.log(patternWeights[i]));
        }
        int siteCount = alignment.getSiteCount();
        multinomialTestStatistic -= siteCount*Math.log(siteCount);


    }

    public double getMultinomialTestStatistic(){
        return multinomialTestStatistic;

    }




    /**
     * Loggable implementation *
     */

    @Override
    public void init(final PrintStream out) throws Exception {
        String sID = getID();
        if (sID == null) {
            sID = "";
        }
        out.print(sID + ".mean\t" + sID + ".variance\t" + sID + ".coefficientOfVariation\t");
    }


    @Override
    public void log(final int nSample, final PrintStream out) {

        out.print(multinomialTestStatistic + "\t");
    }


    @Override
    public void close(final PrintStream out) {
        // nothing to do
    }


}
