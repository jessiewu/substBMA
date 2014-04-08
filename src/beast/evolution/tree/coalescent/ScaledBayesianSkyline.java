package beast.evolution.tree.coalescent;

import beast.core.Input;
import beast.evolution.tree.Scaler;
import beast.math.Binomial;

/**
 * Created by IntelliJ IDEA.
 * User: Jessie Wu
 * Date: 2/12/13
 * Time: 12:23 PM
 * To change this template use File | Settings | File Templates.
 */
public class ScaledBayesianSkyline extends BayesianSkyline {
    public Input<Scaler> scalerInput = new Input<Scaler>(
            "scaler",
            "A component that scales the tree intervals",
            Input.Validate.REQUIRED
    );

    private Scaler scaler;

    public void initAndValidate() throws Exception{
        scaler = scalerInput.get();
        super.initAndValidate();
    }

    public boolean requiresRecalculation(){
        if(scaler.isDirtyCalculation()){
            return true;
        }
        return super.requiresRecalculation();
    }

    @Override
    public double calculateLogP() throws Exception {
        if (!m_bIsPrepared) {
            prepare();
        }

        logP = 0.0;

        double scaleFactor = scaler.getScaleFactor();

        double currentTime = 0.0;

        int groupIndex = 0;
        // int[] groupSizes = getGroupSizes();
        // double[] groupEnds = getGroupHeights();

        int subIndex = 0;

        //ConstantPopulation cp = new ConstantPopulation();// Units.Type.YEARS);

        for (int j = 0; j < intervals.getIntervalCount(); j++) {

            // set the population size to the size of the middle of the current
            // interval
            final double ps = getPopSize(currentTime + (intervals.getInterval(j) / 2.0));
            //cp.setN0(ps);
            if (intervals.getIntervalType(j) == IntervalType.COALESCENT) {
                subIndex += 1;
                if (subIndex >= groupSizes.getValue(groupIndex)) {
                    groupIndex += 1;
                    subIndex = 0;
                }
            }

            logP += calculateIntervalLikelihood(ps, intervals.getInterval(j), currentTime,
                    intervals.getLineageCount(j), intervals.getIntervalType(j),
                    scaleFactor);

            // insert zero-length coalescent intervals
            int diff = intervals.getCoalescentEvents(j) - 1;
            for (int k = 0; k < diff; k++) {
                //cp.setN0(getPopSize(currentTime));
                double fPopSize = getPopSize(currentTime);
                logP += calculateIntervalLikelihood(fPopSize, 0.0,
                        currentTime, intervals.getLineageCount(j) - k - 1,
                        IntervalType.COALESCENT,
                        scaleFactor);
                subIndex += 1;
                if (subIndex >= groupSizes.getValue(groupIndex)) {
                    groupIndex += 1;
                    subIndex = 0;
                }
            }

            currentTime += intervals.getInterval(j);
        }
        return logP;
    }

    public static double calculateIntervalLikelihood(
            double popSize, double width,
            double timeOfPrevCoal, int lineageCount, IntervalType type,
            double scaleFactor) {
        final double timeOfThisCoal = width + timeOfPrevCoal;
        //System.out.println( width +" " + timeOfPrevCoal);
        //final double intervalArea = (timeOfThisCoal - timeOfPrevCoal) / popSize;
        final double intervalArea = width*scaleFactor / popSize;
        //demogFunction.getIntegral(timeOfPrevCoal, timeOfThisCoal);
        //System.out.println(intervalArea+" "+popSize);
        final double kchoose2 = Binomial.choose2(lineageCount);
        double like = -kchoose2 * intervalArea;

        switch (type) {
            case COALESCENT:
                final double demographic = Math.log(popSize);//demogFunction.getLogDemographic(timeOfThisCoal);
                like += -demographic;

                break;
            default:
                break;
        }

        return like;
    }

}
