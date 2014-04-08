package beast.math.distributions;

import beast.core.Distribution;
import beast.core.Input;
import beast.core.parameter.DPValuable;
import beast.core.parameter.RealParameter;
import org.apache.commons.math.distribution.ContinuousDistribution;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Chieh-Hsi Wu
 */
public class FlexibleCompoundDirichletProcess extends ParametricDistribution {
    public Input<DPValuable> dpValuableInput = new Input<DPValuable>(
            "dpVal",
            "reports the counts in each cluster",
            Input.Validate.REQUIRED
    );

    public Input<List<Distribution>> baseDistrInput = new Input<List<Distribution>>(
            "baseDistr",
            "The base distribution of the dirichlet process",
            new ArrayList<Distribution>(),
            Input.Validate.REQUIRED
    );


    public Input<RealParameter> concInput = new Input<RealParameter>(
            "concentrationParameter",
            "The concentration parameter of the Dirichlet process prior",
            Input.Validate.REQUIRED
    );

    private List<Distribution> baseDistributions;
    //double[] alphaPowers;
    RealParameter concentrationParameter;
    DPValuable dpValuable;
    double denominator;
    double storedDenominator;
    double[] gammas;


    public void initAndValidate(){

        concentrationParameter = concInput.get();
        dpValuable = dpValuableInput.get();
        baseDistributions = baseDistrInput.get();

        //Yes I know that we are only counting number of members in the existing clusters
        //So there won't be any clusters of size 0.
        //But for the sake of convenience for later computation I'm going to start from 0.


        initialise();

    }

    public void initialise(){
        //Yes I know that we are only counting number of memebers is the existing clusters
        //So there won't be any clusters of size 0.
        //But for the sake of convenience for later computation I'm going to start from 0.
        refresh();
        gammas = new double[dpValuable.getPointerDimension()+1];
        gammas[0] = Double.NaN;
        gammas[1] = 0;
        for(int i = 2; i < gammas.length; i++){
            gammas[i] = gammas[i-1]+Math.log(i-1);
        }
    }

    public void refresh(){
        //Yes I know that we are only counting number of memebers is the existing clusters
        //So there won't be any clusters of size 0.
        //But for the sake of convenience for later computation I'm going to start from 0.

        denominator = 0;
        int n = dpValuable.getPointerDimension();

        for(int i = 0; i < n;i++){
            denominator += Math.log(concentrationParameter.getValue()+i);
        }

    }
    public double calculateLogP() throws Exception {
        if(requiresRecalculation()){
            refresh();
        }
        int listCount = baseDistributions.size();

        double logP = 0.0;

        for(int i = 0; i < listCount;i++){
            logP += baseDistributions.get(i).calculateLogP();

        }


        int[] counts = dpValuable.getClusterCounts();

        logP+=counts.length*Math.log(concentrationParameter.getValue());
        //System.err.println("flag2: "+logP);

        for(int count: counts){
            logP+=gammas[count];
        }
        //System.err.println("flag3: "+logP);
        logP-=denominator;
        //System.err.println("flag4: "+logP);
        //if(Double.isNaN(logP))
        //    throw new RuntimeException("wrong!!!");
        /*if(logP == Double.NEGATIVE_INFINITY){
             for(int i = 0; i < listCount;i++){
            int dimParam = xLists.get(i).getDimension();
            for(int j = 0; j < dimParam; j ++){
                System.out.print(baseDistributions.get(i).calcLogP((xLists.get(i)).getParameter(j)) +" ");
            }
                 System.out.println();
        }

        }*/
        return logP;


    }
    public boolean requiresRecalculation(){
        if(concentrationParameter.somethingIsDirty()){
            return true;
        }
        for(Distribution paramDistr: baseDistributions){
            if(paramDistr.isDirtyCalculation()){
                return true;
            }

        }
        return false;
    }

    public ContinuousDistribution getDistribution(){
        throw new RuntimeException("Dirichlet process is a discrete distribution.");
    }

    public List<Distribution> getBaseDistributions(){
        return baseDistributions;
    }

    public void store(){
        storedDenominator = denominator;
        super.store();
    }

    public void restore(){
        denominator = storedDenominator;
        super.restore();
    }

    public double getConcParameter(){
        return concentrationParameter.getValue();
    }
}
