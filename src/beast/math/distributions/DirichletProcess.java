package beast.math.distributions;

import org.apache.commons.math.distribution.ContinuousDistribution;
import beast.core.parameter.*;
import beast.core.Input;
import beast.core.Description;
import beast.core.Valuable;

/**
 * @author Chieh-Hsi Wu
 */
@Description("This class that the dirichlet process.")
public class DirichletProcess extends ParametricDistribution{

    public Input<DPValuable> dpValuableInput = new Input<DPValuable>(
            "dpVal",
            "reports the counts in each cluster",
            Input.Validate.REQUIRED
    );

    public Input<ParametricDistribution> baseDistrInput = new Input<ParametricDistribution>(
            "baseDistr",
            "The base distribution of the dirichlet process",
            Input.Validate.REQUIRED
    );
    

    public Input<RealParameter> alphaInput = new Input<RealParameter>(
            "alpha",
            "The concentration parameter of the Dirichlet process prior",
            Input.Validate.REQUIRED
    );

    private ParametricDistribution baseDistribution;
    RealParameter alpha;
    DPValuable dpValuable;
    double denominator;
    double storedDenominator;
    double[] gammas;


    public void initAndValidate(){

        alpha = alphaInput.get();
        dpValuable = dpValuableInput.get();
        baseDistribution = baseDistrInput.get();

        //Yes I know that we are only counting number of memebers is the existing clusters
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
            denominator += Math.log(alpha.getValue()+i);
        }

    }
    public double calcLogP(Valuable xList) throws Exception {
        if(requiresRecalculation()){
            refresh();
        }
        int dimParam = xList.getDimension();
        double logP = 0.0;
        for(int i = 0; i < dimParam; i ++){
            //System.out.println(getID()+" "+((ParameterList)xList).getParameter(i));
            logP+= baseDistribution.calcLogP(((ParameterList)xList).getParameter(i));
        }
        //System.out.println("flag1: "+logP);

        int[] counts = dpValuable.getClusterCounts();

        logP+=counts.length*Math.log(alpha.getValue());
        //System.err.println("flag2: "+logP+" counts: "+counts.length+" a: "+alpha.getValue());

        for(int count: counts){
            //System.err.println(count+" "+gammas[count]);
            logP+=gammas[count];
        }
        //System.err.println("flag3: "+logP);
        logP-=denominator;

        //if(Double.isNaN(logP))
        //    throw new RuntimeException("wrong!!!");
        return logP;


    }

    public void store(){
        storedDenominator = denominator;
        super.store();
    }

    public void restore(){
        denominator = storedDenominator;
        super.restore();
    }

    public boolean requiresRecalculation(){
        return alpha.somethingIsDirty() || baseDistribution.isDirtyCalculation();
    }

    public ContinuousDistribution getDistribution(){
        throw new RuntimeException("Dirichlet process is a discrete distribution.");
    }

    public ParametricDistribution getBaseDistribution(){
        return baseDistribution;
    }

    public double getConcParameter(){
        return alpha.getValue();
    }


}
