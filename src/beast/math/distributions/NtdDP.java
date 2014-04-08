package beast.math.distributions;

import beast.core.Input;
import beast.core.Description;
import beast.core.parameter.ParameterList;

/**
 * @author Chieh-Hsi Wu
 */
@Description("Dirichlet process for model parameters of nucletoide subsitution model.")
public class NtdDP extends DirichletProcess{
    public Input<ParametricDistribution> paramBaseDistrInput = new Input<ParametricDistribution>(
            "paramBaseDistr",
            "The base distribution of the dirichlet process",
            Input.Validate.REQUIRED
    );
    public Input<ParametricDistribution> modelBaseDistrInput = new Input<ParametricDistribution>(
            "modelBaseDistr",
            "The base distribution of the dirichlet process",
            Input.Validate.REQUIRED
    );
    public Input<ParametricDistribution> freqBaseDistrInput = new Input<ParametricDistribution>(
            "freqBaseDistr",
            "The base distribution of the dirichlet process",
            Input.Validate.REQUIRED
    );


    private ParametricDistribution paramBaseDistr;
    private ParametricDistribution modelBaseDistr;
    private ParametricDistribution freqBaseDistr;

    public NtdDP(){
        baseDistrInput.setRule(Input.Validate.OPTIONAL);
    }
    public void initAndValidate(){
        paramBaseDistr = paramBaseDistrInput.get();
        modelBaseDistr = modelBaseDistrInput.get();
        freqBaseDistr = freqBaseDistrInput.get();

        alpha = alphaInput.get();
        dpValuable = dpValuableInput.get();
        initialise();

    }

    public double calcLogP(
            ParameterList paramList,
            ParameterList modelList,
            ParameterList freqList) throws Exception {
        if(requiresRecalculation()){
            refresh();
            //System.out.println("Refreshing!");
        }
        int dimParam = paramList.getDimension();
        double logP = 0.0;
        for(int i = 0; i < dimParam; i ++){
            logP += paramBaseDistr.calcLogP(paramList.getParameter(i));
        }
        //System.out.println("flag1: "+logP);
        int dimModel = modelList.getDimension();
        for(int i = 0; i < dimModel; i ++){
            logP += modelBaseDistr.calcLogP(modelList.getParameter(i));
        }
        //System.out.println("flag2: "+logP);
        int dimFreq = freqList.getDimension();
        for(int i = 0; i < dimFreq; i ++){
            logP += freqBaseDistr.calcLogP(freqList.getParameter(i));
        }
        //System.out.println("flag3: "+logP);

        int[] counts = dpValuable.getClusterCounts();

        logP+=counts.length*Math.log(alpha.getValue());
        //System.out.println("flag4: "+logP+" counts.length: "+counts.length+" alpha.getValue():"+alpha.getValue());

        for(int count: counts){
            logP+=gammas[count];
        }
        //System.out.println("flag5: "+logP);


        logP-=denominator;



        //System.out.println("flag6: "+logP);
        //if(Double.isNaN(logP))
        //    throw new RuntimeException("wrong!!!");
        return logP;


    }

    public boolean requiresRecalculation(){
        return alpha.somethingIsDirty()
                || paramBaseDistr.isDirtyCalculation()
                || modelBaseDistr.isDirtyCalculation()
                || freqBaseDistr.isDirtyCalculation();
    }
    public ParametricDistribution getBaseDistribution(){
        return null;
    }
    public ParametricDistribution getParamBaseDistr(){
        return paramBaseDistr;
    }

    public ParametricDistribution getModelBaseDistr(){
        return modelBaseDistr;
    }

    public ParametricDistribution getFreqBaseDistr(){
        return freqBaseDistr;
    }

}
