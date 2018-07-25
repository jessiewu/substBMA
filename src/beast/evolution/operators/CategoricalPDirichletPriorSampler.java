package beast.evolution.operators;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.ParameterList;
import beast.core.parameter.RealParameter;
import beast.math.distributions.DirichletDistribution;

/**
 * @author Chieh-Hsi Wu
 */
@Description("This is the Gibb sampler for categorial distribution assuming a Dirichlet prior.")
public class CategoricalPDirichletPriorSampler extends Operator {
    public Input<RealParameter> categoricalProbsInput = new Input<RealParameter>(
            "categorialProbs",
            "The probabailities for the categories.",
            Input.Validate.REQUIRED
    );
    public Input<RealParameter> dirichletPriorCountsInput = new Input<RealParameter>(
            "dirichletPriorCounts",
            "The counts parameter of the Dirichlet distribution",
            Input.Validate.REQUIRED
    );
    public Input<ParameterList> xInput = new Input<ParameterList>(
            "x",
            "The list of values that have an common categorical distribution.",
            Input.Validate.REQUIRED
    );

    public Input<Double> offsetInput = new Input<Double>(
            "offset",
            "The offset of the categorical distribution",
            0.0
    );

    private RealParameter dirichletPriorCounts;
    private ParameterList x;
    private double offset;

    public void initAndValidate(){
        dirichletPriorCounts = dirichletPriorCountsInput.get();
        x = xInput.get();
        offset = offsetInput.get();

    }

    public double proposal(){
        Double[] newCounts = new Double[dirichletPriorCounts.getDimension()];
        for(int i = 0; i < newCounts.length; i++){
            newCounts[i] = dirichletPriorCounts.getValue();
        }
        int sampleSize = x.getDimension();

        for(int i = 0; i < sampleSize; i++){
            newCounts[(int)(double)(x.getParameter(i).getValue()-offset)]++;
        }
        try{
            DirichletDistribution dirichlet = new DirichletDistribution(newCounts);
            // offset has been handled above, so not use sampler in QuietRealParameter#getSample
            Double[][] draw = dirichlet.sample(1);
            RealParameter probs = categoricalProbsInput.get(this);
            for(int i = 0 ; i < draw[0].length; i++){
                probs.setValue(i,draw[0][i]);
            }
        }catch(Exception e){
            throw new RuntimeException(e);

        }

        return Double.POSITIVE_INFINITY;



    }


}
