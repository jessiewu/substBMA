package beast.evolution.operators;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.ParameterList;
import beast.core.parameter.RealParameter;
import beast.math.distributions.MultivariateNormal;
import beast.math.matrixAlgebra1.Matrix;

/**
 * @author Chieh-Hsi Wu
 */
@Description("This is a sampler for the mean vector of a multivariate normal distribution when a Multivariate normal prior is applied.")
public class MeanMVNSampler extends Operator{
    public Input<RealParameter> priorMeanInput = new Input<RealParameter>(
            "priorMean",
            "The mean vector of the multivariate normal prior on the mean vector",
            Input.Validate.REQUIRED
    );

    public Input<RealParameter> priorPrecisionInput = new Input<RealParameter>(
            "priorPrecision",
            "The precision matrix (as one-dimensional vector) of the multivariate normal prior on the mean vector",
            Input.Validate.REQUIRED
    );

    public Input<RealParameter> meanInput = new Input<RealParameter>(
            "mean",
            "The mean vector of the multivariate normal to be sampled",
            Input.Validate.REQUIRED
    );

    public Input<RealParameter> precisionInput = new Input<RealParameter>(
            "precision",
            "The precision matrix of the multivariate normal to be sampled",
            Input.Validate.REQUIRED
    );

    public Input<ParameterList> xInput = new Input<ParameterList>(
            "x",
            "The observations assumed to have a multivariate normal distribution.",
            Input.Validate.REQUIRED
    );

    private double[] priorMean;
    private double[][] priorPrecision;
    private ParameterList x;
    private RealParameter mean;
    private int dim;

    public void initAndValidate(){

        RealParameter priorMeanParam = priorMeanInput.get();
        RealParameter priorPrecisionParam = priorPrecisionInput.get();
        dim = priorMeanParam.getDimension();
        if(dim != (int)Math.sqrt(priorPrecisionParam.getDimension())){
            throw new RuntimeException("The rows and column counts of the precision matrix" +
                    " should equal to the dimension of the mean vector.");
        }
        priorMean = new double[dim];
        for(int i = 0;i < dim;i++){
            priorMean[i] = priorMeanParam.getValue(i);

        }



        priorPrecision = new double[dim][dim];

        int k = 0;

        for(int i = 0; i< dim; i++){
            for(int j = 0; j < dim; j++){
                priorPrecision[i][j] = priorPrecisionParam.getValue(k++);
            }
        }
        x = xInput.get();



    }

    public double proposal(){

        //System.out.println("flag1");
        int obsCount = x.getDimension();

        double[] obsMean = new double[dim];
        for(int i = 0; i < dim; i++){
            for(int j = 0; j < obsCount; j++){
                //System.out.println(i+" "+j);
                obsMean[i] += x.getValue(j,i);
            }
            obsMean[i] = obsMean[i]/obsCount;
        }

        //System.out.println("flag2");
        double[][] precision = new double[dim][dim];
        RealParameter precParam = precisionInput.get();
        int k = 0;
        for(int i = 0; i < dim; i++){
            for(int j = 0; j < dim; j++){
                precision[i][j] = precParam.getValue(k++);
            }
        }

        double[][] newPrecision = new double[dim][dim];
        for(int i = 0; i < dim; i++){
            for(int j = 0; j < dim; j++){
                newPrecision[i][j] = priorPrecision[i][j]+obsCount*precision[i][j];
                //System.out.print(newPrecision[i][j]+" ");
            }
            //System.out.println();
        }

        //System.out.println("flag3");
        double[] tmp = new double[dim];
        for(int i = 0; i < dim;i++){
            for(k = 0; k < dim; k++){
                tmp[i] += obsCount*precision[i][k]*obsMean[k]+priorPrecision[i][k]*priorMean[k];
                //System.out.println(precision[i][k]+" "+obsMean[k]);
            }
            //System.out.print(newMean[i]+" ");

        }

        //System.out.println();
        //System.out.println("flag4");
        double[] newMean = new double[dim];
        double[][] newCovariance = (new Matrix(newPrecision)).inverse().toComponents();
        for(int i = 0; i < dim; i++){
            for(k = 0; k < dim; k++){
                newMean[i] += newCovariance[i][k]*tmp[k];
            }
            //System.out.print(newMean2[i]+" ");
        }
        //System.out.println();

        double[] draw = MultivariateNormal.nextMultivariateNormalVariance(newMean,newCovariance);

        mean = meanInput.get(this);
        for(int i = 0; i < draw.length; i++){
            mean.setValue(i,draw[i]);
        }

        //System.out.println("flag5");
        return Double.POSITIVE_INFINITY;
    }

    public static void main(String[] args){
        try{
            RealParameter priorMean = new RealParameter(new Double[]{1.0, 2.0, 3.0});
            RealParameter priorPrec = new RealParameter(new Double[]{
                    5.022036183614738469316, -0.044273610915457527193, -0.198225030689662129468,
                    -0.044273610915457527193,  1.002193556177174871280, -0.058360668934012191467,
                    -0.198225030689662129468, -0.058360668934012184528,  2.011430641363627369600
            });
            RealParameter mean = new RealParameter(new Double[]{0.0,0.0,0.0});
            RealParameter prec = new RealParameter(new Double[]{
                    5.11210762331838619588, -0.26905829596412556004, -1.25560538116591935420,
                    -0.26905829596412561555, 2.64573991031390098883, -0.98654708520179368314,
                    -1.25560538116591935420, -0.98654708520179357212, 2.06278026905829614535
            });

            RealParameter x1 = new RealParameter(new Double[]{1.0,3.0,1.0});
            RealParameter x2 = new RealParameter(new Double[]{2.0,1.0,5.0});
            RealParameter x3 = new RealParameter(new Double[]{2.0,4.0,2.0});
            ParameterList x = new ParameterList();
            x.initByName(
                    "parameter", x1,
                    "parameter", x2,
                    "parameter", x3
            );

            MeanMVNSampler sampler = new MeanMVNSampler();
            sampler.initByName(
                    "priorMean", priorMean,
                    "priorPrecision", priorPrec,
                    "mean", mean,
                    "precision", prec,
                    "x", x
            );

            sampler.proposal();


            /*1.493206204232796 2.570970475304471 2.650196067131461

            by hand: 1.493206204232796 2.570970475304471 2.650196067131461 */

        }catch(Exception e){
            throw new RuntimeException(e);

        }
    }



}
