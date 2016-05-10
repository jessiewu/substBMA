package beast.math.distributions;

import beast.core.Description;
import beast.core.Function;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.util.Randomizer;

/**
 * @author Chieh-Hsi Wu
 */
@Description("Extends Dirichlet class to add the sample method.")
public class DirichletDistribution extends Dirichlet{

    public Input<RealParameter> scaleInput = new Input<RealParameter>(
            "scale",
            "A parameter that scales the alpha parameter",
            new RealParameter(new Double[]{1.0})
    );

    RealParameter alpha;
    RealParameter scale;



    public DirichletDistribution() throws Exception{
        alpha = new RealParameter(new Double[]{Double.NaN,Double.NaN,Double.NaN,Double.NaN});
        scale = scaleInput.get();

    }

    public DirichletDistribution(Double[] alpha) throws Exception{
        this.alpha = new RealParameter(alpha);
        scale = scaleInput.get();
    }

    /*public void setAlpha(Double[] alphaValues){
        for(int i = 0; i < alphaValues.length;i++){
            ((QuietRealParameter)alpha).setValueQuietly(i,alphaValues[i]);
        }
    } */



    public void initAndValidate() {
        //System.out.println("This is initiated!");
        //System.out.println(m_alpha.get().getID()+" "+m_alpha.get().getValue(0));

        alpha = alphaInput.get();
        //System.out.println(alpha.getID()+" " +alpha.getValue(0));
        scale = scaleInput.get();
	}

    public double getScaleValue(){
        return scale.getValue();
    }

    @Override
    public Double[][] sample(int size){
        double scaleVal = scale.getValue();
        int dim = alpha.getDimension();
        Double[][] samples = new Double[size][];
        try{

            for(int i =0; i < size;i++){
                Double[] dirichletSample = new Double[dim];
                double sum = 0.0;
                for(int j  =0; j < dim;j++){
                    dirichletSample[j] = Randomizer.nextGamma(alpha.getValue(j)*scaleVal,1.0);
                    sum += dirichletSample[j];
                }
                for(int j = 0; j < dim;j++){
                    dirichletSample[j] = dirichletSample[j]/sum;
                }
                samples[i] = dirichletSample;

            }
        }catch(Exception e){
            throw new RuntimeException(e);

        }
        return samples;

    }



    public static Double[] nextDirichletScale(Double[] alpha, double scaleVal){

        Double[] sample = new Double[alpha.length];
        double sum = 0.0;
        for(int j  =0; j < sample.length;j++){
            sample[j] = Randomizer.nextGamma(alpha[j]*scaleVal,1.0);
            sum += sample[j];
        }

        for(int j = 0; j < sample.length;j++){
            sample[j] = sample[j]/sum;
        }
        return sample;

    }


    public static double logPDF(Double[] proposal, Double[] alpha, double scaleVal){
        //Put alpha values in an array.
        double[] fAlpha = new double[alpha.length];
        for(int i = 0; i < fAlpha.length; i++){
            fAlpha[i] = alpha[i]*scaleVal;
        }



        double fLogP = 0;
        double fSumAlpha = 0;
        for (int i = 0; i < proposal.length; i++) {
            double fX = proposal[i];
            fLogP += (fAlpha[i]-1) * Math.log(fX);
            fLogP -= org.apache.commons.math.special.Gamma.logGamma(fAlpha[i]);
            fSumAlpha += fAlpha[i];
        }
        fLogP += org.apache.commons.math.special.Gamma.logGamma(fSumAlpha);
        return fLogP;


    }

    public void printDetails(){
        System.err.println("alpha:"+alpha);
        System.err.println("scale:"+scale);
    }

    public static void  main(String[] args){
        try{
            RealParameter alpha = new RealParameter(new Double[]{3.0,2.0,1.0,3.0});
            DirichletDistribution dirichlet = new DirichletDistribution();
            dirichlet.initByName(
                    "alpha", alpha
            );

             Double[][] samples = dirichlet.sample(1000);

            for(int i = 0; i < 1000; i++){
                RealParameter sample = new RealParameter(samples[i]);
                //sample.log(0,System.out);
                //System.out.println();

            }

            //fitted alpha from VGAM
            //3.9830039 1.9792060 0.9378492 2.9612180

            samples = new Double[1000][];
            for(int i = 0 ; i < 1000; i++){
                //System.out.println((i+1)+": ");
                RealParameter sample = new RealParameter(nextDirichletScale(new Double[]{1.5,1.0,0.5,1.5},0.02));
                sample.log(0,System.out);
                System.out.println();
            }



            RealParameter alpha2 = new RealParameter(new Double[]{0.2548157632156566,0.1485647835002092,0.38521176312917504,0.016041715566902426});
            alpha2.setID("XXX");
            DirichletDistribution dirichlet2 = new DirichletDistribution();
            dirichlet2.initByName(
                    "alpha", alpha2
            );
            RealParameter x = new RealParameter(new Double[]{1.0,	4.178667426785708E-27,	4.0458954268572387E-76,	1.567575228801415E-55	});
            System.out.println(dirichlet2.calcLogP(x));
        }catch(Exception e){
            throw new RuntimeException(e);

        }

    }


    @Override
    public double calcLogP(Function pX) {

        double scaleVal = scale.getValue();
        Double [] fAlpha = alpha.getValues();
        //System.out.println(alpha.getID());
        for(int i = 0; i < pX.getDimension(); i++){
            if(pX.getArrayValue(i) == 0.0 || pX.getArrayValue(i) == 1.0){
                return Double.NEGATIVE_INFINITY;
            }
        }


        if (alpha.getDimension() != pX.getDimension()) {
            throw new RuntimeException("Dimensions of alpha and x should be the same, but dim(alpha)=" + alpha.getDimension()
                    + " and dim(x)=" + pX.getDimension());
        }
        for(int i = 0; i < fAlpha.length; i++){
            fAlpha[i] = fAlpha[i]*scaleVal;
            //System.out.println(alpha.getValue(i)+" "+scaleVal);
        }
        double fLogP = 0;
        double fSumAlpha = 0;
        for (int i = 0; i < pX.getDimension(); i++) {
            double fX = pX.getArrayValue(i);
            fLogP += (fAlpha[i]-1) * Math.log(fX);
            //System.out.println((fAlpha[i]-1) +" "+ Math.log(fX));
            fLogP -= org.apache.commons.math.special.Gamma.logGamma(fAlpha[i]);
            fSumAlpha += fAlpha[i];
        }
        fLogP += org.apache.commons.math.special.Gamma.logGamma(fSumAlpha);
        //System.out.println(fLogP);
        return fLogP;
    }

}
