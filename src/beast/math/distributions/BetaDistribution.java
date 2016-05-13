package beast.math.distributions;

import beast.core.Function;
import beast.core.parameter.RealParameter;
import beast.util.Randomizer;

/**
 * Created by IntelliJ IDEA.
 * User: cwu080
 * Date: 2/07/13
 * Time: 12:05 PM
 * To change this template use File | Settings | File Templates.
 */
public class BetaDistribution extends DirichletDistribution {

    public BetaDistribution() throws Exception{
        alpha = null;
        scale = scaleInput.get();

    }

    public BetaDistribution(Double[] alpha) throws Exception{
        if(alpha.length != 2)
            throw new ParameterLengthException();
        this.alpha = new RealParameter(alpha);
        scale = scaleInput.get();

    }


    public void initAndValidate() {
        alpha = alphaInput.get();
        if(alpha.getDimension() != 2){
            throw new ParameterLengthException();
        }
        scale = scaleInput.get();
	}

    @Override
    public Double[][] sample(int size){
        double scaleVal = scale.getValue();
        Double[][] samples = new Double[size][1];
        try{

            for(int i =0; i < size;i++){
                Double[] dirichletSample = new Double[2];
                double sum = 0.0;
                for(int j  =0; j < 2;j++){
                    dirichletSample[j] = Randomizer.nextGamma(alpha.getValue(j) * scaleVal, 1.0);
                    sum += dirichletSample[j];
                }

                dirichletSample[0] = dirichletSample[0]/sum;
                samples[i][0] = dirichletSample[0];

            }
        }catch(Exception e){
            throw new RuntimeException(e);

        }
        return samples;

    }

     public static double nextBetaScale(double alpha, double beta, double scaleVal){

        Double[] sample = new Double[2];

        if(alpha == 0.0){
            return 0.0;
        }
        sample[0] = Randomizer.nextGamma(alpha * scaleVal, 1.0);
        if(beta == 0.0){
            sample[1] = 0.0;
        } else{
            sample[1] = Randomizer.nextGamma(beta * scaleVal, 1.0);
        }



        return sample[0]/(sample[0] + sample[1]);

    }


    public static double nextBetaScale(Double[] alpha, double scaleVal){

        Double[] sample = new Double[alpha.length];
        double sum = 0.0;
        for(int j  = 0; j < sample.length;j++){
            sample[j] = Randomizer.nextGamma(alpha[j] * scaleVal, 1.0);
            sum += sample[j];
        }

        return sample[0]/sum;

    }

    public static double logPDF(double proposal, double alpha, double scaleVal){
        //Put alpha values in an array.
        double[] fAlpha = new double[2];
        fAlpha[0] = alpha*scaleVal;
        fAlpha[1] = (1.0 - alpha)*scaleVal;



        double fLogP = 0;
        double fSumAlpha = 0;
        double[] proposalVec = new double[]{proposal, 1.0 - proposal};

        for (int i = 0; i < proposalVec.length; i++) {
            double fX = proposalVec[i];
            fLogP += (fAlpha[i]-1) * Math.log(fX);
            //System.out.println(fLogP+" "+fX+" "+Math.log(fX));
            fLogP -= org.apache.commons.math.special.Gamma.logGamma(fAlpha[i]);
            //System.out.println(fLogP);
            fSumAlpha += fAlpha[i];
        }

        fLogP += org.apache.commons.math.special.Gamma.logGamma(fSumAlpha);
        return fLogP;


    }


    class ParameterLengthException extends RuntimeException{
        public String getMessage(){
            return  "The length of the alpha vector of a Beta distribution must be two.";

        }

    }

    @Override
    public double calcLogP(Function pX) {

        double scaleVal = scale.getValue();
        Double [] fAlpha = alpha.getValues();

        if(pX.getArrayValue() == 1.0 || pX.getArrayValue()== 0.0){
            return Double.NEGATIVE_INFINITY;
        }

        if(pX.getArrayValue() > 1.0 || pX.getArrayValue() < 0.0){
            throw new RuntimeException("Value must be between 0.0 and 1.0: "+pX.getArrayValue());
        }


        /*if (alpha.getDimension() != pX.getDimension()) {
            throw new Exception("Dimensions of alpha and x should be the same, but dim(alpha)=" + alpha.getDimension()
                    + " and dim(x)=" + pX.getDimension());
        }*/
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


    public static void  main(String[] args){
        try{
            RealParameter alpha = new RealParameter(new Double[]{3.0,2.0});
            BetaDistribution beta = new BetaDistribution();
            beta.initByName(
                    "alpha", alpha
            );

             Double[][] samples = beta.sample(1000);

            for(int i = 0; i < 1000; i++){
                RealParameter sample = new RealParameter(samples[i]);
                //sample.log(0,System.out);
                //System.out.println();

            }

            //fitted alpha from VGAM
            //2.893453 1.898905

            samples = new Double[1000][];
            for(int i = 0 ; i < 1000; i++){
                //System.out.println((i+1)+": ");
                RealParameter sample = new RealParameter(new Double[]{nextBetaScale(new Double[]{15.0,10.0},0.2)});
                //sample.log(0,System.out);
                //System.out.println();
            }
            //Estimated counts from VGAM:2.864905 1.866268




            RealParameter alpha2 = new RealParameter(new Double[]{2.5e-10,7.5e-10});
            alpha2.setID("XXX");
            DirichletDistribution dirichlet2 = new DirichletDistribution();
            dirichlet2.initByName(
                    "alpha", alpha2
            );
            RealParameter x = new RealParameter(new Double[]{0.1,0.9});
            System.out.println(dirichlet2.calcLogP(x));
            //R computation: 0.8440593815828779
        }catch(Exception e){
            throw new RuntimeException(e);

        }

    }
}
