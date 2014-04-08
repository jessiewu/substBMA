package beast.evolution.substitutionmodel;

import beast.core.Description;
import beast.core.parameter.QuietRealParameter;
import beast.evolution.tree.Node;
import beast.core.parameter.RealParameter;

import java.util.Arrays;

/**
 * @author Chieh-Hsi Wu
 */
@Description("This class switches between analytic and numerical computation depends on the substitution model selected.")
public class SwitchingNtdBMA extends NtdBMA{
    private double freqA;
    private double freqC;
    private double freqG;
    private double freqT;
    protected boolean updateTN = true;
    //private boolean storedUpdateTN;
    protected boolean updateHKY = true;
    //protected boolean storedUpdateHKY;

    public SwitchingNtdBMA(){
        super();
    }

    public SwitchingNtdBMA(
            QuietRealParameter logKappa,
            QuietRealParameter logTN,
            QuietRealParameter logAC,
            QuietRealParameter logAT,
            QuietRealParameter logGC,
            //RealParameter logGT,
            QuietRealParameter modelChoose,
            QuietRealParameter frequencies){
        super(logKappa,logTN, logAC, logAT,
                logGC, modelChoose, frequencies);
    }

    public void initialize(
            QuietRealParameter logKappa,
            QuietRealParameter logTN,
            QuietRealParameter logAC,
            QuietRealParameter logAT,
            QuietRealParameter logGC,
            //RealParameter logGT,
            QuietRealParameter modelChoose,
            QuietRealParameter frequencies){
        updateTN = true;
        updateHKY = true;
        super.initialize(logKappa,logTN, logAC, logAT,
                logGC, modelChoose, frequencies);

    }


    public void getTransitionProbabilities(Node node, double fStartTime, double fEndTime, double fRate, double[] matrix) {
        //System.out.println(fRate+" "+fStartTime+" "+fEndTime );
        if(fRate < 0.0){
            for(int j = 0; j < matrix.length; j++){
                matrix[j] = 0.0;
            }
            return;
        }

        if(modelChoose.getValue() == GTR){
            getGTRTransitionProbabilities(node, fStartTime, fEndTime, fRate, matrix);
        }else if(modelChoose.getValue() == TN93){
            getTN93TransitionProbabilities(node, fStartTime, fEndTime, fRate, matrix);

        }else if(modelChoose.getValue() == JC69){
            getJC69TransitionProbabilities(node, fStartTime, fEndTime, fRate, matrix);

        }else{
            getHKY85TransitionProbabilities(node, fStartTime, fEndTime, fRate, matrix);

        }


        double sum= 0.0;

        for(int i = 0; i < NTD_PAIRS;i++){
            sum+= matrix[i];
        }
        if(Double.isNaN(sum)){
            for(int j = 0; j < matrix.length; j++){
                matrix[j] = 0.0;
            }

        }else if((sum-4.0) > 0.005){
            for(int i = 0; i < matrix.length; i++){
                matrix[i] = 0.0;
            }
        }else{
            for(int i = 0; i < matrix.length; i++){
                if(matrix[i] < 0.0){
                    matrix[i] = 0.0;
                }else if(matrix[i] > 1.0){
                    matrix[i] = 1.0;
                }
            }

        }

        /*System.out.println("start: "+fStartTime+" end: "+fEndTime+" fRate: "+fRate);
            printDetails();

        for(int i = 0; i < matrix.length;i++){
            System.out.print(matrix[i]+" ");
        }


        System.out.println();*/
        /*if(matrix[0]>1.5){
            System.out.println("start: "+fStartTime+" end: "+fEndTime+" fRate: "+fRate);
            printDetails();



        }



        System.out.println("Id #: "+getIDNumber()+" ID: "+getID());
        System.out.println("modelChoose "+modelChoose.getValue());
        System.out.println("start: "+fStartTime+" end: "+fEndTime+" fRate: "+fRate);


        */
    }


    public void getJC69TransitionProbabilities(
            Node node,
            double fStartTime,
            double fEndTime,
            double fRate,
            double[] matrix) {
        //System.err.println("matrix: "+matrix.length);
        double distance = (fStartTime - fEndTime) * fRate;

        int k = 0;
        for(int i = 0; i < STATE_COUNT;i++){
            for(int j = 0; j < STATE_COUNT;j++){
                if(i==j){
                    matrix[k++] = (1.0+3.0*Math.exp(-4.0/3.0*distance))/4.0;
                }else{
                    matrix[k++] = (1.0-Math.exp(-4.0/3.0*distance))/4.0;
                }

            }
        }
    }



    private double beta, A_R, A_Y;
    private double tab1A, tab2A, tab3A;
    private double tab1C, tab2C, tab3C;
    private double tab1G, tab2G, tab3G;
    private double tab1T, tab2T, tab3T;

    public void getHKY85TransitionProbabilities(Node node, double fStartTime, double fEndTime, double fRate, double[] matrix) {
        //System.out.println("id number: "+getIDNumber());
      	//System.err.println("hello?");
        double distance = (fStartTime - fEndTime) * fRate;
        synchronized (this) {
            if (updateHKY) {
                setupHKY();
                //System.out.println(" Set up hky");
            }
        }

        final double xx = beta * distance;
        final double bbR = Math.exp(xx * A_R);
        final double bbY = Math.exp(xx * A_Y);

        final double aa = Math.exp(xx);
        final double oneminusa = 1 - aa;

        final double t1Aaa = (tab1A * aa);
        matrix[0] = freqA + t1Aaa + (tab2A * bbR);

        matrix[1] = freqC * oneminusa;
        final double t1Gaa = (tab1G * aa);
        matrix[2] = freqG + t1Gaa - (tab3G * bbR);
        matrix[3] = freqT * oneminusa;

        matrix[4] = freqA * oneminusa;
        final double t1Caa = (tab1C * aa);
        matrix[5] = freqC + t1Caa + (tab2C * bbY);
        matrix[6] = freqG * oneminusa;
        final double t1Taa = (tab1T * aa);
        matrix[7] = freqT + t1Taa - (tab3T * bbY);

        matrix[8] = freqA + t1Aaa - (tab3A * bbR);
        matrix[9] = matrix[1];
        matrix[10] = freqG + t1Gaa + (tab2G * bbR);
        matrix[11] = matrix[3];

        matrix[12] = matrix[4];
        matrix[13] = freqC + t1Caa - (tab3C * bbY);
        matrix[14] = matrix[6];
        matrix[15] = freqT + t1Taa + (tab2T * bbY);
    }

    protected void setupHKY() {
        if(INDICATORS[getCurrModel()][F81_INDEX] == PRESENT){
            freqA = frequencies.getValue(A);
            freqC = frequencies.getValue(C);
            freqG = frequencies.getValue(G);
            freqT = frequencies.getValue(T);
        }else{
            freqA = UNIF_DIST[0];
            freqC = UNIF_DIST[1];
            freqG = UNIF_DIST[2];
            freqT = UNIF_DIST[3];
        }


        final double freqR = freqA + freqG;
        final double freqY = freqC + freqT;

        // small speed up - reduce calculations. Comments show original code

        // (C+T) / (A+G)
        final double r1 = (1 / freqR) - 1;
        tab1A = freqA * r1;

        tab3A = freqA / freqR;
        tab2A = 1 - tab3A;        // (freqR-freqA)/freqR;

        final double r2 = 1 / r1; // ((1 / freqY) - 1);
        tab1C = freqC * r2;

        tab3C = freqC / freqY;
        tab2C = 1 - tab3C;       // (freqY-freqC)/freqY; assert  tab2C + tab3C == 1.0;

        tab1G = freqG * r1;
        tab3G = tab2A;            // 1 - tab3A; // freqG/freqR;
        tab2G = tab3A;            // 1 - tab3G; // (freqR-freqG)/freqR;

        tab1T = freqT * r2;

        tab3T = tab2C;            // 1 - tab3C;  // freqT/freqY;
        tab2T = tab3C;            // 1 - tab3T; // (freqY-freqT)/freqY; //assert tab2T + tab3T == 1.0 ;

        final double k = Math.exp(logKappa.getValue()*INDICATORS[getCurrModel()][K80_INDEX]);
        beta = -1.0 / (2.0 * (freqR * freqY + k * (freqA * freqG + freqC * freqT)));

        A_R = 1.0 + freqR * (k - 1);
        A_Y = 1.0 + freqY * (k - 1);

        updateHKY = false;
    }

    /**
     * Used for precalculations
     */

    private double p1a;
    private double p0a;
    private double p3b;
    private double p2b;
    private double a;
    private double b;
    private double p1aa;
    private double p0aa;
    private double p3bb;
    private double p2bb;
    private double p1aIsa;
    private double p0aIsa;
    private double p3bIsb;
    private double p2bIsb;
    private double k1g;
    private double k1a;
    private double k2t;
    private double k2c;
    private double subrateScale;
    public void getTN93TransitionProbabilities(Node node, double fStartTime, double fEndTime, double fRate, double[] matrix) {
        double distance = (fStartTime - fEndTime) * fRate;
        synchronized (this) {
            if (updateTN || updateHKY) {
                setupTN();
            }
        }

        distance /= subrateScale;

        double[] q = {
                0, k1g, frequencies.getValue(C), frequencies.getValue(T),
                k1a, 0, frequencies.getValue(C), frequencies.getValue(T),
                frequencies.getValue(A), frequencies.getValue(G), 0, k2t,
                frequencies.getValue(A), frequencies.getValue(G), k2c, 0
        };

        q[0] = -(q[1] + q[2] + q[3]);
        q[5] = -(q[4] + q[6] + q[7]);
        q[10] = -(q[8] + q[9] + q[11]);
        q[15] = -(q[12] + q[13] + q[14]);

        double[] fa0 = {
                1 + q[0] - p1aa, q[1] + p1aa, q[2], q[3],
                q[4] + p0aa, 1 + q[5] - p0aa, q[6], q[7],
                q[8], q[9], 1 + q[10] - p3bb, q[11] + p3bb,
                q[12], q[13], q[14] + p2bb, 1 + q[15] - p2bb
        };


        double[] fa1 = {
                -q[0] + p1aIsa, -q[1] - p1aIsa, -q[2], -q[3],
                -q[4] - p0aIsa, -q[5] + p0aIsa, -q[6], -q[7],
                -q[8], -q[9], -q[10] + p3bIsb, -q[11] - p3bIsb,
                -q[12], -q[13], -q[14] - p2bIsb, -q[15] + p2bIsb};

        double et = Math.exp(-distance);

        for (int k = 0; k < 16; ++k) {
            fa1[k] = fa1[k] * et + fa0[k];
        }

        final double eta = Math.exp(distance * a);
        final double etb = Math.exp(distance * b);

        double za = eta / (a * (1 + a));
        double zb = etb / (b * (1 + b));
        double u0 = p1a * za;
        double u1 = p0a * za;
        double u2 = p3b * zb;
        double u3 = p2b * zb;

        fa1[0] += u0;
        fa1[1] -= u0;
        fa1[4] -= u1;
        fa1[5] += u1;

        fa1[10] += u2;
        fa1[11] -= u2;
        fa1[14] -= u3;
        fa1[15] += u3;

        // transpose 2 middle rows and columns
        matrix[0] = fa1[0];
        matrix[1] = fa1[2];
        matrix[2] = fa1[1];
        matrix[3] = fa1[3];
        matrix[4] = fa1[8];
        matrix[5] = fa1[10];
        matrix[6] = fa1[9];
        matrix[7] = fa1[11];
        matrix[8] = fa1[4];
        matrix[9] = fa1[6];
        matrix[10] = fa1[5];
        matrix[11] = fa1[7];
        matrix[12] = fa1[12];
        matrix[13] = fa1[14];
        matrix[14] = fa1[13];
        matrix[15] = fa1[15];

        //System.arraycopy(fa1, 0, matrix, 0, 16);
    }


    private void setupTN() {
        freqA = frequencies.getValue(A);
        freqC = frequencies.getValue(C);
        freqG = frequencies.getValue(G);
        freqT = frequencies.getValue(T);

        double freqR = freqA + freqG;
		double freqY = freqC + freqT;

        double k1 = Math.exp(logKappa.getValue());
        double k2 = Math.exp(logKappa.getValue()+logTN.getValue());

        //System.out.println(getIDNumber() +" Using " + k1 + " " + k2);
        // A hack until I get right this boundary case. gives results accurate to 1e-8 in the P matrix
        // so should be OK even like this.
        if (k1 == 1) {
            k1 += 1E-10;
        }
        if (k2 == 1) {
            k2 += 1e-10;
        }

        double l1 = k1 * k1 * freqR + k1 * (2 * freqY - 1) - freqY;
        double l2 = k2 * k2 * freqY + k2 * (2 * freqR - 1) - freqR;

        p1a = freqG * l1;
        p0a = freqA * l1;
        p3b = freqT * l2;
        p2b = freqC * l2;

        a = -(k1 * freqR + freqY);
        b = -(k2 * freqY + freqR);

        p1aa = p1a / a;
        p0aa = p0a / a;
        p3bb = p3b / b;
        p2bb = p2b / b;

        p1aIsa = p1a / (1 + a);
        p0aIsa = p0a / (1 + a);
        p3bIsb = p3b / (1 + b);
        p2bIsb = p2b / (1 + b);

        k1g = k1 * freqG;
        k1a = k1 * freqA;
        k2t = k2 * freqT;
        k2c = k2 * freqC;

        subrateScale = 2 * (k1 * freqA * freqG + k2 * freqC * freqT + freqR * freqY);
        // updateMatrix = true;
        updateTN = false;
    }

    public void getGTRTransitionProbabilities(Node node, double fStartTime, double fEndTime, double fRate, double[] matrix) {
        double distance = (fStartTime - fEndTime) * fRate;

        int i, j, k;
        double temp;

        // this must be synchronized to avoid being called simultaneously by
        // two different likelihood threads - AJD
        synchronized (this) {
            if (updateMatrix) {
                //System.out.println("UPDATE MATRIX");
            	setupRelativeRates();
            	setupRateMatrix();

            	updateMatrix = false;
            }
        }

        if (!wellConditioned) {
            Arrays.fill(matrix, 0.0);
            return;
        }
        double[][] iexp = new double[STATE_COUNT][STATE_COUNT];//popiexp();




        for (i = 0; i < STATE_COUNT; i++) {
            temp = Math.exp(distance * Eval[i]);
            for (j = 0; j < STATE_COUNT; j++) {
                iexp[i][j] = Ievc[i][j] * temp;
            }
        }


        int u = 0;
        for (i = 0; i < STATE_COUNT; i++) {
            for (j = 0; j < STATE_COUNT; j++) {
                temp = 0.0;
                for (k = 0; k < STATE_COUNT; k++) {
                    temp += Evec[i][k] * iexp[k][j];
                }
                if (temp < 0.0)
                    matrix[u] = minProb;
                else
                    matrix[u] = temp;
                u++;
            }
        }

    } // getTransitionProbabilities

    public boolean requiresRecalculation(){
        boolean recalculate = false;
        //System.err.println("model "+getCurrModel());
        //System.err.println("frequencies.somethingIsDirty() "+frequencies.somethingIsDirty());

        if(modelChoose.somethingIsDirty()){
            if(modelChoose.getValue() > JC69 && modelChoose.getValue() < GTR){
                updateHKY = true;
            }
             recalculate = true;
        }
        if(frequencies.somethingIsDirty()){
            updateHKY = true;
            updateTN = true;
            if(INDICATORS[getCurrModel()][F81_INDEX] == PRESENT){
                recalculate = true;
            }

        }
        if(logKappa.somethingIsDirty()){
            updateHKY = true;
            updateTN = true;
            if(INDICATORS[getCurrModel()][K80_INDEX] == PRESENT){
                recalculate = true;
            }

        }

        if(logTN.somethingIsDirty()){
            updateTN = true;
            if(INDICATORS[getCurrModel()][TN_INDEX] == PRESENT){
                recalculate = true;
            }

        }

        if(logAC.somethingIsDirty() ||
                logAT.somethingIsDirty() ||
                logGC.somethingIsDirty()&&
               INDICATORS[getCurrModel()][GTR_INDEX] == PRESENT){

            recalculate = true;

        }

        if(recalculate){
            updateMatrix = true;
        }
        /*System.out.println("updateHKY: "+updateHKY);
        System.out.println("recalculate: "+recalculate+" "+getCurrModel());
        System.out.println("log kappa: "+logKappa.getValue());
        System.out.println("log tn: "+logTN.getValue());
        System.out.println("log ac: "+logAC.getValue());
        System.out.println("log at: "+logAT.getValue());
        System.out.println("log at: "+logGC.getValue());
        System.out.println("log at: "+modelChoose.getValue());
        System.out.println("freqs: "+frequencies.getValue(0)+" "+frequencies.getValue(1)+" "+
            frequencies.getValue(2)+" "+frequencies.getValue(3));*/
        return recalculate;
    }

    /*public void store(){
        storedUpdateHKY = updateHKY;
        storedUpdateTN = updateTN;
        super.store();
    }*/

    public void restore(){
        //System.out.println("storedUpdateHKY: "+storedUpdateHKY);
        updateHKY = true;
         //storedUpdateHKY;
        updateTN = true;
                // storedUpdateTN;
        super.restore();
        //updateMatrix= true;
    }


    public void setUpdateMatrix(boolean update){
        updateMatrix = update;
        updateHKY = update;
        updateTN = update;

    }


}
