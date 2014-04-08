package beast.evolution.substitutionmodel;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.QuietRealParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.datatype.DataType;
import beast.evolution.datatype.Nucleotide;
import beast.evolution.substitutionmodel.DefaultEigenSystem;
import beast.evolution.substitutionmodel.EigenDecomposition;
import beast.evolution.substitutionmodel.EigenSystem;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.Node;
import beast.math.MachineAccuracy;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.Property;

import java.util.Arrays;

/**
 * @author Chieh-Hsi Wu
 */
@Description("This class augments the parameter space of GTR to allow selection of nucleotide substitution models.")
public class NtdBMA extends SubstitutionModel.Base{




    public Input<QuietRealParameter> logKappaInput = new Input<QuietRealParameter>("logKappa", "parameter representing log of HKY kappa parameter", Input.Validate.REQUIRED);
    public Input<QuietRealParameter> logTNInput = new Input<QuietRealParameter>("logTN", "parameter representing log of TN parameter", Input.Validate.REQUIRED);
    public Input<QuietRealParameter> logACInput = new Input<QuietRealParameter>("logAC", "parameter representing log of AC parameter", Input.Validate.REQUIRED);
    public Input<QuietRealParameter> logATInput = new Input<QuietRealParameter>("logAT", "parameter representing log of AT parameter", Input.Validate.REQUIRED);
    public Input<QuietRealParameter> logGCInput = new Input<QuietRealParameter>("logGC", "parameter representing log of GC parameter", Input.Validate.REQUIRED);
    //public Input<RealParameter> logGTInput = new Input<RealParameter>("logGT", "parameter representing log of GT parameter", Input.Validate.REQUIRED);
    public Input<QuietRealParameter> modelChooseInput = new Input<QuietRealParameter>("modelChoose", "Integer presenting the model", Input.Validate.REQUIRED);
    public Input<QuietRealParameter> freqInput = new Input<QuietRealParameter>("frequenciesParameter", "Stationary frequencies the model", Input.Validate.REQUIRED);

    QuietRealParameter logKappa;
    QuietRealParameter logTN;
    QuietRealParameter logAC;
    QuietRealParameter logAT;
    QuietRealParameter logGC;
    //public RealParameter logGT;
    QuietRealParameter modelChoose;
    QuietRealParameter frequencies;

    public static final int NTD_PAIRS = 16;


    public static final int STATE_COUNT = 4;
    public static final int RATE_COUNT = 6;
    public static final int A = 0;
    public static final int C = 1;
    public static final int G = 2;
    public static final int T = 3;



    public static final int ABSENT = 0;
    public static final int PRESENT = 1;

    public static final int K80_INDEX = 0;
    public static final int F81_INDEX = 1;
    public static final int TN_INDEX = 2;
    public static final int GTR_INDEX = 3;

    public static final Double[] UNIF_DIST = {1.0/STATE_COUNT,1.0/STATE_COUNT,1.0/STATE_COUNT,1.0/STATE_COUNT};

    public static final int JC69 = 0;
    public static final int K80 = 1;
    public static final int F81 = 2;
    public static final int HKY85 = 3;
    public static final int TN93 = 4;
    public static final int GTR = 5;

    public static final int[][] INDICATORS = {
            {ABSENT, ABSENT, ABSENT,ABSENT},
            {PRESENT, ABSENT, ABSENT,ABSENT},
            {ABSENT, PRESENT,ABSENT,ABSENT},
            {PRESENT, PRESENT,ABSENT,ABSENT},
            {PRESENT, PRESENT, PRESENT,ABSENT},
            {PRESENT, PRESENT, PRESENT,PRESENT}
    };


    private int IDNumber = -1;
    public void setIDNumber(int idNum){
        IDNumber = idNum;
    }

    public int getIDNumber(){
        return IDNumber;
    }

    double [][] m_rateMatrix;
    protected double[] relativeRates;
    protected double[] storedRelativeRates;

    private int[] ordr;
    private double[] evali;

    public NtdBMA(){frequenciesInput.setRule(Input.Validate.OPTIONAL);}

    public NtdBMA(
            QuietRealParameter logKappa,
            QuietRealParameter logTN,
            QuietRealParameter logAC,
            QuietRealParameter logAT,
            QuietRealParameter logGC,
            //RealParameter logGT,
            QuietRealParameter modelChoose,
            QuietRealParameter frequencies){


        frequenciesInput.setRule(Input.Validate.OPTIONAL);

        initialize(logKappa,logTN,logAC,logAT,logGC, modelChoose, frequencies);


    }


    @Override
    public void initAndValidate() throws Exception {
        initialize(
                logKappaInput.get(),
                logTNInput.get(),
                logACInput.get(),
                logATInput.get(),
                logGCInput.get(),
                modelChooseInput.get(),
                freqInput.get());


        //q = new double[STATE_COUNT][STATE_COUNT];
    } // initAndValidate

    public void initialize(QuietRealParameter logKappa,
            QuietRealParameter logTN,
            QuietRealParameter logAC,
            QuietRealParameter logAT,
            QuietRealParameter logGC,
            //RealParameter logGT,
            QuietRealParameter modelChoose,
            QuietRealParameter frequencies){
        this.logKappa = logKappa;
        this.logTN = logTN;
        this.logAC = logAC;
        this.logAT = logAT;
        this.logGC = logGC;
        //this.logGT = logGT;
        this.frequencies = frequencies;
        this.modelChoose = modelChoose;
        if(modelChoose.getUpper() > (double)GTR || modelChoose.getLower() < (double)JC69){
            //throw new RuntimeException("Provided model choose value is "+modelChoose.getValue()+".\n The value of model choose needs to be between " + JC69 + " and " + GTR + "inclusive, " +
                    //"where "+ JC69 + " and " + GTR +" represents JC and GTR repectively");
        }
        updateMatrix = true;

        //eigenSystem = new DefaultEigenSystem(STATE_COUNT);
        m_rateMatrix = new double[STATE_COUNT][STATE_COUNT];
        relativeRates = new double[RATE_COUNT];
        storedRelativeRates = new double[RATE_COUNT];
        ordr = new int[STATE_COUNT];
        evali = new double[STATE_COUNT];
        initialiseEigen();

    }





    protected void setupRelativeRates() {

        //rate AG value
    	relativeRates[1] = 1.0;


        //rate CT value
        relativeRates[4] = Math.exp(INDICATORS[getCurrModel()][TN_INDEX]*logTN.getValue());

        //rate AC value
        relativeRates[0] = Math.exp(
                0.0-INDICATORS[getCurrModel()][K80_INDEX]*logKappa.getValue()+
                INDICATORS[getCurrModel()][GTR_INDEX]*logAC.getValue());

        //rate AT value
        relativeRates[2] = Math.exp(
                -INDICATORS[getCurrModel()][K80_INDEX]*logKappa.getValue()+
                        INDICATORS[getCurrModel()][GTR_INDEX]*logAT.getValue());

        //rate GC value
        relativeRates[3] = Math.exp(
                -INDICATORS[getCurrModel()][K80_INDEX]*logKappa.getValue()+
                        INDICATORS[getCurrModel()][GTR_INDEX]*logGC.getValue());

        //rate GT value
        relativeRates[5] = Math.exp(-INDICATORS[getCurrModel()][K80_INDEX]*logKappa.getValue());


        /*System.out.println("AC: "+relativeRates[0]);
        System.out.println("AG: "+relativeRates[1]);
        System.out.println("AT: "+relativeRates[2]);
        System.out.println("GC: "+relativeRates[3]);
        System.out.println("CT: "+relativeRates[4]);
        System.out.println("GT: "+relativeRates[5]);
        System.out.println("indicators: "+INDICATORS[getCurrModel()][GTR_INDEX]);

        System.out.println(logKappa.getValue());

        System.out.println(logTN.getValue());
        System.out.println(logAT.getValue());
        System.out.println(logAC.getValue());
        System.out.println(logGC.getValue());
        System.out.println(frequencies.getValue(0)+" "+
        frequencies.getValue(1)+" "+
        frequencies.getValue(2)+" "+
        frequencies.getValue(3));
        System.out.println(modelChoose.getValue());*/
    }

    /** sets up rate matrix **/
    protected void setupRateMatrix() {
    	Double [] fFreqs;

        if(INDICATORS[getCurrModel()][F81_INDEX] == PRESENT){
            fFreqs = frequencies.getValues();
        }else{
            fFreqs = UNIF_DIST;
        }
        //System.err.println("getCurrModel(): "+getCurrModel());
        //System.err.println("freq0: "+fFreqs[0]);


        int i, j, k = 0;

        // Set the instantaneous rate matrix
        for (i = 0; i < STATE_COUNT; i++) {
            m_rateMatrix[i][i] = 0;

            for (j = i + 1; j < STATE_COUNT; j++) {
                m_rateMatrix[i][j] = relativeRates[k] * fFreqs[j];
                m_rateMatrix[j][i] = relativeRates[k] * fFreqs[i];
                k += 1;
            }
        }


        // set up diagonal
        for (i = 0; i < STATE_COUNT; i++) {
            double fSum = 0.0;
            for (j = 0; j < STATE_COUNT; j++) {
                if (i != j)
                    fSum += m_rateMatrix[i][j];
            }
            m_rateMatrix[i][i] = -fSum;
        }
        // normalise rate matrix to one expected substitution per unit time
        double fSubst = 0.0;
        for (i = 0; i < STATE_COUNT; i++)
            fSubst += -m_rateMatrix[i][i] * fFreqs[i];



       /*System.err.println("Part 2");
       for(i = 0; i < STATE_COUNT; i++){

                System.err.print(m_rateMatrix[i][i]+" ");

            System.err.println();
        }
        System.err.println("fSubst: "+fSubst);*/
        for (i = 0; i < STATE_COUNT; i++) {
            for (j = 0; j < STATE_COUNT; j++) {
            	m_rateMatrix[i][j] = m_rateMatrix[i][j] / fSubst;
                //System.err.print(m_rateMatrix[i][j]+" ");
            }
            //System.err.println();
        }

        /*for (i = 0; i < STATE_COUNT; i++) {
            for (j = 0; j < STATE_COUNT; j++) {
            	m_rateMatrix[i][j] = UsefulShit.truncate(m_rateMatrix[i][j] / fSubst,10);
            }
        }*/
        try{
            elmhes(m_rateMatrix, ordr, STATE_COUNT);
            eltran(m_rateMatrix, Evec, ordr, STATE_COUNT);
            hqr2(STATE_COUNT, 1, STATE_COUNT, m_rateMatrix, Evec, Eval, evali);
            luinverse(Evec, Ievc, STATE_COUNT);


            // Check for valid decomposition
            for (i = 0; i < STATE_COUNT; i++) {
                if (Double.isNaN(Eval[i]) || Double.isInfinite(Eval[i]) ) {
                    wellConditioned = false;
                    return;
                }
            }
        }catch(ArithmeticException e){
            wellConditioned = false;
            return;

        }catch(IllegalArgumentException e){
            wellConditioned = false;
            return;

        }

        updateMatrix = false;
        //wellConditioned = true;

	} // setupRateMatrix

    public boolean requiresRecalculation(){
        boolean recalculate = false;
        //System.err.println("model "+getCurrModel());
        //System.err.println("frequencies.somethingIsDirty() "+frequencies.somethingIsDirty());

        if(modelChoose.somethingIsDirty()){
            recalculate = true;
        }else if(frequencies.somethingIsDirty() &&
                INDICATORS[getCurrModel()][F81_INDEX] == PRESENT){
            //System.err.println(frequencies);
            recalculate = true;

        }else if(logKappa.somethingIsDirty() &&
                INDICATORS[getCurrModel()][K80_INDEX] == PRESENT){

            recalculate = true;

        }else if(logTN.somethingIsDirty() &&
                INDICATORS[getCurrModel()][TN_INDEX] == PRESENT){

            recalculate = true;

        }else if(logAC.somethingIsDirty() &&
               INDICATORS[getCurrModel()][GTR_INDEX] == PRESENT){

            recalculate = true;

        }else if(logAT.somethingIsDirty() &&
                INDICATORS[getCurrModel()][GTR_INDEX] == PRESENT){

            recalculate = true;

        }else if(logGC.somethingIsDirty() &&
                INDICATORS[getCurrModel()][GTR_INDEX] == PRESENT){

            recalculate = true;

        }/*else if(logGT.somethingIsDirty() &&
                INDICATORS[getCurrModel()][GTR_INDEX] == PRESENT){

            recalculate = true;

        }*/
        if(recalculate){
            updateMatrix = true;
        }
        return recalculate;
    }

    public void setUpdateMatrix(boolean update){
        updateMatrix = update;

    }

    protected int getCurrModel(){
        return (int)((double)modelChoose.getValue());

    }

    @Override
    public double[] getFrequencies() {
        Double[] temp;
        if(INDICATORS[getCurrModel()][F81_INDEX] == PRESENT){
            //System.out.println("estimate freqs");
            temp =  frequencies.getValues();
        }else{

            temp =  UNIF_DIST;
        }

        double[] freqs = new double[temp.length];
        for(int i = 0; i < freqs.length;i++){
            freqs[i] = temp[i];
        }
        return freqs;
    }

        @Override
    public void getTransitionProbabilities(Node node, double fStartTime, double fEndTime, double fRate, double[] matrix) {
    	//System.err.println("Get probs: "+updateMatrix);
            //System.out.println("Get probs: "+frequencies);
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
            System.err.println("THIS IS BOTHERSOME");
            System.err.println("distance: "+distance);
            printDetails();
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

    public QuietRealParameter getLogKappa(){
        return logKappa;
    }

    public QuietRealParameter getLogTN(){
        return logTN;
    }

    public QuietRealParameter getLogAC(){
        return logAC;
    }

    public QuietRealParameter getLogAT(){
        return logAT;
    }

    public QuietRealParameter getLogGC(){
        return logGC;
    }

    public QuietRealParameter getModelChoose(){
        return modelChoose;
    }

    public QuietRealParameter getFreqs(){
        return frequencies;
    }

    protected double[] Eval;
    protected double[] storedEval;
    protected double[][] Evec;
    protected double[][] storedEvec;
    protected double[][] Ievc;
    protected double[][] storedIevc;

    /**
     * allocate memory for the Eigen routines
     */
    protected void initialiseEigen() {

        Eval = new double[STATE_COUNT];
        Evec = new double[STATE_COUNT][STATE_COUNT];
        Ievc = new double[STATE_COUNT][STATE_COUNT];

        storedEval = new double[STATE_COUNT];
        storedEvec = new double[STATE_COUNT][STATE_COUNT];
        storedIevc = new double[STATE_COUNT][STATE_COUNT];



        updateMatrix = true;
    }



    /**
     * Restore the additional stored state
     */
    @Override

    public void restore() {

        updateMatrix = storedUpdateMatrix;
        wellConditioned = storedWellConditioned;

        double[] tmp1 = storedEval;
        storedEval = Eval;
        Eval = tmp1;

        double[][] tmp2 = storedIevc;
        storedIevc = Ievc;
        Ievc = tmp2;

        tmp2 = storedEvec;
        storedEvec = Evec;
        Evec = tmp2;

        // To restore all this stuff just swap the pointers...
        double[] tmp4 = storedRelativeRates;
        storedRelativeRates = relativeRates;
        relativeRates = tmp4;

        super.restore();

    }

    public void store() {

        storedUpdateMatrix = updateMatrix;

//        if(updateMatrix)
//            System.err.println("Storing updatable state!");

        storedWellConditioned = wellConditioned;



        // Inherited

        System.arraycopy(Eval, 0, storedEval, 0, STATE_COUNT);
        for(int i = 0; i < Ievc.length;i++){
            System.arraycopy(Ievc[i], 0, storedIevc[i], 0, STATE_COUNT);
        }
        for(int i = 0;i < Evec.length; i++){
            System.arraycopy(Evec[i], 0, storedEvec[i], 0, STATE_COUNT);
        }

        super.store();
    }

	@Override
	public boolean canHandleDataType(DataType dataType) throws Exception {
		if (dataType instanceof Nucleotide) {
			return true;
		}
		throw new Exception("Can only handle nucleotide data");
	}
    /**
     * This function returns the Eigen vectors.
     *
     * @return the array
     */
    @Override
    public EigenDecomposition getEigenDecomposition(Node node) {

        EigenSystem eigenSystem =  new DefaultEigenSystem(STATE_COUNT);
        synchronized (this) {
            if (updateMatrix) {
            	setupRelativeRates();
                setupRateMatrix();

                updateMatrix = false;
            }
        }
        //EigenDecomposition eigenDecomposition = eigenSystem.decomposeMatrix(m_rateMatrix);
        return eigenSystem.decomposeMatrix(m_rateMatrix);
    }

    protected boolean updateMatrix = true;
    protected boolean storedUpdateMatrix = true;
    protected boolean wellConditioned = true;
    private boolean storedWellConditioned;
    protected static final double minProb = Property.DEFAULT.tolerance();

    private void elmhes(double[][] a, int[] ordr, int n) {
        int m, j, i;
        double y, x;

        for (i = 0; i < n; i++) {
            ordr[i] = 0;
        }
        for (m = 2; m < n; m++) {
            x = 0.0;
            i = m;
            for (j = m; j <= n; j++) {
                if (Math.abs(a[j - 1][m - 2]) > Math.abs(x)) {
                    x = a[j - 1][m - 2];
                    i = j;
                }
            }
            ordr[m - 1] = i;
            if (i != m) {
                for (j = m - 2; j < n; j++) {
                    y = a[i - 1][j];
                    a[i - 1][j] = a[m - 1][j];
                    a[m - 1][j] = y;
                }
                for (j = 0; j < n; j++) {
                    y = a[j][i - 1];
                    a[j][i - 1] = a[j][m - 1];
                    a[j][m - 1] = y;
                }
            }
            if (x != 0.0) {
                for (i = m; i < n; i++) {
                    y = a[i][m - 2];
                    if (y != 0.0) {
                        y /= x;
                        a[i][m - 2] = y;
                        for (j = m - 1; j < n; j++) {
                            a[i][j] -= y * a[m - 1][j];
                        }
                        for (j = 0; j < n; j++) {
                            a[j][m - 1] += y * a[j][i];
                        }
                    }
                }
            }
        }
    }

    // Helper variables for mcdiv
    private double cr, ci;

    private void mcdiv(double ar, double ai, double br, double bi) {
        double s, ars, ais, brs, bis;

        s = Math.abs(br) + Math.abs(bi);
        ars = ar / s;
        ais = ai / s;
        brs = br / s;
        bis = bi / s;
        s = brs * brs + bis * bis;
        cr = (ars * brs + ais * bis) / s;
        ci = (ais * brs - ars * bis) / s;
    }

    void hqr2(int n, int low, int hgh, double[][] h, double[][] zz,
              double[] wr, double[] wi) throws ArithmeticException {
        int i, j, k, l = 0, m, en, na, itn, its;
        double p = 0, q = 0, r = 0, s = 0, t, w, x = 0, y, ra, sa, vi, vr, z = 0, norm, tst1, tst2;
        boolean notLast;


        norm = 0.0;
        k = 1;
        /* store isolated roots and compute matrix norm */
        for (i = 0; i < n; i++) {
            for (j = k - 1; j < n; j++) {
                norm += Math.abs(h[i][j]);
            }
            k = i + 1;
            if (i + 1 < low || i + 1 > hgh) {
                wr[i] = h[i][i];
                wi[i] = 0.0;
            }
        }
        en = hgh;
        t = 0.0;
        itn = n * 30;
        while (en >= low) {    /* search for next eigenvalues */
            its = 0;
            na = en - 1;
            while (en >= 1) {
                /* look for single small sub-diagonal element */
                boolean fullLoop = true;
                for (l = en; l > low; l--) {
                    s = Math.abs(h[l - 2][l - 2]) + Math.abs(h[l - 1][l - 1]);
                    if (s == 0.0) {
                        s = norm;
                    }
                    tst1 = s;
                    tst2 = tst1 + Math.abs(h[l - 1][l - 2]);
                    if (tst2 == tst1) {
                        fullLoop = false;
                        break;
                    }
                }
                if (fullLoop) {
                    l = low;
                }

                x = h[en - 1][en - 1];    /* form shift */
                if (l == en || l == na) {
                    break;
                }
                if (itn == 0) {
                    /* eigenvalues have not converged */
                    System.out.println("WARNING: Eigenvalues have not converged.");
                    printDetails();
                    throw new ArithmeticException();
                }
                y = h[na - 1][na - 1];
                w = h[en - 1][na - 1] * h[na - 1][en - 1];
                /* form exceptional shift */
                if (its == 10 || its == 20) {
                    t += x;
                    for (i = low - 1; i < en; i++) {
                        h[i][i] -= x;
                    }
                    s = Math.abs(h[en - 1][na - 1]) + Math.abs(h[na - 1][en - 3]);
                    x = 0.75 * s;
                    y = x;
                    w = -0.4375 * s * s;
                }
                its++;
                itn--;
                /* look for two consecutive small sub-diagonal elements */
                for (m = en - 2; m >= l; m--) {
                    z = h[m - 1][m - 1];
                    r = x - z;
                    s = y - z;
                    p = (r * s - w) / h[m][m - 1] + h[m - 1][m];
                    q = h[m][m] - z - r - s;
                    r = h[m + 1][m];
                    s = Math.abs(p) + Math.abs(q) + Math.abs(r);
                    p /= s;
                    q /= s;
                    r /= s;
                    if (m == l) {
                        break;
                    }
                    tst1 = Math.abs(p) * (Math.abs(h[m - 2][m - 2]) + Math.abs(z) + Math.abs(h[m][m]));
                    tst2 = tst1 + Math.abs(h[m - 1][m - 2]) * (Math.abs(q) + Math.abs(r));
                    if (tst2 == tst1) {
                        break;
                    }
                }
                for (i = m + 2; i <= en; i++) {
                    h[i - 1][i - 3] = 0.0;
                    if (i != m + 2) {
                        h[i - 1][i - 4] = 0.0;
                    }
                }
                for (k = m; k <= na; k++) {
                    notLast = k != na;
                    if (k != m) {
                        p = h[k - 1][k - 2];
                        q = h[k][k - 2];
                        r = 0.0;
                        if (notLast) {
                            r = h[k + 1][k - 2];
                        }
                        x = Math.abs(p) + Math.abs(q) + Math.abs(r);
                        if (x != 0.0) {
                            p /= x;
                            q /= x;
                            r /= x;
                        }
                    }
                    if (x != 0.0) {
                        if (p < 0.0) {    /* sign */
                            s = -Math.sqrt(p * p + q * q + r * r);
                        } else {
                            s = Math.sqrt(p * p + q * q + r * r);
                        }
                        if (k != m) {
                            h[k - 1][k - 2] = -s * x;
                        } else if (l != m) {
                            h[k - 1][k - 2] = -h[k - 1][k - 2];
                        }
                        p += s;
                        x = p / s;
                        y = q / s;
                        z = r / s;
                        q /= p;
                        r /= p;
                        if (!notLast) {
                            for (j = k - 1; j < n; j++) {    /* row modification */
                                p = h[k - 1][j] + q * h[k][j];
                                h[k - 1][j] -= p * x;
                                h[k][j] -= p * y;
                            }
                            j = (en < (k + 3)) ? en : (k + 3); /* min */
                            for (i = 0; i < j; i++) {    /* column modification */
                                p = x * h[i][k - 1] + y * h[i][k];
                                h[i][k - 1] -= p;
                                h[i][k] -= p * q;
                            }
                            /* accumulate transformations */
                            for (i = low - 1; i < hgh; i++) {
                                p = x * zz[i][k - 1] + y * zz[i][k];
                                zz[i][k - 1] -= p;
                                zz[i][k] -= p * q;
                            }
                        } else {
                            for (j = k - 1; j < n; j++) {    /* row modification */
                                p = h[k - 1][j] + q * h[k][j] + r * h[k + 1][j];
                                h[k - 1][j] -= p * x;
                                h[k][j] -= p * y;
                                h[k + 1][j] -= p * z;
                            }
                            j = (en < (k + 3)) ? en : (k + 3); /* min */
                            for (i = 0; i < j; i++) {    /* column modification */
                                p = x * h[i][k - 1] + y * h[i][k] + z * h[i][k + 1];
                                h[i][k - 1] -= p;
                                h[i][k] -= p * q;
                                h[i][k + 1] -= p * r;
                            }
                            /* accumulate transformations */
                            for (i = low - 1; i < hgh; i++) {
                                p = x * zz[i][k - 1] + y * zz[i][k] +
                                        z * zz[i][k + 1];
                                zz[i][k - 1] -= p;
                                zz[i][k] -= p * q;
                                zz[i][k + 1] -= p * r;
                            }
                        }
                    }
                }                 /* for k */
            }                     /* while infinite loop */
            if (l == en) {                 /* one root found */
                h[en - 1][en - 1] = x + t;
                wr[en - 1] = h[en - 1][en - 1];
                wi[en - 1] = 0.0;
                en = na;
                continue;
            }
            y = h[na - 1][na - 1];
            w = h[en - 1][na - 1] * h[na - 1][en - 1];
            p = (y - x) / 2.0;
            q = p * p + w;
            z = Math.sqrt(Math.abs(q));
            h[en - 1][en - 1] = x + t;
            x = h[en - 1][en - 1];
            h[na - 1][na - 1] = y + t;
            if (q >= 0.0) {     /* real pair */
                if (p < 0.0) {    /* sign */
                    z = p - Math.abs(z);
                } else {
                    z = p + Math.abs(z);
                }
                wr[na - 1] = x + z;
                wr[en - 1] = wr[na - 1];
                if (z != 0.0) {
                    wr[en - 1] = x - w / z;
                }
                wi[na - 1] = 0.0;
                wi[en - 1] = 0.0;
                x = h[en - 1][na - 1];
                s = Math.abs(x) + Math.abs(z);
                p = x / s;
                q = z / s;
                r = Math.sqrt(p * p + q * q);
                p /= r;
                q /= r;
                for (j = na - 1; j < n; j++) {    /* row modification */
                    z = h[na - 1][j];
                    h[na - 1][j] = q * z + p * h[en - 1][j];
                    h[en - 1][j] = q * h[en - 1][j] - p * z;
                }
                for (i = 0; i < en; i++) {    /* column modification */
                    z = h[i][na - 1];
                    h[i][na - 1] = q * z + p * h[i][en - 1];
                    h[i][en - 1] = q * h[i][en - 1] - p * z;
                }
                /* accumulate transformations */
                for (i = low - 1; i < hgh; i++) {
                    z = zz[i][na - 1];
                    zz[i][na - 1] = q * z + p * zz[i][en - 1];
                    zz[i][en - 1] = q * zz[i][en - 1] - p * z;
                }
            } else {    /* complex pair */
                wr[na - 1] = x + p;
                wr[en - 1] = x + p;
                wi[na - 1] = z;
                wi[en - 1] = -z;
            }
            en -= 2;
        } /* while en >= low */
        /* backsubstitute to find vectors of upper triangular form */
        if (norm != 0.0) {
            for (en = n; en >= 1; en--) {
                p = wr[en - 1];
                q = wi[en - 1];
                na = en - 1;
                if (q == 0.0) {/* real vector */
                    m = en;
                    h[en - 1][en - 1] = 1.0;
                    if (na != 0) {
                        for (i = en - 2; i >= 0; i--) {
                            w = h[i][i] - p;
                            r = 0.0;
                            for (j = m - 1; j < en; j++) {
                                r += h[i][j] * h[j][en - 1];
                            }
                            if (wi[i] < 0.0) {
                                z = w;
                                s = r;
                            } else {
                                m = i + 1;
                                if (wi[i] == 0.0) {
                                    t = w;
                                    if (t == 0.0) {
                                        tst1 = norm;
                                        t = tst1;
                                        do {
                                            t = 0.01 * t;
                                            tst2 = norm + t;
                                        }
                                        while (tst2 > tst1);
                                    }
                                    h[i][en - 1] = -(r / t);
                                } else {    /* solve real equations */
                                    x = h[i][i + 1];
                                    y = h[i + 1][i];
                                    q = (wr[i] - p) * (wr[i] - p) + wi[i] * wi[i];
                                    t = (x * s - z * r) / q;
                                    h[i][en - 1] = t;
                                    if (Math.abs(x) > Math.abs(z))
                                        h[i + 1][en - 1] = (-r - w * t) / x;
                                    else
                                        h[i + 1][en - 1] = (-s - y * t) / z;
                                }
                                /* overflow control */
                                t = Math.abs(h[i][en - 1]);
                                if (t != 0.0) {
                                    tst1 = t;
                                    tst2 = tst1 + 1.0 / tst1;
                                    if (tst2 <= tst1) {
                                        for (j = i; j < en; j++) {
                                            h[j][en - 1] /= t;
                                        }
                                    }
                                }
                            }
                        }
                    }
                } else if (q > 0.0) {
                    m = na;
                    if (Math.abs(h[en - 1][na - 1]) > Math.abs(h[na - 1][en - 1])) {
                        h[na - 1][na - 1] = q / h[en - 1][na - 1];
                        h[na - 1][en - 1] = (p - h[en - 1][en - 1]) / h[en - 1][na - 1];
                    } else {
                        mcdiv(0.0, -h[na - 1][en - 1], h[na - 1][na - 1] - p, q);
                        h[na - 1][na - 1] = cr;
                        h[na - 1][en - 1] = ci;
                    }
                    h[en - 1][na - 1] = 0.0;
                    h[en - 1][en - 1] = 1.0;
                    if (en != 2) {
                        for (i = en - 3; i >= 0; i--) {
                            w = h[i][i] - p;
                            ra = 0.0;
                            sa = 0.0;
                            for (j = m - 1; j < en; j++) {
                                ra += h[i][j] * h[j][na - 1];
                                sa += h[i][j] * h[j][en - 1];
                            }
                            if (wi[i] < 0.0) {
                                z = w;
                                r = ra;
                                s = sa;
                            } else {
                                m = i + 1;
                                if (wi[i] == 0.0) {
                                    mcdiv(-ra, -sa, w, q);
                                    h[i][na - 1] = cr;
                                    h[i][en - 1] = ci;
                                } else {    /* solve complex equations */
                                    x = h[i][i + 1];
                                    y = h[i + 1][i];
                                    vr = (wr[i] - p) * (wr[i] - p);
                                    vr = vr + wi[i] * wi[i] - q * q;
                                    vi = (wr[i] - p) * 2.0 * q;
                                    if (vr == 0.0 && vi == 0.0) {
                                        tst1 = norm * (Math.abs(w) + Math.abs(q) + Math.abs(x) +
                                                Math.abs(y) + Math.abs(z));
                                        vr = tst1;
                                        do {
                                            vr = 0.01 * vr;
                                            tst2 = tst1 + vr;
                                        }
                                        while (tst2 > tst1);
                                    }
                                    mcdiv(x * r - z * ra + q * sa, x * s - z * sa - q * ra, vr, vi);
                                    h[i][na - 1] = cr;
                                    h[i][en - 1] = ci;
                                    if (Math.abs(x) > Math.abs(z) + Math.abs(q)) {
                                        h[i + 1]
                                                [na - 1] = (q * h[i][en - 1] -
                                                w * h[i][na - 1] - ra) / x;
                                        h[i + 1][en - 1] = (-sa - w * h[i][en - 1] -
                                                q * h[i][na - 1]) / x;
                                    } else {
                                        mcdiv(-r - y * h[i][na - 1], -s - y * h[i][en - 1], z, q);
                                        h[i + 1][na - 1] = cr;
                                        h[i + 1][en - 1] = ci;
                                    }
                                }
                                /* overflow control */
                                t = (Math.abs(h[i][na - 1]) > Math.abs(h[i][en - 1])) ?
                                        Math.abs(h[i][na - 1]) : Math.abs(h[i][en - 1]);
                                if (t != 0.0) {
                                    tst1 = t;
                                    tst2 = tst1 + 1.0 / tst1;
                                    if (tst2 <= tst1) {
                                        for (j = i; j < en; j++) {
                                            h[j][na - 1] /= t;
                                            h[j][en - 1] /= t;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            /* end back substitution. vectors of isolated roots */
            for (i = 0; i < n; i++) {
                if (i + 1 < low || i + 1 > hgh) {
                    for (j = i; j < n; j++) {
                        zz[i][j] = h[i][j];
                    }
                }
            }
            /* multiply by transformation matrix to give vectors of
                               * original full matrix. */
            for (j = n - 1; j >= low - 1; j--) {
                m = ((j + 1) < hgh) ? (j + 1) : hgh; /* min */
                for (i = low - 1; i < hgh; i++) {
                    z = 0.0;
                    for (k = low - 1; k < m; k++) {
                        z += zz[i][k] * h[k][j];
                    }
                    zz[i][j] = z;
                }
            }
        }
    }

    private void eltran(double[][] a, double[][] zz, int[] ordr, int n) {
        int i, j, m;

        for (i = 0; i < n; i++) {
            for (j = i + 1; j < n; j++) {
                zz[i][j] = 0.0;
                zz[j][i] = 0.0;
            }
            zz[i][i] = 1.0;
        }
        if (n <= 2) {
            return;
        }
        for (m = n - 1; m >= 2; m--) {
            for (i = m; i < n; i++) {
                zz[i][m - 1] = a[i][m - 2];
            }
            i = ordr[m - 1];
            if (i != m) {
                for (j = m - 1; j < n; j++) {
                    zz[m - 1][j] = zz[i - 1][j];
                    zz[i - 1][j] = 0.0;
                }
                zz[i - 1][m - 1] = 1.0;
            }
        }
    }

    void luinverse(double[][] inmat, double[][] imtrx, int size) throws IllegalArgumentException {
        int i, j, k, l, maxi = 0, idx, ix, jx;
        double sum, tmp, maxb, aw;
        int[] index;
        double[] wk;
        double[][] omtrx;


        index = new int[size];
        omtrx = new double[size][size];

        /* copy inmat to omtrx */
        for (i = 0; i < size; i++) {
            for (j = 0; j < size; j++) {
                omtrx[i][j] = inmat[i][j];
            }
        }

        wk = new double[size];
        aw = 1.0;
        for (i = 0; i < size; i++) {
            maxb = 0.0;
            for (j = 0; j < size; j++) {
                if (Math.abs(omtrx[i][j]) > maxb) {
                    maxb = Math.abs(omtrx[i][j]);
                }
            }
            if (maxb == 0.0) {
                /* Singular matrix */
                System.out.println("Singular matrix encountered");
                throw new IllegalArgumentException("Singular matrix");
            }
            wk[i] = 1.0 / maxb;
        }
        for (j = 0; j < size; j++) {
            for (i = 0; i < j; i++) {
                sum = omtrx[i][j];
                for (k = 0; k < i; k++) {
                    sum -= omtrx[i][k] * omtrx[k][j];
                }
                omtrx[i][j] = sum;
            }
            maxb = 0.0;
            for (i = j; i < size; i++) {
                sum = omtrx[i][j];
                for (k = 0; k < j; k++) {
                    sum -= omtrx[i][k] * omtrx[k][j];
                }
                omtrx[i][j] = sum;
                tmp = wk[i] * Math.abs(sum);
                if (tmp >= maxb) {
                    maxb = tmp;
                    maxi = i;
                }
            }
            if (j != maxi) {
                for (k = 0; k < size; k++) {
                    tmp = omtrx[maxi][k];
                    omtrx[maxi][k] = omtrx[j][k];
                    omtrx[j][k] = tmp;
                }
                aw = -aw;
                wk[maxi] = wk[j];
            }
            index[j] = maxi;
            if (omtrx[j][j] == 0.0) {
                omtrx[j][j] = MachineAccuracy.EPSILON;
            }
            if (j != size - 1) {
                tmp = 1.0 / omtrx[j][j];
                for (i = j + 1; i < size; i++) {
                    omtrx[i][j] *= tmp;
                }
            }
        }
        for (jx = 0; jx < size; jx++) {
            for (ix = 0; ix < size; ix++) {
                wk[ix] = 0.0;
            }
            wk[jx] = 1.0;
            l = -1;
            for (i = 0; i < size; i++) {
                idx = index[i];
                sum = wk[idx];
                wk[idx] = wk[i];
                if (l != -1) {
                    for (j = l; j < i; j++) {
                        sum -= omtrx[i][j] * wk[j];
                    }
                } else if (sum != 0.0) {
                    l = i;
                }
                wk[i] = sum;
            }
            for (i = size - 1; i >= 0; i--) {
                sum = wk[i];
                for (j = i + 1; j < size; j++) {
                    sum -= omtrx[i][j] * wk[j];
                }
                wk[i] = sum / omtrx[i][i];
            }
            for (ix = 0; ix < size; ix++) {
                imtrx[ix][jx] = wk[ix];
            }
        }
        wk = null;
        index = null;
        omtrx = null;
    }

    public void printDetails(){
            System.out.println("modelChoose: "+modelChoose.getValue());
            System.out.println("logKappa: "+logKappa.getValue());
            System.out.println("logTN: "+logTN.getValue());
            System.out.println("logAC: "+logAC.getValue());
            System.out.println("logAT: "+logAT.getValue());
            System.out.println("logGC: "+logGC.getValue());
            System.out.println("frequences: "+
                    frequencies.getValue(0)+ " "+
                    frequencies.getValue(1)+ " "+
                    frequencies.getValue(2)+ " "+
                    frequencies.getValue(3)+ " "
            );

        }

    protected void makeAccept() {
    	super.accept();
    }
}
