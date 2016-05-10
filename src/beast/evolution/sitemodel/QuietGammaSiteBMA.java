package beast.evolution.sitemodel;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.QuietRealParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.Node;
import org.apache.commons.math.distribution.GammaDistribution;
import org.apache.commons.math.distribution.GammaDistributionImpl;

/**
 * @author Chieh-Hsi Wu
 */
@Description("This class facilitates Bayesian variable selection for parameters in the Gamma site model.")
public class QuietGammaSiteBMA extends QuietSiteModel {
    public Input<QuietRealParameter> modelChoiceInput = new Input<QuietRealParameter>(
            "modelChoice",
            "And indicator that represents the current Gamma site model",
            Input.Validate.REQUIRED
    );

    public Input<Boolean> invPrLogitInput = new Input<Boolean>(
            "invPrLogit",
            "Whether invPr has logit transformation.",
            Input.Validate.REQUIRED
    );

    private boolean invPrLogit;
    private QuietRealParameter modelChoice;
    public void initAndValidate() {
        this.modelChoice = modelChoiceInput.get();
        this.invPrLogit = invPrLogitInput.get();
        gammaCatCount = gammaCategoryCount.get();
        categoryCount = gammaCategoryCount.get();
        super.initAndValidate();

        addCondition(modelChoice);
    }

    private IntegerParameter indicator;

    public static final int SHAPE_INDEX = 0;
    public static final int INVAR_INDEX = 1;
    public static final int PRESENT = 1;
    public static final int ABSENT = 0;
    public static final int[][] INDICATORS = {
            {ABSENT, ABSENT,},
            {ABSENT, PRESENT},
            {PRESENT, ABSENT},
            {PRESENT, PRESENT}
        };


    public QuietGammaSiteBMA(){

    }



    public QuietGammaSiteBMA(
            SubstitutionModel substModel,
            QuietRealParameter muParameter,
            QuietRealParameter modelChoice,
            boolean invPr)throws Exception{

        this(substModel,
                muParameter,
                null,
                null,
                false,
                1,
                modelChoice,
                invPr);


    }

    public QuietGammaSiteBMA(SubstitutionModel substModel,
                             QuietRealParameter muParameter,
                             RealParameter shapeParameter,
                             RealParameter invarParameter,
                             boolean useBeast1StyleGamma,
                             int gammaCategoryCount,
                             QuietRealParameter modelChoice,
                             boolean invPrLogit) throws Exception{
        //m_bPropInvariantIsCategory = false;


        this.modelChoice = modelChoice;
        substitutionModel = (SubstitutionModel.Base)substModel;

        this.useBeast1StyleGamma = useBeast1StyleGamma;

        if (muParameter == null) {
            this.muParameter = new RealParameter("1.0");
        }else{
            this.muParameter = muParameter;
        }
        //System.out.println(muParameter.getLower());
        muParameter.setBounds(Math.max(muParameter.getLower(), 0.0), Math.min(muParameter.getUpper(), Double.POSITIVE_INFINITY));


        this.shapeParameter = shapeParameter;
        if (shapeParameter != null) {
            // The quantile calculator fails when the shape parameter goes much below
            // 1E-3 so we have put a hard lower bound on it. If this is not there then
            // the category rates can go to 0 and cause a -Inf likelihood (whilst this
            // is not a problem as the state will be rejected, it could mask other issues
            // and this seems the better approach.
            shapeParameter.setBounds(Math.max(shapeParameter.getLower(), 1.0E-3), Math.min(shapeParameter.getUpper(), 1.0E3));
        }

        this.invPrLogit = invPrLogit;
        if (invarParameter == null) {

            this.invarParameter = new RealParameter("0.0");
            this.invarParameter.setBounds(Math.max(0.0, this.invarParameter.getLower()), Math.min(1.0, this.invarParameter.getUpper()));

        }else{
            this.invarParameter = invarParameter;

        }



        if (/*invarParameter != null && */(this.invarParameter.getValue() < 0 || this.invarParameter.getValue() > 1) && !invPrLogit) {
            throw new Exception("proportion invariant should be between 0 and 1: "+this.invarParameter.getValue());
        }

        gammaCatCount = gammaCategoryCount;
        categoryCount = gammaCategoryCount;
        //System.out.println("gammaCategoryCount: "+gammaCategoryCount);
        //System.out.println("categoryCount: "+categoryCount);
        //System.out.println("initiate");
        refresh();

        addCondition(muParameter);
        addCondition(invarParameter);
        addCondition(shapeParameter);
        addCondition(modelChoice);


    }



    @Override
	protected void refresh() {
        if (shapeParameter != null) {
            //categoryCount = gammaCategoryCount.get();
            if (categoryCount < 1) {
            	//throw new RuntimeException("SiteModel: Invalid category count (" + categoryCount + ") Setting category count to 1");
               	categoryCount = 1;
            }

            // The quantile calculator fails when the shape parameter goes much below
            // 1E-3 so we have put a hard lower bound on it. If this is not there then
            // the category rates can go to 0 and cause a -Inf likelihood (whilst this
            // is not a problem as the state will be rejected, it could mask other issues
            // and this seems the better approach.
            shapeParameter.setBounds(0.0, Double.POSITIVE_INFINITY);
        } else {
            categoryCount = 1;
        }
        //System.out.println("m_bPropInvariantIsCategory: "+m_bPropInvariantIsCategory);
        //if (invarParameter.getValue() > 0) {
        	if (hasPropInvariantCategory) {
        		categoryCount += 1;
        	}
            //invarParameter.setBounds(0.0, 1.0);
        //}

        categoryRates = new double[categoryCount];
        categoryProportions = new double[categoryCount];
        calculateCategoryRates(null);
        //ratesKnown = false;
	}

    public void setRatesKnown(boolean known){
        ratesKnown = known;
    }


    protected boolean requiresRecalculation() {
        /*if(modelChoice.somethingIsDirty()){
            ratesKnown = false;
            return true;
        }
        return super.requiresRecalculation();*/
       // we only get here if something is dirty in its inputs
        boolean recalculate = false;
        if(substitutionModel.isDirtyCalculation()){
            recalculate = true;

        }else if(modelChoice.somethingIsDirty()){


            recalculate = true;
        }else if(shapeParameter.somethingIsDirty()){
            //if(INDICATORS[getCurrModel()][SHAPE_INDEX] == PRESENT){
                recalculate = true;
            //}
        }else if(invarParameter.somethingIsDirty()){
            //if(INDICATORS[getCurrModel()][INVAR_INDEX] == PRESENT){
                recalculate = true;
            //}
        }

        if(recalculate){
            ratesKnown = false;
        }
        //return recalculate;

        //System.out.println("recalculate: "+recalculate);
        return recalculate||muParameter.somethingIsDirty();


    }



    /*protected boolean requiresRecalculation() {
        // do explicit check whether any of the non-substitution model parameters changed
        if (categoryCount > 1) {
            if (shapeParameter != null && shapeParameter.somethingIsDirty() ||
                    muParameter.somethingIsDirty() ||
                    invarParameter.somethingIsDirty()) {
                ratesKnown = false;
                return true;
            }

        } else if (muParameter.somethingIsDirty() || !m_bPropInvariantIsCategory && invarParameter.somethingIsDirty()) {
                ratesKnown = false;

            return true;
        }

//    	ratesKnown = false;
        // we only get here if something is dirty in its inputs, so always return true
        return substitutionModel.isDirtyCalculation();
    }    */

    private int getCurrModel(){
        return (int)(double)modelChoice.getValue();

    }


 /**
     * discretization of gamma distribution with equal proportions in each
     * category
     */
    @Override
    protected void calculateCategoryRates(Node node) {

        double propVariable = 1.0;
        int cat = 0;

        //if (/*invarParameter != null && */invarParameter.getValue() > 0 ) {
            //System.out.println("----------------");
            double pr;
            if(invPrLogit){
                pr = 1.0/(1.0 + Math.exp(-invarParameter.getValue()));
                //System.out.println(pr);
            }else{
                pr = invarParameter.getValue();
            }

            if (hasPropInvariantCategory) {
                categoryRates[0] = 0.0;
                //System.out.println(getCurrModel()+" "+INVAR_INDEX);
                categoryProportions[0] = pr*INDICATORS[getCurrModel()][INVAR_INDEX];
                //System.out.println(invarParameter.getValue()+" "+INDICATORS[getCurrModel()][INVAR_INDEX]);
                //System.out.println("categoryProportions[0]: "+categoryProportions[0]);
            }

            //System.out.println(invarParameter.getID()+" " +invarParameter.getValue()+" "+pr);
            propVariable = 1.0 - pr*INDICATORS[getCurrModel()][INVAR_INDEX];
            if (hasPropInvariantCategory) {
                cat = 1;
            }
        //}

        //System.out.println("categoryProportions[0]: "+categoryProportions[0]);

        if (INDICATORS[getCurrModel()][SHAPE_INDEX] == PRESENT) {

            final double a = shapeParameter.getValue();
            double mean = 0.0;
            final int gammaCatCount = categoryCount - cat;
            //System.out.println("a: "+a);
            final GammaDistribution g = new GammaDistributionImpl(a, 1.0 / a);
            //System.out.println("gammaCatCount:"+gammaCatCount);
            //if(gammaCatCount == 3){
                //throw new RuntimeException("");
            //}
            for (int i = 0; i < gammaCatCount; i++) {
                try {
                    // RRB: alternative implementation that seems equally good in
                    // the first 5 significant digits, but uses a standard distribution object
                    if(a < 1e-3){

                        categoryRates[i + cat] = Double.NEGATIVE_INFINITY;
                    }else if(a > 1e10){
                        categoryRates[i + cat] = 1.0;
                    }else if (useBeast1StyleGamma) {
                        categoryRates[i + cat] = GammaDistributionQuantile((2.0 * i + 1.0) / (2.0 * gammaCatCount), a, 1.0 / a);
                	} else {
                		categoryRates[i + cat] = g.inverseCumulativeProbability((2.0 * i + 1.0) / (2.0 * gammaCatCount));
                	}

                } catch (Exception e) {
                    e.printStackTrace();
                    System.err.println("Something went wrong with the gamma distribution calculation");
                    System.exit(-1);
                }
                mean += categoryRates[i + cat];

                categoryProportions[i + cat] = propVariable / gammaCatCount;
            }

            if(a >= 1e-3 ){

                mean = (propVariable * mean) / gammaCatCount;

                for (int i = 0; i < gammaCatCount; i++) {

                    categoryRates[i + cat] /= mean;

                }
            }
        } else {

            int gammaCatCount = categoryCount - cat;
            //System.out.println("Hi!");
            for(int i = cat; i < categoryRates.length;i++){
                categoryRates[i] = 1.0 / propVariable/gammaCatCount;

                categoryProportions[i] = propVariable/gammaCatCount;
            }
        }
        /*System.out.println("-------------------------------");
        System.out.println("ID: "+getID());
        System.out.println(modelChoice);
        System.out.print (invarParameter.getValue()*INDICATORS[getCurrModel()][INVAR_INDEX]+" ");
        System.out.println(getID());
        System.out.println("alpha: "+shapeParameter.getValue());
        System.out.println("invPr: "+invarParameter.getValue());
        System.out.println("siteModel: "+modelChoice.getValue());
        System.out.println("rate: "+muParameter.getValue());
        for(int i = 0; i < categoryRates.length;i++){
            System.out.print(categoryRates[i]+" ");
        }
        System.out.println();

        System.out.println(invarParameter.getValue());
        for(int i = 0; i < categoryProportions.length;i++){

            System.out.print(categoryProportions[i]+" ");
        }
        System.out.println();
         System.out.println("----------------"); */


        ratesKnown = true;
    }

    @Override
    public double getProportionInvariant() {
        //if (invarParameter == null) {
        //	return 0;
        //}
        //System.out.println(getCurrModel()+" "+INDICATORS[getCurrModel()]+" "+INDICATORS[getCurrModel()][INVAR_INDEX]);
        double pr = invarParameter.getValue();
        /*if(Math.abs(pr - 6.383340011058765) < 1e-15){
            System.out.println(pr*INDICATORS[getCurrModel()][INVAR_INDEX]+" "+invarParameter.getValue());

        }*/
        if(invPrLogit){
            pr = 1.0/(1.0+Math.exp(-pr));
        }
        //System.out.println("pr: "+pr);
        return pr*INDICATORS[getCurrModel()][INVAR_INDEX];
    }


    public void setMuValueQuietly(double muValue){
        ((QuietRealParameter)muParameter).setValueQuietly(0,muValue);
    }

    public void setShapeValueQuietly(double shapeValue){
        ((QuietRealParameter)shapeParameter).setValueQuietly(0,shapeValue);
    }

    public void setInvPrValueQuietly(double invPrValue){
        ((QuietRealParameter)invarParameter).setValueQuietly(0,invPrValue);
    }


    public void setModelChoiceQuietly(double modelChoiceVal){
        modelChoice.setValueQuietly(0,modelChoiceVal);
    }

    public double getMuValue(){
        return muParameter.getValue();
    }

    public double getAlphaValue(){
        return shapeParameter.getValue();
    }

    public double getInvPrValue(){
        return invarParameter.getValue();
    }

    public double getModelChoiceValue(){
        return modelChoice.getValue();
    }

    public int getRateID(){
        return ((QuietRealParameter)muParameter).getIDNumber();
    }





}
