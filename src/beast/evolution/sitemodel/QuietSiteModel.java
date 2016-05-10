package beast.evolution.sitemodel;

import beast.core.Description;
import beast.core.StateNode;
import beast.core.parameter.QuietRealParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.substitutionmodel.SwitchingNtdBMA;

import java.util.ArrayList;

/**
 * @author Chieh-Hsi Wu
 */
@Description("This site model is used to be created internally on the fly during the MCMC.")
public class QuietSiteModel extends SiteModel {
    @Override
    public void initAndValidate() {
        //System.out.println(getID()+": "+m_pSubstModel.get());
        substitutionModel = (SubstitutionModel.Base) substModelInput.get();

    	useBeast1StyleGamma = true; // useBeast1StyleGammaInput.get();
        muParameter = muParameterInput.get();
        if (muParameter == null) {
            muParameter = new RealParameter("1.0");
        }
        shapeParameter = shapeParameterInput.get();
        invarParameter = invarParameterInput.get();
        if (invarParameter == null) {
            invarParameter = new RealParameter("0.0");
            invarParameter.setBounds(Math.max(0.0, invarParameter.getLower()), Math.min(1.0, invarParameter.getUpper()));
        }

        //if (muParameter != null) {
        muParameter.setBounds(Math.max(muParameter.getLower(), 0.0), Math.min(muParameter.getUpper(), Double.POSITIVE_INFINITY));
        //}
        if (shapeParameter != null) {
            // The quantile calculator fails when the shape parameter goes much below
            // 1E-3 so we have put a hard lower bound on it. If this is not there then
            // the category rates can go to 0 and cause a -Inf likelihood (whilst this
            // is not a problem as the state will be rejected, it could mask other issues
            // and this seems the better approach.
            shapeParameter.setBounds(Math.max(shapeParameter.getLower(), 1.0E-3), Math.min(shapeParameter.getUpper(), 1.0E3));
        }


        if (/*invarParameter != null && */(invarParameter.getValue() < 0 || invarParameter.getValue() > 1)) {
            throw new RuntimeException("proportion invariant should be between 0 and 1");
        }
        refresh();

        addCondition(muParameterInput);
        addCondition(invarParameterInput);
        addCondition(shapeParameterInput);
    }


    public void makeAccept(){
        super.accept();
    }


    int gammaCatCount;
    SubstitutionModel.Base substitutionModel;

    public QuietSiteModel(){}

    public QuietSiteModel(SubstitutionModel substModel,
                          QuietRealParameter muParameter) {
        this(substModel,
                muParameter,
                null,
                null,
                false,
                1);

    }

    public QuietSiteModel(SubstitutionModel substModel,
                          QuietRealParameter muParameter,
                          RealParameter shapeParameter,
                          RealParameter invarParameter,
                          boolean useBeast1StyleGamma,
                          int gammaCategoryCount) {
        substitutionModel = (SubstitutionModel.Base)substModel;

        this.useBeast1StyleGamma = useBeast1StyleGamma;

        if (muParameter == null) {
            this.muParameter = new RealParameter("1.0");
        }else{
            this.muParameter = muParameter;
        }
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

        if (invarParameter == null) {

            this.invarParameter = new RealParameter("0.0");
            this.invarParameter.setBounds(Math.max(0.0, this.invarParameter.getLower()), Math.min(1.0, this.invarParameter.getUpper()));

        }else{
            this.invarParameter = invarParameter;

        }



        if (/*invarParameter != null && */(this.invarParameter.getValue() < 0 || this.invarParameter.getValue() > 1)) {
            throw new RuntimeException("proportion invariant should be between 0 and 1: "+this.invarParameter.getValue());
        }
        gammaCatCount = gammaCategoryCount;
        refresh();

        addCondition(muParameter);
        addCondition(invarParameter);
        addCondition(shapeParameter);

    }


    @Override
    protected void refresh() {
        if (shapeParameter != null) {
            categoryCount = gammaCatCount;
            if (categoryCount < 1) {
                System.out.println("SiteModel: Invalid category count (" + categoryCount + ") Setting category count to 1");
                categoryCount = 1;
            }
        } else {
            categoryCount = 1;
        }

        if (/*invarParameter != null && */invarParameter.getValue() > 0) {
            if (hasPropInvariantCategory) {
                categoryCount += 1;
            }
        }

        categoryRates = new double[categoryCount];
        categoryProportions = new double[categoryCount];
        calculateCategoryRates(null);
        //ratesKnown = false;
    }

    /**
     * add item to the list *
     * @param stateNode
     */
    public void addCondition(final StateNode stateNode) {
        if (stateNode == null) return;

        if (conditions == null) conditions = new ArrayList<String>();

        conditions.add(stateNode.getID());
    }

    @Override
    public SubstitutionModel.Base getSubstitutionModel() {
        return substitutionModel;
    }

    public QuietRealParameter getRateParameter(){
        return (QuietRealParameter)muParameter;
    }

    @Override
    protected boolean requiresRecalculation() {
        // do explicit check whether any of the non-substitution model parameters changed
        if (categoryCount > 1) {
            if (shapeParameter != null && shapeParameter.somethingIsDirty() ||
                    muParameter.somethingIsDirty() ||
                    invarParameter.somethingIsDirty()) {
                ratesKnown = false;
                return true;
            }

        } else if (muParameter.somethingIsDirty() || !hasPropInvariantCategory && invarParameter.somethingIsDirty()) {
                ratesKnown = false;

            return true;
        }

//    	ratesKnown = false;
        // we only get here if something is dirty in its inputs, so always return true
        return substitutionModel.isDirtyCalculation();
    }

    /*public boolean requiresRecalculation(){
        System.out.println(getRateParameter().getIDNumber()+ " "+((SwitchingNtdBMA)substitutionModel).getIDNumber());
        return super.requiresRecalculation();
    } */

    public void printDetails(){
        System.out.println(muParameter);
        System.out.println(shapeParameter);
        System.out.println(invarParameter);
    }
}
