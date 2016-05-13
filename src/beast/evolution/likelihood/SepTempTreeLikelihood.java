package beast.evolution.likelihood;

import beast.core.Description;
import beast.evolution.alignment.Alignment;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.sitemodel.DPNtdRateSepSiteModel;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.substitutionmodel.NtdBMA;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Tree;

/**
 * @author Chieh-Hsi Wu
 */
@Description("The tree likelihood used by the operators when using independent DPPs to the substitution model and rate partitioning.")
public class SepTempTreeLikelihood extends TempTreeLikelihood{
    public Input<DPNtdRateSepSiteModel> dpNtdRateSepSiteModelInput = new Input<DPNtdRateSepSiteModel>(
            "NtdRateSepSiteModel",
            "A site model that records separate clustering of the substitution model and rates",
            Input.Validate.REQUIRED
    );

    private DPNtdRateSepSiteModel dpNtdRateSepSiteModel;

    public SepTempTreeLikelihood(){}

    public SepTempTreeLikelihood(Alignment alignment,
                              Tree tree,
                              boolean useAmbiguities,
                              SiteModel siteModel,
                              BranchRateModel.Base branchRateModel,
                              DPNtdRateSepSiteModel dpNtdRateSepSiteModel){
        super(alignment, tree, useAmbiguities, siteModel, branchRateModel);
        this.dpNtdRateSepSiteModel = dpNtdRateSepSiteModel;

    }
    public void initAndValidate() {
        dpNtdRateSepSiteModel = dpNtdRateSepSiteModelInput.get();
        super.initAndValidate();
    }

    public double calculateLogP(
            RealParameter rateParameter,
            int siteIndex){

        try{


            //Retrive the substModel of the given site
            SubstitutionModel substModel = dpNtdRateSepSiteModel.getModel(siteIndex);
            //Set the substitution model of the siteModel to the substitution model retrieved
            siteModelInput.get().substModelInput.setValue(substModel,siteModelInput.get());
            //Set the rate to the specified rate
            siteModelInput.get().getRateParameter().setValueQuietly(0,rateParameter.getValue());
            //Retrieve the pattern
            int iPat = alignment.getPatternIndex(siteIndex);
            //System.out.println("recompute:");
            logP = treeLiks[iPat].calculateLogP();
        }catch(Exception e){
            throw new RuntimeException(e);

        }
        return logP;

    }

    public double calculateLogP(
            RealParameter modelParameters,
            RealParameter modelCode,
            RealParameter freqs,
            int siteIndex){
        try{

            //Retrive the rate of the given site
            RealParameter muParameter = dpNtdRateSepSiteModel.getRate(siteIndex);
            //Set the mutational of the siteModel to the retrieved rate
            siteModelInput.get().getRateParameter().setValueQuietly(0,muParameter.getValue());

            //Set the subsitution model parameters to the specified values.
            if(substModel instanceof NtdBMA){
                ((NtdBMA)substModel).getLogKappa().setValueQuietly(0,modelParameters.getValue(0));
                ((NtdBMA)substModel).getLogTN().setValueQuietly(0,modelParameters.getValue(1));
                ((NtdBMA)substModel).getLogAC().setValueQuietly(0,modelParameters.getValue(2));
                ((NtdBMA)substModel).getLogAT().setValueQuietly(0,modelParameters.getValue(3));
                ((NtdBMA)substModel).getLogGC().setValueQuietly(0,modelParameters.getValue(4));
                ((NtdBMA)substModel).getModelChoose().setValueQuietly(0,modelCode.getValue());
                ((NtdBMA)substModel).getFreqs().setValueQuietly(0,freqs.getValue(0));
                ((NtdBMA)substModel).getFreqs().setValueQuietly(1,freqs.getValue(1));
                ((NtdBMA)substModel).getFreqs().setValueQuietly(2,freqs.getValue(2));
                ((NtdBMA)substModel).getFreqs().setValueQuietly(3,freqs.getValue(3));

            }else{
                throw new RuntimeException("Need NtdBMA");
            }
            ((NtdBMA)substModel).setUpdateMatrix(true);
            siteModelInput.get().substModelInput.setValue(substModel,siteModelInput.get());

            int iPat = alignment.getPatternIndex(siteIndex);
            logP = treeLiks[iPat].calculateLogP();
        }catch(Exception e){
            throw new RuntimeException(e);

        }
        return logP;
    }



}
