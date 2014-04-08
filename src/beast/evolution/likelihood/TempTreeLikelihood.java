package beast.evolution.likelihood;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.AlignmentSubset;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.sitemodel.QuietSiteModel;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.NtdBMA;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.Tree;

import java.util.List;
import java.util.Random;

/**
 * @author Chieh-Hsi Wu
 */
@Description("Used to compute the tree likelihood for operators that samples the partitioning space. Used when there is only one partitioning.")
public class TempTreeLikelihood extends Distribution {
    public Input<Alignment> dataInput = new Input<Alignment>(
            "data",
            "sequence data for the beast.tree",
            Input.Validate.REQUIRED
    );

    public Input<QuietSiteModel> siteModelInput = new Input<QuietSiteModel>(
            "siteModel",
            "Models the evolution of a site in an alignment",
            Input.Validate.REQUIRED
    );

    public Input<Tree> treeInput = new Input<Tree>(
            "tree",
            "phylogenetic beast.tree with sequence data in the leafs",
            Input.Validate.REQUIRED
    );

    public Input<BranchRateModel.Base> branchRateModelInput = new Input<BranchRateModel.Base>(
            "branchRateModel",
            "A model describing the rates on the branches of the beast.tree."
    );

    public Input<Boolean> useAmbiguitiesInput = new Input<Boolean>(
            "useAmbiguities",
            "flag to indicate leafs that sites containing ambigue states should be handled instead of ignored (the default)",
            false
    );

    protected RealParameter defaultMu;
    protected Alignment alignment;
    protected TempSiteTreeLikelihood[] treeLiks;
    protected SubstitutionModel substModel;

    public TempTreeLikelihood(){}

    public TempTreeLikelihood(Alignment alignment,
                              Tree tree,
                              boolean useAmbiguities,
                              SiteModel siteModel,
                              BranchRateModel.Base branchRateModel){
        int siteCount = alignment.getSiteCount();
        int patternCount = alignment.getPatternCount();
        treeLiks = new TempSiteTreeLikelihood[patternCount];

        int[] firstPatternOccur = new int[patternCount];
        for(int iPat = 0; iPat < firstPatternOccur.length; iPat++){
            firstPatternOccur[iPat] = -1;
        }
        for(int iSite = 0; iSite < siteCount; iSite++){

            int iPat = alignment.getPatternIndex(iSite);
            if(firstPatternOccur[iPat] == -1){
                firstPatternOccur[iPat] = iSite;
            }
        }
        for(int i = 0; i < firstPatternOccur.length; i++){
            AlignmentSubset sitePattern = new AlignmentSubset(alignment,firstPatternOccur[i]);
            TempSiteTreeLikelihood treeLik = new TempSiteTreeLikelihood(
                    sitePattern,
                    tree,
                    useAmbiguities,
                    siteModel,
                    branchRateModel

            );

            treeLiks[i] = treeLik;
        }

        substModel = siteModelInput.get().getSubstitutionModel();

    }
    
    public void initAndValidate() throws Exception{
        defaultMu = new RealParameter(new Double[]{1.0});
        alignment = dataInput.get();
        int siteCount = alignment.getSiteCount();
        int patternCount = alignment.getPatternCount();
        treeLiks = new TempSiteTreeLikelihood[patternCount];

        //
        int[] firstPatternOccur = new int[patternCount];
        for(int iPat = 0; iPat < firstPatternOccur.length; iPat++){
            firstPatternOccur[iPat] = -1;
        }
        for(int iSite = 0; iSite < siteCount; iSite++){

            int iPat = alignment.getPatternIndex(iSite);
            if(firstPatternOccur[iPat] == -1){
                firstPatternOccur[iPat] = iSite;
            }
        }


        for(int i = 0; i < firstPatternOccur.length; i++){
            AlignmentSubset sitePattern = new AlignmentSubset(alignment,firstPatternOccur[i]);
            TempSiteTreeLikelihood treeLik = new TempSiteTreeLikelihood(
                    sitePattern,
                    treeInput.get(),
                    useAmbiguitiesInput.get(),
                    siteModelInput.get(),
                    branchRateModelInput.get()

            );

            treeLiks[i] = treeLik;
        }

        substModel = siteModelInput.get().getSubstitutionModel();
        

    }

    public void setSubstModelParameter(
            RealParameter modelParameters,
            RealParameter modelCode,
            RealParameter freqs){
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
    }

    public void setRateParameter(RealParameter rateParameter){
        siteModelInput.get().getRateParameter().setValueQuietly(0,rateParameter.getValue());
    }

    public double calculateLogP(int site){
        try{
            int iPat = alignment.getPatternIndex(site);
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
            int site){
        try{

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
            /*NtdBMA ntdBMA = DPNtdBMA.createNtdBMA(modelParameters,modelCode,freqs);

            SiteModel siteModel = new SiteModel();
            siteModel.initByName(
                    "substModel", ntdBMA
            );*/
            int iPat = alignment.getPatternIndex(site);
            //treeLiks[iPat].setSiteModel(siteModel);
            ((NtdBMA)substModel).setUpdateMatrix(true);
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
            RealParameter rate,
            int site){
        try{

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
            /*NtdBMA ntdBMA = DPNtdBMA.createNtdBMA(modelParameters,modelCode,freqs);

            SiteModel siteModel = new SiteModel();
            siteModel.initByName(
                    "substModel", ntdBMA
            );*/
            ((NtdBMA)substModel).setUpdateMatrix(true);
            siteModelInput.get().getRateParameter().setValueQuietly(0,rate.getValue());
            int iPat = alignment.getPatternIndex(site);         
            logP = treeLiks[iPat].calculateLogP();

        }catch(Exception e){
            throw new RuntimeException(e);

        }
        //System.out.println("site: "+site);
        //System.out.println("logP: "+logP);
        return logP;
    }


    public double calculateLogP(
            RealParameter rateParameter,
            int site){
        try{
            siteModelInput.get().getRateParameter().setValueQuietly(0,rateParameter.getValue());
            int iPat = alignment.getPatternIndex(site);
            logP = treeLiks[iPat].calculateLogP();
        }catch(Exception e){
            throw new RuntimeException(e);

        }
        return logP;
    }


 

    public List<String> getConditions(){
        return null;
    }

    public List<String> getArguments(){
        return null;

    }

    public boolean requiresRecalculation(){
        return false;
    }

    public void sample(State state, Random random){
        throw new RuntimeException("Not yet implemented as it doesn't make much sense to do so in this case");
    }





}
