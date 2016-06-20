package beast.evolution.likelihood;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.alignment.AlignmentSubset;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.sitemodel.DPMultiAlignSiteModel;
import beast.evolution.sitemodel.DummySiteModel;
import beast.evolution.sitemodel.QuietSiteModel;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.NtdBMA;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Chieh-Hsi Wu
 */
@Description("Temporary likelihood that can take up multiple alignments.")
public class MultiAlignmentTempTreeLikelihood extends TempTreeLikelihood {
    public Input<List<Alignment>> alignmentsInput = new Input<List<Alignment>>(
            "data",
            "sequence data for the beast.tree",
            new ArrayList<Alignment>(),
            Input.Validate.REQUIRED
    );

    public Input<List<DummySiteModel>> siteModelsInput = new Input<List<DummySiteModel>>(
            "siteModel",
            "Site model that has information on the substitution process.",
            new ArrayList<DummySiteModel>(),
            Input.Validate.REQUIRED
    );

    public Input<DPMultiAlignSiteModel> dpMultiAlignSiteModelInput = new Input<DPMultiAlignSiteModel>(
            "dpMultiAlignSiteModel",
            "The site model that accommodates multiple alignments.",
            Input.Validate.REQUIRED
    );

    public MultiAlignmentTempTreeLikelihood(){
        dataInput.setRule(Input.Validate.OPTIONAL);
        siteModelInput.setRule(Input.Validate.OPTIONAL);
    }

    List<DummySiteModel> siteModels;
    List<Alignment> alignments;
    DPMultiAlignSiteModel dpMultiAlignSiteModel;
    int[] prevPatternCount;
    int[] siteIndexWithinAlignment;
    public void initAndValidate() {
        defaultMu = new RealParameter(new Double[]{1.0});
        alignments = alignmentsInput.get();
        siteModels = siteModelsInput.get();
        if(alignments.size() != siteModels.size()){
            throw new RuntimeException("The number of alignments("+alignments.size()+
                    ") provided != the number of site models ("+siteModels.size()+").");
        }

        //int siteCount = 0;
        int patternCount = 0;
        prevPatternCount = new int[alignments.size()];
        for(int i = 1; i < prevPatternCount.length; i++){
            prevPatternCount[i] = prevPatternCount[i - 1]+ alignments.get(i - 1).getPatternCount();
        }
        patternCount = prevPatternCount[prevPatternCount.length - 1]+ alignments.get(alignments.size() - 1).getPatternCount();


        treeLiks = new TempSiteTreeLikelihood[patternCount];


        //
        int[] firstPatternOccur = new int[patternCount];
        for(int iPat = 0; iPat < firstPatternOccur.length; iPat++){
            firstPatternOccur[iPat] = -1;
        }

        int prevPatternCount = 0;
        int[] patternAlignmentMap = new int[patternCount];
        for(int iAlign = 0; iAlign < alignments.size(); iAlign++){
            Alignment alignment = alignments.get(iAlign);

            for(int iSite = 0; iSite < alignment.getSiteCount(); iSite++){
                int iPat = alignment.getPatternIndex(iSite);
                if(firstPatternOccur[prevPatternCount+iPat] == -1){
                    firstPatternOccur[prevPatternCount+iPat] = iSite;
                    patternAlignmentMap[prevPatternCount+iPat] = iAlign;
                }

            }

            prevPatternCount += alignment.getPatternCount();
        }

        //for(int iAlign = 0; iAlign < alignments.size(); iAlign++){
            //System.out.println(iAlign+" patternCount: "+alignments.get(iAlign).getPatternCount());

            BranchRateModel.Base brModel = branchRateModelInput.get() != null
                ? branchRateModelInput.get()
                : trueLikelihoodInput.get().branchRateModelInput.get();

            for(int i = 0; i < firstPatternOccur.length; i++){
                //System.out.println("i: "+i+" "+firstPatternOccur[i]);
                AlignmentSubset sitePattern = new AlignmentSubset(alignments.get(patternAlignmentMap[i]),firstPatternOccur[i]);
                TempSiteTreeLikelihood treeLik = new TempSiteTreeLikelihood(
                        sitePattern,
                        treeInput.get(),
                        useAmbiguitiesInput.get(),
                        siteModels.get(patternAlignmentMap[i]),
                        brModel
                );
                treeLiks[i] = treeLik;
            }
        //}
        dpMultiAlignSiteModel = dpMultiAlignSiteModelInput.get();



        int[] prevAlignEndIndex = new int[alignments.size()];
        prevAlignEndIndex[0] = 0;
        for(int i = 1; i < prevAlignEndIndex.length; i++){
            prevAlignEndIndex[i] = prevAlignEndIndex[i - 1] + alignments.get(i - 1).getSiteCount();
        }
        int siteCount = prevAlignEndIndex[prevAlignEndIndex.length - 1] +
                alignments.get(alignments.size() - 1).getSiteCount();

        siteIndexWithinAlignment = new int[siteCount];
        for(int i = 0; i < siteIndexWithinAlignment.length; i++){
            siteIndexWithinAlignment[i] = i - prevAlignEndIndex[dpMultiAlignSiteModel.getAlignmentIndex(i)];

        }

    }

    public double calculateLogP(
            RealParameter modelParameters,
            RealParameter modelCode,
            RealParameter freqs,
            int site){
        throw new RuntimeException("No proper implementation yet");
    }

    public double calculateLogP(
            RealParameter modelParameters,
            RealParameter modelCode,
            RealParameter freqs,
            RealParameter rate,
            int site){
        throw new RuntimeException("No proper implementation yet");
    }

    public double calculateLogP(
            RealParameter rateParameter,
            int site){
        int alignmentIndex = dpMultiAlignSiteModel.getAlignmentIndex(site);
        try{
            //System.out.println(alignmentIndex);
            siteModels.get(alignmentIndex).getRateParameter().setValueQuietly(0,rateParameter.getValue());
            //siteModels.get(alignmentIndex).getRateParameter().setValueQuietly(0,5.8197453245394814E-5);

            int iPat = alignments.get(alignmentIndex).getPatternIndex(siteIndexWithinAlignment[site]);
            //System.out.println("site: "+site+" "+iPat);
            //System.out.println("site: "+site+" "+alignmentIndex);
            logP = treeLiks[iPat+ prevPatternCount[alignmentIndex]].calculateLogP();
        }catch(Exception e){
            throw new RuntimeException(e);

        }
        return logP;
    }

}
