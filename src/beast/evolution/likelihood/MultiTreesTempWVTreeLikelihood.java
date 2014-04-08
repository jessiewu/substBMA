package beast.evolution.likelihood;


import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.Alignment;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.sitemodel.DummySiteModel;
import beast.evolution.tree.Tree;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * @author Chieh-Hsi Wu
 */
@Description("Weight variable tree likelihood that can handle multiple trees.")
public class MultiTreesTempWVTreeLikelihood extends Distribution {
    public Input<List<Alignment>> alignmentsInput = new Input<List<Alignment>>(
            "alignment",
            "A list of alignments used as data.",
            new ArrayList<Alignment>(),
            Input.Validate.REQUIRED
    );

    public Input<List<Tree>> treesInput = new Input<List<Tree>>(
            "tree",
            "A list of trees representing the underlying phylogenies that generated the data.",
            new ArrayList<Tree>(),
            Input.Validate.REQUIRED
    );
    public Input<Boolean> useAmbiguitiesInput = new Input<Boolean>(
            "useAmbiguities",
            "flag to indicate leafs that sites containing ambigue states should be handled instead of ignored (the default)",
            false
    );

    public Input<BranchRateModel.Base> branchRateModelInput = new Input<BranchRateModel.Base>(
            "branchRateModel",
            "A model describing the rates on the branches of the beast.tree."
    );

    public Input<List<DummySiteModel>> siteModelsInput = new Input<List<DummySiteModel>>(
            "siteModel",
            "Models the evolution of a site in an alignment",
            Input.Validate.REQUIRED
    );

    public Input<DPMultiTreesTreeLikelihood> dpMultiTreesTreeLikelihoodInput = new Input<DPMultiTreesTreeLikelihood>(
            "dpMultiTreesTreeLikelihood",
            "The tree likelihood that facilitates nucleotide substitution model partition by site with multiple alignment/trees.",
            Input.Validate.REQUIRED
    );


    private List<Alignment> alignments;
    private List<Tree> trees;
    private TempWVTreeLikelihood[] tempTreeLikelihoods;
    private DPMultiTreesTreeLikelihood dpMultiTreesTreeLikelihood;
    private int[] alignmentStartingIndex;
    public void initAndValidate(){
        alignments = alignmentsInput.get();
        trees = treesInput.get();
        tempTreeLikelihoods = new TempWVTreeLikelihood[alignments.size()];
        for(int i = 0; i < tempTreeLikelihoods.length; i++){
            tempTreeLikelihoods[i] = new TempWVTreeLikelihood(
                    alignments.get(i).getWeights(),
                    alignments.get(i),
                    trees.get(i),
                    useAmbiguitiesInput.get(),
                    siteModelsInput.get().get(i),
                    branchRateModelInput.get()

            );
        }
        dpMultiTreesTreeLikelihood = dpMultiTreesTreeLikelihoodInput.get();
        alignmentStartingIndex = new int[alignments.size()];
        alignmentStartingIndex[0] = 0;
        for(int i = 1; i < alignmentStartingIndex.length; i++){
            alignmentStartingIndex[i] = alignments.get(i - 1).getSiteCount()+alignmentStartingIndex[0];
        }

    }


    public double[] calculateLogP (
            RealParameter modelParameters,
            RealParameter modelCode,
            RealParameter freqs,
            RealParameter rate,
            int[] sites) throws Exception{
        double[] siteLogP = new double[sites.length];


        //Sort out the weights
        int[] alignmentWeights = new int[alignments.size()];
        int[] siteAlignmentIndex = new int[alignments.size()];
        int[][] multiAlignPatternWeights = setupWeights(sites, alignmentWeights, siteAlignmentIndex);

        for(int i = 0; i <alignmentWeights.length; i++){
            if(alignmentWeights[i] > 0){
                tempTreeLikelihoods[i].setPatternWeights(multiAlignPatternWeights[i]);
                tempTreeLikelihoods[i].calculateLogP(
                        modelParameters,
                        modelCode,
                        freqs,
                        rate
                );

            }
        }

        for(int i = 0; i < sites.length;i++){
            siteLogP[i] = tempTreeLikelihoods[siteAlignmentIndex[i]].getCurrSiteLikelihood(sites[i] - alignmentStartingIndex[siteAlignmentIndex[i]]);
        }
        return siteLogP;
    }


    public double[] calculateLogP (
            RealParameter modelParameters,
            RealParameter modelCode,
            RealParameter freqs,
            int[] sites) throws Exception{
        double[] siteLogP = new double[sites.length];


        //Sort out the weights
        int[] alignmentWeights = new int[alignments.size()];
        int[] siteAlignmentIndex = new int[alignments.size()];
        int[][] multiAlignPatternWeights = setupWeights(sites, alignmentWeights, siteAlignmentIndex);

        for(int i = 0; i <alignmentWeights.length; i++){
            if(alignmentWeights[i] > 0){
                tempTreeLikelihoods[i].setPatternWeights(multiAlignPatternWeights[i]);
                tempTreeLikelihoods[i].calculateLogP(
                        modelParameters,
                        modelCode,
                        freqs
                );

            }
        }

        for(int i = 0; i < sites.length;i++){
            siteLogP[i] = tempTreeLikelihoods[siteAlignmentIndex[i]].getCurrSiteLikelihood(sites[i] - alignmentStartingIndex[siteAlignmentIndex[i]]);
        }
        return siteLogP;
    }

    public double[] calculateLogP (
            RealParameter rate,
            int[] sites) throws Exception{
        double[] siteLogP = new double[sites.length];


        //Sort out the weights
        int[] alignmentWeights = new int[alignments.size()];
        int[] siteAlignmentIndex = new int[alignments.size()];
        int[][] multiAlignPatternWeights = setupWeights(sites, alignmentWeights, siteAlignmentIndex);

        for(int i = 0; i <alignmentWeights.length; i++){
            if(alignmentWeights[i] > 0){
                tempTreeLikelihoods[i].setPatternWeights(multiAlignPatternWeights[i]);
                tempTreeLikelihoods[i].calculateLogP(rate);
            }
        }

        for(int i = 0; i < sites.length;i++){
            siteLogP[i] = tempTreeLikelihoods[siteAlignmentIndex[i]].getCurrSiteLikelihood(sites[i] - alignmentStartingIndex[siteAlignmentIndex[i]]);
        }
        return siteLogP;
    }

    private int[][] setupWeights(int[] sites, int[] alignmentWeights, int[] siteAlignmentIndex){
        int[][] multiAlignPatternWeights = new int[alignments.size()][];
        for(int i = 0; i < multiAlignPatternWeights.length; i++){
            multiAlignPatternWeights[i] = new int[alignments.get(i).getPatternCount()];
        }

        for(int i = 0; i < sites.length; i++){
            siteAlignmentIndex[i] = dpMultiTreesTreeLikelihood.getAlignmentIndexBySite(sites[i]);
            multiAlignPatternWeights[siteAlignmentIndex[i]][sites[i] - alignmentStartingIndex[siteAlignmentIndex[i]]]++;
            alignmentWeights[siteAlignmentIndex[i]]++;
        }
        return multiAlignPatternWeights;

    }

    public List<String> getConditions(){
        return null;
    }

    public List<String> getArguments(){
        return null;
    }

    public void sample(State state, Random random){
        throw new RuntimeException("Not yet implemented as it doesn't make much sense to do so in this case");
    }
}
