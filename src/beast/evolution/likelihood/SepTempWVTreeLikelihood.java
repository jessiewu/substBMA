package beast.evolution.likelihood;

import beast.core.Description;
import beast.core.parameter.QuietRealParameter;
import beast.evolution.sitemodel.DPNtdRateSepSiteModel;
import beast.evolution.sitemodel.DummySiteModel;
import beast.evolution.substitutionmodel.SwitchingNtdBMA;
import beast.evolution.alignment.Alignment;
import beast.evolution.tree.Tree;
import beast.core.Input;
import beast.core.parameter.RealParameter;

import java.util.HashMap;
import java.util.ArrayList;

/**
 * @author Chieh-Hsi Wu
 */
@Description("The tree likelihood used by the operators when using independent DPPs to the substitution model and rate partitioning. It allows pattern weights to be modified.")
public class SepTempWVTreeLikelihood extends TempWVTreeLikelihood{
    public Input<DPNtdRateSepSiteModel> dpNtdRateSepSiteModelInput = new Input<DPNtdRateSepSiteModel>(
            "NtdRateSepSiteModel",
            "A site model that records separate clustering of the substitution model and rates",
            Input.Validate.REQUIRED
    );

    protected DPNtdRateSepSiteModel dpNtdRateSepSiteModel;
    protected SwitchingNtdBMA substModel;

    public void initAndValidate() throws Exception{
        super.initAndValidate();
        dpNtdRateSepSiteModel = dpNtdRateSepSiteModelInput.get();
        substModel = (SwitchingNtdBMA)m_siteModel.getSubstitutionModel();
    }

    public double[] calculateLogP(
            RealParameter modelParameters,
            RealParameter modelCode,
            RealParameter freqs,
            int[] siteIndexP,
            int except){

        int[] siteIndex = new int[siteIndexP.length - 1];
        int k = 0;
        for(int i = 0; i < siteIndexP.length;i++){
            if(siteIndexP[i] != except){
                siteIndex[k++] =  siteIndexP[i];
            }
        }
        return calculateLogP(modelParameters,modelCode,freqs,siteIndex);

    }


    public double[] calculateLogP(
            RealParameter modelParameters,
            RealParameter modelCode,
            RealParameter freqs,
            int[] siteIndex){
        
        Alignment alignment = m_data.get();
        double[] logPs = new double[siteIndex.length];
        setModelParameterVals(modelParameters,modelCode,freqs);

        try{
            //Stores a list of rates
            ArrayList<QuietRealParameter> rates = new ArrayList<QuietRealParameter>();
            //Stores the ID number of the substModel and pattern weights
            HashMap<Integer,int[]> clusterWeightsMap = new HashMap<Integer,int[]>();
            //Stores the ID number of the substModel and the sites with that model
            HashMap<Integer, ArrayList<Integer>> clusterSitesMap = new HashMap<Integer, ArrayList<Integer>>();

            int patternIndex;
            int rateID;

            for(int i = 0; i < siteIndex.length;i++){

                patternIndex = alignment.getPatternIndex(siteIndex[i]);
                QuietRealParameter rate = dpNtdRateSepSiteModel.getRate(siteIndex[i]);
                rateID = rate.getIDNumber();

                if(clusterWeightsMap.containsKey(rateID)){
                    //increment the weight of the pattern at site siteIndex[i] by 1
                    clusterWeightsMap.get(rateID)[patternIndex]++;
                    //add the site to the cluster
                    clusterSitesMap.get(rateID).add(i);

                }else{
                    //Create a new pattern weights vector
                    int[] patternWeights = new int[alignment.getPatternCount()];
                    //Increment the weight of the pattern at site siteIndex[i]
                    patternWeights[patternIndex]++;
                    //Add the new pattern weight vector to the HashMap
                    clusterWeightsMap.put(rateID,patternWeights);

                    //Add the new substModel to the array list
                    rates.add(rate);

                    //Create a new array list to store sites
                    ArrayList<Integer> clusterSites = new ArrayList<Integer>();
                    //Add site index to the corresponding list
                    clusterSites.add(i);
                    //Put the list in the HashMap
                    clusterSitesMap.put(rateID,clusterSites);
                }

            }

            int rateCount = rates.size();
            int site;

            for(int i = 0; i < rateCount;i++){

                m_siteModel.m_pSubstModel.setValue(substModel, m_siteModel);
                ((DummySiteModel)m_siteModel).getRateParameter().setValueQuietly(0,rates.get(i).getValue());
                rateID = rates.get(i).getIDNumber();
                setPatternWeights(clusterWeightsMap.get(rateID));
                calculateLogP();

                ArrayList<Integer> arrayIndex = clusterSitesMap.get(rateID);
                for(Integer index:arrayIndex){
                    site = siteIndex[index];
                    logPs[index] = m_fPatternLogLikelihoods[m_data.get().getPatternIndex(site)] ;
                }
            }

        }catch(Exception e){
            throw new RuntimeException(e);

        }
        return logPs;
    }

    public void setModelParameterVals(
            RealParameter modelParameters,
            RealParameter modelCode,
            RealParameter freqs){
        substModel.getLogKappa().setValueQuietly(0,modelParameters.getValue(0));
        substModel.getLogTN().setValueQuietly(0,modelParameters.getValue(1));
        substModel.getLogAC().setValueQuietly(0,modelParameters.getValue(2));
        substModel.getLogAT().setValueQuietly(0,modelParameters.getValue(3));
        substModel.getLogGC().setValueQuietly(0,modelParameters.getValue(4));
        substModel.getModelChoose().setValueQuietly(0,modelCode.getValue());
        substModel.getFreqs().setValueQuietly(0,freqs.getValue(0));
        substModel.getFreqs().setValueQuietly(1,freqs.getValue(1));
        substModel.getFreqs().setValueQuietly(2,freqs.getValue(2));
        substModel.getFreqs().setValueQuietly(3,freqs.getValue(3));
        substModel.setUpdateMatrix(true);
    }

    public double[] calculateLogP(
            RealParameter rateParameter,
            int[] siteIndexP,
            int except)throws Exception{
        int[] siteIndex = new int[siteIndexP.length - 1];
        int k = 0;
        for(int i = 0; i < siteIndexP.length;i++){
            if(siteIndexP[i] != except){
                siteIndex[k++] =  siteIndexP[i];
            }
        }
        return calculateLogP(rateParameter,siteIndex);

    }


    public double[] calculateLogP (
            RealParameter rateParameter,
            int[] siteIndex) throws Exception{
        Alignment alignment = m_data.get();

        double[] logPs = new double[siteIndex.length];



            //Stores a list of subst models
            ArrayList<SwitchingNtdBMA> ntdBMAMap = new ArrayList<SwitchingNtdBMA>();
            //Stores the ID number of the substModel and pattern weights
            HashMap<Integer,int[]> clusterWeightsMap = new HashMap<Integer,int[]>();            
            //Stores the ID number of the substModel and the sites with that model
            HashMap<Integer, ArrayList<Integer>> clusterSitesMap = new HashMap<Integer, ArrayList<Integer>>();

            int patternIndex;
            int modelID;

            for(int i = 0; i < siteIndex.length;i++){

                patternIndex = alignment.getPatternIndex(siteIndex[i]);
                //System.out.println((SwitchingNtdBMA)dpNtdRateSepSiteModel.getModel(siteIndex[i]));
                SwitchingNtdBMA ntdBMA = ((SwitchingNtdBMA)dpNtdRateSepSiteModel.getModel(siteIndex[i]));
                modelID = ntdBMA.getIDNumber();
                
                if(clusterWeightsMap.containsKey(modelID)){
                    //increment the weight of the pattern at site siteIndex[i] by 1
                    clusterWeightsMap.get(modelID)[patternIndex]++;
                    //add the site to the cluster
                    clusterSitesMap.get(modelID).add(i);
                    
                }else{
                    //Create a new pattern weights vector
                    int[] patternWeights = new int[alignment.getPatternCount()];
                    //Increment the weight of the pattern at site siteIndex[i]
                    patternWeights[patternIndex]++;
                    //Add the new pattern weight vector to the HashMap
                    clusterWeightsMap.put(modelID,patternWeights);

                    //Add the new substModel to the array list
                    ntdBMAMap.add(ntdBMA);

                    //Create a new array list to store sites
                    ArrayList<Integer> clusterSites = new ArrayList<Integer>();
                    //Add site index to the corresponding list
                    clusterSites.add(i);
                    //Put the list in the HashMap
                    clusterSitesMap.put(modelID,clusterSites);
                }
                
            }

            int substModelCount = ntdBMAMap.size();
            int site;

            for(int i = 0; i < substModelCount;i++){

                m_siteModel.m_pSubstModel.setValue(ntdBMAMap.get(i), m_siteModel);
                ((DummySiteModel)m_siteModel).getRateParameter().setValueQuietly(0,rateParameter.getValue());
                modelID = ntdBMAMap.get(i).getIDNumber();
                setPatternWeights(clusterWeightsMap.get(modelID));
                calculateLogP();

                ArrayList<Integer> arrayIndex = clusterSitesMap.get(modelID);
                for(Integer index:arrayIndex){
                    site = siteIndex[index];
                    logPs[index] = m_fPatternLogLikelihoods[m_data.get().getPatternIndex(site)] ;
                }
            }


        return logPs;

    }

    public double calculateLogP() throws Exception {
        m_substitutionModel = m_siteModel.m_pSubstModel.get();
        /*boolean[] unmasked = m_likelihoodCore.getUnmasked();
        for(int i = 0; i < patternWeights.length; i++){
            System.out.print(unmasked[i]+" ");
        }
            System.out.println();   */

        Tree tree = m_tree.get();
        //System.err.println("m_nHasDirt: "+m_nHasDirt);
        m_nHasDirt = Tree.IS_FILTHY;
       	traverse(tree.getRoot());


        calcLogP();

        m_nScale++;
        /*if (logP > 0 || (m_likelihoodCore.getUseScaling() && m_nScale > X)) {
            //System.err.println("Switch off scaling");
            m_likelihoodCore.setUseScaling(1.0);
            m_likelihoodCore.unstore();
            m_nHasDirt = Tree.IS_FILTHY;
            X *= 2;
           	traverse(tree.getRoot());
            calcLogP();
            //return logP;
        } else if (logP == Double.NEGATIVE_INFINITY && m_fScale < 10) { // && !m_likelihoodCore.getUseScaling()) {
        	m_nScale = 0;
        	m_fScale *= 1.01;
            System.err.println("Turning on scaling to prevent numeric instability " + m_fScale);
            m_likelihoodCore.setUseScaling(m_fScale);
            m_likelihoodCore.unstore();
            m_nHasDirt = Tree.IS_FILTHY;
           	traverse(tree.getRoot());
            calcLogP();
            //System.err.println(getID()+": "+logP);
            //return logP;
        }*/

        //System.err.println(getID()+": "+logP);
        return logP;
    }



}
