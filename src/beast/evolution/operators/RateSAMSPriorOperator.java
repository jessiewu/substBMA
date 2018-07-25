package beast.evolution.operators;

import beast.core.Description;
import beast.core.Operator;
import beast.core.Input;
import beast.core.parameter.*;
import beast.evolution.sitemodel.DPNtdRateSepSiteModel;
import beast.math.distributions.DirichletProcess;
import beast.math.distributions.ParametricDistribution;
import beast.evolution.likelihood.*;
import beast.util.Randomizer;

import java.util.HashMap;

/**
 * @author Chieh-Hsi Wu
 */
@Description("Operator that splits and merges categories of substitution rate.")
public class RateSAMSPriorOperator extends Operator {
    public Input<DPPointer> ratePointersInput = new Input<DPPointer>(
            "ratesPointers",
            "array which points a set of unique parameter values",
            Input.Validate.REQUIRED
    );
    public Input<ParameterList> rateListInput = new Input<ParameterList>(
            "ratesList",
            "points at which the density is calculated",
            Input.Validate.REQUIRED
    );


    public Input<DirichletProcess> dpInput = new Input<DirichletProcess>(
            "dirichletProcess",
            "An object of Dirichlet Process",
            Input.Validate.REQUIRED
    );

    public Input<TempWVTreeLikelihood> tempLikelihoodInput = new Input<TempWVTreeLikelihood>(
            "tempLikelihood",
            "The temporary likelihood given the data at site i",
            Input.Validate.REQUIRED
    );

    public Input<DPValuable> dpValuableInput = new Input<DPValuable>(
            "dpVal",
            "reports the counts in each cluster",
            Input.Validate.REQUIRED
    );

    public Input<DPTreeLikelihood> dpTreeLikelihoodInput = new Input<DPTreeLikelihood>(
            "dpTreeLik",
            "Tree likelihood that handle DPP",
            Input.Validate.REQUIRED
    );



    private DPPointer ratePointers;
    private ParameterList rateList;
    private int pointerCount;
    private ParametricDistribution baseDistr;
    private DirichletProcess dp;
    private TempWVTreeLikelihood tempLikelihood;
    private DPTreeLikelihood dpTreeLikelihood;
    public void initAndValidate(){
        ratePointers = ratePointersInput.get();
        pointerCount = ratePointers.getDimension();
        rateList = rateListInput.get();
        dp = dpInput.get();
        baseDistr = dp.getBaseDistribution();
        tempLikelihood = tempLikelihoodInput.get();
        dpTreeLikelihood = dpTreeLikelihoodInput.get();
    }

    
    public double proposal(){
        double logq = 0.0;
        //Pick two indcies at random
        int index1 = Randomizer.nextInt(pointerCount);
        int index2 = index1;
        while(index2 == index1){
            index2 = Randomizer.nextInt(pointerCount);
        }


        int clusterIndex1 = ratePointers.indexInList(index1,rateList);
        int clusterIndex2 = ratePointers.indexInList(index2,rateList);

        //If the randomly draw sites are from the same cluster, perform a split-move.
        if(clusterIndex1 == clusterIndex2){

            int[] clusterSites = dpValuableInput.get().getClusterSites(clusterIndex1);

            double temp = split(index1, index2,clusterIndex1,clusterSites);
            //System.out.println("split: "+temp);
            logq += temp;


        }else{
            //If the the two randomly drawn sites are not from the same cluster, perform a merge-move.

            int[] cluster1Sites = dpValuableInput.get().getClusterSites(clusterIndex1);
            int[] cluster2Sites = dpValuableInput.get().getClusterSites(clusterIndex2);

            //logq = merge(index1, index2,clusterIndex1,clusterIndex2,cluster1Sites,cluster2Sites);
            double temp = merge(
                    index1,
                    index2,
                    clusterIndex1,
                    clusterIndex2,
                    cluster1Sites,
                    cluster2Sites
            );

            //System.out.println("merge: "+temp);
            logq = temp;

        }
        //System.out.println("logq: "+logq);
        return logq;
    }

    public double split(int index1, int index2, int clusterIndex, int[] initClusterSites){
        try{
            double logqSplit = 0.0;

            //Create a parameter by sampling from the prior
            //Double[][] sampleVal = baseDistr.sample(1,rateList.getParameter(clusterIndex).getUpper(),rateList.getParameter(clusterIndex).getLower());
            // keep the original code, and not use sampler in QuietRealParameter#getSample
            Double[][] sampleVal = baseDistr.sample(1);

            QuietRealParameter newParam = new QuietRealParameter(sampleVal[0]);
            newParam.setLower(rateList.getParameter(clusterIndex).getLower());
            newParam.setUpper(rateList.getParameter(clusterIndex).getUpper());


            //Remove the index 1 and index 2 from the cluster
            int[] clusterSites = new int[initClusterSites.length -2];
            int k = 0;
            for(int i = 0 ; i < initClusterSites.length;i++){
                if(initClusterSites[i] != index1 && initClusterSites[i] != index2){
                    clusterSites[k++] = initClusterSites[i];
                }
            }


            //Form a new cluster with index 1
            //ratePointers.point(index1,newParam);

            int[] sitesInCluster1 = new int[initClusterSites.length];
            sitesInCluster1[0] = index1;
            //Shuffle the cluster_-{index_1,index_2} to obtain a random permutation
            Randomizer.shuffle(clusterSites);

            //Create the weight vector of site patterns according to the order of the shuffled index.
            /*int[] tempWeights = new int[tempLikelihood.m_data.get().getPatternCount()];
            int patIndex;
            for(int i = 0; i < clusterSites.length; i++){
                patIndex = tempLikelihood.m_data.get().getPatternIndex(clusterSites[i]);
                tempWeights[patIndex] = 1;
            }

            tempLikelihood.setPatternWeights(tempWeights);*/
            tempLikelihood.setupPatternWeightsFromSites(clusterSites);

            //Site log likelihoods in the order of the shuffled sites
            double[] logLik1 = tempLikelihood.calculateLogP(newParam,clusterSites);
            double[] logLik2 = new double[clusterSites.length];
            for(int i = 0; i < logLik2.length; i++){
                //logLik2[i] = dpTreeLikelihood.getSiteLogLikelihood(clusterIndex,clusterSites[i]);
                logLik2[i] = getSiteLogLikelihood(
                        rateList.getParameter(clusterIndex),
                        clusterIndex,
                        clusterSites[i]
                );
            }



            double[] lik1 = new double[logLik1.length];
            double[] lik2 = new double[logLik2.length];

            double minLog;
            //scale it so it may be more accurate
            for(int i = 0; i < logLik1.length; i++){
                minLog = Math.min(logLik1[i],logLik2[i]);
                if(minLog == logLik1[i]){
                    lik1[i] = 1.0;
                    lik2[i] = Math.exp(logLik2[i] - minLog);
                }else{
                    lik1[i] = Math.exp(logLik1[i] - minLog);
                    lik2[i] = 1.0;
                }

            }

            /*for(int i = 0; i < logLik1.length; i++){
                lik1[i] = Math.exp(logLik1[i]);
                lik2[i] = Math.exp(logLik2[i]);
            }*/
            /*for(int i = 0; i < clusterSites.length;i++){
                System.out.println("clusterSites: "+clusterSites[i]);

            }
            System.out.println("index 1: "+index1+" index2: "+index2);*/

            int cluster1Count = 1;
            int cluster2Count = 1;

            //Assign members of the existing cluster (except for indice 1 and 2) randomly
            //to the existing and the new cluster
            double psi1, psi2, newClusterProb, draw;

            for(int i = 0;i < clusterSites.length; i++){

                psi1 = cluster1Count*lik1[i];
                psi2 = cluster2Count*lik2[i];
                //System.out.println(lik1[i]+" "+lik2[i]);
                newClusterProb = psi1/(psi1+psi2);
                draw = Randomizer.nextDouble();
                if(draw < newClusterProb){
                    //System.out.println("in new cluster: "+clusterSites[i]);
                    //ratePointers.point(clusterSites[i],newParam);
                    sitesInCluster1[cluster1Count] = clusterSites[i];
                    logqSplit += Math.log(newClusterProb);
                    cluster1Count++;
                }else{
                    logqSplit += Math.log(1.0-newClusterProb);
                    cluster2Count++;
                }
                
            }
            //System.out.println("logqSplit: "+logqSplit);
            logqSplit += baseDistr.logDensity(sampleVal[0][0]);

            if(-logqSplit != Double.NEGATIVE_INFINITY){
                //Perform a split
                rateList.splitParameter(clusterIndex,newParam);
                for(int i = 0; i < cluster1Count; i++){
                    ratePointers.point(sitesInCluster1[i],newParam);
                }
            }

            return -logqSplit;


        }catch(Exception e){
            throw new RuntimeException(e);
        }
    }

    public double merge(
            int index1,
            int index2,
            int clusterIndex1,
            int clusterIndex2,
            int[] cluster1Sites,
            int[] cluster2Sites){

        double logqMerge = 0.0;


        HashMap<Integer,Integer> siteMap = new HashMap<Integer, Integer>();

        //The value of the merged cluster will have that of cluster 2 before the merge.
        QuietRealParameter mergedParameter = rateList.getParameter(clusterIndex2);

        //Create a vector that combines the site indices of the two clusters
        int[] mergedClusterSites = new int[cluster1Sites.length+cluster2Sites.length-2];

        int k = 0;
        for(int i = 0; i < cluster1Sites.length;i++){
            //Point every member in cluster 1 to cluster 2
            //ratePointers.point(cluster1Sites[i],mergedParameter);
            if(cluster1Sites[i] != index1){
                // For all members that are not index 1,
                // record the cluster in which they have been before the merge,
                // and assign them to the combined vector.
                siteMap.put(cluster1Sites[i],clusterIndex1);
                mergedClusterSites[k++] = cluster1Sites[i];
            }
        }


        for(int i = 0; i < cluster2Sites.length;i++){
            //All members in cluster 2 remains in cluster2 so no new pointer assignments
            if(cluster2Sites[i] != index2){
                // For all members that are not index 2,
                // record the cluster in which they have been before the merge,
                // and assign them to the combined vector.
                siteMap.put(cluster2Sites[i],clusterIndex2);
                mergedClusterSites[k++] = cluster2Sites[i];
            }
        }



        
        try{

            // Create a weight vector of patterns to inform the temporary tree likelihood
            // which set of pattern likelihoods are to be computed.
            /*int[] tempWeights = new int[tempLikelihood.m_data.get().getPatternCount()];
            for(int i = 0; i < cluster1Sites.length; i++){
                int patIndex = tempLikelihood.m_data.get().getPatternIndex(cluster1Sites[i]);
                tempWeights[patIndex] = 1;
            }
            tempLikelihood.setPatternWeights(tempWeights);*/

            tempLikelihood.setupPatternWeightsFromSites(cluster1Sites);

            double[] cluster1SitesCluster2ParamLogLik = tempLikelihood.calculateLogP(
                    rateList.getParameter(clusterIndex2),
                    cluster1Sites,
                    index1
            );

            //tempWeights = dpTreeLikelihood.getClusterWeights(clusterIndex2);
            /*tempWeights = new int[tempLikelihood.m_data.get().getPatternCount()];
            for(int i = 0; i < cluster2Sites.length; i++){
                int patIndex = tempLikelihood.m_data.get().getPatternIndex(cluster2Sites[i]);
                tempWeights[patIndex] = 1;
            }
            tempLikelihood.setPatternWeights(tempWeights);*/
            tempLikelihood.setupPatternWeightsFromSites(cluster2Sites);


            double[] cluster2SitesCluster1ParamLogLik = tempLikelihood.calculateLogP(
                    rateList.getParameter(clusterIndex1),
                    cluster2Sites,
                    index2
            );

            //System.out.println("populate logLik1:");
            double[] logLik1 = new double[mergedClusterSites.length];
            for(int i = 0; i < (cluster1Sites.length-1); i++){
                //System.out.println(clusterIndex1+" "+mergedClusterSites[i]);
                //logLik1[i] = dpTreeLikelihood.getSiteLogLikelihood(clusterIndex1,mergedClusterSites[i]);
                logLik1[i] = getSiteLogLikelihood(
                        rateList.getParameter(clusterIndex1),
                        clusterIndex1,
                        mergedClusterSites[i]
                );
            }
            System.arraycopy(cluster2SitesCluster1ParamLogLik,0,logLik1,cluster1Sites.length-1,cluster2SitesCluster1ParamLogLik.length);

            double[] logLik2 = new double[mergedClusterSites.length];
            System.arraycopy(cluster1SitesCluster2ParamLogLik,0,logLik2,0,cluster1SitesCluster2ParamLogLik.length);

            //System.out.println("populate logLik2:");
            for(int i = cluster1SitesCluster2ParamLogLik.length; i < logLik2.length; i++){
                //System.out.println(clusterIndex2+" "+mergedClusterSites[i-cluster1SitesCluster2ParamLogLik.length]);
                //logLik2[i] = dpTreeLikelihood.getSiteLogLikelihood(clusterIndex2,mergedClusterSites[i]);
                logLik2[i] = getSiteLogLikelihood(
                        rateList.getParameter(clusterIndex2),
                        clusterIndex2,
                        mergedClusterSites[i]
                );
            }



            double[] lik1 = new double[logLik1.length];
            double[] lik2 = new double[logLik2.length];

            //scale it so it may be more accurate
            double minLog;
            for(int i = 0; i < logLik1.length; i++){
                minLog = Math.min(logLik1[i],logLik2[i]);
                if(minLog == logLik1[i]){
                    lik1[i] = 1.0;
                    lik2[i] = Math.exp(logLik2[i] - minLog);
                }else{
                    lik1[i] = Math.exp(logLik1[i] - minLog);
                    lik2[i] = 1.0;
                }

            }

            /*for(int i = 0; i < logLik1.length;i++){
                lik1[i] = Math.exp(logLik1[i]);
                lik2[i] = Math.exp(logLik2[i]);
            }*/

            //Create a set of indices for random permutation
            int[] shuffle = new int[mergedClusterSites.length];
            for(int i = 0; i < shuffle.length;i++){
                shuffle[i] = i;
            }
            Randomizer.shuffle(shuffle);

            int cluster1Count = 1;
            int cluster2Count = 1;
            int cluster;
            double psi1, psi2, cluster1Prob;
            for(int i = 0; i < mergedClusterSites.length;i++){
                cluster = siteMap.get(mergedClusterSites[shuffle[i]]);
                psi1 = cluster1Count*lik1[shuffle[i]];
                psi2 = cluster2Count*lik2[shuffle[i]];



                cluster1Prob = psi1/(psi1+psi2);
                if(cluster == clusterIndex1){
                    logqMerge += Math.log(cluster1Prob);
                    cluster1Count++;

                }else if(cluster == clusterIndex2){
                    logqMerge += Math.log(1-cluster1Prob);
                    cluster2Count++;

                }else{
                    throw new RuntimeException("Something is wrong.");
                }

            }
            logqMerge += baseDistr.calcLogP(rateList.getParameter(clusterIndex1));

        }catch(Exception e){
            throw new RuntimeException(e);
        }

        if(logqMerge != Double.NEGATIVE_INFINITY){
            rateList.mergeParameter(clusterIndex1,clusterIndex2);
            for(int i = 0; i < cluster1Sites.length;i++){
                //Point every member in cluster 1 to cluster 2
                ratePointers.point(cluster1Sites[i],mergedParameter);

            }
        }

        return logqMerge;

    }

    public double getSiteLogLikelihood(
            QuietRealParameter parameter,
            int clusterIndex,
            int siteIndex){

        double siteLogLik;
        if(dpTreeLikelihood instanceof DPSepTreeLikelihood){
            siteLogLik = ((DPSepTreeLikelihood)dpTreeLikelihood).getSiteLogLikelihood(
                    DPNtdRateSepSiteModel.RATES,
                    parameter.getIDNumber(),
                    siteIndex
            );
            //System.out.println("1: "+siteIndex+" "+siteLogLik);
        }else if(dpTreeLikelihood instanceof SlowDPSepTreeLikelihood){
            siteLogLik = ((SlowDPSepTreeLikelihood)dpTreeLikelihood).getSiteLogLikelihood(
                    DPNtdRateSepSiteModel.RATES,
                    parameter.getIDNumber(),
                    siteIndex
            );
       }else{
           siteLogLik =  dpTreeLikelihood.getSiteLogLikelihood(clusterIndex,siteIndex);
           //System.out.println("2: "+siteIndex+" "+siteLogLik);
       }
       return siteLogLik;
    }

    public void testCorrectness(
            int i,
            int cluster,
            int clusterIndex1,
            int clusterIndex2,
            int[] shuffle,
            int[] mergedClusterSites,
            double[] lik1,
            double[] lik2) throws Exception{

        //int[] tempWeights = new int[tempLikelihood.m_data.get().getPatternCount()];
        //tempWeights[tempLikelihood.m_data.get().getPatternIndex(mergedClusterSites[shuffle[i]])] = 1;
        //tempLikelihood.setupPatternWeights(tempWeights);
        tempLikelihood.setupPatternWeightsFromSites(new int[]{mergedClusterSites[shuffle[i]]});
        double temp1 = Math.exp(tempLikelihood.calculateLogP(rateList.getParameter(clusterIndex1)));
        double temp2 = Math.exp(tempLikelihood.calculateLogP(rateList.getParameter(clusterIndex2)));
        if(temp1 != lik1[shuffle[i]] || temp2 != lik2[shuffle[i]]){
            System.out.println("temp1");
            System.out.println("shuffle_i: "+shuffle[i]);
            System.out.println("mergedClusterSites[shuffle]: "+mergedClusterSites[shuffle[i]]);
            System.out.println("cluster: "+cluster);
            System.out.println(+mergedClusterSites.length+" "+lik1.length);
            for(int j = 0; j < lik1.length;j++){
                System.out.println("merged lik1: "+mergedClusterSites[j]+" "+lik1[j]);
            }
            for(int j = 0; j < lik2.length;j++){
                System.out.println("merged lik2: "+mergedClusterSites[j]+" "+lik2[j]);
            }
            throw new RuntimeException(temp1+" "+lik1[shuffle[i]]+" "+temp2+" "+lik2[shuffle[i]]);

        }

    }

}
