package beast.evolution.operators;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.*;
import beast.evolution.likelihood.DPTreeLikelihood;
import beast.evolution.likelihood.TempWVTreeLikelihood;
import beast.math.distributions.CompoundDirichletProcess;
import beast.math.distributions.ParametricDistribution;
import beast.util.Randomizer;

import java.util.HashMap;
import java.util.List;

/**
 *@author Chieh-Hsi Wu
 */
@Description("Operator that splits and merges categories of substitution model and rate. The substitution model and rate shares the same partition scheme.")
public class NtdBMARateSAMSPriorOperator extends Operator {

    public Input<DPPointer> paramPointersInput = new Input<DPPointer>(
            "parameterPointers",
            "array which points a set of unique parameter values",
            Input.Validate.REQUIRED
    );
    public Input<ParameterList> paramListInput = new Input<ParameterList>(
            "parameterList",
            "points at which the density is calculated",
            Input.Validate.REQUIRED
    );

    public Input<DPPointer> modelPointersInput = new Input<DPPointer>(
            "modelPointers",
            "array which points a set of unique parameter values",
            Input.Validate.REQUIRED
    );
    public Input<ParameterList> modelListInput = new Input<ParameterList>(
            "modelList",
            "points at which the density is calculated",
            Input.Validate.REQUIRED
    );

    public Input<DPPointer> freqsPointersInput = new Input<DPPointer>(
            "freqPointers",
            "array which points a set of unique parameter values",
            Input.Validate.REQUIRED
    );
    public Input<ParameterList> freqsListInput = new Input<ParameterList>(
            "freqsList",
            "points at which the density is calculated",
            Input.Validate.REQUIRED
    );

    public Input<DPPointer> ratesPointersInput = new Input<DPPointer>(
            "ratesPointers",
            "array which points a set of unique parameter values",
            Input.Validate.REQUIRED
    );
    public Input<ParameterList> ratesListInput = new Input<ParameterList>(
            "ratesList",
            "points at which the density is calculated",
            Input.Validate.REQUIRED
    );


    public Input<CompoundDirichletProcess> dpInput = new Input<CompoundDirichletProcess>(
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



    private DPPointer paramPointers;
    private DPPointer modelPointers;
    private DPPointer freqsPointers;
    private DPPointer ratesPointers;
    private ParameterList paramList;
    private ParameterList modelList;
    private ParameterList freqsList;
    private ParameterList ratesList;
    private int pointerCount;
    private ParametricDistribution paramBaseDistr;
    private ParametricDistribution modelBaseDistr;
    private ParametricDistribution freqsBaseDistr;
    private ParametricDistribution ratesBaseDistr;
    private CompoundDirichletProcess dp;
    private TempWVTreeLikelihood tempLikelihood;
    private DPTreeLikelihood dpTreeLikelihood;
    HashMap<Double,double[]> modelNetworkMap= new HashMap<Double,double[]>();

    public void initAndValidate(){
        paramList = paramListInput.get();
        modelList = modelListInput.get();
        freqsList = freqsListInput.get();
        ratesList = ratesListInput.get();

        paramPointers = paramPointersInput.get();
        modelPointers = modelPointersInput.get();
        freqsPointers = freqsPointersInput.get();
        ratesPointers = ratesPointersInput.get();

        pointerCount = paramPointers.getDimension();



        dp = dpInput.get();
        List<ParametricDistribution> distrs = dp.getBaseDistributions();
        paramBaseDistr = distrs.get(0);
        modelBaseDistr = distrs.get(1);
        freqsBaseDistr = distrs.get(2);
        ratesBaseDistr = distrs.get(3);
        tempLikelihood = tempLikelihoodInput.get();
        dpTreeLikelihood = dpTreeLikelihoodInput.get();

        modelNetworkMap.put(1.0,new double[]{3.0});
        modelNetworkMap.put(2.0,new double[]{3.0});
        modelNetworkMap.put(3.0,new double[]{1.0,2.0,4.0});
        modelNetworkMap.put(4.0,new double[]{3.0,5.0});
        modelNetworkMap.put(5.0,new double[]{4.0});
        //System.out.println("is null? "+(modelNetworkMap.get(5.0) == null));

    }


    public double proposal(){
        double logq = 0.0;
        //Pick two indcies at random
        int index1 = Randomizer.nextInt(pointerCount);
        int index2 = index1;
        while(index2 == index1){
            index2 = Randomizer.nextInt(pointerCount);
        }


        int clusterIndex1 = paramPointers.indexInList(index1,paramList);
        int clusterIndex2 = paramPointers.indexInList(index2,paramList);

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
        return logq;
    }

    public QuietRealParameter getSample(ParametricDistribution distr, double upper, double lower) throws Exception{

        Double[][] sampleVals = distr.sample(1);

        QuietRealParameter sampleParameter = new QuietRealParameter(sampleVals[0]);
        sampleParameter.setUpper(upper);
        sampleParameter.setLower(lower);

        return sampleParameter;
    }

    public double split(int index1, int index2, int clusterIndex, int[] initClusterSites){
        try{
            double logqSplit = 0.0;



            //Create a parameter by sampling from the prior
            QuietRealParameter newParam = getSample(paramBaseDistr, paramList.getUpper(), paramList.getLower());
            QuietRealParameter newModel = getSample(modelBaseDistr, modelList.getUpper(), modelList.getLower());
            QuietRealParameter newFreqs = getSample(freqsBaseDistr, freqsList.getUpper(), freqsList.getLower());
            QuietRealParameter newRates = getSample(ratesBaseDistr, ratesList.getUpper(), ratesList.getLower());





            //Remove the index 1 and index 2 from the cluster
            int[] clusterSites = new int[initClusterSites.length -2];
            int k = 0;
            for(int i = 0 ; i < initClusterSites.length;i++){
                if(initClusterSites[i] != index1 && initClusterSites[i] != index2){
                    clusterSites[k++] = initClusterSites[i];
                }
            }
            //Form a new cluster with index 1
            paramPointers.point(index1,newParam);
            modelPointers.point(index1,newModel);
            freqsPointers.point(index1,newFreqs);
            ratesPointers.point(index1,newRates);

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
            double[] logLik1 = tempLikelihood.calculateLogP(
                    newParam,
                    newModel,
                    newFreqs,
                    newRates,
                    clusterSites
            );
            double[] logLik2 = new double[clusterSites.length];
            for(int i = 0; i < logLik2.length; i++){
                logLik2[i] = dpTreeLikelihood.getSiteLogLikelihood(clusterIndex,clusterSites[i]);
            }



            double[] lik1 = new double[logLik1.length];
            double[] lik2 = new double[logLik2.length];

            double minLog;
            //scale it so it may be more accurate
            /*for(int i = 0; i < logLik1.length; i++){
                minLog = Math.min(logLik1[i],logLik2[i]);
                if(minLog == logLik1[i]){
                    lik1[i] = 1.0;
                    lik2[i] = Math.exp(logLik2[i] - minLog);
                }else{
                    lik1[i] = Math.exp(logLik1[i] - minLog);
                    lik2[i] = 1.0;
                }

            }*/

            for(int i = 0; i < logLik1.length; i++){
                lik1[i] = Math.exp(logLik1[i]);
                lik2[i] = Math.exp(logLik2[i]);
                //System.out.println(lik1[i]+" "+lik2[i]);
            }
            /*for(int i = 0; i < clusterSites.length;i++){
                System.out.println("clusterSites: "+clusterSites[i]);

            }
            System.out.println("index 1: "+index1+" index2: "+index2);*/

            int cluster1Count = 1;
            int cluster2Count = 1;

            int[] sitesInCluster1 = new int[initClusterSites.length];
            sitesInCluster1[0] = index1;

            //Assign members of the existing cluster (except for indice 1 and 2) randomly
            //to the existing and the new cluster
            double psi1, psi2, newClusterProb, draw;
            for(int i = 0;i < clusterSites.length; i++){



                psi1 = cluster1Count*lik1[i];
                psi2 = cluster2Count*lik2[i];
                newClusterProb = psi1/(psi1+psi2);
                draw = Randomizer.nextDouble();
                if(draw < newClusterProb){

                    //System.out.println("in new cluster: "+clusterSites[i]);
                    sitesInCluster1[cluster1Count] = clusterSites[i];
                    //paramPointers.point(clusterSites[i],newParam);
                    //modelPointers.point(clusterSites[i],newModel);
                    //freqsPointers.point(clusterSites[i],newFreqs);
                    //ratesPointers.point(clusterSites[i],newRates);
                    logqSplit += Math.log(newClusterProb);
                    cluster1Count++;
                }else{
                    logqSplit += Math.log(1.0-newClusterProb);
                    cluster2Count++;
                }

            }

            logqSplit += paramBaseDistr.calcLogP(newParam)
                    + modelBaseDistr.calcLogP(newModel)
                    + freqsBaseDistr.calcLogP(newFreqs)+
                    + ratesBaseDistr.calcLogP(newRates);

            //Perform a split
            paramList  = paramListInput.get(this);
            modelList  = modelListInput.get(this);
            freqsList  = freqsListInput.get(this);
            ratesList  = ratesListInput.get(this);
            paramPointers = paramPointersInput.get(this);
            modelPointers = modelPointersInput.get(this);
            freqsPointers = freqsPointersInput.get(this);
            ratesPointers = ratesPointersInput.get(this);

            paramList.splitParameter(clusterIndex,newParam);
            modelList.splitParameter(clusterIndex,newModel);
            freqsList.splitParameter(clusterIndex,newFreqs);
            ratesList.splitParameter(clusterIndex,newRates);
            //Form a new cluster with index 1
                paramPointers = paramPointersInput.get(this);
                modelPointers = modelPointersInput.get(this);
                freqsPointers = freqsPointersInput.get(this);
                ratesPointers = ratesPointersInput.get(this);
                for(int i = 0 ; i < cluster1Count ;i++){
                    paramPointers.point(sitesInCluster1[i],newParam);
                    modelPointers.point(sitesInCluster1[i],newModel);
                    freqsPointers.point(sitesInCluster1[i],newFreqs);
                    ratesPointers.point(sitesInCluster1[i],newRates);

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
        QuietRealParameter mergedParam = paramList.getParameter(clusterIndex2);
        QuietRealParameter mergedModel = modelList.getParameter(clusterIndex2);
        QuietRealParameter mergedFreqs = freqsList.getParameter(clusterIndex2);
        QuietRealParameter mergedRates = ratesList.getParameter(clusterIndex2);

        //Create a vector that combines the site indices of the two clusters
        int[] mergedClusterSites = new int[cluster1Sites.length+cluster2Sites.length-2];

        int k = 0;
        for(int i = 0; i < cluster1Sites.length;i++){
            //Point every member in cluster 1 to cluster 2
            //paramPointers.point(cluster1Sites[i],mergedParam);
            //modelPointers.point(cluster1Sites[i],mergedModel);
            //freqsPointers.point(cluster1Sites[i],mergedFreqs);
            //ratesPointers.point(cluster1Sites[i],mergedRates);

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
            //int[] tempWeights = dpTreeLikelihood.getClusterWeights(clusterIndex1);
            /*int[] tempWeights = new int[tempLikelihood.m_data.get().getPatternCount()];
            for(int i = 0; i < cluster1Sites.length; i++){
                int patIndex = tempLikelihood.m_data.get().getPatternIndex(cluster1Sites[i]);
                tempWeights[patIndex] = 1;
            }
            tempLikelihood.setPatternWeights(tempWeights);
            double[] cluster1SitesCluster2ParamLogLik = tempLikelihood.calculateLogP(
                    mergedParam,
                    mergedModel,
                    mergedFreqs,
                    mergedRates,
                    cluster1Sites,
                    index1
            ); */
            k = 0;
            int[] sCluster1Sites = new int[cluster1Sites.length - 1];
            for(int i = 0; i < cluster1Sites.length; i++){
                if(cluster1Sites[i] != index1){
                    sCluster1Sites[k++] = cluster1Sites[i];
                }
            }
            tempLikelihood.setupPatternWeightsFromSites(sCluster1Sites);
            double[] cluster1SitesCluster2ParamLogLik = tempLikelihood.calculateLogP(
                    mergedParam,
                    mergedModel,
                    mergedFreqs,
                    mergedRates,
                    sCluster1Sites
            );

            //tempWeights = dpTreeLikelihood.getClusterWeights(clusterIndex2);
            /*tempWeights = new int[tempLikelihood.m_data.get().getPatternCount()];
            for(int i = 0; i < cluster2Sites.length; i++){
                int patIndex = tempLikelihood.m_data.get().getPatternIndex(cluster2Sites[i]);
                tempWeights[patIndex] = 1;
            }
            tempLikelihood.setPatternWeights(tempWeights);
            RealParameter removedParam = paramList.getParameter(clusterIndex1);
            RealParameter removedModel = modelList.getParameter(clusterIndex1);
            RealParameter removedFreqs = freqsList.getParameter(clusterIndex1);
            RealParameter removedRates = ratesList.getParameter(clusterIndex1);
            double[] cluster2SitesCluster1ParamLogLik = tempLikelihood.calculateLogP(
                    removedParam,
                    removedModel,
                    removedFreqs,
                    removedRates,
                    cluster2Sites,
                    index2
            );*/
            k = 0;
            int[] sCluster2Sites = new int[cluster2Sites.length - 1];
            for(int i = 0; i < cluster2Sites.length; i++){
                if(cluster2Sites[i] != index2){
                    sCluster2Sites[k++] = cluster2Sites[i];
                }
            }
            tempLikelihood.setupPatternWeightsFromSites(sCluster2Sites);
            RealParameter removedParam = paramList.getParameter(clusterIndex1);
            RealParameter removedModel = modelList.getParameter(clusterIndex1);
            RealParameter removedFreqs = freqsList.getParameter(clusterIndex1);
            RealParameter removedRates = ratesList.getParameter(clusterIndex1);
            double[] cluster2SitesCluster1ParamLogLik = tempLikelihood.calculateLogP(
                    removedParam,
                    removedModel,
                    removedFreqs,
                    removedRates,
                    sCluster2Sites
            );


            //System.out.println("populate logLik1:");
            double[] logLik1 = new double[mergedClusterSites.length];
            for(int i = 0; i < (cluster1Sites.length-1); i++){
                //System.out.println(clusterIndex1+" "+mergedClusterSites[i]);

                 logLik1[i] = dpTreeLikelihood.getSiteLogLikelihood(clusterIndex1,mergedClusterSites[i]);
            }
            System.arraycopy(cluster2SitesCluster1ParamLogLik,0,logLik1,cluster1Sites.length-1,cluster2SitesCluster1ParamLogLik.length);

            double[] logLik2 = new double[mergedClusterSites.length];
            System.arraycopy(cluster1SitesCluster2ParamLogLik,0,logLik2,0,cluster1SitesCluster2ParamLogLik.length);

            //System.out.println("populate logLik2:");
            for(int i = cluster1SitesCluster2ParamLogLik.length; i < logLik2.length; i++){
                //System.out.println(clusterIndex2+" "+mergedClusterSites[i-cluster1SitesCluster2ParamLogLik.length]);
                logLik2[i] = dpTreeLikelihood.getSiteLogLikelihood(clusterIndex2,mergedClusterSites[i]);
            }



            double[] lik1 = new double[logLik1.length];
            double[] lik2 = new double[logLik2.length];

            //scale it so it may be more accuate
            double minLog;
            /*for(int i = 0; i < logLik1.length; i++){
                minLog = Math.min(logLik1[i],logLik2[i]);
                if(minLog == logLik1[i]){
                    lik1[i] = 1.0;
                    lik2[i] = Math.exp(logLik2[i] - minLog);
                }else{
                    lik1[i] = Math.exp(logLik1[i] - minLog);
                    lik2[i] = 1.0;
                }

            }*/

            for(int i = 0; i < logLik1.length;i++){

                lik1[i] = Math.exp(logLik1[i]);
                lik2[i] = Math.exp(logLik2[i]);
                //System.out.println(lik1[i]+" "+lik2[i]);
            }

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

                 /*testCorrectness(i,cluster,
                        clusterIndex1,clusterIndex2,shuffle, mergedClusterSites,
                         lik1,lik2);*/

                cluster1Prob = psi1/(psi1+psi2);
                //System.out.println(cluster1Prob);
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

            logqMerge += paramBaseDistr.calcLogP(removedParam)
                    + modelBaseDistr.calcLogP(removedModel)
                    + freqsBaseDistr.calcLogP(removedFreqs)
                    + ratesBaseDistr.calcLogP(removedRates);
            if(logqMerge > Double.NEGATIVE_INFINITY){
                paramList.mergeParameter(clusterIndex1,clusterIndex2);
                modelList.mergeParameter(clusterIndex1,clusterIndex2);
                freqsList.mergeParameter(clusterIndex1,clusterIndex2);
                ratesList.mergeParameter(clusterIndex1,clusterIndex2);
                for(int i = 0; i < cluster1Sites.length;i++){
                    //Point every member in cluster 1 to cluster 2
                    paramPointers.point(cluster1Sites[i],mergedParam);
                    modelPointers.point(cluster1Sites[i],mergedModel);
                    freqsPointers.point(cluster1Sites[i],mergedFreqs);
                    ratesPointers.point(cluster1Sites[i],mergedRates);


                }
            }
        }catch(Exception e){
            throw new RuntimeException(e);
        }

        return logqMerge;

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

        int[] tempWeights = new int[tempLikelihood.dataInput.get().getPatternCount()];
        tempWeights[tempLikelihood.dataInput.get().getPatternIndex(mergedClusterSites[shuffle[i]])] = 1;
        tempLikelihood.setPatternWeights(tempWeights);
        double temp1 = Math.exp(tempLikelihood.calculateLogP(
                paramList.getParameter(clusterIndex1),
                modelList.getParameter(clusterIndex1),
                freqsList.getParameter(clusterIndex1),
                ratesList.getParameter(clusterIndex1))
        );
        double temp2 = Math.exp(tempLikelihood.calculateLogP(
                paramList.getParameter(clusterIndex2),
                modelList.getParameter(clusterIndex2),
                freqsList.getParameter(clusterIndex2),
                ratesList.getParameter(clusterIndex2))
        );
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
