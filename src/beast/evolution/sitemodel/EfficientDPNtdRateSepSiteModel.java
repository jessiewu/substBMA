package beast.evolution.sitemodel;

import beast.core.Description;

import java.util.ArrayList;

/**
 * @author Chieh-Hsi Wu
 */
@Description("An incomplete attempt to improve the efficiency of computation.")
public class EfficientDPNtdRateSepSiteModel extends DPNtdRateSepSiteModel{

    public void updateMap(int siteIndex, int ntdBMACluster, int rateCluster) throws Exception{
       updatedSites.add(siteIndex);
        /*for(int i = 0; i < clusterMap.length;i++){
            System.out.print("NTDBMA: ");
            for(int j = 0; j < clusterMap[i].length;j++){
                System.out.print(clusterMap[i][j]+" ");
            }
            System.out.println();
        }

        for(int i = 0; i < storedClusterMap.length;i++){
            System.out.print("stored NTDBMA: ");
            for(int j = 0; j < storedClusterMap[i].length;j++){
                System.out.print(storedClusterMap[i][j]+" ");
            }
            System.out.println();
        }*/

        mappingChanged = true;

        //Update the weight matrix
        moveWeight(
                clusterMap[NTDBMA][siteIndex],  //current substModel
                clusterMap[RATES][siteIndex],   //current rate
                ntdBMACluster,                  //potentially new substModel
                rateCluster,                    //potentially new rate
                1
        );

        int prevNtdBMAId = clusterMap[NTDBMA][siteIndex];
        int prevRateId = clusterMap[RATES][siteIndex];

        //Update mapping
        clusterMap[NTDBMA][siteIndex] = ntdBMACluster;
        clusterMap[RATES][siteIndex] = rateCluster;

        //Stored map
        storedClusterMap[NTDBMA][siteIndex] = prevNtdBMAId;
        storedClusterMap[RATES][siteIndex] = prevRateId;

        //System.out.println("flag1: "+siteModelWeights[prevNtdBMAId][prevRateId]);
        //Remove site model if it has zero pattern weight
        if(siteModelWeights[prevNtdBMAId][prevRateId] == 0){
            //System.out.println("remove");
            //System.out.println("siteModels: "+siteModels.size());
            //System.out.println(siteModelsMatrix[prevNtdBMAId][prevRateId]);
            siteModels.remove(siteModels.indexOf(siteModelsMatrix[prevNtdBMAId][prevRateId]));
            siteModelsMatrix[prevNtdBMAId][prevRateId] = null;
        }

        //Add site model if this is a new substModel and rate combination
        if(siteModelsMatrix[clusterMap[NTDBMA][siteIndex]][clusterMap[RATES][siteIndex]] == null){
            //System.out.println("add: "+clusterMap[NTDBMA][siteIndex]+" "+clusterMap[RATES][siteIndex]);
            QuietSiteModel siteModel =
                    new QuietSiteModel(dpNtdBMA.getModel(dpNtdBMA.getCurrCluster(siteIndex)),
                    ratesPointers.getParameter(siteIndex));

            siteModelsMatrix[clusterMap[NTDBMA][siteIndex]][clusterMap[RATES][siteIndex]] = siteModel;
            siteModels.add(siteModel);
        }

        /*for(int i = 0; i < clusterMap.length;i++){
            System.out.print("RATE: ");
            for(int j = 0; j < clusterMap[i].length;j++){
                System.out.print(clusterMap[i][j]+" ");
            }
            System.out.println();
        }
        for(int i = 0; i < storedClusterMap.length;i++){
            System.out.print("stored RATE: ");
            for(int j = 0; j < storedClusterMap[i].length;j++){
                System.out.print(storedClusterMap[i][j]+" ");
            }
            System.out.println();
        }*/




    }

    ArrayList<Integer> updatedSites;
    public void store(){
        updatedSites = new ArrayList<Integer>();
        changedInput = -1;
        /*for(int i = 0; i < storedClusterMap.length;i++){
            System.arraycopy(clusterMap[i],0,storedClusterMap[i],0,clusterMap[i].length);
        }
        for(int i = 0; i < storedSiteModelWeights.length;i++){
            System.arraycopy(siteModelWeights[i],0,storedSiteModelWeights[i],0,siteModelWeights[i].length);
        }

        for(int i = 0; i < storedSiteModelsMatrix.length;i++){
            System.arraycopy(siteModelsMatrix[i],0,storedSiteModelsMatrix[i],0,siteModelsMatrix[i].length);
        }*/
        //for(int i = 0; i < clusterMap.length; i++){
        //    System.arraycopy(clusterMap[i],0,storedClusterMap[i],0,clusterMap[i].length);
        //}
        super.store();
    }

    public void backUp(){
        for(int i = 0; i < storedSiteModelWeights.length;i++){
            System.arraycopy(siteModelWeights[i],0,storedSiteModelWeights[i],0,siteModelWeights[i].length);
        }

        for(int i = 0; i < storedSiteModelsMatrix.length;i++){
            System.arraycopy(siteModelsMatrix[i],0,storedSiteModelsMatrix[i],0,siteModelsMatrix[i].length);
        }
    }


    int affectedCluster, remainedCluster;

    void handleMerge(int changedInput) throws Exception{
        backUp();
        /*if(changedInput == NTDBMA){

            affectedCluster = dpNtdBMA.getRemovedModelIDNumber();
            remainedCluster = dpNtdBMA.getDirtyModelIDNumber();
            System.arraycopy(siteModelWeights[affectedCluster],0,storedSiteModelWeights[affectedCluster],0,siteModelWeights[affectedCluster].length);
            System.arraycopy(siteModelWeights[remainedCluster],0,storedSiteModelWeights[remainedCluster],0,siteModelWeights[remainedCluster].length);

        }else if(changedInput == RATES){
            //The rate cluster that takes up the member of the removed rate cluster
            affectedCluster = ratesList.getRemovedIDNumber();
            remainedCluster = ratesList.getDirtyParameterIDNumber();
            for(int i = 0; i < siteModelWeights.length;i++){
                storedSiteModelWeights[i][affectedCluster] = siteModelWeights[i][affectedCluster];
            }
            for(int i = 0; i < siteModelWeights.length;i++){
                storedSiteModelWeights[i][remainedCluster] = siteModelWeights[i][remainedCluster];
            }

        }else{
            throw new RuntimeException("Can only add clusters for either ntdBMA model or rates, changedInput: "+changedInput);
        } */

        super.handleMerge(changedInput);

    }

    protected void handleSplit(int changedInput) throws Exception{
        backUp();
        /*int ntdBMAIDNum;
        int muIDNum;
        if(changedInput == NTDBMA){

            affectedCluster = dpNtdBMA.getLastModel().getIDNumber();
            remainedCluster = dpNtdBMA.getDirtyModelIDNumber();
            System.arraycopy(siteModelWeights[affectedCluster],0,storedSiteModelWeights[affectedCluster],0,siteModelWeights[affectedCluster].length);
            System.arraycopy(siteModelWeights[remainedCluster],0,storedSiteModelWeights[remainedCluster],0,siteModelWeights[remainedCluster].length);

            //Retreive the new NtdBMA model
            ntdBMAIDNum = dpNtdBMA.getLastModel().getIDNumber();
            //Retrieve all the (dirty) sites that uses the newly added model
            clusterSites = dpValSubst.getClusterSites(dpNtdBMA.getDimension()-1);

            for(int dirtySite: clusterSites){
                //Get the rate of the each of the dirty sites
                muIDNum = ratesPointers.getParameter(dirtySite).getIDNumber();
                //Update the substModel and rate combination
                updateMap(dirtySite,ntdBMAIDNum,muIDNum);
            }

        }else if(changedInput == RATES){
            remainedCluster = ratesList.getDirtyParameterIDNumber();
            affectedCluster = ratesList.getLastParameter().getIDNumber();
            for(int i = 0; i < siteModelWeights.length;i++){
                storedSiteModelWeights[i][affectedCluster] = siteModelWeights[i][affectedCluster];
            }
            for(int i = 0; i < siteModelWeights.length;i++){
                storedSiteModelWeights[i][remainedCluster] = siteModelWeights[i][remainedCluster];
            }

            muIDNum = ratesList.getLastParameter().getIDNumber();
            clusterSites = dpValRate.getClusterSites(ratesList.getDimension()-1);
            for(int dirtySite:clusterSites){
                ntdBMAIDNum = dpNtdBMA.getModel(dpNtdBMA.getCurrCluster(dirtySite)).getIDNumber();
                updateMap(dirtySite,ntdBMAIDNum,muIDNum);
            }
        }else{
            throw new RuntimeException("Can only add clusters for either ntdBMA model or rates, changedInput: "+changedInput);
        }*/
        super.handleSplit(changedInput);

    }

    public void restore(){
        if(mappingChanged){
            for(int updatedSite:updatedSites){
                clusterMap[NTDBMA][updatedSite] = storedClusterMap[NTDBMA][updatedSite];
                clusterMap[RATES][updatedSite] = storedClusterMap[RATES][updatedSite];

            }
            /*if(changedInput == NTDBMA){
                int[] tmp1 = siteModelWeights[affectedCluster];
                siteModelWeights[affectedCluster] = storedSiteModelWeights[affectedCluster];
                storedSiteModelWeights[affectedCluster] = tmp1;

                int[] tmp2 = siteModelWeights[remainedCluster];
                siteModelWeights[remainedCluster] = storedSiteModelWeights[remainedCluster];
                storedSiteModelWeights[remainedCluster] = tmp2;

            }else{
                for(int i = 0; i < siteModelWeights.length;i++){
                    siteModelWeights[i][affectedCluster] = storedSiteModelWeights[i][affectedCluster];
                }
                for(int i = 0; i < siteModelWeights.length;i++){
                    siteModelWeights[i][remainedCluster] = storedSiteModelWeights[i][remainedCluster];
                }


            }*/
            int[][] tmp1 = siteModelWeights;
            siteModelWeights = storedSiteModelWeights;
            storedSiteModelWeights = tmp1;

            QuietSiteModel[][] tmp2 = siteModelsMatrix;
            siteModelsMatrix = storedSiteModelsMatrix;
            storedSiteModelsMatrix = tmp2;
        }

        super.restore();
    }

}
