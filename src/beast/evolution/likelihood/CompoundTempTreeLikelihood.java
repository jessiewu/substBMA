package beast.evolution.likelihood;

import beast.app.BeastMCMC;
import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.State;
import beast.core.parameter.QuietRealParameter;
import beast.core.parameter.RealParameter;
import beast.core.util.CompoundDistribution;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.Executors;
import java.util.concurrent.RejectedExecutionException;

/**
 * @author Chieh-Hsi Wu
 */
@Description("This is a collection of temporary tree likelihoods that attempts to facilitate multithreading computation.")
public class CompoundTempTreeLikelihood extends CompoundDistribution {
    public Input<List<TempTreeLikelihood>> tempTreeLiksInput = new Input<List<TempTreeLikelihood>>(
            "tempTreeLikelihood",
            "A list of TempTreeLikelihood objects",
            new ArrayList<TempTreeLikelihood>(),
            Input.Validate.REQUIRED
    );

    boolean useThreads;

    private List<TempTreeLikelihood> tempTreeLiks;
    public void initAndValidate() {
        tempTreeLiks = tempTreeLiksInput.get();
        useThreads = useThreadsInput.get() && (BeastMCMC.m_nThreads > 1);


    }

    class CoreRunnable implements Runnable {
        TempTreeLikelihood distr;
        int site;

        /*CoreRunnable(TempTreeLikelihood core, int site) {
            distr = core;
            this.site = site;
        }*/

        int start, end;

        CoreRunnable(int start, int end, int site) {
            this.site = site;
            this.start = start;
            this.end = end + 1;
        }

        public void run() {
            try {
                //System.out.println("started: "+this);
                for(int i = start; i < end; i++){
                    tempTreeLiks.get(i).calculateLogP(site);  //todo

                }
                //System.out.println("stopped: "+this);

            } catch (Exception e) {
                System.err.println("Something went wrong in a calculation of " + distr.getID());
                e.printStackTrace();
                System.exit(0);
            }

            m_nCountDown.countDown();
        }

    } // CoreRunnable

    CountDownLatch m_nCountDown;

    private double[] calculateLogPUsingThreads(int site) throws Exception {

        try {

            //int nrOfDirtyDistrs = tempTreeLiks.size();
            /*for (Distribution dists : tempTreeLiks) {
                if (dists.isDirtyCalculation()) {
                    nrOfDirtyDistrs++;
                }
            }*/
            //m_nCountDown = new CountDownLatch(nrOfDirtyDistrs);
            m_nCountDown = new CountDownLatch(BeastMCMC.m_nThreads);
            // kick off the threads
            int batchSize = tempTreeLiks.size()/ BeastMCMC.m_nThreads;
            int start, end;
            for (int i = 0; i < BeastMCMC.m_nThreads; i++) {
                //if (dists.isDirtyCalculation()) {    //todo
                start = batchSize*i;
                end = batchSize*(i + 1) - 1;
                if(i == BeastMCMC.m_nThreads -1){
                    end += tempTreeLiks.size()% BeastMCMC.m_nThreads;
                }

                CoreRunnable coreRunnable = new CoreRunnable(start, end,site);

                BeastMCMC.g_exec.execute(coreRunnable);
                //}

            }

            m_nCountDown.await();
            double logPs[] = new double[tempTreeLiks.size()];
            for (int i = 0; i < logPs.length; i++) {
                logPs[i] = tempTreeLiks.get(i).getCurrentLogP();
            }
            return logPs;
        } catch (RejectedExecutionException e) {
            useThreads = false;
            System.err.println("Stop using threads: " + e.getMessage());
            // refresh thread pool
            BeastMCMC.g_exec = Executors.newFixedThreadPool(BeastMCMC.m_nThreads);
            return calculateLogP(site);
        }
    }


    public double[] calculateLogP(int site) throws Exception{
        double[] logPs = new double[tempTreeLiks.size()];

        if (useThreads) {
            logPs = calculateLogPUsingThreads(site);
        } else {
            for (int i = 0; i < logPs.length; i++) {
                //if (tempTreeLiks.get(i).isDirtyCalculation()) {
                    logPs[i] = tempTreeLiks.get(i).calculateLogP(site);
                //} else {
                //    logP += tempTreeLiks.get(i).getCurrentLogP();
                //}
                if (Double.isInfinite(logPs[i]) || Double.isNaN(logPs[i])) {
                    return logPs;
                }
            }
        }
        return logPs;

    }


    public double[] calculateLogP(
            RealParameter[] modelParameterList,
            RealParameter[] modelCodeList,
            RealParameter[] freqsList,
            int site){

        double[] logPs;
        try{
            for(int i = 0;i < tempTreeLiks.size(); i++){
                tempTreeLiks.get(i).setSubstModelParameter(
                        modelParameterList[i],
                        modelCodeList[i],
                        freqsList[i]);
            }

            logPs =  calculateLogP(site);
        }catch(Exception e){
            throw new RuntimeException(e);
        }
        return logPs;
    }


    public double[] calculateLogP(
            QuietRealParameter[] modelParameterList,
            QuietRealParameter[] modelCodeList,
            QuietRealParameter[] freqsList,
            QuietRealParameter[] rateList,
            int site){
        double[] logPs;
        try{
            for(int i = 0;i < tempTreeLiks.size(); i++){
                tempTreeLiks.get(i).setSubstModelParameter(
                        modelParameterList[i],
                        modelCodeList[i],
                        freqsList[i]);
                tempTreeLiks.get(i).setRateParameter(rateList[i]);
            }


            logPs =  calculateLogP(site);
        }catch(Exception e){
            throw new RuntimeException(e);

        }
        return logPs;
    }


    public double[] calculateLogP(
            RealParameter[] rateList,
            int site){
        double[] logPs;
        try{
            for(int i = 0;i < tempTreeLiks.size(); i++){
                tempTreeLiks.get(i).setRateParameter(rateList[i]);
            }
            logPs =  calculateLogP(site);
        }catch(Exception e){
            throw new RuntimeException(e);

        }
        return logPs;
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
