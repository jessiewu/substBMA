package beast.evolution.likelihood;

/**
 * @author Chieh-Hsi Wu
 */
public class WVLikelihoodCore extends BeerLikelihoodCore{
    protected boolean[] unmasked;
    private boolean[] storedUnmasked;
    private double m_fScalingThreshold = 1.0E-100;
    public WVLikelihoodCore(int nStateCount, boolean[] unmasked) {
		super(nStateCount);
        this.unmasked = unmasked;
        storedUnmasked = new boolean[unmasked.length];
	}

    public void setUnmasked(int patID, boolean unmasked){
        this.unmasked[patID] = unmasked;
    }

    public void setUnmasked(boolean[] unmasked){
        this.unmasked = unmasked;
    }

    public boolean[] getUnmasked(){
        return unmasked;
    }

    public boolean getUnmasked(int i){
        return unmasked[i];
    }

    public void initialize(int nNodeCount, int nPatternCount, int nMatrixCount, boolean bIntegrateCategories, boolean useAmbiguities){
        super.initialize(nNodeCount, nPatternCount, nMatrixCount, bIntegrateCategories, useAmbiguities);

        if(unmasked.length != m_nPatterns){
            throw new RuntimeException("The length of the unmasked ("+unmasked.length+
                    ") array needs to be the same as the number of patterns("+m_nPatterns+").");
        }
    }

	/**
	 * Calculates partial likelihoods at a node when both children have states.
	 */
	protected void calculateStatesStatesPruning(int[] iStates1, double[] fMatrices1,
												int[] iStates2, double[] fMatrices2,
												double[] fPartials3){
		int v = 0;

		for (int l = 0; l < m_nMatrices; l++) {

			for (int k = 0; k < m_nPatterns; k++) {
                if(unmasked[k]){
				    int state1 = iStates1[k];
				    int state2 = iStates2[k];

				    int w = l * m_nMatrixSize;

                    if (state1 < m_nStates && state2 < m_nStates) {

				    	for (int i = 0; i < m_nStates; i++) {

				    		fPartials3[v] = fMatrices1[w + state1] * fMatrices2[w + state2];

				    		v++;
				    		w += m_nStates;
				    	}

				    } else if (state1 < m_nStates) {
				    	// child 2 has a gap or unknown state so treat it as unknown

				    	for (int i = 0; i < m_nStates; i++) {

				    		fPartials3[v] = fMatrices1[w + state1];

				    		v++;
				    		w += m_nStates;
				    	}
				    } else if (state2 < m_nStates) {
				    	// child 2 has a gap or unknown state so treat it as unknown

				    	for (int i = 0; i < m_nStates; i++) {

				    		fPartials3[v] = fMatrices2[w + state2];

				    		v++;
				    		w += m_nStates;
				    	}
				    } else {
				    	// both children have a gap or unknown state so set partials to 1

				    	for (int j = 0; j < m_nStates; j++) {
				    		fPartials3[v] = 1.0;
				    		v++;
				    	}
				    }
                }else{
                    v += m_nStates;
                }
			}
		}
	}

    /**
	 * Calculates partial likelihoods at a node when one child has states and one has partials.
	 */
	protected void calculateStatesPartialsPruning(	int[] iStates1, double[] fMatrices1,
													double[] fPartials2, double[] fMatrices2,
													double[] fPartials3){

		double sum, tmp;

		int u = 0;
		int v = 0;

		for (int l = 0; l < m_nMatrices; l++) {
			for (int k = 0; k < m_nPatterns; k++) {
                if(unmasked[k]){
				    int state1 = iStates1[k];

                    int w = l * m_nMatrixSize;

				    if (state1 < m_nStates) {


				    	for (int i = 0; i < m_nStates; i++) {

				    		tmp = fMatrices1[w + state1];

				    		sum = 0.0;
				    		for (int j = 0; j < m_nStates; j++) {
				    			sum += fMatrices2[w] * fPartials2[v + j];
				    			w++;
				    		}

				    		fPartials3[u] = tmp * sum;
				    		u++;
				    	}

				    	v += m_nStates;
				    } else {
				    	// Child 1 has a gap or unknown state so don't use it

				    	for (int i = 0; i < m_nStates; i++) {

				    		sum = 0.0;
				    		for (int j = 0; j < m_nStates; j++) {
				    			sum += fMatrices2[w] * fPartials2[v + j];
				    			w++;
				    		}

				    		fPartials3[u] = sum;
				    		u++;
				    	}

				    	v += m_nStates;
				    }
                }else{
                    u += m_nStates;
                    v += m_nStates;
                }
			}
		}
	}

	/**
	 * Calculates partial likelihoods at a node when both children have partials.
	 */
	protected void calculatePartialsPartialsPruning(double[] fPartials1, double[] fMatrices1,
													double[] fPartials2, double[] fMatrices2,
													double[] fPartials3)
	{
		double sum1, sum2;

		int u = 0;
		int v = 0;

		for (int l = 0; l < m_nMatrices; l++) {

			for (int k = 0; k < m_nPatterns; k++) {
                if(unmasked[k]){
                    int w = l * m_nMatrixSize;

				    for (int i = 0; i < m_nStates; i++) {

				    	sum1 = sum2 = 0.0;

				    	for (int j = 0; j < m_nStates; j++) {
				    		sum1 += fMatrices1[w] * fPartials1[v + j];
				    		sum2 += fMatrices2[w] * fPartials2[v + j];

				    		w++;
				    	}

				    	fPartials3[u] = sum1 * sum2;
				    	u++;
				    }
				    v += m_nStates;
                }else{
                    v += m_nStates;
                    u += m_nStates;
                }
			}
		}
	}



	/**
	 * Calculates partial likelihoods at a node when both children have states.
	 */
	protected void calculateStatesStatesPruning(int[] iStates1, double[] fMatrices1,
												int[] iStates2, double[] fMatrices2,
												double[] fPartials3, int[] iMatrixMap){
		int v = 0;

		for (int k = 0; k < m_nPatterns; k++) {
            if(unmasked[k]){
			    int state1 = iStates1[k];
			    int state2 = iStates2[k];

			    int w = iMatrixMap[k] * m_nMatrixSize;

			    if (state1 < m_nStates && state2 < m_nStates) {

			    	for (int i = 0; i < m_nStates; i++) {

			    		fPartials3[v] = fMatrices1[w + state1] * fMatrices2[w + state2];

			    		v++;
			    		w += m_nStates;
			    	}

			    } else if (state1 < m_nStates) {
			    	// child 2 has a gap or unknown state so treat it as unknown

			    	for (int i = 0; i < m_nStates; i++) {

			    		fPartials3[v] = fMatrices1[w + state1];

			    		v++;
			    		w += m_nStates;
			    	}
			    } else if (state2 < m_nStates) {
			    	// child 2 has a gap or unknown state so treat it as unknown

			    	for (int i = 0; i < m_nStates; i++) {

			    		fPartials3[v] = fMatrices2[w + state2];

			    		v++;
			    		w += m_nStates;
			    	}
			    } else {
			    	// both children have a gap or unknown state so set partials to 1

			    	for (int j = 0; j < m_nStates; j++) {
			    		fPartials3[v] = 1.0;
			    		v++;
			    	}
			    }
            }else{
                v += m_nStates;
            }
		}
	}

	/**
	 * Calculates partial likelihoods at a node when one child has states and one has partials.
	 */
	protected void calculateStatesPartialsPruning(	int[] iStates1, double[] fMatrices1,
													double[] fPartials2, double[] fMatrices2,
													double[] fPartials3, int[] iMatrixMap){

		double sum, tmp;

		int u = 0;
		int v = 0;

		for (int k = 0; k < m_nPatterns; k++) {
            if(unmasked[k]){
			    int state1 = iStates1[k];

			    int w = iMatrixMap[k] * m_nMatrixSize;

			    if (state1 < m_nStates) {

			    	for (int i = 0; i < m_nStates; i++) {

			    		tmp = fMatrices1[w + state1];

			    		sum = 0.0;
			    		for (int j = 0; j < m_nStates; j++) {
			    			sum += fMatrices2[w] * fPartials2[v + j];
			    			w++;
			    		}

			    		fPartials3[u] = tmp * sum;
			    		u++;
			    	}

			    	v += m_nStates;
			    } else {
			    	// Child 1 has a gap or unknown state so don't use it

			    	for (int i = 0; i < m_nStates; i++) {

			    		sum = 0.0;
			    		for (int j = 0; j < m_nStates; j++) {
			    			sum += fMatrices2[w] * fPartials2[v + j];
			    			w++;
			    		}

			    		fPartials3[u] = sum;
			    		u++;
			    	}

			    	v += m_nStates;
			    }
            }else{
                u += m_nStates;
                v += m_nStates;

            }
		}
	}


	/**
	 * Calculates partial likelihoods at a node when both children have partials.
	 */
	protected void calculatePartialsPartialsPruning(double[] fPartials1, double[] fMatrices1,
													double[] fPartials2, double[] fMatrices2,
													double[] fPartials3, int[] iMatrixMap){
        //super.calculatePartialsPartialsPruning(fPartials1,fMatrices1,fPartials2,fMatrices2,fPartials3,iMatrixMap);
		double sum1, sum2;

		int u = 0;
		int v = 0;

		for (int k = 0; k < m_nPatterns; k++) {
            if(unmasked[k]){
			    int w = iMatrixMap[k] * m_nMatrixSize;

			    for (int i = 0; i < m_nStates; i++) {

			    	sum1 = sum2 = 0.0;

			    	for (int j = 0; j < m_nStates; j++) {
			    		sum1 += fMatrices1[w] * fPartials1[v + j];
			    		sum2 += fMatrices2[w] * fPartials2[v + j];

			    		w++;
			    	}

			    	fPartials3[u] = sum1 * sum2;
			    	u++;
			    }
			    v += m_nStates;
            }else{
                u += m_nStates;
                v += m_nStates;

            }
		}
	}

	/**
	 * Integrates partials across categories.
     * @param fInPartials the array of partials to be integrated
	 * @param fProportions the proportions of sites in each category
	 * @param fOutPartials an array into which the partials will go
	 */
	protected void calculateIntegratePartials(double[] fInPartials, double[] fProportions, double[] fOutPartials)
	{

		int u = 0;
		int v = 0;
		for (int k = 0; k < m_nPatterns; k++) {
            if(unmasked[k]){
			    for (int i = 0; i < m_nStates; i++) {

			    	fOutPartials[u] = fInPartials[v] * fProportions[0];

                    u++;
			    	v++;
			    }
            }else{
                u+=m_nStates;
                v+=m_nStates;
            }
		}


		for (int l = 1; l < m_nMatrices; l++) {
			u = 0;

			for (int k = 0; k < m_nPatterns; k++) {
                if(unmasked[k]){
				    for (int i = 0; i < m_nStates; i++) {

				    	fOutPartials[u] += fInPartials[v] * fProportions[l];
                        u++;
				    	v++;
				    }
                }else{
                    u += m_nStates;
                    v += m_nStates;
                }
			}
		}
	}

	/**
	 * Calculates pattern log likelihoods at a node.
	 * @param fPartials the partials used to calculate the likelihoods
	 * @param fFrequencies an array of state frequencies
	 * @param fOutLogLikelihoods an array into which the likelihoods will go
	 */
	public void calculateLogLikelihoods(double[] fPartials, double[] fFrequencies, double[] fOutLogLikelihoods){
        int v = 0;
		for (int k = 0; k < m_nPatterns; k++) {
            if(unmasked[k]){
                double sum = 0.0;
			    for (int i = 0; i < m_nStates; i++) {

			    	sum += fFrequencies[i] * fPartials[v];
                    //System.out.println("sum: "+sum+" log sum: "+Math.log(sum)+" fPartials[v]: "+fPartials[v]+" log fPartials[v]: "+Math.log(fPartials[v]));
			    	v++;
			    }
                /*System.out.println("calculateLogLikelihoods, sum: "+sum +" scaling factor: "+ getLogScalingFactor(k));

                    System.out.println("fPartials.length: "+fPartials.length);
                    for(int i = 0; i < m_fPartials.length;i++){
                        for(int j = 0; j < m_fPartials[i].length; j++){
                            if(m_fPartials[i][j] != null){
                                for(int l = 0; l < m_fPartials[i][j].length;l++){
                                    System.out.println(m_fPartials[i][j][l]+" ");
                                }
                            }

                        //System.out.print("partials "+i+" :"+fPartials[i]+" ");
                        }
                    }
                System.out.println();
                if(sum == 0){

                    throw new RuntimeException("");
                }*/
                fOutLogLikelihoods[k] = Math.log(sum) + getLogScalingFactor(k);

            }else{
                v += m_nStates;

            }
		}
	}


    /**
     * Scale the partials at a given node. This uses a scaling suggested by Ziheng Yang in
     * Yang (2000) J. Mol. Evol. 51: 423-432
     * <p/>
     * This function looks over the partial likelihoods for each state at each pattern
     * and finds the largest. If this is less than the scalingThreshold (currently set
     * to 1E-40) then it rescales the partials for that pattern by dividing by this number
     * (i.e., normalizing to between 0, 1). It then stores the log of this scaling.
     * This is called for every internal node after the partials are calculated so provides
     * most of the performance hit. Ziheng suggests only doing this on a proportion of nodes
     * but this sounded like a headache to organize (and he doesn't use the threshold idea
     * which improves the performance quite a bit).
     *
     * @param iNodeIndex
     */
    protected void scalePartials(int iNodeIndex) {
        int u = 0;

        for (int i = 0; i < m_nPatterns; i++) {
            if(unmasked[i]){
                double scaleFactor = 0.0;
                int v = u;
                for (int k = 0; k < m_nMatrices; k++) {
                    for (int j = 0; j < m_nStates; j++) {
                        if (m_fPartials[m_iCurrentPartials[iNodeIndex]][iNodeIndex][v] > scaleFactor) {
                            scaleFactor = m_fPartials[m_iCurrentPartials[iNodeIndex]][iNodeIndex][v];
                        }
                        v++;
                    }
                    v += (m_nPatterns - 1) * m_nStates;
                }
                //System.out.println("pattern "+i+" scale factor: "+scaleFactor);
                if (scaleFactor < m_fScalingThreshold && scaleFactor > 0.0) {

                    v = u;
                    for (int k = 0; k < m_nMatrices; k++) {
                        for (int j = 0; j < m_nStates; j++) {
                            m_fPartials[m_iCurrentPartials[iNodeIndex]][iNodeIndex][v] /= scaleFactor;
                            v++;
                        }
                        v += (m_nPatterns - 1) * m_nStates;
                    }
                    m_fScalingFactors[m_iCurrentPartials[iNodeIndex]][iNodeIndex][i] = Math.log(scaleFactor);

                } else {
                    m_fScalingFactors[m_iCurrentPartials[iNodeIndex]][iNodeIndex][i] = 0.0;
                }

            }
            u += m_nStates;

        }
    }

/*	@Override
    public void calcRootPsuedoRootPartials(double[] fFrequencies, int iNode, double [] fPseudoPartials) {
		int u = 0;
		double [] fInPartials = m_fPartials[m_iCurrentPartials[iNode]][iNode];
		for (int k = 0; k < m_nPatterns; k++) {
            if(unmasked[k]){
			    for (int l = 0; l < m_nMatrices; l++) {
			    	for (int i = 0; i < m_nStates; i++) {
			    		fPseudoPartials[u] = fInPartials[u] * fFrequencies[i];
			    		u++;
			    	}
			    }
            }else{
                u+=m_nMatrices*m_nStates;
            }
		}
    }


	@Override
    public void calcNodePsuedoRootPartials(double[] fInPseudoPartials, int iNode, double [] fOutPseudoPartials) {
		double [] fPartials = m_fPartials[m_iCurrentPartials[iNode]][iNode];
		double [] fOldPartials = m_fPartials[m_iStoredPartials[iNode]][iNode];
		//int nMaxK = m_nPatterns * m_nMatrices * m_nStates;
        int nMaxJ = m_nMatrices * m_nStates;
        int k = 0;
        for(int i = 0; i < m_nPatterns; i++){
            if(unmasked[i]){
		        for (int j = 0; j < nMaxJ; j++) {
			        fOutPseudoPartials[k] = fInPseudoPartials[k] * fPartials[k] / fOldPartials[k];
                    k++;
		        }
            }
        }
	}


	@Override
    public void calcPsuedoRootPartials(double [] fParentPseudoPartials, int iNode, double [] fPseudoPartials) {
		int v = 0;
		int u = 0;
		double [] fMatrices = m_fMatrices[m_iCurrentMatrices[iNode]][iNode];
		for (int k = 0; k < m_nPatterns; k++) {
            if(unmasked[k]){
			    for (int l = 0; l < m_nMatrices; l++) {
			    	for (int i = 0; i < m_nStates; i++) {
			    		int w = 0;
			    		double fSum = 0;
			    		for (int j = 0; j < m_nStates; j++) {
			    		      fSum += fParentPseudoPartials[u+j] * fMatrices[w + i];
			    		      w+=m_nStates;
			    		}
			    		fPseudoPartials[v] = fSum;
			    		v++;
//			    		int w = l * m_nMatrixSize;
//			    		double fSum = 0;
//			    		for (int j = 0; j < m_nStates; j++) {
//			    		      fSum += fParentPseudoPartials[u+j] * fMatrices[w+j];
//			    		}
//			    		fPseudoPartials[v] = fSum;
//			    		v++;
			    	}
			    	u += m_nStates;
			    }
            }
		}
    }


    @Override
    void integratePartialsP(double [] fInPartials, double [] fProportions, double [] m_fRootPartials) {
		int nMaxK = m_nPatterns * m_nStates;
        int k = 0;
        for(int i = 0; i < m_nPatterns;i++){
            for(int j = 0; j < m_nStates; j++){
                m_fRootPartials[k] = fInPartials[k] * fProportions[0];
                k++;
            }
        }


        for(int l = 1; l < m_nMatrices; l++){
            int n = nMaxK * l;
            k=0;
            for(int i = 0; i < m_nPatterns;i++){
                if(unmasked[i]){
                    for(int j = 0; j < m_nStates; j++){
                        m_fRootPartials[k] += fInPartials[n+k] * fProportions[l];
                    }
                }
            }
        }

    } // integratePartials */

	/**
	 * Calculates pattern log likelihoods at a node.
	 * @param fPartials the partials used to calculate the likelihoods
	 * @param fOutLogLikelihoods an array into which the likelihoods will go
	 */
    /*@Override
	public void calculateLogLikelihoodsP(double[] fPartials,double[] fOutLogLikelihoods)
	{
        int v = 0;
		for (int k = 0; k < m_nPatterns; k++) {
            if(unmasked[k]){
                double sum = 0.0;
			    for (int i = 0; i < m_nStates; i++) {
				    sum += fPartials[v];
				    v++;
			    }
                fOutLogLikelihoods[k] = Math.log(sum) + getLogScalingFactor(k);
            }else{
                v += m_nStates;
            }
		}
	} */

    public void store(){
        System.arraycopy(unmasked, 0, storedUnmasked, 0, unmasked.length);
        super.store();
    }

    public void restore(){
        boolean[] tmp = unmasked;
        unmasked = storedUnmasked;
        storedUnmasked = tmp;
        super.restore();
    }

    public void restoreUnmaskedOnly(){
        //System.err.println("RestoreUnmasked: "+storedUnmasked[39]);
        boolean[] tmp = unmasked;
        unmasked = storedUnmasked;
        storedUnmasked = tmp;
        //System.err.println("RestoreUnmasked: "+unmasked[39]);
    }

}
