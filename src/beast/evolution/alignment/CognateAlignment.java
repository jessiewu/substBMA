package beast.evolution.alignment;

import beast.core.Description;
import beast.core.Input;
import beast.evolution.datatype.DataType;
import beast.evolution.datatype.LanguageIntegerData;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * @author Chieh-Hsi Wu
 */
@Description("This class implements the alignment for cognate data.")
public class CognateAlignment extends Alignment {
    public Input<LanguageIntegerData> dataTypeInput = new Input<LanguageIntegerData>(
            "languageIntegerDataType",
            "The language integer data type congnates.",
            Input.Validate.REQUIRED
    );

    public Input<TaxonSet> taxaInput = new Input<TaxonSet>(
            "taxa",
            "The taxa in the data set.",
            Input.Validate.REQUIRED
    );

    public Input<List<Cognate>> cognatesInput = new Input<List<Cognate>>(
            "cognate",
            "The cognate of the word of interest",
            new ArrayList<Cognate>(),
            Input.Validate.REQUIRED
    );

    public CognateAlignment(){
        m_pSequences.setRule(Input.Validate.OPTIONAL);
        m_sDataType.setRule(Input.Validate.OPTIONAL);
    }


    private List<Cognate> cognates;
    public void initAndValidate(){
        m_dataType = dataTypeInput.get();
        m_sTaxaNames = taxaInput.get().asStringList();
        m_nMaxStateCount = m_dataType.getStateCount();
        cognates = cognatesInput.get();

        int firstWordCount = cognates.get(0).getWordCount();
        for(int i = 1; i < cognates.size(); i++){
            if(cognates.get(i).getWordCount() != firstWordCount){
                throw new RuntimeException("The number of words/taxa should be the same for all cognates. " +
                        "The first cognate has "+firstWordCount+" words, while the "+(i+1)+" cognate has "+cognates.get(i).getWordCount());
            }
        }

        calcPatterns();




    }

    /**
     * calculate patterns from sequence data
     * *
     */
    protected void calcPatterns() {
        int nTaxa = cognates.get(0).getWordCount();
        int nSites = cognates.size();

        // convert data to transposed int array
        int[][] nData = new int[nSites][nTaxa];
        for (int i = 0; i < nSites; i++) {
            nData[i] = cognates.get(i).getWords(m_dataType);
        }

        // sort data
        SiteComparator comparator = new SiteComparator();
        Arrays.sort(nData, comparator);

        // count patterns in sorted data
        int nPatterns = 1;
        int[] weights = new int[nSites];
        weights[0] = 1;
        for (int i = 1; i < nSites; i++) {
            if (comparator.compare(nData[i - 1], nData[i]) != 0) {
                nPatterns++;
                nData[nPatterns - 1] = nData[i];
            }
            weights[nPatterns - 1]++;
        }

        // reserve memory for patterns
        m_nWeight = new int[nPatterns];
        m_nPatterns = new int[nPatterns][nTaxa];
        for (int i = 0; i < nPatterns; i++) {
            m_nWeight[i] = weights[i];
            m_nPatterns[i] = nData[i];
        }

        // find patterns for the sites
        m_nPatternIndex = new int[nSites];
        for (int i = 0; i < nSites; i++) {
        	int [] sites = cognates.get(i).getWords(m_dataType);
            m_nPatternIndex[i] = Arrays.binarySearch(m_nPatterns, sites, comparator);
        }

        // determine maximum state count
        // Usually, the state count is equal for all sites,
        // though for SnAP analysis, this is typically not the case.
        m_nMaxStateCount = m_dataType.getStateCount();
        // report some statistics

        System.out.println(getNrTaxa() + " taxa");
        System.out.println(getSiteCount() + " sites");
        System.out.println(getPatternCount() + " patterns");


        if (m_bStripInvariantSites.get()) {
            // don't add patterns that are invariant, e.g. all gaps
        	System.err.print("Stripping invariant sites");
	        for (int i = 0; i < nPatterns; i++) {
	        	int [] nPattern = m_nPatterns[i];
	        	int iValue = nPattern[0];
	        	boolean bIsInvariant = true;
	        	for (int k = 1; k < nPattern.length; k++) {
	        		if (nPattern[k] != iValue) {
	        			bIsInvariant = false;
	        			break;
	        		}
	        	}
	        	if (bIsInvariant) {
	                m_nWeight[i] = 0;
	                System.err.print(" <" + iValue+"> ");
	        	}
	        }
	        System.err.println();
        }


    } // calcPatterns


}
