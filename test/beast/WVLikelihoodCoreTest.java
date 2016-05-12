package beast;

import beast.core.MCMCNodeFactory;
import beast.evolution.alignment.Sequence;
import beast.util.TreeParser;
import junit.framework.TestCase;
import beast.evolution.alignment.Alignment;
import beast.evolution.tree.Tree;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.HKY;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.likelihood.NewWVTreeLikelihood;
import beast.evolution.likelihood.TreeLikelihood;
import beast.core.parameter.RealParameter;


/**
 * @author Chieh-Hsi Wu
 */
public class WVLikelihoodCoreTest extends TestCase {

    public static SiteModel getDefaultSiteModel() throws Exception{
        RealParameter frequencies = new RealParameter(new Double[]{0.25, 0.25, 0.25, 0.25});
        Frequencies freqs = new Frequencies();
        freqs.initByName(
                "frequencies", frequencies,
                "estimate", false
        );


        HKY hky = new HKY();
        hky.initByName(
                    "kappa", "1.0",
                    "frequencies",freqs
            );
        SiteModel siteModel = new SiteModel();
        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", 1, "substModel", hky);

        return siteModel;

    }

    /*
     * Sanity checking!
     */
    public void testWVLikelihoodCore0(){
        try{
            Alignment data = getAlignment();
            Tree tree = getTree(data);
            SiteModel siteModel = getDefaultSiteModel();

            int[] patternWeights = data.getWeights();
            NewWVTreeLikelihood newTreeLik = new NewWVTreeLikelihood(patternWeights);
            newTreeLik.initByName(
                    "data", data,
                    "tree", tree,
                    "siteModel", siteModel
            );

            TreeLikelihood treeLik = new TreeLikelihood();
            treeLik.initByName(
                    "data", data,
                    "tree", tree,
                    "siteModel", siteModel
            );
        }catch(Exception e){
            throw new RuntimeException(e);
        }
    }

    /* AAA
     * AAA
     * AAA
     * AAA
     * AAC
     * AGA
     * -1.7210602650853448 221
     * -4.811580369169091 8
     * -5.4467194806410095 1
     */
    public void testWVLikelihoodCore1(){
        try{

            Alignment data = getAlignment();
            Tree tree = getTree(data);


		    SiteModel siteModel = getDefaultSiteModel();

            int[] weights1 = new int[data.getPatternCount()];
            weights1[0] = 1;
            weights1[1] = 1;

            NewWVTreeLikelihood newTreeLik = new NewWVTreeLikelihood(weights1);
            newTreeLik.initByName(
                    "data", data,
                    "tree", tree,
                    "siteModel", siteModel
            );

            assertEquals(newTreeLik.calculateLogP(),-6.532640634254436,1e-10);

            newTreeLik.store();
            newTreeLik.addWeight(2, 1);
            MCMCNodeFactory.checkDirtiness(newTreeLik);
            assertEquals(newTreeLik.calculateLogP(),-11.97936011489545,1e-10);

            tree.setEverythingDirty(true);
            MCMCNodeFactory.checkDirtiness(newTreeLik);
            assertEquals(newTreeLik.calculateLogP(),-11.97936011489545,1e-10);
        }catch(Exception e){
            throw new RuntimeException(e);
        }
    }

    /* AAA
     * AAA
     * AAA
     * AAA
     * AAC
     * AGA
     * -1.7210602650853448 221
     * -4.811580369169091 8
     * -5.4467194806410095 1
     */
    public void testWVLikelihoodCore2(){
        try{

            Alignment data = getAlignment();
            Tree tree = getTree(data);


		    SiteModel siteModel = getDefaultSiteModel();

            int[] weights1 = new int[data.getPatternCount()];
            weights1[0] = 1;
            weights1[1] = 1;

            NewWVTreeLikelihood newTreeLik = new NewWVTreeLikelihood(weights1);
            newTreeLik.initByName(
                    "data", data,
                    "tree", tree,
                    "siteModel", siteModel
            );

            assertEquals(newTreeLik.calculateLogP(),-6.532640634254436,1e-10);

            newTreeLik.store();
            newTreeLik.addWeight(2, 1);
            MCMCNodeFactory.checkDirtiness(newTreeLik);
            assertEquals(newTreeLik.calculateLogP(),-11.97936011489545,1e-10);
            newTreeLik.restore();


            tree.setEverythingDirty(true);
            MCMCNodeFactory.checkDirtiness(newTreeLik);
            assertEquals(newTreeLik.calculateLogP(),-6.532640634254436,1e-10);
            System.out.println(newTreeLik.calculateLogP());
        }catch(Exception e){
            throw new RuntimeException(e);
        }
    }

    /* AAA
     * AAA
     * AAA
     * AAA
     * AAC
     * AGA
     * -1.7210602650853448 221
     * -4.811580369169091 8
     * -5.4467194806410095 1
     */
    public void testWVLikelihoodCore3(){
        try{

            Alignment data = getAlignment();
            Tree tree = getTree(data);


		    SiteModel siteModel = getDefaultSiteModel();

            int[] weights1 = new int[data.getPatternCount()];
            weights1[0] = 1;
            weights1[1] = 1;
            weights1[2] = 1;

            NewWVTreeLikelihood newTreeLik = new NewWVTreeLikelihood(weights1);
            newTreeLik.initByName(
                    "data", data,
                    "tree", tree,
                    "siteModel", siteModel
            );

            assertEquals(newTreeLik.calculateLogP(),-11.97936011489545,1e-10);

            newTreeLik.store();
            newTreeLik.removeWeight(1, 1);
            MCMCNodeFactory.checkDirtiness(newTreeLik);
            assertEquals(newTreeLik.calculateLogP(),-7.167779745726355,1e-10);

            tree.setEverythingDirty(true);
            MCMCNodeFactory.checkDirtiness(newTreeLik);
            assertEquals(newTreeLik.calculateLogP(),-7.167779745726355,1e-10);
            System.out.println(newTreeLik.calculateLogP());
        }catch(Exception e){
            throw new RuntimeException(e);
        }
    }

    /* AAA
     * AAA
     * AAA
     * AAA
     * AAC
     * AGA
     * -1.7210602650853448 221
     * -4.811580369169091 8
     * -5.4467194806410095 1
     */
    public void testWVLikelihoodCore4(){
        try{

            Alignment data = getAlignment();
            Tree tree = getTree(data);


		    SiteModel siteModel = getDefaultSiteModel();

            int[] weights1 = new int[data.getPatternCount()];
            weights1[0] = 1;
            weights1[1] = 1;
            weights1[2] = 1;

            NewWVTreeLikelihood newTreeLik = new NewWVTreeLikelihood(weights1);
            newTreeLik.initByName(
                    "data", data,
                    "tree", tree,
                    "siteModel", siteModel
            );

            assertEquals(newTreeLik.calculateLogP(),-11.97936011489545,1e-10);

            newTreeLik.store();
            newTreeLik.removeWeight(1, 1);
            MCMCNodeFactory.checkDirtiness(newTreeLik);
            assertEquals(newTreeLik.calculateLogP(),-7.167779745726355,1e-10);

            newTreeLik.restore();
            tree.setEverythingDirty(true);
            MCMCNodeFactory.checkDirtiness(newTreeLik);
            assertEquals(newTreeLik.calculateLogP(),-11.97936011489545,1e-10);
            System.out.println(newTreeLik.calculateLogP());
        }catch(Exception e){
            throw new RuntimeException(e);
        }
    }


    /* AAA
     * AAA
     * AAA
     * AAA
     * AAC
     * AGA
     * -1.7210602650853448 221
     * -4.811580369169091 8
     * -5.4467194806410095 1
     */
    public void testWVLikelihoodCore5(){
        try{

            Alignment data = getAlignment();
            Tree tree = getTree(data);


		    SiteModel siteModel = getDefaultSiteModel();

            int[] weights1 = new int[data.getPatternCount()];
            weights1[0] = 1;
            weights1[1] = 1;
            weights1[2] = 1;

            NewWVTreeLikelihood newTreeLik = new NewWVTreeLikelihood(weights1);
            newTreeLik.initByName(
                    "data", data,
                    "tree", tree,
                    "siteModel", siteModel
            );

            assertEquals(newTreeLik.calculateLogP(),-11.97936011489545,1e-10);

            newTreeLik.store();
            newTreeLik.removeWeight(1, 1);
            MCMCNodeFactory.checkDirtiness(newTreeLik);
            assertEquals(newTreeLik.calculateLogP(),-7.167779745726355,1e-10);

            newTreeLik.restore();
            assertEquals(newTreeLik.calculateLogP(),-11.97936011489545,1e-10);

            newTreeLik.store();
            newTreeLik.removeWeight(2, 1);
            MCMCNodeFactory.checkDirtiness(newTreeLik);
            assertEquals(newTreeLik.calculateLogP(),-6.532640634254436,1e-10);

            newTreeLik.store();
            newTreeLik.addWeight(0, 2);
            MCMCNodeFactory.checkDirtiness(newTreeLik);
            assertEquals(newTreeLik.calculateLogP(),-9.974761164425125,1e-10);


            newTreeLik.store();
            newTreeLik.addWeight(1, 1);
            MCMCNodeFactory.checkDirtiness(newTreeLik);
            assertEquals(newTreeLik.calculateLogP(),-14.78634153359422,1e-10);

            newTreeLik.restore();
            assertEquals(newTreeLik.calculateLogP(),-9.974761164425125,1e-10);


            newTreeLik.addWeight(1, 0);
            MCMCNodeFactory.checkDirtiness(newTreeLik);
            assertEquals(newTreeLik.calculateLogP(),-9.974761164425125,1e-10);

            newTreeLik.addWeight(2, 0);
            MCMCNodeFactory.checkDirtiness(newTreeLik);
            assertEquals(newTreeLik.calculateLogP(),-9.974761164425125,1e-10);
        }catch(Exception e){
            throw new RuntimeException(e);
        }
    }
    
    static public Alignment getAlignment() throws Exception {
        Sequence human = new Sequence("human", "AGAAATATGTCTGATAAAAGAGTTACTTTGATAGAGTAAATAATAGGAGCTTAAACCCCCTTATTTCTACTAGGACTATGAGAATCGAACCCATCCCTGAGAATCCAAAATTCTCCGTGCCACCTATCACACCCCATCCTAAGTAAGGTCAGCTAAATAAGCTATCGGGCCCATACCCCGAAAATGTTGGTTATACCCTTCCCGTACTAAGAAATTTAGGTTAAATACAGACCAAGAGCCTTCAAAGCCCTCAGTAAGTTG-CAATACTTAATTTCTGTAAGGACTGCAAAACCCCACTCTGCATCAACTGAACGCAAATCAGCCACTTTAATTAAGCTAAGCCCTTCTAGACCAATGGGACTTAAACCCACAAACACTTAGTTAACAGCTAAGCACCCTAATCAAC-TGGCTTCAATCTAAAGCCCCGGCAGG-TTTGAAGCTGCTTCTTCGAATTTGCAATTCAATATGAAAA-TCACCTCGGAGCTTGGTAAAAAGAGGCCTAACCCCTGTCTTTAGATTTACAGTCCAATGCTTCA-CTCAGCCATTTTACCACAAAAAAGGAAGGAATCGAACCCCCCAAAGCTGGTTTCAAGCCAACCCCATGGCCTCCATGACTTTTTCAAAAGGTATTAGAAAAACCATTTCATAACTTTGTCAAAGTTAAATTATAGGCT-AAATCCTATATATCTTA-CACTGTAAAGCTAACTTAGCATTAACCTTTTAAGTTAAAGATTAAGAGAACCAACACCTCTTTACAGTGA");
        Sequence chimp = new Sequence("chimp", "AGAAATATGTCTGATAAAAGAATTACTTTGATAGAGTAAATAATAGGAGTTCAAATCCCCTTATTTCTACTAGGACTATAAGAATCGAACTCATCCCTGAGAATCCAAAATTCTCCGTGCCACCTATCACACCCCATCCTAAGTAAGGTCAGCTAAATAAGCTATCGGGCCCATACCCCGAAAATGTTGGTTACACCCTTCCCGTACTAAGAAATTTAGGTTAAGCACAGACCAAGAGCCTTCAAAGCCCTCAGCAAGTTA-CAATACTTAATTTCTGTAAGGACTGCAAAACCCCACTCTGCATCAACTGAACGCAAATCAGCCACTTTAATTAAGCTAAGCCCTTCTAGATTAATGGGACTTAAACCCACAAACATTTAGTTAACAGCTAAACACCCTAATCAAC-TGGCTTCAATCTAAAGCCCCGGCAGG-TTTGAAGCTGCTTCTTCGAATTTGCAATTCAATATGAAAA-TCACCTCAGAGCTTGGTAAAAAGAGGCTTAACCCCTGTCTTTAGATTTACAGTCCAATGCTTCA-CTCAGCCATTTTACCACAAAAAAGGAAGGAATCGAACCCCCTAAAGCTGGTTTCAAGCCAACCCCATGACCTCCATGACTTTTTCAAAAGATATTAGAAAAACTATTTCATAACTTTGTCAAAGTTAAATTACAGGTT-AACCCCCGTATATCTTA-CACTGTAAAGCTAACCTAGCATTAACCTTTTAAGTTAAAGATTAAGAGGACCGACACCTCTTTACAGTGA");
        Sequence bonobo = new Sequence("bonobo", "AGAAATATGTCTGATAAAAGAATTACTTTGATAGAGTAAATAATAGGAGTTTAAATCCCCTTATTTCTACTAGGACTATGAGAGTCGAACCCATCCCTGAGAATCCAAAATTCTCCGTGCCACCTATCACACCCCATCCTAAGTAAGGTCAGCTAAATAAGCTATCGGGCCCATACCCCGAAAATGTTGGTTATACCCTTCCCGTACTAAGAAATTTAGGTTAAACACAGACCAAGAGCCTTCAAAGCTCTCAGTAAGTTA-CAATACTTAATTTCTGTAAGGACTGCAAAACCCCACTCTGCATCAACTGAACGCAAATCAGCCACTTTAATTAAGCTAAGCCCTTCTAGATTAATGGGACTTAAACCCACAAACATTTAGTTAACAGCTAAACACCCTAATCAGC-TGGCTTCAATCTAAAGCCCCGGCAGG-TTTGAAGCTGCTTCTTTGAATTTGCAATTCAATATGAAAA-TCACCTCAGAGCTTGGTAAAAAGAGGCTTAACCCCTGTCTTTAGATTTACAGTCCAATGCTTCA-CTCAGCCATTTTACCACAAAAAAGGAAGGAATCGAACCCCCTAAAGCTGGTTTCAAGCCAACCCCATGACCCCCATGACTTTTTCAAAAGATATTAGAAAAACTATTTCATAACTTTGTCAAAGTTAAATTACAGGTT-AAACCCCGTATATCTTA-CACTGTAAAGCTAACCTAGCATTAACCTTTTAAGTTAAAGATTAAGAGGACCAACACCTCTTTACAGTGA");
        Sequence gorilla = new Sequence("gorilla", "AGAAATATGTCTGATAAAAGAGTTACTTTGATAGAGTAAATAATAGAGGTTTAAACCCCCTTATTTCTACTAGGACTATGAGAATTGAACCCATCCCTGAGAATCCAAAATTCTCCGTGCCACCTGTCACACCCCATCCTAAGTAAGGTCAGCTAAATAAGCTATCGGGCCCATACCCCGAAAATGTTGGTCACATCCTTCCCGTACTAAGAAATTTAGGTTAAACATAGACCAAGAGCCTTCAAAGCCCTTAGTAAGTTA-CAACACTTAATTTCTGTAAGGACTGCAAAACCCTACTCTGCATCAACTGAACGCAAATCAGCCACTTTAATTAAGCTAAGCCCTTCTAGATCAATGGGACTCAAACCCACAAACATTTAGTTAACAGCTAAACACCCTAGTCAAC-TGGCTTCAATCTAAAGCCCCGGCAGG-TTTGAAGCTGCTTCTTCGAATTTGCAATTCAATATGAAAT-TCACCTCGGAGCTTGGTAAAAAGAGGCCCAGCCTCTGTCTTTAGATTTACAGTCCAATGCCTTA-CTCAGCCATTTTACCACAAAAAAGGAAGGAATCGAACCCCCCAAAGCTGGTTTCAAGCCAACCCCATGACCTTCATGACTTTTTCAAAAGATATTAGAAAAACTATTTCATAACTTTGTCAAGGTTAAATTACGGGTT-AAACCCCGTATATCTTA-CACTGTAAAGCTAACCTAGCGTTAACCTTTTAAGTTAAAGATTAAGAGTATCGGCACCTCTTTGCAGTGA");
        Sequence orangutan = new Sequence("orangutan", "AGAAATATGTCTGACAAAAGAGTTACTTTGATAGAGTAAAAAATAGAGGTCTAAATCCCCTTATTTCTACTAGGACTATGGGAATTGAACCCACCCCTGAGAATCCAAAATTCTCCGTGCCACCCATCACACCCCATCCTAAGTAAGGTCAGCTAAATAAGCTATCGGGCCCATACCCCGAAAATGTTGGTTACACCCTTCCCGTACTAAGAAATTTAGGTTA--CACAGACCAAGAGCCTTCAAAGCCCTCAGCAAGTCA-CAGCACTTAATTTCTGTAAGGACTGCAAAACCCCACTTTGCATCAACTGAGCGCAAATCAGCCACTTTAATTAAGCTAAGCCCTCCTAGACCGATGGGACTTAAACCCACAAACATTTAGTTAACAGCTAAACACCCTAGTCAAT-TGGCTTCAGTCCAAAGCCCCGGCAGGCCTTAAAGCTGCTCCTTCGAATTTGCAATTCAACATGACAA-TCACCTCAGGGCTTGGTAAAAAGAGGTCTGACCCCTGTTCTTAGATTTACAGCCTAATGCCTTAACTCGGCCATTTTACCGCAAAAAAGGAAGGAATCGAACCTCCTAAAGCTGGTTTCAAGCCAACCCCATAACCCCCATGACTTTTTCAAAAGGTACTAGAAAAACCATTTCGTAACTTTGTCAAAGTTAAATTACAGGTC-AGACCCTGTGTATCTTA-CATTGCAAAGCTAACCTAGCATTAACCTTTTAAGTTAAAGACTAAGAGAACCAGCCTCTCTTTGCAATGA");
        Sequence siamang = new Sequence("siamang", "AGAAATACGTCTGACGAAAGAGTTACTTTGATAGAGTAAATAACAGGGGTTTAAATCCCCTTATTTCTACTAGAACCATAGGAGTCGAACCCATCCTTGAGAATCCAAAACTCTCCGTGCCACCCGTCGCACCCTGTTCTAAGTAAGGTCAGCTAAATAAGCTATCGGGCCCATACCCCGAAAATGTTGGTTATACCCTTCCCATACTAAGAAATTTAGGTTAAACACAGACCAAGAGCCTTCAAAGCCCTCAGTAAGTTAACAAAACTTAATTTCTGCAAGGGCTGCAAAACCCTACTTTGCATCAACCGAACGCAAATCAGCCACTTTAATTAAGCTAAGCCCTTCTAGATCGATGGGACTTAAACCCATAAAAATTTAGTTAACAGCTAAACACCCTAAACAACCTGGCTTCAATCTAAAGCCCCGGCAGA-GTTGAAGCTGCTTCTTTGAACTTGCAATTCAACGTGAAAAATCACTTCGGAGCTTGGCAAAAAGAGGTTTCACCTCTGTCCTTAGATTTACAGTCTAATGCTTTA-CTCAGCCACTTTACCACAAAAAAGGAAGGAATCGAACCCTCTAAAACCGGTTTCAAGCCAGCCCCATAACCTTTATGACTTTTTCAAAAGATATTAGAAAAACTATTTCATAACTTTGTCAAAGTTAAATCACAGGTCCAAACCCCGTATATCTTATCACTGTAGAGCTAGACCAGCATTAACCTTTTAAGTTAAAGACTAAGAGAACTACCGCCTCTTTACAGTGA");

        Alignment data = new Alignment();
        data.initByName("sequence", human, "sequence", chimp, "sequence", bonobo, "sequence", gorilla, "sequence", orangutan, "sequence", siamang,
                "dataType", "nucleotide"
        );
        return data;
    }
    
    static public Tree getTree(Alignment data) throws Exception {
        TreeParser tree = new TreeParser();
        tree.initByName("taxa", data,
                "newick", "((((human:0.024003,(chimp:0.010772,bonobo:0.010772):0.013231):0.012035,gorilla:0.036038):0.033087000000000005,orangutan:0.069125):0.030456999999999998,siamang:0.099582);");
        return tree;
    }
}
