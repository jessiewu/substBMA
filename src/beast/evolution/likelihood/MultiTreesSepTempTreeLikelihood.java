package beast.evolution.likelihood;

import beast.core.Input;
import beast.evolution.alignment.Alignment;
import beast.evolution.tree.Tree;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: cwu080
 * Date: 5/05/13
 * Time: 9:37 PM
 * To change this template use File | Settings | File Templates.
 */
public class MultiTreesSepTempTreeLikelihood  extends SepTempTreeLikelihood{
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


    private SepTempTreeLikelihood[] sepTempTreeLikelihoods;
    private List<Alignment> alignments;
    private List<Tree> trees;
    public void initAndValidate(){
        alignments = alignmentsInput.get();
        trees = treesInput.get();

        sepTempTreeLikelihoods = new SepTempTreeLikelihood[alignments.size()];



    }
}
