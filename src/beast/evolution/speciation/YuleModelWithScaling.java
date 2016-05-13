package beast.evolution.speciation;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.Scaler;
import beast.evolution.tree.Tree;

/**
 * @author Chieh-Hsi
 */
@Description("This yule model takes an extra parameter so that node heights are scaled correctly. " +
        "This is used when we cannot estimate time and substitution rate(s) separately.")
public class YuleModelWithScaling extends YuleModel {
    public Input<Scaler> scalerInput = new Input<Scaler>("scaler", "Some thing that scales the tree branches");
    public Input<RealParameter> scalerParameterInput = new Input<RealParameter>("scalerParameter", "Some parameter that scales the tree branches", Input.Validate.XOR, scalerInput);

    private Scaler scaler;
    private RealParameter scalerParameter;
    public void initAndValidate(){
        scaler = scalerInput.get();
        scalerParameter = scalerParameterInput.get();
    }

    protected double calculateTreeLogLikelihood(final Tree tree, final double rho, final double a) {
        final int taxonCount = tree.getLeafNodeCount();
        final double r = birthDiffRateParameterInput.get().getValue();

        double logL = logTreeProbability(taxonCount, r, rho, a);
        double scaleFactor;
        if(scalerParameter == null){
            scaleFactor = scaler.getScaleFactor();
        }else{
            scaleFactor = scalerParameter.getValue();
        }
        final Node [] nodes = tree.getNodesAsArray();
        for (int i = taxonCount; i < nodes.length; i++) {
            assert ( ! nodes[i].isLeaf() );
            logL += calcLogNodeProbability(nodes[i], r, rho, a, taxonCount,scaleFactor);
        }

        return logL;
    }
    protected double calcLogNodeProbability(Node node, double r, double rho, double a, int taxonCount, double scaleFactor) {
        final double height = node.getHeight()*scaleFactor;
        final double mrh = -r * height;

        if( ! conditionalOnRoot ) {
            final double z = Math.log(rho + ((1 - rho) - a) * Math.exp(mrh));
            double l = -2 * z + mrh;

            if(node.isRoot()) {
                l += mrh - z;
            }
            return l;
        } else {
            double l;
            if( !node.isRoot() ) {
                final double z = Math.log(1 - a * Math.exp(mrh));
                l = -2 * z + mrh;
            } else {
                // Root dependent coefficient from each internal node
                final double ca = 1 - a;
                final double emrh = Math.exp(-mrh);
                if( emrh != 1.0 ) {
                  l = (taxonCount - 2) * Math.log(r * ca * (1 + ca /(emrh - 1)));
                } else {  // use exp(x)-1 = x for x near 0
                  l = (taxonCount - 2) * Math.log(ca * (r + ca/height));
                }
            }
            return l;
        }
    }

    public boolean requiresRecalculation(){
        if(scalerParameter == null){
            return super.requiresRecalculation()||scaler.isDirtyCalculation();
        }else{
            return super.requiresRecalculation()||scalerParameter.somethingIsDirty();
        }
    }
}
