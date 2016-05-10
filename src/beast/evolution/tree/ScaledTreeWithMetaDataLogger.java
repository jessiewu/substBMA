package beast.evolution.tree;

import beast.core.*;
import beast.core.parameter.Parameter;
import beast.core.parameter.RealParameter;
import beast.evolution.branchratemodel.BranchRateModel;

import java.io.PrintStream;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.List;

/**
 * @author Chieh-Hsi Wu
 */
@Description("Scales the tree with the provided scaler and then log it with metadata.")
public class ScaledTreeWithMetaDataLogger  extends TreeWithMetaDataLogger implements Loggable {
    public Input<Scaler> m_scaler = new Input<Scaler>("scaler","To scale the branch lengths by some number.", Input.Validate.REQUIRED);
    public Input<RealParameter> scalerParameterInput = new Input<RealParameter>("scalerParameter", "Some parameter that scales the tree branches", Input.Validate.XOR, m_scaler);
    Scaler scaler;
    RealParameter scaleParameter;
    private double scaleFactor;

    private DecimalFormat df;

	@Override
	public void initAndValidate() {
		super.initAndValidate();
        scaler = m_scaler.get();
        scaleParameter = scalerParameterInput.get();

        int dp = decimalPlacesInput.get();

        if (dp < 0) {
            df = null;
        } else {
            // just new DecimalFormat("#.######") (with dp time '#' after the decimal)
            df = new DecimalFormat("#."+new String(new char[dp]).replace('\0', '#'));
            df.setRoundingMode(RoundingMode.HALF_UP);
        }
	}

    @Override
	public void log(int nSample, PrintStream out) {
		// make sure we get the current version of the inputs
        Tree tree = (Tree) treeInput.get().getCurrent();
        List<Function> metadata = parameterInput.get();
        for (int i = 0; i<metadata.size(); i++) {
            if (metadata.get(i) instanceof StateNode)
                metadata.set(i, ((StateNode) metadata.get(i)).getCurrent());
        }
        BranchRateModel.Base branchRateModel = clockModelInput.get();
        // write out the log tree with meta data
        out.print("tree STATE_" + nSample + " = ");
		tree.getRoot().sort();
        if(scaler != null){
            scaleFactor = scaler.getScaleFactor();
        }else{
            scaleFactor = scaleParameter.getValue();
        }
		out.print(toNewick(tree.getRoot(), metadata, branchRateModel, scaleFactor));
        //out.print(tree.getRoot().toShortNewick(false));
        out.print(";");
	}


    /**
     * Appends a double to the given StringBuffer, formatting it using
     * the private DecimalFormat instance, if the input 'dp' has been
     * given a non-negative integer, otherwise just uses default
     * formatting.
     * @param buf
     * @param d
     */
    private void appendDouble(StringBuffer buf, double d) {
        if (df == null) {
            buf.append(d);
        } else {
            buf.append(df.format(d));
        }
    }

    String toNewick(Node node, List<Function> metadataList, BranchRateModel.Base branchRateModel, double scaleFactor) {
        StringBuffer buf = new StringBuffer();
        if (node.getLeft() != null) {
            buf.append("(");
            buf.append(toNewick(node.getLeft(), metadataList, branchRateModel, scaleFactor));
            if (node.getRight() != null) {
                buf.append(',');
                buf.append(toNewick(node.getRight(), metadataList, branchRateModel, scaleFactor));
            }
            buf.append(")");
        } else {
            buf.append(node.labelNr + 1);
        }

        if (someMetaDataNeedsLogging) {
            buf.append("[&");
            if (metadataList.size() > 0) {
                for (Function metadata : metadataList) {
                    buf.append(((BEASTObject)metadata).getID());
                    buf.append('=');
                    if (metadata instanceof Parameter<?>) {
                        Parameter<?> p = (Parameter<?>) metadata;
                        int dim = p.getMinorDimension1();
                        if (dim > 1) {
                            buf.append('{');
                            for (int i = 0; i < dim; i++) {
                                buf.append(p.getMatrixValue(node.labelNr, i));
                                if (i < dim - 1) {
                                    buf.append(',');
                                }
                            }
                            buf.append('}');
                        } else {
                            buf.append(metadata.getArrayValue(node.labelNr));
                        }
                    } else {
                        buf.append(metadata.getArrayValue(node.labelNr));
                    }
                    if (metadataList.indexOf(metadata) < metadataList.size() - 1) {
                        buf.append(",");
                    }
                }
                if (branchRateModel != null) {
                    buf.append(",");
                }
            }
            if (branchRateModel != null) {
                buf.append("rate=");
                appendDouble(buf, branchRateModel.getRateForBranch(node));
            }
            buf.append(']');
        }
        buf.append(":");
        if (substitutions) {
            appendDouble(buf, node.getLength() * branchRateModel.getRateForBranch(node));
        } else {
            appendDouble(buf, node.getLength());
        }
        return buf.toString();
    }

}
