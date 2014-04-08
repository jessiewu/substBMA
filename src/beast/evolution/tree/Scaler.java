package beast.evolution.tree;

import beast.core.CalculationNode;
import beast.core.Description;

/**
 * @author Chieh-Hsi Wu
 */
@Description("This class is used to scale things.")
public abstract class Scaler extends CalculationNode {
    public abstract double getScaleFactor();
}
