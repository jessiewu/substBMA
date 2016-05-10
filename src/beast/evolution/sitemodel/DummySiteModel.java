package beast.evolution.sitemodel;

import beast.core.Description;
import beast.core.parameter.QuietRealParameter;

/**
 * @author Chieh-Hsi Wu
 */
@Description("This site model is used for computing tree likelihoods for operator move. " +
        "It is used by DP samplers so that we can use the same site model object over and over again which is required by the move.")
public class DummySiteModel extends SiteModel {
    public void initAndValidate() {

        super.initAndValidate();
        if(!(muParameter instanceof QuietRealParameter) ){
            throw new RuntimeException("Quiet mu parameter required");
        }

    }
    public QuietRealParameter getRateParameter(){
        return (QuietRealParameter)muParameter;
    }

    public void setShapeParameterValueQuietly(double shape){
        ((QuietRealParameter)shapeParameter).setValueQuietly(0,shape);
    }
}
