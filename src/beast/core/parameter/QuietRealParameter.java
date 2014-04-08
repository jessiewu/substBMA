package beast.core.parameter;

import beast.core.Description;

/**
 * @author Chieh-Hsi Wu
 */
@Description("A real parameter class that allows values to be set quietly. It's convenience is beyond imagination.")
public class QuietRealParameter extends RealParameter{
    private int idNumber = -1;
    public QuietRealParameter(){

    }

    public void initAndValidate() throws Exception{
        super.initAndValidate();
    }

    public QuietRealParameter(Double value) throws Exception{
        super(new Double[]{value});

    }



    public QuietRealParameter(Double[] values) throws Exception{
        super(values);

    }


    public void setValueQuietly(int dim, Double value){
        values[dim] = value;
        m_bIsDirty[dim] = true;
        m_nLastDirty = dim;
    }

    public void setIDNumber(int idNumber){
        this.idNumber = idNumber;
    }

    public int getIDNumber(){
        return idNumber;
    }

    @Override
    public QuietRealParameter copy() {
        try {
            @SuppressWarnings("unchecked") final
            QuietRealParameter copy = (QuietRealParameter) this.clone();
            copy.values = values.clone();//new Boolean[values.length];
            copy.m_bIsDirty = new boolean[values.length];
            return copy;
        } catch (Exception e) {
            e.printStackTrace();
        }
        return null;
    }



}