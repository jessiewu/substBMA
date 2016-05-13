package beast.evolution.datatype;

import beast.core.Description;
import beast.core.Input;
import beast.util.HeapSort;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * @author Chieh-Hsi Wu
 */
@Description("Data type for languages coded as numbers.")
public class LanguageIntegerData extends DataType.Base{
    public Input<String> symbolsInput = new Input<String>(
            "symbols",
            "A series of comma delimited symbols used to code the data.",
            Input.Validate.REQUIRED
    );

    public Input<String> ambiguitiesInput = new Input<String>(
            "ambiguity",
            "Ambiguous state defined by the user",
            Input.Validate.REQUIRED
    );

    private HashMap<String,Integer> codeMap;




    public void initAndValidate(){
        String[] symbols = symbolsInput.get().split("\\s+");
        String[] ambiguities = ambiguitiesInput.get().split("\\s+");
        mapCodeToStateSet = new int[symbols.length+ambiguities.length+1][];


        codeMap = new HashMap<String, Integer>();
        int k = 0;
        for(int i = 0;i < symbols.length;i++){
            if(!codeMap.containsKey(symbols[i])){
                mapCodeToStateSet[k] = new int[]{k};
                codeMap.put(symbols[i],k++);

            }else{
                throw new RuntimeException("Error: Repeated symbols in code map!");
            }
        }
        stateCount = codeMap.size();


        for(int i = 0; i < ambiguities.length; i++){
            String[] ambiguityStr = ambiguities[i].trim().split("_");
            int[] ambiguity = new int[ambiguityStr.length];
            for(int j = 0; j < ambiguity.length; j++){
                ambiguity[j] = Integer.parseInt(ambiguityStr[j]);
            }
            heapSort(ambiguity);
            String ambiguityKey = ""+ambiguity[0];
            for(int j = 1; j < ambiguity.length;j++){
                ambiguityKey+="_"+ambiguity[j];

            }
            if(!codeMap.containsKey(ambiguityKey)){
                int[] stateSet = new int[ambiguity.length];
                for(int j = 0; j < stateSet.length; j++){
                    System.out.println(codeMap.get(""+ambiguity[j])+" "+ambiguity[j]);
                    stateSet[j] = codeMap.get(""+ambiguity[j]);
                }
                mapCodeToStateSet[k] = stateSet;
                codeMap.put(ambiguityKey,k++);
            }else{
                throw new RuntimeException("Error: Repeated ambiguity specified!");
            }

        }
        codeMap.put("?",k);
        codeMap.put("-",k);
        mapCodeToStateSet[k] = new int[symbols.length];
        for(int i = 0;i < symbols.length; i++){
            mapCodeToStateSet[k][i] = i;
        }

        for(int i = 0; i < mapCodeToStateSet.length; i++){
            System.out.print(i+": ");
            for(int j = 0; j < mapCodeToStateSet[i].length; j++){
                System.out.print(mapCodeToStateSet[i][j]+", ");
            }
            System.out.println();
        }

        //m_sCodeMap = symbolsInput.get();

    }

    public static void heapSort(int[] array) {

        int temp;
        int j, n = array.length;

        // turn input array into a heap
        for (j = n / 2; j > 0; j--) {
            adjust(array, j, n);
        }

        // remove largest elements and put them at the end
        // of the unsorted region until you are finished
        for (j = n - 1; j > 0; j--) {
            temp = array[0];
            array[0] = array[j];
            array[j] = temp;
            adjust(array, 1, j);
        }
    }

    private static void adjust(int[] array, int lower, int upper) {

        int j, k;
        int temp;

        j = lower;
        k = lower * 2;

        while (k <= upper) {
            if ((k < upper) && (array[k - 1] < array[k])) {
                k += 1;
            }
            if (array[j - 1] < array[k - 1]) {
                temp = array[j - 1];
                array[j - 1] = array[k - 1];
                array[k - 1] = temp;
            }
            j = k;
            k *= 2;
        }
    }



    public int string2code(String s){
        String key;
        if(s.contains("_")){
            String[] ambiguityStr = s.trim().split("_");
            int[] ambiguity = new int[ambiguityStr.length];
            for(int j = 0; j < ambiguity.length; j++){
                ambiguity[j] = Integer.parseInt(ambiguityStr[j]);
            }
            heapSort(ambiguity);
            String ambiguityKey = ""+ambiguity[0];
            for(int j = 1; j < ambiguity.length;j++){
                ambiguityKey+="_"+ambiguity[j];
            }
            key = ambiguityKey;
        }else{
            key = s;

        }
        if(codeMap.containsKey(key)){
            return codeMap.get(key);

        }else{
            throw new RuntimeException("The symbol, "+ key+",is not specified.");
        }

    }


    @Override
    public String getTypeDescription() {
        return "Language integer data";
    }
}
