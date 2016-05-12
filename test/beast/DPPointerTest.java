package beast;

import beast.core.parameter.*;
import beast.core.State;
//import com.sun.servicetag.SystemEnvironment;
import junit.framework.TestCase;

/**
 * @author Chieh-Hsi Wu
 */
public class DPPointerTest extends TestCase {

    interface Instance {
        public void setup();
        public DPPointer getDPPointer();
        public ParameterList getParameterList();

    }

    Instance test0 = new Instance(){
        ParameterList paramList;
        DPPointer pointer;
        public void setup(){
            try{

                QuietRealParameter parameter1 = new QuietRealParameter(new Double[]{0.0});
                QuietRealParameter parameter2 = new QuietRealParameter(new Double[]{1.0});
                QuietRealParameter parameter3 = new QuietRealParameter(new Double[]{2.0});
                paramList = new ParameterList();
                paramList.initByName(
                        "parameter", parameter1,
                        "parameter", parameter2,
                        "parameter", parameter3
                );
                paramList.setID("parameterList");
                IntegerParameter initAssign = new IntegerParameter(new Integer[]{0,2,1,1,2});
                pointer = new DPPointer();
                pointer.initByName(
                        "uniqueParameter", parameter1,
                        "uniqueParameter", parameter2,
                        "uniqueParameter", parameter3,
                        "initialAssignment", initAssign
                );
                pointer.setID("pointer");
            }catch(Exception e){
                throw new RuntimeException(e);
            }
        }
        public DPPointer getDPPointer(){
            return pointer;
        }

        public ParameterList getParameterList(){
            return paramList;
        }


        

    };



    Instance[] tests = new Instance[]{test0};
    public void testDDPointer(){
        try{
        for(Instance test: tests){
            test.setup();
            DPPointer pointer = test.getDPPointer();
            ParameterList paramList =  test.getParameterList();
            State state = new State();

        state.initByName(
                "stateNode", pointer,
                "stateNode", paramList
        );
        state.initialise();

            operationEx1(pointer,paramList);
        }
        }catch(Exception e){
            throw new RuntimeException(e);

        }

    }



    public void operationEx1(DPPointer pointer,ParameterList paramList){
        try{

            QuietRealParameter newVal = new QuietRealParameter(new Double[]{3.0});
            int index = 3;
            int[] assignment = new int[]{0,2,1,1,2};
            RealParameter[] expectedPointer= new RealParameter[assignment.length];
            RealParameter[] expectedList = new RealParameter[assignment.length];
            operation1(pointer, paramList, index, newVal,assignment,expectedPointer,expectedList);

            for(int i = 0; i < expectedPointer.length; i++)
                assertTrue(pointer.sameParameter(i,expectedPointer[i]));

            for(int i = 0; i < paramList.getDimension(); i++)
                assertTrue(paramList.getParameter(i)==expectedList[i]);

            index = 1;
            int existingValIndex = 2;
            assignment = new int[]{0,2,1,3,2};
            operation2(pointer, paramList, index, existingValIndex,assignment,expectedPointer);
            for(int i = 0; i < expectedPointer.length; i++){
                //System.err.println(expectedPointer[i]);
                assertTrue(pointer.sameParameter(i,expectedPointer[i]));
            }

            for(int i = 0; i < paramList.getDimension(); i++)
                assertTrue(paramList.getParameter(i)==expectedList[i]);

            index = 0;
            newVal = new QuietRealParameter(new Double[]{4.0});
            assignment = new int[]{0,1,1,3,2};
            operation3(pointer, paramList, index, newVal,assignment,expectedPointer,expectedList);

            for(int i = 0; i < expectedPointer.length; i++){
                assertTrue(pointer.sameParameter(i,expectedPointer[i]));
            }

            for(int i = 0; i < paramList.getDimension(); i++){
                assertTrue(paramList.getParameter(i) == expectedList[i]);
            }

            index = 3;
            existingValIndex = 4;
            assignment = new int[]{3,0,0,2,1};
            operation4(pointer, paramList, index, existingValIndex,assignment,expectedPointer,expectedList);

            for(int i = 0; i < expectedPointer.length; i++){
                assertTrue(pointer.sameParameter(i,expectedPointer[i]));
            }

            for(int i = 0; i < paramList.getDimension(); i++)
                assertTrue(paramList.getParameter(i)==expectedList[i]);

            newVal = new QuietRealParameter(new Double[]{5.0});
            index = 4;
            //assignment = new int[]{3,0,0,1,1};
            operation5(pointer, paramList, index, newVal);

            for(int i = 0; i < expectedPointer.length; i++){
                assertTrue(pointer.sameParameter(i,expectedPointer[i]));
            }

            for(int i = 0; i < paramList.getDimension(); i++)
                assertTrue(paramList.getParameter(i)==expectedList[i]);


            index = 1;
            existingValIndex = 4;
            //assignment = new int[]{3,0,0,1,1};
            operation6(pointer, paramList, index, existingValIndex, expectedPointer);

            for(int i = 0; i < expectedPointer.length; i++)
                assertTrue(pointer.sameParameter(i,expectedPointer[i]));


            for(int i = 0; i < paramList.getDimension(); i++)
                assertTrue(paramList.getParameter(i)==expectedList[i]);

            newVal = new QuietRealParameter(new Double[]{6.0});
            index = 0;
            //assignment = new int[]{3,0,0,1,1};
            operation7(pointer, paramList, index, newVal);

            for(int i = 0; i < expectedPointer.length; i++)
                assertTrue(pointer.sameParameter(i,expectedPointer[i]));


            for(int i = 0; i < paramList.getDimension(); i++)
                assertTrue(paramList.getParameter(i)==expectedList[i]);


            index = 0;
            existingValIndex = 3;
            //assignment = new int[]{3,0,0,1,1};
            operation8(pointer, paramList, index, existingValIndex, expectedPointer);

            for(int i = 0; i < expectedPointer.length; i++)
                assertTrue(pointer.sameParameter(i,expectedPointer[i]));


            for(int i = 0; i < paramList.getDimension(); i++)
                assertTrue(paramList.getParameter(i)==expectedList[i]);


        }catch(Exception e){
            throw new RuntimeException(e);

        }

    }

    private void operation1(
            DPPointer pointer,
            ParameterList paramList,
            int index,
            QuietRealParameter newVal,
            int[] assignment,
            RealParameter[] expectedPointer,
            RealParameter[] expectedList) throws Exception{


        for(int i = 0; i < expectedPointer.length; i++)
            expectedPointer[i] = paramList.getParameter(assignment[i]);
        expectedPointer[index] = newVal;

        for(int i = 0; i < paramList.getDimension(); i++)
            expectedList[i] = paramList.getParameter(i);
        expectedList[paramList.getDimension()] = newVal;

        pointer.point(index,newVal);
        paramList.addParameter(newVal);

        pointer.setEverythingDirty(false);
        paramList.setEverythingDirty(false);

    }
    
    private void operation2(
            DPPointer pointer,
            ParameterList paramList,
            int index,
            int existingValIndex,
            int[] assignment,
            RealParameter[] expectedPointer) throws Exception{

        for(int i = 0; i < pointer.getDimension();i++){
            expectedPointer[i] = paramList.getParameter(assignment[i]);
        }        

        int listIndex = paramList.indexOf(expectedPointer[existingValIndex]);
        expectedPointer[index] = expectedPointer[existingValIndex];

        pointer.point(index,paramList.getParameter(listIndex));

        pointer.setEverythingDirty(false);
        paramList.setEverythingDirty(false);


    }

    private void operation3(
            DPPointer pointers,
            ParameterList paramList,
            int index,
            QuietRealParameter newVal,
            int[] assignment,
            RealParameter[] expectedPointer,
            RealParameter[] expectedList) throws Exception{


        for(int i = 0; i < expectedPointer.length; i++)
            expectedPointer[i] = paramList.getParameter(assignment[i]);
        expectedPointer[index] = newVal;


        int counter = 0;
        for(int i = 0; i < paramList.getDimension(); i++){
            if(assignment[index] != i){
                expectedList[counter++] = paramList.getParameter(i);
            }
        }
        expectedList[counter] = newVal;

        DPValuable dpValuable = new DPValuable();
        dpValuable.initByName(
                "paramList", paramList,
                "pointers", pointers
        );

        pointers.point(index,newVal);
        paramList.addParameter(newVal);
        dpValuable.requiresRecalculation();
        int[] counts = dpValuable.getClusterCounts();

        for(int i = 0; i < counts.length;i++){
            if(counts[i] == 0)
                paramList.removeParameter(i);
        }

        pointers.setEverythingDirty(false);
        paramList.setEverythingDirty(false);
    }

    private void operation4(
            DPPointer pointers,
            ParameterList paramList,
            int index,
            int existingValIndex,
            int[] assignment,
            RealParameter[] expectedPointer,
            RealParameter[] expectedList) throws Exception{

        for(int i = 0; i < pointers.getDimension();i++){
            expectedPointer[i] = paramList.getParameter(assignment[i]);
        }

        int listIndex = paramList.indexOf(expectedPointer[existingValIndex]);
        expectedPointer[index] = expectedPointer[existingValIndex];
        int counter = 0;
        for(int i = 0; i < paramList.getDimension(); i++){
            if(assignment[index] != i){
                
                expectedList[counter++] = paramList.getParameter(i);
            }
        }

        DPValuable dpValuable = new DPValuable();
        dpValuable.initByName(
                "paramList", paramList,
                "pointers", pointers
        );

        pointers.point(index,paramList.getParameter(listIndex));
        dpValuable.requiresRecalculation();
        int[] counts = dpValuable.getClusterCounts();
        for(int i = 0; i < counts.length;i++){
            if(counts[i] == 0)
                paramList.removeParameter(i);
        }
        pointers.setEverythingDirty(false);
        paramList.setEverythingDirty(false);


    }

    private void operation5(
            DPPointer pointer,
            ParameterList paramList,
            int index,
            QuietRealParameter newVal) throws Exception{
      
        pointer.point(index,newVal);
        paramList.addParameter(newVal);

        pointer.restore();
        paramList.restore();
        
        pointer.setEverythingDirty(false);
        paramList.setEverythingDirty(false);

    }


    private void operation6(
            DPPointer pointer,
            ParameterList paramList,
            int index,
            int existingValIndex,
            RealParameter[] expectedPointer) throws Exception{

        int listIndex = paramList.indexOf(expectedPointer[existingValIndex]);

        pointer.point(index,paramList.getParameter(listIndex));

        pointer.restore();

        pointer.setEverythingDirty(false);
        paramList.setEverythingDirty(false);


    }

    private void operation7(
            DPPointer pointers,
            ParameterList paramList,
            int index,
            QuietRealParameter newVal) throws Exception{



        DPValuable dpValuable = new DPValuable();
        dpValuable.initByName(
                "paramList", paramList,
                "pointers", pointers
        );

        pointers.point(index,newVal);
        paramList.addParameter(newVal);
        dpValuable.requiresRecalculation();
        int[] counts = dpValuable.getClusterCounts();
        for(int i = 0; i < counts.length;i++){
            if(counts[i] == 0)
                paramList.removeParameter(i);
        }

        pointers.restore();
        paramList.restore();
        dpValuable.restore();

        pointers.setEverythingDirty(false);
        paramList.setEverythingDirty(false);
    }

    private void operation8(
            DPPointer pointers,
            ParameterList paramList,
            int index,
            int existingValIndex,
            RealParameter[] expectedPointer) throws Exception{


        int listIndex = paramList.indexOf(expectedPointer[existingValIndex]);

        DPValuable dpValuable = new DPValuable();
        dpValuable.initByName(
                "paramList", paramList,
                "pointers", pointers
        );

        pointers.point(index,paramList.getParameter(listIndex));
        dpValuable.requiresRecalculation();
        int[] counts = dpValuable.getClusterCounts();
        for(int i = 0; i < counts.length;i++){
            if(counts[i] == 0)
                paramList.removeParameter(i);
        }

        pointers.restore();
        paramList.restore();
        dpValuable.restore();

        pointers.setEverythingDirty(false);
        paramList.setEverythingDirty(false);
    }
}
