package beast.evolution.operators.util;

import beast.core.BEASTObject;
import beast.core.Input;
import beast.core.Description;

/**
 * @author Chieh-Hsi Wu
 */
@Description("This class represents a vertex in a network.")
public class Vertex extends BEASTObject {
    public static final String ID_NUM = "idNum";
    public static final String NEIGHBOURS = "neighbours";
    public Input<Integer> idNum = new Input<Integer>(ID_NUM, "Unique id number of the vertex.");
    public Input<String> neighbours = new Input<String>(NEIGHBOURS, "A string of the id number of adjacent nodes separated by a while space");
         private int id;
    private int[] neighbourIds;
         public void initAndValidate() {
        id = idNum.get();
        String[] neighboursStr = neighbours.get().split("\\s+");
        neighbourIds = new int[neighboursStr.length];
        for(int i = 0; i < neighboursStr.length;i++){
            neighbourIds[i] = Integer.parseInt(neighboursStr[i]);
        }
         }

         public int getId(){
        return id;
    }
         public int[] getNeighbours(){
        return neighbourIds;
    }
         public String toString(){
        return "Id: "+id+ ", neighbours: "+neighbours.get() ;
    }
}
