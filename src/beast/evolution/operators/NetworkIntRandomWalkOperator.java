package beast.evolution.operators;

import beast.core.*;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.operators.util.Vertex;
import beast.evolution.substitutionmodel.NtdBMA;
import beast.util.Randomizer;

import java.util.ArrayList;

/**
 * @author Chieh-Hsi Wu
 *
 */
@Description("Performs a random walk on a network which a single connected component.")
public class NetworkIntRandomWalkOperator extends Operator {

    public Input<Integer> offsetInput = new Input<Integer>(
            "offset",
            "Code off set to zero so that coding can start from another number.",
            0
    );


    protected int offset;
    public Input<ArrayList<Vertex>> vertices = new Input<ArrayList<Vertex>>("vertex", "Unique id number of the vertex.", new ArrayList<Vertex>());
    public Input<RealParameter> parameterInput = new Input<RealParameter>("parameter", "Vertex indicator which indicates the current vertex in the network.",Input.Validate.REQUIRED);
    double[][] logHastingsRatios;
    boolean[][] permittedRoutes;
    int[][] neighbours;
    int vertexCount = -1;


    public void initAndValidate() throws Exception {
        offset = offsetInput.get();
        ArrayList<Vertex> vertices = this.vertices.get();
        vertexCount = vertices.size();
        for(Vertex vertex:vertices){
            System.err.println(vertex);
        }
        
        //hasSingleComponent(vertices);
        neighbours = new int[vertices.size()][];

        for(int i = 0; i < neighbours.length;i++){
            neighbours[i] = vertices.get(i).getNeighbours();
        }


        permittedRoutes = new boolean[vertices.size()][vertices.size()];
        logHastingsRatios = new double[vertices.size()][vertices.size()];
        markPermittedRoutes();
        computeHastingsRatio();
    }

    private void markPermittedRoutes(){
        for(int i = 0; i < vertexCount;i++){
            for(int j = 0; j < neighbours[i].length; j++){
                permittedRoutes[i][neighbours[i][j]]= true;
            }
        }
    }

    private void computeHastingsRatio(){
        for(int i = 0; i < vertexCount; i++){
            for(int j = i+1; j < vertexCount;j++){
                if(permittedRoutes[i][j] != permittedRoutes[j][i]){
                    System.err.println("Route specified is not symmetrical");
                    System.err.println((permittedRoutes[i][j]) ? "Can ":"Cannot "+"get from "+i+" to "+ j);
                    System.err.println((permittedRoutes[j][i]) ? "Can ":"Cannot "+"get from "+j+" to "+ i);
                    throw new RuntimeException();
                }else{
                    logHastingsRatios[i][j] = Math.log(neighbours[i].length)-Math.log(neighbours[j].length);
                    logHastingsRatios[j][i] = Math.log(neighbours[j].length)-Math.log(neighbours[i].length);

                }
            }
        }
        /*for(int i =0; i < vertexCount; i++){
            for(int j = 0; j < vertexCount; j++){
                System.out.print(logHastingsRatios[i][j]+ " ");
            }
            System.out.println();
        }*/

    }

    public double proposal(){

        RealParameter parameter = parameterInput.get(this);
        int currVertex = (int)(double)parameter.getValue() - offset;
        int nextVertex = neighbours[currVertex][Randomizer.nextInt(neighbours[currVertex].length)];
        //System.out.println("currVertex: "+currVertex+ ", nextVertex: "+nextVertex);
        parameter.setValue(0,(double)nextVertex+offset);
        //System.out.println("currVal: "+this.parameter.get());
        //System.out.println(logHastingsRatios[currVertex][nextVertex]);
        return logHastingsRatios[currVertex][nextVertex];

    }

    
    public void hasSingleComponent(ArrayList<Vertex> vertices){

        ArrayList<Vertex> compQ = new ArrayList<Vertex>();
        compQ.add(vertices.get(0));
        ArrayList<Integer> markedQ = new ArrayList<Integer>();
        markedQ.add(vertices.get(0).getId());

        while(compQ.size() >0){

            Vertex v = compQ.remove(0);
            int[] neighbours = v.getNeighbours();

            for(int j = 0; j < neighbours.length; j++){
                System.err.println("j: "+j);

                if(markedQ.indexOf(neighbours[j]) == -1){

                    markedQ.add(neighbours[j]);
                    compQ.add(vertices.get(neighbours[j]));

                }
            }

        }

        if(markedQ.size() != vertices.size()){
            throw new RuntimeException("There are disjointed components in the network: "+markedQ.size()+" "+vertices.size());
        }
    }





    /** basic implementation of a SubstitutionModel bringing together relevant super class**/



    public static void main(String[] args){
        NetworkIntRandomWalkOperator op = new NetworkIntRandomWalkOperator();
        try{
            op.testHasSingleComponent();
        }catch(Exception e){
            throw new RuntimeException(e);
        }
    }
    public void testHasSingleComponent() throws Exception{
        Vertex jc = new Vertex();
        jc.initByName(
                Vertex.ID_NUM, NtdBMA.JC69,
                Vertex.NEIGHBOURS, "1 2"
        );

        Vertex k80 = new Vertex();
        k80.initByName(
                Vertex.ID_NUM, NtdBMA.K80,
                Vertex.NEIGHBOURS, "0 3"
        );

        Vertex f81 = new Vertex();
        f81.initByName(
                Vertex.ID_NUM, NtdBMA.F81,
                Vertex.NEIGHBOURS, "0 3"
        );

        Vertex hky85 = new Vertex();
        hky85.initByName(
                Vertex.ID_NUM, NtdBMA.HKY85,
                Vertex.NEIGHBOURS, "1 2 4"
        );

        Vertex tn93 = new Vertex();
        tn93.initByName(
                Vertex.ID_NUM, NtdBMA.TN93,
                Vertex.NEIGHBOURS, "4 5"
        );

        Vertex gtr = new Vertex();
        gtr.initByName(
                Vertex.ID_NUM, NtdBMA.GTR,
                Vertex.NEIGHBOURS, "4"
        );

        ArrayList<Vertex> vertices = new ArrayList<Vertex>();
        vertices.add(jc);
        vertices.add(k80);
        vertices.add(f81);
        vertices.add(hky85);
        vertices.add(tn93);
        vertices.add(gtr);

        for(int i = 0; i < vertices.size(); i++){
            System.out.println(vertices.get(i));
        }

        hasSingleComponent(vertices);


    }


}
