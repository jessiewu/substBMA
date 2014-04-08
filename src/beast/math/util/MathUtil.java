package beast.math.util;


import beast.util.Randomizer;

/**
 * Methods from MathUtils in BEAST1
 */
public class MathUtil  {
    /**
 * Returns sqrt(a^2 + b^2) without under/overflow.
 */
    public static double hypot(double a, double b) {
        double r;

	    if (Math.abs(a) > Math.abs(b)) {
		    r = b/a;
		    r = Math.abs(a)*Math.sqrt(1+r*r);
	    } else if (b != 0) {
		    r = a/b;
		    r = Math.abs(b)*Math.sqrt(1+r*r);
	    } else {
		    r = 0.0;
	    }
	    return r;
    }

    public static int[] sample(int size, int[] vec, boolean replace){
        int[] temp = new int[size];
        if(replace){

            for(int i = 0; i < size; i++){

                temp[i] = vec[Randomizer.nextInt(vec.length)];


            }

        }else{
            //System.out.println("Without replace");
            int sampledSiteIndex;
            int[] temp2 = new int[vec.length];
            System.arraycopy(vec, 0, temp2, 0, vec.length);
            /*for(int tempSite:temp2){
                System.out.print(tempSite+" ");
            }
        System.out.println();*/
            for(int i = 0; i < size; i++){
                sampledSiteIndex  = Randomizer.nextInt(vec.length - i);
                //System.out.println("sampledSiteIndex: "+sampledSiteIndex+" replaced index: "+(vec.length - i - 1));
                int siteTemp = temp2[vec.length - i - 1];
                temp2[vec.length - i - 1] = temp2[sampledSiteIndex];
                temp2[sampledSiteIndex] = siteTemp;
                /*System.out.print("Round "+(i+1)+": ");
                for(int tempSite:temp2){
                    System.out.print(tempSite+" ");
                }
                System.out.println();   */
            }

            System.arraycopy(temp2,vec.length-size,temp,0,size);

        }
         /*for(int tempSite:temp){
             System.out.print(tempSite+" ");
         }
        System.out.println();*/
         return temp;

    }


    public static int[] sample(int size, int[] vec, boolean replace, int[] indicesSampled){
        int[] temp = new int[size];
        if(replace){

            for(int i = 0; i < size; i++){

                temp[i] = vec[Randomizer.nextInt(vec.length)];


            }

        }else{
            int[] allIndices = new int[vec.length];
            for(int i = 0;i < allIndices.length;i++){
                allIndices[i] = i;
            }
            //System.out.println("Without replace");
            int sampledSiteIndex;
            int[] temp2 = new int[vec.length];
            System.arraycopy(vec, 0, temp2, 0, vec.length);
            /*for(int tempSite:temp2){
                System.out.print(tempSite+" ");
            }
        System.out.println();*/
            for(int i = 0; i < size; i++){
                sampledSiteIndex  = Randomizer.nextInt(vec.length - i);

                //System.out.println("sampledSiteIndex: "+sampledSiteIndex+" replaced index: "+(vec.length - i - 1));
                int siteTemp = temp2[vec.length - i - 1];
                temp2[vec.length - i - 1] = temp2[sampledSiteIndex];
                temp2[sampledSiteIndex] = siteTemp;

                int tempSiteIndex = allIndices[vec.length - i - 1];
                allIndices[vec.length - i - 1] = allIndices[sampledSiteIndex];
                allIndices[sampledSiteIndex] = tempSiteIndex;
                /*System.out.print("Round "+(i+1)+": ");
                for(int tempSite:temp2){
                    System.out.print(tempSite+" ");
                }
                System.out.println();   */
                /*System.out.println(allIndices.length);
                for(int j = 0; j < allIndices.length;j++){
                    System.out.print(allIndices[j]+" ");

                }
                System.out.println(); */

            }

            System.arraycopy(temp2,vec.length-size,temp,0,size);
            System.arraycopy(allIndices,allIndices.length-size,indicesSampled,0,size);

            /*System.out.println();
            for(int i = 0; i < size;i++){
                System.out.print(indicesSampled[i]+" ");

            }
            System.out.println();   */

        }
         /*for(int tempSite:temp){
             System.out.print(tempSite+" ");
         }
        System.out.println();*/
         return temp;

    }


   
}
