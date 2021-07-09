  
package beast.evolution.substitutionmodel;
import beast.evolution.datatype.DataType;
import beast.evolution.datatype.IntegerData;
import beast.evolution.tree.Node;

import java.util.Arrays;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;


public class BD extends SubstitutionModel.Base {
	public BD() {
        // this is added to avoid a parsing error inherited from superclass because frequencies are not provided.
        frequenciesInput.setRule(Validate.OPTIONAL);
        try {
            // this call will be made twice when constructed from XML
            // but this ensures that the object is validly constructed for testing purposes.
            initAndValidate();
        } catch (Exception e) {
            e.printStackTrace();
            throw new RuntimeException("initAndValidate() call failed when constructing BD()");
        }
    }
	public Input<RealParameter> nstates = new Input<RealParameter>("nstates", "number of states to consider", Validate.REQUIRED);
	//public static final int nstates = 30;
    @Override
    public double[] getFrequencies() {
        return null;
    }
    
    public int getStateCount() {
        nrOfStates = 30;
        return nrOfStates;
    }
    protected int nrOfStates = 30;
    
    
    protected static double[] binom_array = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.0, 0.0, 0.0, -0.0, 0.0, -0.0, 0.0, -0.0, 0.0, -0.0, 0.0, -0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 3.0, 3.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 4.0, 6.0, 4.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 5.0, 10.0, 10.0, 5.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 6.0, 15.0, 20.0, 15.0, 6.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 7.0, 21.0, 35.0, 35.0, 21.0, 7.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 8.0, 28.0, 56.0, 70.0, 56.0, 28.0, 8.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 9.0, 36.0, 84.0, 126.0, 126.0, 84.0, 36.0, 9.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 10.0, 45.0, 120.0, 210.0, 252.0, 210.0, 120.0, 45.0, 10.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 11.0, 55.0, 165.0, 330.0, 462.0, 462.0, 330.0, 165.0, 55.0, 11.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 12.0, 66.0, 220.0, 495.0, 792.0, 924.0, 792.0, 495.0, 220.0, 66.0, 12.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 13.0, 78.0, 286.0, 715.0, 1287.0, 1716.0, 1716.0, 1287.0, 715.0, 286.0, 78.0, 13.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 14.0, 91.0, 364.0, 1001.0, 2002.0, 3003.0, 3432.0, 3003.0, 2002.0, 1001.0, 364.0, 91.0, 14.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 15.0, 105.0, 455.0, 1365.0, 3003.0, 5005.0, 6435.0, 6435.0, 5005.0, 3003.0, 1365.0, 455.0, 105.0, 15.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 16.0, 120.0, 560.0, 1820.0, 4368.0, 8008.0, 11440.0, 12870.0, 11440.0, 8008.0, 4368.0, 1820.0, 560.0, 120.0, 16.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 17.0, 136.0, 680.0, 2380.0, 6188.0, 12376.0, 19448.0, 24310.0, 24310.0, 19448.0, 12376.0, 6188.0, 2380.0, 680.0, 136.0, 17.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 18.0, 153.0, 816.0, 3060.0, 8568.0, 18564.0, 31824.0, 43758.0, 48620.0, 43758.0, 31824.0, 18564.0, 8568.0, 3060.0, 816.0, 153.0, 18.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 19.0, 171.0, 969.0, 3876.0, 11628.0, 27132.0, 50388.0, 75582.0, 92378.0, 92378.0, 75582.0, 50388.0, 27132.0, 11628.0, 3876.0, 969.0, 171.0, 19.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 20.0, 190.0, 1140.0, 4845.0, 15504.0, 38760.0, 77520.0, 125970.0, 167960.0, 184756.0, 167960.0, 125970.0, 77520.0, 38760.0, 15504.0, 4845.0, 1140.0, 190.0, 20.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 21.0, 210.0, 1330.0, 5985.0, 20349.0, 54264.0, 116280.0, 203490.0, 293930.0, 352716.0, 352716.0, 293930.0, 203490.0, 116280.0, 54264.0, 20349.0, 5985.0, 1330.0, 210.0, 21.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 22.0, 231.0, 1540.0, 7315.0, 26334.0, 74613.0, 170544.0, 319770.0, 497420.0, 646646.0, 705432.0, 646646.0, 497420.0, 319770.0, 170544.0, 74613.0, 26334.0, 7315.0, 1540.0, 231.0, 22.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 23.0, 253.0, 1771.0, 8855.0, 33649.0, 100947.0, 245157.0, 490314.0, 817190.0, 1144066.0, 1352078.0, 1352078.0, 1144066.0, 817190.0, 490314.0, 245157.0, 100947.0, 33649.0, 8855.0, 1771.0, 253.0, 23.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 24.0, 276.0, 2024.0, 10626.0, 42504.0, 134596.0, 346104.0, 735471.0, 1307504.0, 1961256.0, 2496144.0, 2704156.0, 2496144.0, 1961256.0, 1307504.0, 735471.0, 346104.0, 134596.0, 42504.0, 10626.0, 2024.0, 276.0, 24.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 25.0, 300.0, 2300.0, 12650.0, 53130.0, 177100.0, 480700.0, 1081575.0, 2042975.0, 3268760.0, 4457400.0, 5200300.0, 5200300.0, 4457400.0, 3268760.0, 2042975.0, 1081575.0, 480700.0, 177100.0, 53130.0, 12650.0, 2300.0, 300.0, 25.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 26.0, 325.0, 2600.0, 14950.0, 65780.0, 230230.0, 657800.0, 1562275.0, 3124550.0, 5311735.0, 7726160.0, 9657700.0, 10400600.0, 9657700.0, 7726160.0, 5311735.0, 3124550.0, 1562275.0, 657800.0, 230230.0, 65780.0, 14950.0, 2600.0, 325.0, 26.0, 1.0, 0.0, 0.0, 0.0, 1.0, 27.0, 351.0, 2925.0, 17550.0, 80730.0, 296010.0, 888030.0, 2220075.0, 4686825.0, 8436285.0, 13037895.0, 17383860.0, 20058300.0, 20058300.0, 17383860.0, 13037895.0, 8436285.0, 4686825.0, 2220075.0, 888030.0, 296010.0, 80730.0, 17550.0, 2925.0, 351.0, 27.0, 1.0, 0.0, 0.0, 1.0, 28.0, 378.0, 3276.0, 20475.0, 98280.0, 376740.0, 1184040.0, 3108105.0, 6906900.0, 13123110.0, 21474180.0, 30421755.0, 37442160.0, 40116600.0, 37442160.0, 30421755.0, 21474180.0, 13123110.0, 6906900.0, 3108105.0, 1184040.0, 376740.0, 98280.0, 20475.0, 3276.0, 378.0, 28.0, 1.0, 0.0, 1.0, 29.0, 406.0, 3654.0, 23751.0, 118755.0, 475020.0, 1560780.0, 4292145.0, 10015005.0, 20030010.0, 34597290.0, 51895935.0, 67863915.0, 77558760.0, 77558760.0, 67863915.0, 51895935.0, 34597290.0, 20030010.0, 10015005.0, 4292145.0, 1560780.0, 475020.0, 118755.0, 23751.0, 3654.0, 406.0, 29.0, 1.0}; 
    
    public static double binomi(int n, int k) {

    	return binom_array[n*30+k];
    }
	public static double bd_prob(int child, int ancestor, double bd_rate, double distance) {
		double p = 0;
		int j;
		double p_j;
		int range = Math.min(child, ancestor) + 1;
		if (distance <=1) {
			for (j = 1; j < range; ++j ){
				p_j = binomi(ancestor, j) * binomi(child - 1, j - 1) * Math.pow(bd_rate * distance, child + ancestor -2*j);
				p += p_j;
			}
			p = p * Math.pow(bd_rate / (1 + distance * bd_rate), child + ancestor);
			//if (Math.pow(bd_rate / (1 + distance * bd_rate), child + ancestor) == 0) {
				//p = 0;
			//}	
		}
		else {
			for (j = 1; j < range; ++j ){
				p_j = binomi(ancestor, j) * binomi(child - 1, j - 1) * Math.pow(bd_rate * distance, -2 * j);
				p += p_j;
			}
			p = p * Math.pow(bd_rate*distance / (1 + distance * bd_rate), child + ancestor);
			if (Math.pow(bd_rate*distance / (1 + distance * bd_rate), child + ancestor) == 0) {
				//p = 0;
			}
		}

		return p; 
	}
	
	protected boolean checkTransitionMatrix(double[] matrix) {
		double sum = 0;
		int i, j;
		int index;
		for (i = 0; i < nrOfStates ; ++i) {
			for (j = 0; j < nrOfStates ; ++j) {
				index = i * nrOfStates + j;
				sum = sum + matrix[index]; 
			}
			if (sum > 1.01 |sum < 0.95) {
				//System.out.println("current index:" + i);
				//System.out.println(sum);
				return true;
			}
		sum = 0;
		}		
		return true;
		
	}

	 
	@Override
	public void getTransitionProbabilities(Node node, double startTime, double endTime, double rate, double[] matrix) {
		// TODO Auto-generated method stub
		//assume birth rate = death rate = 1
		//System.out.println("TRANSITIONPROB");
		double bd_rate = 1;
		int index;
		int i, j;
		double prob;
		double distance = (startTime - endTime) * rate;
		//System.out.println("D" + distance);
		for (i = 0; i < nrOfStates ; ++i) {
			for (j = 0; j < nrOfStates ; ++j) {
				index = i * nrOfStates + j;
				if (i == 0) {
					if (j == 0) {
						prob = 1;
						//matrix[index] = 1;
					}
					else {
						prob = 0;
						//matrix[index] = 0;
					}
				} else if(j == 0){
					prob =  Math.pow((bd_rate * distance) / (1 + bd_rate * distance), i);
					//System.out.println("j == 0, " + prob);
					//matrix[index] = Math.pow((bd_rate * distance) / (1 + bd_rate * distance), i);
					
				} else if (i == 1) {
					prob = Math.pow(distance, j - 1) / Math.pow((1 + distance), j + 1);
					//System.out.println("i == 1, " + prob);
					//matrix[index] = Math.pow(distance, j - 1) / Math.pow((1 + distance), j + 1);
				} else {
					prob = bd_prob(j, i, bd_rate, distance);

					//matrix[index] = bd_prob(j, i, bd_rate, distance);
					//System.out.println("else, " + prob);
				}
			matrix[index] = prob;
			}
		}
		//assert checkTransitionMatrix(matrix):"Transition Matrix does not sum up to 1.";

	}

	@Override
	public EigenDecomposition getEigenDecomposition(Node node) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public boolean canHandleDataType(DataType dataType) {
		// TODO Auto-generated method stub
		return dataType instanceof IntegerData;
	}

}