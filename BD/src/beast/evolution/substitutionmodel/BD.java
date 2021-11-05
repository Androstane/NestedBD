  
package beast.evolution.substitutionmodel;
import beast.evolution.datatype.DataType;
import beast.evolution.datatype.IntegerData;
import beast.evolution.tree.Node;

import java.util.Arrays;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;


public class BD extends SubstitutionModel.Base {
	public Input<RealParameter> nstate = new Input<RealParameter>("nstate", "same as what in BD model", Validate.REQUIRED);
	protected static int nrOfStates;
	protected static double[] binom_array;

	public BD() {
        // this is added to avoid a parsing error inherited from superclass because frequencies are not provided.
		frequenciesInput.setRule(Validate.OPTIONAL);
        try {
            // this call will be made twice when constructed from XML
            // but this ensures that the object is validly constructed for testing purposes.
            //System.out.println("Class Constructor");
        	//System.out.println(nstate.get().getValue());
            //initAndValidate();
        } catch (Exception e) {
            e.printStackTrace();
            throw new RuntimeException("initAndValidate() call failed when constructing BD()");
        }
    }
	
	@Override
    public void initAndValidate() {
		int i, j, index;
        super.initAndValidate();
        //System.out.println("BD");
        //System.out.println(nstate.get().getValue());
        if (nstate.get() == null) {
            throw new IllegalArgumentException("number of states to consider is required");
        }
        nrOfStates = (int)Math.round(nstate.get().getValue());
        //nrOfStates = 30;
        binom_array = new double[nrOfStates * nrOfStates];
        for (i = 0; i < nrOfStates ; ++i) {
			for (j = 0; j < nrOfStates ; ++j) {
				index = i * nrOfStates + j;
				binom_array[index] = binomialCoeff(i,j);
			}
        }
        //trMatrix = new double[(nrOfStates - 1) * (nrOfStates - 1)];
    }
	

	//public static final int nstates = 30;
    @Override
    public double[] getFrequencies() {
        return null;
    }
    
    public int getStateCount() {
        return nrOfStates;
    }
    //protected int nrOfStates = 30;
    
    
   
    public static double binomi(int n, int k) {
    	return binom_array[n*nrOfStates+k];
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
	public int binomialCoeff(int n, int k) {
	    int res = 1;
	 
	    // Since C(n, k) = C(n, n-k)
	    if (k > n - k)
	        k = n - k;
	 
	    // Calculate value of
	    // [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]
	    for (int i = 0; i < k; ++i) {
	        res *= (n - i);
	        res /= (i + 1);
	    }
	 
	    return res;
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