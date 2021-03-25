package beast.evolution.substitutionmodel;
import beast.evolution.datatype.DataType;
import beast.evolution.datatype.IntegerData;
import beast.evolution.tree.Node;
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
    
    public static long binomi(int n, int k) {
        if ((n == k) || (k == 0))
            return 1;
        else
            return binomi(n - 1, k) + binomi(n - 1, k - 1);
    }
	public static double bd_prob(int child, int ancestor, double bd_rate, double distance) {
		double p = 0;
		int j;
		double p_j;
		int range = Math.min(child, ancestor) + 1;
		for (j = 1; j < range; ++j ){
			p_j = binomi(ancestor, j) * binomi(child - 1, j - 1) * Math.pow(bd_rate * distance, -2 * j);
			p += p_j;
		}
		p *= Math.pow(distance * bd_rate / (1 + distance * bd_rate), child + ancestor);
		return p; 
	}

	 
	@Override
	public void getTransitionProbabilities(Node node, double startTime, double endTime, double rate, double[] matrix) {
		// TODO Auto-generated method stub
		//assume birth rate = death rate = 1
		int bd_rate = 1;
		int index;
		int i, j;
		double distance = (startTime - endTime) * rate;
		for (i = 0; i < nrOfStates -1 ; ++i) {
			for (j = 0; j < nrOfStates -1 ; ++j) {
				index = i * nrOfStates + j;
				if (i == 0) {
					if (j == 0) {
						matrix[index] = 1;
					}
					else {
						matrix[index] = 0;
					}
				} else if(j == 0){
					matrix[index] = Math.pow((bd_rate * distance) / (1 + bd_rate * distance), i);
					
				} else if (i == 1) {
					matrix[index] = Math.pow(distance, j - 1) / Math.pow((1 + distance), j + 1);
				} else {
					matrix[index] = bd_prob(j, i, bd_rate, distance);
				}
			}
		}
		
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
