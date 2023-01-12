package beast.evolution.errormodel;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.datatype.DataType;
import beast.evolution.datatype.IntegerData;

@Description("Negative Binomial error model with parameters as success probabilities")
public class NormalErrorModel extends ErrorModel {
	final public Input<RealParameter> vrInput = new Input<>("vr","variance residue");
	final public Input<RealParameter> vfInput = new Input<>("vf","varaince_factor");
	final public Input<RealParameter> errorprob = new Input<>("p","probability of the inferred copy number is not correct");
	public Input<RealParameter> nstate = new Input<RealParameter>("nstate", "same as what in BD model", Validate.REQUIRED);
	protected static int nrOfStates;
	protected static double[] binom_array;
	private RealParameter v_r;

	@Override
	public void initAndValidate() {
		// init base
		super.initAndValidate();
        v_r = vrInput.get();
		if (updateMatrix) {
			setupErrorMatrix();
			updateMatrix = false;
		}
		if (nstate.get() == null) {
            throw new IllegalArgumentException("number of states to consider is required");
        }
	}

	@Override
	public void setupErrorMatrix() {
		nrOfStates = (int)Math.round(nstate.get().getValue());
		if (errorMatrix == null) {
			errorMatrix = new double[nrOfStates][nrOfStates];
		}
		for (int observedState = 0; observedState < nrOfStates; observedState++){
			double sum_ = 0;
			
			for (int trueState = 0; trueState < nrOfStates; trueState++){
				// rows are observed states X, columns are true states Y
				double temp = getProbability(observedState, trueState);
            	//System.out.println(observedState);
            	//System.out.println(trueState);
            	//System.out.println(getProbability(observedState, trueState));
				double p = 1.0;
				if (errorprob.get() != null){
					 p = errorprob.get().getValue();}
				if (observedState == trueState){
					temp = temp * p + 1.0 - p;
				}
				else {
					temp = temp * (1.0-p);
				}
				errorMatrix[observedState][trueState] = temp;
				sum_ = sum_ + temp;
			}
			//normalize to make sure probability summarize to 1
			//System.out.println(observedState);
			for (int trueState = 0; trueState < nrOfStates; trueState++){
				// rows are observed states X, columns are true states Y
				//System.out.println(observedState);
				errorMatrix[observedState][trueState] = errorMatrix[observedState][trueState]/sum_;
			}
		}
	}
	
	 public static double pdf(double x) {
	        return Math.exp(-x*x / 2) / Math.sqrt(2 * Math.PI);
	    }

	 
    // return pdf(x, mu, signma) = Gaussian pdf with mean mu and stddev sigma
    public static double pdf(double x, double mu, double sigma) {
    	// System.out.println(pdf((x - mu) / sigma) / sigma);
        return pdf((x - mu) / sigma) / sigma;
    }
	
    
	@Override
	public double getProbability(int observedState, int trueState) {
		double vf_val = 0;
		// TODO Auto-generated method stub
        if (vfInput.get() != null) {
        	vf_val = vfInput.get().getValue();
        }
		double variance = v_r.getValue() + observedState * vf_val;
		//System.out.println(variance);
		double prob;
		if (variance == 0) {
			if (observedState == trueState) {
				return 1;
			}
			else {
				return 0;
			}
		}
		prob = pdf(trueState, observedState, variance);
		//System.out.println(prob);
		return prob;
	}

	@Override
	public double[] getProbabilities(int observedState) {
		if (updateMatrix) {
			setupErrorMatrix();
			updateMatrix = false;
		}
		return errorMatrix[observedState];
	}

	@Override
	public boolean canHandleDataType(DataType dataType) {
		return dataType instanceof IntegerData;
	}

	@Override
	public double[] getProbabilities(int observedState, double w, int totalread) {
		// TODO Auto-generated method stub
		return null;
	}


}
