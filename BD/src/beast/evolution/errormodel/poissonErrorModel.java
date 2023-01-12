package beast.evolution.errormodel;
import org.apache.commons.math.distribution.PoissonDistributionImpl;
import java.lang.Math;
import java.util.Arrays;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.datatype.DataType;
import beast.evolution.datatype.IntegerData;
import beast.math.distributions.Poisson;

public class poissonErrorModel extends ErrorModel{
	public Input<RealParameter> nstate = new Input<RealParameter>("nstate", "same as what in BD model", Validate.REQUIRED);
	final public Input<RealParameter> smoothing = new Input<>("smoothing","smoothing parameter");
	final public Input<RealParameter> ploidy = new Input<>("ploidy","ploidy level of genome");
	protected static int nrOfStates;

	@Override
	public void initAndValidate(){
		super.initAndValidate();
		if (nstate.get() == null) {
            throw new IllegalArgumentException("number of states to consider is required");
        }
	};
	public double getProbability(int nread, int copynumber, double w, double ploidy, int totalread) {
		// TODO Auto-generated method stub
		if (nread == 0) {
			if (copynumber == 0) {
				return 1.0;
			}
			return 0.0;
		}

		//mu: expected number of reads at current bin
		//convert to integer for poisson distribution 
		//assume expected number of reads to be 1 when copy number is zero to allow sample errors 
		double mu = totalread * copynumber * w/ploidy;
		int mu_int = (int)mu;
		if (mu_int == 0) {

			mu_int = 1;
		}
		double p = new PoissonDistributionImpl(mu_int).probability(nread);
		return p;

	}
	@Override
	public double[] getProbabilities(int observedState, double w, int trueState) {
		nrOfStates = (int)Math.round(nstate.get().getValue());
        double[] prob = new double[nrOfStates];
        double max_ = 0;
        //double max_copy = 0;
        double plevel = 2.5;
        if (ploidy.get() != null){
        	plevel = ploidy.get().getValue();
        }
        for (int i = 0; i < nrOfStates; i++) {
            prob[i] = getProbability(observedState, i, w, plevel, trueState);
            if (prob[i] > max_){
            	max_ = prob[i];
            }
        }
        //most likely outlier in data, assume the readcount is not informative 
        if (max_ == 0) {
        	for (int i = 0; i < nrOfStates; i++) {
                prob[i] = 1.0/nrOfStates;
                //System.out.println(1/nrOfStates);
                //System.out.println(prob[i]);
        	}
        	//System.out.println(Arrays.toString(prob));
        }
        //smoothing/unsmoothing
        //System.out.print(observedState);
    	//System.out.println(Arrays.toString(prob));
        double s = 1.0;
        if (smoothing.get() != null) {
        	s = smoothing.get().getValue();
        }
    	double sum_ = 0.0;
    	//normalize probability 
        if (s != 1.0) {
        	for (int i = 0; i < nrOfStates; i++) {
                prob[i] = Math.pow(prob[i], s);
                sum_  += prob[i];
        	}
        }
        else {
        	for (int i = 0; i < nrOfStates; i++) {
                sum_  += prob[i];
        	}
        }
    	for (int i = 0; i < nrOfStates; i++) {
            prob[i] = prob[i]/sum_;
    	}
    	//System.out.println(Arrays.toString(prob));
        return prob;
    }
	
	@Override
	public boolean canHandleDataType(DataType dataType) {
		// TODO Auto-generated method stub
		return dataType instanceof IntegerData;
	}


	@Override
	public double getProbability(int observedState, int trueState) {
		throw new IllegalArgumentException("weight is required for read count data, another overloaded version of this method should be selected instead");
	}
	@Override
	public double[] getProbabilities(int observedState) {
		throw new IllegalArgumentException("Probabilities:weight is required for read count data, another overloaded version of this method should be selected instead");
	}




}
