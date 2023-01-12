package beast.evolution.errormodel;
import org.apache.commons.math.special.Gamma;
import java.lang.Math;
import java.util.Arrays;

import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.datatype.DataType;
import beast.evolution.datatype.IntegerData;

public class readcountErrorModel extends ErrorModel{
	public Input<RealParameter> nstate = new Input<RealParameter>("nstate", "same as what in BD model", Validate.REQUIRED);
	final public Input<RealParameter> theta = new Input<>("theta","inverse dispersion of negative binomail distribution", Input.Validate.REQUIRED);
	final public Input<RealParameter> smoothing = new Input<>("smoothing","inverse dispersion of negative binomail distribution");
	final public Input<RealParameter> ploidy = new Input<>("ploidy","ploidy level of genome");
	protected static int nrOfStates;

	@Override
	public void initAndValidate(){
		super.initAndValidate();
		if (nstate.get() == null) {
            throw new IllegalArgumentException("number of states to consider is required");
        }
		if (theta.get() == null) {
			throw new IllegalArgumentException("inverse dispersion of negative binomail distribution is required");
		}
		if (theta.get().getValue() == 0) {
			throw new IllegalArgumentException("inverse dispersion cannot be zero");
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
		double theta_val = theta.get().getValue();
		//expected read count for the region 
		//assume expected number of reads to be 1 when copy number is zero to allow sample errors 
		double mu = totalread * copynumber * w/ploidy;
		if (mu == 0.0) {
			mu = 1.0;
		}
		double val1 = Gamma.logGamma(nread + theta_val) - Gamma.logGamma(nread + 1) - Gamma.logGamma(theta_val); 
		double val2 = theta_val * Math.log(theta_val/(theta_val + mu));
		double val3 = nread * Math.log(mu/(theta_val + mu));
		//if (val1 > 1) {
			//System.out.println("val1:" + val1);
			//System.out.println(Gamma.logGamma(nread + theta_val) - Gamma.logGamma(nread + 1) - Gamma.logGamma(theta_val));
		//System.out.println("nreads:" + nread + " ,copy number " + copynumber + ", expected reads" + mu + ", prob:" + Math.exp(val1 + val2 + val3));
		//}
		//System.out.println(val2 + val3);
		//if (Double.isNaN(Math.exp(val1 + val2 + val3))){
			//System.out.println("nreads:" + nread + " ,copy number " + copynumber + ", expected reads" + mu + ", prob:" + Math.exp(val1 + val2 + val3));
		//}
		double p = Math.exp(val1 + val2 + val3);
		//if (p == 0 && copynumber != 0){
			//System.out.println("nreads:" + nread + " ,copy number " + copynumber + ", expected reads" + mu + ", prob:" + Math.exp(val1 + val2 + val3));
		//}
		return p;
	}
	@Override
	public double[] getProbabilities(int observedState, double w, int trueState) {
		nrOfStates = (int)Math.round(nstate.get().getValue());
        double[] prob = new double[nrOfStates];
        double max_ = 0;
        //double max_copy = 0;
        double plevel = 2.0;
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
        double s = 2.0;
        if (smoothing.get() != null) {
        	s = smoothing.get().getValue();
        }
        //if (s != 1.0) {
        	double sum_ = 0.0;
        	for (int i = 0; i < nrOfStates; i++) {
                prob[i] = Math.pow(prob[i], s);
                sum_  += prob[i];
        	}
        	for (int i = 0; i < nrOfStates; i++) {
                prob[i] = prob[i]/sum_;
        	}
        //}

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
		throw new IllegalArgumentException("weight is required for read count data, another overloaded version of this method should be selected instead");
	}




}
