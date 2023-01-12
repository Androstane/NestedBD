package beast.evolution.likelihood;

import java.lang.reflect.Array;
import java.util.Arrays;
import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.evolution.alignment.Alignment;
import beast.evolution.datatype.IntegerData;
import beast.evolution.tree.Node;
import beast.evolution.errormodel.ErrorModel;
import beast.evolution.errormodel.poissonErrorModel;
import beast.evolution.errormodel.readcountErrorModel;
@Description("Tree likelihood calculation using DiploidOriginLikelihood with error models")
public class DiploidOriginLikelihoodWithError extends DiploidOriginLikelihood {

	final public Input<beast.evolution.errormodel.ErrorModel> errorModelInput = new Input<>("errorModel", "error model to use for leaf partials");
	final public Input<String> weight = new Input<>("weight", "weight for segment, required for readcountErrorModel");
	protected ErrorModel errorModel;
	
	protected boolean useTipLikelihoods = true;

	@Override
	public void initAndValidate() {
		// get error model
		errorModel = errorModelInput.get();
		if (errorModel == null){
			throw new IllegalArgumentException("Error Model is required");
		}	
		super.m_useTipLikelihoods.setValue(useTipLikelihoods, this);
		super.implementationInput.setValue("beast.evolution.likelihood.DiploidOriginLikelihood", this);
		super.initAndValidate();
	}

	protected int getTaxonIndex(String taxon, Alignment data) {
		int taxonIndex = data.getTaxonIndex(taxon);
		if (taxonIndex == -1) {
			if (taxon.startsWith("'") || taxon.startsWith("\"")) {
				taxonIndex = data.getTaxonIndex(taxon.substring(1, taxon.length() - 1));
			}
			if (taxonIndex == -1) {
				throw new RuntimeException("Could not find sequence " + taxon + " in the alignment");
			}
		}
		return taxonIndex;
	}

	protected double[] getLeafPartials(Node node) {
		Alignment data = dataInput.get();
		int nrOfStates = (int)Math.round(nstates.get().getValue());
		int nrOfPatterns = data.getPatternCount();
		int nrOfSites = data.getSiteCount();
		int t = getTaxonIndex(node.getID(), data); // taxon index
		//System.out.print();
		//System.out.println(node.getID());
		int totalread = data.getStateCounts().get(t);
		int weight_mode = 0;
		//System.out.println(TotalCount.get(t));
		double[] partials = new double[nrOfPatterns * nrOfStates];
		double pattern_weight = 0;
		int counter = 0;
		double [] w = new double[nrOfPatterns];
		if (errorModel instanceof readcountErrorModel || errorModel instanceof poissonErrorModel) {
			//equal weights, use only sampled bins 
			if (weight.get() == null){
				for (int i = 0; i < nrOfPatterns; i ++) {
					w[i] = 1.0/nrOfSites * data.getPatternWeight(i);
					//Assume equal weight of each pattern
				}
			}
			else {
				String w_string = weight.get();
				String[] weight_array = w_string.split(",");
				//equal weights, use all bins
				if (weight_array.length == 1) {
					weight_mode = 1;
					for (int i = 0; i < nrOfPatterns; i ++) {
						double total_site = Double.valueOf(weight_array[0]);
						w[i] = 1.0/ total_site * data.getPatternWeight(i);	
						}
				}
				//input weights, use all bins
				else {
					weight_mode = 2;
					for (int i = 0; i < nrOfPatterns; i ++) {
						for (int j = 0; j < data.getPatternWeight(i); j++) {
							pattern_weight += Double.valueOf(weight_array[counter]);
							counter += 1;
						}
						w[i] = pattern_weight;
						//System.out.println(pattern_weight);
						
						pattern_weight = 0;
					}
				}
			}
		}
		//int t = getTaxonIndex(node.getID(), data); // taxon index
		//System.out.print(t);
		int i = 0;
		if (errorModel instanceof readcountErrorModel || errorModel instanceof poissonErrorModel) {
			if (weight_mode == 0 | totalread == -1){
				totalread = 0;
				for (int p = 0; p < nrOfPatterns; p++) {
					totalread += data.getPattern(i, p);
					}
				}
			}
		for (int p = 0; p < nrOfPatterns; p++) {
			int state = data.getPattern(t, p);
			double[] tipLikelihoods;
			if (errorModel instanceof readcountErrorModel || errorModel instanceof poissonErrorModel) {
				tipLikelihoods = errorModel.getProbabilities(state, w[p], totalread);
				//System.out.println(state);
				//System.out.println(Arrays.toString(tipLikelihoods));
			}
			else {
				tipLikelihoods = errorModel.getProbabilities(state);
				
			}
			for (int s = 0; s < nrOfStates; s++) {
				partials[i] = tipLikelihoods[s];
				i++;
			}
		}
		return partials;
	}

	@Override
	protected void setPartials(Node node, int nrOfPatterns) {
		if (node.isLeaf()) {
			//System.out.println(node.getNr());
			double[] partials = getLeafPartials(node);
			likelihoodCore.setNodePartials(node.getNr(), partials);
		} else {
			setPartials(node.getChild(0), nrOfPatterns);
			setPartials(node.getChild(1), nrOfPatterns);
		}
	}

}
