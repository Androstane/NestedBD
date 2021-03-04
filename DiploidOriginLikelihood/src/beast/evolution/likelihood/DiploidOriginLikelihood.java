package beast.evolution.likelihood;
import beast.core.Description;
import beast.core.Input;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.alignment.Alignment;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.branchratemodel.StrictClockModel;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;
import java.util.Collections;


import java.util.*;

public class DiploidOriginLikelihood extends TreeLikelihood {
	public Input<RealParameter> origtime = new Input<RealParameter>("origtime", "time between diploid and tree root");
    protected int nrofStates = dataInput.get().getMaxStateCount();
    protected int nrofPattern = dataInput.get().getPatternCount();
    protected static double bd_rate = 1.0;
    
    /* Assumes there IS a branch rate model as opposed to traverse() */
    @Override
    protected int traverse(final Node node) {

        int update = (node.isDirty() | hasDirt);

        final int nodeIndex = node.getNr();

        final double branchRate = branchRateModel.getRateForBranch(node);
        final double branchTime = node.getLength() * branchRate;

        // First update the transition probability matrix(ices) for this branch
        //if (!node.isRoot() && (update != Tree.IS_CLEAN || branchTime != m_StoredBranchLengths[nodeIndex])) {
        if (!node.isRoot() && (update != Tree.IS_CLEAN || branchTime != m_branchLengths[nodeIndex])) {
            m_branchLengths[nodeIndex] = branchTime;
            final Node parent = node.getParent();
            likelihoodCore.setNodeMatrixForUpdate(nodeIndex);
            for (int i = 0; i < m_siteModel.getCategoryCount(); i++) {
                final double jointBranchRate = m_siteModel.getRateForCategory(i, node) * branchRate;
                substitutionModel.getTransitionProbabilities(node, parent.getHeight(), node.getHeight(), jointBranchRate, probabilities);
                //System.out.println(node.getNr() + " " + Arrays.toString(m_fProbabilities));
                likelihoodCore.setNodeMatrix(nodeIndex, i, probabilities);
            }
            update |= Tree.IS_DIRTY;
        }

        // If the node is internal, update the partial likelihoods.
        if (!node.isLeaf()) {

            // Traverse down the two child nodes
            final Node child1 = node.getLeft(); //Two children
            final int update1 = traverse(child1);

            final Node child2 = node.getRight();
            final int update2 = traverse(child2);

            // If either child node was updated then update this node too
            if (update1 != Tree.IS_CLEAN || update2 != Tree.IS_CLEAN) {

                final int childNum1 = child1.getNr();
                final int childNum2 = child2.getNr();

                likelihoodCore.setNodePartialsForUpdate(nodeIndex);
                update |= (update1 | update2);
                if (update >= Tree.IS_FILTHY) {
                    likelihoodCore.setNodeStatesForUpdate(nodeIndex);
                }

                if (m_siteModel.integrateAcrossCategories()) {
                    likelihoodCore.calculatePartials(childNum1, childNum2, nodeIndex);
                } else {
                    throw new RuntimeException("Error TreeLikelihood 201: Site categories not supported");
                    //m_pLikelihoodCore->calculatePartials(childNum1, childNum2, nodeNum, siteCategories);
                }

                if (node.isRoot()) {
                    // No parent this is the root of the beast.tree -
                    // calculate the pattern likelihoods

                	
                    final double[] proportions = m_siteModel.getCategoryProportions(node);
                    likelihoodCore.integratePartials(node.getNr(), proportions, m_fRootPartials);

                    if (constantPattern != null) { // && !SiteModel.g_bUseOriginal) {
                        proportionInvariant = m_siteModel.getProportionInvariant();
                        // some portion of sites is invariant, so adjust root partials for this
                        for (final int i : constantPattern) {
                            m_fRootPartials[i] += proportionInvariant;
                        }
                    }
                    DipOrigLikelihoods(m_fRootPartials, patternLogLikelihoods);
                }

            }
        }
        return update;
    } // traverseWithBRM
    
    protected long binomi(int n, int k) {
        if ((n == k) || (k == 0))
            return 1;
        else
            return binomi(n - 1, k) + binomi(n - 1, k - 1);
    }
	
	protected double bd_prob(int child, int ancestor, double bd_rate, double distance) {
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

    
    protected double DipOrigProb(int child, double distance) {
    	double prob;
    	if (child == 0){
    		prob = Math.pow(distance * bd_rate/(1 + distance * bd_rate), 2);
    	}
    	else {
    		prob = bd_prob(child, 2, bd_rate, distance);
    	}
    	return prob;
    }
    
    
    public void DipOrigLikelihoods(double[] partials, double[] outLogLikelihoods) {
    	double distance = origtime.get().getValue();
        int v = 0;
        for (int k = 0; k < nrofPattern; k++) {
            double max_prob = 0.0;
            for (int i = 0; i < nrofStates; i++) {
            	max_prob = Double.max(max_prob, partials[v]* DipOrigProb(i, distance));
                v++;
            }
            outLogLikelihoods[k] = Math.log(max_prob) + likelihoodCore.getLogScalingFactor(k);
        }
    }
}
