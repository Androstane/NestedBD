package beast.evolution.likelihood;
import beast.core.Description;
import beast.core.Input;
import beast.core.State;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.alignment.Alignment;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.branchratemodel.StrictClockModel;
import beast.evolution.likelihood.TreeLikelihood.Scaling;
import beast.evolution.sitemodel.SiteModel;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;

import java.util.*;

public class DiploidOriginLikelihood extends TreeLikelihood {
	public Input<RealParameter> origtime = new Input<RealParameter>("origtime", "time between diploid and tree root");
	public Input<RealParameter> nstates = new Input<RealParameter>("nstates", "same as what in BD model", Validate.REQUIRED);
    //protected int nrofStates = dataInput.get().getMaxStateCount();
    //protected int nrofPattern = dataInput.get().getPatternCount();
	//protected int nrofStates = 30;
	protected int nrofStates;
    protected static double bd_rate = 0.01;
    protected boolean recalc = true;
    
    /* override to avoid beagle*/
    @Override
    public void initAndValidate() {
    	// System.out.println("Likelihood class");
    	// System.out.println(nstates.get().getValue());
    	nrofStates = (int)Math.round(nstates.get().getValue());
        // sanity check: alignment should have same #taxa as tree
        if (dataInput.get().getTaxonCount() != treeInput.get().getLeafNodeCount()) {
            throw new IllegalArgumentException("The number of nodes in the tree does not match the number of sequences");
        }
        int nodeCount = treeInput.get().getNodeCount();
        if (!(siteModelInput.get() instanceof SiteModel.Base)) {
        	throw new IllegalArgumentException("siteModel input should be of type SiteModel.Base");
        }
        // System.out.println(nstates.get().getValue());
        if (nstates.get() == null) {
            throw new IllegalArgumentException("number of states to consider is required");
        }
        m_siteModel = (SiteModel.Base) siteModelInput.get();
        //System.out.println(m_siteModel);
        m_siteModel.setDataType(dataInput.get().getDataType());
        //System.out.println("DataType");
        //System.out.println(dataInput.get().getDataType());
        substitutionModel = m_siteModel.substModelInput.get();
        //System.out.println(substitutionModel);
        if (branchRateModelInput.get() != null) {
            branchRateModel = branchRateModelInput.get();
        } else {
            branchRateModel = new StrictClockModel();
        }
        m_branchLengths = new double[nodeCount];
        storedBranchLengths = new double[nodeCount];
        int stateCount = (int)Math.round(nstates.get().getValue());
        //int stateCount = dataInput.get().getMaxStateCount();
        //int stateCount = 30;
        int patterns = dataInput.get().getPatternCount();
        if (stateCount == 4) {
            likelihoodCore = new BeerLikelihoodCore4();
        } else {
            likelihoodCore = new BeerLikelihoodCore(stateCount);
        }
                String className = getClass().getSimpleName();

        Alignment alignment = dataInput.get();

        Log.info.println(className + "(" + getID() + ") uses " + likelihoodCore.getClass().getSimpleName());
        Log.info.println("  " + alignment.toString(true));
        // print startup messages via Log.print*

        proportionInvariant = m_siteModel.getProportionInvariant();
        m_siteModel.setPropInvariantIsCategory(false);
        if (proportionInvariant > 0) {
            calcConstantPatternIndices(patterns, stateCount);
        }

        initCore();

        patternLogLikelihoods = new double[patterns];
        m_fRootPartials = new double[patterns * stateCount];
        matrixSize = (stateCount + 1) * (stateCount + 1);
        probabilities = new double[(stateCount+1) * (stateCount+1)];
        Arrays.fill(probabilities, 1.0);
        if (dataInput.get().isAscertained) {
            useAscertainedSitePatterns = true;
        }
        
        
    }
    @Override
    protected boolean requiresRecalculation() {
    	//final TreeInterface tree = treeInput.get();
        //System.out.println(origtime.get().somethingIsDirty());
        //System.out.println(traverse(tree.getRoot()) != Tree.IS_CLEAN);
        return (origtime.get().somethingIsDirty() | recalc);
        	
        //return true;
    }
    
    
   
    
    /**
     * Calculate the log likelihood of the current state.
     *
     * @return the log likelihood.
     */
    double m_fScale = 1.01;
    int m_nScale = 0;
    int X = 100;

    @Override
    public double calculateLogP() {
        if (beagle != null) {
            logP = beagle.calculateLogP();
            return logP;
        }
        final TreeInterface tree = treeInput.get();

        try {
        	if (traverse(tree.getRoot()) != Tree.IS_CLEAN )
        		calcLogP();
        }
        catch (ArithmeticException e) {
        	System.out.println("exception occured");
        	return Double.NEGATIVE_INFINITY;
        }
        m_nScale++;
        if (logP > 0 || (likelihoodCore.getUseScaling() && m_nScale > X)) {
//            System.err.println("Switch off scaling");
//            m_likelihoodCore.setUseScaling(1.0);
//            m_likelihoodCore.unstore();
//            m_nHasDirt = Tree.IS_FILTHY;
//            X *= 2;
//            traverse(tree.getRoot());
//            calcLogP();
//            return logP;
        } else if (logP == Double.NEGATIVE_INFINITY && m_fScale < 10 && !scaling.get().equals(Scaling.none)) { // && !m_likelihoodCore.getUseScaling()) {
            m_nScale = 0;
            m_fScale *= 1.01;
            Log.warning.println("Turning on scaling to prevent numeric instability " + m_fScale);
            likelihoodCore.setUseScaling(m_fScale);
            likelihoodCore.unstore();
            hasDirt = Tree.IS_FILTHY;
            traverse(tree.getRoot());
            calcLogP();
            return logP;
        }
        return logP;
    }

    
    void calcLogP() {
        logP = 0.0;
        if (useAscertainedSitePatterns) {
            final double ascertainmentCorrection = dataInput.get().getAscertainmentCorrection(patternLogLikelihoods);
            for (int i = 0; i < dataInput.get().getPatternCount(); i++) {
                logP += (patternLogLikelihoods[i] - ascertainmentCorrection) * dataInput.get().getPatternWeight(i);
            }
        } else {
            for (int i = 0; i< dataInput.get().getPatternCount(); i++) {
                logP += patternLogLikelihoods[i] * dataInput.get().getPatternWeight(i);
                //System.out.println(patternLogLikelihoods[i]);
            }
        }
    }
    
    /**
     * set leaf partials in likelihood core *
     */

    // for testing

    
    /* Assumes there IS a branch rate model as opposed to traverse() 
     * compute the diploid origin likelihood instead of likelihood
     */
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
                //if (node.getNr() == 29) {
                	//System.out.println(node.getNr() + " " + Arrays.toString(probabilities));
                //}
                
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
                	//System.out.println(nodeIndex);
                    likelihoodCore.calculatePartials(childNum1, childNum2, nodeIndex);
                } else {
                    throw new RuntimeException("Error TreeLikelihood 201: Site categories not supported");
                    //m_pLikelihoodCore->calculatePartials(childNum1, childNum2, nodeNum, siteCategories);
                }

                if (node.isRoot()) {
                	//System.out.println("root");
                    // No parent this is the root of the beast.tree -
                    // calculate the pattern likelihoods
                    final double[] proportions = m_siteModel.getCategoryProportions(node);
                    //System.out.println(node.getNr() + " Proportions" + Arrays.toString(proportions));
                    likelihoodCore.integratePartials(node.getNr(), proportions, m_fRootPartials);
                    //System.out.println(node.getNr() + " m_RootPartials " + Arrays.toString(m_fRootPartials));
                    if (constantPattern != null) { // && !SiteModel.g_bUseOriginal) {
                        proportionInvariant = m_siteModel.getProportionInvariant();
                        // some portion of sites is invariant, so adjust root partials for this
                        for (final int i : constantPattern) {
                            m_fRootPartials[i] += proportionInvariant;
                        }
                    }

                	double distance = origtime.get().getValue();
                	//System.out.println("DD" + distance);
                	substitutionModel.getTransitionProbabilities(node, distance, 0.0,  1.0, probabilities);
                	//System.out.println("probabilitys" + Arrays.toString(probabilities));
                    DipOrigLikelihoods(m_fRootPartials, probabilities, patternLogLikelihoods);
                }

            }
        }
        if (update == 0) {
        	recalc = false;
        }
        return update;
    } // traverseWithBRM
    

    
    
    public void DipOrigLikelihoods(double[] partials, double[] transition_prob, double[] outLogLikelihoods) {
    	
        int v = 0;
        //System.out.println("Where");
        //double distance = origtime.get().getValue();
        //System.out.println("DDD" + distance);
    	//System.out.println(dataInput.get().getPatternCount());
        for (int k = 0; k <dataInput.get().getPatternCount(); k++) {
            //double max_prob =  Double.NEGATIVE_INFINITY;
            double max_prob =  0;
            for (int i = 0; i < nrofStates; i++) {
            	//System.out.println(DipOrigProb(i, distance));
            	//System.out.println(transition_prob[62+i]);
            	//max_prob = Double.max(max_prob, partials[v]* transition_prob[2 * nrofStates +i]);
            	max_prob = max_prob + partials[v]* transition_prob[2 * nrofStates +i];
               if (Double.isNaN(partials[v])){
                	System.out.println(k);
                }
            	//if (k == dataInput.get().getPatternCount()-1) {
            		//System.out.println(partials[v]);}
                v++;
            }

            outLogLikelihoods[k] = Math.log(max_prob) + likelihoodCore.getLogScalingFactor(k);
        }
    }
}
