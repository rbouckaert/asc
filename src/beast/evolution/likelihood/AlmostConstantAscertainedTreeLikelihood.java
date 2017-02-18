package beast.evolution.likelihood;

import java.util.Arrays;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.util.Log;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

//TODO: deal with missing data

@Description("Tree likelihood that ascertains for sites being almost-constant, that is, constant or deviate from constant on at most "
		+ "k sites (where k can be specified to be a number between 0 and at most half the number of taxa).")
public class AlmostConstantAscertainedTreeLikelihood extends TreeLikelihood {
	public Input<Integer> baseInput = new Input<>("maxDeviation", "ascertains on all sites that are constants or "
			+ "have at most 'maxDeviation' deviations from constant sites. maxDeviation must be a positive "
			+ "number not more than half the number of taxa.", 1);
	public Input<GenericTreeLikelihood> treeLikelihoodInput = new Input<>("treelikelihood", "tree likelihood object that calculates the uncorrected likelihood of the tree", Validate.REQUIRED);
	
	AlmostConstantAscertainedBeerLikelihoodCore [] ascertainmentCore;
	double [][] m_fAscRootPartials;
	double [][] ascPatternLogLikelihoods;
	int stateCount;
	
	GenericTreeLikelihood treelikelihood;
	
	@Override
	public void initAndValidate() {
		if (baseInput.get() < 0) {
			throw new IllegalArgumentException("maxDeviation cannot be negative, should at least 0");
		}
		if (baseInput.get() * 2 >= treeInput.get().getLeafNodeCount()) {
			throw new IllegalArgumentException("maxDeviation should be smaller than half the number of taxa in the tree");
		}

		stateCount = dataInput.get().getMaxStateCount();
		m_fAscRootPartials = new double[stateCount][(2) * (baseInput.get() + 1)];
		ascPatternLogLikelihoods = new double[stateCount][baseInput.get() + 1];
		treelikelihood = treeLikelihoodInput.get();

		boolean org = !Boolean.valueOf(System.getProperty("java.only")); 
		if (org) {
        	System.setProperty("java.only", "true");
        }

		super.initAndValidate();
		
		if (org) {
        	System.setProperty("java.only", "false");
        }
		
		if (useAscertainedSitePatterns) {
			Log.warning("WARNING: Using both ascertained patterns and ascertainment on almst-constant sites. "
					+ "If these overlap, the likelihood calculation can be incorrect.");
		}

	}

	@Override
	protected void initCore() {
		super.initCore();
		ascertainmentCore = new AlmostConstantAscertainedBeerLikelihoodCore[stateCount];
        final int nodeCount = treeInput.get().getNodeCount();
		for (int i = 0; i < stateCount; i++) {
			ascertainmentCore[i] = new AlmostConstantAscertainedBeerLikelihoodCore(2, baseInput.get());
	        ascertainmentCore[i].initialize(
	                nodeCount,
	                1,
	                m_siteModel.getCategoryCount(),
	                true, m_useAmbiguities.get()
	        );
		}
        setStates(treeInput.get().getRoot());
        hasDirt = Tree.IS_FILTHY;
        final int extNodeCount = nodeCount / 2 + 1;
        final int intNodeCount = nodeCount / 2;
        for (int i = 0; i < intNodeCount; i++) {
    		for (int k = 0; k < stateCount; k++) {
    			ascertainmentCore[k].createNodePartials(extNodeCount + i);
    		}
        }
	}
	
	protected void setStates(Node node) {
        if (node.isLeaf()) {
            int[] states = new int[1];
            for (int i = 0; i < stateCount; i++) {
            	states[0] = 0;
                ascertainmentCore[i].setNodeStates(node.getNr(), states);
            }
        } else {
            setStates(node.getLeft());
            setStates(node.getRight());
        }
    }

	
	@Override
    void calcLogP() {
        logP = treelikelihood.getCurrentLogP();
        double ascertainmentCorrection = calcAscertainmentCorrection();
//        if (useAscertainedSitePatterns) {
//            ascertainmentCorrection += dataInput.get().getAscertainmentCorrection(patternLogLikelihoods);
//        }
        for (int i = 0; i < dataInput.get().getPatternCount(); i++) {
            logP -=  ascertainmentCorrection * dataInput.get().getPatternWeight(i);
        }
    }

	@Override
	int traverse(final Node node) {

        int update = (node.isDirty() | hasDirt);

        final int nodeIndex = node.getNr();

        final double branchRate = branchRateModel.getRateForBranch(node);
        final double branchTime = node.getLength() * branchRate;

        // First update the transition probability matrix(ices) for this branch
        //if (!node.isRoot() && (update != Tree.IS_CLEAN || branchTime != m_StoredBranchLengths[nodeIndex])) {
        if (!node.isRoot() && (update != Tree.IS_CLEAN || branchTime != m_branchLengths[nodeIndex])) {
            m_branchLengths[nodeIndex] = branchTime;
            final Node parent = node.getParent();
            //likelihoodCore.setNodeMatrixForUpdate(nodeIndex);
            for (int k = 0; k < stateCount; k++) {
            	ascertainmentCore[k].setNodeMatrixForUpdate(nodeIndex);
            }
            
            for (int i = 0; i < m_siteModel.getCategoryCount(); i++) {
                final double jointBranchRate = m_siteModel.getRateForCategory(i, node) * branchRate;
                substitutionModel.getTransitionProbabilities(node, parent.getHeight(), node.getHeight(), jointBranchRate, probabilities);
                //System.out.println(node.getNr() + " " + Arrays.toString(m_fProbabilities));
                //likelihoodCore.setNodeMatrix(nodeIndex, i, probabilities);
                
                for (int k = 0; k < stateCount; k++) {
                	ascertainmentCore[k].setNodeMatrix(nodeIndex, i, collapseTransProbs(probabilities, k));
                }
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

                //likelihoodCore.setNodePartialsForUpdate(nodeIndex);
                for (int k = 0; k < stateCount; k++) {
                	ascertainmentCore[k].setNodePartialsForUpdate(nodeIndex);
                }
                update |= (update1 | update2);
                if (update >= Tree.IS_FILTHY) {
                    //likelihoodCore.setNodeStatesForUpdate(nodeIndex);
                    for (int k = 0; k < stateCount; k++) {
                    	ascertainmentCore[k].setNodePartialsForUpdate(nodeIndex);
                    }
                }

                if (m_siteModel.integrateAcrossCategories()) {
                    //likelihoodCore.calculatePartials(childNum1, childNum2, nodeIndex);
                    for (int k = 0; k < stateCount; k++) {
                    	ascertainmentCore[k].calculatePartials(childNum1, childNum2, nodeIndex);
                    }
                } else {
                    throw new RuntimeException("Error TreeLikelihood 201: Site categories not supported");
                    //m_pLikelihoodCore->calculatePartials(childNum1, childNum2, nodeNum, siteCategories);
                }

                if (node.isRoot()) {
                    // No parent this is the root of the beast.tree -
                    // calculate the pattern likelihoods
                    final double[] frequencies = //m_pFreqs.get().
                            substitutionModel.getFrequencies();

                    final double[] proportions = m_siteModel.getCategoryProportions(node);
                    //likelihoodCore.integratePartials(node.getNr(), proportions, m_fRootPartials);
                    for (int k = 0; k < stateCount; k++) {
                    	ascertainmentCore[k].integratePartials(node.getNr(), proportions, m_fAscRootPartials[k]);
                    }

                    if (constantPattern != null) { // && !SiteModel.g_bUseOriginal) {
                        proportionInvariant = m_siteModel.getProportionInvariant();
                        // some portion of sites is invariant, so adjust root partials for this
                        for (final int i : constantPattern) {
                            m_fRootPartials[i] += proportionInvariant;
                        }
                    }

                    //likelihoodCore.calculateLogLikelihoods(m_fRootPartials, frequencies, patternLogLikelihoods);
                    for (int k = 0; k < stateCount; k++) {
                    	ascertainmentCore[k].calculateLogLikelihoods(m_fAscRootPartials[k], collapseFreqs(frequencies, k), ascPatternLogLikelihoods[k]);
                    }
                }

            }
        }
        return update;
    } // traverse
    
	private double[] collapseTransProbs(double[] probabilities, int k) {
		double [] collapsed = new double[4];
        Arrays.fill(collapsed, 1.0);
        collapsed[0] = probabilities[k * stateCount + k];
        collapsed[1] = 1.0 - collapsed[0];
        double p = -collapsed[0];
        for (int i = 0; i < stateCount; i++) {
        	p += probabilities[i*stateCount + k];
        }
        collapsed[2] = p;
        collapsed[3] = 1.0 - p;
        return collapsed;
	}

	private double[] collapseFreqs(double[] frequencies, int k) {
		double [] collapsed = new double[2];
		collapsed[0] = frequencies[k];
		collapsed[1] = 1.0 - frequencies[k];
		return collapsed;
	}

	double calcAscertainmentCorrection() {
		final double maxDeviationPlusOne = baseInput.get() + 1;
		double max = ascPatternLogLikelihoods[0][0];
		for (int i = 0; i < stateCount; i++) {
			for (int k = 0; k < maxDeviationPlusOne; k++) {
				max = Math.max(max, ascPatternLogLikelihoods[i][k]);
			}
		}
		
		double p = 0;
		for (int i = 0; i < stateCount; i++) {
			for (int k = 0; k < maxDeviationPlusOne; k++) {
				p += Math.exp(ascPatternLogLikelihoods[i][k] - max);
			}
		}
		p = max + Math.log(p);
		return p;
	}

}
