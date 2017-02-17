package beast.evolution.likelihood;

import org.junit.Test;

import beast.core.MCMC;
import beast.util.XMLParser;
import beast.util.XMLParserException;
import junit.framework.TestCase;

public class AlmostConstantAscertainmentTreeLikelihoodTest extends TestCase {

	@Test
	public void testBinaryCase() throws XMLParserException {
		// Start likelihood: -13.563960708066837 
		String xml = "<beast namespace='beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood' version='2.4'>\n" + 
				"\n" + 
				"<data id='test' dataType='binary'>\n" + 
				"    <sequence id='seq_0' taxon='A' totalcount='2' value='01000'/>\n" + 
				"    <sequence id='seq_1' taxon='B' totalcount='2' value='00100'/>\n" + 
				"    <sequence id='seq_2' taxon='C' totalcount='2' value='00010'/>\n" + 
				"    <sequence id='seq_3' taxon='D' totalcount='2' value='00001'/>\n" + 
				"</data>\n" + 
				"\n" + 
				"<run id='mcmc' spec='MCMC' chainLength='1'>\n" + 
				"    <state id='state' storeEvery='5000'>\n" + 
				"	    <tree name='stateNode' id='Tree.t:tree' IsLabelledNewick='true' spec='beast.util.TreeParser' newick='((A:0.25,B:0.25):0.25,(C:0.25,D:0.25):0.25);'>\n" + 
				"            <taxa idref='test'/>\n" + 
				"        </tree>\n" + 
				"    </state>\n" + 
				"\n" + 
				"\n" + 
				"            <distribution id='treeLikelihood.test' spec='TreeLikelihood' tree='@Tree.t:tree'>\n" + 
				"                <data idref='test'/>\n" + 
				"                <siteModel id='SiteModel.s:test' spec='SiteModel' gammaCategoryCount='1' mutationRate='1.0' shape='1.0'>\n" + 
				"                    <parameter id='proportionInvariant.s:test' estimate='false' lower='0.0' name='proportionInvariant' upper='1.0'>0.0</parameter>\n" + 
				"                    <substModel id='CTMC.s:test' spec='GeneralSubstitutionModel'>\n" + 
				"                        <parameter id='rates.s:test' dimension='2' estimate='false' lower='0.0' name='rates'>1.0 1.0</parameter>\n" + 
				"                        <frequencies id='estimatedFreqs.s:test' spec='Frequencies' frequencies='0.5 0.5'/>\n" + 
				"                    </substModel>\n" + 
				"                </siteModel>\n" + 
				"                <branchRateModel id='StrictClock.c:clock' spec='beast.evolution.branchratemodel.StrictClockModel' clock.rate='1.0'/>\n" + 
				"            </distribution>\n" + 
				"\n" + 
				"    <logger id='screenlog' logEvery='1000'>\n" + 
				"        <log idref='treeLikelihood.test'/>\n" + 
				"    </logger>\n" + 
				"</run>\n" + 
				"\n" + 
				"</beast>\n"; 
		XMLParser parser = new XMLParser();
		MCMC mcmc = (MCMC) parser.parseFragment(xml, true);
		GenericTreeLikelihood treelikelihood = (GenericTreeLikelihood) mcmc.posteriorInput.get();
		double trueP = mcmc.robustlyCalcPosterior(treelikelihood);
		
		AlmostConstantAscertainedTreeLikelihood acaLikelihood = new AlmostConstantAscertainedTreeLikelihood();
		acaLikelihood.initByName("treelikelihood", treelikelihood, "maxDeviation", 1, 
				"data", treelikelihood.dataInput.get(), "tree", treelikelihood.treeInput.get(), 
				"siteModel", treelikelihood.siteModelInput.get(),
				"branchRateModel", treelikelihood.branchRateModelInput.get());
		acaLikelihood.requiresRecalculation();
		acaLikelihood.calculateLogP();
		double aca = acaLikelihood.calcAscertainmentCorrection();
		assertEquals(trueP, aca, 1e-10);
	}

	@Test
	public void testTernaryCase() {
		
	}

	@Test
	public void testFourStateCase() {
		
	}
}
