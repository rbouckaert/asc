package beast.evolution.likelihood;

import org.junit.Test;

import beast.core.MCMC;
import beast.util.XMLParser;
import beast.util.XMLParserException;
import junit.framework.TestCase;

public class AlmostConstantAscertainmentTreeLikelihoodTest extends TestCase {
/*
	@Test
	public void testBinaryCase() throws XMLParserException {
		String xml = "<beast namespace='beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood' version='2.4'>\n" + 
				"\n" + 
				"<data id='test' dataType='binary'>\n" + 
// for maxDeviation = 0
//				"    <sequence id='seq_0' taxon='A' totalcount='2' value='10'/>\n" + 
//				"    <sequence id='seq_1' taxon='B' totalcount='2' value='10'/>\n" + 
//				"    <sequence id='seq_2' taxon='C' totalcount='2' value='10'/>\n" + 
//				"    <sequence id='seq_3' taxon='D' totalcount='2' value='10'/>\n" + 
//				"    <sequence id='seq_4' taxon='E' totalcount='2' value='10'/>\n" + 

// for nrOfTaxa = 4, maxDeviation = 1
//				"    <sequence id='seq_0' taxon='A' totalcount='2' value='1010000111'/>\n" + 
//				"    <sequence id='seq_1' taxon='B' totalcount='2' value='1001001011'/>\n" + 
//				"    <sequence id='seq_2' taxon='C' totalcount='2' value='1000101101'/>\n" + 
//				"    <sequence id='seq_3' taxon='D' totalcount='2' value='1000011110'/>\n" + 

//for nrOfTaxa = 5, maxDeviation = 1
				"    <sequence id='seq_0' taxon='A' totalcount='2' value='101000001111'/>\n" + 
				"    <sequence id='seq_1' taxon='B' totalcount='2' value='100100010111'/>\n" + 
				"    <sequence id='seq_2' taxon='C' totalcount='2' value='100010011011'/>\n" + 
				"    <sequence id='seq_3' taxon='D' totalcount='2' value='100001011101'/>\n" + 
				"    <sequence id='seq_4' taxon='E' totalcount='2' value='100000111110'/>\n" + 
				"</data>\n" + 
				"\n" + 
				"<run id='mcmc' spec='MCMC' chainLength='1'>\n" + 
				"    <state id='state' storeEvery='5000'>\n" + 
				"	    <tree name='stateNode' id='Tree.t:tree' IsLabelledNewick='true' spec='beast.util.TreeParser' newick='((A:0.25,B:0.25):0.25,(C:0.2,(D:0.1,E:0.1):0.1):0.3);'>\n" + 
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
				"                        <frequencies id='estimatedFreqs.s:test' spec='Frequencies' frequencies='0.6 0.4'/>\n" + 
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
		System.setProperty("java.only", "true");
		MCMC mcmc = (MCMC) parser.parseFragment(xml, true);
		TreeLikelihood treelikelihood = (TreeLikelihood) mcmc.posteriorInput.get();
		mcmc.robustlyCalcPosterior(treelikelihood);
		double [] p = treelikelihood.patternLogLikelihoods;
		double trueP = 0;
		for (double d : p) {
			trueP += Math.exp(d);
		}
		trueP = Math.log(trueP);
		
		{   // constant sites only
			AlmostConstantAscertainedTreeLikelihood acaLikelihood = new AlmostConstantAscertainedTreeLikelihood();
			acaLikelihood.initByName("treelikelihood", treelikelihood, "maxDeviation", 0, 
					"data", treelikelihood.dataInput.get(), "tree", treelikelihood.treeInput.get(), 
					"siteModel", treelikelihood.siteModelInput.get(),
					"branchRateModel", treelikelihood.branchRateModelInput.get());
			acaLikelihood.requiresRecalculation();
			acaLikelihood.calculateLogP();
			double aca = acaLikelihood.calcAscertainmentCorrection();
			// 4 taxa case, equal freqs: 
			// assertEquals(-1.1536566517135067, aca, 1e-10);
			// 5 taxa case, equal freqs
			// assertEquals(-1.2668460989842023, aca, 1e-10);
			// 5 taxa case, unequal freqs
			assertEquals(-1.2129088256740665, aca, 1e-10);
		}
		
		{   // constant sites + singletons
			AlmostConstantAscertainedTreeLikelihood acaLikelihood = new AlmostConstantAscertainedTreeLikelihood();
			acaLikelihood.initByName("treelikelihood", treelikelihood, "maxDeviation", 1, 
					"data", treelikelihood.dataInput.get(), "tree", treelikelihood.treeInput.get(), 
					"siteModel", treelikelihood.siteModelInput.get(),
					"branchRateModel", treelikelihood.branchRateModelInput.get());
			acaLikelihood.requiresRecalculation();
			acaLikelihood.calculateLogP();
			double aca = acaLikelihood.calcAscertainmentCorrection();
			// 4 taxa case, equal freqs: 
			// assertEquals(-0.31083933173222555, aca, 1e-10);
			// 5 taxa case, equal freqs
			// assertEquals(-0.4719578439004392, aca, 1e-10);
			// 5 taxa case, unequal freqs
			assertEquals(-0.45231882251944283, aca, 1e-10);
		}
		
		{   // constant sites + singletons + doubles = all sites for nrOfTaxa = 5
			AlmostConstantAscertainedTreeLikelihood acaLikelihood = new AlmostConstantAscertainedTreeLikelihood();
			acaLikelihood.initByName("treelikelihood", treelikelihood, "maxDeviation", 2, 
					"data", treelikelihood.dataInput.get(), "tree", treelikelihood.treeInput.get(), 
					"siteModel", treelikelihood.siteModelInput.get(),
					"branchRateModel", treelikelihood.branchRateModelInput.get());
			acaLikelihood.requiresRecalculation();
			acaLikelihood.calculateLogP();
			double aca = acaLikelihood.calcAscertainmentCorrection();
			// 4 taxa case: assertEquals(-0.31083933173222555, aca, 1e-10);
			assertEquals(0, aca, 1e-10);
		}
	}
*/
	
	@Test
	public void testTernaryCase() throws XMLParserException {
		String xml = "<beast namespace='beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood' version='2.4'>\n" + 
				"\n" + 
				"<data id='test'>\n" +
				"    <userDataType spec='beast.evolution.datatype.UserDataType' states='3' codelength='1' codeMap='0=0,1=1,2=2'/>\n" + 
// for maxDeviation = 0
				"    <sequence id='seq_0' taxon='A' totalcount='3' value='210'/>\n" + 
				"    <sequence id='seq_1' taxon='B' totalcount='3' value='210'/>\n" + 
//				"    <sequence id='seq_2' taxon='C' totalcount='3' value='210'/>\n" + 
//				"    <sequence id='seq_3' taxon='D' totalcount='3' value='210'/>\n" + 
//				"    <sequence id='seq_4' taxon='E' totalcount='3' value='210'/>\n" + 

// for nrOfTaxa = 4, maxDeviation = 1
//				"    <sequence id='seq_0' taxon='A' totalcount='3' value='1010000111'/>\n" + 
//				"    <sequence id='seq_1' taxon='B' totalcount='3' value='1001001011'/>\n" + 
//				"    <sequence id='seq_2' taxon='C' totalcount='3' value='1000101101'/>\n" + 
//				"    <sequence id='seq_3' taxon='D' totalcount='3' value='1000011110'/>\n" + 

//for nrOfTaxa = 5, maxDeviation = 1
//				"    <sequence id='seq_0' taxon='A' totalcount='3' value='101000001111'/>\n" + 
//				"    <sequence id='seq_1' taxon='B' totalcount='3' value='100100010111'/>\n" + 
//				"    <sequence id='seq_2' taxon='C' totalcount='3' value='100010011011'/>\n" + 
//				"    <sequence id='seq_3' taxon='D' totalcount='3' value='100001011101'/>\n" + 
//				"    <sequence id='seq_4' taxon='E' totalcount='3' value='100000111110'/>\n" + 
				"</data>\n" + 
				"\n" + 
				"<run id='mcmc' spec='MCMC' chainLength='1'>\n" + 
				"    <state id='state' storeEvery='5000'>\n" + 
				"	    <tree name='stateNode' id='Tree.t:tree' IsLabelledNewick='true' spec='beast.util.TreeParser' " +
//				"newick='((A:0.25,B:0.25):0.25,(C:0.2,(D:0.1,E:0.1):0.1):0.3);'>\n" + 
				"newick='(A:0.25,B:0.25);'>\n" + 
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
				"                        <parameter id='rates.s:test' dimension='2' estimate='false' lower='0.0' name='rates'>1.0 1.0 1.0 1.0 1.0 1.0</parameter>\n" + 
//				"                        <frequencies id='estimatedFreqs.s:test' spec='Frequencies' frequencies='0.5 0.3 0.2'/>\n" + 
				"                        <frequencies id='estimatedFreqs.s:test' spec='Frequencies' frequencies='0.33333333 0.33333333 0.33333333'/>\n" + 
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
		System.setProperty("java.only", "true");
		MCMC mcmc = (MCMC) parser.parseFragment(xml, true);
		TreeLikelihood treelikelihood = (TreeLikelihood) mcmc.posteriorInput.get();
		mcmc.robustlyCalcPosterior(treelikelihood);
		double [] p = treelikelihood.patternLogLikelihoods;
		double trueP = 0;
		for (double d : p) {
			trueP += Math.exp(d);
		}
		trueP = Math.log(trueP);

		{   // constant sites only
			AlmostConstantAscertainedTreeLikelihood acaLikelihood = new AlmostConstantAscertainedTreeLikelihood();
			acaLikelihood.initByName("treelikelihood", treelikelihood, "maxDeviation", 0, 
					"data", treelikelihood.dataInput.get(), "tree", treelikelihood.treeInput.get(), 
					"siteModel", treelikelihood.siteModelInput.get(),
					"branchRateModel", treelikelihood.branchRateModelInput.get());
			acaLikelihood.requiresRecalculation();
			acaLikelihood.calculateLogP();
			double aca = acaLikelihood.calcAscertainmentCorrection();
			//assertEquals(-1.3326164017084265, aca, 1e-10);
			//assertEquals(-0.4203872506896104, aca, 1e-10);
			assertEquals(-0.43348755548863116, aca, 1e-10);
			
			
		}
		
	}

	@Test
	public void testFourStateCase() {
		
	}
}
