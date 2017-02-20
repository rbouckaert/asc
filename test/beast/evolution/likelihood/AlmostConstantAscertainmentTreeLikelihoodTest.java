package beast.evolution.likelihood;

import org.junit.Test;

import beast.core.MCMC;
import beast.util.XMLParser;
import beast.util.XMLParserException;
import junit.framework.TestCase;

public class AlmostConstantAscertainmentTreeLikelihoodTest extends TestCase {

	@Test
	public void testBinaryCase() throws XMLParserException {
		String xml = "<beast namespace='beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood' version='2.4'>\n" + 
				"\n" + 
				"<data id='test' dataType='binary'>\n" + 
// for maxDeviation = 0
				"    <sequence id='seq_0' taxon='A' totalcount='2' value='10'/>\n" + 
				"    <sequence id='seq_1' taxon='B' totalcount='2' value='10'/>\n" + 
				"    <sequence id='seq_2' taxon='C' totalcount='2' value='10'/>\n" + 
				"    <sequence id='seq_3' taxon='D' totalcount='2' value='10'/>\n" + 
				"    <sequence id='seq_4' taxon='E' totalcount='2' value='10'/>\n" + 

// for nrOfTaxa = 4, maxDeviation = 1
//				"    <sequence id='seq_0' taxon='A' totalcount='2' value='1010000111'/>\n" + 
//				"    <sequence id='seq_1' taxon='B' totalcount='2' value='1001001011'/>\n" + 
//				"    <sequence id='seq_2' taxon='C' totalcount='2' value='1000101101'/>\n" + 
//				"    <sequence id='seq_3' taxon='D' totalcount='2' value='1000011110'/>\n" + 

//for nrOfTaxa = 5, maxDeviation = 1
//				"    <sequence id='seq_0' taxon='A' totalcount='2' value='101000001111'/>\n" + 
//				"    <sequence id='seq_1' taxon='B' totalcount='2' value='100100010111'/>\n" + 
//				"    <sequence id='seq_2' taxon='C' totalcount='2' value='100010011011'/>\n" + 
//				"    <sequence id='seq_3' taxon='D' totalcount='2' value='100001011101'/>\n" + 
//				"    <sequence id='seq_4' taxon='E' totalcount='2' value='100000111110'/>\n" + 
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
	
	@Test
	public void testTernaryCase() throws XMLParserException {
		String xml = "<beast namespace='beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood' version='2.4'>\n" + 
				"\n" + 
				"<data id='test'>\n" +
				"    <userDataType spec='beast.evolution.datatype.UserDataType' states='3' codelength='1' codeMap='0=0,1=1,2=2'/>\n" + 
// for maxDeviation = 0
//				"    <sequence id='seq_0' taxon='A' totalcount='3' value='210'/>\n" + 
//				"    <sequence id='seq_1' taxon='B' totalcount='3' value='210'/>\n" + 
//				"    <sequence id='seq_2' taxon='C' totalcount='3' value='210'/>\n" + 
//				"    <sequence id='seq_3' taxon='D' totalcount='3' value='210'/>\n" + 
//				"    <sequence id='seq_4' taxon='E' totalcount='3' value='210'/>\n" + 

//for nrOfTaxa = 5, maxDeviation = 1
//				"    <sequence id='seq_0' taxon='A' totalcount='3' value='210 2100000000 2011111111 0122222222'/>\n" + 
//				"    <sequence id='seq_1' taxon='B' totalcount='3' value='210 0021000000 1120111111 2201222222'/>\n" + 
//				"    <sequence id='seq_2' taxon='C' totalcount='3' value='210 0000210000 1111201111 2222012222'/>\n" + 
//				"    <sequence id='seq_3' taxon='D' totalcount='3' value='210 0000002100 1111112011 2222220122'/>\n" + 
//				"    <sequence id='seq_4' taxon='E' totalcount='3' value='210 0000000021 1111111120 2222222201'/>\n" + 
//for nrOfTaxa = 5, maxDeviation = 2
				"<sequence id='seq_0' taxon='A' totalcount='3' value='000000000000000000000000000000000000000000000000000111111111111111111111111111111111111111111111111111222222222222222222222222222222222222222222222222222'/>\n" +
				"<sequence id='seq_1' taxon='B' totalcount='3' value='000000000000000000000111111111111111222222222222222000000000000000111111111111111111111222222222222222000000000000000111111111111111222222222222222222222'/>\n" +
				"<sequence id='seq_2' taxon='C' totalcount='3' value='000000000111111222222000000111111222000000111222222000000111111222000000111111111222222000111111222222000000111222222000111111222222000000111111222222222'/>\n" +
				"<sequence id='seq_3' taxon='D' totalcount='3' value='000111222000112000122000112001112012000122012001222000112001112012001112000111222011122012011122011222000122012001222012011122011222001222011222000111222'/>\n" +
				"<sequence id='seq_4' taxon='E' totalcount='3' value='012012012012010012002012010010121012012002012022012012010010121012010121012012012101212012101212212012012002012022012012101212212012022012212012012012012'/>\n" +
				"</data>\n" + 
				"\n" + 
				"<run id='mcmc' spec='MCMC' chainLength='1'>\n" + 
				"    <state id='state' storeEvery='5000'>\n" + 
				"	    <tree name='stateNode' id='Tree.t:tree' IsLabelledNewick='true' spec='beast.util.TreeParser' " +
				"newick='((A:0.25,B:0.25):0.25,(C:0.2,(D:0.1,E:0.1):0.1):0.3);'>\n" + 
//				"newick='(A:0.25,B:0.25);'>\n" + 
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
				"                        <frequencies id='estimatedFreqs.s:test' spec='Frequencies' frequencies='0.5 0.3 0.2'/>\n" + 
//				"                        <frequencies id='estimatedFreqs.s:test' spec='Frequencies' frequencies='0.33333333 0.33333333 0.33333333'/>\n" + 
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
			// 5 taxa, equal freqs
			// assertEquals(-1.4046981265141947, aca, 1e-10);
			// 5 taxa, unequal freqs
			assertEquals(-1.3326164017084265, aca, 1e-10);
		}
		
		{   // include singletons and constant sites
			AlmostConstantAscertainedTreeLikelihood acaLikelihood = new AlmostConstantAscertainedTreeLikelihood();
			acaLikelihood.initByName("treelikelihood", treelikelihood, "maxDeviation", 1, 
					"data", treelikelihood.dataInput.get(), "tree", treelikelihood.treeInput.get(), 
					"siteModel", treelikelihood.siteModelInput.get(),
					"branchRateModel", treelikelihood.branchRateModelInput.get());
			acaLikelihood.requiresRecalculation();
			acaLikelihood.calculateLogP();
			double aca = acaLikelihood.calcAscertainmentCorrection();
			// 5 taxa, equal freqs
			// assertEquals(-0.629694078440168, aca, 1e-10);
			// 5 taxa, equal unfreqs
			assertEquals(-0.588342183595726, aca, 1e-10);			
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
			assertEquals(-0.07992305457106781, aca, 1e-10);
		}
	}

	
	@Test
	public void testFourStateCase() throws XMLParserException {
		String xml = "<beast namespace='beast.core:beast.evolution.alignment:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.nuc:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood' version='2.4'>\n" + 
				"\n" + 
				"<data id='test'>\n" +
				"    <userDataType spec='beast.evolution.datatype.UserDataType' states='4' codelength='1' codeMap='0=0,1=1,2=2,3=3'/>\n" + 
// for maxDeviation = 0
//				"    <sequence id='seq_0' taxon='A' totalcount='4' value='3210'/>\n" + 
//				"    <sequence id='seq_1' taxon='B' totalcount='4' value='3210'/>\n" + 
//				"    <sequence id='seq_2' taxon='C' totalcount='4' value='3210'/>\n" + 
//				"    <sequence id='seq_3' taxon='D' totalcount='4' value='3210'/>\n" + 
//				"    <sequence id='seq_4' taxon='E' totalcount='4' value='3210'/>\n" + 

//for nrOfTaxa = 5, maxDeviation = 1
//"    <sequence id='seq_0' taxon='A' totalcount='4' value='0000000000000000111111111111111122222222222222223333333333333333'/>\n" +
//"    <sequence id='seq_1' taxon='B' totalcount='4' value='0000000000112233001111111111223300112222222222330011223333333333'/>\n" +
//"    <sequence id='seq_2' taxon='C' totalcount='4' value='0000000123010203010111111123121302120122222223230313230123333333'/>\n" +
//"    <sequence id='seq_3' taxon='D' totalcount='4' value='0000123000010203011011112311121302122201222232230313233330123333'/>\n" +
//"    <sequence id='seq_4' taxon='E' totalcount='4' value='0123000000010203011101231111121302122222012322230313233333330123'/>\n" +
	
//for nrOfTaxa = 5, maxDeviation = 2
"    <sequence id='seq_0' taxon='A' totalcount='4' value='0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111122222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222223333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333333'/>\n" +
"    <sequence id='seq_1' taxon='B' totalcount='4' value='0000000000000000000000000000000000000000111111111111111111111122222222222222222222223333333333333333333333000000000000000000000011111111111111111111111111111111111111112222222222222222222222333333333333333333333300000000000000000000001111111111111111111111222222222222222222222222222222222222222233333333333333333333330000000000000000000000111111111111111111111122222222222222222222223333333333333333333333333333333333333333'/>\n" +
"    <sequence id='seq_2' taxon='C' totalcount='4' value='0000000000000000111111112222222233333333000000001111111122233300000000111222222223330000000011122233333333000000001111111122233300000000111111111111111122222222333333330001111111122222222333000111111112223333333300000000111222222223330001111111122222222333000000001111111122222222222222223333333300011122222222333333330000000011122233333333000111111112223333333300011122222222333333330000000011111111222222223333333333333333'/>\n" +
"    <sequence id='seq_3' taxon='D' totalcount='4' value='0000111122223333000011230000122300001233000011230011112301201300001223012001222230230000123301302300123333000011230011112301201300111123000011112222333301111223011112330120111122301122223123013011112331230112333300001223012001222230230120111122301122223123001222230112222300001111222233330122223302312301222233012233330000123301302300123333013011112331230112333302312301222233012233330012333301123333012233330000111122223333'/>\n" +
"    <sequence id='seq_4' taxon='E' totalcount='4' value='0123012301230123012301000123002001230003012301000101231101201301230020012022012320230123000301302303330123012301000101231101201301012311012301230123012310123121101231130121012312121201232123013101231131233133012301230020012022012320230121012312121201232123022012322120123201230123012301232201232302312322012323332301230123000301302303330123013101231131233133012302312322012323332301230333012331330123332301230123012301230123'/>\n" +
				"</data>\n" + 
				"\n" + 
				"<run id='mcmc' spec='MCMC' chainLength='1'>\n" + 
				"    <state id='state' storeEvery='5000'>\n" + 
				"	    <tree name='stateNode' id='Tree.t:tree' IsLabelledNewick='true' spec='beast.util.TreeParser' " +
				"newick='((A:0.25,B:0.25):0.25,(C:0.2,(D:0.1,E:0.1):0.1):0.3);'>\n" + 
//				"newick='(A:0.25,B:0.25);'>\n" + 
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
				"                        <parameter id='rates.s:test' dimension='12' estimate='false' lower='0.0' name='rates'>1.0</parameter>\n" + 
				"                        <frequencies id='estimatedFreqs.s:test' spec='Frequencies' frequencies='0.5 0.3 0.1 0.1'/>\n" + 
//				"                        <frequencies id='estimatedFreqs.s:test' spec='Frequencies' frequencies='0.33333333 0.33333333 0.33333333'/>\n" + 
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
			// 5 taxa, unequal freqs
			assertEquals(-1.3287452003443339, aca, 1e-10);
		}

		{   // include singletons and constant sites
			AlmostConstantAscertainedTreeLikelihood acaLikelihood = new AlmostConstantAscertainedTreeLikelihood();
			acaLikelihood.initByName("treelikelihood", treelikelihood, "maxDeviation", 1, 
					"data", treelikelihood.dataInput.get(), "tree", treelikelihood.treeInput.get(), 
					"siteModel", treelikelihood.siteModelInput.get(),
					"branchRateModel", treelikelihood.branchRateModelInput.get());
			acaLikelihood.requiresRecalculation();
			acaLikelihood.calculateLogP();
			double aca = acaLikelihood.calcAscertainmentCorrection();
			// 5 taxa, equal unfreqs
			assertEquals(-0.6039289340612545, aca, 1e-10);			
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
			assertEquals(-0.09642213012889077, aca, 1e-10);
		}
	}
}
