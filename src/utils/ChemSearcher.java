package utils;


/**
 * @author Djoumbou Feunang, Yannick
 * 
 *
 */

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.graph.CycleFinder;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.smiles.smarts.SMARTSQueryTool;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

public class ChemSearcher {

	/**
	 * 
	 * @param molecule
	 *            : The molecule of interest
	 * @param smartsPattern
	 *            : The SMARTS pattern to be search for in the molecule
	 * @return : A list containing a list with the indexes of each corresponding
	 *         atom match for each occurrence of the pattern
	 * @throws Exception
	 */

	public List<List<Integer>> getMatchingAtoms(IAtomContainer molecule,
			SMARTSQueryTool smartsPattern) throws Exception {
		//IChemObjectBuilder builder = SilentChemObjectBuilder.getInstance();
		//boolean occurrences = smartsPattern.matches(molecule);
		List<List<Integer>> indexes = smartsPattern.getMatchingAtoms();
//		for (int i = 0; i < indexes.size(); i++) {
//			List<Integer> currentAtomIndexes = indexes.get(i);
//			System.out.println("\nMatch nr." + (i + 1) + ": " + (currentAtomIndexes));
//			for (int j = 0; j < currentAtomIndexes.size(); j++) {
//				System.out.println(currentAtomIndexes.get(j) + ": "
//						+ molecule.getAtom((int) currentAtomIndexes.get(j)));
//			}
//		}
		return indexes;
	}

	/**
	 * 
	 * @return : A HashMap with the SMARTS expressions for functional groups and
	 *         structural patterns
	 * @throws Exception
	 */
	public static LinkedHashMap<String, String> getFingerprintPatterns() throws Exception {
		LinkedHashMap<String, String> queries = new LinkedHashMap<String, String>();
		queries.put("carboxyl", "[#8;A;X2H1,X1-][#6]([#6,#1;A])=O");
		queries.put("hydroxyl", "[#6][OX2H1]");
		queries.put("aromatic", "[*;a]");
		queries.put("sulfuric acid",
				"[SX4](=[OX1])(=[OX1])([$([OX2H]),$([OX1-])])[$([OX2H]),$([OX1-])]");
		queries.put("carbonyl", "[#6]-[#6;X3](-[*,#1;O,C,#1,N,X])=[O;v2X1]");
		queries.put("aldehyde", "[#6;X3H1](-[#6])=[O;v2X1]");
		queries.put("aryl aldehyde", "[#6;X3H1](-[#6;a])=[O;v2X1]");
		queries.put("indole",
				"[#7;R1]-1-[#6;R1]=[#6;R1]-[#6]-2=[#6]-1-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-2");
		queries.put("imidazole", "[c;R1]1[c;R1][n;R1][c;R1][n;R1]1");
		queries.put("halide", "[F,Cl,Br,I]");
		queries.put("acyl chloride", "[Cl;X1][#6;A;X3;$([R0][#6]),$([H1R0])]=[O;X1]");
		queries.put("nitrile", "C#N");
		queries.put(
				"organophosphate",
				"[#8;A;X2;$([OX2][#6])][#15;A;X4D4]([$([OX2H]),$([OX1-]),$([OX2][#6])])([$([OX2H]),$([OX1-]),$([OX2][#6])])=[O;X1]");
		queries.put("arylphosphoester",
				"[#8;A;X1-,X2H1,X2C][P;X4]([#8;A;X1-,X2H1,X2C])(=[O;X1])c:c");
		queries.put("alkenyl",
				"[CX3;$([H2]),$([H1][#6]),$(C([#6])[#6])]=[CX3;$([H2]),$([H1][#6]),$(C([#6])[#6])]");
		queries.put("alkyne", "[CX2]#[CX2]");
		queries.put("allene", "[CX3]=[CX2]=[CX3]");
		queries.put("furan", "[c;R1]1[c;R1][c;R1][o;R1][c;R1]1");
		queries.put("dialkylether",
				"[OX2]([CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])])[CX4;!$(C([OX2])[O,S,#7,#15])]");
		queries.put("dialkylthioether",
				"[SX2]([CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])])[CX4;!$(C([OX2])[O,S,#7,#15])]");
		queries.put("alkylarylether", "[OX2](c)[CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])]");
		queries.put("diarylether", "[c][OX2][c]");
		queries.put("alkylarylthioether",
				"[SX2](c)[CX4;!$(C([OX2])[O,S,#7,#15,F,Cl,Br,I])]");
		queries.put("diarylthioether", "[c][SX2][c]");
		queries.put("aliphatic alcohol", "[#6;A;X4H1][#8;A;H1X2]");
		queries.put(
				"hydroxylamine",
				"[$([#8;v2;H1]-[#7;v3](-[#6])-[*,#1;C,c,#1]),$([#8;v1-]-[#7;v4+](-[*,#1;C,c,#1])=[#6]),$([#7;A;X3;$([H2]),$([H1][#6]),$([H0]([#6])[#6]);!$(NC=[O,N,S])][#8;A;X2;$([H1]),$(O[#6;!$(C=[N,O,S])])])]");
		// queries.put("hydroxylamine","[NX3H2,$([H1][#6]),$([H0]([#6])[#6]);!$(NC=[O,N,S])][OX2;$([H1]),$(O[#6;!$(C=[N,O,S])])]");
		queries.put("phenol",
				"[$([#8;A;H1X2][#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1),$([#8;A;H1X2][c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1)]");
		queries.put("primary carboxamide", "[CX3;$([R0][#6]),$([H1R0])](=[OX1])[NX3H2]");
		queries.put("secondary carboxamide",
				"[CX3;$([R0][#6]),$([H1R0])](=[OX1])[#7X3H1][#6;!$(C=[O,N,S])]");
		queries.put("tertiary carboxamide",
				"[CX3;$([R0][#6]),$([H1R0])](=[OX1])[#7X3H0]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])]");
		queries.put("aryl halide", "[F,Cl,Br,I][c]");
		queries.put("C ONS bond", "[#6]~[#7,#8,#16]");
		queries.put(
				"13-dicarbonyl",
				"[*,#1;*,#1;O,C,#1,N,X]-[#6;X3](=[O;X1])-[#6;H2X4]-[#6;X3](-[*,#1;*,#1;O,C,#1,N,X])=[O;X1]");
		queries.put("quaternary aliph ammonium",
				"[NX4H0+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]");
		queries.put("quaternary arom ammonium",
				"[NX4H0+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]");
		queries.put("o-quinone",
				"[O;X1]=[#6;R1]-1-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]-1=[O;X1]");
		queries.put("p-quinone",
				"[O;X1]=[#6;R1]-1-[#6;R1]=[#6;R1]-[#6;R1](=O)-[#6;R1]=[#6;R1]-1");
		queries.put("alkylguanidine",
				"[#6;X4][#7;A;v3X3,v4X4+][#6;X3R0]([#7;A;v3X3,v4X4+])=[#7;A;v3X2,v4X3+]"); // add
																							// this
																							// class
		queries.put("arylguanidine",
				"[#6;a][#7;A;v3X3,v4X4+][#6;X3R0]([#7;A;v3X3,v4X4+])=[#7;A;v3X2,v4X3+]"); // add
																							// this
																							// class
		queries.put("nitroarene",
				"[$([NX3](=O)=O),$([NX3+](=O)[O-])]-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1"); // add
																											// this
																											// class
		queries.put("nitrosoarene",
				"[O;X1]=[#7;X2]-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1"); // add
																						// this
																						// class
		queries.put("sulfide", "[#6]S[#6]");
		queries.put("sulfone",
				"[$([SX4](=[OX1])(=[OX1])([#6])[#6]),$([SX4+2]([OX1-])([OX1-])([#6])[#6])]");
		queries.put("sulfoxide",
				"[$([SX3](=[OX1])([#6])[#6]),$([SX3+]([OX1-])([#6])[#6])]");
		queries.put("thiodisulfinate", "[#6][S;X3]([#16;X2])=[O;X1]");
		queries.put(
				"thioamide",
				"[#7;X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])][$([CX3;!R][#6]),$([CX3H;!R])]=[S;X1]");
		queries.put(
				"amidoxime",
				"[#7;A;X3+0,X4+;!$([N]~[!#6]);!$([N]*~[#7,#8,#15,#16])][#6;X3]([#6,#1;A])=[#7;X2]-[#8;H1X2]");
		queries.put("carboxamidine",
				"[#7;A;X3;!$(NC=[O,S])][#6;A;X3;$([CH]),$([C][#6])]=[#7;A;X2;!$(NC=[O,S])]");
		queries.put(
				"n-hydroxyguanidine",
				"[#7;A;v3X3,v4X4+][#6;X3](=[#7;A;v3X2,v4X3+]/[#8;H1X2])[#7;A;v3X3,v4X4+]([#6,#1;A])[#6,#1;A]");
		queries.put("arylhydrazine", "[#6;a]-[#7H1]-[#7H1]-[#6;a]");
		queries.put("thiol", "[#6]-[#16;H1X2]");
		queries.put("alkylthiol", "[SX2H][CX4;!$(C([SX2H])~[O,S,#7,#15])]");
		queries.put("acyloxymethyl ester",
				"[C;X4H2]([#8;X2]-[#6;X3]([#6,#1;A])=[O;X1])[#8;X2]-[#6;X3]([#6,#1;A])=[O;X1]");
		queries.put("1-(acyloxy)ethyl ester",
				"[C;X4H1]([#6;H3X4])([#8;X2]-[#6;X3]([#6,#1;A])=[O;X1])[#8;X2]-[#6;X3]([#6,#1;A])=[O;X1]");
		queries.put("alkyl carbamate",
				"[#6]-[#8;X2]-[#6;X3](=[O;X1])-[#7;X3](-[#6;X4])[#6,#1;A]");
		queries.put(
				"(acyloxy)alkyl carbamate",
				"[C;X4H1]([#6,#1;A])([#8]-[#6;X3]([#6,#1;A])=[O;X1])[#8;X2]-[#6;X3](=[O;X1])-[#7;X3](-[#6;X4])[#6,#1;A]");
		queries.put("epoxide", "[OX2r3]1[#6r3][#6r3]1");
		queries.put("aryl ester", "[#6;a]-[#8;X2]-[#6;X3]([#6,#1;A])=[O;X1]");
		queries.put("purine", "[c;R1]1[n;R1]c2[c;R1][n;R1][c;R1][n;R1]c2[n;R1]1");
		queries.put("pyrimidine", "[c;R1]1[c;R1][n;R1][c;R1][n;R1][c;R1]1");
		queries.put("thiocyanate", "[NX2]=[CX2]=[SX1]");
		queries.put("isothiocyanate", "[SX2][CX2]#[NX1]");
		queries.put("alpha,beta-unsaturated system",
				"[CX3]=[CX3][$([CX3]=[O,N,S]),$(C#[N]),$([S,P]=[OX1]),$([NX3]=O),$([NX3+](=O)[O-])]"); // Michael's
																										// acceptor
		queries.put("ketene", "[CX3]=[CX2]=[OX1]");
		queries.put("allylic alcohol", "[#8;H0X2]-[#6]([#6,#1;A])-[#6]=[#6]");
		queries.put(
				"CH-acidic",
				"[$([CX4;!$([H0]);!$(C[!#6;!$([P,S]=O);!$(N(~O)~O)])][$([CX3]=[O,N,S]),$(C#[N]),$([S,P]=[OX1]),$([NX3]=O),$([NX3+](=O)[O-]);!$(*[S,O,N;H1,H2]);!$([*+0][S,O;X1-])]),$([CX4;!$([H0])]1[CX3]=[CX3][CX3]=[CX3]1)]");
		queries.put(
				"CH-acidic strong",
				"[CX4;!$([H0]);!$(C[!#6;!$([P,S]=O);!$(N(~O)~O)])]([$([CX3]=[O,N,S]),$(C#[N]),$([S,P]=[OX1]),$([NX3]=O),$([NX3+](=O)[O-]);!$(*[S,O,N;H1,H2]);!$([*+0][S,O;X1-])])[$([CX3]=[O,N,S]),$(C#[N]),$([S,P]=[OX1]),$([NX3]=O),$([NX3+](=O)[O-]);!$(*[S,O,N;H1,H2]);!$([*+0][S,O;X1-])]");
		queries.put("N-aryl_Np-hydroxyguanidine",
				"[#6;a][#7;A;v3X3,v4X4+]([#6,#1;A])[#6;X3]([#7;A;v3X3,v4X4+])=[#7;A;v3X2,v4X3+][#8;H1X2]");
		queries.put("primary aliph amine",
				"[NX3H2+0,NX4H3+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]");
		queries.put("secondary aliph amine",
				"[NX3H1+0,NX4H2+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]");
		queries.put("tertiary aliph amine",
				"[NX3H0+0,NX4H1+;!$([N][!C]);!$([N]*~[#7,#8,#15,#16])]");
		queries.put("primary arom amine", "[NX3H2+0,NX4H3+]c");
		queries.put("secondary arom amine",
				"[NX3H1+0,NX4H2+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]");
		queries.put("tertiary arom amine",
				"[NX3H0+0,NX4H1+;!$([N][!c]);!$([N]*~[#7,#8,#15,#16])]");
		queries.put("peroxo group", "[OX2D2][OX2D2]");
		queries.put("oximether",
				"[NX2](=[CX3;$([CH2]),$([CH][#6]),$([C]([#6])[#6])])[OX2][#6;!$(C=[#7,#8])]");
		queries.put(
				"thioamide S-oxide derivative",
				"[#7;X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])][$([CX3;!R][#6]),$([CX3H;!R])]=[S;X2+][#8;X1-]");
		queries.put(
				"thioamide SS-dioxide derivative",
				"[#8;X1-][S;X3](=[O;X1])[$([CX3;!R][#6]),$([CX3H;!R])]=[#7;NX2H1,$([#7][#6;!$(C=[O,N,S])])]");
		queries.put("arylsulfate", "[#6;a][S;X4]([#8;A;X1-,X2H1])(=[O;X1])=[O;X1]");
		queries.put("primary alcohol", "[OX2H][CX4H2;!$(C([OX2H])[O,S,#7,#15])]");
		queries.put("secondary alcohol", "[OX2H][CX4H1;!$(C([OX2H])[O,S,#7,#15])]");
		queries.put("tertiary alcohol", "[OX2H][CX4D4;!$(C([OX2H])[O,S,#7,#15])]");
		queries.put("n-nitroso",
				"[#7;A;X3+0,X4+;!$([N]~[!#6]);!$([N]*~[#7,#8,#15,#16])][#7;X2]=[O;X1]");
		queries.put("c-nitroso", "[#6]-[#7;X2]=[O;X1]");
		queries.put("dithiocarbamic acid ester",
				"[#6]-[#16;X2]-[#6;X3](=[S;X1])-[#7;X3]([#6,#1;A])[#6,#1;A]");
		queries.put("dithiocarbamic acid",
				"[#16;X2H1]-[#6;X3](=[S;X1])-[#7;X3]([#6,#1;A])[#6,#1;A]");
		queries.put(
				"sulfonamide",
				"[#7;X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])][#16;A;X4;$([H1]),$([H0][#6])](=[O;X1])=[O;X1]");
		queries.put(
				"steroid",
				"[#6;R1]-,=1-,=[#6;R1]-,=[#6]-,=2-,=[#6;R1]-,=[#6;R1]-,=[#6]-3-,=[#6]-,=4-,=[#6;R1]-,=[#6;R1]-,=[#6;R1]-,=[#6;R1]-,=[#6]-,=4-,=[#6;R1]-,=[#6;R1]-,=[#6]-3-,=[#6]-,=2-,=[#6;R1]-,=1");
		queries.put("imidazolidine", "N1CCNC1");
		queries.put("alkylarylsulfoxide",
				"[$([SX3](=[OX1])([#6])[#6;a]),$([SX3+]([OX1-])([#6])[#6;a])]");
		queries.put("alkylarylsulfone", "[$([SX4](=[OX1])(=[OX1])([#6])[#6;a])]");
		queries.put("thiophene", "C1=CSC=C1");
		queries.put("thiophene s-oxide", "C1=C[S+]([OX1-])C=C1");
		queries.put("1,3-thiazole", "S1C=CN=C1"); // SMARTCyp - a 2D-method for
													// Prediction of Cytochrome
													// P450 Mediated Drug
													// Metabolism_Supplement
		queries.put("1,2-thiazole", "S1C=CC=N1");
		queries.put("quinoline", "[$(C1=CC=C2N=CC=CC2=C1),$(c1ccc2ncccc2c1)]");
		queries.put("pyridazine", "C1=CC=NN=C1");
		queries.put("1,3-dioxolane", "C1OCOC1");
		queries.put("xanthine", "[$(Oc1nc(O)c2[nH]cnc2n1),$(O=C1NC2=C(NC=N2)C(=O)N1)]");
		queries.put(
				"biphenyl",
				"[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1](-[#6;R1]=[#6;R1]-1)-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1");
		queries.put(
				"1,2-aminoalcohol",
				"[#8;X2H1][C;X4]([*,#1;#1,CX4,$([cX3](:*):*),$([cX2+](:*):*)])([*,#1;#1,CX4,$([cX3](:*):*),$([cX2+](:*):*)])[C;X4]([*,#1;#1,CX4,$([cX3](:*):*),$([cX2+](:*):*)])([*,#1;#1,CX4,$([cX3](:*):*),$([cX2+](:*):*)])[#7;X3]([H])-[*,#1;#1,CX4,$([cX3](:*):*),$([cX2+](:*):*)]");
		queries.put("organosulfur compound", "[#6]~[#16]");
		queries.put("secondary aliphatic/aromatic amine",
				"[#7;v3X3H1]([#6;A;X4])-[$([cX3](:*):*),$([cX2+](:*):*)]");
		queries.put(
				"propargyl-type 13-dipolar organic compound",
				"[$([#6,#7,#8;A;-][N+]#C),$([#6-]=[N+]=[#6,#7,#8;A]),$([#6,#7,#8;A;+][#7]=[#6-]),$([#6]-[#7]=[#6,#7,#8;A])].[#6]");
		queries.put(
				"diphenylmethane",
				"[$([*,#1;!$(c1ccccc1)]C([*,#1;!$(c1ccccc1)])(!@-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1)!@-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1),$(*=[#6](!@-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1)!@-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1)]");
		queries.put(
				"phenol ether",
				"[#6;A;X4;!$(C([SX2])[O,S,#7,#15])][#8;X2]!@-[#6;R1]-1=[#6;R1]-[#6;R1]=[#6;R1]-[#6;R1]=[#6;R1]-1");
		queries.put("heteroaromatic", "[a;!c]"); // 107 bits
		queries.put("nitro", "[#6]-[#7;X3+](-[#8;X1-])=[O;X1]");
		queries.put("azo", "[#6]-[#7;X2]=[#7;X2]-[#6]");
		queries.put(
				"hydroxamid_acid",
				"[CX3;$([H0][#6]),$([H1])](=[OX1])[#7X3;$([H1]),$([H0][#6;!$(C=[O,N,S])])][$([OX2H]),$([OX1-])]");
		queries.put(
				"hydroxamid_acid_ester",
				"[CX3;$([H0][#6]),$([H1])](=[OX1])[#7X3;$([H1]),$([H0][#6;!$(C=[O,N,S])])][OX2][#6;!$(C=[O,N,S])]");
		queries.put("pyridine", "C1=CC=NC=C1");
		queries.put("thiophenol", "[#16;H1X2]-[c;R1]1[c;R1][c;R1][c;R1][c;R1][c;R1]1");
		queries.put("organohalide", "[#6][F,Cl,Br,I]");
		queries.put("hydrazine_derivative", "[#6][NX3H1][NX3H2]");
		queries.put("carboxylic ester",
				"[CX3;$([R0][#6]),$([H1R0])](=[OX1])[OX2][#6;!$(C=[O,N,S])]");
		queries.put(
				"steroid",
				"[#6,#8,#7,#16]-,=1-,=[#6,#8,#7,#16]-,=[#6,#8,#7,#16;A;R2]-,=2-,=[#6,#8,#7,#16]-,=[#6,#8,#7,#16;A;R1]-,=[#6,#8,#7,#16]~3~[#6,#8,#7,#16;A;R2]~4~[#6,#8,#7,#16;A;R1]~[#6,#8,#7,#16]~[#6,#8,#7,#16]~[#6,#8,#7,#16]~[#6,#8,#7,#16;A]~4~[#6,#8,#7,#16]~[#6,#8,#7,#16;A;R1]~[#6,#8,#7,#16;A;R2]~3-,=[#6,#8,#7,#16]-,=2-,=[#6,#8,#7,#16]-,=1");
		queries.put("phosphate monoester",
				"[#6]-[#8;X2]P([#8;A;X2H1,X1-])([#8;A;X2H1,X1-])=[O;X1]");
		queries.put("3-acylpyruvate",
				"[#8-]-[#6](=[O;X1])-[#6;X3](=[O;X1])-[#6;H2X4]-[#6;X3]([#6,#1;A])=[O;X1]");
		queries.put("n-acyl-aromatic alpha-amino-acid",
				"[#6;a]-[#6][#6;A;H1X4;@]([#7;A;H1X3][#6]([#6,#1;A])=O)[#6]([#8;A;X2H1,X1-])=O");
		queries.put(
				"n-acyl-aliphatic alpha-amino-acid",
				"[#6;A;X4;CX4H3,$([CX4]-[C;A])][#6;A;H1X4;@]([#7;A;H1X3][#6]([#6,#1;A])=O)[#6]([#8;A;X2H1,X1-])=O");
		queries.put("nitrosodialkylamine",
				"[#6;X4][#7;A;X3+0,X4+;!$([N]~[!#6]);!$([N]*!@[#7,#8,#15,#16])]([#6;X4])[#7]=[O;X1]");
		queries.put("vinyl alcohol", "[#8;X2H1]-[#6]=[#6]");
		queries.put("nitroimine", "[#8;X1-]-[#7;X3+](=[O;X1])-[#7]=[#6;X3]");
		queries.put("1,2-oxazole", "[#8]-1-[#6]=[#6]-[#6]=[#7]-1");
		queries.put("1,3-oxazole", "[#8]-1-[#6]=[#6]-[#7]=[#6]-1");
		queries.put("glycerol", "[#8][#6;A;H2X4][#6;A;H1X4]([#8])[#6;A;H2X4][#8]");
		queries.put("peptide",
				"[#7]-[#6](-[$([#1,*])])-[#6](=O)-[#7]-[#6](-[$([#1,*])])-[#6]([#8,#7;A])=O");
		queries.put(
				"coenzyme a",
				"[#6;R0]-,=[#6;R0]-,=[#6;R0][#6;A;X4R0;H1,H2][#6;H2X4R0]-[#6;R0](=[O;R0])-[#16;R0]-[#6;R0]-[#6;R0]-[#7;R0]-,=[#6;R0](-,=[#8])-[#6;R0]-[#6;R0]-[#7;R0]-,=[#6;R0](-,=[#8;R0])[#6;A;H1X4]([#8])[C;R0]([#6;R0])([#6;R0])[#6;R0]-[#8;R0][P;R0]([!#1!#6;O,$([O-])])(=[O;R0])[#8;R0][P;R0]([!#1!#6;O,$([O-])])(=[O;R0])[#8]-[#6]-[#6]-1-[#8]-[#6](-[#6](-[#8])-[#6]-1-[#8]P([!#1!#6;O,$([O-])])([!#1!#6;O,$([O-])])=O)-n1cnc2c(-[#7])ncnc12");
		queries.put("fatty acyl",
				"[#6]!@[#6]!@[#6]-[#6](-[OX2H1,$(O-[#6]),OX1-,NX3,SX2])=O");
		queries.put(
				"2-N-linked ribose deriv.",
				"[$([#7;R1]-[#6]-1-[#8]-[#6](-[#6]-[#8])-[#6]=[#6]-1),$([#7;R1]-[#6]-1=[#6]-[#6]-[#6](-[#6]-[#8])-[#8]-1),$([#7;R1]-[#6]-1-,=[#6]-[#6]-,=[#6](-[#6]-[#8])-[#8]-1)]");
		queries.put("aromatic aocohol", "[#6;a][#6;A;X4;H2][#8;H1X2]");
		queries.put("glycerol-3-phosphate",
				"[#8][#6;A;H2X4][#6;A;H1X4]([#8])[#6;A;H2X4][#8]P([#8])([#8])=O");

		return queries;
	}
	
	/**
	 * 
	 * @param mole
	 *            : A molecule of interest
	 * @param queries
	 *            : A HashMap with the functional groups and patterns to detect,
	 *            with their SMARTS patterns
	 * @return : A list containing atom-based list of the bits (1 or 0) that
	 *         describe for every functional group/pattern, whether the atom of
	 *         interest is part of a match.
	 * @throws Exception
	 */

	public ArrayList<ArrayList<Integer>> generateClassyfireAtomFingeprint(
			IAtomContainer mole, LinkedHashMap<String, String> queries) throws Exception {
		ArrayList<ArrayList<Integer>> results = new ArrayList<ArrayList<Integer>>();
		ArrayList<ArrayList<Integer>> allBits = new ArrayList<ArrayList<Integer>>();
		//LinkedHashMap<Integer, int[]> atomFingerprints = new LinkedHashMap<Integer, int[]>();

		for (Map.Entry<String, String> item : queries.entrySet()) {
			// System.out.println(item.getKey());
			IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();
			SMARTSQueryTool smartsPattern = new SMARTSQueryTool(item.getValue(), builder);
			List<List<Integer>> matchingAtoms = getMatchingAtoms(mole, smartsPattern);
			ArrayList<Integer> bits = new ArrayList<Integer>();

			if (matchingAtoms.size() > 0) {
//				System.out.println(item.getKey());
//				System.out.println("matching_atoms: " + matchingAtoms);
				for (int i = 0; i < matchingAtoms.size(); i++) {
					// System.out.println(matching_atoms.get(i));
					List<Integer> indexes = matchingAtoms.get(i);
					for (int j = 0; j < indexes.size(); j++) {
						bits.add(indexes.get(j));
						// results.add(bits);
					}
				}
			}
//			if (bits.size() > 0) {
//				System.out.println("Bits: " + bits);
//			}

			allBits.add(new ArrayList<Integer>(new HashSet<Integer>((bits)))); // to
																				// remove
																				// duplicates
		}

//		System.out.println(allBits);
//		System.out.println(allBits.size());

		for (int k = 0; k < mole.getAtomCount(); k++) {
			ArrayList<Integer> nilArraylist = new ArrayList<Integer>();
			for (int a = 0; a < queries.size(); a++) {
				nilArraylist.add(0);
			}
			results.add(nilArraylist);
		}

		for (int l = 0; l < allBits.size(); l++) {
			if (allBits.get(l).size() > 0) {
				for (int m = 0; m < allBits.get(l).size(); m++) {
					int index = allBits.get(l).get(m);
					results.get(index).set(l, 1);
				}
			}
		}

		return results;
	}

	/**
	 * 
	 * @param sdfInput
	 *            : a SDF file
	 * @param queries
	 *            : a HashMap with fingerprint patterns and their SMARTS
	 *            expressions
	 * @return An array list of fingerprints for every atom of each molecule in
	 *         the SDF file
	 * @throws Exception
	 */

	public ArrayList<ArrayList<ArrayList<Integer>>> generateSerialAtomfingerprintToArraylist(
			File sdfInput, LinkedHashMap<String, String> queries) throws Exception {
//		String absolutePath = sdfInput.getAbsolutePath();
//		String destinationDirectory = absolutePath.substring(0,
//				absolutePath.lastIndexOf(File.separator));
//
//		String[] substrings = sdfInput.getName().split("/");
//		String name = substrings[substrings.length - 1].split(".sdf")[0];
//		File sdf_output = new File(destinationDirectory + "/atom_fingerprint_" + name
//				+ ".tsv");
//		FileOutputStream sdfOutputOs = new FileOutputStream(sdf_output);
//		BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(sdfOutputOs));
//
//		StringBuilder patternNames = new StringBuilder();
//		for (Map.Entry<String, String> item : queries.entrySet()) {
//			patternNames.append(item.getKey()).append("\t");
//		}
//		System.out.println(patternNames);
//		bw.write("molecule ID\tatom ID\t" + patternNames + "\n"); // add header

		IteratingSDFReader reader = new IteratingSDFReader(new FileInputStream(sdfInput),
				DefaultChemObjectBuilder.getInstance());
		ArrayList<ArrayList<ArrayList<Integer>>> totalBits = new ArrayList<ArrayList<ArrayList<Integer>>>();
		ElectronDonation model = ElectronDonation.cdk();
		CycleFinder cycles = Cycles.cdkAromaticSet();
		Aromaticity aromaticity = new Aromaticity(model, cycles);

		int i = 0;
		while (reader.hasNext()) {
			i = i + 1;
			IAtomContainer molecule = (IAtomContainer) reader.next();
			AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule); // 1
			aromaticity.apply(molecule); // 1
			// Kekulization.kekulize(molecule); // 2

			 CDKHydrogenAdder adder =
			 CDKHydrogenAdder.getInstance(molecule.getBuilder()); // 3
			 adder.addImplicitHydrogens(molecule); // 3
			// CDKHydrogenAdder.getInstance(builder).addImplicitHydrogens(molecule);
			// //http://efficientbits.blogspot.ca/2013/12/new-smiles-behaviour-parsing-cdk-154.html
			
//			System.out.println("\n\n\n\n----------------------> " + molecule.getProperty(CDKConstants.TITLE));
			ArrayList<ArrayList<Integer>> atom_fingerprint = generateClassyfireAtomFingeprint(
					molecule, queries);
			totalBits.add(atom_fingerprint);
		}

		for (int j = 0; j < totalBits.size(); j++) {
//			System.out.println("Number of atoms of molecule " + (j + 1) + ": "
//					+ totalBits.get(j).size());
			for (int k = 0; k < totalBits.get(j).size(); k++) {
				ArrayList<Integer> bits = totalBits.get(j).get(k);
				StringBuilder fingerprintBits = new StringBuilder();

				for (int bit : bits) {
					fingerprintBits.append(bit).append("\t"); // separating
																// contents
																// using commas
				}
//				bw.write((j + 1)
//						+ "\t"
//						+ (k + 1)
//						+ "\t"
//						+ fingerprintBits.deleteCharAt((fingerprintBits.length() - 1))
//								.toString() + "\n");
//				System.out.println((j + 1) + "\t" + (k + 1) + "\t" + fingerprintBits);
				// sdf_output.

			}
		}
//		bw.close();
		return totalBits;

	}



	/**
	 * 
	 * @param sdfInput
	 *            : a SDF file
	 * @param queries
	 *            : a HashMap with fingerprint patterns and their SMARTS
	 *            expressions
	 * @return : A string with the atom-based fingerprints of every molecule. A
	 *         .csv file is also saved.
	 * @throws Exception
	 */

	public String saveSerialAtomFingerprinterToCSV(File sdfInput,
			LinkedHashMap<String, String> queries) throws Exception {
		String absolutePath = sdfInput.getAbsolutePath();
		String destinationDirectory = absolutePath.substring(0,
				absolutePath.lastIndexOf(File.separator));
		String[] substrings = sdfInput.getName().split("/");
		String name = substrings[substrings.length - 1].split(".sdf")[0];
		String output = destinationDirectory + "/" + name + "_atom_fingerprint" + ".csv";
		File sdfOutput = new File(output);
		FileOutputStream sdfOutputOs = new FileOutputStream(sdfOutput);
		BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(sdfOutputOs));

		StringBuilder patternNames = new StringBuilder();
		for (Map.Entry<String, String> item : queries.entrySet()) {
			patternNames.append(item.getKey()).append(",");
		}
		// System.out.println(pattern_names);
		// bw.write("molecule ID\tatom ID\t" + pattern_names + "\n"); //add
		// header

		IteratingSDFReader reader = new IteratingSDFReader(new FileInputStream(sdfInput),
				DefaultChemObjectBuilder.getInstance());
		ArrayList<ArrayList<ArrayList<Integer>>> totalBits = new ArrayList<ArrayList<ArrayList<Integer>>>();
		ElectronDonation model = ElectronDonation.cdk();
		CycleFinder cycles = Cycles.cdkAromaticSet();
		Aromaticity aromaticity = new Aromaticity(model, cycles);

		int i = 0;
		while (reader.hasNext()) {
			i = i + 1;
			IAtomContainer molecule = (IAtomContainer) reader.next();
			AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule); // 1
			aromaticity.apply(molecule); // 1
			ArrayList<ArrayList<Integer>> atomFingerprint = generateClassyfireAtomFingeprint(
					molecule, queries);
			totalBits.add(atomFingerprint);
		}

		for (int j = 0; j < totalBits.size(); j++) {
			for (int k = 0; k < totalBits.get(j).size(); k++) {
				ArrayList<Integer> bits = totalBits.get(j).get(k);
				StringBuilder fingerprintBits = new StringBuilder();

				for (int bit : bits) {
					fingerprintBits.append(bit).append(","); // separating
																// contents
																// using commas
				}
				bw.write((j + 1)
						+ ","
						+ (k + 1)
						+ ","
						+ fingerprintBits.deleteCharAt((fingerprintBits.length() - 1))
								.toString() + "\n");
			}
		}
		bw.close();
		return output;
	}

}