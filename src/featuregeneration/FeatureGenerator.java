/**
 * 
 */
package featuregeneration;

/**
 * This class implements several functions to generate the different features of molecules,
 * including separate functions for atomic and molecular features generation, 
 * and a function used to read and preprocess the molecules from a chemical file.
 * @author Zheng Shi
 *
 */
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Properties;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.atomtype.SybylAtomTypeMatcher;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.NoSuchAtomTypeException;
import org.openscience.cdk.geometry.surface.NumericalSurface;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.graph.PathTools;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.io.listener.PropertiesListener;
import org.openscience.cdk.qsar.DescriptorValue;
import org.openscience.cdk.qsar.descriptors.atomic.AtomDegreeDescriptor;
import org.openscience.cdk.qsar.descriptors.atomic.AtomHybridizationDescriptor;
import org.openscience.cdk.qsar.descriptors.atomic.AtomValenceDescriptor;
import org.openscience.cdk.qsar.descriptors.atomic.EffectiveAtomPolarizabilityDescriptor;
import org.openscience.cdk.qsar.descriptors.atomic.PartialSigmaChargeDescriptor;
import org.openscience.cdk.qsar.descriptors.atomic.PartialTChargeMMFF94Descriptor;
import org.openscience.cdk.qsar.descriptors.atomic.PiElectronegativityDescriptor;
import org.openscience.cdk.qsar.descriptors.atomic.SigmaElectronegativityDescriptor;
import org.openscience.cdk.qsar.descriptors.atomic.StabilizationPlusChargeDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.ALOGPDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.APolDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.HBondAcceptorCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.HBondDonorCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.MomentOfInertiaDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.RotatableBondsCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.TPSADescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.WeightDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.XLogPDescriptor;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import utils.ChemSearcher;
import weka.core.Attribute;
import weka.core.FastVector;
import weka.core.Instance;
import weka.core.Instances;


public class FeatureGenerator {
	/**
	 * This functions reads molecules from a chemical file and preprocess the molecules.
	 * 
	 * @param pathToInputFile
	 *            : path of the file which contains the input molecules
	 * @return an IAtomContainerSet contains all the input molecules
	 * @throws FileNotFoundException
	 * @throws CDKException 
	 */
	public IAtomContainerSet readFile(String pathToInputFile)
			throws FileNotFoundException, CDKException {
		IChemObjectBuilder bldr = SilentChemObjectBuilder.getInstance();
		IteratingSDFReader sdfr = new IteratingSDFReader(new FileReader(pathToInputFile),
				bldr);
		Properties prop = new Properties();
		prop.setProperty("ForceReadAs3DCoordinates", "true");
		PropertiesListener listener = new PropertiesListener(prop);
		sdfr.addChemObjectIOListener(listener);
		sdfr.customizeJob();
		IAtomContainerSet MOLS = DefaultChemObjectBuilder.getInstance().newInstance(
				IAtomContainerSet.class);

		while (sdfr.hasNext())
				MOLS.addAtomContainer( preprocessContainer(sdfr.next()));
		return MOLS;

	}
	/**
	 * This function applies some preprocessing operations, such as setting the
	 * flag of atoms from aromatic rings to "ISAROMATIC", and kelulizing
	 * molecules.
	 * 
	 * @param molecule
	 *            : A molecule of interest
	 * @return : A processed molecule (AtomContainer)
	 * @throws CDKException
	 */
	public IAtomContainer preprocessContainer(IAtomContainer molecule)
			throws CDKException {
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
		CDKHydrogenAdder.getInstance(molecule.getBuilder())
				.addImplicitHydrogens(molecule);
	    Aromaticity aromaticity = new Aromaticity( ElectronDonation.cdk(), Cycles.cdkAromaticSet());
		for (IBond bond : molecule.bonds()) {
			if (bond.getFlag(CDKConstants.SINGLE_OR_DOUBLE)) {
				bond.setFlag(CDKConstants.ISAROMATIC, true);
				bond.getAtom(0).setFlag(CDKConstants.ISAROMATIC, true);
				bond.getAtom(1).setFlag(CDKConstants.ISAROMATIC, true);
			} 
		}
		aromaticity.apply(molecule);
		return molecule;
	}
	/**
	 * This function is used to generate the molecular features for the
	 * molecules contained in the input file.
	 * 
	 * @param pathToInputFile
	 *            : path of the file which contains the input molecules
	 * @return an arraylist which contains the molecular features.
	 * @throws IOException
	 * @throws CDKException
	 */
	public ArrayList<String> generateMolecularFeatures(IAtomContainerSet set, String pathToInputFile)
			throws IOException, CDKException {
		StringBuffer sb = new StringBuffer();
		int i;
		for (i = 0; i < 11; i++) {
			sb.append(i + 1).append(",");
		}
		sb.deleteCharAt(sb.length() - 1);
		ArrayList<String> result = new ArrayList<String>();
		int l=set.getAtomContainerCount();
		for(i=0;i<l;i++) {
			IAtomContainer container = set.getAtomContainer(i);
			AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(container);
			sb.setLength(0);
			try {
				String s = generateMolecularFeatures(container);
				result.add(s);
			}
			catch (CDKException e) {
				System.err.println("Compound nr." + (i+1) + ": Title = " + container.getProperty(CDKConstants.TITLE) + "\n"+ "CDKException.");
			}
			catch (java.lang.NullPointerException je) {
				System.err.println("Compound nr." + (i+1) + ": Title = " + container.getProperty(CDKConstants.TITLE) + "\n"+ "java.lang.NullPointerException.");
			}
		}
		return result;
	}

	/**
	 * This function generates molecular features for a single molecule
	 * 
	 * @param molecule
	 * @return String representation of molecular features where each feature is separated by a comma
	 * @throws CDKException
	 */
	public String generateMolecularFeatures(IAtomContainer molecule) throws CDKException {
		IAtomContainer mol = molecule;
		String molecular;
		String []res=new String[9];
		ALOGPDescriptor Alogp = new ALOGPDescriptor();
		res[0] =  Alogp.calculate(mol).getValue().toString().split(",")[0];
		APolDescriptor Apol = new APolDescriptor();
		res[1] = Apol.calculate(mol).getValue().toString();
		HBondAcceptorCountDescriptor Hb = new HBondAcceptorCountDescriptor();
		res[2] = Hb.calculate(mol).getValue().toString();
		HBondDonorCountDescriptor HB = new HBondDonorCountDescriptor();
		res[3] = HB.calculate(mol).getValue().toString();
		MomentOfInertiaDescriptor Mo = new MomentOfInertiaDescriptor();
		res[4] = Mo.calculate(mol).getValue().toString().split(",")[6];
		RotatableBondsCountDescriptor Ro = new RotatableBondsCountDescriptor();
		res[5] = Ro.calculate(mol).getValue().toString();
		TPSADescriptor Tp = new TPSADescriptor();
		res[6] = Tp.calculate(mol).getValue().toString();
		WeightDescriptor We = new WeightDescriptor();
		res[7]= We.calculate(mol).getValue().toString();
		XLogPDescriptor Xl = new XLogPDescriptor();
		res[8]= Xl.calculate(mol).getValue().toString();
		StringBuffer sb=new StringBuffer();
		for(int i=0;i<res.length;i++)
			sb.append(res[i]).append(",");
		NumericalSurface ns = new NumericalSurface(mol);
		
		try{
			ns.calculateSurface();
		}catch (java.lang.IllegalArgumentException je){
			System.out.println(je.getMessage());
		}
	
		double ASA = ns.getTotalSurfaceArea();
		molecular = sb.toString() + String.valueOf(ASA);
		return molecular;
	}
	
	
	/**
	 * This function is used to generate the atomic features for the molecules
	 * contained in the input file.
	 * 
	 * @param input_path
	 *            : path of the file which contains the input molecules
	 * @return String representation of atomic features where each feature is separated by a comma
	 * @throws CDKException
	 * @throws IOException
	 */
	public ArrayList<ArrayList<String>> generateAtomicFeatures(IAtomContainerSet set)
			throws CDKException, IOException {
		ArrayList<ArrayList<String>> result = new ArrayList<ArrayList<String>>();
		int i,j,length;
		int l=set.getAtomContainerCount();
		for(i=0;i<l;i++) {
			IAtomContainer container = set.getAtomContainer(i);
			CDKHydrogenAdder adder = CDKHydrogenAdder.getInstance(container.getBuilder());
			AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(container);
			adder.addImplicitHydrogens(container);
			length = container.getAtomCount();
			ArrayList<String> list = new ArrayList<String>();
			try {
				String[] atomicFeature = generateAtomicFeatures(container);
				String[] atomicType = generateAtomType(container);
				String[] path1Feature = generatePath(container, 1);
				String[] path2Feature = generatePath(container, 2);

				for ( j = 0; j < length; j++) {
					StringBuffer sb=new StringBuffer();
					sb.append(atomicFeature[j]).append(",");
					sb.append(atomicType[j]).append(",");
					sb.append(path1Feature[j] ).append(",");
					sb.append(path2Feature[j]);
					list.add(sb.toString());
				}
				result.add(list);

			}
			catch (CDKException e) {
				System.err.println("Compound nr." + (i+1) + ": Title = " + container.getProperty(CDKConstants.TITLE) + "\n"+ "CDKException: " + e.getMessage());
			}
			
		}
		return result;
	}

	/**
	 * definition for all the adopted SYBYL atom types
	 */
	public static final String[]	str	= { "C.1", "C.2", "C.3", "C.ar", "C.cat", "N.1",
			"N.2", "N.3", "N.4", "N.ar", "N.am", "N.pl3", "O.2", "O.3", "O.co2", "S.2",
			"S.3", "S.O", "S.O2", "P.3", "F", "Cl", "Br", "I" };

	/**
	 * This function is used to generate the atomic type for a single molecule.
	 * 
	 * @param molecule
	 * @return String representation of atomic type of one molecule 
	 * where each feature is separated by a comma
	 * @throws CDKException
	 */
	public String[] generateAtomType(IAtomContainer molecule) throws CDKException {
		IChemObjectBuilder bldr = SilentChemObjectBuilder.getInstance();
		SybylAtomTypeMatcher satm = SybylAtomTypeMatcher.getInstance(bldr);
		IAtomContainer mol = molecule;
		Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
		aromaticity.apply(mol);
		int length = mol.getAtomCount();
		int i;
		int[] symbol = new int[24];
		String[] molecular = new String[length];
		for (i = 0; i < length; i++) {
			int position = 0;
			IAtom a = mol.getAtom(i);
			IAtomType Atom = satm.findMatchingAtomType(mol, a);
			String type = Atom.getAtomTypeName();
			for (int k = 0; k < 24; k++)
				symbol[k] = 0;
			for (int k = 0; k < 24; k++) {
				if (type.equalsIgnoreCase(str[k]))
					position = k;
			}
			symbol[position] = 1;
			StringBuffer strbuf = new StringBuffer();
			for (int n = 0; n < 24; n++) {
				strbuf.append(symbol[n]).append(',');
			}
			strbuf.deleteCharAt(strbuf.length() - 1);
			molecular[i] = strbuf.toString();
		}
		return molecular;

	}

	/**
	 * This function is used to generate the environmental features for a single
	 * molecule.
	 * 
	 * @param molecule
	 * @param depth
	 *            : the bond length of an atom to the centered atom.
	 * @return String representation of the environmental features for a single
	 *         molecule where each feature is separated by a comma
	 * @throws CDKException
	 */
	public String[] generatePath(IAtomContainer molecule, int depth) throws CDKException {
		IChemObjectBuilder bldr = SilentChemObjectBuilder.getInstance();
		IAtomContainer mol = molecule;
		Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
		aromaticity.apply(mol);
		int dep = depth;
		SybylAtomTypeMatcher satm = SybylAtomTypeMatcher.getInstance(bldr);
		int length = mol.getAtomCount();
		int i = 0;
		int o;
		int[] symbol = new int[24];
		String[] molecular = new String[length];
		for (i = 0; i < length; i++) {
			IAtom a = mol.getAtom(i);
			for (int k = 0; k < 24; k++)
				symbol[k] = 0;
			int bond_1 = 0;
			int bond_2 = 0;
			int bond_3 = 0;
			int bond_4 = 0;
			List<List<IAtom>> b = PathTools.getPathsOfLength(mol, a, dep);
			List<IAtom> c;
			IAtom d, e;
			int m = b.size();
			for (int j = 0; j < m; j++) {
				int position = 0;
				c = b.get(j);
				d = c.get(c.size() - 1);
				e = c.get(c.size() - 2);
				IAtomType s_1 = satm.findMatchingAtomType(mol, d);
				String type_1 = s_1.getAtomTypeName();
				for (int k = 0; k < 24; k++) {
					if (type_1.equalsIgnoreCase(str[k]))
						position = k;
				}
				symbol[position] = 1;
				IBond bond = mol.getBond(d, e);
				IBond.Order order = bond.getOrder();
				o = order.numeric();
				boolean flag = bond.getFlag(CDKConstants.ISAROMATIC);
				if (flag)
					bond_4 = 1;
				else if (o == 1)
					bond_1 = 1;
				else if (o == 2)
					bond_2 = 1;
				else
					bond_3 = 1;
			}
			StringBuffer strbuf = new StringBuffer();
			for (int n = 0; n < 24; n++) {
				strbuf.append(symbol[n]).append(',');
			}
			strbuf.append(bond_1).append(",");
			strbuf.append(bond_2).append(",");
			strbuf.append(bond_3).append(",");
			strbuf.append(bond_4);
			molecular[i] = strbuf.toString();
		}
		return molecular;
	}

	/**
	 * This function is used to generate the atomic features for a single
	 * molecule
	 * 
	 * @param molecule
	 * @return String representation of the atomic features for a single
	 *         molecule where each feature is separated by a comma
	 * @throws CDKException, NoSuchAtomTypeException
	 */
	public String[] generateAtomicFeatures(IAtomContainer molecule) throws CDKException, NoSuchAtomTypeException {
		IAtomContainer mol = molecule;
		AtomDegreeDescriptor degree = new AtomDegreeDescriptor();
		AtomHybridizationDescriptor hy = new AtomHybridizationDescriptor();
		AtomValenceDescriptor va = new AtomValenceDescriptor();
		EffectiveAtomPolarizabilityDescriptor ep = new EffectiveAtomPolarizabilityDescriptor();
		PartialSigmaChargeDescriptor psc = new PartialSigmaChargeDescriptor();
		PartialTChargeMMFF94Descriptor ptc = new PartialTChargeMMFF94Descriptor();
		PiElectronegativityDescriptor pen = new PiElectronegativityDescriptor();
		SigmaElectronegativityDescriptor sen = new SigmaElectronegativityDescriptor();
		StabilizationPlusChargeDescriptor spc = new StabilizationPlusChargeDescriptor();
		int length = mol.getAtomCount();
		int i;
		String[] atomics = new String[length];
		for (i = 0; i < length; i++) {
			IAtom a = mol.getAtom(i);
			String []res=new String[9];
			DescriptorValue d = degree.calculate(a, mol);
			res[0] = d.getValue().toString();
			DescriptorValue h = hy.calculate(a, mol);
			res[1] = h.getValue().toString();
			DescriptorValue v = va.calculate(a, mol);
			res[2] = v.getValue().toString();
			DescriptorValue e = ep.calculate(a, mol);
			res[3] = e.getValue().toString();
			DescriptorValue ps = psc.calculate(a, mol);
			res[4] = ps.getValue().toString();
			DescriptorValue pt = ptc.calculate(a, mol);
			res[5] = pt.getValue().toString();
			DescriptorValue pe = pen.calculate(a, mol);
			res[6]= pe.getValue().toString();
			DescriptorValue se = sen.calculate(a, mol);
			res[7] = se.getValue().toString();
			DescriptorValue sp = spc.calculate(a, mol);
			res[8] = sp.getValue().toString();
			StringBuffer sb=new StringBuffer();
			for(int j=0;j<res.length;j++)
				sb.append(res[j]).append(",");
			sb.deleteCharAt(sb.length() - 1);
			atomics[i] = sb.toString();

		}
		return atomics;
	}

	/**
	 * This function is used to generate all the features for the molecules
	 * contained in the input file
	 * 
	 * @param molf
	 *            : the molecular features of all the molecules
	 * @param atomf
	 *            : the atomic features of all the molecules
	 * @param pathToInputFile
	 *            : path of the file which contains the input molecules
	 * @return the path of the output csv file which contains all the features
	 *         for the molecules
	 * @throws Exception
	 */
	/*public String generateAllFeatures(ArrayList<String> molf,
			ArrayList<ArrayList<String>> atomf, String pathToInputFile) throws Exception {
		ArrayList<ArrayList<String>> result = new ArrayList<ArrayList<String>>();
		File sdfFile = new File(pathToInputFile);
		
		Aromaticity aromaticity = new Aromaticity(ElectronDonation.cdk(), Cycles.all());
		String pathToOutputFile = pathToInputFile.split(".sdf")[0] + ".csv";
		FileOutputStream fos = new FileOutputStream(pathToOutputFile);
		OutputStreamWriter osw = new OutputStreamWriter(fos);
		StringBuffer sb = new StringBuffer();
		int i, j, k;
		for (i = 0; i < 453; i++) {
			sb.append(i + 1).append(",");
		}
		sb.deleteCharAt(sb.length() - 1);
		osw.write(sb.toString() + "\r\n");
		ChemSearcher csearcher = new ChemSearcher();
		LinkedHashMap<String, String> queries = csearcher.getFingerprintPatterns();
		ArrayList<ArrayList<ArrayList<Integer>>> temp = csearcher
				.generateSerialAtomfingerprintToArraylist(sdfFile, queries);
		ArrayList<ArrayList<String>> fp = new ArrayList<ArrayList<String>>();
		int l_1 = temp.size();
		ArrayList<ArrayList<Integer>> array = new ArrayList<ArrayList<Integer>>();
		ArrayList<Integer> list1 = new ArrayList<Integer>();
		for (i = 0; i < l_1; i++) {
			ArrayList<String> list = new ArrayList<String>();
			array = temp.get(i);
			int l_2 = array.size();
			for (j = 0; j < l_2; j++) {
				list1 = array.get(j);
				StringBuffer sb1 = new StringBuffer();
				for (k = 0; k < list1.size(); k++)
					sb1.append(list1.get(k)).append(",");
				sb1.deleteCharAt(sb1.length() - 1);
				list.add(sb1.toString());
			}
			fp.add(list);

		}
		int l = molf.size();
		for (i = 0; i < l; i++) {
			String mol = molf.get(i);
			ArrayList<String> list2 = new ArrayList<String>();
			int l_3 = atomf.get(i).size();
			for (j = 0; j < l_3; j++) {
				for (k = 0; k < l_3; k++) {
					if (k > j) {
						String s1 = atomf.get(i).get(j);
						String s2 = fp.get(i).get(j);
						String s3 = atomf.get(i).get(k);
						String s4 = fp.get(i).get(k);
						String s5 = s1 + "," + s3 + "," + s2 + "," + s4 + "," + mol + ","
								+ "?";
						list2.add(s5);
						osw.write(s5 + "\r\n");
					}
				}

			}
			result.add(list2);
		}
		osw.close();
		fos.close();
		// System.out.println("finish.");
		return pathToOutputFile;
	}*/
	public Instances generateAllFeatures(ArrayList<String> molf,
			ArrayList<ArrayList<String>> atomf, String pathToInputFile) throws Exception {
		ArrayList<ArrayList<String>> result = new ArrayList<ArrayList<String>>();
		File sdfFile = new File(pathToInputFile);
		int i, j, k,m;
		FastVector fvWekaAttributes = new FastVector(453);
		for(i=0;i<452;i++)
		{
			Attribute Attribute = new Attribute(String.valueOf(i));
			fvWekaAttributes.addElement(Attribute);
		}
		 FastVector fvClassVal = new FastVector(1);
		 fvClassVal.addElement("?");
		 Attribute ClassAttribute = new Attribute("Label", fvClassVal);
		 fvWekaAttributes.addElement(ClassAttribute);
		// Create an empty training set
		 Instances isTrainingSet = new Instances("Rel", fvWekaAttributes, 100000);
		 // Set class index
		 isTrainingSet.setClassIndex(452);
	
		ChemSearcher csearcher = new ChemSearcher();
		LinkedHashMap<String, String> queries = ChemSearcher.getFingerprintPatterns();
		ArrayList<ArrayList<ArrayList<Integer>>> temp = csearcher
				.generateSerialAtomfingerprintToArraylist(sdfFile, queries);
		ArrayList<ArrayList<String>> fp = new ArrayList<ArrayList<String>>();
		int l_1 = temp.size();
		ArrayList<ArrayList<Integer>> array = new ArrayList<ArrayList<Integer>>();
		ArrayList<Integer> list1 = new ArrayList<Integer>();
		for (i = 0; i < l_1; i++) {
			ArrayList<String> list = new ArrayList<String>();
			array = temp.get(i);
			int l_2 = array.size();
			for (j = 0; j < l_2; j++) {
				list1 = array.get(j);
				StringBuffer sb1 = new StringBuffer();
				for (k = 0; k < list1.size(); k++)
					sb1.append(list1.get(k)).append(",");
				sb1.deleteCharAt(sb1.length() - 1);
				list.add(sb1.toString());
			}
			fp.add(list);

		}
		int l = molf.size();
		for (i = 0; i < l; i++) {
			String mol = molf.get(i);
			ArrayList<String> list2 = new ArrayList<String>();
			int l_3 = atomf.get(i).size();
			for (j = 0; j < l_3; j++) {
				for (k = 0; k < l_3; k++) {
					if (k > j) {
						String s1 = atomf.get(i).get(j);
						String s2 = fp.get(i).get(j);
						String s3 = atomf.get(i).get(k);
						String s4 = fp.get(i).get(k);
						String s5 = s1 + "," + s3 + "," + s2 + "," + s4 + "," + mol + ","
								+ "?";;
						String []ss=s5.split(",");
						// Create the instance
						int len=ss.length;
						Instance iExample = new Instance(len); 
						 for(m=0;m<len;m++)
						 {
							 if(ss[m].equals("?"))
								 iExample.setValue((Attribute)fvWekaAttributes.elementAt(m),ss[m]);
							 else
								 iExample.setValue((Attribute)fvWekaAttributes.elementAt(m), Double.parseDouble(ss[m]));
						 }
						
						 
						 // add the instance
						 isTrainingSet.add(iExample);
					}
				}

			}
			result.add(list2);
		}
		return isTrainingSet;
	}


	
	
}