package somprediction;


import java.util.ArrayList;
import java.util.LinkedHashMap;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;

import weka.core.Instances;
import featuregeneration.FeatureGenerator;

public class PredictionRunner {

	/**
	 * predict the SOMs for the molecules with a specific enzyme and write out the output
	 * @param input_file
	 * 					:path of the input file
	 * @param model
	 * 				:path of the model file
	 * @param enzyme
	 * 				:name of the enzyme
	 * @return A hashmap which stores the prediction result,
	 * 		    where the key is the input molecule and the value is the corresponding predicted SOMs.
	 * @throws Exception
	 */
	public LinkedHashMap<String, LinkedHashMap<IAtomContainer, ArrayList<Integer>>> runningprediction(String input_file,String model,String enzyme) throws Exception  
	{

		FeatureGenerator fg = new FeatureGenerator();
		IAtomContainerSet MOLS = fg.readFile(input_file);
		ArrayList<String> molf = fg.generateMolecularFeatures(MOLS,input_file);
		ArrayList<ArrayList<String>> atomf = fg.generateAtomicFeatures(MOLS);
		Instances ins = fg.generateAllFeatures(molf, atomf, input_file);
		String model_name="models/"+enzyme+"_SOM.model";
		SOMPredictor sp=new SOMPredictor();
		LinkedHashMap<String, LinkedHashMap<IAtomContainer, ArrayList<Integer>>> result = sp
					.makePrediction(input_file, ins,model_name, MOLS);
		return result;
	}

}
