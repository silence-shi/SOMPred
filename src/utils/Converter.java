package utils;


/**
 * @author Zheng Shi
 *
 */

import java.io.File;

import weka.core.Instances;
import weka.core.converters.ArffSaver;
import weka.core.converters.CSVLoader;

public class Converter {
	/**
	 * convert .csv to arff file
	 * 
	 * @param input_path
	 *            : the path of the input csv file
	 * @return the path of the output arff file
	 * @throws Exception
	 */
	public String convert(String input_path) throws Exception {
		System.out.println("start");
		// load CSV
		CSVLoader loader = new CSVLoader();
//		File inFile = new File(input_path);
		loader.setSource(new File(input_path));
		String[] opts = new String[] {  "-N last"};
				
		System.out.println("Setting loader options");
		loader.setOptions(opts);
		
		for(String s : loader.getOptions()){
			System.out.println(s);
		}
			
		System.out.println("Getting Dataset");
		Instances data = loader.getDataSet();
		String output_path = input_path.split(".csv")[0] + ".arff";
		// String output_path = "D:/metabolism/data/new data/HMDB.arff";
		// save ARFF
		
		ArffSaver saver = new ArffSaver();
		System.out.println("setting instances");
		saver.setInstances(data);
		saver.setFile(new File(output_path));
		// saver.setDestinatio n(new File(output_path));
		System.out.println("Saving arff");
		saver.writeBatch();
		 System.out.println("finished the conversion");
		return output_path;
	}

}
