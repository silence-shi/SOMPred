package somprediction;



/**
 * @author Zheng
 *
 */

import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.io.SDFWriter;
import weka.classifiers.Classifier;
import weka.core.Instances;
import weka.filters.Filter;
import weka.filters.unsupervised.attribute.Add;
import weka.filters.unsupervised.attribute.Remove;
import weka.filters.unsupervised.attribute.Standardize;

public class SOMPredictor {
	/**
	 * predict the SOMs for the molecules with a specific enzyme
	 * 
	 * @param pathToInputFile
	 *            : the path of the file which contains all the generated
	 *            features
	 * @param model
	 *            : the path of the pretrained model for the enzyme
	 * @param molecules
	 *            : all the molecules
	 * @return the predicted SOMs for the molecules
	 * @throws Exception
	 */

	public LinkedHashMap<String, LinkedHashMap<IAtomContainer, ArrayList<Integer>>> makePrediction(
			String pathToInputFile, Instances ins, String model, IAtomContainerSet molecules)
			throws Exception {
		LinkedHashMap<String, LinkedHashMap<IAtomContainer, ArrayList<Integer>>> result = new LinkedHashMap<String, LinkedHashMap<IAtomContainer, ArrayList<Integer>>>();
		LinkedHashMap<IAtomContainer, ArrayList<Integer>> map = new LinkedHashMap<IAtomContainer, ArrayList<Integer>>();
		String[] ss = model.split(".model")[0].split("/");
		String s = ss[ss.length - 1].split("_")[0];
		ArrayList<Double> ret = new ArrayList<Double>();
		String s4 = pathToInputFile.split(".sdf")[0] + "_SOM_prediction_" + s + ".sdf";
		File sdfOoutput = new File(s4);
		SDFWriter writer = new SDFWriter(new FileWriter(sdfOoutput));
		Instances unlabeled =ins;
		unlabeled.setClassIndex(unlabeled.numAttributes() - 1);
		Standardize filter = new Standardize();
		filter.setInputFormat(unlabeled);
		Instances newUnlabeled = Filter.useFilter(unlabeled, filter);
		Remove remove = new Remove();
		String[] opts = new String[] { "-R", "last" };
		remove.setOptions(opts);
		remove.setInputFormat(newUnlabeled);
		Instances newData = Filter.useFilter(newUnlabeled, remove);
		Add add = new Add();
		String[] opts1 = new String[] { "-T", "NOM", "-L", "0.0,-1.0,1.0", "-C", "last" };
		add.setOptions(opts1);
		add.setInputFormat(newData);
		Instances filteredData = Filter.useFilter(newData, add);
		filteredData.setClassIndex(filteredData.numAttributes() - 1);
		Classifier cls = (Classifier) weka.core.SerializationHelper.read(model);
		int i, j, length, count = 0;
		IAtomContainer temp;
		Instances labeled = new Instances(filteredData);
		for (i = 0; i < filteredData.numInstances(); i++) {
			double clsLabel = cls.classifyInstance(filteredData.instance(i));
			String ss1 = filteredData.classAttribute().value((int) clsLabel);
			ret.add(Double.parseDouble(ss1));
			labeled.instance(i).setClassValue(clsLabel);
		}
		int l = molecules.getAtomContainerCount();
		for (i = 0; i < l; i++) {
			temp = molecules.getAtomContainer(i);
			length = temp.getAtomCount();
			ArrayList<Double> list = new ArrayList<Double>();
			for (j = 0; j < length * (length - 1) / 2; j++) {
				list.add(ret.get(count));
				count = count + 1;
			}
			ArrayList<Integer> pre = new ArrayList<Integer>();
			pre = predictSOMsByVote(list, length);
			double[] score= rankedSitesByVote( list, length);
			ArrayList<Integer> rank_index= getIndexForSitesByVote(list, length);
			map.put(temp, pre);
			StringBuffer sb = new StringBuffer();
			for (j = 0; j < pre.size(); j++) {
				sb.append(String.valueOf(pre.get(j) + 1)).append(" ");
			}
			temp.setProperty("Prediction_SOM_"+s, sb.toString());
			sb.setLength(0);
			for (j = 0; j < rank_index.size(); j++) {
				sb.append(String.valueOf(rank_index.get(j) + 1)).append(" ");
			}
			temp.setProperty("Ranking_Sites_"+s, sb.toString());
			sb.setLength(0);
			for (j = 0; j < score.length; j++) {
				sb.append(String.valueOf(score[rank_index.get(j)])).append(" ");
			}
			temp.setProperty("Prediction_Scores_"+s, sb.toString());
			writer.write(temp);
			
		}
		writer.close();
		result.put(s, map);
		return result;
	}	
	
	/**
	 * predicted SOMs with voting preference aggregation
	 * 
	 * @param list
	 *            : predicted preference relation
	 * @param l
	 *            : length of a molecule
	 * @return predicted SOMs for one molecule
	 */
	public ArrayList<Integer> predictSOMsByVote(ArrayList<Double> list, int l) {

		ArrayList<Integer> result = new ArrayList<Integer>();
		double[][] a = new double[l][l];
		int i, j, count = 0;
		double sum;
		for (i = 0; i < l - 1; i++) {
			for (j = 0; j < l; j++) {
				if (i == j)
					a[i][j] = 0;
				else {
					if (j > i) {
						a[i][j] = list.get(count);
						count = count + 1;
					}
					else {
						a[i][j] = 0 - a[j][i];
					}
				}
			}
		}
		for (j = 0; j < l; j++)
			a[l - 1][j] = 0 - a[j][l - 1];
		for (i = 0; i < l; i++) {
			sum = 0;
			for (j = 0; j < l; j++) {
				sum = sum + a[i][j];
			}
			if (sum > 0) {

				result.add(i);
			}
		}
		return result;
	}
	/**
	 * 
	 * @param list
	 *            : predicted preference relation
	 * @param l
	 *            : length of a molecule
	 * @return each site with its score for one molecule
	 */
	public double[] rankedSitesByVote(ArrayList<Double> list, int l)
	{
		double[]result= new double[l];
		double[][] a = new double[l][l];
		int i, j, count = 0;
		double sum;
		for (i = 0; i < l - 1; i++) {
			for (j = 0; j < l; j++) {
				if (i == j)
					a[i][j] = 0;
				else {
					if (j > i) {
						a[i][j] = list.get(count);
						count = count + 1;
					}
					else {
						a[i][j] = 0 - a[j][i];
					}
				}
			}
		}
		for (j = 0; j < l; j++)
			a[l - 1][j] = 0 - a[j][l - 1];
		for (i = 0; i < l; i++) {
			sum = 0;
			for (j = 0; j < l; j++) {
				sum = sum + a[i][j];
			}
			result[i]=sum;
		}
		
		return result;
	}
	/**
	 * 
	 * @param list
	 *            : predicted preference relation
	 * @param l
	 *            : length of a molecule
	 * @return the ranked index for each site using its score for one molecule
	 */
	public ArrayList<Integer> getIndexForSitesByVote(ArrayList<Double> list, int l)
	{
		double []score=rankedSitesByVote( list,  l);
		double []temp=Arrays.copyOf(score, score.length);
		Arrays.sort(temp);
		ArrayList<Integer> index=new ArrayList<Integer> ();
		int i,j;
		boolean []visited=new boolean[l];
		for(i=l-1;i>=0;i--)
		{
			for(j=0;j<l;j++)
			{
				if(!visited[j]&&score[j]==temp[i])
				{
					index.add(j);
					visited[j]=true;
				}
			}
		}
		return index;
	}
}
