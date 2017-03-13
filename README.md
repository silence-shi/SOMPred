# SOMPred
A computational tool for CYP450 sites of metabolism prediction

To run the SOM predictor to get the desired result you want, write a script as below:
	PredictionRunner pr=new PredictionRunner();
	LinkedHashMap<String, LinkedHashMap<IAtomContainer, ArrayList<Integer>>> result=pr.runningprediction
				(path,model_name,cyp_name);
where the path is a string path name for your input file, the model is the string name for the prediction model,
and the cyp_name refers to the name of the CYP.
For example:
	PredictionRunner pr=new PredictionRunner();
	LinkedHashMap<String, LinkedHashMap<IAtomContainer, ArrayList<Integer>>> result=pr.runningprediction
				("\Test_sample\example.sdf","CYP1A2_SOM.model","CYP1A2");
