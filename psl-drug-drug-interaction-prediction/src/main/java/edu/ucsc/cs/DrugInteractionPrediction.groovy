package edu.ucsc.cs;

import java.util.Collections;
import java.util.Iterator;
import java.util.Random;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


import org.linqs.psl.application.inference.MPEInference
import org.linqs.psl.application.learning.weight.maxlikelihood.MaxLikelihoodMPE
import org.linqs.psl.config.*
import org.linqs.psl.core.*
import org.linqs.psl.core.inference.*
import org.linqs.psl.database.DataStore
import org.linqs.psl.database.Database
import org.linqs.psl.database.DatabasePopulator
import org.linqs.psl.database.Partition
import org.linqs.psl.database.rdbms.RDBMSDataStore
import org.linqs.psl.database.rdbms.driver.H2DatabaseDriver
import org.linqs.psl.database.rdbms.driver.H2DatabaseDriver.Type
import org.linqs.psl.utils.evaluation.result.*
import org.linqs.psl.utils.evaluation.statistics.RankingScore
import org.linqs.psl.utils.evaluation.statistics.SimpleRankingComparator
import org.linqs.psl.groovy.*
import org.linqs.psl.model.atom.GroundAtom;
import org.linqs.psl.model.atom.QueryAtom
import org.linqs.psl.model.term.*;
import org.linqs.psl.model.rule.*;
import org.linqs.psl.model.rule.arithmetic.UnweightedArithmeticRule;
import org.linqs.psl.model.weight.*;
import org.linqs.psl.model.formula.*;
import org.linqs.psl.model.function.*;
import org.linqs.psl.model.kernel.*;
import org.linqs.psl.model.predicate.Predicate
import org.linqs.psl.database.*
import org.linqs.psl.utils.dataloading.InserterUtils;


import com.google.common.collect.Iterables

class DrugInteractionPrediction {
	Logger log = LoggerFactory.getLogger(this.class);

	/**
	 * DEFINES k-fold cross-validation
	 * stratified cross-validation refers to the process of rearranging
	 * the data as to ensure each fold is a good representative of the whole.
	 */
	class DrugInteractionFold {
		public int cvFold;
		public int wlFold;

		public Partition cvTruth;
		public Partition cvTest;
		public Partition cvTrain;

		public Partition wlTruth;
		public Partition wlTest;
		public Partition wlTrain;

		public Partition wlSim;
		public Partition cvSim;

		DrugInteractionFold(Integer cvFold, Integer wlFold){
			this.cvFold = cvFold;
			this.wlFold = wlFold;
		}

	}

	/**
	 * DEFINES experiment parameter definition
	 */
	class DrugInteractionExperiment {
		public ConfigManager cm;
		public ConfigBundle cb;
		public double initialWeight;
		public int numFolds;
		public int numDrugs;
		public String interaction_type;
		public String experiment_name;
		public boolean doWeightLearning;
		public boolean createNewDataStore;
		public String base_dir;
		public String interactions_dir;
		public String similarities_dir;
		public int blockingType;
		public int blockingNum;

		Set<Predicate> closedPredicatesInference;
		Set<Predicate> closedPredicatesTruth;

	}

	/**
	 * DEFINES resulting parameters
	 * auc : Area under ROC Curve
	 * aupr : Area under PR Curve
	 * method: First ensure that the 5 folds are all stratified
	 * (e.g. they all are equally distributed with respect to the originak sample distribution).
	 * Then, after training some classifier, calculate the AUC the
	 * classifier's prediction achieves for each of the 5 folds.
	 */
	class DrugInteractionEvalResults {
    	// Creating the variables to save the results of each fold
	    public double[] AUC_Folds;
	    public double[] AUPR_P_Folds;
	    public double[] AUPR_N_Folds;

	    DrugInteractionEvalResults(Integer folds){
	    	this.AUC_Folds = new double[folds+1];
	    	this.AUPR_P_Folds = new double[folds+1];
	    	this.AUPR_N_Folds = new double[folds+1];
    	}
    }


	/**
	 * DEFINES config information for the experiment
	 * returns `config` object with all information
	 */
    def setupConfig(numFolds, numDrugs, experimentName, blocking_k, blocking_type, interaction_type, createDS){
    	DrugInteractionExperiment config = new DrugInteractionExperiment();
    	config.cm = ConfigManager.getManager();
    	config.cb = config.cm.getBundle("fakhraei_sridhar_bioinformatics");

    	config.initialWeight = 5.0;
    	config.numFolds = numFolds;
    	config.numDrugs = numDrugs;
    	config.doWeightLearning = true;
    	config.createNewDataStore = createDS;

    	config.blockingNum = blocking_k;
    	config.interaction_type = interaction_type;
    	config.blockingType = blocking_type;
    	config.experiment_name = experimentName;


    	config.base_dir = 'data'+java.io.File.separator;
    	config.interactions_dir = config.base_dir + config.experiment_name + java.io.File.separator + config.interaction_type + java.io.File.separator;
    	config.similarities_dir = config.blockingType + '_' + config.blockingNum.toString() + java.io.File.separator;

    	return config

    }

	/**
	 * DEFINES datastore location
	 */
    def setupDataStore(config) {
    	String dbpath = "/tmp/psl_fsbio_collective";
    	DataStore data = new RDBMSDataStore(new H2DatabaseDriver(Type.Disk, dbpath, config.createNewDataStore), config.cb);

    	return data
    }

	/**
	 * DEFINES PSL rules
	 * rules are based on similarity condition,
	 * that is, if two drugs are similar and one other drug is similar to one of previous pair,
	 * then all drugs are similar and hence interaction is possible
	 */
    def defineModel(config, data, m){

    	m.add predicate : "ATCSimilarity" , types:[ConstantType.UniqueID, ConstantType.UniqueID]
    	m.add predicate : "SideEffectSimilarity" , types:[ConstantType.UniqueID, ConstantType.UniqueID]
    	m.add predicate : "GOSimilarity" , types:[ConstantType.UniqueID, ConstantType.UniqueID]
    	m.add predicate : "ligandSimilarity" , types:[ConstantType.UniqueID, ConstantType.UniqueID]
    	m.add predicate : "chemicalSimilarity" , types:[ConstantType.UniqueID, ConstantType.UniqueID]
    	m.add predicate : "seqSimilarity" , types:[ConstantType.UniqueID, ConstantType.UniqueID]
    	m.add predicate : "distSimilarity" , types:[ConstantType.UniqueID, ConstantType.UniqueID]
    	m.add predicate : "interacts" , types:[ConstantType.UniqueID, ConstantType.UniqueID]
    	m.add predicate: "validInteraction" , types:[ConstantType.UniqueID, ConstantType.UniqueID]

    	m.add rule : (ATCSimilarity(D1, D2) & interacts(D1, D3) & validInteraction(D1, D3) & validInteraction(D2, D3) & (D2 - D3) & (D1 - D2))>>interacts(D2, D3) , weight:config.initialWeight
		m.add rule : (SideEffectSimilarity(D1, D2) & interacts(D1, D3) & validInteraction(D1, D3) & validInteraction(D2, D3) & (D2 - D3) & (D1 - D2))>>interacts(D2, D3) , weight:config.initialWeight
		m.add rule : (GOSimilarity(D1, D2) & interacts(D1, D3) & validInteraction(D1, D3) & validInteraction(D2, D3) & (D2 - D3) & (D1 - D2))>>interacts(D2, D3) , weight:config.initialWeight
		m.add rule : (ligandSimilarity(D1, D2) & interacts(D1, D3) & validInteraction(D1, D3) & validInteraction(D2, D3) & (D2 - D3) & (D1 - D2))>>interacts(D2, D3) , weight:config.initialWeight
		m.add rule : (chemicalSimilarity(D1, D2) & interacts(D1, D3) & validInteraction(D1, D3) & validInteraction(D2, D3) & (D2 - D3) & (D1 - D2))>>interacts(D2, D3) , weight:config.initialWeight
		m.add rule : (seqSimilarity(D1, D2) & interacts(D1, D3) & validInteraction(D1, D3) & validInteraction(D2, D3) & (D2 - D3) & (D1 - D2))>>interacts(D2, D3) , weight:config.initialWeight
		m.add rule : (distSimilarity(D1, D2) & interacts(D1, D3) & validInteraction(D1, D3) & validInteraction(D2, D3) & (D2 - D3) & (D1 - D2))>>interacts(D2, D3) , weight:config.initialWeight

		//prior
		m.add rule : validInteraction(D1,D2) >> ~interacts(D1,D2),  weight : config.initialWeight

		m.add rule : "interacts(D1, D2) = interacts(D2, D1) ."

		config.closedPredicatesInference = [validInteraction, ATCSimilarity, distSimilarity, seqSimilarity, ligandSimilarity, GOSimilarity, SideEffectSimilarity, chemicalSimilarity];
		config.closedPredicatesTruth = [interacts];


	}


	/**
	 * DEFINES loads data, creates partitions
	 */
	def loadData(data,config) {

	    // Loading the data
	    // ================

	    // Creating the partition to read the data
	    Partition readSimilarities =  data.getPartition("all_sims");

	    def insert;
	    def simDir = config.base_dir + config.experiment_name + java.io.File.separator;

	    for (Predicate p : [ATCSimilarity, distSimilarity, seqSimilarity, ligandSimilarity, GOSimilarity, SideEffectSimilarity, chemicalSimilarity])
	    {
	    	log.debug("Reading " + p.getName());
	    	insert = data.getInserter(p, readSimilarities)
	    	InserterUtils.loadDelimitedDataTruth(insert, simDir +p.getName().toLowerCase()+".csv")
	    }
	}


	/**
	 * DEFINES similarity data loading,
	 * reads triad target similarity data
	 * wlfold : k-folds for weight learning
	 * cvfold : k-folds for cross-validation
	 */
	def loadSimilarityData(config, data, df){
		def cvfold = df.cvFold;
		def wlfold = df.wlFold

		log.debug("\n-------------------");
		def timeNow = new Date();
		log.debug("Fold "+ cvfold +" Start: "+timeNow);

		def current_sim_dir = config.interactions_dir + cvfold + java.io.File.separator + config.similarities_dir;

		if (!config.createNewDataStore){
			//clearing the partitions in case already in datastore
			data.deletePartition(data.getPartition(df.wlSim));
			data.deletePartition(data.getPartition(df.cvSim));
		}

		df.wlSim = data.getPartition("wlsim" + wlfold)
		df.cvSim = data.getPartition("cvsim" + cvfold)

		// Reading triad target similarities from file
		for (Predicate p : [ATCSimilarity, distSimilarity, seqSimilarity, ligandSimilarity, GOSimilarity, SideEffectSimilarity, chemicalSimilarity])
		{
			log.debug("Reading " + p.getName());
			def insert = data.getInserter(p, df.cvSim)
			InserterUtils.loadDelimitedDataTruth(insert, current_sim_dir+p.getName().toLowerCase()+"_cv.csv")
			insert = data.getInserter(p, df.wlSim)
			InserterUtils.loadDelimitedDataTruth(insert, current_sim_dir+p.getName().toLowerCase()+"_wl.csv")

		}

	}

	/**
	 * DEFINES labelling
	 */
	def resetDataForFold(config, data, df){
		def cvFold = df.cvFold;
		def wlFold = df.wlFold;

		df.wlSim = data.getPartition("wlsim" + wlFold);
		df.cvSim = data.getPartition("cvsim" + cvFold);

		df.cvTruth =  data.getPartition("cvlabels" + cvFold); // Labels for the Cross-validation hold-outs
		df.cvTest =  data.getPartition("cvwrite" + cvFold); // Partition to write the predictions in
		df.cvTrain =  data.getPartition("cvread" + cvFold); // Observed training data for the training (i.e., all data minus hold-outs)

		df.wlTruth =  data.getPartition("wllabels" + wlFold); // Labels for the weight Learning
		df.wlTrain =  data.getPartition("wlread"  + wlFold); // Training data for Weight Learning
		df.wlTest =  data.getPartition("wlwrite" + wlFold); // Partition to write prediction in for Weight Learning

	}

	/**
	 * DEFINES
	 *
	 */
	def blockSimilarityData(config, data, df){
		def cvfold = df.cvFold;
		def wlfold = df.wlFold;

		TriadBlocking triadBlock = new TriadBlocking();

		log.debug("\n-------------------");
		def timeNow = new Date();
		log.debug("Fold "+ cvfold +" Start: "+timeNow);

		if (!config.createNewDataStore){
			//clearing the partitions in case already in datastore
			data.deletePartition(data.getPartition(df.wlSim));
			data.deletePartition(data.getPartition(df.cvSim));
		}

		df.cvSim = data.getPartition('cvsim' + cvfold);
		df.wlSim = data.getPartition('wlsim' + wlfold);

		Partition readSimilarities = data.getPartition("all_sims");
		Database getSimilarities = data.getDatabase(data.getPartition("simDummy" + cvfold), readSimilarities);

		//Dummy dbs for wl train, wl test, cv train, cv test for more involved blocking techniques - DS

		Database wlTrainLinks = data.getDatabase(data.getPartition("wlTrainDummy" + wlfold ), df.wlTrain);
		Database wlTestLinks = data.getDatabase(data.getPartition("wltestdummy" + wlfold), df.wlTruth);

		Database cvTrainLinks = data.getDatabase(data.getPartition("cvTrainDummy" + cvfold ), df.cvTrain);
		Database cvTestLinks = data.getDatabase(data.getPartition("cvtestdummy" + cvfold), df.cvTruth);

		Set<GroundAtom> wlTrainingLinks = Queries.getAllAtoms(wlTrainLinks, interacts);
		Set<GroundAtom> cvTrainingLinks = Queries.getAllAtoms(cvTrainLinks, interacts);
		Set<GroundAtom> wlTestingLinks = Queries.getAllAtoms(wlTestLinks, interacts);
		Set<GroundAtom> cvTestingLinks = Queries.getAllAtoms(cvTestLinks, interacts);

		wlTestingLinks.removeAll(cvTestingLinks);

		// Reading triad target similarities from file
		for (Predicate p : [ATCSimilarity, distSimilarity, seqSimilarity, ligandSimilarity, GOSimilarity, SideEffectSimilarity, chemicalSimilarity])
		{
			log.debug("Working on " + p.getName());

			Set<GroundAtom> allSimilarities = Queries.getAllAtoms(getSimilarities, p); //SF

			def k = config.blockingNum;
			def cvBlockedSim, wlBlockedSim;
			switch(config.blockingType){
				case 1:
					cvBlockedSim = triadBlock.observedLinksInTraining(allSimilarities, cvTrainingLinks, 5, k );
					wlBlockedSim = triadBlock.observedLinksInTraining(allSimilarities, wlTrainingLinks, 5, k );
					break;
				case 2:
					cvBlockedSim = triadBlock.allLinksInTesting(allSimilarities, cvTestingLinks, 5, k );
					wlBlockedSim = triadBlock.allLinksInTesting(allSimilarities, wlTestingLinks, 5, k );
					break;
				case 3:
					cvBlockedSim = triadBlock.testTermsWithObservedLinks(allSimilarities, cvTrainingLinks, cvTestingLinks, 5, k );
					wlBlockedSim = triadBlock.testTermsWithObservedLinks(allSimilarities, wlTrainingLinks, wlTestingLinks, 5, k );
					break;

				case 4:
					cvBlockedSim = triadBlock.observedTrainWithTestTerm(allSimilarities, cvTrainingLinks, cvTestingLinks, 5, k );
					wlBlockedSim = triadBlock.observedTrainWithTestTerm(allSimilarities, wlTrainingLinks, wlTestingLinks, 5, k );
					break;
				case 5:
					cvBlockedSim = triadBlock.knnBlockingSet(true,k,0, triadBlock.GetQueryTerms(allSimilarities, 0), allSimilarities);
					wlBlockedSim = triadBlock.knnBlockingSet(true,k,0, triadBlock.GetQueryTerms(allSimilarities, 0), allSimilarities);
			}

			// Inserting them into the partition - SF
			def cvInsert = data.getInserter(p, df.cvSim)
			def wlInsert = data.getInserter(p, df.wlSim)
			Iterator<GroundAtom> itrSimilarities = cvBlockedSim.iterator();
			while (itrSimilarities.hasNext()){
				GroundAtom itemSimilarity = itrSimilarities.next();
				cvInsert.insertValue(itemSimilarity.getValue(), itemSimilarity.getArguments());
			}

			itrSimilarities = wlBlockedSim.iterator();
			while (itrSimilarities.hasNext()){
				GroundAtom itemSimilarity = itrSimilarities.next();
				wlInsert.insertValue(itemSimilarity.getValue(), itemSimilarity.getArguments());
			}

		}
		getSimilarities.close(); //SF

		wlTrainLinks.close();
		wlTestLinks.close();
		cvTrainLinks.close();
		cvTestLinks.close();

	}

	/**
	 * DEFINES interaction data loading, going through folds data
	 */
	def loadInteractionsData(config, data, df, isSupervisedMode){
		def cvfold = df.cvFold;
		def wlfold = df.wlFold;
		def interactions_file = 'interacts.csv'
		def interactions_ids = 'interactsids.csv'
		def interactions_positive = 'interacts_positives.csv'
		def interactions_negative = 'interacts_negatives.csv'

		def cvlabels = 'cvlabels' + cvfold
		def cvwrite = 'cvwrite' + cvfold
		def cvread = 'cvread' + cvfold

		def wllabels = 'wllabels' + wlfold
		def wlwrite = 'wlwrite' + wlfold
		def wlread = 'wlread' + wlfold

		log.debug("\n-------------------");
		def timeNow = new Date();
		log.debug("Fold "+ cvfold +" Start: "+timeNow);

		if (!config.createNewDataStore){
		//clearing from partition if it already exists
			data.deletePartition(data.getPartition(cvlabels));
			data.deletePartition(data.getPartition(cvwrite));
			data.deletePartition(data.getPartition(cvread));
			data.deletePartition(data.getPartition(wlwrite));
			data.deletePartition(data.getPartition(wlread));
			data.deletePartition(data.getPartition(wllabels));
		}

		// Setting up the partitions for final model
		df.cvTruth =  data.getPartition(cvlabels); // Labels for the Cross-validation hold-outs
		df.cvTest =  data.getPartition(cvwrite); // Partition to write the predictions in
		df.cvTrain =  data.getPartition(cvread); // Observed training data for the training (i.e., all data minus hold-outs)

		//Setting up the the partitions for weight learning
		df.wlTruth =  data.getPartition(wllabels); // Labels for the weight Learning
		df.wlTrain =  data.getPartition(wlread); // Training data for Weight Learning
		df.wlTest =  data.getPartition(wlwrite); // Partition to write prediction in for Weight Learning

		// Setting up the inserters
		def insertWLTrain = data.getInserter(interacts, df.wlTrain);
		def insertWLLabels = data.getInserter(interacts, df.wlTruth);
		def insertWLValid = data.getInserter(validInteraction, df.wlTrain);

		def insertCVValid = data.getInserter(validInteraction, df.cvTrain);
		def insertCVTrain = data.getInserter(interacts, df.cvTrain);
		def insertCVLabels = data.getInserter(interacts, df.cvTruth);

		def insertWLTest = data.getInserter(interacts, df.wlTest);
		def insertCVTest = data.getInserter(interacts, df.cvTest);

		// Reading the interactions and setting the data for current fold

		log.debug("\nReading INTERACTS files for fold "+ cvfold +" ");

		// Reading all the other folds as training data
		for (int j=1;j<=config.numFolds;j++)
		{
			def current_interactions_dir = config.interactions_dir + j + java.io.File.separator;
			log.debug(current_interactions_dir);
			if ((j!=cvfold) && (j!=wlfold))
			{

				InserterUtils.loadDelimitedData(insertWLValid, current_interactions_dir + interactions_ids);
				InserterUtils.loadDelimitedData(insertCVValid, current_interactions_dir + interactions_ids);

				if(isSupervisedMode){

					InserterUtils.loadDelimitedDataTruth(insertCVTrain, current_interactions_dir + interactions_file);
					InserterUtils.loadDelimitedDataTruth(insertWLTrain, current_interactions_dir + interactions_file);

				}else{

					InserterUtils.loadDelimitedDataTruth(insertCVTrain, current_interactions_dir + interactions_positive);
					InserterUtils.loadDelimitedDataTruth(insertWLTrain, current_interactions_dir + interactions_positive);

					//for semi-supervised setting
					InserterUtils.loadDelimitedData(insertWLTest, current_interactions_dir + interactions_negative);
					InserterUtils.loadDelimitedData(insertCVTest, current_interactions_dir + interactions_negative);
				}
			}
		}

		def train_interactions_dir = config.interactions_dir + wlfold + java.io.File.separator;

		InserterUtils.loadDelimitedData(insertWLValid, train_interactions_dir + interactions_ids);
		InserterUtils.loadDelimitedData(insertCVValid, train_interactions_dir + interactions_ids);

		// Adding the weight learning held-out to the final training data
		// Reading the weight learning held-out as the the labels - use positives file if in semi supervised mode

		if(isSupervisedMode){

			InserterUtils.loadDelimitedDataTruth(insertCVTrain, train_interactions_dir + interactions_file);

			InserterUtils.loadDelimitedDataTruth(insertWLLabels, train_interactions_dir + interactions_file);

			InserterUtils.loadDelimitedData(insertWLTest, train_interactions_dir + interactions_ids);



		}else{
			InserterUtils.loadDelimitedDataTruth(insertCVTrain, train_interactions_dir + interactions_positive);
			InserterUtils.loadDelimitedData(insertCVTest, train_interactions_dir + interactions_negative);

			InserterUtils.loadDelimitedDataTruth(insertWLLabels, train_interactions_dir + interactions_positive);
			InserterUtils.loadDelimitedData(insertWLTest, train_interactions_dir + interactions_ids);
		}

		// Reading the cross-validation held out as labels for weight learning and also into the ignored_interacts_DT.
		// Weight learning will not be able to see these labels because of they will be ignored and the body of their rules will be 0.

		def holdout_interactions_dir = config.interactions_dir + cvfold + java.io.File.separator;

		InserterUtils.loadDelimitedData(insertCVValid, holdout_interactions_dir + interactions_ids);

		// Reading the cross-validation held out as labels for the final model
		InserterUtils.loadDelimitedDataTruth(insertCVLabels, holdout_interactions_dir + interactions_file);

		//Insert drug pair ids into write partitions
		InserterUtils.loadDelimitedData(insertCVTest, holdout_interactions_dir + interactions_ids);
	}


	/**
	 * DEFINES weight learning process
	 * uses MaxLikelihoodMPE algo for it
	 */
	def learnWeights(m, data, config, df, learnerOption){

	    /// Weight Learning
	    // ===============
	    def fold = df.cvFold;

	    def timeNow = new Date();
	    log.debug("Fold "+ fold +" Weight Learning: " + timeNow);

	    if(fold > 1){
            for(Rule r : m.getRules()){
            	if ((r instanceof WeightedRule)){
            		Weight w = new PositiveWeight(config.initialWeight);
                	r.setWeight(w);
            	}
            }
        }

	    Database dbWLTrain = data.getDatabase(df.wlTest, config.closedPredicatesInference, df.wlTrain, df.wlSim);
	    Database dbWLLabels = data.getDatabase(df.wlTruth, config.closedPredicatesTruth);


	    def wLearn = new MaxLikelihoodMPE(m,dbWLTrain,dbWLLabels,config.cb);
	    wLearn.learn();
	    wLearn.close();



	    dbWLTrain.close();
	    dbWLLabels.close();

	    // Printing the new weights
	    println m;
	}

	// Inferring
	// =========
	/**
	* DEFINES inference on data
	* uses MPE Inference algorithm
	*/
	def runInference(m,data,config,df, testPartition, evidencePartition, simPartition) {
		def fold = df.cvFold.toString();
		def timeNow = new Date();
		log.debug("Fold "+ fold +" Inferring: " + timeNow);

		Database dbCVTrain = data.getDatabase(testPartition, config.closedPredicatesInference , evidencePartition, simPartition);

		MPEInference mpe = new MPEInference(m, dbCVTrain, config.cb);
		mpe.mpeInference();
		mpe.close();
		mpe.finalize();
		dbCVTrain.close();

		timeNow = new Date();
		log.debug("Fold "+ fold +" End: "+timeNow);
		log.debug("-------------------");
	}

	/**
	 * DEFINES evaluation of result
	 *	uses AUPRC, NegAUPRC & Area-under-ROC
	 */
	def evaluateResults(data,config,df,de, truthPartition, predictionPartition, writeToFile){

	    //Begin Evaluate
	    //==============

	    def fold = df.cvFold;

	    System.out.println("Evaluating fold: " + fold.toString());

	    def labelsDB = data.getDatabase(truthPartition, config.closedPredicatesTruth)
	    Database predictionsDB = data.getDatabase(data.getPartition("dummy_cv_preds_"+fold), predictionPartition)
	    def comparator = new SimpleRankingComparator(predictionsDB);
	    comparator.setBaseline(labelsDB);

	    if(writeToFile){
	    	OutputWriter writer = new OutputWriter(predictionsDB, interacts, 'collective_'+config.interaction_type, "allFolds", "psl");
	    	writer.outputToFile(true);
	    }

	    // Choosing what metrics to report
	    def metrics = [RankingScore.AUPRC, RankingScore.NegAUPRC,  RankingScore.AreaROC]
	    double [] score = new double[metrics.size()]

	    try {
	    	for (int i = 0; i < metrics.size(); i++) {
	    		comparator.setRankingScore(metrics.get(i))
	    		score[i] = comparator.compare(interacts)
	    	}
	      //Storing the performance values of the current fold
	      de.AUPR_P_Folds[fold]=score[0];
	      de.AUPR_N_Folds[fold]=score[1];
	      de.AUC_Folds[fold]=score[2];

	      System.out.println("Area under positive-class PR curve: " + score[0])
	      System.out.println("Area under negetive-class PR curve: " + score[1])
	      System.out.println("Area under ROC curve: " + score[2])
	      System.out.println("-------------------");
	      } catch (ArrayIndexOutOfBoundsException e) {
	      	System.out.println("No evaluation data! Terminating!");
	      }

	      predictionsDB.close();
	      labelsDB.close();
	  }

	/**
	 * DEFINES output print on `/output/..` files
	 */
	def outputResultsCSV(interactionType, experimentName, array, resultName){
		BufferedWriter writer = null;
		String resultsFile = "./output/psl/" + resultName + '_' + interactionType + '_' + experimentName;

		try {
			writer = new BufferedWriter(new FileWriter(resultsFile));
			StringBuilder output = new StringBuilder();
			for (double d : array){
				output.append(d+',');
			}
			writer.append(output.toString());
			writer.flush();

		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	/**
	 * DEFINES processing config file and running experiment using it
	 */
	def runExperiment(numDrugs, numFolds, experimentName, blockingNum, blockingType, interactionType, createDS){

		def config = this.setupConfig(numFolds, numDrugs, experimentName, blockingNum, blockingType, interactionType, createDS);
		def data = this.setupDataStore(config);
		PSLModel m = new PSLModel(this, data);

		this.defineModel(config, data, m);
		DrugInteractionEvalResults de = new DrugInteractionEvalResults(config.numFolds);

		if(config.createNewDataStore){
			this.loadData(data, config);
		}

		for(int fold = 1; fold <= config.numFolds; fold++){
			def wl = (fold % config.numFolds) + 1
			DrugInteractionFold df = new DrugInteractionFold(fold, wl);

			if(config.createNewDataStore){
				this.loadInteractionsData(config,data, df, true);
				this.blockSimilarityData(config, data, df);
			}
			else{
				this.resetDataForFold(config, data, df);
			}

			if(config.doWeightLearning){
				this.learnWeights(m,data,config,df, 1);
			}
			this.runInference(m,data,config,df, df.cvTest, df.cvTrain, df.cvSim);

			this.evaluateResults(data,config,df,de, df.cvTruth, df.cvTest, false);

		}

		data.close();

		return de;

	}

	/**
	 * DEFINES calling all functions
	 */
	static void main(args){

		String[] experiments = ['all_dataset2', 'all_dataset1']
		String[] interactionTypes = ['all'];

		for(String exp: experiments){
			def numDrugs = 315;

			if(exp == "all_dataset1"){
				interactionTypes = ["crd", "ncrd"];
				numDrugs = 807;
			}

			for(String it: interactionTypes){
				def blockingNum = 15;
				def blockingType = 5;


				if (it == "crd"){
					blockingType = 3;
				}


				System.out.println("Working on interaction type " + it + " with blocking type " + blockingType);
				def dip = new DrugInteractionPrediction();
				def evaluation = dip.runExperiment(numDrugs, 10, exp, blockingNum, blockingType, it, true);

				dip.outputResultsCSV(it, exp, evaluation.AUC_Folds, "auc");
				dip.outputResultsCSV(it, exp, evaluation.AUPR_N_Folds, "auprNeg");
				dip.outputResultsCSV(it, exp, evaluation.AUPR_P_Folds, "auprPos");

			}
		}
	}
}
