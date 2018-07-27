package pt.uminho.ceb.biosystems.merlin.local.alignments.core;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Observable;
import java.util.Observer;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;

import javax.xml.bind.JAXBContext;

import org.apache.axis2.i18n.ProjectResourceBundle;
import org.biojava.nbio.core.sequence.template.AbstractSequence;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import pt.uminho.ceb.biosystems.merlin.bioapis.externalAPI.ncbi.CreateGenomeFile;
import pt.uminho.ceb.biosystems.merlin.local.alignments.core.ModelMerge.BlastAlignment;
import pt.uminho.ceb.biosystems.merlin.local.alignments.core.ModelMerge.ModelAlignments;
import pt.uminho.ceb.biosystems.merlin.utilities.Enumerators.AlignmentPurpose;
import pt.uminho.ceb.biosystems.merlin.utilities.Enumerators.AlignmentScoreType;
import pt.uminho.ceb.biosystems.merlin.utilities.Enumerators.Method;
import pt.uminho.ceb.biosystems.merlin.utilities.blast.ncbi_blastparser.BlastOutput;
import pt.uminho.ceb.biosystems.merlin.utilities.containers.capsules.AlignmentCapsule;
import pt.uminho.ceb.biosystems.merlin.utilities.io.FileUtils;

/**
 * @author ODias
 *
 */
public class RunSimilaritySearch extends Observable implements Observer {

	private Map<String, AbstractSequence<?>> staticGenesSet;
	private AtomicBoolean cancel;
	private AtomicInteger counter;
	private AtomicInteger querySize;
	private double similarity_threshold;
	private Method method;
	private ConcurrentHashMap<String, AbstractSequence<?>> querySequences;
	private List<String> annotatedGenes;
	private ConcurrentLinkedQueue<String> sequencesWithoutSimilarities;
	private String ec_number;
	private Map<String, Set<Integer>> modules;
	private Map<String, Set<String>> closestOrthologs;
	private int referenceTaxonomyScore;
	private Map<String, Integer> kegg_taxonomy_scores;
	private double referenceTaxonomyThreshold;
	private boolean compareToFullGenome;
	private AlignmentScoreType alignmentScoreType;
//	private String tcdbFastaFilePath;
	private String subjectFastaFilePath;
	private boolean gapsIdentification;
	private String workspaceTaxonomyFolderPath;
	
	final static Logger logger = LoggerFactory.getLogger(RunSimilaritySearch.class);


	/**
	 * Run similarity searches constructor.
	 * 
	 * @param dbAccess
	 * @param staticGenesSet
	 * @param minimum_number_of_helices
	 * @param similarity_threshold
	 * @param method
	 * @param querySequences
	 * @param cancel
	 * @param querySize
	 * @param counter
	 * @param project_id
	 * @param alignmentScoreType
	 * @throws Exception
	 */
	public RunSimilaritySearch(Map<String, AbstractSequence<?>> staticGenesSet, double similarity_threshold, Method method, ConcurrentHashMap<String, AbstractSequence<?>> querySequences, 
			AtomicBoolean cancel, AtomicInteger querySize, AtomicInteger counter, AlignmentScoreType alignmentScoreType) throws Exception {

		this.setCounter(counter);
		this.setQuerySize(querySize);
		this.setCancel(cancel);
		this.staticGenesSet = staticGenesSet;
		this.similarity_threshold = similarity_threshold;
		this.method = method; 
		this.querySequences = querySequences;
		this.sequencesWithoutSimilarities = null;
		this.alignmentScoreType = alignmentScoreType;
		
		this.gapsIdentification = false;
	}
	
	///////////////////////////////////
	
	/**
	 * Run the transport similarity searches.
	 * If some threshold parameters were null, this method use the default values.
	 * 
	 * Default values: evalueThreshold(1E-6), bitScoreThreshold(50), queryCoverageThreshold(0.80)
	 * 
	 * @param isTransportersSearch
	 * @param eValueThreshold 
	 * @param bitScoreThreshold
	 * @param queryCoverageThreshold
	 * @return
	 * @throws Exception
	 */
	public ConcurrentLinkedQueue<AlignmentCapsule> runBlastSearch(boolean isTransportersSearch, Double eValueThreshold, Double bitScoreThreshold, Double queryCoverageThreshold) throws Exception {
		
		List<Thread> threads = new ArrayList<Thread>();
//		ConcurrentLinkedQueue<String> queryArray = new ConcurrentLinkedQueue<String>(this.querySequences.keySet());
		int numberOfCores = Runtime.getRuntime().availableProcessors();
		//int numberOfCores = new Double(Runtime.getRuntime().availableProcessors()*1.5).intValue();

		if(this.querySequences.keySet().size()<numberOfCores)
			numberOfCores=this.querySequences.keySet().size();

		this.querySize.set(new Integer(this.querySequences.size()));
		setChanged();
		notifyObservers();
		
		//Distribute querySequences into fastaFiles
		
		logger.debug("Writting query sequences temporary fasta files... ");
		
		List<String> queryFilesPaths = new ArrayList<>();
		List<Map<String,AbstractSequence<?>>> queriesSubSetList = new ArrayList<>();
		
		String path = this.workspaceTaxonomyFolderPath.concat("/queryBlast");
		
		File f = new File (path);
		if(!f.exists())
			f.mkdir();
		
		CreateGenomeFile.buildSubFastaFiles(path, this.querySequences, queriesSubSetList, queryFilesPaths, numberOfCores);
		
//		int batch_size= this.querySequences.size()/numberOfCores;
//		
//		Map<String, AbstractSequence<?>> queriesSubSet = new HashMap<>();
//		List<String> queryFilesPaths = new ArrayList<>();
//		List<Map<String,AbstractSequence<?>>> queriesSubSetList = new ArrayList<>();
//		
//		String path = this.currentTempFolderDirectory.concat("queryBlastSubFasta_");
//		String fastaFileName;
//		
//		int c=0;
//		for (String query : this.querySequences.keySet()) {
//			
//			queriesSubSet.put(query, this.querySequences.get(query));
//
//			if ((c+1)%batch_size==0 && ((c+1)/batch_size < numberOfCores)) {
//				
//				fastaFileName = path.concat(Integer.toString((c+1)/batch_size)).concat("_of_").
//						concat(Integer.toString(numberOfCores)).concat(FileExtensions.PROTEIN_FAA.getExtension());
//				
//				CreateGenomeFile.buildFastaFile(fastaFileName, queriesSubSet);
//				queryFilesPaths.add(fastaFileName);
//				queriesSubSetList.add(queriesSubSet);
//				
//				queriesSubSet = new HashMap<>();
//			}
//			c++;
//		}
//		
//		fastaFileName = path.concat(Integer.toString(numberOfCores)).concat("_of_").
//				concat(Integer.toString(numberOfCores)).concat(FileExtensions.PROTEIN_FAA.getExtension());
//		
//		CreateGenomeFile.buildFastaFile(fastaFileName, queriesSubSet);
//		queriesSubSetList.add(queriesSubSet);
//		queryFilesPaths.add(fastaFileName);
		
		
		ConcurrentLinkedQueue<AlignmentCapsule> alignmentContainerSet = new ConcurrentLinkedQueue<>();
		JAXBContext jc = JAXBContext.newInstance(BlastOutput.class);
		
		for(int i=0; i<numberOfCores; i++) {
			
			ModelAlignments blastAlign;
			
			if(eValueThreshold!=null && bitScoreThreshold!=null && queryCoverageThreshold!=null){
				
				blastAlign = new BlastAlignment(queryFilesPaths.get(i), this.subjectFastaFilePath, queriesSubSetList.get(i), 
						this.similarity_threshold, eValueThreshold, bitScoreThreshold, queryCoverageThreshold, isTransportersSearch, this.cancel, alignmentContainerSet, jc);
			}
			else{
				blastAlign= new BlastAlignment(queryFilesPaths.get(i), this.subjectFastaFilePath, queriesSubSetList.get(i), 
					this.similarity_threshold, isTransportersSearch, this.cancel, alignmentContainerSet, jc);
				
				if(eValueThreshold!=null)
					((BlastAlignment) blastAlign).setEvalueThreshold(eValueThreshold);
				if(bitScoreThreshold!=null)
					((BlastAlignment) blastAlign).setBitScoreThreshold(bitScoreThreshold);
				if(queryCoverageThreshold!=null)
					((BlastAlignment) blastAlign).setQueryCoverageThreshold(queryCoverageThreshold);	
			}
			
			((BlastAlignment) blastAlign).addObserver(this); 

			Thread thread = new Thread(blastAlign);
			threads.add(thread);
			thread.start();
		}
		
		for(Thread thread :threads)
			thread.join();
		
		return alignmentContainerSet;
	}
	///////////////////////////////
	
	
	/**
	 * Run the transport similarity searches.
	 * 
	 * @param allSequences
	 * @throws Exception
	 */
	public ConcurrentLinkedQueue<AlignmentCapsule> runTransportSearch(Map <String, Double> querySpecificThreshold) throws Exception {
		
		List<Thread> threads = new ArrayList<Thread>();
		ConcurrentLinkedQueue<String> queryArray = new ConcurrentLinkedQueue<String>(this.querySequences.keySet());
		int numberOfCores = Runtime.getRuntime().availableProcessors();
		//int numberOfCores = new Double(Runtime.getRuntime().availableProcessors()*1.5).intValue();

		if(this.querySequences.keySet().size()<numberOfCores)
			numberOfCores=this.querySequences.keySet().size();

		this.querySize.set(new Integer(this.querySequences.size()));
		setChanged();
		notifyObservers();

		ConcurrentLinkedQueue<AlignmentCapsule> alignmentContainerSet = new ConcurrentLinkedQueue<>();
		
		for(int i=0; i<numberOfCores; i++) {
			
			Runnable lc	= new PairwiseSequenceAlignement(this.method, this.querySequences, this.staticGenesSet, queryArray, this.similarity_threshold,
				this.counter, this.cancel, AlignmentPurpose.TRANSPORT, this.alignmentScoreType, alignmentContainerSet);
			
			((PairwiseSequenceAlignement) lc).setQuerySpecificThreshold(querySpecificThreshold);
			((PairwiseSequenceAlignement) lc).addObserver(this); 
			Thread thread = new Thread(lc);
			threads.add(thread);
			thread.start();
		}

		for(Thread thread :threads)
			thread.join();
		
		return alignmentContainerSet;
	}


	/**
	 * @throws Exception
	 */
	public ConcurrentLinkedQueue<AlignmentCapsule> run_OrthologGapsSearch(Map<String, List<String>> sequenceIdsSet, ConcurrentLinkedQueue<AlignmentCapsule> alignmentContainerSet) throws Exception {

		boolean recursive = false;
		
		ConcurrentHashMap<String, AbstractSequence<?>> all_sequences = new ConcurrentHashMap<>(querySequences);
		
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if(all_sequences.keySet().size()>0) {

//			this.setAlreadyProcessed(false);
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			List<Thread> threads = new ArrayList<Thread>();
			ConcurrentLinkedQueue<String> queryArray = new ConcurrentLinkedQueue<String>(querySequences.keySet());

			this.querySize.set(new Integer(all_sequences.size()));
			setChanged();
			notifyObservers();

			Map<String, AbstractSequence<?>> ecNumberAnnotations = new HashMap<>();
			ecNumberAnnotations.putAll(this.staticGenesSet);
			
			if(this.sequencesWithoutSimilarities==null) {

				if(this.annotatedGenes!= null && !this.annotatedGenes.isEmpty()) 
					ecNumberAnnotations.keySet().retainAll(this.annotatedGenes);
					
				if(!recursive) {

					this.sequencesWithoutSimilarities = new ConcurrentLinkedQueue<String>();
					this.sequencesWithoutSimilarities.addAll(queryArray);
				}
			}
			else  {
				recursive = true;
				queryArray.retainAll(this.sequencesWithoutSimilarities);
			}
			
//			System.out.println("EC NUBMER ANNOT--->"+ecNumberAnnotations.keySet());
//			System.out.println("EC NUBMER ANNOT--->"+ecNumberAnnotations.size());

			int numberOfCores = Runtime.getRuntime().availableProcessors();

			if(queryArray.size()<numberOfCores)
				numberOfCores=queryArray.size();
			
			if(this.method.equals(Method.Blast)){// && !isGapsIdentification()){
				
				//Distribute querySequences into fastaFiles
				logger.info("Writting query sequences temporary fasta files... ");
				
				List<String> queryFilesPaths = new ArrayList<>();
				List<Map<String,AbstractSequence<?>>> queriesSubSetList = new ArrayList<>();
				
				String path = this.workspaceTaxonomyFolderPath.concat("/queryBlast");
				
				File f = new File (path);
				if(!f.exists())
					f.mkdir();
				
				CreateGenomeFile.buildSubFastaFiles(path, all_sequences, queriesSubSetList, queryFilesPaths, numberOfCores);
				//Distribute querySequences into fastaFiles
				
				//Subject Fasta File
				CreateGenomeFile.buildFastaFile(this.subjectFastaFilePath, ecNumberAnnotations);
				//Subject Fasta File
				
				JAXBContext jc = JAXBContext.newInstance(BlastOutput.class);
				
				logger.info("Starting Blast homology searches... ");
				
				for(int i=0; i<numberOfCores; i++) {
					
					ModelAlignments blastAlign	= new BlastAlignment(queryFilesPaths.get(i), this.subjectFastaFilePath, queriesSubSetList.get(i), 
							this.similarity_threshold, false, this.cancel, alignmentContainerSet, jc);
					
					((BlastAlignment) blastAlign).setEc_number(this.ec_number);
					((BlastAlignment) blastAlign).setModules(this.modules);
					((BlastAlignment) blastAlign).setClosestOrthologs(this.closestOrthologs);
					((BlastAlignment) blastAlign).setSequencesWithoutSimilarities(this.sequencesWithoutSimilarities);
					
					((BlastAlignment) blastAlign).setReferenceTaxonomyScore(this.referenceTaxonomyScore);
					((BlastAlignment) blastAlign).setKegg_taxonomy_scores(this.kegg_taxonomy_scores);
					((BlastAlignment) blastAlign).setReferenceTaxonomyThreshold(this.referenceTaxonomyThreshold);
					((BlastAlignment) blastAlign).setSequenceIdsSet(sequenceIdsSet);
					((BlastAlignment) blastAlign).setBlastPurpose(AlignmentPurpose.ORTHOLOGS);
					
					((BlastAlignment) blastAlign).addObserver(this); 

					Thread thread = new Thread(blastAlign);
					threads.add(thread);
					thread.start();
				}
			}		
			else{
				
				logger.info("Starting pairwise sequence alignements... ");

				for(int i=0; i<numberOfCores; i++) {

					Runnable lc	= new PairwiseSequenceAlignement(method, all_sequences, ecNumberAnnotations, queryArray,
							similarity_threshold, this.counter, this.cancel, AlignmentPurpose.ORTHOLOGS, this.alignmentScoreType, 
							alignmentContainerSet);

					((PairwiseSequenceAlignement) lc).setSequencesWithoutSimilarities(this.sequencesWithoutSimilarities);
					((PairwiseSequenceAlignement) lc).setEc_number(this.ec_number);
					((PairwiseSequenceAlignement) lc).setModules(this.modules);
					((PairwiseSequenceAlignement) lc).setClosestOrthologs(this.closestOrthologs);
					((PairwiseSequenceAlignement) lc).setReferenceTaxonomyScore(this.referenceTaxonomyScore);
					((PairwiseSequenceAlignement) lc).setKegg_taxonomy_scores(this.kegg_taxonomy_scores);
					((PairwiseSequenceAlignement) lc).setReferenceTaxonomyThreshold(this.referenceTaxonomyThreshold);
					((PairwiseSequenceAlignement) lc).setSequenceIdsSet(sequenceIdsSet);

					((PairwiseSequenceAlignement) lc).addObserver(this); 
					Thread thread = new Thread(lc);
					threads.add(thread);
					thread.start();
				}
			}

			for(Thread thread :threads)
				thread.join();
		}
		
		return alignmentContainerSet;
	}
	
	/**
	 * @param sequenceIdsSet
	 * @param alignmentContainerSet
	 * @return
	 * @throws Exception
	 */
	public ConcurrentLinkedQueue<AlignmentCapsule> run_OrthologsSearch(Map<String, List<String>> sequenceIdsSet, ConcurrentLinkedQueue<AlignmentCapsule> alignmentContainerSet) throws Exception {
		
		Boolean recursive = false;
		
		ConcurrentHashMap<String, AbstractSequence<?>> all_sequences = new ConcurrentHashMap<>(querySequences);

		if(all_sequences.keySet().size()>0) {

			if(this.sequencesWithoutSimilarities!=null)
				recursive = true;

			this.run_OrthologGapsSearch(sequenceIdsSet, alignmentContainerSet);


			if(this.compareToFullGenome && !recursive && this.sequencesWithoutSimilarities!=null && !this.sequencesWithoutSimilarities.isEmpty())
				this.run_OrthologsSearch(sequenceIdsSet, alignmentContainerSet);

		}
		return alignmentContainerSet;
	}

	

//	/**
//	 * @return the alreadyProcessed
//	 */
//	public boolean isAlreadyProcessed() {
//		return alreadyProcessed;
//	}
//
//	/**
//	 * @param alreadyProcessed the alreadyProcessed to set
//	 */
//	public void setAlreadyProcessed(boolean alreadyProcessed) {
//		this.alreadyProcessed = alreadyProcessed;
//	}
//
//	/**
//	 * @return the processed
//	 */
//	public boolean isProcessed() {
//		return processed;
//	}
//
//	/**
//	 * @param processed the processed to set
//	 */
//	public void setProcessed(boolean processed) {
//		this.processed = processed;
//	}

	/**
	 * @param counter the counter to set
	 */
	public void setCounter(AtomicInteger counter) {
		this.counter = counter;
	}

	/**
	 * @param querySize the querySize to set
	 */
	public void setQuerySize(AtomicInteger querySize) {
		this.querySize = querySize;
	}

	/**
	 * @return the querySize
	 */
	public AtomicInteger getQuerySize() {
		return querySize;
	}


	/**
	 * @param cancel the cancel to set
	 */
	public void setCancel(AtomicBoolean cancel) {
		this.cancel = cancel;
	}

	/**
	 * @param annotatedGenes
	 */
	public void setAnnotatedGenes(List<String> annotatedGenes) {

		this.annotatedGenes = annotatedGenes;		
	}

	/**
	 * @return
	 */
	public List<String> getAnnotatedGenes() {

		return this.annotatedGenes;		
	}

	/**
	 * @return the ec_number
	 */
	public String getEc_number() {
		return ec_number;
	}

	/**
	 * @param ec_number the ec_number to set
	 */
	public void setEc_number(String ec_number) {
		this.ec_number = ec_number;
	}

	/**
	 * @return the idModule
	 */
	public Map<String, Set<Integer>> getModules() {
		return modules;
	}

	/**
	 * @param genes_ko_modules the idModule to set
	 */
	public void setModules(Map<String, Set<Integer>> modules) {
		this.modules = modules;
	}

	/**
	 * @return the sequencesWithoutSimilarities
	 */
	public ConcurrentLinkedQueue<String> getSequencesWithoutSimilarities() {
		return sequencesWithoutSimilarities;
	}

	/**
	 * @param sequencesWithoutSimilarities the sequencesWithoutSimilarities to set
	 */
	public void setSequencesWithoutSimilarities(
			ConcurrentLinkedQueue<String> sequencesWithoutSimilarities) {
		this.sequencesWithoutSimilarities = sequencesWithoutSimilarities;
	}

	/**
	 * @param revertMapFromSet
	 */
	public void setClosestOrthologs(Map<String, Set<String>> closestOrthologs) {

		this.closestOrthologs = closestOrthologs;
	}

	public void setReferenceTaxonomyScore(int referenceTaxonomyScore) {
		this.referenceTaxonomyScore = referenceTaxonomyScore;

	}

	public void setKegg_taxonomy_scores(Map<String, Integer> kegg_taxonomy_scores) {

		this.kegg_taxonomy_scores = kegg_taxonomy_scores;
	}

	public double getReferenceTaxonomyThreshold() {
		return referenceTaxonomyThreshold;
	}

	public void setReferenceTaxonomyThreshold(double referenceTaxonomyThreshold) {
		this.referenceTaxonomyThreshold = referenceTaxonomyThreshold;
	}
	
	
	/**
	 * @return
	 */
	public String getSubjectFastaFilePath() {
		return this.subjectFastaFilePath;
	}

	/**
	 * @param subjectFastaFilePath
	 */
	public void setSubjectFastaFilePath(String subjectFastaFilePath) {
		this.subjectFastaFilePath = subjectFastaFilePath;
	}

//	/**
//	 * @return the tcdbFastaFilePath
//	 */
//	public String getTcdbFastaFilePath() {
//		return tcdbFastaFilePath;
//	}
//
//	/**
//	 * @param tcdbFastaFilePath the tcdbFastaFilePath to set
//	 */
//	public void setTcdbFastaFilePath(String tcdbFastaFilePath) {
//		this.tcdbFastaFilePath = tcdbFastaFilePath;
//	}

	public boolean isCompareToFullGenome() {
		return compareToFullGenome;
	}

	public void setCompareToFullGenome(boolean compareToFullGenome) {
		this.compareToFullGenome = compareToFullGenome;
	}
	
	/**
	 * @return the gapsIdentification
	 */
	public boolean isGapsIdentification() {
		return gapsIdentification;
	}

	/**
	 * @param gapsIdentification the gapsIdentification to set
	 */
	public void setGapsIdentification(boolean gapsIdentification) {
		this.gapsIdentification = gapsIdentification;
	}

	public String getWorkspaceTaxonomyFolderPath() {
		return workspaceTaxonomyFolderPath;
	}

	public void setWorkspaceTaxonomyFolderPath(String workspaceTaxonomyFolderPath) {
		this.workspaceTaxonomyFolderPath = workspaceTaxonomyFolderPath;
	}

	@Override
	public void update(Observable arg0, Object arg1) {

		setChanged();
		notifyObservers();
	}

}
