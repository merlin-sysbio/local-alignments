package pt.uminho.sysbio.common.local.alignments.core;

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

import org.biojava.nbio.core.sequence.template.AbstractSequence;

import pt.uminho.sysbio.merIin.utilities.capsules.AlignmentCapsule;
import pt.uminho.sysbio.merlin.utilities.Enumerators.AlignmentPurpose;
import pt.uminho.sysbio.merlin.utilities.Enumerators.AlignmentScoreType;
import pt.uminho.sysbio.merlin.utilities.Enumerators.Method;

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
	private Map<String, AbstractSequence<?>> querySequences;
	private List<String> annotatedGenes;
	private ConcurrentLinkedQueue<String> sequencesWithoutSimilarities;
	private String ec_number;
	private Map<String, Set<String>> modules;
	private Map<String, Set<String>> closestOrthologs;
	private int referenceTaxonomyScore;
	private Map<String, Integer> kegg_taxonomy_scores;
	private double referenceTaxonomyThreshold;
	private boolean compareToFullGenome;
	private AlignmentScoreType alignmentScoreType;


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
	public RunSimilaritySearch(Map<String, AbstractSequence<?>> staticGenesSet, double similarity_threshold,  
			Method method, ConcurrentHashMap<String, AbstractSequence<?>> querySequences, 
			AtomicBoolean cancel, AtomicInteger querySize, AtomicInteger counter, 
			AlignmentScoreType alignmentScoreType) throws Exception {

		this.setCounter(counter);
		this.setQuerySize(querySize);
		this.setCancel(cancel);
		this.staticGenesSet = staticGenesSet;
		this.similarity_threshold = similarity_threshold;
		this.method = method; 
		this.querySequences = querySequences;
		this.sequencesWithoutSimilarities = null;
		this.alignmentScoreType = alignmentScoreType;
	}

	/**
	 * Run the transport similarity searches.
	 * 
	 * @param allSequences
	 * @param querySpecificThreshold
	 * @throws Exception
	 */
	public ConcurrentLinkedQueue<AlignmentCapsule> runTransportSearch(ConcurrentHashMap<String, AbstractSequence<?>> allSequences, Map <String, Double> querySpecificThreshold) throws Exception {

			List<Thread> threads = new ArrayList<Thread>();
			ConcurrentLinkedQueue<String> queryArray = new ConcurrentLinkedQueue<String>(allSequences.keySet());
			int numberOfCores = Runtime.getRuntime().availableProcessors();
			//int numberOfCores = new Double(Runtime.getRuntime().availableProcessors()*1.5).intValue();

			if(allSequences.keySet().size()<numberOfCores)
				numberOfCores=allSequences.keySet().size();

			System.out.println("number Of threads: "+numberOfCores);

			this.querySize.set(new Integer(allSequences.size()));
			setChanged();
			notifyObservers();

			ConcurrentLinkedQueue<AlignmentCapsule> alignmentContainerSet = new ConcurrentLinkedQueue<>();
			
			for(int i=0; i<numberOfCores; i++) {

				Runnable lc	= new PairwiseSequenceAlignement(this.method, allSequences, this.staticGenesSet, queryArray, this.similarity_threshold,
					this.counter, this.cancel, AlignmentPurpose.TRANSPORT, this.alignmentScoreType, alignmentContainerSet);

				((PairwiseSequenceAlignement) lc).addObserver(this); 
				Thread thread = new Thread(lc);
				threads.add(thread);
				System.out.println("Start "+i);
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

			int numberOfCores = Runtime.getRuntime().availableProcessors();

			if(queryArray.size()<numberOfCores)
				numberOfCores=queryArray.size();

			System.out.println("number Of threads: "+numberOfCores);
			
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
				System.out.println("Start "+i);
				thread.start();
			}

			for(Thread thread :threads)
				thread.join();

		}
		return alignmentContainerSet;
	}
	
	public ConcurrentLinkedQueue<AlignmentCapsule> run_OrthologsSearch(Map<String, List<String>> sequenceIdsSet, ConcurrentLinkedQueue<AlignmentCapsule> alignmentContainerSet) throws Exception {

		boolean recursive = false;
		
		ConcurrentHashMap<String, AbstractSequence<?>> all_sequences = new ConcurrentHashMap<>(querySequences);
		
		System.out.println("3.1 all_sequences map size ---->" + all_sequences.size());

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

			int numberOfCores = Runtime.getRuntime().availableProcessors();

			if(queryArray.size()<numberOfCores)
				numberOfCores=queryArray.size();

			System.out.println("number Of threads: "+numberOfCores);

			for(int i=0; i<numberOfCores; i++) {

				Runnable lc	= new PairwiseSequenceAlignement(method, all_sequences, ecNumberAnnotations, queryArray,
						similarity_threshold, this.counter, this.cancel, AlignmentPurpose.ORTHOLOGS, this.alignmentScoreType, alignmentContainerSet);
				
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
				System.out.println("Start "+i);
				thread.start();
			}

			for(Thread thread :threads)
				thread.join();

			if(this.compareToFullGenome && !recursive && this.sequencesWithoutSimilarities!=null && !this.sequencesWithoutSimilarities.isEmpty())
				this.run_OrthologsSearch(sequenceIdsSet, alignmentContainerSet);

		}
		System.out.println("3.1 all_sequences map size ---->" + all_sequences.size());
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
	public Map<String, Set<String>> getModules() {
		return modules;
	}

	/**
	 * @param genes_ko_modules the idModule to set
	 */
	public void setModules(Map<String, Set<String>> modules) {
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

	public boolean isCompareToFullGenome() {
		return compareToFullGenome;
	}

	public void setCompareToFullGenome(boolean compareToFullGenome) {
		this.compareToFullGenome = compareToFullGenome;
	}
	
	@Override
	public void update(Observable arg0, Object arg1) {

		setChanged();
		notifyObservers();
	}

}
