package pt.uminho.sysbio.common.local.alignments.core;
/**
 * 
 */

import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Observable;
import java.util.Observer;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;

import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.io.FastaReader;
import org.biojava.nbio.core.sequence.io.GenericFastaHeaderParser;
import org.biojava.nbio.core.sequence.io.ProteinSequenceCreator;
import org.biojava.nbio.core.sequence.template.AbstractSequence;

import pt.uminho.sysbio.common.database.connector.databaseAPI.TransportersAPI;
import pt.uminho.sysbio.common.database.connector.datatypes.Connection;
import pt.uminho.sysbio.common.database.connector.datatypes.DatabaseAccess;
import pt.uminho.sysbio.common.local.alignments.core.Enumerators.AlignmentPurpose;
import pt.uminho.sysbio.common.local.alignments.core.Enumerators.AlignmentScoreType;
import pt.uminho.sysbio.common.local.alignments.core.Enumerators.Method;
import pt.uminho.sysbio.merlin.utilities.DatabaseProgressStatus;

/**
 * @author ODias
 *
 */
public class RunSimilaritySearch extends Observable implements Observer {

	private DatabaseAccess dbAccess;
	private Map<String, AbstractSequence<?>> staticGenesSet;
	private boolean alreadyProcessed, processed;
	private AtomicBoolean cancel;
	private AtomicInteger counter;
	private AtomicInteger querySize;
	private int minimum_number_of_helices;
	private double similarity_threshold;
	private Method method;
	private Map<String, AbstractSequence<?>> querySequences;
	private Set<String> processedGenes;
	private List<String> annotatedGenes;
	private ConcurrentLinkedQueue<String> sequencesWithoutSimilarities;
	private String ec_number;
	private Map<String, Set<String>> modules;
	private int project_id = -1;
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
	 * @param genome_dir
	 * @param project_id
	 * @param alignmentScoreType
	 * @param idLocus
	 * @throws Exception
	 */
	public RunSimilaritySearch(DatabaseAccess dba, Map<String, AbstractSequence<?>> staticGenesSet, int minimum_number_of_helices,
			double similarity_threshold, Method method, File genome_dir, int project_id, AlignmentScoreType alignmentScoreType, Map<String, String> idLocus) throws Exception {

		Map<String, AbstractSequence<?>> querySequences = new HashMap<>();
		if(genome_dir.isDirectory()) {

			for(File genome_file:genome_dir.listFiles()) {

				if(genome_file.isFile()) {

					FastaReader<ProteinSequence, AminoAcidCompound> fastaReader = new FastaReader<ProteinSequence,AminoAcidCompound>(genome_file, 
							new GenericFastaHeaderParser<ProteinSequence,AminoAcidCompound>(), 
							new ProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet()));
					
					Map<String, AbstractSequence<?>> genome_map  =  new HashMap<>();
					genome_map.putAll(fastaReader.process());

					querySequences.putAll(genome_map);
				}
			}
		}

		this.setCounter(new AtomicInteger(0));
		this.setQuerySize(new AtomicInteger(0));
		this.setCancel(new AtomicBoolean(false));
		this.dbAccess = dba;
		this.staticGenesSet = staticGenesSet;
		this.minimum_number_of_helices = minimum_number_of_helices;
		this.similarity_threshold = similarity_threshold;
		this.method = method; 
		this.querySequences = querySequences;
		this.sequencesWithoutSimilarities = null;
		this.project_id = project_id;
		this.alignmentScoreType = alignmentScoreType;
	}

	/**
	 * Run similarity searches constructor.
	 * 
	 * @param dbAccess
	 * @param staticGenesSet
	 * @param minimum_number_of_helices
	 * @param similarity_threshold
	 * @param method
	 * @param genome
	 * @param project_id
	 * @param alignmentScoreType
	 * @param idLocus
	 * @throws Exception
	 */
	public RunSimilaritySearch(DatabaseAccess dba, Map<String, AbstractSequence<?>> staticGenesSet, int minimum_number_of_helices, double similarity_threshold, 
			Method method, Map<String, AbstractSequence<?>> genome,
			int project_id, AlignmentScoreType alignmentScoreType, Map<String, String> idLocus) throws Exception {

		this.setCounter(new AtomicInteger(0));
		this.setQuerySize(new AtomicInteger(0));
		this.setCancel(new AtomicBoolean(false));
		this.dbAccess = dba;
		this.staticGenesSet = staticGenesSet;
		this.minimum_number_of_helices = minimum_number_of_helices;
		this.similarity_threshold = similarity_threshold;
		this.method = method; 
		this.querySequences = genome;
		this.sequencesWithoutSimilarities = null;
		this.project_id = project_id;
		this.alignmentScoreType = alignmentScoreType;
	}

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
	public RunSimilaritySearch(DatabaseAccess dba, Map<String, AbstractSequence<?>> staticGenesSet, int minimum_number_of_helices, double similarity_threshold,  
			Method method, Map<String, AbstractSequence<?>> querySequences, 
			AtomicBoolean cancel, AtomicInteger querySize, AtomicInteger counter, 
			int project_id, AlignmentScoreType alignmentScoreType) throws Exception {

		this.setCounter(counter);
		this.setQuerySize(querySize);
		this.setCancel(cancel);
		this.dbAccess = dba;
		this.staticGenesSet = staticGenesSet;
		this.minimum_number_of_helices = minimum_number_of_helices;
		this.similarity_threshold = similarity_threshold;
		this.method = method; 
		this.querySequences = querySequences;
		this.sequencesWithoutSimilarities = null;
		this.project_id = project_id;
		this.alignmentScoreType = alignmentScoreType;
	}

	/**
	 * Run similarity searches constructor.
	 *  
	 * @param dbAccess
	 * @param staticGenesSet
	 * @param similarity_threshold
	 * @param method
	 * @param orthologs
	 * @param cancel
	 * @param querySize
	 * @param counter
	 * @param d
	 * @param alignment
	 */
	public RunSimilaritySearch(DatabaseAccess dbAccess, Map<String, AbstractSequence<?>> staticGenesSet, double similarity_threshold, Method method, 
			ConcurrentHashMap<String, AbstractSequence<?>> orthologs, 
			AtomicBoolean cancel, AtomicInteger querySize, AtomicInteger counter, 
			int project_id, AlignmentScoreType alignmentScoreType) throws Exception {
	
		this.setCounter(counter);
		this.setQuerySize(querySize);
		this.setCancel(cancel);
		this.dbAccess = dbAccess;
		this.staticGenesSet = staticGenesSet;
		this.similarity_threshold = similarity_threshold;
		this.method = method; 
		this.querySequences = orthologs;
		this.sequencesWithoutSimilarities = null;
		this.project_id = project_id;
		this.alignmentScoreType = alignmentScoreType;
	}
	

	/**
	 * Run the transport similarity searches.
	 * 
	 * @param transmembraneGenes 
	 * @param allSequences 
	 * @throws Exception
	 */
	public void runTransportSearch(Map<String, Integer> transmembraneGenes, ConcurrentHashMap<String, AbstractSequence<?>> allSequences) throws Exception {

		Connection conn = new Connection(dbAccess);

		this.processedGenes = TransportersAPI.retrieveProcessedTransportAlignmentGenes(conn);

		this.setProcessed(false);
		ConcurrentHashMap<String, String> locus_ids = new ConcurrentHashMap<String, String>();

		for(String sequence_id : new HashSet<String>(allSequences.keySet())) {

			String status = null;

			if(transmembraneGenes.get(sequence_id) >= minimum_number_of_helices && !this.processedGenes.contains(sequence_id)) {

				int seqLength = allSequences.get(sequence_id).getLength();

				String matrix;
				if(seqLength<35)
					matrix="pam30";
				else if(seqLength<50)
					matrix="pam70";
				else if(seqLength<85)
					matrix="blosum80";
				else
					matrix="blosum62";

				status = DatabaseProgressStatus.PROCESSING.toString();
				TransportersAPI.loadTransportAlignmentsGenes(sequence_id, matrix, transmembraneGenes.get(sequence_id), conn, locus_ids, status, this.project_id);
				processedGenes.add(sequence_id);
			}
			else {

				status = DatabaseProgressStatus.PROCESSED.toString();
				TransportersAPI.loadTransportAlignmentsGenes(sequence_id, null, transmembraneGenes.get(sequence_id),conn, locus_ids, status, this.project_id);
				allSequences.remove(sequence_id);
				transmembraneGenes.remove(sequence_id);
			}
		}
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if(transmembraneGenes.keySet().size()>0) {

			this.setAlreadyProcessed(false);
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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

			for(int i=0; i<numberOfCores; i++) {

				Runnable lc	= new PairwiseSequenceAlignement(this.method, allSequences, this.staticGenesSet, queryArray, this.dbAccess, this.similarity_threshold, 
						transmembraneGenes, locus_ids, this.counter, this.cancel, AlignmentPurpose.TRANSPORT, this.alignmentScoreType);

				((PairwiseSequenceAlignement) lc).addObserver(this); 
				Thread thread = new Thread(lc);
				threads.add(thread);
				System.out.println("Start "+i);
				thread.start();
			}

			for(Thread thread :threads)
				thread.join();
		}
		else {

			this.setAlreadyProcessed(true);
		}
		this.setProcessed(true);
		conn.closeConnection();
	}


	/**
	 * @throws Exception
	 */
	public void run_OrthologsSearch() throws Exception {

		boolean recursive = false;

		ConcurrentHashMap<String, AbstractSequence<?>> all_sequences = new ConcurrentHashMap<>(querySequences);

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if(all_sequences.keySet().size()>0) {

			this.setAlreadyProcessed(false);
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

				Runnable lc	= new PairwiseSequenceAlignement(method, all_sequences, ecNumberAnnotations, queryArray,  dbAccess, 
						similarity_threshold, null, null, this.counter, this.cancel, AlignmentPurpose.ORTHOLOGS, this.alignmentScoreType);
				
				((PairwiseSequenceAlignement) lc).setSequencesWithoutSimilarities(this.sequencesWithoutSimilarities);
				((PairwiseSequenceAlignement) lc).setEc_number(this.ec_number);
				((PairwiseSequenceAlignement) lc).setModules(this.modules);
				((PairwiseSequenceAlignement) lc).setClosestOrthologs(this.closestOrthologs);
				((PairwiseSequenceAlignement) lc).setReferenceTaxonomyScore(this.referenceTaxonomyScore);
				((PairwiseSequenceAlignement) lc).setKegg_taxonomy_scores(this.kegg_taxonomy_scores);
				((PairwiseSequenceAlignement) lc).setReferenceTaxonomyThreshold(this.referenceTaxonomyThreshold);

				((PairwiseSequenceAlignement) lc).addObserver(this); 
				Thread thread = new Thread(lc);
				threads.add(thread);
				System.out.println("Start "+i);
				thread.start();
			}

			for(Thread thread :threads)
				thread.join();

			if(this.compareToFullGenome && !recursive && this.sequencesWithoutSimilarities!=null && !this.sequencesWithoutSimilarities.isEmpty())
				this.run_OrthologsSearch();

		}
		else {

			this.setAlreadyProcessed(true);
		}
		this.setProcessed(true);
	}

	/**
	 * 
	 * 
	 * @param url
	 * @return
	 * @throws Exception
	 */
	public static Map<String, AbstractSequence<?>> setTCDB(String url) throws Exception {

		InputStream tcdbInputStream = (new URL(url)).openStream();
		BufferedReader br= new BufferedReader(new InputStreamReader(tcdbInputStream));
		StringBuilder sb = new StringBuilder();
		String line;

		while ((line = br.readLine()) != null)
			sb.append(line+"\n");

			String theString = sb.toString().replace("</p>", "").replace("<p>", "").replace(">gnl|TC-DB|xxxxxx 3.A.1.205.14 \ndsfgdfg", "");
		byte[] bytes = theString.getBytes("utf-8");
		tcdbInputStream =  new ByteArrayInputStream(bytes);

		FastaReader<ProteinSequence,AminoAcidCompound> fastaReader = new FastaReader<ProteinSequence,AminoAcidCompound>(
				tcdbInputStream, 
				//tcdbFile,
				new GenericFastaHeaderParser<ProteinSequence,AminoAcidCompound>(), 
				new ProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet()));
		
		Map<String, AbstractSequence<?>> tcdb  =  new HashMap<>();
		 tcdb.putAll(fastaReader.process());
		return tcdb;
	}

	/**
	 * @return the alreadyProcessed
	 */
	public boolean isAlreadyProcessed() {
		return alreadyProcessed;
	}

	/**
	 * @param alreadyProcessed the alreadyProcessed to set
	 */
	public void setAlreadyProcessed(boolean alreadyProcessed) {
		this.alreadyProcessed = alreadyProcessed;
	}

	/**
	 * @return the processed
	 */
	public boolean isProcessed() {
		return processed;
	}

	/**
	 * @param processed the processed to set
	 */
	public void setProcessed(boolean processed) {
		this.processed = processed;
	}

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

	public Set<String> getProcessedGenes() {
		return processedGenes;
	}

	public void setProcessedGenes(Set<String> processedGenes) {
		this.processedGenes = processedGenes;
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
