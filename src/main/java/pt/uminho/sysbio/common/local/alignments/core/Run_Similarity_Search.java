package pt.uminho.sysbio.common.local.alignments.core;
/**
 * 
 */

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
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

import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava3.core.sequence.io.FastaReader;
import org.biojava3.core.sequence.io.GenericFastaHeaderParser;
import org.biojava3.core.sequence.io.ProteinSequenceCreator;

import pt.uminho.sysbio.common.bioapis.externalAPI.ncbi.NcbiEFetchSequenceStub_API;
import pt.uminho.sysbio.common.bioapis.externalAPI.ncbi.NcbiAPI;
import pt.uminho.sysbio.common.database.connector.datatypes.MySQLMultiThread;
import pt.uminho.sysbio.merlin.utilities.DatabaseProgressStatus;

/**
 * @author ODias
 *
 */
public class Run_Similarity_Search extends Observable implements Observer {

	private MySQLMultiThread msqlmt;
	private Map<String, ProteinSequence> staticGenesSet;
	private boolean alreadyProcessed, processed;
	private AtomicBoolean cancel;
	private AtomicInteger counter;
	private AtomicInteger querySize;
	private List<File> tmhmmFiles;
	private int minimum_number_of_helices;
	private double similarity_threshold;
	private Method method;
	private boolean isNCBI;
	private Map<String, ProteinSequence> querySequences;
	private boolean isUseProxy, isUseAuthentication;
	private String host, port, user, pass;
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


	/**
	 * @param msqlmt
	 * @param staticGenesSet
	 * @param tmhmm_file_dir
	 * @param minimum_number_of_helices
	 * @param similarity_threshold
	 * @param method
	 * @param isNCBI
	 * @param genome_dir
	 * @param project_id
	 * @throws Exception
	 */
	public Run_Similarity_Search(MySQLMultiThread msqlmt, Map<String, ProteinSequence> staticGenesSet, File tmhmm_file_dir, int minimum_number_of_helices,
			double similarity_threshold, Method method, boolean isNCBI, File genome_dir, int project_id) throws Exception {

		List<File> tmhmmFiles = new ArrayList<File>();
		if(tmhmm_file_dir.isDirectory()) {

			for(File tmhmm_file:tmhmm_file_dir.listFiles()) {

				if(tmhmm_file.isFile()) {

					tmhmmFiles.add(tmhmm_file);
				}
			}
		}

		Map<String, ProteinSequence> querySequences = new HashMap<String, ProteinSequence>();
		if(genome_dir.isDirectory()) {

			for(File genome_file:genome_dir.listFiles()) {

				if(genome_file.isFile()) {

					FastaReader<ProteinSequence, AminoAcidCompound> fastaReader = new FastaReader<ProteinSequence,AminoAcidCompound>(genome_file, 
							new GenericFastaHeaderParser<ProteinSequence,AminoAcidCompound>(), 
							new ProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet()));
					Map<String, ProteinSequence> genome_map  = fastaReader.process();

					querySequences.putAll(genome_map);
				}
			}
		}

		this.setCounter(new AtomicInteger(0));
		this.setQuerySize(new AtomicInteger(0));
		this.setCancel(new AtomicBoolean(false));
		this.msqlmt = msqlmt;
		this.staticGenesSet = staticGenesSet;
		this.tmhmmFiles = tmhmmFiles;
		this.minimum_number_of_helices = minimum_number_of_helices;
		this.similarity_threshold = similarity_threshold;
		this.method = method; 
		this.isNCBI = isNCBI;
		this.querySequences = querySequences;
		this.isUseProxy = false;
		this.isUseAuthentication = false;
		this.sequencesWithoutSimilarities = null;
		this.project_id = project_id;
	}

	/**
	 * @param msqlmt
	 * @param staticGenesSet
	 * @param tmhmm_file_dir
	 * @param minimum_number_of_helices
	 * @param similarity_threshold
	 * @param method
	 * @param isNCBI
	 * @param genome
	 * @param project_id
	 * @throws Exception
	 */
	public Run_Similarity_Search(MySQLMultiThread msqlmt, Map<String, ProteinSequence> staticGenesSet, File tmhmm_file_dir, int minimum_number_of_helices,
			double similarity_threshold, Method method, boolean isNCBI, Map<String, ProteinSequence> genome, int project_id) throws Exception {

		List<File> tmhmmFiles = new ArrayList<File>();
		if(tmhmm_file_dir.isDirectory()) {

			for(File tmhmm_file:tmhmm_file_dir.listFiles()) {

				if(tmhmm_file.isFile()) {

					tmhmmFiles.add(tmhmm_file);
				}
			}
		}

		this.setCounter(new AtomicInteger(0));
		this.setQuerySize(new AtomicInteger(0));
		this.setCancel(new AtomicBoolean(false));
		this.msqlmt = msqlmt;
		this.staticGenesSet = staticGenesSet;
		this.tmhmmFiles = tmhmmFiles;
		this.minimum_number_of_helices = minimum_number_of_helices;
		this.similarity_threshold = similarity_threshold;
		this.method = method; 
		this.isNCBI = isNCBI;
		this.querySequences = genome;
		this.isUseProxy = false;
		this.isUseAuthentication = false;
		this.sequencesWithoutSimilarities = null;
		this.project_id = project_id;
	}

	/**
	 * @param msqlmt
	 * @param staticGenesSet
	 * @param tmhmmFiles
	 * @param minimum_number_of_helices
	 * @param similarity_threshold
	 * @param method
	 * @param isNCBI
	 * @param querySequences
	 * @param cancel
	 * @param querySize
	 * @param counter
	 * @param project_id
	 * @throws Exception
	 */
	public Run_Similarity_Search(MySQLMultiThread msqlmt, Map<String, ProteinSequence> staticGenesSet, List<File> tmhmmFiles, int minimum_number_of_helices,
			double similarity_threshold, Method method, boolean isNCBI, Map<String, ProteinSequence> querySequences, AtomicBoolean cancel, 
			AtomicInteger querySize, AtomicInteger counter, int project_id) throws Exception {

		this.setCounter(counter);
		this.setQuerySize(querySize);
		this.setCancel(cancel);
		this.msqlmt = msqlmt;
		this.staticGenesSet = staticGenesSet;
		this.tmhmmFiles = tmhmmFiles;
		this.minimum_number_of_helices = minimum_number_of_helices;
		this.similarity_threshold = similarity_threshold;
		this.method = method; 
		this.isNCBI = isNCBI;
		this.querySequences = querySequences;
		this.isUseProxy = false;
		this.isUseAuthentication = false;
		this.sequencesWithoutSimilarities = null;
		this.project_id = project_id;
	}

	/**
	 * @param msqlmt
	 * @param staticGenesSet
	 * @param similarity_threshold
	 * @param method
	 * @param querySequences
	 * @param cancel
	 * @param querySize
	 * @param counter
	 * @param project_id
	 */
	public Run_Similarity_Search(MySQLMultiThread msqlmt, Map<String, ProteinSequence> staticGenesSet, double similarity_threshold,
			Method method, ConcurrentHashMap<String, ProteinSequence> querySequences, AtomicBoolean cancel, AtomicInteger querySize, AtomicInteger counter, int project_id) {

		this.setCounter(counter);
		this.setQuerySize(querySize);
		this.setCancel(cancel);
		this.msqlmt = msqlmt;
		this.staticGenesSet = staticGenesSet;
		this.tmhmmFiles = null;
		this.minimum_number_of_helices = -1;
		this.similarity_threshold = similarity_threshold;
		this.method = method; 
		this.querySequences = querySequences;
		this.isUseProxy = false;
		this.isUseAuthentication = false;
		this.sequencesWithoutSimilarities = null;
		this.project_id = project_id;
	}

	/**
	 * @throws Exception
	 */
	public void run_TransportSearch() throws Exception {

		this.retrieveProcessedGenes();
		Connection conn = msqlmt.openConnection();
		this.setProcessed(false);
		ConcurrentHashMap<String, String> locus_ids = new ConcurrentHashMap<String, String>();
		Map<String, Integer> tmhmm_genes = new ConcurrentHashMap<String, Integer>();
		Map<String, ProteinSequence> new_sequences = new HashMap<String, ProteinSequence>();
		ConcurrentHashMap<String, ProteinSequence> all_sequences = new ConcurrentHashMap<String, ProteinSequence>();
		Map<String, Integer> new_tmhmm_genes = new HashMap<String, Integer>();

		NcbiEFetchSequenceStub_API fetchStub = null;

		if(isNCBI) {

			fetchStub = new NcbiEFetchSequenceStub_API(50);
		}

		for(File tmhmm_file:tmhmmFiles) {

			if(tmhmm_file.isFile()) {

				if(isNCBI) {

					tmhmm_genes.putAll(NcbiAPI.readTMHMMGenbankNCBI(tmhmm_file, 0));
					Map<String, String> idLocus = fetchStub.getLocusFromID(tmhmm_genes.keySet(),100);

					for (String id : idLocus.keySet()) {

						new_tmhmm_genes.put(idLocus.get(id), tmhmm_genes.get(id));

						if(querySequences.containsKey(id)) {

							new_sequences.put(idLocus.get(id), querySequences.get(id));
						}
						else {

							throw new IOException("Wrong TMHMM file");
						}
					}
					tmhmm_genes = new_tmhmm_genes;
					all_sequences = new ConcurrentHashMap<String, ProteinSequence>(new_sequences);
				}
				else {

					tmhmm_genes.putAll(NcbiAPI.readTMHMMGenbank(tmhmm_file, 0));

					for (String id : querySequences.keySet()) {

						if(tmhmm_genes.containsKey(id)) {

							all_sequences.put(id,querySequences.get(id));
						}
					}
				}
			}
		}

		if(tmhmm_genes.size()==0) {

			throw new IOException ("Verify tmhmm files path!");
		}

		for(String locus_tag:new HashSet<String>(all_sequences.keySet())) {

			if(tmhmm_genes.get(locus_tag) >= minimum_number_of_helices && !this.processedGenes.contains(locus_tag)) {

				int seqLength = all_sequences.get(locus_tag).getLength();

				String matrix;
				if(seqLength<35) { matrix="pam30";}
				else if(seqLength<50) { matrix="pam70";}
				else if(seqLength<85) { matrix="blosum80";}
				else { matrix="blosum62";}
				this.load_locus_tag(locus_tag, matrix, tmhmm_genes.get(locus_tag),conn, locus_ids);
				processedGenes.add(locus_tag);
			}
			else {

				this.load_locus_tag(locus_tag, null, tmhmm_genes.get(locus_tag),conn, locus_ids);
				all_sequences.remove(locus_tag);
				tmhmm_genes.remove(locus_tag);
			}
		}
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if(tmhmm_genes.keySet().size()>0) {

			this.setAlreadyProcessed(false);
			/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			List<Thread> threads = new ArrayList<Thread>();
			ConcurrentLinkedQueue<String> queryArray = new ConcurrentLinkedQueue<String>(all_sequences.keySet());
			int numberOfCores = Runtime.getRuntime().availableProcessors();
			//int numberOfCores = new Double(Runtime.getRuntime().availableProcessors()*1.5).intValue();

			if(all_sequences.keySet().size()<numberOfCores) {

				numberOfCores=all_sequences.keySet().size();
			}
			System.out.println("number Of threads: "+numberOfCores);

			this.querySize.set(new Integer(all_sequences.size()));
			setChanged();
			notifyObservers();

			for(int i=0; i<numberOfCores; i++) {

				Runnable lc	= new PairwiseSequenceAlignement(method, all_sequences, this.staticGenesSet, queryArray, msqlmt, 
						similarity_threshold,tmhmm_genes, locus_ids, this.counter, this.cancel, AlignmentPurpose.Transport);

				((PairwiseSequenceAlignement) lc).addObserver(this); 
				Thread thread = new Thread(lc);
				threads.add(thread);
				System.out.println("Start "+i);
				thread.start();
			}

			for(Thread thread :threads) {

				thread.join();
			}
		}
		else {

			this.setAlreadyProcessed(true);
		}
		this.setProcessed(true);
		conn.close();
	}


	/**
	 * @throws Exception
	 */
	public void run_OrthologsSearch() throws Exception {

		Connection conn = msqlmt.openConnection();

		boolean recursive = false;

		ConcurrentHashMap<String, ProteinSequence> all_sequences = new ConcurrentHashMap<String, ProteinSequence>(querySequences);

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

			Map<String, ProteinSequence> ec_number_annotations = new HashMap<String, ProteinSequence>();
			ec_number_annotations.putAll(this.staticGenesSet);

			if(this.sequencesWithoutSimilarities==null) {

				if(this.annotatedGenes!= null && !this.annotatedGenes.isEmpty()) 					
					ec_number_annotations.keySet().retainAll(this.annotatedGenes);

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

			if(queryArray.size()<numberOfCores) {

				numberOfCores=queryArray.size();
			}
			System.out.println("number Of threads: "+numberOfCores);

			for(int i=0; i<numberOfCores; i++) {

				Runnable lc	= new PairwiseSequenceAlignement(method, all_sequences, ec_number_annotations, queryArray,  msqlmt, 
						similarity_threshold, null, null, this.counter, this.cancel, AlignmentPurpose.Orthologs);
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

			for(Thread thread :threads) {

				thread.join();
			}

			if(this.compareToFullGenome && !recursive && this.sequencesWithoutSimilarities!=null && !this.sequencesWithoutSimilarities.isEmpty())
				this.run_OrthologsSearch();

		}
		else {

			this.setAlreadyProcessed(true);
		}
		this.setProcessed(true);
		conn.close();

	}

	/**
	 * @param url
	 * @return
	 * @throws Exception
	 */
	public static Map<String, ProteinSequence> set_TCDB(String url) throws Exception {

		InputStream tcdbInputStream = (new URL(url)).openStream();
		BufferedReader br= new BufferedReader(new InputStreamReader(tcdbInputStream));
		StringBuilder sb = new StringBuilder();
		String line;

		while ((line = br.readLine()) != null) {

			sb.append(line+"\n");
		} 
		String theString = sb.toString().replace("</p>", "").replace("<p>", "").replace(">gnl|TC-DB|xxxxxx 3.A.1.205.14 \ndsfgdfg", "");
		byte[] bytes = theString.getBytes("utf-8");
		tcdbInputStream =  new ByteArrayInputStream(bytes);

		//File tcdbFile = new File("C:/tcdb.txt");
		FastaReader<ProteinSequence,AminoAcidCompound> fastaReader
		= new FastaReader<ProteinSequence,AminoAcidCompound>(
				tcdbInputStream, 
				//tcdbFile,
				new GenericFastaHeaderParser<ProteinSequence,AminoAcidCompound>(), 
				new ProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet()));
		Map<String, ProteinSequence> tcdb  = fastaReader.process();
		return tcdb;
	}

	/**
	 * @param queryMap 
	 * @throws SQLException 
	 */
	public void removeLoadedCandidates(Set<String> candidates, Map<String, ProteinSequence> queryMap) throws SQLException {

		Connection conn = this.msqlmt.openConnection();
		Statement stmt;

		try {

			stmt = conn.createStatement();
			ResultSet rs = stmt.executeQuery("SELECT DISTINCT(locus_tag) FROM sw_reports"); //LEFT JOIN sw_similarities ON sw_reports_id=sw_reports.id WHERE similarity>"+threshold);

			while(rs.next()) {

				candidates.remove(rs.getString(1));
				queryMap.remove(rs.getString(1));
			}
			stmt.close();
			this.msqlmt.closeConnection(conn);
		}
		catch (SQLException e) {e.printStackTrace();}
	}

	/**
	 * @param msmt 
	 * @param dataSource
	 * @throws IOException 
	 * @throws SQLException 
	 */
	public static void makeFile(String output_file_name, MySQLMultiThread msqlmt) throws IOException, SQLException {

		Connection conn = msqlmt.openConnection();
		FileWriter fstream = new FileWriter(output_file_name);
		BufferedWriter out = new BufferedWriter(fstream);

		Statement stmt = conn.createStatement();

		ResultSet rs = stmt.executeQuery("SELECT * FROM sw_reports "
				+"LEFT JOIN sw_similarities ON sw_reports.id=sw_similarities.sw_reports_id "
				+"LEFT JOIN sw_hits ON sw_hits.id=sw_similarities.sw_hits_id "
				+"ORDER BY sw_reports.locus_tag, similarity DESC ");

		out.write("locus tag\tsimilarity\thomologue ID\tTCDB ID\tnumber of helices\n");
		String locus="";
		while(rs.next()) {

			if(!locus.equals(rs.getString(1)) && rs.getString(8)!=null) {

				locus=rs.getString(1);
			}

			if(rs.getString(8)!=null) {

				out.write(rs.getString(2)+"\t"+rs.getString(8)+"\t"+rs.getString(10)+"\t"+rs.getString(11)+"\t"+rs.getString(5)+"\n");
			}
		}
		//Close the output stream
		out.close();
		msqlmt.closeConnection(conn);
	}	


	/**
	 * @author ODias
	 * available methods for alignment
	 */
	public static enum Method {

		SmithWaterman,
		NeedlemanWunsch
	}

	/**
	 * @author ODias
	 * available methods for alignment
	 */
	public static enum AlignmentPurpose {

		Transport,
		Orthologs
	}

	/**
	 * @param locus_tag
	 * @param matrix
	 * @param locus_ids 
	 * @return
	 * @throws SQLException
	 */
	private String load_locus_tag(String locus_tag, String matrix, int tmd, Connection conn, ConcurrentHashMap<String,String> locus_ids) throws SQLException {

		String result = null;
		if(locus_ids.contains(locus_tag)) {

			result=locus_ids.get(locus_tag);
		}
		else {

			String status = null;
			if(tmd<minimum_number_of_helices) {

				status = DatabaseProgressStatus.PROCESSED.toString();
			}
			else {

				status = DatabaseProgressStatus.PROCESSING.toString();
			}


			java.sql.Date sqlToday = new java.sql.Date((new java.util.Date()).getTime());
			Statement stmt = conn.createStatement();
			ResultSet rs = stmt.executeQuery("SELECT id FROM sw_reports WHERE locus_tag='"+locus_tag+"'");

			if(!rs.next()) {

				stmt.execute("INSERT INTO sw_reports (locus_tag, date, matrix, number_TMD, project_id, status) " +
						"VALUES ('"+locus_tag+"','"+sqlToday+"','"+matrix+"','"+tmd+"',"+this.project_id+",'"+status+"')");
				rs = stmt.executeQuery("SELECT LAST_INSERT_ID()");
				rs.next();
			}
			result = rs.getString(1);
			rs.close();
			stmt=null;
			locus_ids.put(locus_tag,result);
		}
		return result;
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
	 * @return the isUseProxy
	 */
	public boolean isUseProxy() {
		return isUseProxy;
	}

	/**
	 * @param isUseProxy the isUseProxy to set
	 */
	public void setUseProxy(boolean isUseProxy) {
		this.isUseProxy = isUseProxy;
	}

	/**
	 * @return the isUseAuthentication
	 */
	public boolean isUseAuthentication() {
		return isUseAuthentication;
	}

	/**
	 * @param isUseAuthentication the isUseAuthentication to set
	 */
	public void setUseAuthentication(boolean isUseAuthentication) {
		this.isUseAuthentication = isUseAuthentication;
	}

	/**
	 * @return the host
	 */
	public String getHost() {
		return host;
	}

	/**
	 * @param host the host to set
	 */
	public void setHost(String host) {
		this.host = host;
	}

	/**
	 * @return the port
	 */
	public String getPort() {
		return port;
	}

	/**
	 * @param port the port to set
	 */
	public void setPort(String port) {
		this.port = port;
	}

	/**
	 * @return the user
	 */
	public String getUser() {
		return user;
	}

	/**
	 * @param user the user to set
	 */
	public void setUser(String user) {
		this.user = user;
	}

	/**
	 * @return the pass
	 */
	public String getPass() {
		return pass;
	}

	/**
	 * @param pass the pass to set
	 */
	public void setPass(String pass) {
		this.pass = pass;
	}

	/**
	 * @param cancel the cancel to set
	 */
	public void setCancel(AtomicBoolean cancel) {
		this.cancel = cancel;
	}

	@Override
	public void update(Observable arg0, Object arg1) {

		setChanged();
		notifyObservers();
	}

	public Set<String> getProcessedGenes() {
		return processedGenes;
	}

	public void setProcessedGenes(Set<String> processedGenes) {
		this.processedGenes = processedGenes;
	}

	//	/**
	//	 * @param args
	//	 * @throws Exception 
	//	 * @throws IOException 
	//	 * @throws MalformedURLException 
	//	 */
	//	public static void main(String[] args) throws MalformedURLException, IOException, Exception {
	//
	//		MySQLMultiThread msqlmt = new MySQLMultiThread("localhost","3306","test_transporters","root","");
	//		boolean isncbi=true;
	//		Run_Similarity_Search run_Similarity_Search = new Run_Similarity_Search(
	//				msqlmt,
	//				Run_Similarity_Search.set_TCDB("http://www.tcdb.org/public/tcdb"), 
	//				new File("../transport_systems/test/tmhmm"), //dir tmhmm
	//				1,
	//				0.10,
	//				Method.SmithWaterman,
	//				isncbi,
	//				new File("../transport_systems/test") // dir genome
	//				);
	//
	//		run_Similarity_Search.run_TransportSearch();
	//
	//		String similaritiesFile = "../transport_systems/test/tmhmm/simfiletest"; // ficheiro necessário para próxima etapa
	//		Run_Similarity_Search.makeFile(similaritiesFile,msqlmt);
	//	}

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

	/**
	 * @throws SQLException
	 */
	public void retrieveProcessedGenes() throws SQLException{

		this.processedGenes = new HashSet<String>();

		Connection conn = this.msqlmt.openConnection();
		Statement statement = conn.createStatement();

		ResultSet rs = statement.executeQuery("SELECT locus_tag FROM genes WHERE status <> 'PROCESSING'");

		while(rs.next()) {

			processedGenes.add(rs.getString(1));
		}
		statement.close();
		this.msqlmt.closeConnection(conn);
	}

}
