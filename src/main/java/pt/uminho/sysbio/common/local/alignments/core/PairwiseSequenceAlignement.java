package pt.uminho.sysbio.common.local.alignments.core;

import java.io.InputStreamReader;
import java.io.Reader;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Observable;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;

import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava.nbio.alignment.Alignments.PairwiseSequenceScorerType;
import org.biojava.nbio.alignment.NeedlemanWunsch;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.alignment.SmithWaterman;
import org.biojava.nbio.alignment.template.AbstractPairwiseSequenceAligner;
import org.biojava.nbio.alignment.template.GapPenalty;
import org.biojava.nbio.alignment.template.PairwiseSequenceAligner;
import org.biojava.nbio.core.alignment.matrices.SimpleSubstitutionMatrix;
import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.template.CompoundSet;

import pt.uminho.ceb.biosystems.mew.utilities.datastructures.pair.Pair;
import pt.uminho.sysbio.common.database.connector.databaseAPI.HomologyAPI;
import pt.uminho.sysbio.common.database.connector.databaseAPI.ModelAPI;
import pt.uminho.sysbio.common.database.connector.databaseAPI.TransportersAPI;
import pt.uminho.sysbio.common.database.connector.datatypes.Connection;
import pt.uminho.sysbio.common.database.connector.datatypes.DatabaseAccess;
import pt.uminho.sysbio.common.database.connector.datatypes.Enumerators.DatabaseType;
import pt.uminho.sysbio.common.local.alignments.core.Enumerators.AlignmentPurpose;
import pt.uminho.sysbio.common.local.alignments.core.Enumerators.AlignmentScoreType;
import pt.uminho.sysbio.common.local.alignments.core.Enumerators.Matrix;
import pt.uminho.sysbio.common.local.alignments.core.Enumerators.Method;
import pt.uminho.sysbio.common.local.alignments.core.datatype.AlignmentContainer;
import pt.uminho.sysbio.merlin.utilities.DatabaseProgressStatus;

/**
 * @author ODias
 *
 */
public class PairwiseSequenceAlignement extends Observable implements Runnable {

	private static int _MIN_RESIDUES = 100;
	private static double _SCORE_INCREMENT = 0.2, _LENGTH_DIVISION = 2;
	private ConcurrentLinkedQueue<String> queryArray;
	private Map<String, ProteinSequence> staticSubjectMap; 
	private Map<String,Integer> numberOfHelicesMap;
	private ConcurrentHashMap<String, ProteinSequence> concurrentQueryMap;
	private DatabaseAccess databaseAccess;
	private Method method;
	private double threshold;
	private ConcurrentHashMap<String, String> locus_ids;
	private AtomicInteger counter;
	private AtomicBoolean cancel;
	private AlignmentPurpose alignmentPurpose;
	private ConcurrentLinkedQueue<String> sequencesWithoutSimilarities;
	private Map<String, Set<String>> modules;
	private String ec_number;
	private Map<String, Set<String>> closestOrthologs;
	private Map<String, Integer> kegg_taxonomy_scores;
	private int referenceTaxonomyScore;
	private double referenceTaxonomyThreshold;
	private ConcurrentLinkedQueue<AlignmentContainer> alignmentContainerSet;
	private String ko;
	private AlignmentScoreType alignmentScoreType;
	private double minAlignedResidues;

	/**
	 * Perform multiple sequence alignments.
	 * 
	 * @param method
	 * @param concurrentQueryMap
	 * @param staticSubjectMap
	 * @param queryArray
	 * @param dba
	 * @param threshold
	 * @param numberOfHelicesMap
	 * @param locus_ids
	 * @param counter
	 * @param cancel
	 * @param alignmentPurpose
	 * @param alignmentScoreType
	 */
	public PairwiseSequenceAlignement(Method method, ConcurrentHashMap<String, ProteinSequence> concurrentQueryMap, Map<String, ProteinSequence> staticSubjectMap, 
			ConcurrentLinkedQueue<String> queryArray, DatabaseAccess dba, double  threshold, Map<String, Integer> numberOfHelicesMap, 
			ConcurrentHashMap<String, String> locus_ids, AtomicInteger counter, AtomicBoolean cancel, AlignmentPurpose alignmentPurpose, AlignmentScoreType alignmentScoreType) {

		this.method=method;
		this.queryArray = queryArray;
		this.staticSubjectMap = staticSubjectMap;
		this.concurrentQueryMap = concurrentQueryMap;
		this.databaseAccess = dba;
		this.threshold=threshold;
		this.numberOfHelicesMap=numberOfHelicesMap;
		this.locus_ids = locus_ids;
		this.setCounter(counter);
		this.cancel = cancel;
		this.alignmentPurpose = alignmentPurpose;
		this.alignmentScoreType = alignmentScoreType; 
	}

	/**
	 * Perform multiple sequence alignments.
	 * 
	 * @param method
	 * @param concurrentQueryMap
	 * @param staticSubjectMap
	 * @param queryArray
	 * @param dba
	 * @param threshold
	 * @param numberOfHelicesMap
	 * @param locus_ids
	 * @param alignmentPurpose
	 * @param alignmentScoreType
	 */
	public PairwiseSequenceAlignement(Method method, ConcurrentHashMap<String, ProteinSequence> concurrentQueryMap, 
			Map<String, ProteinSequence> staticSubjectMap, ConcurrentLinkedQueue<String> queryArray, 
			DatabaseAccess dba, double  threshold, Map<String, Integer> numberOfHelicesMap, 
			ConcurrentHashMap<String, String> locus_ids, AlignmentPurpose alignmentPurpose, AlignmentScoreType alignmentScoreType) {

		this.method=method;
		this.queryArray = queryArray;
		this.staticSubjectMap = staticSubjectMap;
		this.concurrentQueryMap = concurrentQueryMap;
		this.databaseAccess = dba;
		this.threshold=threshold;
		this.numberOfHelicesMap=numberOfHelicesMap;
		this.locus_ids = locus_ids;
		this.counter = new AtomicInteger(0);
		this.cancel = new AtomicBoolean(false);
		this.alignmentPurpose = alignmentPurpose;
		this.alignmentScoreType = alignmentScoreType; 
	}

	@Override
	public void run() {

		final int size = this.queryArray.size();

		while(this.queryArray.size()>0 && !this.cancel.get()) {

			String query = this.queryArray.poll();
			try {

				if(this.alignmentPurpose.equals(AlignmentPurpose.TRANSPORT)) {

					Connection conn = new Connection(this.databaseAccess);
					this.getSimilarityTransport(query,conn);
					conn.closeConnection();
				}
				else if(this.alignmentPurpose.equals(AlignmentPurpose.ORTHOLOGS)) {

					Connection conn = new Connection(this.databaseAccess);
					this.getSimilarityOrthologs(query,conn);
					conn.closeConnection();
				}
				else if(this.alignmentPurpose.equals(AlignmentPurpose.OTHER)) {

					this.getSimilarityOther(query);
				}
			}
			catch (Exception e) {

				this.queryArray.add(query);
				e.printStackTrace();
			}
			catch (OutOfMemoryError oue) {

				this.queryArray.add(query);
				oue.printStackTrace();
			}
			System.gc();

			this.counter.incrementAndGet();
			setChanged();
			notifyObservers();
		}

		if(this.cancel.get())
			this.counter.set(size);

		setChanged();
		notifyObservers();
	}


	/**
	 * @param query
	 * @param conn
	 * @throws SQLException 
	 */
	private void getSimilarityTransport(String query, Connection conn) throws SQLException {

		ProteinSequence querySequence= this.concurrentQueryMap.get(query);
		int seqLength = querySequence.getLength();
		Matrix matrix;
		short gapOpenPenalty=10, gapExtensionPenalty=1;

		//		if(seqLength<35){matrix=Matrix.PAM30;gapOpenPenalty=9;}
		//		else if(seqLength<50){matrix=Matrix.PAM70;}
		//		else 
		if(seqLength<85)
			matrix=Matrix.BLOSUM80;
		else
			matrix=Matrix.BLOSUM62;

		AbstractPairwiseSequenceAligner<ProteinSequence,AminoAcidCompound> alignmentMethod;
		GapPenalty gp = new SimpleGapPenalty(gapOpenPenalty ,gapExtensionPenalty);
		CompoundSet<AminoAcidCompound> aa = new AminoAcidCompoundSet();

		Reader rd = new InputStreamReader(getClass().getClassLoader().getResourceAsStream(matrix.getPath()));
		SubstitutionMatrix<AminoAcidCompound> sb = new SimpleSubstitutionMatrix<AminoAcidCompound>(aa, rd,matrix.getPath());

		double helicesDependentSimilarity=this.threshold;
		ProteinSequence tcdbAAsequence=null;

		for(String tcdbRecord: this.staticSubjectMap.keySet()) {

			if(!this.cancel.get()) {

				tcdbAAsequence = this.staticSubjectMap.get(tcdbRecord);

				if(this.method.equals(Method.SmithWaterman))
					alignmentMethod=new SmithWaterman<ProteinSequence,AminoAcidCompound>(querySequence, tcdbAAsequence, gp, sb);
				else
					alignmentMethod=new NeedlemanWunsch<ProteinSequence,AminoAcidCompound>(querySequence, tcdbAAsequence, gp, sb);

				double alignmentScore = alignmentMethod.getSimilarity(); //(((double)alignmentMethod.getScore()-alignmentMethod.getMinScore())/(alignmentMethod.getMaxScore()-alignmentMethod.getMinScore()))
				double similarityScore = ((double)alignmentMethod.getPair().getNumSimilars()/alignmentMethod.getPair().getLength());
				double identityScore = ((double)alignmentMethod.getPair().getNumIdenticals()/alignmentMethod.getPair().getLength());

				double score = -1;
				if(this.method.equals(PairwiseSequenceScorerType.LOCAL))
					score = alignmentScore;
				else if(this.method.equals(PairwiseSequenceScorerType.LOCAL_IDENTITIES))
					score = identityScore;
				else if(this.method.equals(PairwiseSequenceScorerType.LOCAL_SIMILARITIES))
					score = similarityScore;

				if(this.numberOfHelicesMap.get(query)>5)
					helicesDependentSimilarity=(1-((this.numberOfHelicesMap.get(query)/2-2)*0.1))*this.threshold;

				if(score>helicesDependentSimilarity && score>(this.threshold/2)) {

					if(this.sequencesWithoutSimilarities!=null && this.sequencesWithoutSimilarities.contains(query))
						this.sequencesWithoutSimilarities.remove(query);

					String[] similarityData = new String[6];

					similarityData[0]= new StringTokenizer(query," ").nextToken();

					StringTokenizer st = new StringTokenizer(this.staticSubjectMap.get(tcdbRecord).getOriginalHeader(),"|");
					st.nextToken();
					st.nextToken();
					similarityData[1]= st.nextToken().toUpperCase();
					similarityData[2]= st.nextToken().split(" ")[0].toUpperCase();

					similarityData[3]= alignmentScore+"";

					similarityData[4]= matrix.toString();

					similarityData[5]= this.numberOfHelicesMap.get(query).toString();

					this.loadTransportData(similarityData, conn);

					similarityData=null;
				}
				alignmentMethod=null;
			}
		}

		if(!this.cancel.get())
			TransportersAPI.setProcessed(this.locus_ids.get(query), conn, DatabaseProgressStatus.PROCESSED.toString());
	}

	/**
	 * @param data
	 * @param conn
	 * @throws SQLException
	 */
	public void loadTransportData(String[] data, Connection conn) throws SQLException {

		String idLT = this.locus_ids.get(data[0]);
		Statement stmt = conn.createStatement();
		ResultSet rs = stmt.executeQuery("SELECT id FROM sw_hits WHERE tcdb_id='"+data[2]+"' AND acc='"+data[1]+"'");

		if(!rs.next()) {

			stmt.execute("INSERT INTO sw_hits (acc,tcdb_id) VALUES ('"+data[1]+"', '"+data[2]+"')");
			rs = stmt.executeQuery("SELECT LAST_INSERT_ID();");
			rs.next();
		}

		String idHIT = rs.getString(1);

		rs = stmt.executeQuery("SELECT * FROM sw_similarities WHERE sw_report_id="+idLT+" AND sw_hit_id="+idHIT+"");

		if(!rs.next()) {

			stmt.execute("INSERT INTO sw_similarities (sw_report_id,sw_hit_id,similarity) VALUES("+idLT+","+idHIT+","+data[3]+")");
		}

		stmt.close();
		stmt=null;
		rs.close();
	}

	/**
	 * @param query_id
	 * @param conn
	 * @throws SQLException 
	 */
	private void getSimilarityOrthologs(String query, Connection conn) throws SQLException {

		try {

			String [] query_array = query.split(":"); 

			String query_org = query_array [0].trim();
			String queryLocus = query_array[1].trim();

			List<String> locusTags = ModelAPI.checkDatabase(conn, queryLocus);

			if(locusTags.isEmpty()) {

				ProteinSequence querySequence= this.concurrentQueryMap.get(query);
				int seqLength = querySequence.getLength();
				Matrix matrix;
				short gapOpenPenalty=10, gapExtensionPenalty=1;

				if(seqLength<85){matrix=Matrix.BLOSUM80;}
				else{matrix=Matrix.BLOSUM62;}

				AbstractPairwiseSequenceAligner<ProteinSequence,AminoAcidCompound> alignmentMethod;
				GapPenalty gp = new SimpleGapPenalty(gapOpenPenalty ,gapExtensionPenalty);
				CompoundSet<AminoAcidCompound> aa = new AminoAcidCompoundSet();

				Reader rd = new InputStreamReader(getClass().getClassLoader().getResourceAsStream(matrix.getPath()));
				SubstitutionMatrix<AminoAcidCompound> sb = new SimpleSubstitutionMatrix<AminoAcidCompound>(aa, rd,matrix.getPath());

				ProteinSequence genomeAAsequence=null;

				double threshold = this.threshold;

				if(this.kegg_taxonomy_scores.get(query_org)>=this.referenceTaxonomyScore) {

					System.out.println("Using reference taxonomy:\t"+this.referenceTaxonomyScore+"\tthreshold\t"+this.referenceTaxonomyThreshold);
					threshold = this.referenceTaxonomyThreshold;
				}

				if(querySequence.getLength()>0) {

					for(String genome: this.staticSubjectMap.keySet()) {

						if(!this.cancel.get()) {

							genomeAAsequence = this.staticSubjectMap.get(genome);

							if(this.method.equals(Method.SmithWaterman))
								alignmentMethod=new SmithWaterman<ProteinSequence,AminoAcidCompound>(querySequence, genomeAAsequence, gp, sb);
							else
								alignmentMethod=new NeedlemanWunsch<ProteinSequence,AminoAcidCompound>(querySequence, genomeAAsequence, gp, sb);

							double alignmentScore = alignmentMethod.getSimilarity(); //(((double)alignmentMethod.getScore()-alignmentMethod.getMinScore())/(alignmentMethod.getMaxScore()-alignmentMethod.getMinScore()))
							double similarityScore = ((double)alignmentMethod.getPair().getNumSimilars()/alignmentMethod.getPair().getLength());
							double identityScore = ((double)alignmentMethod.getPair().getNumIdenticals()/alignmentMethod.getPair().getLength());

							double score = -1;
							if(this.method.equals(PairwiseSequenceScorerType.LOCAL))
								score = alignmentScore;
							else if(this.method.equals(PairwiseSequenceScorerType.LOCAL_IDENTITIES))
								score = identityScore;
							else if(this.method.equals(PairwiseSequenceScorerType.LOCAL_SIMILARITIES))
								score = similarityScore;

							//TODO replace String array by Container

							if(score > threshold) {

								//no_similarity = false;

								this.sequencesWithoutSimilarities.remove(query);

								String[] similarityData = new String[4];

								similarityData[0]= query;
								similarityData[1]= this.staticSubjectMap.get(genome).getOriginalHeader();
								similarityData[2]= alignmentScore+"";
								similarityData[3]= matrix.toString();

								PairwiseSequenceAlignement.loadOrthologsData(similarityData, conn, ec_number, closestOrthologs, modules, this.databaseAccess.get_database_type());

								similarityData=null;
							}
							alignmentMethod=null;
						}
					}
				}
			}
			else {

				this.sequencesWithoutSimilarities.remove(query);

				for(String locus : locusTags) {

					String[] similarityData = new String[4];

					similarityData[0]= query;
					similarityData[1]= locus;
					similarityData[2]= null;
					similarityData[3]= null;

					PairwiseSequenceAlignement.loadOrthologsData(similarityData, conn, ec_number, closestOrthologs, modules, this.databaseAccess.get_database_type());
				}
			}
		}
		catch (Exception e) {

			e.printStackTrace();
			System.err.println(query);
			System.err.println();
			System.out.println();
		}
	}

	/**
	 * @param query
	 * @param conn
	 * @throws SQLException
	 */
	private void getSimilarityOther(String query) throws SQLException {

		try {

			String [] query_array = query.split(":"); 
			String query_org = query_array [0].trim();
			//String queryLocus = query_array[1].trim();

			ProteinSequence querySequence= this.concurrentQueryMap.get(query);
			int seqLength = querySequence.getLength();
			
			//parameters
			Matrix matrix;
			short gapOpenPenalty=10, gapExtensionPenalty=1;

			if(seqLength<85){matrix=Matrix.BLOSUM80;}
			else{matrix=Matrix.BLOSUM62;}

			GapPenalty gp = new SimpleGapPenalty(gapOpenPenalty ,gapExtensionPenalty);
			CompoundSet<AminoAcidCompound> aa = new AminoAcidCompoundSet();
			Reader rd = new InputStreamReader(getClass().getClassLoader().getResourceAsStream(matrix.getPath()));
			SubstitutionMatrix<AminoAcidCompound> sb = new SimpleSubstitutionMatrix<AminoAcidCompound>(aa, rd,matrix.getPath());

			ProteinSequence genomeAAsequence=null;

			double workingThreshold = this.threshold;

			if(this.kegg_taxonomy_scores.containsKey(query_org) && 
					this.kegg_taxonomy_scores.get(query_org)>=this.referenceTaxonomyScore) {

				System.out.println("Using reference taxonomy:\t"+this.referenceTaxonomyScore+"\tthreshold\t"+this.referenceTaxonomyThreshold);
				workingThreshold = this.referenceTaxonomyThreshold;
			}

			//alignment
			PairwiseSequenceAlignerType alignmentType = PairwiseSequenceAlignerType.LOCAL;
			if(this.method.equals(Method.NeedlemanWunsch))
				alignmentType = PairwiseSequenceAlignerType.GLOBAL;

			if(querySequence.getLength()>0) {

				for(String genome: this.staticSubjectMap.keySet()) {

					double verifyingThreshold = workingThreshold;

					try {

						genomeAAsequence = this.staticSubjectMap.get(genome);

						PairwiseSequenceAligner<ProteinSequence, AminoAcidCompound> result = Alignments.getPairwiseAligner(querySequence, genomeAAsequence, alignmentType, gp, sb);

						double alignmentScore = result.getSimilarity();
						double similarityScore = ((double) result.getPair().getNumSimilars()/result.getPair().getLength());
						double identityScore = ((double)result.getPair().getNumIdenticals()/result.getPair().getLength());

						double score = -1;
						boolean go = false;

						if(this.alignmentScoreType.equals(AlignmentScoreType.ALIGNMENT)) {

							score = alignmentScore;
							go = score > verifyingThreshold;
						}
						else if(this.alignmentScoreType.equals(AlignmentScoreType.IDENTITY)) {

							double minResidues = PairwiseSequenceAlignement._MIN_RESIDUES;
							score = identityScore;

							if(result.getPair().getLength()>=minResidues) {

								go = score > verifyingThreshold;
							}
							else {

								while (result.getPair().getLength()<minResidues && minResidues>this.minAlignedResidues) {

									minResidues = minResidues/PairwiseSequenceAlignement._LENGTH_DIVISION;;
									verifyingThreshold += PairwiseSequenceAlignement._SCORE_INCREMENT;
									
									if(verifyingThreshold>1)
										verifyingThreshold=1;
								}

								if(result.getPair().getLength()> this.minAlignedResidues)
									go = score > verifyingThreshold;
							}
						}
						else if(this.alignmentScoreType.equals(AlignmentScoreType.SIMILARITY)) {

							score = similarityScore;
							go = score > verifyingThreshold;
						}

						if(go) {

							AlignmentContainer alignmentContainer = new AlignmentContainer(query, this.staticSubjectMap.get(genome).getOriginalHeader(), ko, 
									score, result.getMaxScore(), result.getMinScore(), result.getScore(), 
									result.getPair().getNumIdenticals(), result.getPair().getNumSimilars(), 
									result.getPair().getLength(), result.getPair().getQuery().getLength(), 
									result.getPair().getTarget().getLength(), matrix.toString(), this.method, this.alignmentScoreType);
							alignmentContainerSet.add(alignmentContainer);
						}
					} 
					catch (OutOfMemoryError ee) {

						System.err.println("query "+query);
						System.err.println("query "+querySequence);
						System.out.println("genome "+genomeAAsequence);
						ee.printStackTrace();
					}
				}
			}
		}
		catch (Exception e) {

			e.printStackTrace();
			System.err.println(query);
			System.err.println();
			System.out.println();
		}
	}

	//	/**
	//	 * @param conn
	//	 * @param query
	//	 * @throws SQLException 
	//	 */
	//	private void addNoSimilarityOrtholog(Connection conn, String query) throws SQLException {
	//
	//		Statement stmt = conn.createStatement();
	//		stmt.execute("INSERT INTO orthology (entry_id) VALUES ('"+query+"')");
	//		stmt.close();
	//		stmt=null;
	//
	//	}

	/**
	 * @param data
	 * @param conn
	 * @param ec_number
	 * @param closestOrthologs
	 * @param modules
	 * @param databaseType
	 * @throws SQLException
	 */
	public static void loadOrthologsData(String[] data, Connection conn, String ec_number, 
			Map<String, Set<String>> closestOrthologs, Map<String, Set<String>> modules, DatabaseType databaseType) throws SQLException {

		String queryLocus = data[0].split(":")[1];
		//FIXME
		System.out.println();
		System.out.println(data[1]);
		System.out.println(data[2]);
		System.out.println(data[3]);
		System.out.println(closestOrthologs.get(data[0]));
		//System.out.println(this.modules.get(data[0]));

		Statement statement = conn.createStatement();

		synchronized (conn) {

			ResultSet rs = null;
			Pair<String, String> geneNames = HomologyAPI.getGeneLocusFromHomologyData(data[1], statement);
			String idGene = ModelAPI.loadGene(geneNames, data[1], null, statement, "KO",databaseType);


			//			if(data[1]!=null) {
			//
			//				rs = stmt.executeQuery("SELECT idgene FROM gene WHERE sequence_id ='"+data[1]+"'");
			//
			//				if(!rs.next()) {
			//
			//					//					rs = stmt.executeQuery("SELECT idchromosome FROM chromosome WHERE name ='unknown'");
			//					//
			//					//					if(!rs.next()) {
			//					//
			//					//						stmt.execute("INSERT INTO chromosome (name) VALUES ('unknown')");
			//					//						rs = stmt.executeQuery("SELECT LAST_INSERT_ID();");
			//					//						rs.next();
			//					//					}
			//					//					String idchromosome = rs.getString(1);
			//
			//					stmt.execute("INSERT INTO gene (locusTag,"
			//							//							+ " chromosome_idchromosome,"
			//							+ " sequence_id, "
			//							+ " origin) VALUES ('"+data[1]+"',"
			//							//									+ " "+idchromosome+","
			//							+ "'"+data[1]+"',"
			//							+ "'KO')");
			//					rs = stmt.executeQuery("SELECT LAST_INSERT_ID();");
			//					rs.next();
			//				}
			//				idGene = rs.getString(1);
			//			}

			if(idGene != null) {

				for (String ortholog : closestOrthologs.get(data[0])) {

					rs = statement.executeQuery("SELECT id FROM orthology WHERE entry_id ='"+ortholog+"' AND (locus_id is null OR locus_id = '')");

					String orthology_id = "";

					if(rs.next()) {

						orthology_id = rs.getString(1);
						statement.execute("UPDATE orthology SET locus_id = '"+queryLocus+"' WHERE entry_id = '"+ortholog+"';");
					}
					else {

						rs = statement.executeQuery("SELECT id FROM orthology WHERE entry_id ='"+ortholog+"' AND locus_id ='"+queryLocus+"';");

						if(!rs.next()) {

							statement.execute("INSERT INTO orthology (entry_id, locus_id) VALUES ('"+ortholog+"', '"+queryLocus+"')");
							rs = statement.executeQuery("SELECT LAST_INSERT_ID();");
							rs.next();
						}
						orthology_id = rs.getString(1);
					}

					rs = statement.executeQuery("SELECT * FROM gene_has_orthology WHERE gene_idgene='"+idGene+"' AND orthology_id='"+orthology_id+"'");

					if(rs.next()) {

						System.out.println("Entry exists!! "+idGene+"\t"+orthology_id);
					}
					else {

						statement.execute("INSERT INTO gene_has_orthology (gene_idgene,orthology_id, similarity) VALUES("+idGene+","+orthology_id+", "+data[2]+" )");
					}

					rs = statement.executeQuery("SELECT protein_idprotein FROM enzyme WHERE ecnumber='"+ec_number+"'");
					rs.next();
					int protein_idprotein = rs.getInt(1);

					rs = statement.executeQuery("SELECT module_id, note FROM subunit WHERE gene_idgene='"+idGene+"' AND enzyme_ecnumber = '"+ec_number+"'");

					List<String> modules_ids = new ArrayList<String>();
					boolean exists = false, noModules=true;

					String note = "unannotated";

					while(rs.next()) {

						exists = true;

						if(rs.getInt(1)>0) {

							noModules = false;
							modules_ids.add(rs.getString(1));
						}

						if(rs.getString(2)!=null && !rs.getString(2).equalsIgnoreCase("null"))
							note = rs.getString(2);
						else
							note = "";
					}

					for(String module_id : modules.get(ortholog)) {

						if(modules_ids.contains(module_id)) {

						}
						else {

							if(exists) {

								if(noModules) {

									statement.execute("UPDATE subunit SET module_id = "+module_id+" WHERE gene_idgene = '"+idGene+"' AND enzyme_ecnumber = '"+ec_number+"'");
									noModules = false;
									modules_ids.add(module_id);
								}
							}
							else {

								statement.execute("INSERT INTO subunit (module_id, gene_idgene, enzyme_ecnumber, enzyme_protein_idprotein, note) " +
										"VALUES("+module_id+", "+idGene+", '"+ec_number+"', "+protein_idprotein+", '"+note+"')");
								//exists = true;
							}

						}
					}
				}
			}

			rs.close();
		} 
		statement.close();
		statement=null;
	}


	/**
	 * @return the counter
	 */
	public AtomicInteger getCounter() {

		return counter;
	}

	/**
	 * @param counter the counter to set
	 */
	public void setCounter(AtomicInteger counter) {
		this.counter = counter;
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
	 * @return the idModule
	 */
	public Map<String, Set<String>> getModules() {
		return modules;
	}

	/**
	 * @param idModule the idModule to set
	 */
	public void setModules(Map<String, Set<String>> modules) {
		this.modules = modules;
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
	 * @param closestOrthologs
	 */
	public void setClosestOrthologs(Map<String, Set<String>> closestOrthologs) {

		this.closestOrthologs = closestOrthologs;
	}

	public void setKegg_taxonomy_scores(Map<String, Integer> kegg_taxonomy_scores) {

		this.kegg_taxonomy_scores = kegg_taxonomy_scores;
	}

	public void setReferenceTaxonomyScore(int referenceTaxonomyScore) {

		this.referenceTaxonomyScore = referenceTaxonomyScore;
	}

	public void setReferenceTaxonomyThreshold(double referenceTaxonomyThreshold) {

		this.referenceTaxonomyThreshold = referenceTaxonomyThreshold;
	}

	public void setAlignmentContainerSet(ConcurrentLinkedQueue<AlignmentContainer> alignmentContainerSet) {

		this.alignmentContainerSet = alignmentContainerSet;
	}


	/**
	 * @param ko
	 */
	public void setKO(String ko) {

		this.ko = ko;		
	}

	/**
	 * @return the minAlignedResidues
	 */
	public double getMinAlignedResidues() {
		return minAlignedResidues;
	}

	/**
	 * @param minAlignedResidues the minAlignedResidues to set
	 */
	public void setMinAlignedResidues(double minAlignedResidues) {
		this.minAlignedResidues = minAlignedResidues;
	}

}