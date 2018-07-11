package pt.uminho.ceb.biosystems.merlin.local.alignments.core;

import java.io.InputStreamReader;
import java.io.Reader;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
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
import org.biojava.nbio.alignment.NeedlemanWunsch;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.alignment.SmithWaterman;
import org.biojava.nbio.alignment.template.AbstractPairwiseSequenceAligner;
import org.biojava.nbio.alignment.template.GapPenalty;
import org.biojava.nbio.alignment.template.PairwiseSequenceAligner;
import org.biojava.nbio.core.alignment.matrices.SimpleSubstitutionMatrix;
import org.biojava.nbio.core.alignment.template.AlignedSequence;
import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.template.AbstractSequence;
import org.biojava.nbio.core.sequence.template.CompoundSet;

import pt.uminho.ceb.biosystems.merlin.utilities.Enumerators.AlignmentPurpose;
import pt.uminho.ceb.biosystems.merlin.utilities.Enumerators.AlignmentScoreType;
import pt.uminho.ceb.biosystems.merlin.utilities.Enumerators.Matrix;
import pt.uminho.ceb.biosystems.merlin.utilities.Enumerators.Method;
import pt.uminho.ceb.biosystems.merlin.utilities.containers.capsules.AlignmentCapsule;


/**
 * @author ODias
 *
 */
public class PairwiseSequenceAlignement extends Observable implements Runnable {

	private static int _MIN_RESIDUES = 100;
	private static double _SCORE_INCREMENT = 0.2, _LENGTH_DIVISION = 2;
	private ConcurrentLinkedQueue<String> queryArray;
	private Map<String, AbstractSequence<?>> staticSubjectMap; 
	private ConcurrentHashMap<String, AbstractSequence<?>> concurrentQueryMap;
	private Method method;
	private double threshold;
	private AtomicInteger counter;
	private AtomicBoolean cancel;
	private AlignmentPurpose alignmentPurpose;
	private ConcurrentLinkedQueue<String> sequencesWithoutSimilarities;
	private Map<String, Set<Integer>> modules;
	private String ec_number;
	private Map<String, Set<String>> closestOrthologs;
	private Map<String, Integer> kegg_taxonomy_scores;
	private int referenceTaxonomyScore;
	private double referenceTaxonomyThreshold;
	private ConcurrentLinkedQueue<AlignmentCapsule> alignmentContainerSet;
	private String ko;
	private AlignmentScoreType alignmentScoreType;
	private double minAlignedResidues;
	private Map <String, Double> querySpecificThreshold;
	private Map<String, List<String>> sequenceIdsSet;
	
	/**
	 * Perform multiple sequence alignments.
	 * 
	 * @param method
	 * @param concurrentQueryMap
	 * @param staticSubjectMap
	 * @param queryArray
	 * @param threshold
	 * @param counter
	 * @param cancel
	 * @param alignmentPurpose
	 * @param alignmentScoreType
	 * @param alignmentContainerSet
	 */
	public PairwiseSequenceAlignement(Method method, ConcurrentHashMap<String, AbstractSequence<?>> concurrentQueryMap, Map<String, AbstractSequence<?>> staticSubjectMap, 
			ConcurrentLinkedQueue<String> queryArray, double  threshold,
			AtomicInteger counter, AtomicBoolean cancel, AlignmentPurpose alignmentPurpose, AlignmentScoreType alignmentScoreType,
			ConcurrentLinkedQueue<AlignmentCapsule> alignmentContainerSet) {

		this.method=method;
		this.queryArray = queryArray;
		this.staticSubjectMap = staticSubjectMap;
		this.concurrentQueryMap = concurrentQueryMap;
		this.threshold=threshold;
		this.setCounter(counter);
		this.cancel = cancel;
		this.alignmentPurpose = alignmentPurpose;
		this.alignmentScoreType = alignmentScoreType;
		this.alignmentContainerSet = alignmentContainerSet;
		this.querySpecificThreshold = new HashMap<>();
	}
	
	@Override
	public void run() {

		final int size = this.queryArray.size();
		
		while(this.queryArray.size()>0 && !this.cancel.get()) {

			String query = this.queryArray.poll();
			try {
				
				System.out.println(this.queryArray.size());
				System.out.println("RUNNING...");
				
				this.getSimilarity(query);
//				if(this.alignmentPurpose.equals(AlignmentPurpose.TRANSPORT)) 
//					this.alignmentContainerSet.addAll(new ArrayList<>(this.getSimilarityTransport(query)));
//				
//				else if(this.alignmentPurpose.equals(AlignmentPurpose.ORTHOLOGS))
//					this.alignmentContainerSet.addAll(new ArrayList<>(this.getSimilarityOrthologs(query)));
//				
//				else if(this.alignmentPurpose.equals(AlignmentPurpose.OTHER)) 
//					this.getSimilarityOther(query);
				
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
	 * @return
	 */
	private void getSimilarity(String query){

		try {
			double workingThreshold = this.threshold;
			
			String [] query_array; 
			String query_org = "";
			String queryLocus = "";

			if(query.contains(":")) {
				query_array = query.split(":"); 
				query_org = query_array [0].trim();
				queryLocus = query_array[1].trim();
				
				if(this.kegg_taxonomy_scores.containsKey(query_org) && this.kegg_taxonomy_scores.get(query_org)>=this.referenceTaxonomyScore)			
					workingThreshold = this.referenceTaxonomyThreshold;
			}
			
			if(!this.alignmentPurpose.equals(AlignmentPurpose.ORTHOLOGS) || (!sequenceIdsSet.containsKey(queryLocus) || sequenceIdsSet.get(queryLocus).isEmpty())) {

				if(this.alignmentPurpose.equals(AlignmentPurpose.TRANSPORT)) {
					if (this.querySpecificThreshold.containsKey(query)) 
						workingThreshold = this.querySpecificThreshold.get(query);
				}
				
				AbstractSequence<?> querySequence= this.concurrentQueryMap.get(query);
				int seqLength = querySequence.getLength();
				Matrix matrix;
				short gapOpenPenalty=10, gapExtensionPenalty=1;
				
				if(seqLength<85){matrix=Matrix.BLOSUM80;}
				else{matrix=Matrix.BLOSUM62;}
				
				//Alignment
				PairwiseSequenceAlignerType alignmentType = null;
				if(this.alignmentPurpose.equals(AlignmentPurpose.OTHER)) {

					alignmentType = PairwiseSequenceAlignerType.LOCAL;
					if(this.method.equals(Method.NeedlemanWunsch))
						alignmentType = PairwiseSequenceAlignerType.GLOBAL;	
				}

				PairwiseSequenceAligner<ProteinSequence,AminoAcidCompound> alignmentMethod;
				GapPenalty gp = new SimpleGapPenalty(gapOpenPenalty ,gapExtensionPenalty);
				CompoundSet<AminoAcidCompound> aa = new AminoAcidCompoundSet();

				Reader rd = new InputStreamReader(getClass().getClassLoader().getResourceAsStream(matrix.getPath()));
				SubstitutionMatrix<AminoAcidCompound> sb = new SimpleSubstitutionMatrix<AminoAcidCompound>(aa, rd,matrix.getPath());

				AbstractSequence<?> subjectSequence=null;

				if(seqLength>0) {

					for(String gene: this.staticSubjectMap.keySet()) {
						
						if(!this.cancel.get()) {

							try {
								
								subjectSequence = this.staticSubjectMap.get(gene);
								
								double verifyingThreshold = workingThreshold;

								if(this.alignmentPurpose.equals(AlignmentPurpose.OTHER)) {

									alignmentMethod = Alignments.getPairwiseAligner(new ProteinSequence(querySequence.getSequenceAsString())
													, new ProteinSequence(subjectSequence.getSequenceAsString()), alignmentType, gp, sb);
								}
								else {

									if(this.method.equals(Method.SmithWaterman))
										alignmentMethod=new SmithWaterman<ProteinSequence,AminoAcidCompound>(new ProteinSequence(querySequence.getSequenceAsString()), new ProteinSequence(subjectSequence.getSequenceAsString()), gp, sb);
									else
										alignmentMethod=new NeedlemanWunsch<ProteinSequence,AminoAcidCompound>(new ProteinSequence(querySequence.getSequenceAsString()), new ProteinSequence(subjectSequence.getSequenceAsString()), gp, sb);
								}

								double alignmentScore = alignmentMethod.getSimilarity(); //(((double)alignmentMethod.getScore()-alignmentMethod.getMinScore())/(alignmentMethod.getMaxScore()-alignmentMethod.getMinScore()))
								double similarityScore = ((double)alignmentMethod.getPair().getNumSimilars()/alignmentMethod.getPair().getLength());
								double identityScore = alignmentMethod.getPair().getPercentageOfIdentity();//((double)alignmentMethod.getPair().getNumIdenticals()/alignmentMethod.getPair().getLength());

								
								double score = -1;
								boolean go = false;

								if(this.alignmentScoreType.equals(AlignmentScoreType.ALIGNMENT)) {
									score = alignmentScore;
									go = score > verifyingThreshold;
								}
								else if(this.alignmentScoreType.equals(AlignmentScoreType.IDENTITY)) {

									if(this.alignmentPurpose.equals(AlignmentPurpose.OTHER)) {
										double minResidues = PairwiseSequenceAlignement._MIN_RESIDUES;
										score = identityScore;

										if(alignmentMethod.getPair().getLength()>=minResidues) {

											go = score > verifyingThreshold;
										}
										else {

											while (alignmentMethod.getPair().getLength()<minResidues && minResidues>this.minAlignedResidues) {

												minResidues = minResidues/PairwiseSequenceAlignement._LENGTH_DIVISION;;
												verifyingThreshold += PairwiseSequenceAlignement._SCORE_INCREMENT;

												if(verifyingThreshold>1)
													verifyingThreshold=1;
											}

											if(alignmentMethod.getPair().getLength()> this.minAlignedResidues)
												go = score > verifyingThreshold;
										}
									}
									else {
										score = identityScore;
										go = score > verifyingThreshold;
									}
								}
								else if(this.alignmentScoreType.equals(AlignmentScoreType.SIMILARITY)) {
									score = similarityScore;
									go = score > verifyingThreshold;
								}
								
								if(go) {

									if(this.sequencesWithoutSimilarities!=null && this.sequencesWithoutSimilarities.contains(query))
										this.sequencesWithoutSimilarities.remove(query);

									double coverage1 = alignmentMethod.getPair().getAlignedSequence(1).getCoverage();
									double coverage2 = alignmentMethod.getPair().getAlignedSequence(2).getCoverage();
									
									String alignQuery = query;
									String tcdbID = "";
									String target = subjectSequence.getOriginalHeader();
									
									if(this.alignmentPurpose.equals(AlignmentPurpose.TRANSPORT)) {

										alignQuery = new StringTokenizer(query," ").nextToken();

										StringTokenizer st = new StringTokenizer(subjectSequence.getOriginalHeader(),"|");
										st.nextToken();
										st.nextToken();

										target = st.nextToken().toUpperCase();
										tcdbID = st.nextToken().split(" ")[0].toUpperCase();
									}
									
									AlignmentCapsule container = new AlignmentCapsule(alignQuery, target, ko, score,
											alignmentMethod.getMaxScore(), alignmentMethod.getMinScore(), alignmentMethod.getScore(), 
											alignmentMethod.getPair().getNumIdenticals(), alignmentMethod.getPair().getNumSimilars(), 
											alignmentMethod.getPair().getLength(), querySequence.getLength(), subjectSequence.getLength(),
											coverage1, coverage2, alignmentMethod.getPair().getAlignedSequence(1).getNumGapPositions(), 
											alignmentMethod.getPair().getAlignedSequence(2).getNumGapPositions(),
											matrix.toString(), this.method, this.alignmentScoreType);
									
									container.setTcdbID(tcdbID);
									container.setEcNumber(ec_number);
									container.setClosestOrthologues(closestOrthologs);
									container.setModules(modules);

									this.alignmentContainerSet.add(container);
								}
								alignmentMethod=null;
							}
							catch (OutOfMemoryError ee) {

								System.err.println("query "+query);
								System.err.println("query "+querySequence);
								ee.printStackTrace();
							}
						}
					}
				}
			}
			else {
				if(this.sequencesWithoutSimilarities!=null && this.sequencesWithoutSimilarities.contains(query))
					this.sequencesWithoutSimilarities.remove(query);;
			}
		}
		catch (Exception e) {

			e.printStackTrace();
			System.err.println(query);
			System.err.println();
		}
	}


	/**
	 * @param query
	 * @param conn
	 * @throws SQLException 
	 * @throws CompoundNotFoundException 
	 */
	private List<AlignmentCapsule> getSimilarityTransport(String query) throws SQLException, CompoundNotFoundException {

		List<AlignmentCapsule> similarityData = new ArrayList<>();
		
		AbstractSequence<?> querySequence= this.concurrentQueryMap.get(query);
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

		AbstractSequence<?> tcdbAAsequence=null;

		for(String tcdbRecord: this.staticSubjectMap.keySet()) {

			if(!this.cancel.get()) {

				tcdbAAsequence = this.staticSubjectMap.get(tcdbRecord);
				
				if(this.method.equals(Method.SmithWaterman))
					alignmentMethod=new SmithWaterman<ProteinSequence,AminoAcidCompound>(new ProteinSequence(querySequence.getSequenceAsString()), new ProteinSequence(tcdbAAsequence.getSequenceAsString()), gp, sb);
				else
					alignmentMethod=new NeedlemanWunsch<ProteinSequence,AminoAcidCompound>(new ProteinSequence(querySequence.getSequenceAsString()), new ProteinSequence(tcdbAAsequence.getSequenceAsString()), gp, sb);

				double alignmentScore = alignmentMethod.getSimilarity(); //(((double)alignmentMethod.getScore()-alignmentMethod.getMinScore())/(alignmentMethod.getMaxScore()-alignmentMethod.getMinScore()))
				double similarityScore = ((double)alignmentMethod.getPair().getNumSimilars()/alignmentMethod.getPair().getLength());
				double identityScore = ((double)alignmentMethod.getPair().getNumIdenticals()/alignmentMethod.getPair().getLength());

				double score = -1;
				if(this.alignmentScoreType.equals(AlignmentScoreType.ALIGNMENT))
					score = alignmentScore;
				else if(this.alignmentScoreType.equals(AlignmentScoreType.IDENTITY))
					score = identityScore;
				else if(this.alignmentScoreType.equals(AlignmentScoreType.SIMILARITY))
					score = similarityScore;
				
				double specificThreshold = this.threshold;
				if (this.querySpecificThreshold.containsKey(query))
					specificThreshold = this.querySpecificThreshold.get(query);

				if(score>specificThreshold) {

					if(this.sequencesWithoutSimilarities!=null && this.sequencesWithoutSimilarities.contains(query))
						this.sequencesWithoutSimilarities.remove(query);

					String query2 = new StringTokenizer(query," ").nextToken();

					StringTokenizer st = new StringTokenizer(this.staticSubjectMap.get(tcdbRecord).getOriginalHeader(),"|");
					st.nextToken();
					st.nextToken();
					
					String target = st.nextToken().toUpperCase();
					String tcdbID = st.nextToken().split(" ")[0].toUpperCase();

//					similarityData[3]= alignmentScore+"";

//					similarityData[4]= matrix.toString();

//					similarityData[5]= this.numberOfHelicesMap.get(query).toString();			VERIFICAR SE ESTE FAZ FALTA!!!!!!!!!!!!!!
					
					AlignmentCapsule container = new AlignmentCapsule(query2, target, tcdbID, matrix.toString(), alignmentScore);
					similarityData.add(container);

//					this.loadTransportData(similarityData, conn);

//					similarityData=null;
				}
				alignmentMethod=null;
			}
		}

		return similarityData;
	}



	/**
	 * @param query_id
	 * @param conn
	 * @throws SQLException 
	 */
	private List<AlignmentCapsule> getSimilarityOrthologs(String query) throws SQLException {

		try {
			
			List<AlignmentCapsule> similarityData = new ArrayList<>();
			
			String [] query_array = query.split(":"); 

			String query_org = query_array [0].trim();
			String queryLocus = query_array[1].trim();
			
			if(!sequenceIdsSet.containsKey(queryLocus) || sequenceIdsSet.get(queryLocus).isEmpty()) {
				
				AbstractSequence<?> querySequence= this.concurrentQueryMap.get(query);
				
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

				AbstractSequence<?> genomeAAsequence=null;

				double threshold = this.threshold;
				
				if(this.kegg_taxonomy_scores.get(query_org)>=this.referenceTaxonomyScore) 
					threshold = this.referenceTaxonomyThreshold;
				
				if(seqLength>0) {

					for(String gene: this.staticSubjectMap.keySet()) {
						
						if(!this.cancel.get()) {

							genomeAAsequence = this.staticSubjectMap.get(gene);

							if(this.method.equals(Method.SmithWaterman))
								alignmentMethod=new SmithWaterman<ProteinSequence,AminoAcidCompound>(new ProteinSequence(querySequence.getSequenceAsString()),new ProteinSequence(genomeAAsequence.getSequenceAsString()), gp, sb);
							else
								alignmentMethod=new NeedlemanWunsch<ProteinSequence,AminoAcidCompound>(new ProteinSequence(querySequence.getSequenceAsString()), new ProteinSequence(genomeAAsequence.getSequenceAsString()), gp, sb);
							
							double alignmentScore = alignmentMethod.getSimilarity(); //(((double)alignmentMethod.getScore()-alignmentMethod.getMinScore())/(alignmentMethod.getMaxScore()-alignmentMethod.getMinScore()))
							double similarityScore = ((double)alignmentMethod.getPair().getNumSimilars()/alignmentMethod.getPair().getLength());
							double identityScore = ((double)alignmentMethod.getPair().getNumIdenticals()/alignmentMethod.getPair().getLength());

							double score = -1;
							if(this.alignmentScoreType.equals(AlignmentScoreType.ALIGNMENT))
								score = alignmentScore;
							else if(this.alignmentScoreType.equals(AlignmentScoreType.IDENTITY))
								score = identityScore;
							else if(this.alignmentScoreType.equals(AlignmentScoreType.SIMILARITY))
								score = similarityScore;
							
							
							if(score > threshold) {
								
								this.sequencesWithoutSimilarities.remove(query);
								
								AlignedSequence<ProteinSequence, AminoAcidCompound> a = alignmentMethod.getPair().getAlignedSequence(1);
								AlignedSequence<ProteinSequence, AminoAcidCompound> b = alignmentMethod.getPair().getAlignedSequence(2);
								
								double c1 = a.getCoverage();
								double c2 = b.getCoverage();
								
								AlignmentCapsule container = new AlignmentCapsule(query, c1, c2, this.staticSubjectMap.get(gene).getOriginalHeader(), 
										alignmentScore, ec_number, closestOrthologs, modules);
								similarityData.add(container);
							}
							alignmentMethod=null;
						}
					}
				}
			}
			else 
				this.sequencesWithoutSimilarities.remove(query);
			
			return similarityData;
		}
		catch (Exception e) {

			e.printStackTrace();
		}
		return null;
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

			AbstractSequence<?> querySequence= this.concurrentQueryMap.get(query);
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

			AbstractSequence<?> genomeAAsequence=null;

			double workingThreshold = this.threshold;

			if(this.kegg_taxonomy_scores.containsKey(query_org) && 
					this.kegg_taxonomy_scores.get(query_org)>=this.referenceTaxonomyScore) {

				workingThreshold = this.referenceTaxonomyThreshold;
			}

			//alignment
			PairwiseSequenceAlignerType alignmentType = PairwiseSequenceAlignerType.LOCAL;
			if(this.method.equals(Method.NeedlemanWunsch))
				alignmentType = PairwiseSequenceAlignerType.GLOBAL;

			if(seqLength>0) {

				for(String genome: this.staticSubjectMap.keySet()) {

					double verifyingThreshold = workingThreshold;

					try {

						genomeAAsequence = this.staticSubjectMap.get(genome);

						PairwiseSequenceAligner<ProteinSequence, AminoAcidCompound> result = Alignments.getPairwiseAligner(new ProteinSequence(querySequence.getSequenceAsString()), new ProteinSequence(genomeAAsequence.getSequenceAsString()), alignmentType, gp, sb);

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

							AlignmentCapsule alignmentContainer = new AlignmentCapsule(query, this.staticSubjectMap.get(genome).getOriginalHeader(), ko, 
									score, result.getMaxScore(), result.getMinScore(), result.getScore(), 
									result.getPair().getNumIdenticals(), result.getPair().getNumSimilars(), 
									result.getPair().getLength(), querySequence.getLength(), genomeAAsequence.getLength(),
									result.getPair().getAlignedSequence(1).getCoverage(), result.getPair().getAlignedSequence(2).getCoverage(),
									result.getPair().getAlignedSequence(1).getNumGapPositions(), result.getPair().getAlignedSequence(2).getNumGapPositions(),
									matrix.toString(), this.method, this.alignmentScoreType);
							alignmentContainerSet.add(alignmentContainer);
						}
					} 
					catch (OutOfMemoryError ee) {

						System.err.println("query "+query);
						System.err.println("query "+querySequence);
						ee.printStackTrace();
					}
				}
			}
		}
		catch (Exception e) {

			e.printStackTrace();
			System.err.println(query);
			System.err.println();
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
	public Map<String, Set<Integer>> getModules() {
		return modules;
	}

	/**
	 * @param idModule the idModule to set
	 */
	public void setModules(Map<String, Set<Integer>> modules) {
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

	/**
	 * @return the querySpecificThreshold
	 */
	public Map<String, Double> getQuerySpecificThreshold() {
		return querySpecificThreshold;
	}

	/**
	 * @param querySpecificThreshold the querySpecificThreshold to set
	 */
	public void setQuerySpecificThreshold(Map<String, Double> querySpecificThreshold) {
		this.querySpecificThreshold = querySpecificThreshold;
	}

	/**
	 * @return the sequenceIdsSet
	 */
	public Map<String, List<String>> getSequenceIdsSet() {
		return sequenceIdsSet;
	}

	/**
	 * @param sequenceIdsSet the sequenceIdsSet to set
	 */
	public void setSequenceIdsSet(Map<String, List<String>> sequenceIdsSet) {
		this.sequenceIdsSet = sequenceIdsSet;
	}

//	public ConcurrentHashMap<String, List<AlignmentContainer>> getAlignmentResults() {
//		return alignmentResults;
//	}
//
//	public void setAlignmentResults(ConcurrentHashMap<String, List<AlignmentContainer>> alignmentResults) {
//		this.alignmentResults = alignmentResults;
//	}

}