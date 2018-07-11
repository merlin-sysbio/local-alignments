package pt.uminho.ceb.biosystems.merlin.local.alignments.core.ModelMerge;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Observable;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.atomic.AtomicBoolean;

import javax.xml.bind.JAXBContext;

import org.biojava.nbio.core.sequence.template.AbstractSequence;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import pt.uminho.ceb.biosystems.merlin.local.alignments.core.AlignmentsUtils;
import pt.uminho.ceb.biosystems.merlin.utilities.Enumerators.AlignmentPurpose;
import pt.uminho.ceb.biosystems.merlin.utilities.blast.ncbi_blastparser.BlastIterationData;
import pt.uminho.ceb.biosystems.merlin.utilities.blast.ncbi_blastparser.Hit;
import pt.uminho.ceb.biosystems.merlin.utilities.blast.ncbi_blastparser.NcbiBlastParser;
import pt.uminho.ceb.biosystems.merlin.utilities.containers.capsules.AlignmentCapsule;


/**
 * @author amaromorais
 *
 */
public class BlastAlignment extends Observable implements ModelAlignments{

	private static final double FIXED_THRESHOLD =  10e-6;

	private static final double ALIGNMENT_MIN_SCORE = 0.0;
	private static final double BITSCORE_THRESHOLD = 50;
	private static final double COVERAGE_THRESHOLD = 0.20;
	private static final double ALIGNMENT_QUERY_LEN_THRESHOLD = 0.25;
	private static final double QUERY_HIT_LEN_THRESHOLD = 0.25;

	final static Logger logger = LoggerFactory.getLogger(BlastAlignment.class);


	private NcbiBlastParser blout;
	private ConcurrentLinkedQueue<AlignmentCapsule> alignmentContainerSet;
	private String alignmentMatrix, queryFasta, subjectFasta;
	private boolean isTransportersSearch = false;
	private double threshold;
	private AtomicBoolean cancel; 
	private Map<String,AbstractSequence<?>> querySequences;
	private JAXBContext jc;
	private String ec_number;
	private Map<String,Set<String>> closestOrthologs;
	private Map<String,Set<Integer>> modules;
	private ConcurrentLinkedQueue<String> sequencesWithoutSimilarities;
	private AlignmentPurpose blastPurpose;

	private Map<String, List<String>> sequenceIdsSet;
	private Map<String, Integer> kegg_taxonomy_scores;
	private Integer referenceTaxonomyScore;
	private Double referenceTaxonomyThreshold;


	public BlastAlignment(String queryFasta, String subjectFasta, Map<String,AbstractSequence<?>> querySequences, double treshold,  boolean transportersSearch, AtomicBoolean cancel, ConcurrentLinkedQueue<AlignmentCapsule> alignmentContainerSet, JAXBContext jc){

		this.queryFasta = queryFasta;
		this.subjectFasta = subjectFasta;
		this.threshold = treshold;
		this.isTransportersSearch = transportersSearch;
		this.querySequences = querySequences;
		this.alignmentContainerSet = alignmentContainerSet;
		this.cancel = cancel;
		this.jc = jc;

	}

	@Override
	public void run(){

		if(!this.cancel.get()) {

			try {

				File tcdbfile = new File(subjectFasta);

				String outputFileName = queryFasta.substring(queryFasta.lastIndexOf("/")).replace(".faa", "").concat("_blastReport.xml");
				if(isTransportersSearch)
					outputFileName = outputFileName.replace(".xml", "_transporters.xml");

				File outputFile = new File(tcdbfile.getParent().concat("\\..\\").concat("reports").concat(outputFileName));
				outputFile.getParentFile().mkdirs();
				
				System.out.println("queryFasta---->"+this.queryFasta);
				System.out.println("subjectFasta---->"+this.subjectFasta);
				System.out.println("outputFile.getAbsolutePath---->"+outputFile.getAbsolutePath());

				Process p = Runtime.getRuntime().exec("blastp -query " + this.queryFasta + " -subject " 
						+ this.subjectFasta + " -out " + outputFile.getAbsolutePath() + " -outfmt 5");

				p.waitFor();

				this.blout = new NcbiBlastParser(outputFile, this.jc);
				this.alignmentMatrix = blout.getMatrix();

				buildAlignmentCapsules();

			} catch (IOException | InterruptedException e) {

				e.printStackTrace();
			}
			catch (OutOfMemoryError oue) {

				oue.printStackTrace();
			}
			System.gc();

			setChanged();
			notifyObservers();
		}

		setChanged();
		notifyObservers();
	}


	public void buildAlignmentCapsules(){

		List<BlastIterationData> iterations = this.blout.getResults();

		Map<String, Double> queriesMaxScores = AlignmentsUtils.getSequencesAlignmentMaxScoreMap(querySequences, alignmentMatrix);

		for(BlastIterationData iteration : iterations){

			String queryID = iteration.getQueryDef().trim();
			Integer queryLength = iteration.getQueryLen();

			String [] query_array; 
			String query_org = "";
			String queryLocus = "";

			if(queryID.contains(":")) {
				query_array = queryID.split(":"); 
				query_org = query_array [0].trim();
				queryLocus = query_array[1].trim();
			}
			else if(queryID.contains(" "))
				queryID = new StringTokenizer(queryID," ").nextToken();
			

			if(this.blastPurpose==null || !this.blastPurpose.equals(AlignmentPurpose.ORTHOLOGS) || (!this.sequenceIdsSet.containsKey(queryLocus) || sequenceIdsSet.get(queryLocus).isEmpty())){

				double maxScore = queriesMaxScores.get(queryID);
				double specificThreshold = this.threshold;
				
				if(this.kegg_taxonomy_scores!=null && this.referenceTaxonomyScore!=null && this.referenceTaxonomyThreshold!=null)
					if(this.kegg_taxonomy_scores.get(query_org)>=this.referenceTaxonomyScore) 
						specificThreshold = this.referenceTaxonomyThreshold;

				List<Hit> hits = iteration.getHits();

				if(hits!=null && !hits.isEmpty()){

					for(Hit hit : hits){

						if(!this.cancel.get()){

							try {
								String tcdbID = "";
								String hitNum = hit.getHitNum();
								String target = hit.getHitId();
								
								Integer targetLength = iteration.getHitLength(hitNum);
								Integer alingmentLength = iteration.getHitAlignmentLength(hitNum);

								double alignmentScore = (iteration.getHitScore(hit)-ALIGNMENT_MIN_SCORE)/(maxScore-ALIGNMENT_MIN_SCORE);//alignmentMethod.getSimilarity(); //(((double)alignmentMethod.getScore()-alignmentMethod.getMinScore())/(alignmentMethod.getMaxScore()-alignmentMethod.getMinScore()))
								//double similarityScore = iteration.getPositivesScore(hitNum);
								//double identityScore = iteration.getIdentityScore(hitNum);

								double bitScore = iteration.getHitBitScore(hit);
								double eValue = iteration.getHitEvalue(hit);
								
								double queryCoverage = iteration.getHitQueryCoverage(hitNum);//(double)(alingmentLength-iteration.getHitAlignmentGaps(hitNum))/(double)queryLength;
								double tragetCoverage = iteration.getHiTargetCoverage(hitNum);//(double)(alingmentLength-iteration.getHitAlignmentGaps(hitNum))/(double)targetLength;

								double l1 = (double)queryLength/(double)targetLength;
//								double l2 = (double)alingmentLength/(double)queryLength;
//								double l3 = (double)alingmentLength/(double)targetLength;

								double score = alignmentScore;//-1;

								//				if(this.alignmentScoreType.equals(AlignmentScoreType.ALIGNMENT))
								//					score = alignmentScore;
								//				else if(this.alignmentScoreType.equals(AlignmentScoreType.IDENTITY))
								//					score = identityScore;
								//				else if(this.alignmentScoreType.equals(AlignmentScoreType.SIMILARITY))
								//					score = similarityScore;

								System.out.println(queryID+"\t"+target+"\t"+score+"\t"+specificThreshold+"\t"+iteration.getHitEvalue(hit)+"\t"+iteration.getHitBitScore(hit)
								+"\t"+l1);//)+"\t"+l2+"\t"+l3);

								boolean go = false;
								
								if(isTransportersSearch)
									if(eValue<FIXED_THRESHOLD && bitScore>BITSCORE_THRESHOLD && Math.abs(1-queryCoverage)<=COVERAGE_THRESHOLD)
										go=true;
								else if(blastPurpose.equals(AlignmentPurpose.ORTHOLOGS))
									if(score>specificThreshold)
										go=true;
									
								if(go){
									//									&& Math.abs(1-l1)<=ALIGNMENT_QUERY_LEN_THRESHOLD && Math.abs(1-l2)<=QUERY_HIT_LEN_THRESHOLD){

									if(this.sequencesWithoutSimilarities!=null && this.sequencesWithoutSimilarities.contains(queryID))
										this.sequencesWithoutSimilarities.remove(queryID);

									if(isTransportersSearch){

										String hitdef = hit.getHitDef();

										StringTokenizer st = new StringTokenizer(hitdef,"|");
										st.nextToken();
										st.nextToken();
										target = st.nextToken().toUpperCase().trim();
										tcdbID = st.nextToken().split("\\s+")[0].toUpperCase().trim();
									}

									AlignmentCapsule alignContainer = new AlignmentCapsule(queryID, target, tcdbID, this.alignmentMatrix, score);

									alignContainer.setEvalue(eValue);
									alignContainer.setBitScore(bitScore);
									alignContainer.setAlignmentLength(alingmentLength);
									alignContainer.setQueryLength(queryLength);
									alignContainer.setTargetLength(targetLength);	
									alignContainer.setNumIdenticals(iteration.getHitIdentity(hitNum));
									alignContainer.setNumSimilars(iteration.getHitPositive(hitNum));
									alignContainer.setCoverageQuery(queryCoverage);
									alignContainer.setCoverageTarget(tragetCoverage);

									alignContainer.setEcNumber(this.ec_number);
									alignContainer.setClosestOrthologues(this.closestOrthologs);
									alignContainer.setModules(modules);

									//					alignContainer.setMaxScore(maxScore);
									//					alignContainer.setMinScore(0);
									//					alignContainer.setAlignedScore(alignedScore);

									//					iterationAlignments.add(align);
									this.alignmentContainerSet.add(alignContainer);
								}
							} 
							catch (Exception e) {
								// TODO Auto-generated catch block
								e.printStackTrace();
							}
						}
					}

					//			this.alignments.put(queryID,iterationAlignments);
				}
				else{

					logger.warn(iteration.getIteration().getIterationMessage().concat(" for {}"), queryID);
				}
			}
			else{
				if(this.sequencesWithoutSimilarities!=null && this.sequencesWithoutSimilarities.contains(queryID))
					this.sequencesWithoutSimilarities.remove(queryID);
			}
		}
	}


	public ConcurrentLinkedQueue<AlignmentCapsule> getAlignmentsCapsules(){

		return this.alignmentContainerSet;
	}


	/**
	 * @return
	 */
	public Map<String,List<AlignmentCapsule>> getAlignmentsByQuery(){

		Map<String,List<AlignmentCapsule>> alignmentMap = new HashMap<>();

		for(AlignmentCapsule alignContainer : this.alignmentContainerSet){

			String query = alignContainer.getQuery();

			if(alignmentMap.containsKey(query)){

				alignmentMap.get(query).add(alignContainer);
			}
			else{
				List<AlignmentCapsule> containersList = new ArrayList<>();
				containersList.add(alignContainer);
				alignmentMap.put(query, containersList);
			}
		}

		return alignmentMap;
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
	 * @return the closestOrthologs
	 */
	public Map<String,Set<String>> getClosestOrthologs() {
		return closestOrthologs;
	}

	/**
	 * @param closestOrthologs the closestOrthologs to set
	 */
	public void setClosestOrthologs(Map<String,Set<String>> closestOrthologs) {
		this.closestOrthologs = closestOrthologs;
	}

	/**
	 * @return the modules
	 */
	public Map<String,Set<Integer>> getModules() {
		return modules;
	}

	/**
	 * @param modules the modules to set
	 */
	public void setModules(Map<String,Set<Integer>> modules) {
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
	public void setSequencesWithoutSimilarities(ConcurrentLinkedQueue<String> sequencesWithoutSimilarities) {
		this.sequencesWithoutSimilarities = sequencesWithoutSimilarities;
	}

	/**
	 * @return the blastPurpose
	 */
	public AlignmentPurpose getBlastPurpose() {
		return blastPurpose;
	}

	/**
	 * @param blastPurpose the blastPurpose to set
	 */
	public void setBlastPurpose(AlignmentPurpose blastPurpose) {
		this.blastPurpose = blastPurpose;
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

	/**
	 * @return the kegg_taxonomy_scores
	 */
	public Map<String, Integer> getKegg_taxonomy_scores() {
		return kegg_taxonomy_scores;
	}

	/**
	 * @param kegg_taxonomy_scores the kegg_taxonomy_scores to set
	 */
	public void setKegg_taxonomy_scores(Map<String, Integer> kegg_taxonomy_scores) {
		this.kegg_taxonomy_scores = kegg_taxonomy_scores;
	}

	/**
	 * @return the referenceTaxonomyScore
	 */
	public int getReferenceTaxonomyScore() {
		return referenceTaxonomyScore;
	}

	/**
	 * @param referenceTaxonomyScore the referenceTaxonomyScore to set
	 */
	public void setReferenceTaxonomyScore(int referenceTaxonomyScore) {
		this.referenceTaxonomyScore = referenceTaxonomyScore;
	}

	/**
	 * @return the referenceTaxonomyThreshold
	 */
	public double getReferenceTaxonomyThreshold() {
		return referenceTaxonomyThreshold;
	}

	/**
	 * @param referenceTaxonomyThreshold the referenceTaxonomyThreshold to set
	 */
	public void setReferenceTaxonomyThreshold(double referenceTaxonomyThreshold) {
		this.referenceTaxonomyThreshold = referenceTaxonomyThreshold;
	}

}
