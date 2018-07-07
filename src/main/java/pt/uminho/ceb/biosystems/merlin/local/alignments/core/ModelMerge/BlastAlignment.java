package pt.uminho.ceb.biosystems.merlin.local.alignments.core.ModelMerge;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Observable;
import java.util.StringTokenizer;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;

import javax.xml.bind.JAXBContext;

import org.apache.commons.io.output.ThresholdingOutputStream;
import org.apache.commons.math3.geometry.euclidean.oned.Interval;
import org.apache.commons.math3.geometry.partitioning.Region.Location;
import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.core.sequence.template.AbstractSequence;

import pt.uminho.ceb.biosystems.merlin.local.alignments.core.AlignmentsUtils;
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
	private static final Interval ALIGNMENT_QUERY_LEN_THRESHOLD = new Interval(0.75, 1.25);
	private static final Interval ALIGNMENT_HIT_LEN_THRESHOLD = new Interval(0.75, 1.25);
	
	private NcbiBlastParser blout;
	private ConcurrentLinkedQueue<AlignmentCapsule> alignmentContainerSet;
	private String alignmentMatrix, queryFasta, subjectFasta;
	private boolean isTransportersSearch = false;
	private double threshold;
	private AtomicBoolean cancel; 
	private Map<String,AbstractSequence<?>> querySequences;
	private JAXBContext jc;

	
//	public BlastAlignment(String queryFasta, String subjectFasta, double threshold, AtomicBoolean cancel){
//		
//		this.queryFasta = queryFasta;
//		this.subjectFasta = subjectFasta;
//		this.threshold = threshold;
//		this.cancel = cancel;
//	}
	
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
				
				Process p = Runtime.getRuntime().exec("blastp -query " + this.queryFasta + " -subject " 
						+ this.subjectFasta + " -out " + outputFile.getAbsolutePath() + " -outfmt 5");

				p.waitFor();
				
				this.blout = new NcbiBlastParser(outputFile, this.jc);

				this.alignmentMatrix = blout.getMatrix();
				
				System.out.println("Building alignment capsules....");
				
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
			
			System.out.println(iteration.getQueryID());
			System.out.println(iteration.getQueryDef());
			System.out.println("File path: "+this.queryFasta.substring(this.queryFasta.lastIndexOf("/")));
			
			System.out.println("queryID"+"\t"+"target"+"\t"+"Score"+"\t"+"Threshold"+"\t"+"Evalue(hit)"+"\t"+"BitScore(hit)");


			String queryID = iteration.getQueryDef().trim();
			
			if(queryID.contains(" "))
				queryID = new StringTokenizer(queryID," ").nextToken();
			
			double maxScore = queriesMaxScores.get(queryID);
			
//			List<AlignmentCapsule> iterationAlignments = new ArrayList<>();
			List<Hit> hits = iteration.getHits();

			for(Hit hit : hits){
				
				if(!this.cancel.get()){
					
					try {
						String tcdbID = "";
						String hitNum = hit.getHitNum();
						String target = hit.getHitId();
						//				String targetLength = hit.getHitLen();
						Integer alingmentLength = iteration.getHitAlignmentLength(hitNum);


						double alignmentScore = (iteration.getHitScore(hit)-ALIGNMENT_MIN_SCORE)/(maxScore-ALIGNMENT_MIN_SCORE);//alignmentMethod.getSimilarity(); //(((double)alignmentMethod.getScore()-alignmentMethod.getMinScore())/(alignmentMethod.getMaxScore()-alignmentMethod.getMinScore()))
						//double similarityScore = iteration.getPositivesScore(hitNum);
						//double identityScore = iteration.getIdentityScore(hitNum);
						
						double bitScore = iteration.getHitBitScore(hit);
						double eValue = iteration.getHitEvalue(hit);

						double score = -1;

						if(isTransportersSearch){
							score = alignmentScore;
						}
						else{
							//				if(this.alignmentScoreType.equals(AlignmentScoreType.ALIGNMENT))
							//					score = alignmentScore;
							//				else if(this.alignmentScoreType.equals(AlignmentScoreType.IDENTITY))
							//					score = identityScore;
							//				else if(this.alignmentScoreType.equals(AlignmentScoreType.SIMILARITY))
							//					score = similarityScore;
						}
						
						
						double l1 = (double)alingmentLength/(double)iteration.getQueryLen();
						double l2 = (double)alingmentLength/(double)iteration.getHitLength(hitNum);

						double specificThreshold = this.threshold;
//						if (this.querySpecificThreshold.containsKey(queryID))
//							specificThreshold = this.querySpecificThreshold.get(queryID);
						
						System.out.println(queryID+"\t"+target+"\t"+score+"\t"+specificThreshold+"\t"+iteration.getHitEvalue(hit)+"\t"+iteration.getHitBitScore(hit)
							+"\t"+l1+"\t"+l2);

						
						if(eValue<FIXED_THRESHOLD && bitScore>BITSCORE_THRESHOLD && ALIGNMENT_QUERY_LEN_THRESHOLD.checkPoint(l1, 0).equals(Location.INSIDE) 
								&& ALIGNMENT_HIT_LEN_THRESHOLD.checkPoint(l2, 0).equals(Location.INSIDE)){
//						if(score>specificThreshold) {

							if(isTransportersSearch){
								
								String hitdef = hit.getHitDef();

								StringTokenizer st = new StringTokenizer(hitdef,"|");
								st.nextToken();
								st.nextToken();
								target = st.nextToken().toUpperCase().trim();
								tcdbID = st.nextToken().split("\\s+")[0].toUpperCase().trim();
							}
							
							AlignmentCapsule alignContainer = new AlignmentCapsule(queryID, target, tcdbID, this.alignmentMatrix, score);

							alignContainer.setEvalue(iteration.getHitEvalue(hit));
							alignContainer.setBitScore(iteration.getHitBitScore(hit));
							alignContainer.setAlignmentLength((Integer)iteration.getHitAlignmentLength(hitNum));
							alignContainer.setQueryLength(iteration.getQueryLen());
							alignContainer.setTargetLength(iteration.getHitLength(hitNum));	

							alignContainer.setNumIdenticals(iteration.getHitIdentity(hitNum));
							alignContainer.setNumSimilars(iteration.getHitPositive(hitNum));

							//					alignContainer.setMaxScore(maxScore);
							//					alignContainer.setMinScore(0);
							//					alignContainer.setAlignedScore(alignedScore);

							//					alignContainer.setGapsQuery(gapsQuery);
							//					alignContainer.setGapsTarget(gapsTarget); 

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


	
}