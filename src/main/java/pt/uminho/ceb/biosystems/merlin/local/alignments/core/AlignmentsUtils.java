package pt.uminho.ceb.biosystems.merlin.local.alignments.core;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentLinkedQueue;

import org.biojava.nbio.core.alignment.matrices.SubstitutionMatrixHelper;
import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.template.AbstractSequence;

import pt.uminho.ceb.biosystems.merlin.utilities.containers.capsules.AlignmentCapsule;

public class AlignmentsUtils {
	
	
	/**
	 * @param sequence
	 * @return
	 */
	public static double getSequenceAlignmentMaxScore(AbstractSequence<AminoAcidCompound> sequence, String matrixName){
		
		double maxScore = 0;
		
		SubstitutionMatrix<AminoAcidCompound> matrix = SubstitutionMatrixHelper.getAminoAcidSubstitutionMatrix(matrixName.toLowerCase());
		
		List<AminoAcidCompound> sequenceCompounds = sequence.getAsList();
		
		for(AminoAcidCompound compound : sequenceCompounds){
			
			maxScore += (double) matrix.getValue(compound, compound);
		}
		
		return maxScore;
	}
	
	
	/**
	 * @param sequences
	 * @return
	 */
	public static Map<String,Double> getSequencesAlignmentMaxScoreMap(Map<String,AbstractSequence<?>> sequences, String matrixName){
		
		Map<String,Double> maxScores = new HashMap<>();
		
		for(String sequenceID : sequences.keySet()){
			
			@SuppressWarnings("unchecked")
			AbstractSequence<AminoAcidCompound> sequence = (AbstractSequence<AminoAcidCompound>) sequences.get(sequenceID);
			
			double sequenceMaxScore = getSequenceAlignmentMaxScore(sequence,matrixName);
			
			maxScores.put(sequenceID, sequenceMaxScore);
		}
		
		return maxScores;
		
	}
	
	
	/**
	 * Method to verify if BLAST+ is installed
	 * 
	 * @return
	 */
	public static boolean checkBlastInstalation() {
		
		try {
			Process p = Runtime.getRuntime().exec("blastp -version");
			p.waitFor();
			return true;
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			if(e.getMessage().contains("The system cannot find the file specified")){
				return false;
			}
			else {
				e.printStackTrace();
				return false;
			}
		} catch (InterruptedException e) {
			e.printStackTrace();
			return false;
		}
	}
	
	
	/**
	 * @param alignmentContainerSet
	 * @return
	 */
	public static Map<String,List<AlignmentCapsule>> getAlignmentsByQuery(ConcurrentLinkedQueue<AlignmentCapsule> alignmentContainerSet){

		Map<String,List<AlignmentCapsule>> alignmentMap = new HashMap<>();

		for(AlignmentCapsule alignContainer : alignmentContainerSet){

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
	 * @param alignmentContainerSet
	 * @return
	 */
	public static Map<String,String> getOrthologsGenesMap(ConcurrentLinkedQueue<AlignmentCapsule> alignmentContainerSet){
		
		Map<String,List<AlignmentCapsule>> alignmentsMap = getAlignmentsByQuery(alignmentContainerSet);
		
		Map<String,String> orthologsGenesMap = new HashMap<>();
		
		
		
		
		return orthologsGenesMap;
		
	}
	
	
}
