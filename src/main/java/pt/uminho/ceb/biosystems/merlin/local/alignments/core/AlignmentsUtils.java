package pt.uminho.ceb.biosystems.merlin.local.alignments.core;

import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.biojava.nbio.core.alignment.matrices.SubstitutionMatrixHelper;
import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.template.AbstractSequence;

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
			if(e.getMessage().contains("The system cannot find the file specified")) {
				
				System.out.println("BLAST NOT found\n"+e.getMessage());
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
	
	
}
