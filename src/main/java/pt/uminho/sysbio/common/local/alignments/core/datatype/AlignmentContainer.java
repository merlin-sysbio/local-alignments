package pt.uminho.sysbio.common.local.alignments.core.datatype;

import java.util.Arrays;

import org.biojava.nbio.alignment.Alignments.PairwiseSequenceScorerType;

import pt.uminho.sysbio.common.local.alignments.core.Enumerators.Method;

/**
 * @author Oscar
 *
 */
public class AlignmentContainer {
	
	private String query;
	private String target;
	private PairwiseSequenceScorerType alignmentScoreType;
	private String matrix;
	private String ko;
	private double queryLength, targetLength;
	private double score, maxScore, minScore;
	private Method method;
	private int[][][] scoreMatrix;
	
	
	/**
	 * @param alignmentScoreType
	 * @param query
	 * @param target
	 * @param score
	 * @param ko
	 * @param queryLength
	 * @param targetLength
	 * @param maxScore
	 * @param minScore
	 * @param matrix
	 * @param method
	 */
	public AlignmentContainer(PairwiseSequenceScorerType alignmentScoreType, String query, String target, double score, String ko, int queryLength, int targetLength, double maxScore, double minScore, String matrix, Method method) {
		
		super();
		this.alignmentScoreType = alignmentScoreType;
		this.query = query;
		this.target = target;
		this.matrix = matrix;
		this.ko = ko;
		this.setQueryLength(queryLength);
		this.setTargetLength(targetLength);
		this.setScore(score);
		this.setMaxScore(maxScore);
		this.setMinScore(minScore);
		this.setMethod(method);
	}
	
	/**
	 * @param alignmentScoreType
	 * @param query
	 * @param target
	 * @param score
	 * @param ko
	 * @param queryLength
	 * @param targetLength
	 * @param matrix
	 * @param method
	 */
	public AlignmentContainer(PairwiseSequenceScorerType alignmentScoreType, String query, String target, double score, String ko, int queryLength, int targetLength, String matrix, Method method) {
		
		super();
		this.alignmentScoreType = alignmentScoreType;
		this.query = query;
		this.target = target;
		this.ko = ko;
		this.setQueryLength(queryLength);
		this.setTargetLength(targetLength);
		this.matrix = matrix;
		this.setScore(score);
		this.setMethod(method);
	}
	
	/**
	 * @param alignmentScoreType
	 * @param query
	 * @param target
	 * @param score
	 * @param queryLength
	 * @param targetLength
	 * @param matrix
	 * @param method
	 */
	public AlignmentContainer(PairwiseSequenceScorerType alignmentScoreType, String query, String target, double score, int queryLength, int targetLength, String matrix, Method method) {
		
		super();
		this.alignmentScoreType = alignmentScoreType;
		this.query = query;
		this.target = target;
		this.setQueryLength(queryLength);
		this.setTargetLength(targetLength);
		this.matrix = matrix;
		this.setScore(score);
	}
	
	
	/**
	 * @param alignmentScoreType
	 * @param query
	 * @param target
	 * @param score
	 * @param ko
	 * @param matrix
	 */
	public AlignmentContainer(PairwiseSequenceScorerType alignmentScoreType, String query, String target, double score, String ko, String matrix, Method method) {
		super();
		this.alignmentScoreType = alignmentScoreType;
		this.query = query;
		this.target = target;
		this.score = score;
		this.ko = ko;
		this.matrix = matrix;
		this.setMethod(method);
	}

	/**
	 * @param alignmentScoreType
	 * @param query
	 * @param locusTag
	 * @param score
	 * @param identityScore
	 * @param matrix
	 */
	/**
	 * @param alignmentScoreType
	 * @param query
	 * @param locusTag
	 * @param score
	 * @param identityScore
	 * @param matrix
	 * @param method
	 */
	public AlignmentContainer(PairwiseSequenceScorerType alignmentScoreType, String query, String locusTag, double score, double identityScore, String matrix, Method method) {
		
		super();
		this.alignmentScoreType = alignmentScoreType;
		this.query = query;
		this.target = locusTag;
		this.score = score;
		this.matrix = matrix;
		this.setMethod(method);
	}

	/**
	 * @return the alignmentScoreType
	 */
	public PairwiseSequenceScorerType getAlignmentScoreType() {
		return alignmentScoreType;
	}

	/**
	 * @param alignmentScoreType the alignmentScoreType to set
	 */
	public void setAlignmentScoreType(PairwiseSequenceScorerType alignmentScoreType) {
		this.alignmentScoreType = alignmentScoreType;
	}

	/**
	 * @return the query
	 */
	public String getQuery() {
		return query;
	}

	/**
	 * @param query the query to set
	 */
	public void setQuery(String query) {
		this.query = query;
	}

	/**
	 * @return the locusTag
	 */
	public String getTarget() {
		return target;
	}

	/**
	 * @param locusTag the locusTag to set
	 */
	public void setTarget(String target) {
		this.target = target;
	}


	/**
	 * @return the matrix
	 */
	public String getMatrix() {
		return matrix;
	}

	/**
	 * @param matrix the matrix to set
	 */
	public void setMatrix(String matrix) {
		this.matrix = matrix;
	}

	/**
	 * @return the ko
	 */
	public String getKo() {
		return ko;
	}



	/**
	 * @param ko the ko to set
	 */
	public void setKo(String ko) {
		this.ko = ko;
	}


	/**
	 * @return the score
	 */
	public double getScore() {
		return score;
	}

	/**
	 * @param score the score to set
	 */
	public void setScore(double score) {
		this.score = score;
	}

	/**
	 * @return the maxScore
	 */
	public double getMaxScore() {
		return maxScore;
	}

	/**
	 * @param maxScore the maxScore to set
	 */
	public void setMaxScore(double maxScore) {
		this.maxScore = maxScore;
	}

	/**
	 * @return the minscore
	 */
	public double getMinScore() {
		return minScore;
	}

	/**
	 * @param minscore the minscore to set
	 */
	public void setMinScore(double minscore) {
		this.minScore = minscore;
	}

	/**
	 * @return the queryLenght
	 */
	public double getQueryLength() {
		return queryLength;
	}

	/**
	 * @param queryLenght the queryLenght to set
	 */
	public void setQueryLength(double queryLength) {
		this.queryLength = queryLength;
	}

	/**
	 * @return the targetLength
	 */
	public double getTargetLength() {
		return targetLength;
	}

	/**
	 * @param targetLength the targetLength to set
	 */
	public void setTargetLength(double targetLength) {
		this.targetLength = targetLength;
	}

	/**
	 * @param scoreMatrix
	 */
	public void setScoreMatrix(int[][][] scoreMatrix) {
		
		this.scoreMatrix = scoreMatrix;
	}
	
	/**
	 * @return
	 */
	public int[][][] getScoreMatrix() {
		
		return this.scoreMatrix;
	}

	/**
	 * @return the method
	 */
	public Method getMethod() {
		return method;
	}

	/**
	 * @param method the method to set
	 */
	public void setMethod(Method method) {
		this.method = method;
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return "AlignmentContainer [query=" + query + ", target=" + target + ", alignmentScoreType="
				+ alignmentScoreType + ", matrix=" + matrix + ", ko=" + ko
				+ ", queryLength=" + queryLength + ", targetLength=" + targetLength + ", score=" + score + ", maxScore="
				+ maxScore + ", minScore=" + minScore + ", scoreMatrix=" + Arrays.toString(scoreMatrix) + "]";
	}

	
	
}
