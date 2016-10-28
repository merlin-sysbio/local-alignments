package pt.uminho.sysbio.common.local.alignments.core.datatype;

import java.util.Arrays;

import pt.uminho.sysbio.common.local.alignments.core.Enumerators.AlignmentScoreType;
import pt.uminho.sysbio.common.local.alignments.core.Enumerators.Method;

/**
 * @author Oscar
 *
 */
public class AlignmentContainer {

	private String query, queryLocus;
	private String target, targetLocus;
	private AlignmentScoreType alignmentScoreType;
	private String matrix;
	private String ko;
	private double score, alignedScore, maxScore, minScore;
	private Method method;
	private int[][][] scoreMatrix;
	private int numIdenticals, numSimilars, queryLength, targetLength, alignmentLength;

	/**
	 * Alignment information container.
	 * 
	 * @param query
	 * @param target
	 * @param ko
	 * @param score
	 * @param maxScore
	 * @param minScore
	 * @param alignedScore
	 * @param numIdenticals
	 * @param numSimilars
	 * @param alignmentLength
	 * @param queryLength
	 * @param targetLength
	 * @param matrix
	 * @param method
	 * @param alignmentScoreType
	 */
	public AlignmentContainer(String query, String target, String ko, double score, double maxScore, double minScore, double alignedScore,
			int numIdenticals, int numSimilars, int alignmentLength, int queryLength, int targetLength, 
			String matrix, Method method, AlignmentScoreType alignmentScoreType) {

		super();
		this.setQuery(query);
		this.setTarget(target);
		this.setKo(ko);
		this.setScore(score);
		this.setMaxScore(maxScore);
		this.setMinScore(minScore);
		this.setAlignedScore(alignedScore);
		this.setNumIdenticals(numIdenticals);
		this.setNumSimilars(numSimilars);
		this.setAlignmentLength(alignmentLength);
		this.setQueryLength(queryLength);
		this.setTargetLength(targetLength);
		this.setMethod(method);
		this.setMatrix(matrix);
		this.setAlignmentScoreType(alignmentScoreType);
	}

	/**
	 * @return the alignmentScoreType
	 */
	public AlignmentScoreType getAlignmentScoreType() {
		return alignmentScoreType;
	}

	/**
	 * @param alignmentScoreType the alignmentScoreType to set
	 */
	public void setAlignmentScoreType(AlignmentScoreType alignmentScoreType) {
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
	public int getQueryLength() {
		return queryLength;
	}

	/**
	 * @param queryLenght the queryLenght to set
	 */
	public void setQueryLength(int queryLength) {
		this.queryLength = queryLength;
	}

	/**
	 * @return the targetLength
	 */
	public int getTargetLength() {
		return targetLength;
	}

	/**
	 * @param targetLength the targetLength to set
	 */
	public void setTargetLength(int targetLength) {
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

	/**
	 * @return the queryLocus
	 */
	public String getQueryLocus() {
		return queryLocus;
	}

	/**
	 * @param queryLocus the queryLocus to set
	 */
	public void setQueryLocus(String queryLocus) {
		this.queryLocus = queryLocus;
	}

	/**
	 * @return the targetLocus
	 */
	public String getTargetLocus() {
		return targetLocus;
	}

	/**
	 * @param targetLocus the targetLocus to set
	 */
	public void setTargetLocus(String targetLocus) {
		this.targetLocus = targetLocus;
	}

	/**
	 * @return the numIdenticals
	 */
	public int getNumIdenticals() {
		return numIdenticals;
	}

	/**
	 * @param numIdenticals the numIdenticals to set
	 */
	public void setNumIdenticals(int numIdenticals) {
		this.numIdenticals = numIdenticals;
	}

	/**
	 * @return the numSimilars
	 */
	public int getNumSimilars() {
		return numSimilars;
	}

	/**
	 * @param numSimilars the numSimilars to set
	 */
	public void setNumSimilars(int numSimilars) {
		this.numSimilars = numSimilars;
	}

	/**
	 * @return the alignmentLength
	 */
	public int getAlignmentLength() {
		return alignmentLength;
	}

	/**
	 * @param alignmentLength the alignmentLength to set
	 */
	public void setAlignmentLength(int alignmentLength) {
		this.alignmentLength = alignmentLength;
	}

	/**
	 * @return the alignedScore
	 */
	public double getAlignedScore() {
		return alignedScore;
	}

	/**
	 * @param alignedScore the alignedScore to set
	 */
	public void setAlignedScore(double alignedScore) {
		this.alignedScore = alignedScore;
	}

	/**
	 * Get overall alignment score.
	 * 
	 * @return double alignment score
	 */
	public double getAlignmentScore(){

		return (this.getAlignedScore()-this.getMinScore())/(this.getMaxScore()-this.getMinScore());
	}

	/**
	 * Get overall similarity score.
	 * 
	 * @return double similarity score
	 */
	public double getSimilarityScore(){

		return ((double)this.getNumSimilars()/(double)this.getAlignmentLength());
	}

	/**
	 * Get overall identity score.
	 * 
	 * @return double identity score
	 */
	public double getIdentityScore(){

		return ((double)this.getNumIdenticals()/(double)this.getAlignmentLength());
	}

	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	@Override
	public String toString() {
		return "AlignmentContainer [query=" + query + ", queryLocus=" + queryLocus + ", target=" + target
				+ ", targetLocus=" + targetLocus + ", alignmentScoreType=" + alignmentScoreType + ", matrix=" + matrix
				+ ", ko=" + ko + ", score=" + score + ", alignedScore=" + alignedScore + ", maxScore=" + maxScore
				+ ", minScore=" + minScore + ", method=" + method + ", scoreMatrix=" + Arrays.toString(scoreMatrix)
				+ ", numIdenticals=" + numIdenticals + ", numSimilars=" + numSimilars + ", queryLength=" + queryLength
				+ ", targetLength=" + targetLength + ", alignmentLength=" + alignmentLength  
				+ ", getIdentityScore=" + this.getIdentityScore() + ", getNumSimilars=" + this.getSimilarityScore() 
				+ ", getAlignmentScore=" + this.getAlignmentScore() 
				+ "]";
	}
	
}
