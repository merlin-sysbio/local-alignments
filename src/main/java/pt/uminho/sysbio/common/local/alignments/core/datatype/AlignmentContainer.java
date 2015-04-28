package pt.uminho.sysbio.common.local.alignments.core.datatype;



/**
 * @author Oscar
 *
 */
public class AlignmentContainer {
	
	private String query;
	private String target;
	private double alignmentScore;
	private double similarityScore;
	private double identityScore;
	private String matrix;
	private String ko;
	private int alignmentLength;
	private double identityLength;
	private double similarityLength;
	private double queryLength, targetLength;
	private double score, maxScore, minScore;
	private short[][][] scoreMatrix;
	
	
	/**
	 * @param query
	 * @param target
	 * @param alignmentScore
	 * @param similarityScore
	 * @param identityScore
	 * @param matrix
	 * @param ko
	 * @param alignmentLength
	 * @param identityLength
	 * @param similarityLength
	 * @param queryLength
	 * @param targetLength
	 * @param score
	 * @param maxScore
	 * @param minScore
	 */
	public AlignmentContainer(String query, String target,
			double alignmentScore, double similarityScore,
			double identityScore, String matrix, String ko,
			int alignmentLength, double identityLength,
			double similarityLength, double queryLength, 
			double targetLength, double score, 
			double maxScore, double minScore) {
		super();
		this.query = query;
		this.target = target;
		this.alignmentScore = alignmentScore;
		this.similarityScore = similarityScore;
		this.identityScore = identityScore;
		this.matrix = matrix;
		this.ko = ko;
		this.setAlignmentLength(alignmentLength);
		this.setIdentityLength(identityLength);
		this.setSimilarityLength(similarityLength);
		this.setQueryLength(queryLength);
		this.setTargetLength(targetLength);
		this.setScore(score);
		this.setMaxScore(maxScore);
		this.setMinScore(minScore);
	}
	
	/**
	 * @param query
	 * @param target
	 * @param alignmentScore
	 * @param similarityScore
	 * @param identityScore
	 * @param matrix
	 * @param ko
	 */
	public AlignmentContainer(String query, String target,
			double alignmentScore, double similarityScore,
			double identityScore, String matrix, String ko) {
		super();
		this.query = query;
		this.target = target;
		this.alignmentScore = alignmentScore;
		this.setSimilarityScore(similarityScore);
		this.setIdentityScore(identityScore);
		this.matrix = matrix;
		this.ko = ko;
	}

	/**
	 * @param query
	 * @param locusTag
	 * @param alignmentScore
	 * @param similarityScore
	 * @param identityScore
	 * @param matrix
	 */
	public AlignmentContainer(String query, String locusTag,
			double alignmentScore, double similarityScore,
			double identityScore, String matrix) {
		super();
		this.query = query;
		this.target = locusTag;
		this.alignmentScore = alignmentScore;
		this.setSimilarityScore(similarityScore);
		this.setIdentityScore(identityScore);
		this.matrix = matrix;
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
	 * @return the alignmentScore
	 */
	public double getAlignmentScore() {
		return alignmentScore;
	}

	/**
	 * @param alignmentScore the alignmentScore to set
	 */
	public void setAlignmentScore(double alignmentScore) {
		this.alignmentScore = alignmentScore;
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
	 * @return the similarityScore
	 */
	public double getSimilarityScore() {
		return similarityScore;
	}



	/**
	 * @param similarityScore the similarityScore to set
	 */
	public void setSimilarityScore(double similarityScore) {
		this.similarityScore = similarityScore;
	}



	/**
	 * @return the identityScore
	 */
	public double getIdentityScore() {
		return identityScore;
	}



	/**
	 * @param identityScore the identityScore to set
	 */
	public void setIdentityScore(double identityScore) {
		this.identityScore = identityScore;
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
	 * @return the identityLength
	 */
	public double getIdentityLength() {
		return identityLength;
	}

	/**
	 * @param identityLength the identityLength to set
	 */
	public void setIdentityLength(double identityLength) {
		this.identityLength = identityLength;
	}

	/**
	 * @return the similarityLength
	 */
	public double getSimilarityLength() {
		return similarityLength;
	}

	/**
	 * @param similarityLength the similarityLength to set
	 */
	public void setSimilarityLength(double similarityLength) {
		this.similarityLength = similarityLength;
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

	@Override
	public String toString() {
		return "AlignmentContainer [query=" + query + ", locusTag=" + target
				+ ", alignmentScore=" + alignmentScore + ", similarityScore="
				+ similarityScore + ", identityScore=" + identityScore
				+ ", matrix=" + matrix + ", ko=" + ko + ", alignmentLength="
				+ alignmentLength + ", identityLength=" + identityLength
				+ ", similarityLength=" + similarityLength + ", queryLenght="
				+ queryLength + ", targetLength=" + targetLength + ", score="
				+ score + ", maxScore=" + maxScore + ", minScore=" + minScore
				+ "]";
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
	public void setScoreMatrix(short[][][] scoreMatrix) {
		
		this.scoreMatrix = scoreMatrix;
	}
	
	/**
	 * @return
	 */
	public short[][][] getScoreMatrix() {
		
		return this.scoreMatrix;
	}

}
