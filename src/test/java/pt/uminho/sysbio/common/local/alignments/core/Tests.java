package pt.uminho.sysbio.common.local.alignments.core;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Reader;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava.nbio.alignment.Alignments.PairwiseSequenceScorerType;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.alignment.SmithWaterman;
import org.biojava.nbio.alignment.template.AbstractPairwiseSequenceAligner;
import org.biojava.nbio.alignment.template.GapPenalty;
import org.biojava.nbio.alignment.template.PairwiseSequenceAligner;
import org.biojava.nbio.alignment.template.PairwiseSequenceScorer;
import org.biojava.nbio.core.alignment.matrices.SimpleSubstitutionMatrix;
import org.biojava.nbio.core.alignment.template.AlignedSequence;
import org.biojava.nbio.core.alignment.template.SequencePair;
import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.sequence.ProteinSequence;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompound;
import org.biojava.nbio.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava.nbio.core.sequence.io.FastaReaderHelper;
import org.biojava.nbio.core.sequence.template.CompoundSet;
import org.biojava.nbio.core.sequence.template.Sequence;
import org.junit.Test;

import pt.uminho.sysbio.common.local.alignments.core.Enumerators.Matrix;

public class Tests {


	//@Test
	public void align() throws Exception {


		System.out.println("######################### old sw align #######################");
		String ko = null;
		String uni_id = "A1TX06";//"Q0VNJ6"; // ABO Q0VNJ6

		ProteinSequence	querySequence = FastaReaderHelper.readFastaProteinSequence(new URL(String.format("http://www.uniprot.org/uniprot/%s.fasta", uni_id)).openStream()).get(uni_id);

		Set<String> orthologs = new HashSet<>();
		//orthologs.add("C1B498");
		//		orthologs.add("Q93ZR6"); // ARA TH Q93ZR6
		//		orthologs.add("C1ASK8");
		//		orthologs.add("Q14693");
		//		orthologs.add("O14494");
		//		orthologs.add("A4IWB4");
		//		orthologs.add("D2A465");
		//		orthologs.add("I1G979");
		//		orthologs.add("Q0BKE1");
		//		orthologs.add("F7B3Z1");
		//		orthologs.add("M3X9V6");
		//		orthologs.add("B8XSI9");
		//		orthologs.add("B8XSI7");
		orthologs.add("A1U572");

		for (String uniprot_id : orthologs) {

			ProteinSequence	genomeAAsequence = FastaReaderHelper.readFastaProteinSequence(new URL(String.format("http://www.uniprot.org/uniprot/%s.fasta", uniprot_id)).openStream()).get(uniprot_id);

			Matrix matrix = Matrix.BLOSUM62;
			short gapOpenPenalty=10, gapExtensionPenalty=1;

			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

			//		SubstitutionMatrix<AminoAcidCompound> matrixS = new SimpleSubstitutionMatrix<AminoAcidCompound>();
			//		SequencePair<ProteinSequence, AminoAcidCompound> pair = Alignments.getPairwiseAlignment(querySequence, genomeAAsequence,
			//				PairwiseSequenceAlignerType.LOCAL, new SimpleGapPenalty(), matrixS);
			//		System.out.printf("%n%s vs %s%n%s", pair.getQuery().getAccession(), pair.getTarget().getAccession(), pair, pair.getNumIdenticals());

			//		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


			GapPenalty gp = new SimpleGapPenalty(gapOpenPenalty ,gapExtensionPenalty);
			CompoundSet<AminoAcidCompound> aa = new AminoAcidCompoundSet();
			Reader rd = new InputStreamReader(getClass().getClassLoader().getResourceAsStream(matrix.getPath()));
			SubstitutionMatrix<AminoAcidCompound> sb = new SimpleSubstitutionMatrix<AminoAcidCompound>(aa, rd, matrix.getPath());

			AbstractPairwiseSequenceAligner<Sequence<AminoAcidCompound>,AminoAcidCompound> alignmentMethod 
			= new SmithWaterman<Sequence<AminoAcidCompound>,AminoAcidCompound>((Sequence<AminoAcidCompound>)querySequence, (Sequence<AminoAcidCompound>)genomeAAsequence, gp, sb);

			double	alignmentScore = alignmentMethod.getSimilarity();
			double	similarityScore = ((double)alignmentMethod.getPair().getNumSimilars()/alignmentMethod.getPair().getLength());
			double	identityScore = ((double)alignmentMethod.getPair().getNumIdenticals()/alignmentMethod.getPair().getLength());

			//			System.out.println( alignmentMethod.getSimilarity());
			//			System.out.println(alignmentMethod.getMaxScore());
			//			System.out.println(alignmentMethod.getMinScore());
			//			System.out.println(alignmentMethod.getDistance());
			System.out.println("Alignment Score "+alignmentMethod.getScore());
			System.out.println("Min Score "+alignmentMethod.getMinScore());
			System.out.println("Max Score "+alignmentMethod.getMaxScore());
			System.out.println("Score - MinScore / MaxScore - MinScore "+(((double)alignmentMethod.getScore()-alignmentMethod.getMinScore())/(alignmentMethod.getMaxScore()-alignmentMethod.getMinScore()))+"\n");
			System.out.println("NumSimilars "+alignmentMethod.getPair().getNumSimilars());
			System.out.println("Lenght "+alignmentMethod.getPair().getLength());
			System.out.println("Similarity Score (NumSimilars/Length) "+similarityScore+"\n");
			System.out.println("NumIdenticals "+alignmentMethod.getPair().getNumIdenticals());
			System.out.println("Lenght "+alignmentMethod.getPair().getLength());
			System.out.println("Identity Score (NumIdenticals/Length) "+identityScore+"\n");
			System.out.println("Matrix ");
			FileWriter fileWriter = null;

			fileWriter = new FileWriter("C:/Users/Oscar/Desktop/scoreMatrix.txt");
			fileWriter.write(alignmentMethod.getScoreMatrixAsString().replaceAll(" +", "\t").replace("0", ""));
			fileWriter.close();
			System.out.println();
			//			for(short[][] saa : alignmentMethod.getScoreMatrix()) {
			//				for(short[] sa : saa) {
			//					for(short s : sa) {
			//
			//						System.out.print(s+" ");
			//					}
			//					System.out.println();
			//				}
			//				System.out.println();
			//			}

			//AlignmentContainer alignmentContainer = new AlignmentContainer(PairwiseSequenceScorerType.LOCAL, uni_id, uniprot_id, alignmentScore, matrix.toString(), ko, Method.SmithWaterman);
			//System.out.println(alignmentContainer);

		}
	}
	
	@Test
	public void newAlignBetter() throws MalformedURLException, IOException, Exception{

		System.out.println("######################### new sw align #######################");

		//		String uni_id1 = "Q0VKV8";	// "Q0VNJ6";
		//		String uni_id2 = "C1B498";	// "Q93ZR6";
		//		String uni_id1 = "A1TX06";
		//		String uni_id2 = ("A1U572");
				String uni_id1 = "P49374";
				String uni_id2 = ("P0AGF4");
				ProteinSequence s1 = FastaReaderHelper.readFastaProteinSequence(new URL(String.format("http://www.uniprot.org/uniprot/%s.fasta", uni_id1)).openStream()).get(uni_id1);
				ProteinSequence s2 = FastaReaderHelper.readFastaProteinSequence(new URL(String.format("http://www.uniprot.org/uniprot/%s.fasta", uni_id2)).openStream()).get(uni_id2);

//		ProteinSequence s1 = FastaReaderHelper.readFastaProteinSequence(new FileInputStream(new File("C:/Users/Oscar Dias/Desktop/seq1.txt"))).get("seq1");
//		ProteinSequence s2 = FastaReaderHelper.readFastaProteinSequence(new FileInputStream(new File("C:/Users/Oscar Dias/Desktop/seq2.txt"))).get("seq2");

		Matrix matrix = Matrix.BLOSUM62;
		short gapOpenPenalty=10, gapExtensionPenalty=1;

		GapPenalty gp = new SimpleGapPenalty(gapOpenPenalty, gapExtensionPenalty);
		CompoundSet<AminoAcidCompound> aa = new AminoAcidCompoundSet();
		Reader rd = new InputStreamReader(getClass().getClassLoader().getResourceAsStream(matrix.getPath()));
		SubstitutionMatrix<AminoAcidCompound> sb = new SimpleSubstitutionMatrix<AminoAcidCompound>(aa, rd, matrix.getPath());

		PairwiseSequenceAligner<ProteinSequence, AminoAcidCompound> result =  Alignments.getPairwiseAligner(s1,s2, PairwiseSequenceAlignerType.LOCAL, gp, sb);
		
		double minResidues = 100;
		double score = ((double) result.getPair().getNumIdenticals()/result.getPair().getLength()), threshold = 0.3;
		
		System.out.println("Length " + result.getPair().getLength());
		
		if(result.getPair().getLength()>=minResidues) {
			
			if(score > threshold) {
				
				System.out.println("Length ok. score "+score+"\t"+(score > threshold));
			}
		}
		else {
			
			while (result.getPair().getLength()<minResidues) {
				
				minResidues = minResidues/2;
				threshold +=0.2;
				System.out.println("new min residues "+minResidues);
				System.out.println("new threshold "+threshold);
				
				if(threshold>1)
					threshold=1;
			}
			
			System.out.println("Length not ok "+minResidues+" score "+(score > threshold));
		}
		
		System.out.println("Similarity "+result.getSimilarity());
		System.out.println("Score "+result.getScore());
		System.out.println("Min Score "+result.getMinScore());
		System.out.println("Max Score "+result.getMaxScore());
		System.out.println("Score calc "+ ((double) result.getScore()/result.getMaxScore()));
		System.out.println("Length "+result.getPair().getLength());
		System.out.println("query: "+result.getQuery().getLength()+"\ttarget: "+ result.getTarget().getLength() +"\tdistance: "+result.getDistance());
		System.out.println("Similarity score "+ ((double) result.getPair().getNumSimilars()/result.getPair().getLength()));
		System.out.println("Identity score "+ ((double) result.getPair().getNumIdenticals()/result.getPair().getLength()));
		System.out.println();
	}

	@Test
	public void newAlign() throws MalformedURLException, IOException, Exception{

		System.out.println("######################### new sw align #######################");

		//		String uni_id1 = "Q0VKV8";	// "Q0VNJ6";
		//		String uni_id2 = "C1B498";	// "Q93ZR6";
		//		String uni_id1 = "A1TX06";
		//		String uni_id2 = ("A1U572");
//				String uni_id1 = "P49374";
//				String uni_id2 = ("P0AGF4");
//				ProteinSequence s1 = FastaReaderHelper.readFastaProteinSequence(new URL(String.format("http://www.uniprot.org/uniprot/%s.fasta", uni_id1)).openStream()).get(uni_id1);
//				ProteinSequence s2 = FastaReaderHelper.readFastaProteinSequence(new URL(String.format("http://www.uniprot.org/uniprot/%s.fasta", uni_id2)).openStream()).get(uni_id2);

		ProteinSequence s1 = FastaReaderHelper.readFastaProteinSequence(new FileInputStream(new File("C:/Users/Oscar Dias/Desktop/seq1.txt"))).get("seq1");
		ProteinSequence s2 = FastaReaderHelper.readFastaProteinSequence(new FileInputStream(new File("C:/Users/Oscar Dias/Desktop/seq2.txt"))).get("seq2");

		Matrix matrix = Matrix.BLOSUM62;
		short gapOpenPenalty=10, gapExtensionPenalty=1;

		GapPenalty gp = new SimpleGapPenalty(gapOpenPenalty, gapExtensionPenalty);
		CompoundSet<AminoAcidCompound> aa = new AminoAcidCompoundSet();
		Reader rd = new InputStreamReader(getClass().getClassLoader().getResourceAsStream(matrix.getPath()));
		SubstitutionMatrix<AminoAcidCompound> sb = new SimpleSubstitutionMatrix<AminoAcidCompound>(aa, rd, matrix.getPath());

		SequencePair<ProteinSequence, AminoAcidCompound> pair = Alignments.getPairwiseAlignment(
				s1,
				s2,
				PairwiseSequenceAlignerType.LOCAL, gp, sb);

		System.out.printf("%n%s vs %s%n%s %n%s %n%s %n%s %n", pair.getQuery().getAccession(), pair.getTarget().getAccession(), pair, "NumSimilars "+pair.getNumSimilars(),"NumIdenticals "+pair.getNumIdenticals(),"Lenght "+ pair.getLength());

		for(AlignedSequence<ProteinSequence, AminoAcidCompound> p : pair.getAlignedSequences())
			System.out.println(p.getOriginalSequence().getOriginalHeader()+"\t"+p.getStart()+ " " + p.getEnd());

		System.out.println("Similarity score "+ ((double) pair.getNumSimilars()/pair.getLength()));
		System.out.println("Identity score "+ ((double) pair.getNumIdenticals()/pair.getLength()));
		System.out.println("Length " + pair.getLength());
		System.out.println();
		List<ProteinSequence> lps = new ArrayList<>();
		lps.add(s1);
		lps.add(s2);


		List<PairwiseSequenceScorer<ProteinSequence, AminoAcidCompound>> result =  Alignments.getAllPairsScorers(lps, PairwiseSequenceScorerType.LOCAL, gp, sb);

		for(PairwiseSequenceScorer<ProteinSequence, AminoAcidCompound> res : result) {

			System.out.println("distance "+res.getDistance());
			System.out.println(PairwiseSequenceScorerType.LOCAL+" "+res.getSimilarity());
			System.out.println("Score "+res.getScore());
			System.out.println("Min Score "+res.getMinScore());
			System.out.println("Max Score "+res.getMaxScore());
			System.out.println("Score calc "+ ((double) res.getScore()/res.getMaxScore()));
			System.out.println();
		}

		result =  Alignments.getAllPairsScorers(lps, PairwiseSequenceScorerType.LOCAL_SIMILARITIES, gp, sb);

		for(PairwiseSequenceScorer<ProteinSequence, AminoAcidCompound> res : result) {

			System.out.println("distance "+res.getDistance());
			System.out.println(PairwiseSequenceScorerType.LOCAL_SIMILARITIES+" "+res.getSimilarity());
			System.out.println("Score "+res.getScore());
			System.out.println("Min Score "+res.getMinScore());
			System.out.println("Max Score "+res.getMaxScore());
			System.out.println("Score calc "+ ((double) res.getScore()/res.getMaxScore()));
			System.out.println();
		}

		result =  Alignments.getAllPairsScorers(lps, PairwiseSequenceScorerType.LOCAL_IDENTITIES, gp, sb);

		for(PairwiseSequenceScorer<ProteinSequence, AminoAcidCompound> res : result) {

			System.out.println("distance "+res.getDistance());
			System.out.println(PairwiseSequenceScorerType.LOCAL_IDENTITIES+" "+res.getSimilarity());
			System.out.println("Min Score "+res.getMinScore());
			System.out.println("Score "+res.getScore());
			System.out.println("Max Score "+res.getMaxScore());
			System.out.println("Score calc "+ ((double) res.getScore()/res.getMaxScore()));
			System.out.println();
		}
	}

	//
	//	@Test
	//	public void testBiojava () throws Exception {
	//
	//
	//		String[] id = new String[] {"C1ATM0", "A4IWB4"};
	//
	//		URL uniprotFasta = new URL(String.format("http://www.uniprot.org/uniprot/%s.fasta", id[0]));
	//		ProteinSequence s1 = FastaReaderHelper.readFastaProteinSequence(uniprotFasta.openStream()).get(id[0]);
	//		System.out.printf("id : %s %s%n%s%n", id[0], s1, s1.getOriginalHeader());
	//
	//
	//		uniprotFasta = new URL(String.format("http://www.uniprot.org/uniprot/%s.fasta", id[1]));
	//		ProteinSequence s2 = FastaReaderHelper.readFastaProteinSequence(uniprotFasta.openStream()).get(id[1]);
	//		System.out.printf("id : %s %s%n%s%n", id[1], s2, s2.getOriginalHeader());
	//
	//
	//		{
	//			SubstitutionMatrix<AminoAcidCompound> matrix = new SimpleSubstitutionMatrix<AminoAcidCompound>();
	//			SequencePair<ProteinSequence, AminoAcidCompound> pair = Alignments.getPairwiseAlignment(s1, s2,
	//					PairwiseSequenceAlignerType.LOCAL, new SimpleGapPenalty(), matrix);
	//			System.out.printf("%n%s vs %s%n%s", pair.getQuery().getAccession(), pair.getTarget().getAccession(), pair);
	//		}
	//	}

	//@Test
	//	public void test() throws InterruptedException {
	//
	//		ConcurrentHashMap<String, ProteinSequence> querySequences = new ConcurrentHashMap<>(),
	//				staticGenesSet = new ConcurrentHashMap<>();
	//		Map<String, Set<String>> closestOrthologs = new HashMap<>(); 
	//		Map<String, Integer> kegg_taxonomy_scores = new HashMap<>();
	//		Set<AlignmentContainer> out = new HashSet<AlignmentContainer>();
	//
	//		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//		double similarity_threshold = 0.1;
	//		double referenceTaxonomyThreshold = 0.2;
	//		String ec_number = "3.1.3.4";
	//		int referenceTaxonomyScore = 0;
	//		String ko = null;
	//
	//		querySequences.put("ROP_50370", new ProteinSequence("MLFDQSVLDSAVEDRAPWLIDVVTVVTHSGGTVAAWIVSTVLTAVLLLQDRRREAVLVAG"+
	//				"AMLSGLAVMTALKNLFERQRPPLPDRLVEISSFSFPSGHAMMTAILASVVVAVMLRVVLV"+
	//				"RHVRIALVALLLLYTLAVGLSRVYLGAHWMTDVLAGWAFGALWAAFWILLTRNRPRSART"));
	//
	//		staticGenesSet.put("FTW_0252", new ProteinSequence("MNYFDHKLAKTTAIYGFLTLIIAILSYNFLDIKFATLVHSSELFGTGISTIAAFTSNIFS"+
	//				"PKVWTVITAIATVICIYKHIVKKPSQKLYIMSLSLIMTIIITTIVKVILARYRPEMLLFD"
	//				+ "NHYGFHFFSFKKAYNSMPSGHTALTFAGLLAIANFFEKKYITLIAIIISGLVAVSRIIIL"
	//				+ "DHFISDVIVAAYIGIFTYLWSKAFVESK"));
	//
	//		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	//		List<Thread> threads = new ArrayList<Thread>();
	//		ConcurrentLinkedQueue<String> queryArray = new ConcurrentLinkedQueue<String>(querySequences.keySet());
	//
	//		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	//		int numberOfCores = Runtime.getRuntime().availableProcessors();
	//		if(queryArray.size()<numberOfCores)
	//			numberOfCores=queryArray.size();
	//		System.out.println("number Of threads: "+numberOfCores);
	//
	//		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//
	//
	//		for(int i=0; i<numberOfCores; i++) {
	//
	//			Runnable lc	= new PairwiseSequenceAlignement(Method.SmithWaterman, querySequences, staticGenesSet, queryArray,  null, 
	//					similarity_threshold, null, null, new AtomicInteger(0), new AtomicBoolean(false), AlignmentPurpose.OTHER);
	//			((PairwiseSequenceAlignement) lc).setEc_number(ec_number);
	//			((PairwiseSequenceAlignement) lc).setKO(ko);
	//			((PairwiseSequenceAlignement) lc).setClosestOrthologs(closestOrthologs);
	//			((PairwiseSequenceAlignement) lc).setReferenceTaxonomyScore(referenceTaxonomyScore);
	//			((PairwiseSequenceAlignement) lc).setKegg_taxonomy_scores(kegg_taxonomy_scores);
	//			((PairwiseSequenceAlignement) lc).setReferenceTaxonomyThreshold(referenceTaxonomyThreshold);
	//			((PairwiseSequenceAlignement) lc).setAlignmentContainerSet(out);
	//
	//			Thread thread = new Thread(lc);
	//			threads.add(thread);
	//			System.out.println("Start "+i);
	//			thread.start();
	//		}
	//
	//		for(Thread thread :threads) {
	//
	//			thread.join();
	//		}
	//
	//		System.out.println(out);
	//	}

}
