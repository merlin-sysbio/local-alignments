package pt.uminho.ceb.biosystems.merlin.local.alignments.core.ModelMerge;

import java.util.concurrent.ConcurrentLinkedQueue;

import pt.uminho.ceb.biosystems.merlin.utilities.containers.capsules.AlignmentCapsule;

/**
 * @author amaromorais
 *
 */
public interface ModelAlignments extends Runnable {
	
	/**
	 * 
	 */
	@Override
	public void run();
	
	/**
	 * 
	 */
	public void buildAlignmentCapsules();
	
	/**
	 * 
	 */
	public ConcurrentLinkedQueue<AlignmentCapsule> getAlignmentsCapsules();

	
}
