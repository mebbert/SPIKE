/**
 * 
 */
package spike.datastructures;

import htsjdk.samtools.util.Interval;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;

import org.apache.log4j.Logger;

import spike.datastructures.SamLocusIterator.LocusInfo;

/**
 * @author markebbert
 *
 */
public class LocusInfoQueue {
	
	@SuppressWarnings("unused")
	private static Logger logger = Logger.getLogger(LocusInfoQueue.class);

	private LinkedList<LocusInfo> locusInfoList;

	private static int MAXDISTANCE;
	private static int MAXSIZE;

	
	/**
	 * Instantiate LocusInfoQueue with no maximum size
	 * @throws IOException 
	 */
	public LocusInfoQueue(final int maxDistance) throws IOException{
		this(maxDistance, -1);
	}

	/**
	 * Instantiate LocusInfoQueue with maximum size
	 * @throws IOException 
	 */
	public LocusInfoQueue(final int maxDistance, final int maxSize) {
		this.locusInfoList = new LinkedList<LocusInfo>();
		LocusInfoQueue.MAXDISTANCE = maxDistance;
		LocusInfoQueue.MAXSIZE = maxSize;
	}

	public int size() {
		return this.locusInfoList.size();
	}

	public boolean isEmpty() {
		return this.locusInfoList.isEmpty();
	}

	public Iterator<LocusInfo> iterator() {
		return this.locusInfoList.iterator();
	}

	public boolean remove(int pos) {
		return (this.locusInfoList.remove(pos) != null) ? true : false;
	}

	public void clear() {
		this.locusInfoList.clear();
	}

	/**
	 * Add locusInfo to the back of the queue. Will only
	 * keep up to LocusInfoQueue.MAXSIZE objects
	 * 
	 * @param locus
	 * @return true if successful
	 */
	public boolean add(LocusInfo locus) {
		this.locusInfoList.add(locus);

		if(LocusInfoQueue.MAXSIZE > -1 &&
				locusInfoList.size() > LocusInfoQueue.MAXSIZE){
			this.locusInfoList.remove();
		}
		
		return true;
	}


    /**
     * Borrowed and modified this code from assertOrderedNonOverlapping from
     * IntervalUtil.java in htsjdk.samtools.util. Will determine whether the
     * set of intervals are both ordered and non-overlapping. Intervals that
     * are technically in order, but farther than LocusInfoQueue.MAXDISTANCE
     * away will be considered out of order. Any sequence from a different
     * contig is out of order.
     * 
     * @param intervals
     * @param sequenceDictionary used to determine order of sequences
     */
    public static boolean intervalsOrderedNonOverlapping(
    		ArrayList<Interval> prevIntervals, ArrayList<Interval> currIntervals) {
    	
    	Interval prev1 = prevIntervals.get(0), prev2 = prevIntervals.get(1),
    			curr1 = currIntervals.get(0), curr2 = currIntervals.get(1);
    	
    	if(prev1 == null || curr1 == null){
    		throw new RuntimeException("ERROR: No primary intervals should be null");
    	}
    	
		/* Return true if any of them are ordered and overlapping */
		return (intervalsOrderedNonOverlapping(prev1, curr1)
				|| intervalsOrderedNonOverlapping(prev1, curr2)
				|| intervalsOrderedNonOverlapping(prev2, curr1)
				|| intervalsOrderedNonOverlapping(prev2, curr2));

    }
    
    
    public static boolean intervalsOrderedNonOverlapping(Interval prevInterval,
    		Interval currInterval){
     	   
		if(prevInterval == null || currInterval == null)
			return false;
		
		/* Test if the intervals overlap 
		 * TODO: Do we care if they intersect?
		 */
		if (prevInterval.intersects(currInterval)) {
			return false;
		}
		
		/* Test that the contigs are the same, and positions are ordered.*/
		if (prevInterval.compareTo(currInterval) >= 0) {
			return false;
		}
		
		/* Test that interval is not greater than max distance away */
		if(currInterval.getStart() - prevInterval.getEnd()
				> LocusInfoQueue.MAXDISTANCE){
			return false;
		}
		return true;   	
    }
    
    public String toString(){
    	StringBuilder sb = new StringBuilder();
    	for(LocusInfo locus : locusInfoList){
    		sb.append(locus.toString());
    		sb.append(", PrimaryIntervalMassTuple: ");
    		sb.append(locus.getTopRefIntervalMassTuple().toString());
    		sb.append(", SecondaryIntervalMassTuple: ");
    		sb.append(locus.getIntervalMassTuples().get(1).toString());
    		sb.append("; ");
    	}
    	return sb.toString();
    }
    
 }
