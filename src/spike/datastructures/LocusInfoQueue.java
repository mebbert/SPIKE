/**
 * 
 */
package spike.datastructures;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.SamLocusIterator;
import htsjdk.samtools.util.SamLocusIterator.LocusInfo;
import htsjdk.samtools.util.SamLocusIterator.RecordAndOffset;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

/**
 * @author markebbert
 *
 */
public class LocusInfoQueue {

	private LinkedList<LocusInfo> locusInfoList;
	private LinkedList<Interval> humanGenRefIntervalListByMass;

	private static int MAXSIZE;
	private static int MAXDISTANCE;

	
	/**
	 * Instantiate LocusInfoQueue with no maximum size
	 */
	public LocusInfoQueue(final int maxDistance){
		this(maxDistance, -1);
	}

	/**
	 * Instantiate LocusInfoQueue with maximum size
	 */
	public LocusInfoQueue(final int maxDistance, final int maxSize){
		this.locusInfoList = new LinkedList<LocusInfo>();
		this.humanGenRefIntervalListByMass = new LinkedList<Interval>();
		LocusInfoQueue.MAXDISTANCE = maxDistance;
		LocusInfoQueue.MAXSIZE = maxSize;
	}

	public int size() {
		return this.locusInfoList.size();
	}

	public boolean isEmpty() {
		return this.locusInfoList.isEmpty();
	}

	public Iterator<SamLocusIterator.LocusInfo> iterator() {
		return this.locusInfoList.iterator();
	}

	public boolean remove(int pos) {
		return (this.locusInfoList.remove(pos) != null) ? true : false;
	}

	public void clear() {
		this.locusInfoList.clear();
	}

	/**
	 * Add locusInfo to the back of the queue.
	 * 
	 * @param chromPos
	 * @param depth
	 */
	public Interval add(LocusInfo locusInfo) {
		this.locusInfoList.add(locusInfo);

		Interval interval = generateHumanGenomeRefIntervalByMass(locusInfo);
		if(interval != null) this.humanGenRefIntervalListByMass.add(interval);
		
		if(LocusInfoQueue.MAXSIZE > -1 &&
				locusInfoList.size() > LocusInfoQueue.MAXSIZE){
			this.locusInfoList.remove();
			this.humanGenRefIntervalListByMass.remove();
		}
		
		return interval;
	}
	
	
	/**
	 * Generate a list of Interval objects for each LocusInfo object. The
	 * interval positions will be adjusted for the human reference genome.
	 * The interval at this locus is based on the most common reference
	 * genome position, based on the reference segments origin in the 
	 * reference genome.
	 * 
	 * @return
	 */
	public static Interval generateHumanGenomeRefIntervalByMass(LocusInfo locus){
		
		String hgRefLocation, hgRefContig;
		int hgRefStart, hgRefPos;
		Pattern p = Pattern.compile("(\\s*\\d*)\\[(\\d+)-\\d+\\]");
		Matcher m;
		Interval interval;

		List<RecordAndOffset> raoList = locus.getRecordAndPositions();
		
		ArrayList<String> hgRefPosList = new ArrayList<String>(raoList.size());
		for(RecordAndOffset rao : raoList){
			hgRefLocation = rao.getRecord().getReadName();
			
			/* Ignore reads on negative strand */
			if(rao.getRecord().getMateNegativeStrandFlag()) continue;
			
			m = p.matcher(hgRefLocation);

			hgRefContig = m.group(1);
			hgRefStart = Integer.parseInt(m.group(2));
			
			/* Determine the position for this locus on the actual
			 * human reference genome.
			 */
			hgRefPos = hgRefStart + rao.getOffset();

			hgRefPosList.add(hgRefContig + ":" + hgRefPos);
		}
		
		/* Count the frequency for each value in the lists using lambda functions */
		Map<Object, Long> posFreqs = 
				  hgRefPosList.stream().collect(Collectors.groupingBy(w -> w, Collectors.counting()));
		
		String commonPosStr = (String) LocusInfoQueue.getMostCommonValue(posFreqs);
		if(commonPosStr != null){
			/* This interval will start and end at the same position */
			String[] commonContigPos = commonPosStr.split(":");
			interval = new Interval(commonContigPos[0],
					Integer.parseInt(commonContigPos[1]),
					Integer.parseInt(commonContigPos[1]));
			return interval;
		}
		else{
			/* There were multiple modes of the same height. Just ignore. */
			return null;
		}
	}
	
//	/**
//	 * Generate a list of Interval objects for each LocusInfo object. The
//	 * interval positions will be adjusted for the human reference genome.
//	 * The interval at this locus is based on the most common reference
//	 * genome position, based on the reference segments origin in the 
//	 * reference genome.
//	 * 
//	 * @return
//	 */
//	public List<Interval> generateHumanGenomeRefIntervalListByMass(){
//		ArrayList<Interval> intervalList = new ArrayList<Interval>();
//		
//		String hgRefLocation, hgRefContig;
//		int hgRefStart, hgRefPos;
//		Pattern p = Pattern.compile("(\\s*\\d*)\\[(\\d+)-\\d+\\]");
//		Matcher m;
//		Interval interval;
//
//		for(LocusInfo i : this.locusInfoList){
//			List<RecordAndOffset> raoList = i.getRecordAndPositions();
//			
//			ArrayList<Integer> hgRefPosList = new ArrayList<Integer>(raoList.size());
//			ArrayList<String> hgRefContigList = new ArrayList<String>(raoList.size());
//			for(RecordAndOffset rao : raoList){
//				hgRefLocation = rao.getRecord().getReadName();
//				
//				m = p.matcher(hgRefLocation);
//
//				hgRefContig = m.group(1);
//				hgRefStart = Integer.parseInt(m.group(2));
//				
//				/* Determine the position for this locus on the actual
//				 * human reference genome.
//				 */
//				hgRefPos = hgRefStart + rao.getOffset();
//
//				hgRefPosList.add(hgRefPos);
//				hgRefContigList.add(hgRefContig);
//			}
//			
//			/* Count the frequency for each value in the lists using lambda functions */
//			Map<Object, Long> posFreqs = 
//					  hgRefPosList.stream().collect(Collectors.groupingBy(w -> w, Collectors.counting()));
//			Map<Object, Long> contigFreqs = 
//					  hgRefContigList.stream().collect(Collectors.groupingBy(w -> w, Collectors.counting()));
//			
//			Integer commonPos = (Integer) getMostCommonValue(posFreqs);
//			String commonContig = (String) getMostCommonValue(contigFreqs);
//			if(commonPos == null || commonContig == null){
//				/* This interval will start and end at the same position */
//				interval = new Interval(commonContig, commonPos, commonPos);
//				intervalList.add(interval);
//			}
//			else{
//				/* There were multiple modes of the same height. Just ignore. */
//			}
//				
//		}
//		return intervalList;
//	}
	
    /**
     * Borrowed and modified this code from assertOrderedNonOverlapping from
     * IntervalUtil.java in htsjdk.samtools.util. Will determine whether the
     * set of intervals are both ordered and non-overlapping. Intervals that
     * are technically in order, but a far distance away will be considered out
     * of order. Any sequence from a different contig is out of order.
     * 
     * @param intervals
     * @param sequenceDictionary used to determine order of sequences
     */
    public static boolean intervalsOrderedNonOverlapping(
    		Interval prevInterval, Interval currInterval) {
    	   
		if(prevInterval == null || currInterval == null)
			return true;
		
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

    /**
     * Borrowed and modified this code from assertOrderedNonOverlapping from
     * IntervalUtil.java in htsjdk.samtools.util. Will determine whether the
     * set of intervals are both ordered and non-overlapping. Intervals that
     * are technically in order, but a far distance away will be considered out
     * of order. Any sequence from a different contig is out of order.
     * 
     * This will return the set of intervals that break order. Intervals will
     * be inserted in sets of two, where each pair represent both sides of a
     * given order-breaking boundary.
     *
     * @param intervals
     * @param sequenceDictionary used to determine order of sequences
     */
//    public ArrayList<Interval> intervalsOrderedNonOverlapping() {
//    	
//    	final Iterator<Interval> intervals =
//    			this.humanGenRefIntervalListByMass.iterator();
//    	
//    	ArrayList<Interval> orderBreakingLoci = new ArrayList<Interval>();
//
//        if (!intervals.hasNext()) {
//            return null;
//        }
//
//        Interval prevInterval = intervals.next();
//
//        while (intervals.hasNext()) {
//            final Interval interval = intervals.next();
//            
//            if(interval == null) continue;
//            
//            /* Test if the intervals overlap 
//             * TODO: Do we care if they intersect?
//             */
//            if (prevInterval.intersects(interval)) {
//            	orderBreakingLoci.add(prevInterval);
//            	orderBreakingLoci.add(interval);
//            }
//            
//            /* Test that the contigs are the same, and positions are ordered.*/
//            if (prevInterval.compareTo(interval) >= 0) {
//            	orderBreakingLoci.add(prevInterval);
//            	orderBreakingLoci.add(interval);
//            }
//            
//            /* Test that interval is not greater than max distance away */
//            if(interval.getStart() - prevInterval.getEnd()
//            		> LocusInfoQueue.MAXDISTANCE){
//            	orderBreakingLoci.add(prevInterval);
//            	orderBreakingLoci.add(interval);
//            }
//            prevInterval = interval;
//        }
//        return orderBreakingLoci;
//    }
// 
    
    /**
     * Determine which object has the highest count. If multiple have the
     * highest (i.e., they're equal), return null.
     * 
     * @param freqs
     * @return
     */
    private static Object getMostCommonValue(Map<Object, Long> freqs){
    	Iterator<Object> it = freqs.keySet().iterator();
    	Object o, mostCommon = null;
    	long count, secondMax = 0, max = 0;
    	while(it.hasNext()){
    		o = it.next();
    		count = freqs.get(o);
    		if(count > max){
    			max = count;
    			mostCommon = o;
    		}
    		/* If count ever equals existing max, track it. */
    		else if(count == max){
    			secondMax = count;
    		}
    	}
    	
    	if(max > secondMax){
    		return mostCommon;
    	}
    	else if(max == secondMax){
    		return null;
    	}
    	else{
    		throw new RuntimeException("ERROR: Max count should never be less"
    				+ " than the previous max!");
    	}
    }
	
}
