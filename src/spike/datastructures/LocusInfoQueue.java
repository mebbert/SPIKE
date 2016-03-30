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
//	public Interval add(LocusInfo locusInfo) {
//		this.locusInfoList.add(locusInfo);
//
//		Interval interval = generateHumanGenomeRefIntervalByMass(locusInfo);
//		if(interval != null) this.humanGenRefIntervalListByMass.add(interval);
//		
//		if(LocusInfoQueue.MAXSIZE > -1 &&
//				locusInfoList.size() > LocusInfoQueue.MAXSIZE){
//			this.locusInfoList.remove();
//			this.humanGenRefIntervalListByMass.remove();
//		}
//		
//		return interval;
//	}
	
	
	/**
	 * Generate a list of Interval objects for each LocusInfo object. The
	 * interval positions will be adjusted for the human reference genome.
	 * The interval at this locus is based on the most common reference
	 * genome position, based on the reference segments origin in the 
	 * reference genome. The second will be null if it does not constitute
	 * at least 25% of the mass.
	 * 
	 * @return ArrayList<Interval> of the top two intervals, or null if
	 * there is a tie.
	 */
	public static ArrayList<Interval> generateHumanGenomeRefIntervalByMass(LocusInfo locus){
		
		/* TODO: When calculating interval by mass, I need the reads
		 * that have gaps at this location! They are not being included!
		 * Or are the positions at gaps not getting calculated properly?
		 */
		
		String hgRefLocation, hgRefContig;
		int hgRefStart, hgRefPos;
		Pattern p = Pattern.compile("(\\s*\\d+)\\[(\\d+)-\\d+\\]");
		Matcher m;

		List<RecordAndOffset> raoList = locus.getRecordAndPositions();
		
		ArrayList<String> hgRefPosList = new ArrayList<String>(raoList.size());
		for(RecordAndOffset rao : raoList){
			hgRefLocation = rao.getRecord().getReadName();
			
			/* Ignore reads on negative strand */
			if(rao.getRecord().getReadNegativeStrandFlag()) continue;
			
			m = p.matcher(hgRefLocation);

			if(!m.matches()){
				throw new RuntimeException("ERROR: All read names should match"
						+ "the regular expression. This one failed: "
						+ hgRefLocation);
			}

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
				  hgRefPosList.stream().collect(Collectors.groupingBy(w -> w,
						  Collectors.counting()));
		
		ArrayList<String> commonPosStr = (ArrayList<String>) LocusInfoQueue.getMostCommonValues(posFreqs);
		if(commonPosStr != null){

			ArrayList<Interval> topTwo = new ArrayList<Interval>();

			/* This interval will start and end at the same position */
			String[] commonContigPos1 = commonPosStr.get(0).split(":");
			Interval topInterval = new Interval(commonContigPos1[0],
					Integer.parseInt(commonContigPos1[1]),
					Integer.parseInt(commonContigPos1[1]));
			topTwo.add(topInterval);

			if(commonPosStr.get(1) != null){
				String[] commonContigPos2 = commonPosStr.get(1).split(":");
				Interval secondInterval = new Interval(commonContigPos2[0],
						Integer.parseInt(commonContigPos2[1]),
						Integer.parseInt(commonContigPos2[1]));
				topTwo.add(secondInterval);
			}
			else{
				topTwo.add(null);
			}

			return topTwo;
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
    		
		return intervalsOrderedNonOverlapping(prev1, curr1)
				|| intervalsOrderedNonOverlapping(prev1, curr2)
				|| intervalsOrderedNonOverlapping(prev2, curr1)
				|| intervalsOrderedNonOverlapping(prev2, curr2);

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
     * Determine which object has the highest count. Return the top two IF
     * the second is at least 25% of the mass.
     * 
     * @param freqs
     * @return
     */
    private static ArrayList<String> getMostCommonValues(Map<Object, Long> freqs){
    	Iterator<Object> it = freqs.keySet().iterator();
    	String o, mostCommon = null, secondCommon = null;
    	long count, secondMax = 0, max = 0, totalCount = 0;

    	double secondMaxProportion = 0.25;

    	while(it.hasNext()){
    		o = (String) it.next();
    		count = freqs.get(o);
    		totalCount += count;
    		if(count > max){
    			secondMax = max;
    			max = count;
    			secondCommon = mostCommon;
    			mostCommon = o;
    		}
    	}
    	
    	if(max >= secondMax){
    		ArrayList<String> topTwo = new ArrayList<String>();
    		topTwo.add(mostCommon);
    		
    		if(secondMax / totalCount >= secondMaxProportion){
				topTwo.add(secondCommon);
    		}
    		else{
    			topTwo.add(null);
    		}
    		return topTwo;
    	}
    	else{
    		throw new RuntimeException("ERROR: Max count should never be less"
    				+ " than the previous max!");
    	}
    }
	
}
