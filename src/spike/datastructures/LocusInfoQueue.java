/**
 * 
 */
package spike.datastructures;

import htsjdk.samtools.util.Interval;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import org.apache.log4j.Logger;

import spike.datastructures.SamLocusIterator.LocusInfo;
import spike.datastructures.SamLocusIterator.RecordAndOffset;

/**
 * @author markebbert
 *
 */
public class LocusInfoQueue {
	
	private static Logger logger = Logger.getLogger(LocusInfoQueue.class);

	private LinkedList<LocusInfo> locusInfoList;

	private static int MAXDISTANCE;
	
	File file = new File("locus_stats.txt");

	
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
		LocusInfoQueue.MAXDISTANCE = maxDistance;
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
	 * the the max is below the expected proportion by mass.
	 */
	public static ArrayList<Interval>
			generateHumanGenomeRefIntervalByMass(LocusInfo locus){
		
		/* TODO: When calculating interval by mass, I need the reads
		 * that have gaps at this location! They are not being included!
		 * Or are the positions at gaps not getting calculated properly?
		 */
		
		String hgRefLocation, hgRefContig;
		int hgRefStart, hgRefPos;
		Pattern p = Pattern.compile("(\\s*\\d+)\\[(\\d+)-\\d+\\]");
		Matcher m;

//		List<RecordAndOffset> raoList = locus.getRecordAndPositions();
		
		/* Include reads that "cover" this locus with an insertion when
		 * calculating the mass. This is for situations where properly
		 * aligned reads have an insertion, but improperly aligned reads
		 * do not. This also should not cause us to miss any structural
		 * variants because no insertion within a read will be that big.
		 */
		List<RecordAndOffset> raoWithDelsList =
				locus.getRecordAndPositionsWithDeletions();
		
		
		List<RecordAndOffset> raoList =
				locus.getRecordAndPositions();
		
		int depth = raoList.size();
		int depthWithDels = raoWithDelsList.size();
		
		ArrayList<String> hgRefPosList = new ArrayList<String>(raoWithDelsList.size());
		for(RecordAndOffset rao : raoWithDelsList){
			hgRefLocation = rao.getRecord().getReadName();
			
			/* Ignore reads on negative strand */
			/* All negative strand reads are currently being ignored 
			 * by filter in StructuralVariantCaller
			 */
//			if(rao.getRecord().getReadNegativeStrandFlag()) continue;
			
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
		
		/* Count the frequency for each value in the lists using lambda functions:
		 * http://stackoverflow.com/questions/24119065/rank-arrayliststring-with-int-in-java
		 */
		Map<String, Long> posFreqs = 
				  hgRefPosList.stream().collect(Collectors.groupingBy(w -> w,
						  Collectors.counting()));
		
		if(posFreqs.size() == 0){
			return null;
		}

		return LocusInfoQueue.getMostCommonValues(posFreqs, depthWithDels);
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
    
    
    /**
     * Determine which object has the highest count. Return the top two if
     * 1. The first is at least 50% of the mass
     * 2. The second is least 25% of the mass.
     * 
     * @param freqs
     * @return an ArrayList<Entry<String, Long>> of length two with the
     * top two if both conditions are true. The second value will be null
     * if the second condition is false. Return null if neither are true.
     */
    private static ArrayList<Interval>
    					getMostCommonValues(Map<String, Long> freqs, int totalDepth){
    	
    	ArrayList<Interval> topTwo = null;
    	
		/* Sorting with lambda function in Java 8:
		 * http://stackoverflow.com/questions/8119366/sorting-hashmap-by-values
		 */
		Map<String, Long> freqsSortedByValue = 
			     freqs.entrySet().stream()
			    .sorted(Entry.<String, Long>comparingByValue().reversed())
			    .collect(Collectors.toMap(Entry::getKey, Entry::getValue,
			                              (e1, e2) -> e1, LinkedHashMap::new));
		
		Iterator<String> it = freqsSortedByValue.keySet().iterator();
		String topPosition, secondTopPosition = null;
		Long topCount, secondTopCount = null;

		topPosition = it.next();
		topCount = freqsSortedByValue.get(topPosition);
		
		logger.debug("top pos: " + topPosition);
		if(topPosition.equals("1:10509")){
			logger.debug("here");
		}
		
		if(it.hasNext()){
			secondTopPosition = it.next();
			secondTopCount = freqsSortedByValue.get(secondTopPosition);
		}

		double topProportion = 0.5, secondProportion = 0.25;
		if(topCount / (double) totalDepth >= topProportion){
			topTwo = new ArrayList<Interval>();
			
			String[] topContigPos1 = topPosition.split(":");
			topTwo.add(new Interval(topContigPos1[0], 
					Integer.parseInt(topContigPos1[1]),
					Integer.parseInt(topContigPos1[1])));

			
			if(secondTopPosition == null){
				topTwo.add(null);
			}
			else if(secondTopCount / (double) totalDepth >= secondProportion){

				String[] topContigPos2 = secondTopPosition.split(":");
				topTwo.add(new Interval(topContigPos2[0], 
						Integer.parseInt(topContigPos2[1]),
						Integer.parseInt(topContigPos2[1])));
			}
			return topTwo;
		}
		return null;
	
    }

    
    /**
     * Determine which object has the highest count. Return the top two IF
     * the second is at least 25% of the mass.
     * 
     * @param freqs
     * @return
     */
//    private static ArrayList<String> getMostCommonValues(Map<String, Long> freqs){
//    	Iterator<String> it = freqs.keySet().iterator();
//    	String o, mostCommon = null, secondCommon = null;
//    	long count, secondMax = 0, max = 0, totalCount = 0;
//
//    	double secondMinProportion = 0.25;
//
//    	while(it.hasNext()){
//    		o = it.next();
//    		count = freqs.get(o);
//    		totalCount += count;
//    		if(count > max){
//    			secondMax = max;
//    			max = count;
//    			secondCommon = mostCommon;
//    			mostCommon = o;
//    		}
//    		else if(count > secondMax){
//    			secondMax = count;
//    			secondCommon = o;
//    		}
//    	}
//    	
//    	if(Integer.parseInt(mostCommon.split(":")[1]) > 249000000){
//    		logger.debug("most common: " + mostCommon);
//    		logger.debug("most common count: " + max);
//    		logger.debug("second common: " + secondCommon);
//    		logger.debug("second common count: " + secondMax);
//    		logger.debug("total count: " + totalCount + "\n\n");
//    	}
//    	
//    	if(max >= secondMax){
//    		ArrayList<String> topTwo = new ArrayList<String>();
//    		topTwo.add(mostCommon);
//    		
//    		/* Include the secondMax if it is secondMinProportion of
//    		 * the max (not total count). 
//    		 */
//    		if((secondMax / (double) max) >= secondMinProportion){
//				topTwo.add(secondCommon);
//    		}
//    		else{
//    			topTwo.add(null);
//    		}
//    		return topTwo;
//    	}
//    	else{
//    		throw new RuntimeException("ERROR: Max count should never be less"
//    				+ " than the previous max!");
//    	}
//    }
	
}
