/**
 * 
 */
package spike.datastructures;

import htsjdk.samtools.SAMRecord;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.TreeSet;

/**
 * @author markebbert
 *
 */
public class LocusCoverageQueue {
	
	HashMap<Integer, Integer> positionCoverageMap;
	
	public LocusCoverageQueue(){
		this.positionCoverageMap = new HashMap<Integer, Integer>();
	}

	public int size() {
		return this.positionCoverageMap.size();
	}

	public boolean isEmpty() {
		return this.positionCoverageMap.isEmpty();
	}

	/**
	 * Will determine if the queue contains the chromosome position. The
	 * position should be specified as chr:pos (e.g., 1:1245).
	 * 
	 * @param chromPos
	 * @return
	 */
	public boolean contains(int pos) {
		return this.positionCoverageMap.containsKey(pos);
	}

	public Iterator<Integer> iterator() {
		return this.positionCoverageMap.keySet().iterator();
	}

	public boolean remove(int pos) {
		return (this.positionCoverageMap.remove(pos) != null) ? true : false;
	}

	public void clear() {
		this.positionCoverageMap.clear();
	}

	/**
	 * Add the position and depth.
	 * 
	 * @param chromPos
	 * @param depth
	 */
	public void add(int pos, int depth) {
		this.positionCoverageMap.put(pos, depth);
	}
	
	/**
	 * Add the position and depth. If the position already exists, it will
	 * add 'depth' to the existing depth value.
	 * 
	 * @param rec
	 */
	public void addCoverageForRead(SAMRecord rec){
		
		int start = rec.getAlignmentStart();
		int end = rec.getAlignmentEnd();
		Integer currVal;
		
		for(int i = start; i < end; i++){
			
			if((currVal = this.positionCoverageMap.get(i)) != null){
				this.positionCoverageMap.put(i, currVal++);
			}
			else{
				this.positionCoverageMap.put(i, 1); // Just add depth = 1
			}
		}
	}
	
	/**
	 * Return whether there is a coverage gap at this position. There is
	 * a gap if the depth = 0.
	 * 
	 * @param pos
	 * @return
	 */
	public boolean hasCoverageGap(int pos){
		
		if(!this.positionCoverageMap.containsKey(pos)){
			return true;
		}
		return false;
	}

}
