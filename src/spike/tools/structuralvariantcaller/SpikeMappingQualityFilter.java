/**
 * 
 */
package spike.tools.structuralvariantcaller;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.filter.SamRecordFilter;

/**
 * @author markebbert
 *
 */
public class SpikeMappingQualityFilter implements SamRecordFilter {
    private int lowMappingQualityThreshold = Integer.MIN_VALUE,
    		highMappingQualityThreshold = Integer.MAX_VALUE;

    public SpikeMappingQualityFilter(final int lowMappingQualityThreshold,
    		final int highMappingQualityThreshold) {
        this.lowMappingQualityThreshold = lowMappingQualityThreshold;
        this.highMappingQualityThreshold = highMappingQualityThreshold;
    }

    @Override
    public boolean filterOut(final SAMRecord record) {
    	
    	/* Filter reads where the MAPQ is BETWEEN the low and high
    	 * thresholds. Reads with MAPQ == 0 && MAPQ == 60 are the most
    	 * important.
    	 */
        return record.getMappingQuality() > this.lowMappingQualityThreshold
        		&& record.getMappingQuality() < this.highMappingQualityThreshold;
    }

    @Override
    public boolean filterOut(final SAMRecord first, final SAMRecord second) {
        return filterOut(first) || filterOut(second);
    }
}
