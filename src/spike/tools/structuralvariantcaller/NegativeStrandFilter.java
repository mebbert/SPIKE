/**
 * 
 */
package spike.tools.structuralvariantcaller;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.filter.SamRecordFilter;

/**
 * @author markebbert
 *
 * Filter alignments to the minus strand.
 */
public class NegativeStrandFilter implements SamRecordFilter {


    public NegativeStrandFilter() {
    	return;
    }

    @Override
    public boolean filterOut(final SAMRecord record) {
        return record.getReadNegativeStrandFlag();
    }

    @Override
    public boolean filterOut(final SAMRecord first, final SAMRecord second) {
        return filterOut(first) || filterOut(second);
    }
}
