/**
 * 
 */
package spike.tools.structuralvariantcaller;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

import java.io.File;

import org.apache.log4j.Logger;

import spike.datastructures.ReferenceSegment;

/**
 * @author markebbert
 *
 */
public class StructuralVariantCaller {
	
	private static Logger logger = Logger.getLogger(StructuralVariantCaller.class);

	/**
	 * 
	 */
	public StructuralVariantCaller() {
		return;
	}
	
	public void startWalking(final File samFile, final int minSVSize,
			final int minMapQual, final ValidationStringency vs)
			throws StructuralVariantCallerException{
		
		SamReader reader = openSam(samFile, vs);
        final SAMFileHeader header = reader.getFileHeader();

        if (header.getSortOrder() != SortOrder.coordinate) {
        	throw new StructuralVariantCallerException("Input file " +
        			samFile.getAbsolutePath() + " is not coordinate sorted.");
        }   

        SAMRecordIterator samIt = reader.iterator();
        SAMRecord tmpRec = null, prevRec = null, currRec = null;
        ReferenceSegment prevRecRs = null, currRecRs = null;
        int prevRefRecAlignedStart, currRefRecAlignedStart,
        	prevRefRecAlignedEnd, currRefRecAlignedEnd;


        /* Get the first record */
        if(samIt.hasNext()){
        	prevRec = samIt.next();
        	prevRecRs = new ReferenceSegment(prevRec.getReadName());
        }

        /* Iterate over the records */
        while(samIt.hasNext()){
        	tmpRec = samIt.next();
        	
        	/* only consider reads if they are mapped and have mapping quality
        	 * greater than the minimum.
        	 */
        	if(!tmpRec.getReadUnmappedFlag()
        			&& tmpRec.getMappingQuality() > minMapQual){
        		
        		currRec = tmpRec;
				currRecRs = new ReferenceSegment(currRec.getReadName());
        		/* Calculate the aligned end of the reference segment. This
        		 * is the end position of the reference segment minus the number
        		 * of clipped bases.
        		 */
				prevRefRecAlignedEnd = prevRecRs.getEndPos() -
						(prevRec.getUnclippedEnd() - prevRec.getAlignmentEnd());
				
				currRefRecAlignedStart = currRecRs.getStartPos() +
						(currRec.getAlignmentStart() - currRec.getUnclippedStart());
				
				if(currRec.getAlignmentStart() == 51421496)
					logger.debug("here");
				
				if(currRefRecAlignedStart - prevRefRecAlignedEnd > minSVSize){
					logger.debug("Prev Rec: " + prevRec.toString());
					logger.debug("Curr Rec: " + currRec.toString());
					logger.debug("Diff: " + (currRefRecAlignedStart - prevRefRecAlignedEnd + "\n\n"));
				}

				prevRec = currRec;
				prevRecRs = currRecRs;
        	}
        }
	}
	
    private SamReader openSam(final File samFile, ValidationStringency vs) {
    	
//    	System.setProperty("java.io.tmpdir", "");

		final SamReaderFactory factory =
		          SamReaderFactory.makeDefault()
		              .enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS,
		            		  SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
		              .validationStringency(vs);

        final SamReader reader = factory.open(samFile);
        
        return reader;
    }
    

}
