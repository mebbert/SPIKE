/**
 * 
 */
package spike.tools.structuralvariantcaller;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.TreeSet;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.log4j.Logger;

import spike.datastructures.LocusCoverageQueue;

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
	
	public void startWalking(final File samFile, final File outVCF,
			File sampleRef, final int minSVSize, final int minMapQual,
			final int minDepth, final ValidationStringency vs)
			throws StructuralVariantCallerException, IOException{
		
		SamReader reader = openSam(samFile, vs);
        final SAMFileHeader header = reader.getFileHeader();

        if (header.getSortOrder() != SortOrder.coordinate) {
        	throw new StructuralVariantCallerException("Input file " +
        			samFile.getAbsolutePath() + " is not coordinate sorted.");
        }   

        SAMRecordIterator samIt = reader.iterator();
        SAMRecord tmpRec = null,/* prevRec = null,*/ currRec = null;
//        ReferenceSegment prevRecRs = null, currRecRs = null;
//        int currRefRecAlignedStart, prevRefRecAlignedEnd;
        
        int anticipatedDepth = 250;
        ArrayList<SAMRecord> clippedStarts = new ArrayList<SAMRecord>(anticipatedDepth);
        ArrayList<SAMRecord> clippedEnds = new ArrayList<SAMRecord>(anticipatedDepth);
        
        LocusCoverageQueue runningDepth = new LocusCoverageQueue();
        
        TreeSet<String> samples = new TreeSet<String>();
        for(SAMReadGroupRecord group : header.getReadGroups()){
        	samples.add(group.getSample());
        }

        VCFWriter vw = new VCFWriter(outVCF, sampleRef, samples);
        
        IndexedFastaSequenceFile referenceReader = new IndexedFastaSequenceFile(sampleRef);
        
        boolean potentialSV = false;


        /* Get the first record */
//        if(samIt.hasNext()){
//        	prevRec = samIt.next();
//        	prevRecRs = new ReferenceSegment(prevRec.getReadName());
//        }

        /* Iterate over the records */
        while(samIt.hasNext()){
        	tmpRec = samIt.next();
        	
        	/* only consider reads if they are mapped and have mapping quality
        	 * greater than the minimum.
        	 */
        	if(!tmpRec.getReadUnmappedFlag()
        			&& tmpRec.getMappingQuality() > minMapQual){
        		
        		currRec = tmpRec;
        		
        		if(currRec.getAlignmentStart() % 100000 == 0){
					logger.debug(currRec.toString());
        		}
//				currRecRs = new ReferenceSegment(currRec.getReadName());
				
        		/* Do I ever clear this? */
        		runningDepth.addCoverageForRead(currRec);
				
				/* Store reads where the ends are clipped until we see a read
				 * that is unclipped, UNLESS we are starting to see reads where
				 * the start is clipped. Ignore any reads where both ends are
				 * clipped. Probably misaligned.
				 * 
				 * TODO: Test that there isn't coverage too?
				 */
				if(isStartClipped(currRec) && isEndClipped(currRec)){
					continue;
				}
				else if(isEndClipped(currRec)){
					clippedEnds.add(currRec);
				}
				else if(isStartClipped(currRec)){
					clippedStarts.add(currRec);
					
					/* If we're at the end of the file, see if we have enough
					 * coverage to test for an SV.
					 */
					if(!samIt.hasNext() && clippedEnds.size() >= minDepth){
						potentialSV = true;
					}
				}
				else if(clippedStarts.size() >= minDepth
						|| clippedEnds.size() >= minDepth){ 
					/* We've hit a read that is not clipped at all, but we have
					 * a set of start- and end-clipped reads. Let's look for a
					 * structural variant.
					 */
					
					potentialSV = true;
					
				}
				else if(!clippedStarts.isEmpty() || !clippedEnds.isEmpty()){
					clippedEnds.clear();
					clippedStarts.clear();
				}
				else if(clippedStarts.isEmpty() && clippedEnds.isEmpty()){
					continue;
				}
				else{
					/* If we get here, we want to know why */
					logger.debug(currRec.toString());
					throw new RuntimeException("Please report this bug. We need"
							+ " to know what condition caused it.");
				}
				
				/* Test for SV */
				if(potentialSV){
					/* 'Insertion' in this context is very general. This could
					 * be an inversion or inverted translocation
					 */
					int[] boundaries = testForInsertion(clippedStarts,
							clippedEnds, runningDepth, currRec);

					int start = boundaries[0], end = boundaries[1];

					if(start != -1){
						
						ArrayList<Allele> alleles = new ArrayList<Allele>();
						Allele ref = Allele.create(referenceReader.getSubsequenceAt(currRec.getContig(), start-1, start-1).getBases(), true);
						Allele alt = Allele.create(referenceReader.getSubsequenceAt(currRec.getContig(), start, end).getBases(), false);
						alleles.add(ref);
						alleles.add(alt);
						
						GenotypesContext gc = GenotypesContext.create(new GenotypeBuilder(header.getReadGroups().get(0).getSample(), alleles).make());
						
						VariantContextBuilder vcBuilder = new VariantContextBuilder();
						vcBuilder.chr(currRec.getContig());
						vcBuilder.start(start);
						vcBuilder.stop(start);
						vcBuilder.alleles(alleles);
						vcBuilder.genotypes(gc);
						
						vw.writeVariantToVCF(vcBuilder.make());
						
						logger.debug(
								"\n#################################################\n" +
								"# Found INSERTION between: " + start +
								" and " + end + " #\n" + 
								"#################################################\n");
						System.out.println("");
						
					}
					
					clippedEnds.clear();
					clippedStarts.clear();
					potentialSV = false;
				}

//
//        		/* Calculate the aligned end of the reference segment. This
//        		 * is the end position of the reference segment minus the number
//        		 * of clipped bases.
//        		 */
//				prevRefRecAlignedEnd = prevRecRs.getEndPos() -
//						(prevRec.getUnclippedEnd() - prevRec.getAlignmentEnd());
//				
//				currRefRecAlignedStart = currRecRs.getStartPos() +
//						(currRec.getAlignmentStart() - currRec.getUnclippedStart());
//				
//				if(currRec.getAlignmentStart() == 51421496)
//					logger.debug("here");
//				
//				if(currRefRecAlignedStart - prevRefRecAlignedEnd > minSVSize){
//					logger.debug("Prev Rec: " + prevRec.toString());
//					logger.debug("Curr Rec: " + currRec.toString());
//					logger.debug("Diff: " + (currRefRecAlignedStart - prevRefRecAlignedEnd + "\n\n"));
//				}
//
//				prevRec = currRec;
//				prevRecRs = currRecRs;
        	}
        }
        
		referenceReader.close();
	}
	
    /**
     * 
     * Open a SAM/BAM file for reading and return the SamReader obj
     * 
     * @param samFile
     * @param vs
     * @return SamReader
     */
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
    
    /**
     * Will determine whether a given read was clipped at the beginning
     * 
     * @param rec
     * @return boolean
     */
    private boolean isStartClipped(SAMRecord rec){
    	
    	if(rec.getAlignmentStart() - rec.getUnclippedStart() > 0)
    		return true;
    	return false;
    }
    
    /**
     * Will determine whether a given read was clipped at the end
     * @param rec
     * @return boolean
     */
    private boolean isEndClipped(SAMRecord rec){
    	
    	if(rec.getUnclippedEnd() - rec.getAlignmentEnd() > 0)
    		return true;
    	return false;
    }
    
    /**
     * Assess whether the lists of SAMRecords are bounding an insertion.
     * 
     * @param clippedStarts
     * @param clippedEnds
     */
    private int[] testForInsertion(ArrayList<SAMRecord> clippedStarts,
    		ArrayList<SAMRecord> clippedEnds, LocusCoverageQueue lcq,
    		SAMRecord currRec){
    	
    	double[] leftBoundaries = new double[clippedEnds.size()];
    	double[] rightBoundaries = new double[clippedStarts.size()];

    	/* Loop over the records with clipped ends and determine the left
    	 * boundary on the individual's contig.
    	 */
    	for(int i = 0; i < clippedEnds.size(); i++){
    		leftBoundaries[i] = clippedEnds.get(i).getAlignmentEnd();
    	}
    	
    	double endsMean = 0, startsMean = 0, endsSD = 0, startsSD = 0;
    	
    	if(leftBoundaries.length > 0){
		    DescriptiveStatistics endsStats =
				    new DescriptiveStatistics(leftBoundaries);
		    double endsMin = endsStats.getMin();
		    double endsMax = endsStats.getMax();
		    endsSD = endsStats.getStandardDeviation();
		    endsMean = endsStats.getMean();
		    double endsSkewness = endsStats.getSkewness();
		    
		    if(endsSD <= 5){
			    logger.debug("\nLeft boundary stats: \nRange: " + endsMin + "-" + endsMax
					    + "\nMean (Stand. Dev.): " + endsMean + " (" + endsSD + ")"
					    + "\nSkewness: " + endsSkewness
					    + "\nRecs: " + clippedEnds.toString() + "\n");
			    System.out.println("");
		    }
    	}
    	
    	if(rightBoundaries.length > 0){
		    for(int i = 0; i < clippedStarts.size(); i++){
			    rightBoundaries[i] = clippedStarts.get(i).getAlignmentStart();
		    }

		     DescriptiveStatistics startsStats =
				    new DescriptiveStatistics(rightBoundaries);
			double startsMin = startsStats.getMin();
			double startsMax = startsStats.getMax();
			startsSD = startsStats.getStandardDeviation();
			startsMean = startsStats.getMean();
			double startsSkewness = startsStats.getSkewness();
		    
		    if(startsSD <= 5){
			    logger.debug("\nRight boundary stats: \nRange: " + startsMin + "-" + startsMax
					    + "\nMean (Stand. Dev.): " + startsMean + " (" + startsSD + ")"
					    + "\nSkewness: " + startsSkewness
					    + "\nRecs: " + clippedStarts.toString() + "\n");
			    System.out.println("");
		    }
    	}
  
    	int[] boundaries = new int[]{-2, -2};
    	
    	if(endsMean == 1645061.0){
    		logger.info("here");
    	}
	    /* if we get a startsMean but not an ends mean, start at startsMean
	     * and go backwards until you hit a lot of coverage.
	     */
	    if(startsMean != 0 && endsMean == 0){
		    boundaries = getSVBoundaries(lcq, new Double(startsMean).intValue(),
				    new Double(endsMean).intValue());
	    }
	    /* do the opposite, but go until the start of the next read. */
	    else if(startsMean == 0 && endsMean != 0){
		    boundaries = getSVBoundaries(lcq, new Double(endsMean).intValue(),
				    new Double(currRec.getAlignmentStart()).intValue());
	    }
	    /* if we get both, check coverage across the defined region. */
	    if(startsMean != 0 && endsMean != 0){
		    /* TODO: require running average to be below some value? Look for bowl? */
		    boundaries = getSVBoundaries(lcq, new Double(endsMean).intValue(),
				    new Double(startsMean).intValue());
	    }
		    
		logger.debug("Boundaries: " + boundaries[0] + " to " + boundaries[1]);
		if((boundaries[0] == -1 && boundaries[1] != -1) ||
				(boundaries[0] != -1 && boundaries[1] == -1)){
			throw new RuntimeException("Found one boundary, but not the other!");
		}
    	
    	return boundaries;
    }


    
    /**
     * Calculate the left and right boundaries for the structural variant.
     * 
     * @param lcq
     * @param start
     * @param end
     * @return the left and right boundaries for the variant, or {-1, -1} if not
     * found
     */
    private int [] getSVBoundaries(LocusCoverageQueue lcq, int start, int end){

    	
    	/* TODO: Make this method call boundaries by majority. Genotype? 90/10 */
    	
    	int leftBoundary = -1, rightBoundary = -1;
    	
    	if(start < end){
    		
    		/* go from start to end (incrementing). The first gap is the left
    		 * boundary. Once you have a left boundary and then hit where there
    		 * is coverage, that's the right boundary.
    		 */
    		for(int i = start; i <= end; i++){ // less than or equal to. Include the next reads position
			    if(i - start > 5 && leftBoundary == -1){  /* quit if no coverage gap within 5 bases */
				    return new int[]{leftBoundary, rightBoundary};
			    }
			    else if(lcq.hasCoverageGap(i) && leftBoundary == -1){
				    leftBoundary = i;
			    }
			    else if(leftBoundary != -1 && !lcq.hasCoverageGap(i)){
				    rightBoundary = i - 1; // Take the position before

				    /* return immediately */
				    return new int[]{leftBoundary, rightBoundary};
			    }
    		}
    	}
    	else{
    		for(int i = start; i > end; i--){
			    if(start - i > 5 && rightBoundary == -1){  /* quit if no coverage gap within 5 bases */
				    return new int[]{leftBoundary, rightBoundary};
			    }
			    else if(lcq.hasCoverageGap(i) && rightBoundary == -1){
				    rightBoundary = i;
			    }
			    else if(rightBoundary != -1 && !lcq.hasCoverageGap(i)){
				    leftBoundary = i + 1;
				    
				    /* return immediately */
				    return new int[]{leftBoundary, rightBoundary};
    			}
    		}
    	}
    	
    	/* This will only happen when we start and end have a distance less 
    	 * than 5. e.g., if we have a "wall" of reads clipped at the ends and
    	 * the next read is less than 5 bps away.
    	 */
    	return new int[]{leftBoundary, rightBoundary};
    }
    

}
