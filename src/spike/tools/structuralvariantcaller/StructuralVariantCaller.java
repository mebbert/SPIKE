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
import htsjdk.samtools.filter.MappingQualityFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.SamLocusIterator;
import htsjdk.samtools.util.SamLocusIterator.LocusInfo;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Stack;
import java.util.TreeSet;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.log4j.Logger;

import spike.datastructures.LocusCoverageQueue;
import spike.datastructures.LocusInfoQueue;

/**
 * @author markebbert
 *
 */
public class StructuralVariantCaller /*implements Runnable*/ {
	
	private static Logger logger = Logger.getLogger(StructuralVariantCaller.class);
	private static final int MAX_STD_DEV = 5;
	private static int minDepth;
	private static int minSVSize;
	private static int minMapQual;
	private static ValidationStringency samValidationStringency;
	
//	private int startWalking, endWalking;
	
	private File samFile;

	private LinkedList<Interval> preList, postList,
				deletionPreList, deletionPostList;
	private LinkedList<VariantContext> deletionVarList;

	private IndexedFastaSequenceFile referenceReader;
	private SAMFileHeader header;
	private SamReader reader;

	private VCFWriter vw;

	/**
	 * @throws FileNotFoundException 
	 * 
	 */
	public StructuralVariantCaller(final File samFile, final File outVCF,
			File sampleRef, final int minSVSize, final int minMapQual,
			final int minDepth, final ValidationStringency vs/*, int startWalking,
			int endWalking*/) throws FileNotFoundException {
		
//		this.startWalking = startWalking;
//		this.endWalking = endWalking;
		
		
		StructuralVariantCaller.minSVSize = minSVSize;
		StructuralVariantCaller.minMapQual = minMapQual;
		StructuralVariantCaller.minDepth = minDepth;
		StructuralVariantCaller.samValidationStringency = vs;
		
		this.preList = new LinkedList<Interval>();
		this.postList = new LinkedList<Interval>();
		this.deletionVarList = new LinkedList<VariantContext>();
		
		this.referenceReader =
				new IndexedFastaSequenceFile(sampleRef);

		this.samFile = samFile;
		SamReader reader = openSam(samFile,
				StructuralVariantCaller.samValidationStringency);
		this.header = reader.getFileHeader();
		
		/* Get sample name(s) from the sam/bam file and create
		 * a VCFWriter.
		 */
        TreeSet<String> samples = new TreeSet<String>();
        for(SAMReadGroupRecord group : header.getReadGroups()){
        	samples.add(group.getSample());
        }
		this.vw = new VCFWriter(outVCF, sampleRef, samples);
		
	}

	
	/**
	 * @throws FileNotFoundException
	 */
	public void startWalkingByLocus() throws FileNotFoundException{

		SamLocusIterator sli = new SamLocusIterator(reader);
		

		/* Filter reads with a mapQ < StructuralVariantCaller.minMapQual 
		 * and reads aligning to the negative strand
		 */
		List<SamRecordFilter> samFilters = Arrays.asList(
				new MappingQualityFilter(StructuralVariantCaller.minMapQual),
				new NegativeStrandFilter());

		sli.setSamFilters(samFilters);
		
		Iterator<SamLocusIterator.LocusInfo> iter = sli.iterator();
		



		/* Keep a running set of loci.
		 * Only keep information for this many (MAXSIZE)
		 * loci.
		 */
		final int MAXSIZE = 50;
//		LocusInfoQueue running_liq =
//				new LocusInfoQueue(StructuralVariantCaller.minSVSize, MAXSIZE),
//				possibleVarLiq;
		LocusInfo locus;
		
		/* These locus Intervals are in context of the reference genome.
		 * We use the reference segment's original location (read name)
		 * to determine the locus interval.
		 */
		Interval prevLocusInterval = null, currLocusInterval;

		boolean inOrder;


		VariantContext var;
		
		/* Track positions coverage was lost and gained, as
		 * measured by minDepth.
		 */
		long positionCovLost, positionCovGained, locusCount = 0;
		
		while(iter.hasNext()){

			locus = iter.next();
				
			currLocusInterval =
					LocusInfoQueue.generateHumanGenomeRefIntervalByMass(locus);
			
			/* At every locus, check if a mass of reads are out of order
			 * between the previous and current locus
			 */
			inOrder = LocusInfoQueue.
						intervalsOrderedNonOverlapping(prevLocusInterval,
								currLocusInterval);
		
			/* If out of order, determine whether it's a deletion or
			 * translocation.
			 */
			if(!inOrder){

				/* If prev is before curr, we'll assume it's a deletion, unless
				 * we find it's a translocation, later.
				 */
				if(prevLocusInterval.compareTo(currLocusInterval) <= 0){
					
				}
				else{
					preList.push(prevLocusInterval);
					postList.push(currLocusInterval);
					resolveAllStructuralVariantsInRegion(iter);
					
					if(!preList.isEmpty() || !postList.isEmpty() ||
							!deletionPreList.isEmpty() ||
							!deletionPostList.isEmpty() ||
							!deletionVarList.isEmpty()){
						
						throw new RuntimeException("ERROR: all stacks should be"
								+ " empty. Something is very wrong!\nPreStack: "
								+ preList.toString() + "\n\nPostStack: " 
								+ postList.toString() + "\n\nDeletionPreStack: "
								+ deletionPreList.toString() + "\n\nDeletionPostStack: "
								+ deletionPostList.toString() + "\n\nDeletionVarStack: "
								+ deletionVarList.toString());
					}
				}
		
			}

			/* Every (MAXSIZE/2) loci, check for a mass of reads that are
			 * out of order.
			 */
//			if(locusCount++ % (MAXSIZE/2) == 0){
//				
//				/* If the running LocusInfoQueue is out of order, there
//				 * may be a deletion or translocation.
//				 */
//				ArrayList<Interval> orderBreakingLoci =
//						running_liq.intervalsOrderedNonOverlapping();
//				if(orderBreakingLoci.size() > 0){
//					
//					for(int i = 0; i < orderBreakingLoci.size(); i += 2){
//						orderBreakingLoci.get(i).
//					}
//				}
//			}

			int depth = locus.getRecordAndPositions().size();

			/* If we drop below acceptable coverage, we're looking at a
			 * possible insertion, duplication, or inversion. Go until
			 * we rise above. Duplications are caught here because we're
			 * filtering reads with mapQ < StructuralVariantCaller.minMapQual.
			 * Reads mapped to a duplication will have a mapQ == 0 because the
			 * mapQ is set to 0 if a given read can map to two locations equally
			 * well.
			 */
			if(depth < StructuralVariantCaller.minDepth){
				
				/* Create new LIQ for the low coverage region.
				 * No max size.
				 */
//				possibleVarLiq =
//						new LocusInfoQueue(StructuralVariantCaller.minSVSize);
//				possibleVarLiq.add(locusInfo);

				positionCovLost = locus.getPosition();
				
				/* Keep iterating until we rise above minDepth */
				while(iter.hasNext()){
					locus = iter.next();
					if(depth >= StructuralVariantCaller.minDepth){
						positionCovGained = locus.getPosition();

						logger.debug("Potential insertion/inversion/duplication between: " +
								positionCovLost + " and " + positionCovGained +
								"; Depth: " + depth);

						break;
					}
					else{
//						possibleVarLiq.add(locusInfo);
					}
				}

			}
			else{
				/* Only add the locus if it has acceptable coverage */
//				running_liq.add(locusInfo);
			}
			
			prevLocusInterval = currLocusInterval;

		}

		sli.close();
	}
	

	/**
	 * @param iter
	 */
	private void resolveAllStructuralVariantsInRegion(
			Iterator<SamLocusIterator.LocusInfo> iter){

		
		/* These locus Intervals are in context of the reference genome.
		 * We use the reference segment's original location (read name)
		 * to determine the locus interval.
		 */
		Interval refStartInterval, refEndInterval, varStartInterval,
					varEndInterval, pre, post;

		VariantContext var;
		
		pre = preList.peek();
		post = postList.peek();
		
		/* If prev is before curr, we'll assume it's a deletion, unless
		 * we find it's a translocation, later. We must already be in a
		 * variant where an end order-breaking
		 * boundary exists, so call it a deletion, push it on the stack,
		 * and proceed. 
		 */
		if(pre.getContig().equals(post.getContig())
				&& pre.getEnd() < post.getStart()){

			var = buildVariant(this.referenceReader, this.header,
					pre.getContig(),
					pre.getStart(),
					post.getStart(),
					pre.getStart(),
					pre.getStart());

			deletionVarList.push(var);

			/* Have to pop so we don't get in an infinite loop,
			 * but track these in case they're actually a translocation.
			 */
			deletionPreList.push(preList.pop());
			deletionPostList.push(postList.pop());
			
			/* keep resolving */
			resolveAllStructuralVariantsInRegion(iter);
		}
		else{ // Must be translocation, run until we return to prevLocusInterval
			
			ArrayList<Interval> orderBreakingBoundary =
					walkToNextOrderBreakingBoundary(iter);

			pre = orderBreakingBoundary.get(0);
			post = orderBreakingBoundary.get(1);

			/* If the 'post' interval (the right side of the next
			 * order-breaking boundary) is "in order" with the 
			 * the previous 'pre' interval, then we've found the
			 * end of this variant.
			 */
			for(int i = preList.size() - 1; i >= 0; i--){

				if(LocusInfoQueue.intervalsOrderedNonOverlapping(
								this.preList.get(i), post)){
					
					/* If the deletionPostStack.peek() is lexicographically before
					 * pre, then this must have been another translocation, not a
					 * deletion. The deletionPostStack would actually mark the
					 * beginning of the translocation, because it was the 'post'
					 * position of that boundary.
					 */
					if(deletionPostList.peek().compareTo(pre) <= 0){
						varStartInterval = deletionPostList.pop();
						varEndInterval = pre;
						refStartInterval = deletionPreList.pop();
						refEndInterval = refStartInterval;
						var = buildVariant(this.referenceReader, this.header,
								refStartInterval.getContig(),
								refStartInterval.getStart(), 
								refEndInterval.getStart(),
								varStartInterval.getStart(),
								varEndInterval.getStart());
						this.vw.writeVariantToVCF(var);
						
						deletionVarList.pop(); // get rid of the false deletion

					}
					
					/* Create and write the variant */
					varStartInterval = this.postList.get(i);
					varEndInterval = pre;
					refStartInterval = this.preList.get(i);
					refEndInterval = refStartInterval;
					var = buildVariant(this.referenceReader, this.header,
							refStartInterval.getContig(),
							refStartInterval.getStart(), 
							refEndInterval.getStart(),
							varStartInterval.getStart(),
							varEndInterval.getStart());
					this.vw.writeVariantToVCF(var);
				}
				else{
					preList.push(pre);
					postList.push(post);
					resolveAllStructuralVariantsInRegion(iter);
				}
			}
		}
	}
	
	/**
	 * Walk until we hit another order-breaking boundary and return
	 * the two intervals (positions) that are out of order. Or return
	 * null if we never find one.
	 * 
	 * @param iter
	 * @return
	 */
	private ArrayList<Interval> walkToNextOrderBreakingBoundary(
			Iterator<SamLocusIterator.LocusInfo> iter){

		LocusInfo locus;

		/* These locus Intervals are in context of the reference genome.
		 * We use the reference segment's original location (read name)
		 * to determine the locus interval.
		 */
		Interval prevLocusInterval = null, currLocusInterval;
		
		boolean inOrder;
		
		while(iter.hasNext()){
			locus = iter.next();
			currLocusInterval = 
					LocusInfoQueue.generateHumanGenomeRefIntervalByMass(locus);
			
			/* Check if currLocusInterval is out of order
			 * with prevLocusInterval
			 */
			inOrder = LocusInfoQueue.
					intervalsOrderedNonOverlapping(
							prevLocusInterval, currLocusInterval);
			
			if(!inOrder){
				ArrayList<Interval> orderBreakingIntervals = new ArrayList<Interval>();
				orderBreakingIntervals.add(prevLocusInterval);
				orderBreakingIntervals.add(currLocusInterval);
				return orderBreakingIntervals;
			}
			prevLocusInterval = currLocusInterval;
		}
		
		/* We went to the end of the contig without hitting a new boundary */
		return null;
	}
	
	public void startWalking()
			throws StructuralVariantCallerException, IOException{
		

//        DiskBasedBAMFileIndex dbfi = new DiskBasedBAMFileIndex(new File(samFile.getAbsoluteFile() + ".bai"),
//        		header.getSequenceDictionary());
        
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


        boolean potentialSV = false;


        /* Get the first record */
//        if(samIt.hasNext()){
//        	prevRec = samIt.next();
//        	prevRecRs = new ReferenceSegment(prevRec.getReadName());
//        }

        /* Iterate over the records */
        int count = 0;
        while(samIt.hasNext()){
        	tmpRec = samIt.next();
        	count++;
        	
        	if(count % 10000000 == 0){
        		logger.debug("Assessed " + count + " reads.");
        	}
        	
//        	Integer as = (Integer) tmpRec.getAttribute("AS");
//        	Integer xs = (Integer) tmpRec.getAttribute("XS");
//        	
//        	if(xs != null){
//        		logger.debug("AS: " + as);
//        		logger.debug("XS: " + xs);
//        	}
        	
        	/* only consider reads if they are mapped and have mapping quality
        	 * greater than the minimum.
        	 * 
        	 * TODO: look for any regions where there's a bunch of reads aligning
        	 * to the '-' strand. What would it mean?
        	 */
        	if(!tmpRec.getReadUnmappedFlag()
        			/* && !tmpRec.getReadNegativeStrandFlag() */
        			/* && tmpRec.getMappingQuality() > minMapQual */){
        		
        		currRec = tmpRec;
        		
        		if(currRec.getAlignmentStart() % 10000000 == 0){
					logger.debug(currRec.toString() + " Aligned at start = "
							+ currRec.getAlignmentStart());
        		}
//				currRecRs = new ReferenceSegment(currRec.getReadName());
				
        		/* TODO: Do I ever clear this? */
        		runningDepth.addCoverageForRead(currRec);
				
				/* Store reads where the ends are clipped until we see a read
				 * that is unclipped, UNLESS we start seeing reads where
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

					/* Ignore if the boundaries are not greater than 0 */
					if(start >= 0){
						
						vw.writeVariantToVCF(buildVariant(referenceReader,
								header, currRec.getContig(), start-1, start-1, 
								start-1, end));
						
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
	 * Build a VariantContext object based on the reference genome. Will return
	 * null if an 'N' is found in either the reference or alternate alleles.
	 * 
	 * @param referenceReader
	 * @param header
	 * @param currRec
	 * @param start
	 * @param end
	 * @return
	 */
	private VariantContext buildVariant(IndexedFastaSequenceFile referenceReader,
			SAMFileHeader header, String contig, int refStart, int refEnd,
			int altStart, int altEnd){

		/*
		 * TODO: 
		 * 1. Convert to b37 positions and write a VCF with
		 * those positions?
		 * 2. If #1, get the b37 positions. Use the read names
		 * to get this? Search for sequence in b37?
		 * 3. Get correct genotypes based on percentages (90/10).
		 * This will require figuring out both haplotypes in the
		 * assembly
		 */
		ArrayList<Allele> alleles = new ArrayList<Allele>();
		Allele ref = Allele.create(
				referenceReader.getSubsequenceAt(contig, refStart, refEnd).getBases(), true);
		Allele alt = Allele.create(
				referenceReader.getSubsequenceAt(contig, altStart, altEnd).getBases(), false);
		
		/* Ignore any regions containing "N"s. */
		if(ref.getBaseString().contains("N") || alt.getBaseString().contains("N")){
			return null;
		}
		alleles.add(ref);
		alleles.add(alt);
		
		GenotypesContext gc = GenotypesContext.
				create(new GenotypeBuilder(header.getReadGroups().get(0).getSample(), alleles).make());
		
		return new VariantContextBuilder()
			.chr(contig)
			.start(refStart)
			.stop(refEnd)
			.alleles(alleles)
			.genotypes(gc)
			.make();
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
		    endsSD = endsStats.getStandardDeviation();
		    endsMean = endsStats.getMean();
		    
//		    if(endsSD <= 5){
//			    logger.debug("\nLeft boundary stats: \nRange: " + endsMin + "-" + endsMax
//					    + "\nMean (Stand. Dev.): " + endsMean + " (" + endsSD + ")"
//					    + "\nSkewness: " + endsSkewness
//					    + "\nRecs: " + clippedEnds.toString() + "\n");
//			    System.out.println("");
//		    }
    	}
    	
    	if(rightBoundaries.length > 0){
		    for(int i = 0; i < clippedStarts.size(); i++){
			    rightBoundaries[i] = clippedStarts.get(i).getAlignmentStart();
		    }

		     DescriptiveStatistics startsStats =
				    new DescriptiveStatistics(rightBoundaries);
			startsSD = startsStats.getStandardDeviation();
			startsMean = startsStats.getMean();
		    
//		    if(startsSD <= 5){
//			    logger.debug("\nRight boundary stats: \nRange: " + startsMin + "-" + startsMax
//					    + "\nMean (Stand. Dev.): " + startsMean + " (" + startsSD + ")"
//					    + "\nSkewness: " + startsSkewness
//					    + "\nRecs: " + clippedStarts.toString() + "\n");
//			    System.out.println("");
//		    }
    	}
  
    	int[] boundaries = new int[]{-2, -2};
    	
    	if(endsMean == 1645061.0){
    		logger.info("here");
    	}
	    /* if we get a startsMean but not an ends mean, start at startsMean
	     * and go backwards until you hit a lot of coverage.
	     */
	    if(startsSD < MAX_STD_DEV && startsMean != 0 && endsMean == 0){
		    boundaries = getSVBoundaries(lcq, new Double(startsMean).intValue(),
				    new Double(endsMean).intValue());
	    }
	    /* else, if we have an ends mean, but no starts mean, do the opposite.
	     * but go until the start of the next read.
	     */
	    else if(startsMean == 0 && endsSD < MAX_STD_DEV && endsMean != 0){
		    boundaries = getSVBoundaries(lcq, new Double(endsMean).intValue(),
				    new Double(currRec.getAlignmentStart()).intValue());
	    }
	    /* if we get both, check coverage across the defined region. */
	    if(startsSD < MAX_STD_DEV && startsMean != 0 &&
	    		endsSD < MAX_STD_DEV && endsMean != 0){

		    /* TODO: require running average to be below some value? Look for bowl? */
		    boundaries = getSVBoundaries(lcq, new Double(endsMean).intValue(),
				    new Double(startsMean).intValue());
	    }
		    
	    if(boundaries[0] >= 0 && boundaries[1] >= 0){
			logger.debug("Boundaries: " + boundaries[0] + " to " + boundaries[1]);
	    }
		if((boundaries[0] == -1 && boundaries[1] != -1) ||
				(boundaries[0] != -1 && boundaries[1] == -1)){
			throw new RuntimeException("Found one boundary, but not the other!");
		}
		
		if(boundaries[0] == 267510){
			logger.debug("here");
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

//	@Override
//	public void run() {
//		// TODO Auto-generated method stub
//		
//	}
    

}
