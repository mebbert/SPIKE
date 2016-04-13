/**
 * 
 */
package spike.tools.structuralvariantcaller;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.Interval;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.stream.Collectors;

import org.apache.log4j.Logger;

import spike.datastructures.IntervalMassTuple;
import spike.datastructures.LocusInfoQueue;
import spike.datastructures.SamLocusIterator;
import spike.datastructures.SamLocusIterator.LocusInfo;
import spike.datastructures.SamLocusIterator.RecordAndOffset;
import spike.datastructures.StructuralVariantBoundary;
import spike.tools.utilitybelt.UtilityBelt;

/**
 * @author markebbert
 *
 */
public class StructuralVariantCaller /*implements Runnable*/ {
	
	private static Logger logger = Logger.getLogger(StructuralVariantCaller.class);
	private static int MIN_DEPTH, MIN_SV_SIZE, LOW_MAP_QUAL, HIGH_MAP_QUAL,
						 WINDOW_SIZE, INFO_QUEUE_MAX_SIZE = 20,
						 AVERAGE_POS_QUEUE_MAX_SIZE = 20;
	private static ValidationStringency SAM_VALIDATION_STRINGENCY;
	
	private ArrayList<StructuralVariantBoundary> svBoundaries;
	
	Writer writer;
	
//	private int startWalking, endWalking;
	

//	private LinkedList<Interval> sampPreList, sampPostList, hgPreList, hgPostList,
//				sampDeletionPreList, sampDeletionPostList, hgDeletionPreList,
//				hgDeletionPostList;
	private LinkedList<VariantContext> deletionVarList;
	
	/* Variant list where Integer is the variant's starting position based
	 * on the sample's contig.
	 */
	private TreeMap<Integer, VariantContext> varList;

	private IndexedFastaSequenceFile sampleRefReader, hgRefReader;
	private SAMFileHeader header;
	private SamReader reader;
	

	private VCFWriter vw;

	/**
	 * @throws IOException 
	 * 
	 */
	public StructuralVariantCaller(final File samFile, final File outVCF,
			File sampleRef, File hgRef, final int minSVSize, final int lowMapQual,
			final int highMapQual,
			final int minDepth, final int maxAcceptableCoverage,
			final int maxAcceptableClip,
			final ValidationStringency vs/*, int startWalking,
			int endWalking*/) throws IOException {
		
//		this.startWalking = startWalking;
//		this.endWalking = endWalking;
		
		writer = new BufferedWriter(new OutputStreamWriter(
	              new FileOutputStream("locus_stats.txt"), "utf-8"));
		writer.write("Position\tTotalDepth\tTopHGRefChrom\tTopHGRefPosition\tTopHGRefCount"
				+ "\tSecondTopHGRefChrom\tSecondTopHGRefPosition\tSecondTopHGRefCount");


		svBoundaries = new ArrayList<StructuralVariantBoundary>();
		
		
		StructuralVariantCaller.MIN_SV_SIZE = minSVSize;
		StructuralVariantCaller.LOW_MAP_QUAL = lowMapQual;
		StructuralVariantCaller.HIGH_MAP_QUAL = highMapQual;
		StructuralVariantCaller.MIN_DEPTH = minDepth;
		StructuralVariantCaller.SAM_VALIDATION_STRINGENCY = vs;
		StructuralVariantCaller.WINDOW_SIZE = 5;
		
//		this.sampPreList = new LinkedList<Interval>();
//		this.sampPostList = new LinkedList<Interval>();
//		this.hgPreList = new LinkedList<Interval>();
//		this.hgPostList = new LinkedList<Interval>();
//		this.sampDeletionPreList = new LinkedList<Interval>();
//		this.sampDeletionPostList = new LinkedList<Interval>();
//		this.hgDeletionPreList = new LinkedList<Interval>();
//		this.hgDeletionPostList = new LinkedList<Interval>();
		this.deletionVarList = new LinkedList<VariantContext>();
		this.varList = new TreeMap<Integer, VariantContext>();
		
		this.sampleRefReader =
				new IndexedFastaSequenceFile(sampleRef);
		this.hgRefReader =
				new IndexedFastaSequenceFile(hgRef);

		this.reader = UtilityBelt.openSam(samFile,
				StructuralVariantCaller.SAM_VALIDATION_STRINGENCY);
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
	 * All variants will be represented in a VCF as follows:
	 * 1. Contig name: set by the sample's contig name
	 * 2. Position: is the position on the sample's contig where
	 *              the variant is observed.
	 * 3. Ref: Is the sequence from the HG reference
	 * 4. Alt: The sequence observed on the sample
	 * 5. For translocations, inversions/transversions, and copy gains, the
	 *    INFO field will also contain where in the HG Ref it
	 *    came from.
	 * @throws IOException 
	 */
	public void startWalkingByLocus() throws IOException{

		SamLocusIterator sli = new SamLocusIterator(reader);
		
		/* Set the base quality score cutoff to 0 
		 * so the iterator won't try to validate base qualities. Any
		 * secondary alignment won't have the qualities, and they're
		 * all made up anyway.
		 */
		int baseQualityScoreCutoff = 0;
		sli.setQualityScoreCutoff(baseQualityScoreCutoff);
		

		/* Filter reads with a mapQ < StructuralVariantCaller.minMapQual,
		 * overclipped reads,
		 * and reads aligning to the negative strand
		 */
		List<SamRecordFilter> samFilters = Arrays.asList(
				new SpikeMappingQualityFilter(StructuralVariantCaller.LOW_MAP_QUAL, 
						StructuralVariantCaller.HIGH_MAP_QUAL),
//				new OverclippedReadFilter(maxAcceptableClip, filterSingleEndClips),
				new NegativeStrandFilter());

		sli.setSamFilters(samFilters);
		
		Iterator<LocusInfo> iter = sli.iterator();
		

		/* Keep a running set of loci.
		 * Only keep information for this many (MAXSIZE)
		 * loci.
		 */
		LocusInfoQueue runningLiq =
				new LocusInfoQueue(StructuralVariantCaller.MIN_SV_SIZE, INFO_QUEUE_MAX_SIZE),
				possibleVarLiq;
		LinkedList<Interval> averagePosQueue = new LinkedList<Interval>();

		LocusInfo locus;
		
		/* These locus Intervals are in context of the reference genome.
		 * We use the reference segment's original location (read name)
		 * to determine the locus interval.
		 */
		Interval /*prevHGRefLocusInterval = null,*/ currHGRefLocusInterval,
				prevSampLocusInterval = null, currSampLocusInterval,
				preHGRefLocus, postHGRefLocus, preSampLocus, postSampLocus,
				refVarStart, refVarEnd, sampVarStart, sampVarEnd,
				sampRefLocusIntervalCovLost,
				sampRefLocusIntervalCovGained, hgRefLocusIntervalCovLost,
				hgRefLocusIntervalCovGained;
		ArrayList<IntervalMassTuple> currTopTwoHGRefLociIntervals,
				prevTopTwoHGRefLociIntervals;

		boolean inOrder;
		StructuralVariantBoundary svBoundary;


		VariantContext var;
		
		/* Track positions coverage was lost and gained, as
		 * measured by minDepth.
		 */
		long count = 0;
		
		/* Run until we hit the first position with coverage */
		while(iter.hasNext()){
			locus = iter.next();
			int depth = locus.getRecordAndPositions().size();
			if(depth > 0){
				break;
			}
//			if(locus.getPosition() % 10000 == 0){
//				logger.debug("Processed " + locus.getPosition() + " loci.");
//			}
//			if(locus.getPosition() > 66405){
//				break;
//			}
		}
		
		
		/* Get the first position and set to 'prev' */
		if(iter.hasNext()){
			locus = iter.next();
			prevTopTwoHGRefLociIntervals =
					generateHumanGenomeRefIntervalByMass(locus);
		}
		else{
			return;
		}

		while(iter.hasNext()){

			locus = iter.next();
			
			logger.debug("locus: " + locus.getPosition());
        	
        	if(++count % 10000 == 0){
        		logger.debug("Assessed " + count + " loci.");
        	}
        	
        	if(locus.getPosition() == 36634){
        		logger.debug("here");
        	}

			/* Get depth at this position.
			 */
			int depth = locus.getRecordAndPositions().size();
			
			currTopTwoHGRefLociIntervals =
					generateHumanGenomeRefIntervalByMass(locus);
			

			currSampLocusInterval =
					new Interval(locus.getSequenceName(), locus.getPosition(),
							locus.getPosition());
			
			
//			if(currTopTwoHGRefLociIntervals == null ||
//					prevTopTwoHGRefLociIntervals == null){
//				
//				/* keep prev up with curr while going through
//				 * regions with too much coverage
//				 */
//				prevSampLocusInterval = currSampLocusInterval;
//				prevTopTwoHGRefLociIntervals = currTopTwoHGRefLociIntervals;
//				continue;
//			}
			
			
			/* If we bump into reads from another contig, determine
			 * whether we're seeing a translocation.
			 */
			if(crossesContigs(iter, runningLiq, prevTopTwoHGRefLociIntervals,
					currTopTwoHGRefLociIntervals, prevSampLocusInterval,
					currSampLocusInterval)){
				
				/* record the boundary */
				svBoundaries.add(new StructuralVariantBoundary(
						prevTopTwoHGRefLociIntervals.get(0).interval,
						currTopTwoHGRefLociIntervals.get(0).interval,
						prevSampLocusInterval,
						currSampLocusInterval));
				
				/* Clear the liq and average queue*/
				runningLiq.clear();
				averagePosQueue.clear();
				
				continue;
			}
			
			
			/* add the LocusMassIntervals to the locusInfo object and
			 * add the locus to the queue
			 */
			locus.addRefIntervalMassTuples(currTopTwoHGRefLociIntervals);
			runningLiq.add(locus);
			
			calculateNextSlidingWindowMass(runningLiq, averagePosQueue);

			/* Determine if we've had any major breaks in order */
			svBoundary = checkForBoundaryBreak(averagePosQueue, runningLiq);
			if(svBoundary != null){
				svBoundaries.add(svBoundary);

				/* Clear the liq and average queue*/
				runningLiq.clear();
				averagePosQueue.clear();
			}
			
			
			/* At every locus, check if a mass of reads are out of order
			 * between the previous and current locus. The order is based
			 * on where the sequences came from in the Human Reference Genome
			 * (e.g., b37, hg19, etc.).
			 * 
			 * We're considering the two most frequent ref loci aligning to this
			 * locus. This accounts for low-complexity situations.
			 * 
			 * TODO: We may need to keep track of the last several positions and
			 * not rely solely on the last position.
			 */
//			inOrder = LocusInfoQueue.
//						intervalsOrderedNonOverlapping(prevTopTwoHGRefLociIntervals,
//								currTopTwoHGRefLociIntervals);

		
			/* If out of order, determine whether it's a deletion or
			 * translocation.
			 */
//			if(!inOrder){
				
				/* 'pre' and 'post' refer to the intervals on either side of the
				 * order-breaking boundary. HGRef is based on the human reference
				 * sequence position, and 'Samp' is based on the individual's
				 * actual contig. The 'Samp' is for reporting the variant
				 * according to the individual's position and 'HGRef' is so
				 * we know where in the reference genome everything associates
				 * with.
				 */
//				preHGRefLocus = prevTopTwoHGRefLociIntervals.get(0).interval; // take the first if none were in order
//				postHGRefLocus = currTopTwoHGRefLociIntervals.get(0).interval;
//				
//				preSampLocus = prevSampLocusInterval;
//				postSampLocus = currSampLocusInterval;
//
//				/* If pre is before post (based on original location in human
//				 * reference genome), we'll assume it's a deletion, unless
//				 * we find it's a translocation, later.
//				 * 
//				 * TODO: Keep looping until the next boundary and determine if it was
//				 * a deletion, or something else.
//				 */
//				if(preHGRefLocus.getContig().equals(postHGRefLocus.getContig())
//						&& preHGRefLocus.getEnd() < postHGRefLocus.getStart()){
//					
//					if(postHGRefLocus.getStart() - preHGRefLocus.getEnd()
//							> StructuralVariantCaller.minSVSize){
//						
//
//						/* For deletions, the alt's start and end location will
//						 * be the nucleotide (on the sample's sequence) just
//						 * before the deleted sequence (sampleVarStart).
//						 * 
//						 * The ref allele's start will be at the same location,
//						 * and will end on the last nucleotide of the deleted
//						 * sequence.
//						 */
//						var = buildVariant(this.header,
//								preSampLocus.getContig(), // Samp contig
//								preSampLocus.getEnd(),    // samp locus before the deletion
//								preSampLocus.getEnd(), // also samp locus before the deletion
//								preHGRefLocus.getContig(),// Ref contig
//								preHGRefLocus.getEnd(), // first locus of deletion with HG ref coordinates
//								postHGRefLocus.getStart() - 1 // last locus of deletion
//								);
//							
//						if(var != null){
//							deletionVarList.add(var);
////							pendingVarList.put(preHGRef.getStart(), var);
//							pendingVarList.put(preSampLocus.getStart(), var);
//							
////							hgPreList.add(preHGRefLocus);
////							hgPostList.add(postHGRefLocus);
////							sampPreList.add(preSampLocus);
////							sampPostList.add(postSampLocus);
////
////							hgDeletionPreList.add(preHGRefLocus);
////							hgDeletionPostList.add(postHGRefLocus);	
////							sampDeletionPreList.add(preSampLocus);
////							sampDeletionPostList.add(postSampLocus);	
//						}
//
//					}
//				}
//				else{ // must be a translocation or copy gain
					
//					logger.debug("PreList size: " + sampPreList.size());
					
					/* Loop backwards over the lists to see if we can
					 * resolve any boundaries
					 */
//					for(int i = sampPreList.size() - 1; i >= 0; i--){
//						
//						if(LocusInfoQueue.intervalsOrderedNonOverlapping(
//										this.sampPreList.get(i), postSampLocus)){
//
//							/* See if preList.get(i) was previously treated
//							 * as a deletion. If so, get rid of it
//							 * 
//							 * TODO: Make sure the Interval.class equals method
//							 * is being called appropriately.
//							 */
//							int index = sampDeletionPreList.indexOf(sampPreList.get(i));
//							if(index >= 0){
//								pendingVarList.remove(deletionVarList.get(index).getStart());
//								sampDeletionPreList.remove(index);
//								sampDeletionPostList.remove(index);
//								hgDeletionPreList.remove(index);
//								hgDeletionPostList.remove(index);
//								deletionVarList.remove(index);
//							}
//							
//							
//							/* Since we found the match, now loop forward from
//							 * where we are, resolving all boundaries in between.
//							 */
//							for(int j = i + 1; j < sampPreList.size(); ){
//
//								/* For translocations/copy gains, the ref allele's
//								 * start and end location will
//								 * be the nucleotide (on the HG ref's sequence) just
//								 * before the inserted sequence (refVarStart).
//								 * 
//								 * The alt allele's start will be at the same location,
//								 * and will end on the last nucleotide of the inserted
//								 * sequence.
//								 */
//
//								/* Create and write the variant */
//								refVarStart = this.hgPreList.get(i);
//								refVarEnd = refVarStart;
//								sampVarStart = this.sampPreList.get(i);
//								sampVarEnd = this.sampPreList.get(j);
//								var = buildVariant(this.header,
//										sampVarStart.getContig(),
//										sampVarStart.getStart(),
//										sampVarEnd.getEnd(),
//										refVarStart.getContig(),
//										refVarStart.getEnd(),
//										refVarEnd.getStart()
//										);
////								this.vw.writeVariantToVCF(var);
//								
//								if(var != null){
//									pendingVarList.put(sampVarStart.getStart(), var);
//									
//									/* remove this boundary from pre and post list */
//									sampPreList.remove(i);
//									sampPostList.remove(i);
//									hgPreList.remove(i);
//									hgPostList.remove(i);
//
//									/* don't increment i or j, because we just
//									 * removed the objects at i.
//									 */
//								}
//								else{
//									logger.debug("Var == null");
//									logger.debug("i: " + i);
//									logger.debug("j: " + j);
//								}
//								
//							}
//							
//							
//							break;
//						}
//					}

					/* If we didn't resolve this boundary with any existing,
					 * just add them.
					 */
//					hgPreList.add(preHGRefLocus);
//					hgPostList.add(postHGRefLocus);
//					sampPreList.add(preSampLocus);
//					sampPostList.add(postSampLocus);
//				}
//		
//			}

			/* If we drop below acceptable coverage, we're looking at a
			 * possible insertion, duplication, or inversion. Go until
			 * we rise above. Duplications are caught here because we're
			 * filtering reads with mapQ < StructuralVariantCaller.minMapQual.
			 * Reads mapped to a duplication will have a mapQ == 0 because the
			 * mapQ is set to 0 if a given read can map to two locations equally
			 * well.
			 */
			if(depth < StructuralVariantCaller.MIN_DEPTH){
				
				/* Create new LIQ for the low coverage region.
				 * No max size.
				 */
//				possibleVarLiq =
//						new LocusInfoQueue(StructuralVariantCaller.minSVSize);
//				possibleVarLiq.add(locusInfo);

				/* The position we dropped below acceptable coverage */
				hgRefLocusIntervalCovLost =
						generateHumanGenomeRefIntervalByMass(locus)
						.get(0).interval; // We trust the first in this situation.
				sampRefLocusIntervalCovLost = new Interval(
								locus.getSequenceName(), locus.getPosition(),
								locus.getPosition());

				/* May be null if there are multiple statistical modes */
				if(hgRefLocusIntervalCovLost != null){
//				logger.debug("CovLost (" + locus.getPosition() + "): " + hgRefLocusIntervalCovLost.getStart());
//				logger.debug("depth (" + locus.getPosition() + "): " + depth);


					/* Keep iterating until we rise above minDepth */
					while(iter.hasNext()){
	
						locus = iter.next();
						ArrayList<IntervalMassTuple> topTwo = 
								generateHumanGenomeRefIntervalByMass(locus);

						if(topTwo != null){
							currHGRefLocusInterval = topTwo.get(0).interval;

							depth = locus.getRecordAndPositions().size();
//							logger.debug("depth (" + locus.getPosition() + "): " + depth);
							if(depth >= StructuralVariantCaller.MIN_DEPTH){

								hgRefLocusIntervalCovGained = currHGRefLocusInterval;

								sampRefLocusIntervalCovGained = new Interval(
														locus.getSequenceName(),
														locus.getPosition(),
														locus.getPosition());

								/* If the interval, according to the coverage
								 * gap in the Sample, is smaller than minSVSize,
								 * just move on.
								 */
								if(sampRefLocusIntervalCovGained.getStart() -
									sampRefLocusIntervalCovLost.getStart()
									< StructuralVariantCaller.MIN_SV_SIZE){
									
									/* break if the structural variant wasn't
									 * large enough.
									 */
									break;
								}

								/* Create and write the variant */
								logger.debug("CovGained (" + locus.getPosition() + "): " + hgRefLocusIntervalCovGained.getStart() + "\n");
								
								if(!hgRefLocusIntervalCovLost.getContig()
										.equals(hgRefLocusIntervalCovGained.getContig())){
									throw new RuntimeException("ERROR: HG Ref"
											+ " contigs do not match before and"
											+ " after insertion!");
								}

								/* TODO: also pass in where the insertion came from
								 * on the reference and put it in the info field
								 */
								var = buildVariant(this.header,
										sampRefLocusIntervalCovLost.getContig(), // Sample contig
										sampRefLocusIntervalCovLost.getStart() - 1, // Sample Var starts where coverage was lost
										sampRefLocusIntervalCovGained.getStart() - 1, // Sample Var ends one before where coverage was gained
										hgRefLocusIntervalCovLost.getContig(),  // HGRef contig
										hgRefLocusIntervalCovLost.getEnd() - 1,	   // position in HGRef just before the var
										hgRefLocusIntervalCovLost.getEnd() - 1 // Should theoretically match the position before the var starts
										);

								if(var != null){
//									pendingVarList.put(hgRefLocusIntervalCovLost.getStart() - 1, var);
									varList.put(sampRefLocusIntervalCovLost.getStart() - 1, var);
								}

								/* Write any deletions that remain in the list */
								/* TODO: Make sure this is the appropriate place. Want
								 * to do it once we've seen the next insertion/inversion/duplication
								 */
//								writePendingVars();
								break;
							}
							else{
		//						possibleVarLiq.add(locusInfo);
							}
						}
//						prevHGRefLocusInterval = currHGRefLocusInterval;
						prevTopTwoHGRefLociIntervals = currTopTwoHGRefLociIntervals;
						prevSampLocusInterval = currSampLocusInterval;
					}
				}

			}
			else{
				/* Only add the locus if it has acceptable coverage */
//				running_liq.add(locusInfo);
			}
			
//			prevHGRefLocusInterval = currHGRefLocusInterval;
			prevTopTwoHGRefLociIntervals = currTopTwoHGRefLociIntervals;
			prevSampLocusInterval = currSampLocusInterval;
		}

		for(StructuralVariantBoundary svb : svBoundaries){
			System.out.println(svb);
		}
		sli.close();
		
		resolveVariantsBetween(0, svBoundaries.size());
	}
	
	
	private void resolveVariantsBetween(int beginIndex, int endIndex){
		
		/* Loop over the boundaries to resolve all 
		 * non-deletion variants
		 */
		for(int i = beginIndex; i < endIndex; i++){
			for(int j = i+1; j < endIndex; j++){
				if(LocusInfoQueue.
						intervalsOrderedNonOverlapping(svBoundaries.get(i).hgPre,
								svBoundaries.get(j).hgPost)){
					
					/* create variant and remove boundaries, then resolve all
					 * in between.
					 */
					
					/* If i and j are adjacent boundaries, they must be independent
					 * deletions.
					 */
					if(i+1 != j){
						
						/* If there is one boundary between i and j, it must
						 * be a deletion. If there are multiple boundaries
						 * between i and j, recurse.
						 */
						if(i+1 != j+1){
							resolveVariantsBetween(i+1, j-1);
						}
						else{ // This must be a deletion
							/* Create the variant and advance i to j, which will
							 * increment to j+1 at the end of the loop.
							 */
						}
					}
					i = j;
					break; /* Break from j loop and keep going */
				}
			}
			
			/* If we get through all of 'j' without finding a boundary match,
			 * this must be a deletion.
			 * 
			 * Throw new RuntimeException if the signature does not match
			 * a deletion. Something must be wrong.
			 */
//			buildVariant(header, sampleContig, sampVarStart, sampVarEnd, refContig, refStart, refEnd)
		}
		
		/* Loop backwards over the lists to see if we can
		 * resolve any boundaries
		 */
//		for(int i = sampPreList.size() - 1; i >= 0; i--){
//			
//			if(LocusInfoQueue.intervalsOrderedNonOverlapping(
//							this.sampPreList.get(i), postSampLocus)){
//
//				/* See if preList.get(i) was previously treated
//				 * as a deletion. If so, get rid of it
//				 * 
//				 * TODO: Make sure the Interval.class equals method
//				 * is being called appropriately.
//				 */
//				int index = sampDeletionPreList.indexOf(sampPreList.get(i));
//				if(index >= 0){
//					pendingVarList.remove(deletionVarList.get(index).getStart());
//					sampDeletionPreList.remove(index);
//					sampDeletionPostList.remove(index);
//					hgDeletionPreList.remove(index);
//					hgDeletionPostList.remove(index);
//					deletionVarList.remove(index);
//				}
//				
//				
//				/* Since we found the match, now loop forward from
//				 * where we are, resolving all boundaries in between.
//				 */
//				for(int j = i + 1; j < sampPreList.size(); ){
//
//					/* For translocations/copy gains, the ref allele's
//					 * start and end location will
//					 * be the nucleotide (on the HG ref's sequence) just
//					 * before the inserted sequence (refVarStart).
//					 * 
//					 * The alt allele's start will be at the same location,
//					 * and will end on the last nucleotide of the inserted
//					 * sequence.
//					 */
//
//					/* Create and write the variant */
//					refVarStart = this.hgPreList.get(i);
//					refVarEnd = refVarStart;
//					sampVarStart = this.sampPreList.get(i);
//					sampVarEnd = this.sampPreList.get(j);
//					var = buildVariant(this.header,
//							sampVarStart.getContig(),
//							sampVarStart.getStart(),
//							sampVarEnd.getEnd(),
//							refVarStart.getContig(),
//							refVarStart.getEnd(),
//							refVarEnd.getStart()
//							);
////					this.vw.writeVariantToVCF(var);
//					
//					if(var != null){
//						pendingVarList.put(sampVarStart.getStart(), var);
//						
//						/* remove this boundary from pre and post list */
//						sampPreList.remove(i);
//						sampPostList.remove(i);
//						hgPreList.remove(i);
//						hgPostList.remove(i);
//
//						/* don't increment i or j, because we just
//						 * removed the objects at i.
//						 */
//					}
//					else{
//						logger.debug("Var == null");
//						logger.debug("i: " + i);
//						logger.debug("j: " + j);
//					}
//					
//				}
//				
//				
//				break;
//			}
//		}
	}
	
	
	
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
	 * @throws IOException 
	 */
	private ArrayList<IntervalMassTuple>
			generateHumanGenomeRefIntervalByMass(LocusInfo locus) throws IOException{
		
		/* Include reads that "cover" this locus with a deletion when
		 * calculating the mass. This is for situations where properly
		 * aligned reads have an deletion, but improperly aligned reads
		 * do not.
		 * 
		 * TODO: We will need to take additional precautions to catch
		 * structural variants between the minimum size and the reference
		 * segment (read) size for that portion of the HG Ref.
		 */
		List<RecordAndOffset> raoWithDelsList =
				locus.getRecordAndPositionsWithDeletions();
		Map<String, Integer> hgRefPositionFreqsSorted =
				locus.getHGRefPositionFreqsWithDeletionsSortedByValue();
		
		int depthWithDels = raoWithDelsList.size();
		
		if(depthWithDels == 0){
			return null;
		}

		return getMostCommonValues(locus, hgRefPositionFreqsSorted, depthWithDels);
	}
	
	
	   
    /**
     * Determine which object has the highest count. Return the top two if
     * 1. The first is at least 10% of the mass
     * 2. The second is least 5% of the mass.
     * 
     * @param freqs
     * @return an ArrayList<Entry<String, Long>> of length two with the
     * top two if both conditions are true. The second value will be null
     * if the second condition is false. Return null if neither are true.
     * @throws IOException 
     */
    private ArrayList<IntervalMassTuple>
    					getMostCommonValues(LocusInfo locus,
    							Map<String, Integer> freqsSortedByValue,
    							int totalDepth) throws IOException{
    	
    	ArrayList<IntervalMassTuple> topTwo = null;
		
		Iterator<String> it = freqsSortedByValue.keySet().iterator();
		String topPosition, secondTopPosition = null;
		Integer topCount, secondTopCount = null;
		
//		logger.debug("locus pos: " + locus.getPosition());

		/* Get the top position */
		topPosition = it.next();
		topCount = freqsSortedByValue.get(topPosition);
		
		String secondTopChrom = null, secondTopPos = null;
		
		/* Get the second most frequent position, if it exists. */
		if(it.hasNext()){
			secondTopPosition = it.next();
			secondTopCount = freqsSortedByValue.get(secondTopPosition);
			String[] secondTopToks = secondTopPosition.split(":");
			secondTopChrom = secondTopToks[0];
			secondTopPos = secondTopToks[1];
		}
		
		/* Temporary. Just tracking stats. */
		String[] topPosToks = topPosition.split(":");
		writer.write(locus.getPosition() + "\t" + totalDepth + "\t"
				+ topPosToks[0] + "\t" + topPosToks[1] + "\t" + topCount + "\t"
				+ secondTopChrom + "\t" + secondTopPos + "\t"
				+ "\t" + secondTopCount + "\n");

		/* Define required proportions and calculate them */
		double topRequiredProportion = 0.1, secondRequiredProportion = 0.05,
				topProportion = topCount/(double) totalDepth;

		if(topCount / (double) totalDepth >= topRequiredProportion){
			topTwo = new ArrayList<IntervalMassTuple>();
			
			/* Create IntervalMassTuple */
			String[] topContigPos1 = topPosition.split(":");
			IntervalMassTuple imt = new IntervalMassTuple(
					new Interval(topContigPos1[0], 
						Integer.parseInt(topContigPos1[1]),
						Integer.parseInt(topContigPos1[1])),
					topProportion);

			/* Store it */
			topTwo.add(imt);
			
			/* Do it for the second interval, it necessary */
			double secondProportion;
			if(secondTopPosition == null){
				topTwo.add(null);
			}
			else if((secondProportion = secondTopCount / (double) totalDepth)
						> secondRequiredProportion){
				String[] topContigPos2 = secondTopPosition.split(":");
				IntervalMassTuple imt2 = new IntervalMassTuple(
						new Interval(topContigPos2[0], 
							Integer.parseInt(topContigPos2[1]),
							Integer.parseInt(topContigPos2[1])),
						secondProportion);
				topTwo.add(imt2);

			}
			else{
				topTwo.add(null);
			}
			return topTwo;
		}
		return null;
    }
    
    /**
     * Calculate the sliding window mass by StructuralVariantCaller.windowSize.
     * By the time this method is called, any cross-contig translocations
     * should be caught and accounted for, clearing the runninLiq. So any
     * cross-contig masses in runningLiq should only be blips. Just ignore
     * them.
     * 
     * @param runningLiq
     * @param averagePosQueue
     */
    private void calculateNextSlidingWindowMass(LocusInfoQueue runningLiq,
    		LinkedList<Interval> averagePosQueue){
    	
    	int count = 0, liqSize = runningLiq.size();
    	
    	/* quit if liq is not at least as big as StructuralVariantCaller.windowSize */
    	if(liqSize < StructuralVariantCaller.WINDOW_SIZE) return;
    	
    	Iterator<LocusInfo> it = runningLiq.iterator();
    	LocusInfo locus;
    	double sum = 0;
    	int sumCount = 0;
    	String currentHGRefContig = null;
    	ArrayList<String> observedHGRefContigs = new ArrayList<String>();
    	IntervalMassTuple topRefIntervalMassTuple;
    	while(it.hasNext()){
    		locus = it.next();
    		topRefIntervalMassTuple = locus.getTopRefIntervalMassTuple();
    		if(sum == 0 && topRefIntervalMassTuple != null){
    			observedHGRefContigs.add(topRefIntervalMassTuple.interval.getContig());
    		}
    		
    		/* Only include the last 'StructuralVariantCaller.windowSize' values */
    		if(count++ >= liqSize - StructuralVariantCaller.WINDOW_SIZE){
    			
    			/* We have to determine which HGRef contig we're dealing with.
    			 * Use the most frequent contig up until we start calculating mass.
    			 */
    			if(currentHGRefContig == null){
    				/* Count the frequency for each value in the lists using lambda functions:
    				 * http://stackoverflow.com/questions/24119065/rank-arrayliststring-with-int-in-java
    				 */
    				Map<String, Long> posFreqs = 
    						  observedHGRefContigs.stream().collect(Collectors.groupingBy(w -> w,
    								  Collectors.counting()));
    				
    				/* Sorting with lambda function in Java 8:
    	    		 * http://stackoverflow.com/questions/8119366/sorting-hashmap-by-values
    	    		 */
    				Map<String, Long> posFreqsSortedByValue=
    						posFreqs.entrySet().stream()
    						.sorted(Entry.<String, Long>comparingByValue().reversed())
    						.collect(Collectors.toMap(Entry::getKey, Entry::getValue,
    								(e1, e2) -> e1, LinkedHashMap::new));
    				
    				currentHGRefContig = posFreqsSortedByValue.keySet().iterator().next();
    			}
    			if(topRefIntervalMassTuple != null &&
    					currentHGRefContig.equals(topRefIntervalMassTuple.interval.getContig())){
					sum += locus.getTopRefIntervalMassTuple().interval.getStart();
					sumCount++; // track how many we could actually use in this window
    			}
    		}
    	}
    	
    	if(sumCount >= 0){
    		int avg = new Double(sum/sumCount).intValue();
			averagePosQueue.add(new Interval(currentHGRefContig, avg, avg));
    	}
    	else{
    		averagePosQueue.add(null);
    	}
    	
    	/* Remove front element if too big */
    	if(averagePosQueue.size() > StructuralVariantCaller.AVERAGE_POS_QUEUE_MAX_SIZE){
    		averagePosQueue.remove();
    	}
    }
    
    /**
     * This method will check for boundary breaks using the averagePosQueue.
     * 
     * @param averagePosQueue
     * @return
     */
    private StructuralVariantBoundary checkForBoundaryBreak(
    		LinkedList<Interval> averagePosQueue, LocusInfoQueue runningLiq){
    	
    	/* The break should be clear when comparing the sliding
    	 * window averages that do NOT include the actual break.
    	 * By comparing only values that are
    	 * StructuralVariantCaller.windowSize apart.
    	 */

    	Object[] averages = averagePosQueue.toArray();
    	
    	int j;
    	Interval intervalI, intervalJ;
    	for(int i = 0; i < averages.length - StructuralVariantCaller.WINDOW_SIZE; i++){
			j = i + StructuralVariantCaller.WINDOW_SIZE;
			
			intervalI = (Interval) averages[i];
			intervalJ = (Interval) averages[j];
			if(intervalI == null || intervalJ == null){
				continue;
			}
			else if(!LocusInfoQueue.intervalsOrderedNonOverlapping(intervalI, intervalJ)){
				return getBoundaryBreak(runningLiq);
			}
    	}
    	return null;
    }
    
    
    private StructuralVariantBoundary getBoundaryBreak(LocusInfoQueue runningLiq){
    	Iterator<LocusInfo> it = runningLiq.iterator();
    	LocusInfo prevLocus = null, currLocus;
    	ArrayList<IntervalMassTuple> imts;
    	while(it.hasNext()){
    		currLocus = it.next();
    		
    		if(prevLocus != null){
    			if(!LocusInfoQueue.intervalsOrderedNonOverlapping(
    					prevLocus.getTopRefIntervalMassTuple().interval,
    					currLocus.getTopRefIntervalMassTuple().interval)){
    				
    				Interval hgPre = prevLocus.getTopRefIntervalMassTuple().interval,
    						hgPost = currLocus.getTopRefIntervalMassTuple().interval,
    						sampPre = prevLocus.getLocusInterval(),
    						sampPost = currLocus.getLocusInterval();

    				return new StructuralVariantBoundary(hgPre, hgPost, sampPre, sampPost);
    			}
    			else{
    				
    				/* some may overlap, so check if there's a break at this
    				 * locus. 
    				 * 
    				 * TODO: It's unclear exactly where the boundary is when
    				 * they overlap. 
    				 */
    				imts = currLocus.getIntervalMassTuples();
    				if(imts.get(1) != null
    						&& !LocusInfoQueue.intervalsOrderedNonOverlapping(
    								imts.get(0).interval, imts.get(1).interval)){

						Interval hgPre = prevLocus.getTopRefIntervalMassTuple().interval,
    						hgPost = currLocus.getTopRefIntervalMassTuple().interval,
    						sampPre = prevLocus.getLocusInterval(),
    						sampPost = currLocus.getLocusInterval();
    					return new StructuralVariantBoundary(hgPre, hgPost, sampPre, sampPost);
    				}
    			}
    		}
    		prevLocus = currLocus;
    	}
    	throw new RuntimeException("ERROR: there must be a StructuralVariantBoundary here!");
    }
    
    
    private boolean crossesContigs(Iterator<LocusInfo> iter,
    		LocusInfoQueue runningLiq,
    		ArrayList<IntervalMassTuple> prevTopTwoHGRefLociIntervals,
    		ArrayList<IntervalMassTuple> currTopTwoHGRefLociIntervals,
    		Interval prevSampLocusInterval,
    		Interval currSampLocusInterval) throws IOException{
    	
    	IntervalMassTuple prev1 = prevTopTwoHGRefLociIntervals.get(0),
    			curr1 = currTopTwoHGRefLociIntervals.get(0);
    	
    	if(prev1 == null || curr1 == null){
    		throw new RuntimeException("ERROR: No primary intervals should be null");
    	}	
    	

    	/* if contigs for the top curr and prev intervals don't match,
    	 * see if there is enough evidence from secondary hits to justify
    	 * skipping this as a translocation.
    	 */
    	String curr1Contig = curr1.interval.getContig(),
    			prev1Contig = prev1.interval.getContig();
    	if(!curr1Contig.equals(prev1Contig)){
    		
    		/* Loop until we are convinced whether it is a 
    		 * cross-contig translocation
    		 */
    		LocusInfo locus;
    		int count = 0;
    		while(iter.hasNext() && count < StructuralVariantCaller.WINDOW_SIZE){
    			locus = iter.next();

    			currTopTwoHGRefLociIntervals =
    					generateHumanGenomeRefIntervalByMass(locus);
    			

				currSampLocusInterval =
					new Interval(locus.getSequenceName(), locus.getPosition(),
							locus.getPosition());
    						
				if(currTopTwoHGRefLociIntervals == null){
					
					/* keep prev up with curr while going through
					 * regions with too much coverage
					 */
					prevSampLocusInterval = currSampLocusInterval;
					prevTopTwoHGRefLociIntervals = currTopTwoHGRefLociIntervals;
					continue;
				}

				/* add the LocusMassIntervals to the locusInfo object and
				 * add the locus to the queue
				 */
				locus.addRefIntervalMassTuples(currTopTwoHGRefLociIntervals);
				runningLiq.add(locus);


				/* if we hit a locus where the contig for the top HGRef position
				 * by mass equals the prev1Contig, then assume it's just a blip. 
				 */
    			if(currTopTwoHGRefLociIntervals.get(0).interval.getContig().equals(prev1Contig)){
    				return false;
    			}
    			
    			prevSampLocusInterval = currSampLocusInterval;
				prevTopTwoHGRefLociIntervals = currTopTwoHGRefLociIntervals;
    		}
    		
    		/* We have continued StructuralVariantCaller.windowSize positions
    		 * and we still see the new contig. This must be a cross-contig
    		 * translocation.
    		 */
    		return true;
    	}
    	else{
    		return false;
    	}
    }
	
	
	/**
	 * 
	 */
//	private void writePendingVars(){
//		Iterator<Integer> it = pendingVarList.keySet().iterator();
//		int pos;
//		VariantContext var;
//		while(it.hasNext()){
//			pos = it.next();
//			var = pendingVarList.get(pos);
//			if(var != null){
//			logger.debug("Var: " + pendingVarList.get(pos).getStart());
//			vw.writeVariantToVCF(pendingVarList.get(pos));
//			}
//		}
//		
//		deletionVarList.clear();
//		sampDeletionPreList.clear();
//		sampDeletionPostList.clear();
//		pendingVarList.clear();
//	}
	
	
	/**
	 * Build a VariantContext object based on the sample's genome. Will return
	 * null if an 'N' is found in either the reference or alternate alleles.
	 * 
	 * @param header
	 * @param sampleContig: The sample's contig name to extract the alternate allele from
	 * @param sampVarStart: The position where the alternate allele starts
	 * @param sampVarEnd: The position where the alternate allele ends
	 * @param refContig: The reference contig name to extract the reference allele from
	 * @param refStart: The position where the reference allele starts
	 * @param refEnd: The position where the reference allele ends
	 * @return
	 */
	private VariantContext buildVariant(SAMFileHeader header,
			String sampleContig, int sampVarStart, int sampVarEnd,
			String refContig, int refStart, int refEnd){
		
//		if(refEnd > 249000000 || altEnd > 249000000){
//			logger.debug("Contig: " + contig);
//			logger.debug("\n#########################\nref pos: " + refStart + "-" + refEnd);
//			logger.debug("alt pos: " + altStart + "-" + altEnd + "\n###########################\n");
//			logger.debug("\n#########################\nsampRef pos: " + sampRefStart + "-" + sampRefEnd);
//			logger.debug(referenceReader.getSequence("1").length());
//			logger.debug(referenceReader.getSequenceDictionary().getSequence("1").getSequenceLength());
//		}

		/*
		 * TODO: 
		 * 1. Get correct genotypes based on percentages (90/10).
		 * This will require figuring out both haplotypes in the
		 * assembly
		 */
		ArrayList<Allele> alleles = new ArrayList<Allele>();
		Allele ref = Allele.create(
				hgRefReader.getSubsequenceAt(refContig, refStart, refEnd).getBases(), true);
		Allele alt = Allele.create(
				sampleRefReader.getSubsequenceAt(sampleContig, sampVarStart,
						sampVarEnd).getBases(), false);
		
		/* Ignore any regions containing "N"s. */
//		if(ref.getBaseString().contains("N") || alt.getBaseString().contains("N")){
//		if(ref.getBaseString().contains("N")){
//			logger.debug("ref: " + ref.getBaseString());
//			logger.debug("alt: " + alt.getBaseString());
//			logger.debug("\n#########################\nref pos: " + refStart + "-" + refEnd);
//			logger.debug("alt pos: " + altStart + "-" + altEnd + "\n###########################\n");
//			return null;
//		}
		
		if(refEnd - refStart > 10000 || sampVarEnd - sampVarStart > 10000){
			logger.debug("locus pos: " + sampVarStart);
			logger.debug("Ref size: " + (refEnd-refStart));
			logger.debug("Alt size: " + (sampVarEnd-sampVarStart));
		}

		logger.debug("\n#########################\nref pos: " + refStart + "-" + refEnd);
		logger.debug("alt pos: " + sampVarStart + "-" + sampVarEnd + "\n###########################\n");
		logger.debug("ref: " + ref.getBaseString());
		logger.debug("alt: " + alt.getBaseString());
		alleles.add(ref);
		alleles.add(alt);
		
		GenotypesContext gc = GenotypesContext.
				create(new GenotypeBuilder(header.getReadGroups().get(0).getSample(), alleles).make());
		
		/* The length between start and stop must reflect the length of 
		 * the ref allele, but start needs to reflect where it is in
		 * the sample's genome. The 'stop' will be:
		 * sampVarStart + (refEnd - refStart);
		 */
		return new VariantContextBuilder()
			.chr(sampleContig)
			.start(sampVarStart)
			.stop(sampVarStart + (refEnd - refStart))
			.alleles(alleles)
			.genotypes(gc)
			.make();
	}
	
	
//	
//    /**
//     * 
//     * Open a SAM/BAM file for reading and return the SamReader obj
//     * 
//     * @param samFile
//     * @param vs
//     * @return SamReader
//     */
//    private SamReader openSam(final File samFile, ValidationStringency vs) {
//    	
////    	System.setProperty("java.io.tmpdir", "");
//
//		final SamReaderFactory factory =
//		          SamReaderFactory.makeDefault()
//		              .enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS,
//		            		  SamReaderFactory.Option.VALIDATE_CRC_CHECKSUMS)
//		              .validationStringency(vs);
//
//        final SamReader reader = factory.open(samFile);
//        
//        return reader;
//    }
   

}
