/**
 * 
 */
package spike.tools.structuralvariantcaller;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.TreeMap;
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
	

	private LinkedList<Interval> preList, postList,
				deletionPreList, deletionPostList;
	private LinkedList<VariantContext> deletionVarList;
	private TreeMap<Integer, VariantContext> pendingVarList;

	private IndexedFastaSequenceFile sampleRefReader, hgRefReader;
	private SAMFileHeader header;
	private SamReader reader;
	

	private VCFWriter vw;

	/**
	 * @throws FileNotFoundException 
	 * 
	 */
	public StructuralVariantCaller(final File samFile, final File outVCF,
			File sampleRef, File hgRef, final int minSVSize, final int minMapQual,
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
		this.deletionPreList = new LinkedList<Interval>();
		this.deletionPostList = new LinkedList<Interval>();
		this.deletionVarList = new LinkedList<VariantContext>();
		this.pendingVarList = new TreeMap<Integer, VariantContext>();
		
		this.sampleRefReader =
				new IndexedFastaSequenceFile(sampleRef);
		this.hgRefReader =
				new IndexedFastaSequenceFile(hgRef);

		this.reader = openSam(samFile,
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
		
		/* Set the base quality score cutoff to 0 
		 * so the iterator won't try to validate base qualities. Any
		 * secondary alignment won't have the qualities, and they're
		 * all made up anyway.
		 */
		int baseQualityScoreCutoff = 0;
		sli.setQualityScoreCutoff(baseQualityScoreCutoff);
		

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
		LocusInfoQueue running_liq =
				new LocusInfoQueue(StructuralVariantCaller.minSVSize, MAXSIZE),
				possibleVarLiq;
		LocusInfo locus;
		
		/* These locus Intervals are in context of the reference genome.
		 * We use the reference segment's original location (read name)
		 * to determine the locus interval.
		 */
		Interval /*prevHGRefLocusInterval = null,*/ currHGRefLocusInterval,
				prevSampLocusInterval = null, currSampLocusInterval,
				preHGRef, postHGRef, preSamp, postSamp,
				refStartInterval, refEndInterval, varStartInterval,
				varEndInterval, sampRefLocusIntervalCovLost,
				sampRefLocusIntervalCovGained, hgRefLocusIntervalCovLost,
				hgRefLocusIntervalCovGained;
		ArrayList<Interval> currTopTwoHGRefLociIntervals,
				prevTopTwoHGRefLociIntervals;

		boolean inOrder;


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
		}
		
		/* Get the first position and set to 'prev' */
		if(iter.hasNext()){
			locus = iter.next();
			prevTopTwoHGRefLociIntervals = LocusInfoQueue.generateHumanGenomeRefIntervalByMass(locus);
		}
		else{
			return;
		}

		while(iter.hasNext()){

			locus = iter.next();
        	
        	if(++count % 1000 == 0){
        		logger.debug("Assessed " + count + " loci.");
        	}
			if(locus.getPosition() == 10146){
				logger.debug("here");
			}				
			currTopTwoHGRefLociIntervals =
					LocusInfoQueue.generateHumanGenomeRefIntervalByMass(locus);

			currSampLocusInterval =
					new Interval(locus.getSequenceName(), locus.getPosition(),
							locus.getPosition());
			
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
			inOrder = LocusInfoQueue.
						intervalsOrderedNonOverlapping(prevTopTwoHGRefLociIntervals,
								currTopTwoHGRefLociIntervals);
		
			/* If out of order, determine whether it's a deletion or
			 * translocation.
			 */
			if(!inOrder){
				
				/* 'pre' and 'post' refer to the intervals on either side of the
				 * order-breaking boundary. HGRef is based on the human reference
				 * sequence position, and 'Samp' is based on the individual's
				 * actual contig. The 'Samp' is for reporting the variant
				 * according to the individual's position and 'HGRef' is so
				 * we know where in the reference genome everything associates
				 * with.
				 */
				preHGRef = prevTopTwoHGRefLociIntervals.get(0); // take the first if none were in order
				postHGRef = currTopTwoHGRefLociIntervals.get(0);
				
				preSamp = prevSampLocusInterval;
				postSamp = currSampLocusInterval;

				/* If pre is before post (based on original location in human
				 * reference genome), we'll assume it's a deletion, unless
				 * we find it's a translocation, later.
				 * 
				 * TODO: Keep looping until the next boundary and determine if it was
				 * a deletion, or something else.
				 */
				if(preHGRef.getContig().equals(postHGRef.getContig())
						&& preHGRef.getEnd() < postHGRef.getStart()){
					
					if(postHGRef.getStart() - preHGRef.getEnd()
							> StructuralVariantCaller.minSVSize){
						

						var = buildVariant(this.hgRefReader, this.header,
								preHGRef.getContig(),
								preHGRef.getStart(),
								postHGRef.getStart(),
								preHGRef.getStart(),
								preHGRef.getStart(),
								preSamp.getStart(),
								postSamp.getStart());

						if(var != null){
							deletionVarList.add(var);
//							pendingVarList.put(preHGRef.getStart(), var);
							pendingVarList.put(preSamp.getStart(), var);
							
//							preList.add(preHGRef);
//							postList.add(postHGRef);
							preList.add(preSamp);
							postList.add(postSamp);

//							deletionPreList.add(preHGRef);
//							deletionPostList.add(postHGRef);	
							deletionPreList.add(preSamp);
							deletionPostList.add(postSamp);	
						}

					}
				}
				else{
					
					logger.debug("PreList size: " + preList.size());
					
					/* Loop backwards over the lists to see if we can
					 * resolve any boundaries
					 */
					for(int i = preList.size() - 1; i >= 0; i--){
						
						if(LocusInfoQueue.intervalsOrderedNonOverlapping(
										this.preList.get(i), postHGRef)){

							/* See if preList.get(i) was previously treated
							 * as a deletion. If so, get rid of it
							 * 
							 * TODO: Make sure the Interval.class equals method
							 * is be called appropriately.
							 */
							int index = deletionPreList.indexOf(preList.get(i));
							if(index >= 0){
								pendingVarList.remove(deletionVarList.get(index).getStart());
								deletionPreList.remove(index);
								deletionPostList.remove(index);
								deletionVarList.remove(index);
							}
							
							
							/* Since we found the match, now loop forward from
							 * where we are, resolving all boundaries in between.
							 */
							for(int j = i + 1; j < preList.size(); ){

								/* Create and write the variant */
								refStartInterval = this.preList.get(i);
								refEndInterval = refStartInterval;
								varStartInterval = this.postList.get(i);
								varEndInterval = this.preList.get(j);
								var = buildVariant(this.sampleRefReader, this.header,
										refStartInterval.getContig(),
										refStartInterval.getStart(), 
										refEndInterval.getStart(),
										varStartInterval.getStart(),
										varEndInterval.getStart(),
										null, null);
//								this.vw.writeVariantToVCF(var);
								
								if(var != null){
									pendingVarList.put(refStartInterval.getStart(), var);
									
									/* remove this boundary from pre and post list */
									preList.remove(i);
									postList.remove(i);

									/* don't increment i or j, because we just
									 * removed the objects at i.
									 */
								}
								else{
									logger.debug("Var == null");
									logger.debug("i: " + i);
									logger.debug("j: " + j);
								}
								
							}
							
							
							break;
						}
					}

					/* If we didn't resolve this boundary with any existing,
					 * just add them.
					 */
//					preList.add(preHGRef);
//					postList.add(postHGRef);
					preList.add(preSamp);
					postList.add(postSamp);
				}
		
			}

			/* Get depth at this position */
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


				if(locus.getPosition() == 10398){
					logger.debug("here");
				}	
				/* The position we dropped below acceptable coverage */
				hgRefLocusIntervalCovLost =
						LocusInfoQueue.generateHumanGenomeRefIntervalByMass(locus)
						.get(0); // We trust the first in this situation.
				sampRefLocusIntervalCovLost = new Interval(
								locus.getSequenceName(), locus.getPosition(),
								locus.getPosition());

				/* May be null if there are multiple statistical modes */
				if(hgRefLocusIntervalCovLost != null){
				logger.debug("CovLost (" + locus.getPosition() + "): " + hgRefLocusIntervalCovLost.getStart());
				logger.debug("depth (" + locus.getPosition() + "): " + depth);


					/* Keep iterating until we rise above minDepth */
					while(iter.hasNext()){
	
						locus = iter.next();
						if(locus.getPosition() == 10399){
							logger.debug("here");
						}	
						currHGRefLocusInterval = 
								LocusInfoQueue.generateHumanGenomeRefIntervalByMass(locus)
								.get(0);

						if(currHGRefLocusInterval != null){

							depth = locus.getRecordAndPositions().size();
							logger.debug("depth (" + locus.getPosition() + "): " + depth);
							if(depth >= StructuralVariantCaller.minDepth){

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
									< StructuralVariantCaller.minSVSize){
									
									/* break if the structural variant wasn't
									 * large enough.
									 */
									break;
								}

								/* Create and write the variant */
								logger.debug("CovGained (" + locus.getPosition() + "): " + hgRefLocusIntervalCovGained.getStart() + "\n");
//								var = buildVariant(this.sampleReferenceReader, this.header,
//										hgRefLocusIntervalCovLost.getContig(),
//										hgRefLocusIntervalCovLost.getStart() - 1, // Ref start
//										hgRefLocusIntervalCovLost.getStart() - 1, // Ref end
//										hgRefLocusIntervalCovLost.getStart() - 1,     // Var start
//										hgRefLocusIntervalCovGained.getStart() - 1); // Var end
								var = buildVariant(this.sampleRefReader, this.header,
										sampRefLocusIntervalCovLost.getContig(),
										sampRefLocusIntervalCovLost.getStart() - 1, // Ref start
										sampRefLocusIntervalCovLost.getStart() - 1, // Ref end
										sampRefLocusIntervalCovLost.getStart() - 1,     // Var start
										sampRefLocusIntervalCovGained.getStart() - 1,
										null, null); // Var end

								if(var != null){
//									pendingVarList.put(hgRefLocusIntervalCovLost.getStart() - 1, var);
									pendingVarList.put(sampRefLocusIntervalCovLost.getStart() - 1, var);
								}

								/* Write any deletions that remain in the list */
								/* TODO: Make sure this is the appropriate place. Want
								 * to do it once we've seen the next insertion/inversion/duplication
								 */
								writePendingVars();
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

		sli.close();
	}
	
	
	/**
	 * 
	 */
	private void writePendingVars(){
		Iterator<Integer> it = pendingVarList.keySet().iterator();
		int pos;
		VariantContext var;
		while(it.hasNext()){
			pos = it.next();
			var = pendingVarList.get(pos);
			if(var != null){
			logger.debug("Var: " + pendingVarList.get(pos).getStart());
			vw.writeVariantToVCF(pendingVarList.get(pos));
			}
		}
		
		deletionVarList.clear();
		deletionPreList.clear();
		deletionPostList.clear();
		pendingVarList.clear();
	}
	
	
	/**
	 * Build a VariantContext object based on the sample's genome. Will return
	 * null if an 'N' is found in either the reference or alternate alleles.
	 * 
	 * @param referenceReader
	 * @param header
	 * @param contig
	 * @param refStart
	 * @param refEnd
	 * @param altStart
	 * @param altEnd
	 * @param sampRefStart
	 * @param sampRefEnd
	 * @return
	 */
	private VariantContext buildVariant(IndexedFastaSequenceFile referenceReader,
			SAMFileHeader header, String contig, int refStart, int refEnd,
			int altStart, int altEnd, Integer sampRefStart, Integer sampRefEnd){
		
		if(refEnd > 249000000 || altEnd > 249000000){
			logger.debug("\n#########################\nref pos: " + refStart + "-" + refEnd);
			logger.debug("alt pos: " + altStart + "-" + altEnd + "\n###########################\n");
			logger.debug("\n#########################\nsampRef pos: " + sampRefStart + "-" + sampRefEnd);
		}

		/*
		 * TODO: 
		 * 1. Build a VCF where position references the individual's position
		 * and put the corresponding HGRef positions in the INFO field.
		 * 
		 * 2. Get correct genotypes based on percentages (90/10).
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
//			logger.debug("ref: " + ref.getBaseString());
//			logger.debug("alt: " + alt.getBaseString());
			logger.debug("\n#########################\nref pos: " + refStart + "-" + refEnd);
			logger.debug("alt pos: " + altStart + "-" + altEnd + "\n###########################\n");
			return null;
		}
		alleles.add(ref);
		alleles.add(alt);
		
		GenotypesContext gc = GenotypesContext.
				create(new GenotypeBuilder(header.getReadGroups().get(0).getSample(), alleles).make());
		
		if(sampRefStart == null|| sampRefEnd == null){
			return new VariantContextBuilder()
				.chr(contig)
				.start(refStart)
				.stop(refEnd)
				.alleles(alleles)
				.genotypes(gc)
				.make();
		}
		else{
			return new VariantContextBuilder()
				.chr(contig)
				.start(sampRefStart)
				.stop(sampRefStart + (ref.length() - 1))
				.alleles(alleles)
				.genotypes(gc)
				.make();
		}
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
