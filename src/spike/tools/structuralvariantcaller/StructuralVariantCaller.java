/**
 * 
 */
package spike.tools.structuralvariantcaller;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.filter.MappingQualityFilter;
import htsjdk.samtools.filter.SamRecordFilter;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.Interval;
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

import org.apache.log4j.Logger;

import spike.datastructures.LocusInfoQueue;
import spike.datastructures.OverclippedReadFilter;
import spike.datastructures.SamLocusIterator;
import spike.datastructures.SamLocusIterator.LocusInfo;

/**
 * @author markebbert
 *
 */
public class StructuralVariantCaller /*implements Runnable*/ {
	
	private static Logger logger = Logger.getLogger(StructuralVariantCaller.class);
	private static int minDepth, minSVSize, minMapQual, maxAcceptableCoverage,
						 maxAcceptableClip;
	private static boolean filterSingleEndClips = true;
	private static ValidationStringency samValidationStringency;
	
//	private int startWalking, endWalking;
	

	private LinkedList<Interval> sampPreList, sampPostList, hgPreList, hgPostList,
				sampDeletionPreList, sampDeletionPostList, hgDeletionPreList,
				hgDeletionPostList;
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
			final int minDepth, final int maxAcceptableCoverage,
			final int maxAcceptableClip,
			final ValidationStringency vs/*, int startWalking,
			int endWalking*/) throws FileNotFoundException {
		
//		this.startWalking = startWalking;
//		this.endWalking = endWalking;
		
		
		StructuralVariantCaller.minSVSize = minSVSize;
		StructuralVariantCaller.minMapQual = minMapQual;
		StructuralVariantCaller.minDepth = minDepth;
		StructuralVariantCaller.samValidationStringency = vs;
		StructuralVariantCaller.maxAcceptableCoverage = maxAcceptableCoverage;
		StructuralVariantCaller.maxAcceptableClip = maxAcceptableClip;
		
		this.sampPreList = new LinkedList<Interval>();
		this.sampPostList = new LinkedList<Interval>();
		this.hgPreList = new LinkedList<Interval>();
		this.hgPostList = new LinkedList<Interval>();
		this.sampDeletionPreList = new LinkedList<Interval>();
		this.sampDeletionPostList = new LinkedList<Interval>();
		this.hgDeletionPreList = new LinkedList<Interval>();
		this.hgDeletionPostList = new LinkedList<Interval>();
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
	 * All variants will be represented in a VCF as follows:
	 * 1. Contig name: set by the sample's contig name
	 * 2. Position: is the position on the sample's contig where
	 *              the variant is observed.
	 * 3. Ref: Is the sequence from the HG reference
	 * 4. Alt: The sequence observed on the sample
	 * 5. For translocations, inversions/transversions, and copy gains, the
	 *    INFO field will also contain where in the HG Ref it
	 *    came from.
	 * 
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
		

		/* Filter reads with a mapQ < StructuralVariantCaller.minMapQual,
		 * overclipped reads,
		 * and reads aligning to the negative strand
		 */
		List<SamRecordFilter> samFilters = Arrays.asList(
				new MappingQualityFilter(StructuralVariantCaller.minMapQual),
				new OverclippedReadFilter(maxAcceptableClip, filterSingleEndClips),
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
				refVarStart, refVarEnd, sampVarStart, sampVarEnd,
				sampRefLocusIntervalCovLost,
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
        	
        	if(locus.getPosition() == 104036){
        		logger.debug("here");
        	}

			currTopTwoHGRefLociIntervals = LocusInfoQueue.generateHumanGenomeRefIntervalByMass(locus);

//			logger.debug("prev1: " + prevTopTwoHGRefLociIntervals.get(0));
//			logger.debug("prev2: " + prevTopTwoHGRefLociIntervals.get(1));
//
//			logger.debug("curr1: " + currTopTwoHGRefLociIntervals.get(0));
//			logger.debug("curr2: " + currTopTwoHGRefLociIntervals.get(1));

			if(currTopTwoHGRefLociIntervals == null ||
					prevTopTwoHGRefLociIntervals == null) continue;

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

			/* Get depth at this position. If there is excessive coverage,
			 * we can't trust the data in this region.
			 */
			int depth = locus.getRecordAndPositions().size();

		
			/* If out of order, determine whether it's a deletion or
			 * translocation.
			 */
			if(!inOrder && depth < maxAcceptableCoverage){
				
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
						

						/* For deletions, the alt's start and end location will
						 * be the nucleotide (on the sample's sequence) just
						 * before the deleted sequence (sampleVarStart).
						 * 
						 * The ref allele's start will be at the same location,
						 * and will end on the last nucleotide of the deleted
						 * sequence.
						 */
						var = buildVariant(this.header,
								preSamp.getContig(), // Samp contig
								preSamp.getEnd(),    // samp locus before the deletion
								preSamp.getEnd(), // also samp locus before the deletion
								preHGRef.getContig(),// Ref contig
								preHGRef.getEnd(), // first locus of deletion with HG ref coordinates
								postHGRef.getStart() - 1 // last locus of deletion
								);
							
						if(var != null){
							deletionVarList.add(var);
//							pendingVarList.put(preHGRef.getStart(), var);
							pendingVarList.put(preSamp.getStart(), var);
							
							hgPreList.add(preHGRef);
							hgPostList.add(postHGRef);
							sampPreList.add(preSamp);
							sampPostList.add(postSamp);

							hgDeletionPreList.add(preHGRef);
							hgDeletionPostList.add(postHGRef);	
							sampDeletionPreList.add(preSamp);
							sampDeletionPostList.add(postSamp);	
						}

					}
				}
				else{ // must be a translocation or copy gain
					
					logger.debug("PreList size: " + sampPreList.size());
					
					/* Loop backwards over the lists to see if we can
					 * resolve any boundaries
					 */
					for(int i = sampPreList.size() - 1; i >= 0; i--){
						
						if(LocusInfoQueue.intervalsOrderedNonOverlapping(
										this.sampPreList.get(i), postSamp)){

							/* See if preList.get(i) was previously treated
							 * as a deletion. If so, get rid of it
							 * 
							 * TODO: Make sure the Interval.class equals method
							 * is being called appropriately.
							 */
							int index = sampDeletionPreList.indexOf(sampPreList.get(i));
							if(index >= 0){
								pendingVarList.remove(deletionVarList.get(index).getStart());
								sampDeletionPreList.remove(index);
								sampDeletionPostList.remove(index);
								hgDeletionPreList.remove(index);
								hgDeletionPostList.remove(index);
								deletionVarList.remove(index);
							}
							
							
							/* Since we found the match, now loop forward from
							 * where we are, resolving all boundaries in between.
							 */
							for(int j = i + 1; j < sampPreList.size(); ){

								/* For translocations/copy gains, the ref allele's
								 * start and end location will
								 * be the nucleotide (on the HG ref's sequence) just
								 * before the inserted sequence (refVarStart).
								 * 
								 * The alt allele's start will be at the same location,
								 * and will end on the last nucleotide of the inserted
								 * sequence.
								 */

								/* Create and write the variant */
								refVarStart = this.hgPreList.get(i);
								refVarEnd = refVarStart;
								sampVarStart = this.sampPreList.get(i);
								sampVarEnd = this.sampPreList.get(j);
								var = buildVariant(this.header,
										sampVarStart.getContig(),
										sampVarStart.getStart(),
										sampVarEnd.getEnd(),
										refVarStart.getContig(),
										refVarStart.getEnd(),
										refVarEnd.getStart()
										);
//								this.vw.writeVariantToVCF(var);
								
								if(var != null){
									pendingVarList.put(sampVarStart.getStart(), var);
									
									/* remove this boundary from pre and post list */
									sampPreList.remove(i);
									sampPostList.remove(i);
									hgPreList.remove(i);
									hgPostList.remove(i);

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
					hgPreList.add(preHGRef);
					hgPostList.add(postHGRef);
					sampPreList.add(preSamp);
					sampPostList.add(postSamp);
				}
		
			}

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

				/* The position we dropped below acceptable coverage */
				hgRefLocusIntervalCovLost =
						LocusInfoQueue.generateHumanGenomeRefIntervalByMass(locus)
						.get(0); // We trust the first in this situation.
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
						ArrayList<Interval> topTwo = 
								LocusInfoQueue.generateHumanGenomeRefIntervalByMass(locus);

						if(topTwo != null){
							currHGRefLocusInterval = topTwo.get(0);

							depth = locus.getRecordAndPositions().size();
//							logger.debug("depth (" + locus.getPosition() + "): " + depth);
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
		sampDeletionPreList.clear();
		sampDeletionPostList.clear();
		pendingVarList.clear();
	}
	
	
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
   

}
