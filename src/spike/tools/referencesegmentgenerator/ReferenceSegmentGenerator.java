/**
 * 
 */
package spike.tools.referencesegmentgenerator;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.fastq.FastqWriter;
import htsjdk.samtools.fastq.FastqWriterFactory;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.UnsupportedEncodingException;

import spike.tools.utilitybelt.UtilityBelt;

/**
 * @author markebbert
 *
 */
public class ReferenceSegmentGenerator {

	private FastqWriter acceptableReadsWriter, unAcceptableReadsWriter;

	public ReferenceSegmentGenerator(String fastqOut) throws UnsupportedEncodingException, FileNotFoundException{
		acceptableReadsWriter = new FastqWriterFactory().newWriter(new File(fastqOut));
		unAcceptableReadsWriter = new FastqWriterFactory().newWriter(new File("unacceptable_reads.fastq"));
	}

	public void generate(String refFasta, int segmentLength){
		
		// generate the initial reads, then optimize
//		optimize();
	}
	
	public void optimize(String refFastqOrSamBam, int minMapQ, int incrementLengthBy,
			int maxSegmentLength){
		
		// align the fastqs
		File inputFile = new File(refFastqOrSamBam);
		
		/* Determine if it's a FASTQ or SAM/BAM */
		String extension = "";
		int i = inputFile.getPath().lastIndexOf('.');
		if (i > 0) {
		    extension = inputFile.getPath().substring(i+1);
		}

		if(extension.equalsIgnoreCase(".fastq")){
			
			// perform alignment, then proceed
			// align()
			// optimize(bam)
		}
		else if(extension.equalsIgnoreCase(".bam")
					|| extension.equalsIgnoreCase(".sam")){
//			optimize(bam);
		}

	}
	
	/**
	 * Separate "acceptable" and "unacceptable" reads, accordingly. Continue
	 * optimizing the "unacceptable" reads until done.
	 * 
	 * @param samFile
	 * @param minMapQ
	 * @param incrementLengthBy
	 * @param maxSegmentLength
	 */
	private void optimize(File samFile, int minMapQ, int incrementLengthBy,
			int maxSegmentLength){
		
			/* Sort "acceptable" and "unacceptable" reads, increment
			 * the "unacceptable" and realign. Continue doing it
			 * until done.
			 */
			while(true){
				// sortReads
			}
	}
	
	private void sortReads(File samFile, int minMapQ, int incrementLengthBy,
			int maxSegmentLength){
		
		SamReader reader = UtilityBelt.openSam(samFile, ValidationStringency.SILENT);
		SAMRecordIterator samIt = reader.iterator();
		
		SAMRecord rec;
		while(samIt.hasNext()){
			rec = samIt.next();
			
			if(rec.getMappingQuality() < minMapQ
					&& rec.getReadLength() < maxSegmentLength){

				// increment the read and write to a new fastq
			}
			else{
//				FastqRecord fr = new FastqRecord(seqHeaderPrefix, seqLine, qualHeaderPrefix, qualLine)
//				acceptableReadsWriter.write(rec);
			}
		}
	}
}
