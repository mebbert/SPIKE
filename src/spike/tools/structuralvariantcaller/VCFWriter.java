/**
 * 
 */
package spike.tools.structuralvariantcaller;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFHeaderVersion;
import htsjdk.variant.vcf.VCFUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.LinkedHashSet;
import java.util.Set;

/**
 * @author markebbert
 *
 */
public class VCFWriter {

	private IndexedFastaSequenceFile indexedFastaSequenceFile;
	private VariantContextWriter out;
	private static final int BUFFER_SIZE = 0;

	/**
	 * @throws FileNotFoundException 
	 * 
	 */
	public VCFWriter(File file, File refDict, Set<String> sampleNames) throws FileNotFoundException {
		if(refDict == null){
			throw new RuntimeException("Received a 'null' SAMSequenceDictionary. Something is very wrong!");
		}

		indexedFastaSequenceFile = new IndexedFastaSequenceFile(refDict);
		SAMSequenceDictionary dict = indexedFastaSequenceFile.getSequenceDictionary();

		if(dict == null){
			throw new FileNotFoundException("The reference sequence specified ("
					+ refDict.getAbsolutePath() +
					") does not have the appropriate dictionary file. Please use"
					+ " Picard's CreateSequenceDictionary.jar to generate this file.");
		}

		VariantContextWriterBuilder builder = new VariantContextWriterBuilder()
	       .setReferenceDictionary(dict)
	       .setOption(Options.INDEX_ON_THE_FLY)
	       .setBuffer(BUFFER_SIZE);
		
		this.out = builder
	       .setOutputFile(file)
	       .build();
		
		this.out.writeHeader(this.generateBasicHeader(dict, sampleNames));
	}
	
	/**
	 * Generate a basic header for the VCF
	 * 
	 * @param refDict
	 */
	private VCFHeader generateBasicHeader(SAMSequenceDictionary refDict, Set<String> sampleNames){
		LinkedHashSet<VCFHeaderLine> headerLines = new LinkedHashSet<VCFHeaderLine>();
		
		/* Add the 'fileFormat' header line (must be first) */
		headerLines.add(new VCFHeaderLine(VCFHeaderVersion.VCF4_2.getFormatString(),
				VCFHeaderVersion.VCF4_2.getVersionString()));
		
		/* Format field must have at least one value. The Genotype in this case. */
		headerLines.add(new VCFFormatHeaderLine("GT", 1, VCFHeaderLineType.String, "Genotype"));
		
		/* Create contig header lines */
		headerLines.addAll(VCFUtils.makeContigHeaderLines(refDict, null));
		
		return new VCFHeader(headerLines, sampleNames);
	}
	
	/**
	 * Write a variant to the vcf. This method is buffered, so output may
	 * not be immediate.
	 * 
	 * @param file
	 * @param vp
	 * @param refDict
	 * @param repairHeader
	 * @throws FileNotFoundException
	 */
	public void writeVariantToVCF(VariantContext var){
		this.out.add(var);
	}

}
