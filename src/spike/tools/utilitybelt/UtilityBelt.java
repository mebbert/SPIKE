/**
 * 
 */
package spike.tools.utilitybelt;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;

import java.io.File;
import java.math.BigDecimal;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import net.sourceforge.argparse4j.inf.ArgumentParser;

import org.apache.log4j.Logger;

import spike.datastructures.AlignmentBlock;

/**
 * @author markebbert
 *
 */
public class UtilityBelt {

	public UtilityBelt(){
		return;
	}
	
	 public static String roundDoubleToString(double d) {
	        DecimalFormat df = new DecimalFormat("#.##");
	        return String.valueOf(df.format(d));
    }
	 
    /**
     * @param unrounded
     * @param precision
     * @param roundingMode
     * @return
     */
    public static double round(double unrounded, int precision, int roundingMode)
    {
        BigDecimal bd = new BigDecimal(unrounded);
        BigDecimal rounded = bd.setScale(precision, roundingMode);
        return rounded.doubleValue();
    }
	
	
	/**
	 * Print the error. Then print the usage and help
	 * information and exit
	 * @param e
	 */
	public static void printErrorUsageHelpAndExit(ArgumentParser parser, Logger logger, Exception e){
		System.err.println("\nERROR: " + e.getMessage() + "\n");
//		logger.error(e.getMessage());
		printUsageHelpAndExit(parser);
	}
	
	/**
	 * Print the error. Then print the usage
	 * information and exit
	 * @param e
	 */
	public static void printErrorUsageAndExit(ArgumentParser parser, Logger logger, Exception e){
		System.err.println("\nERROR: " + e.getMessage() + "\n");
		parser.printUsage();
		System.exit(1);
	}
	
	/**
	 * Print only the usage and help information and exit.
	 */
	public static void printUsageHelpAndExit(ArgumentParser parser){
		parser.printUsage();
		parser.printHelp();
		System.exit(1);		
	}
	
    /**
     * Taken from SAMUtils, modified by markebbert to include insertions
     * 
     * Given a Cigar, Returns blocks of the sequence that have been aligned directly to the
     * reference sequence. Note that clipped portions and deleted bases (vs. the reference)
     * are not represented in the alignment blocks.
     *
     * @param cigar The cigar containing the alignment information
     * @param alignmentStart The start (1-based) of the alignment
     * @param cigarTypeName The type of cigar passed - for error logging.
     * @return List of alignment blocks
     */
    public static List<AlignmentBlock> getAlignmentBlocks(final Cigar cigar,
    		final int alignmentStart, final String cigarTypeName,
    		final boolean includeDeletions) {

        if (cigar == null) return Collections.emptyList();

        final List<AlignmentBlock> alignmentBlocks = new ArrayList<AlignmentBlock>();
        int readBase = 1;
        int refBase = alignmentStart;
		int length;

        for (final CigarElement e : cigar.getCigarElements()) {
            switch (e.getOperator()) {
                case H:
                    break; // ignore hard clips
                case P:
                    break; // ignore pads
                case S:
                    readBase += e.getLength();
                    break; // soft clip read bases
                case N:
                    refBase += e.getLength();
                    break;  // reference skip
                case I:
					readBase += e.getLength();
					break;
                case D:
                	if(!includeDeletions){
						refBase += e.getLength();
                	}
                	else{
						length = e.getLength();
						alignmentBlocks.add(new AlignmentBlock(readBase, refBase, length, CigarOperator.D));
						refBase += length;
                	}
                    break;
                case M:
                    length = e.getLength();
                    alignmentBlocks.add(new AlignmentBlock(readBase, refBase, length, CigarOperator.M));
                    readBase += length;
                    refBase += length;
                    break;
                case EQ:
                    length = e.getLength();
                    alignmentBlocks.add(new AlignmentBlock(readBase, refBase, length, CigarOperator.EQ));
                    readBase += length;
                    refBase += length;
                    break;
                case X:
                    length = e.getLength();
                    alignmentBlocks.add(new AlignmentBlock(readBase, refBase, length, CigarOperator.X));
                    readBase += length;
                    refBase += length;
                    break;
                default:
                    throw new IllegalStateException("Case statement didn't deal with " + cigarTypeName + " op: " + e.getOperator());
            }
        }
        return Collections.unmodifiableList(alignmentBlocks);
    }
    
	/**
	 * 
	 * Open a SAM/BAM file for reading and return the SamReader obj
	 * 
	 * @param samFile
	 * @param vs
	 * @return SamReader
	 */
	 public static SamReader openSam(final File samFile, ValidationStringency vs) {
		 
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
