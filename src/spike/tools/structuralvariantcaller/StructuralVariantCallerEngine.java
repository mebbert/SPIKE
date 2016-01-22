/**
 * 
 */
package spike.tools.structuralvariantcaller;

import htsjdk.samtools.ValidationStringency;

import java.io.File;

import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.inf.ArgumentGroup;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;

import org.apache.log4j.Logger;

import spike.Engine;
import spike.tools.utilitybelt.UtilityBelt;

/**
 * @author markebbert
 *
 */
public class StructuralVariantCallerEngine implements Engine {

	private static Logger logger = Logger.getLogger(StructuralVariantCallerEngine.class);
	private static ArgumentParser parser;
	private Namespace parsedArgs;
	
	public StructuralVariantCallerEngine(String[] args) {
		init(args);
	}
	
	private void init(String[] args){
		parser = ArgumentParsers.newArgumentParser("StructuralVariantCaller");
		parser.description("The Structural Variant Caller will identify"
				+ " structural variants using aligned long reads from"
				+ " long-read technologies such as PacBio SMRT or the"
				+ " Oxford Nanopore Technologies' MinIon.");
		parser.defaultHelp(true);
		
		ArgumentGroup ioOptions = parser.addArgumentGroup("input/output arguments");
		ArgumentGroup svcOptions = parser.addArgumentGroup("variant caller arguments");
		
		/* Setup IO options */
		ioOptions
				.addArgument("-i", "--input")
				.dest("SAM")
				.metavar("SAM/BAM")
				.type(String.class)
				.required(true)
				.help("Specify the input file. This can be a SAM or BAM file.");
		
		ioOptions
				.addArgument("-o", "--output")
				.dest("VCF")
				.type(String.class)
				.required(true)
				.help("Specify the output VCF file.");
		
		
		/* Setup SVC options */
		svcOptions
				.addArgument("-s", "--min-structural-variant-size")
				.dest("MIN_SIZE")
				.metavar("SIZE")
				.setDefault(50)
				.type(Integer.class)
				.help("Set the minimum size structural variant to consider.");

		svcOptions
				.addArgument("-m", "--min-mapping-quality")
				.dest("MinMapQual")
				.metavar("QUAL")
				.setDefault(30)
				.type(Integer.class)
				.help("Set the minimum mapping quality to consider a given"
						+ " read.");
			
		svcOptions
				.addArgument("-v", "--validation-stringency")
				.dest("STRINGENCY")
				.setDefault("STRICT")
				.choices("STRICT", "LENIENT", "SILENT")
				.type(String.class)
				.help("Set the validation stringency when parsing the SAM/BAM"
						+ " file. 'STRICT' will throw errors if something "
						+ " is amiss, 'LENIENT' will give warnings but continue,"
						+ " and 'SILENT' will continue and keep our mouth shut.");
		
		
		
		try{
			parsedArgs = parser.parseArgs(args);
		} catch (ArgumentParserException e){
			parser.handleError(e);
			System.exit(1);
		}

	}
	
	public void callVariants(){
		
		String sam = parsedArgs.getString("SAM");
		String vcf = parsedArgs.getString("VCF");

		int minSVSize = parsedArgs.getInt("MIN_SIZE");
		int minMapQual = parsedArgs.getInt("MinMapQual");
		String stringency = parsedArgs.getString("STRINGENCY");
		ValidationStringency vs = null;
		
		if("strict".equalsIgnoreCase(stringency)){
			vs = ValidationStringency.STRICT;
		}
		else if("lenient".equalsIgnoreCase(stringency)){
			vs = ValidationStringency.LENIENT;
		}
		else if("silent".equalsIgnoreCase(stringency)){
			vs = ValidationStringency.SILENT;
		}
		
		// Start doing your thing.
		StructuralVariantCaller svc = new StructuralVariantCaller();
		try {
			svc.startWalking(new File(sam), minSVSize, minMapQual, vs);
		} catch (StructuralVariantCallerException e) {
			UtilityBelt.printErrorUsageHelpAndExit(parser, logger, e);
		}

	}

}
