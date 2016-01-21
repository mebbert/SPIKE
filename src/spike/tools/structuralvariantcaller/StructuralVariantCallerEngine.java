/**
 * 
 */
package spike.tools.structuralvariantcaller;

import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.inf.ArgumentGroup;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;

import org.apache.log4j.Logger;

import spike.Engine;

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
		
		ArgumentGroup ioOptions = parser.addArgumentGroup("input arguments");
		
		ioOptions
					.addArgument("-i", "--input")
					.dest("BAM")
					.metavar("SAM/BAM")
					.required(true)
					.help("Specify the input file. This can be a SAM or BAM file.");
		
		ioOptions
				.addArgument("-o", "--output")
				.dest("VCF")
				.required(true)
				.help("Specify the output VCF file.");
		
		try{
			parsedArgs = parser.parseArgs(args);
		} catch (ArgumentParserException e){
			parser.handleError(e);
			System.exit(1);
		}

	}
	
	public void callVariants(){
		
	}

}
