/**
 * 
 */
package spike;

import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;

import org.apache.log4j.Logger;

import spike.tools.structuralvariantcaller.StructuralVariantCallerEngine;

/**
 * @author markebbert
 *
 */
public class SpikeEngine implements Engine {

	private static Logger logger = Logger.getLogger(SpikeEngine.class);
	private static ArgumentParser parser;

	public SpikeEngine() {
		return;
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		SpikeEngine se = new SpikeEngine();
		
		parser = se.instantiateArgParser();
		Namespace parsedArgs;

		try {
			/* args[0] is the argument for VTCEngine specifying the tool to create */
			if(args.length == 0){
				printUsageHelpAndExit();
			}
			String[] spikeArgs = new String[] {args[0]};
			
			/* Remove args[0] so we can pass the rest to the appropriate tool */
			String[] toolArgs = removeElement(args, 0);

			/* Parse first argument. Must specify the tool. Rest are passed to appropriate tool */
			parsedArgs = parser.parseArgs(spikeArgs);
			
			/* Get string specifying appropriate tool */
			String toolName = (String) parsedArgs.get("ToolName");
			KnownTools tool = se.getTool(toolName);
			if(tool == null){
				throw new ArgumentParserException("Invalid tool specified: " + toolName, parser);
			}
			
			/* Determine which tool was specified and call it */
			if(tool == KnownTools.ALIGNER){
//				AlignerEngine ale = new AlignerEngine(toolArgs);
//				ale.align();
			}
			else if(tool == KnownTools.ASSEMBLER){
//                AssemblerEngine ase = new AssemblerEngine(toolArgs);
//                ase.assemble();
            }
			else if(tool == KnownTools.STRUCTURAL_VARIANT_CALLER){
				StructuralVariantCallerEngine svc = new StructuralVariantCallerEngine(toolArgs);
				svc.callVariants();
			}
			
		} catch (ArgumentParserException e) {
			printErrorUsageHelpAndExit(e);
		} catch (Exception e) {
			logger.error("Caught unexpected exception, something is very wrong!");
			e.printStackTrace();
		}
	}
	
	/**
	 * Remove element from array 'orig' by copying to a new array without the element
	 * @param orig
	 * @param element
	 * @return a new String[] without element
	 */
	private static String[] removeElement(String[] orig, int element){
	    String[] n = new String[orig.length - 1];
	    System.arraycopy(orig, 0, n, 0, element );
	    System.arraycopy(orig, element+1, n, element, orig.length - element-1);
	    return n;
	}
	
	private ArgumentParser instantiateArgParser(){
		ArgumentParser parser = ArgumentParsers.newArgumentParser("VTC", false, "-");
		parser.description("SPIKE was designed to perform alignment, assembly, and structural"
				+ " variant calling using reads from long-read technologies such"
				+ " as PacBio SMRT and Oxford Nanopore Technologies' MinIon.");
		parser.usage("java -jar spike.jar ToolName");
		
		parser.addArgument("ToolName").dest("ToolName")
				.help("Specify the tool to use. Available tools are: " + createToolCommandLineToString());
		
		return parser;
	}
	
	/**
	 * Determine of 'tool' exists in the Variant Tool Chest (VTC). Names are not case-sensitive.
	 * @param tool
	 * @return
	 */
	private KnownTools getTool(String tool){
		for(KnownTools t : KnownTools.values()){
			if(t.permittedCommandsContain(tool)){
				return t;
			}
		}
		return null;
	}
	
	private static void printErrorUsageHelpAndExit(Exception e){
		logger.error(e.getMessage());
		printUsageHelpAndExit();
	}
	
	private static void printUsageHelpAndExit(){
		parser.printUsage();
		parser.printHelp();
		System.exit(1);		
	}
	
	private static String createToolCommandLineToString(){
		StringBuilder sb = new StringBuilder();
		int count = 1;
		for(KnownTools t : KnownTools.values()){
			sb.append("\n" + Integer.toString(count) + ". " + t.toString());
			count++;
		}
		return sb.toString();
	}

}
