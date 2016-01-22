/**
 * 
 */
package spike.tools.utilitybelt;

import java.math.BigDecimal;
import java.text.DecimalFormat;

import net.sourceforge.argparse4j.inf.ArgumentParser;

import org.apache.log4j.Logger;

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
    
}
