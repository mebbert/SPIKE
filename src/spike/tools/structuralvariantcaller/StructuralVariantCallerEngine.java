/**
 * 
 */
package spike.tools.structuralvariantcaller;

import net.sourceforge.argparse4j.inf.ArgumentParser;
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
		
	}
	
	public void callVariants(){
		
	}

}
