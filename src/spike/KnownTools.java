/**
 * 
 */
package spike;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * @author markebbert
 *
 */
public enum KnownTools {

	REFERENCE_SEGMENT_GENERATOR("ReferenceSegementGenerator", "RSG", "Generate"
			+ " and optimize reference segments to align to an individual's"
			+ " assembled genome.",
			new ArrayList<String>(Arrays.asList(new String[]{"RSG", "ReferenceSegementGenerator"}))),
    ALIGNER("Aligner", "AL", "Align long-read technology reads",
			new ArrayList<String>(Arrays.asList(new String[]{"AL", "Aligner"}))),
    ASSEMBLER("Assembler", "AS", "Assemble an individual's genome",
            new ArrayList<String>(Arrays.asList(new String[]{"A", "Assembler"}))),
	STRUCTURAL_VARIANT_CALLER("StructuralVariantCaller", "SVC", "Identify structural variants using long-read technologies",
			new ArrayList<String>(Arrays.asList(new String[]{"SVC", "StructuralVariantCaller"})));
	
	private String name, shortCommand, briefDescription;
	private ArrayList<String> permittedCommands;
	private KnownTools(String name, String shortCommand, String briefDescription, ArrayList<String> permittedCommands){
		this.name = name;
		this.shortCommand = shortCommand;
		this.briefDescription = briefDescription;
		this.permittedCommands = permittedCommands;
	}
	
	public String getName(){
		return this.name;
	}
	
	public String getShortCommand(){
		return this.shortCommand;
	}
	
	public String getBriefDescription(){
		return this.briefDescription;
	}
	
	public ArrayList<String> getPermittedCommands(){
		return this.permittedCommands;
	}
	
	public boolean permittedCommandsContain(String command){
		for(String s : permittedCommands){
			if(s.equalsIgnoreCase(command)){
				return true;
			}
		}
		return false;
	}
	
	@Override
	public String toString(){
		return  getName() + " (" +
				getShortCommand() +
				") -- " + getBriefDescription() +
				".";
	}
}
