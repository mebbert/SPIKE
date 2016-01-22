package spike.datastructures;

import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class ReferenceSegment {
	
	private int chromosome;
	private int startPos;
	private int endPos;

	public ReferenceSegment(String referenceReadName) {
		parseReferenceReadName(referenceReadName);
	}

	/**
	 * Parse the chromosome, start, and end positions out of the
	 * reference segment's name.
	 * 
	 * @param referenceReadName
	 */
	private void parseReferenceReadName(String referenceReadName){
		Pattern p = Pattern.compile("(\\d+)\\[(\\d+)-(\\d+)\\]");
		Matcher m = p.matcher(referenceReadName);

		if (m.find()) {
			this.chromosome = Integer.parseInt(m.group(1));
			this.startPos = Integer.parseInt(m.group(2));
			this.endPos = Integer.parseInt(m.group(3));
		}
	}
	
	public int getChrom(){
		return this.chromosome;
	}
	
	public int getStartPos(){
		return this.startPos;
	}
	
	public int getEndPos(){
		return this.endPos;
	}
}
