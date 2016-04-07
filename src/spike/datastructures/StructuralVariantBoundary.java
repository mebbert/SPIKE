/**
 * 
 */
package spike.datastructures;

import htsjdk.samtools.util.Interval;

/**
 * A StructuralVariantBoundary is defined as a point between
 * two nucleotides where a structural variant begins or ends.
 * The boundary is defined on both the sample's assembled
 * sequence (sampPre and sampPos) and on the human genome
 * reference (hgPre and hgPost). The 'pre' values are the
 * position before the boundary, and the 'post' values are
 * the position after the boundary.
 * 
 * @author markebbert
 *
 */
public class StructuralVariantBoundary {
	
	public final Interval hgPre, hgPost, sampPre, sampPost;

	public StructuralVariantBoundary(Interval hgPre, Interval hgPost,
			Interval sampPre, Interval sampPost){
		this.hgPre = hgPre;
		this.hgPost = hgPost;
		this.sampPre = sampPre;
		this.sampPost = sampPost;
	}
	
	public String toString(){
		
		 return "Ref Pre/post: " + hgPre + "/" + hgPost +
				"\nSamp Pre/post:" + sampPre + "/" + sampPost +
				"\n";
	}
}
