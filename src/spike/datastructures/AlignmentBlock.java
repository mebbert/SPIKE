package spike.datastructures;

import htsjdk.samtools.CigarOperator;

import java.io.Serializable;

/**
 * Represents the contiguous alignment of a subset of read bases to a reference
 * sequence. Simply put an alignment block tells you that read bases from
 * readStart are aligned to the reference (matching or mismatching) from
 * referenceStart for length bases.
 *
 * @author Tim Fennell
 */
public class AlignmentBlock implements Serializable {
    public static final long serialVersionUID = 1L;

    private int readStart;
    private int referenceStart;
    private int length;
    private CigarOperator blockType;

    /** Constructs a new alignment block with the supplied read and ref starts and length. */
    public AlignmentBlock(int readStart, int referenceStart, int length, CigarOperator blockType) {
        this.readStart = readStart;
        this.referenceStart = referenceStart;
        this.length = length;
        this.blockType = blockType;
    }

    /** The first, 1-based, base in the read that is aligned to the reference reference. */
    public int getReadStart() { return readStart; }

    /** The first, 1-based, position in the reference to which the read is aligned. */
    public int getReferenceStart() { return referenceStart; }

    /** The number of contiguous bases aligned to the reference. */
    public int getLength() { return length; }
    
    /** The type of block. */
    public CigarOperator getBlockType(){ return blockType; }
}
