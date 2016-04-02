/*
 * The MIT License
 *
 * Copyright (c) 2015 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

package spike.tools.structuralvariantcaller;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.filter.SamRecordFilter;

import java.util.List;

import org.apache.log4j.Logger;

/**
 * Filters out reads with very few unclipped bases, likely due to the read coming
 * from a foreign organism, e.g. bacterial contamination.
 *
 * Based on GATK's OverclippedReadFilter.
 * 
 * @author modified by markebbert
 */
public class OverclippedReadFilter implements SamRecordFilter {
	
	private static Logger logger = Logger.getLogger(OverclippedReadFilter.class);
	
    // if the number of clipped bases is above this threshold, the read is considered overclipped
    private final int clippedBasesThreshold;
    // if set to true, then reads with at least one clipped end will be filtered; if false, we require both ends to be clipped
    private final boolean filterSingleEndClips;

    public OverclippedReadFilter(final int clippedBasesThreshold, final boolean filterSingleEndClips) {
        if (clippedBasesThreshold < 0) throw new SAMException("clippedBasesThreshold must be non-negative");
        this.clippedBasesThreshold = clippedBasesThreshold;
        this.filterSingleEndClips = filterSingleEndClips;
    }

    @Override
    public boolean filterOut(final SAMRecord record) {
    	int frontClippedLength = 0;
    	int backClippedLength = 0;
        CigarElement element;
        CigarOperator lastOperator = null;

		List<CigarElement> cigarElements = record.getCigar().getCigarElements();
		int nElements = cigarElements.size();
        for(int i = 0; i < nElements; i++){

        	element = cigarElements.get(i);
            if ( element.getOperator() == CigarOperator.S ) {
            	
            	/* If the last operator is null, this must be the front
            	 * clip
            	 */
            	if(lastOperator == null){
            
					frontClippedLength += element.getLength();

					// Treat consecutive S blocks as a single one
            		while(i + 1 < nElements && 
            				cigarElements.get(i+1).getOperator() == CigarOperator.S){
            			element = cigarElements.get(++i);
						frontClippedLength += element.getLength();
            		}
            	}
            	else if(lastOperator != CigarOperator.S){ // This must be the back clip

            		backClippedLength =+ element.getLength();

					// Treat consecutive S blocks as a single one
            		while(i + 1 < nElements &&
            				cigarElements.get(i+1).getOperator() == CigarOperator.S){
            			element = cigarElements.get(++i);
						backClippedLength += element.getLength();
            		}
            	}

            }
            lastOperator = element.getOperator();
        }
        
//        if(frontClippedLength >= clippedBasesThreshold || backClippedLength >= clippedBasesThreshold){
//        	
//        	logger.debug("filtering: " + record.getCigarString());
//        }
        if(filterSingleEndClips){
			return (frontClippedLength >= clippedBasesThreshold 
						|| backClippedLength >= clippedBasesThreshold);
        }
        else{
        	return (frontClippedLength >= clippedBasesThreshold 
						&& backClippedLength >= clippedBasesThreshold);
        }
    }

    @Override
    public boolean filterOut(final SAMRecord first, final SAMRecord second) {
        return filterOut(first) || filterOut(second);
    }
}
