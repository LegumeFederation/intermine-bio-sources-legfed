package org.intermine.bio.dataconversion;

import org.apache.log4j.Logger;

/*
 * Copyright (C) 2015-2016 NCGR
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  See the LICENSE file for more
 * information or http://www.gnu.org/copyleft/lesser.html.
 *
 */

/**
 * Encapsulates a single tab-delimited marker-chromosome file record.
 *
 * primaryIdentifier secondaryIdentifier Type Chromosome Start End Motif
 *
 * @author Sam Hokin, NCGR
 */
public class MarkerChromosomeRecord implements Comparable {

    private static final Logger LOG = Logger.getLogger(MarkerChromosomeRecord.class);

    String primaryIdentifier;
    String secondaryIdentifier;
    String type;
    String chromosome;
    int start;
    int end;
    String motif;

    /**
     * Instantiate from a line from a QTL file. Do nothing if it's a comment.
     */
    public MarkerChromosomeRecord(String line) {
        String[] parts = line.split("\t");
        primaryIdentifier = parts[0];
        secondaryIdentifier = parts[1];
        type = parts[2];
        chromosome = parts[3];
        start = Integer.parseInt(parts[4]);
        end = Integer.parseInt(parts[5]);
        motif = parts[6];
    }

    /**
     * For alpha sorting on chromosome and begin
     */
    public int compareTo(Object o) {
        MarkerChromosomeRecord that = (MarkerChromosomeRecord)o;
        if (this.primaryIdentifier.equals(that.primaryIdentifier) || this.secondaryIdentifier.equals(that.secondaryIdentifier)) {
            return this.start - that.start;
        } else {
            return this.primaryIdentifier.compareTo(that.primaryIdentifier);
        }
    }

}
