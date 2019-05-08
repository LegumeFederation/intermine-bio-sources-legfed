package org.intermine.bio.dataconversion;

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
 * Encapsulates a single tab-delimited LIS data store annot_info file record.
 *
 * Format:
 * pacId locusName transcriptName peptideName Pfam Panther KOG ec KO GO Best-hit-arabi-name arabi-symbol arabi-defline
 * 37170591 Phvul.001G000400 Phvul.001G000400.1 Phvul.001G000400.1.p PF00504 PTHR21649,PTHR21649:SF24 1.10.3.9 K14172 GO:0016020,GO:0009765 AT1G76570.1	Chlorophyll A-B binding family protein
 *
 * Assume that any of the descriptors (Pfam, etc.) can have multiple comma-separated values.
 *
 * @author Sam Hokin
 */
public class AnnotInfoRecord {
    
    String pacId;
    String locusName;
    String transcriptName;
    String peptideName;
    String[] Pfam;    // optional
    String[] Panther; // optional
    String[] KOG;     // optional
    String[] ec;      // optional
    String[] KO;      // optional
    String[] GO;      // optional
    String bestHitAtName;    // optional
    String bestHitAtSymbol;  // optional
    String bestHitAtDefline; // optional

    /**
     * Instantiate from a line from an LIS data store annot_nfo file. Do nothing if it's a header line.
     */
    public AnnotInfoRecord(String line) {
        if (!line.startsWith("#")) {
            try {
                // parse line
                String[] parts = line.split("\t");
                pacId = parts[0];
                locusName = parts[1];
                transcriptName = parts[2];
                peptideName = parts[3];
                if (parts.length>4) Pfam = parts[4].split(",");
                if (parts.length>5) Panther = parts[5].split(",");
                if (parts.length>6) KOG = parts[6].split(",");
                if (parts.length>7) ec = parts[7].split(",");
                if (parts.length>8) KO = parts[8].split(",");
                if (parts.length>9) GO = parts[9].split(",");
                if (parts.length>10) bestHitAtName = parts[10];
                if (parts.length>11) bestHitAtSymbol = parts[11];
                if (parts.length>12) bestHitAtDefline = parts[12];
            } catch (Exception ex) {
                throw new RuntimeException("Error parsing AnnotInfo line:\n" +
                                           line+"\n" +
                                           " pacId="+pacId +
                                           " locusName="+locusName +
                                           " transcriptName="+transcriptName +
                                           " peptideName="+peptideName +
                                           " Pfam="+Pfam +
                                           " Panther="+Panther +
                                           " KOG="+KOG +
                                           " ec="+ec +
                                           " KO="+KO +
                                           " GO="+GO +
                                           " bestHitAtName="+bestHitAtName +
                                           " bestHitAtSymbol="+bestHitAtSymbol +
                                           " bestHitAtDefline="+bestHitAtDefline);
            }
        }

    }
}
