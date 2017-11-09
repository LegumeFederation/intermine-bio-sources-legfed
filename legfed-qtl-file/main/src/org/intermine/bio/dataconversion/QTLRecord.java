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
 * Encapsulates a single tab-delimited QTL file record.
 *
 * QTLName Parent_1 Parent_2 TraitName [PMID] [Journal Year Volume Page Title]
 *
 * @author Sam Hokin, NCGR
 */
public class QTLRecord implements Comparable {

    private static final Logger LOG = Logger.getLogger(QTLRecord.class);

    String qtlName;
    String parent1;
    String parent2;
    String traitName;
    int pubPMID;
    String pubJournal;
    String pubYear;
    String pubVolume;
    String pubPages;
    String pubTitle;

    /**
     * Instantiate from a line from a QTL file. Do nothing if it's a comment.
     */
    public QTLRecord(String line) {
        String[] parts = line.split("\t");
        qtlName = parts[0];
        parent1 = parts[1];
        parent2 = parts[2];
        traitName = parts[3];
        try {
            pubPMID = Integer.parseInt(parts[4]);
        } catch (Exception e) {
            // oh, well, no PMID or not one that's parseable
        }
        try {
            pubJournal = parts[5];
            pubYear = parts[6];
            pubVolume = parts[7];
            pubPages = parts[8];
            pubTitle = parts[9];
        } catch (Exception e) {
            // do nothing, short on journal columns
        }
    }

    /**
     * For alpha sorting on QTL name
     */
    public int compareTo(Object o) {
        QTLRecord that = (QTLRecord)o;
        return this.qtlName.compareTo(that.qtlName);
    }

}
