package org.intermine.bio.dataconversion;

/*
 * Copyright (C) 2015-2019 NCGR
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  See the LICENSE file for more
 * information or http://www.gnu.org/copyleft/lesser.html.
 *
 */

/**
 * Encapsulates a single tab-delimited GWAS file record from a gwas.txt file.
 *
 * 0                    1       2                 3       4       5                   6 
 * CHR                  BP      MARKER            PVAL    BPEND   PHENOTYPE           ONTOLOGY_IDENTIFIER
 * glyma.Wm82.gnm2.Gm02 2380973 glyma.ss107912620 3.27E-7 2380973 Stem diameter, main SOY:0001624
 *
 * @author Sam Hokin
 */
public class GwasRecord {
    String chr;
    int bp;
    String marker;
    double pval;
    int bpend;
    String phenotype;
    String ontologyIdentifier;

    /**
     * Instantiate from a line from an LIS data store gwas.txt file. Do nothing if it's a header line.
     */
    public GwasRecord(String line) {
        String[] parts = line.split("\t");
        chr = parts[0];
        bp = Integer.parseInt(parts[1]);
        marker = parts[2];
        if (parts.length>3 && parts[3].length()>0) pval = Double.parseDouble(parts[3]);
        if (parts.length>4 && parts[4].length()>0) bpend = Integer.parseInt(parts[4]);
        if (parts.length>5 && parts[5].length()>0) phenotype = parts[5];
        if (parts.length>6 && parts[6].length()>0) ontologyIdentifier = parts[6];
    }
}
