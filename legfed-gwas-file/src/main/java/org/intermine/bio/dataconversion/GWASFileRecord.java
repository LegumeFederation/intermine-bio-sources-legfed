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
 * Encapsulates a single tab-delimited GWAS experiment file record.
 *
 * 0          1                   2           3        4          5        6
 * phenotype  ontology_identifier marker      p_value  chromosome start    end
 * Seed oil   SOY:0001668         ss715591649 1.12E-09 Gm05       41780982 41780982
 *
 * @author Sam Hokin, NCGR
 */
public class GWASFileRecord implements Comparable {

    // file record values
    String phenotype;          // 0
    String ontologyIdentifier; // 1
    String marker;             // 2
    double pValue;             // 3
    String chromosome;         // 4
    int start;                 // 5
    int end;                   // 6

    // set this based on start,end
    String type;

    /**
     * Instantiate from a line in the GWAS file.
     */
    public GWASFileRecord(String line) {
        String[] parts = line.split("\t");
        if (parts.length!=7) {
            System.err.println("ERROR: GWASFileRecord input does not have 7 parts:"+line);
            System.exit(1);
        }
        phenotype = parts[0];
        ontologyIdentifier = parts[1];
        marker = parts[2];
        if (parts[3]!=null && parts[3].length()>0) pValue = Double.parseDouble(parts[3]);
        chromosome = parts[4];
        start = Integer.parseInt(parts[5]);
        end = Integer.parseInt(parts[6]);
        // set type based on length
        if (start==end) {
            type = "SNP";
        } else {
            type = "SSR";
        }
    }

    /**
     * For alpha sorting on chromosome and start
     */
    public int compareTo(Object o) {
        GWASFileRecord that = (GWASFileRecord)o;
        if (this.chromosome.equals(that.chromosome)) {
            return this.start - that.start;
        } else {
            return this.chromosome.compareTo(that.chromosome);
        }
    }

    /**
     * For diagnostics
     */
    public String toString() {
        String str = "";
        str += "phenotype="+phenotype;
        str += "ontologyIdentifier="+ontologyIdentifier;
        str += "marker="+marker;
        str += "pValue="+pValue;
        str += "chromosome="+chromosome;
        str += "start="+start;
        str += "end="+end;
        return str;
    }
}
