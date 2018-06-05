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
 * Encapsulates a single tab-delimited GWAS experiment file record.
 *
 * Records:
 * gwas_id gwas_name         other_name gwas_family    gwas_class   trait_SOY_number QTL_type    QTL_category               locus_id locus_name  p_value  LOD  R2   chromosome linkage_group qtl_start qtl_end
 * 1       Seed protein 3-g1 Foobar     Seed protein 3 Seed protein SOY:0001676	     QTL_protein seed_composition_and_yield 119892   ss715637225 3.24E-06 0.87 0.67 Gm20       I	     30743655  30743655
 *
 * Note that the "qtl_start" and "qtl_end" are really the start, end of the associated marker, either a SNP (start=end) or a SSR (start<end). So those are stored with the marker, not the QTL.
 *
 * @author Sam Hokin, NCGR
 */
public class GWASFileRecord implements Comparable {

    private static final Logger LOG = Logger.getLogger(GWASFileRecord.class);

    // file record values
    int gwasId;
    String gwasName;
    String otherName;
    String gwasFamily;
    String gwasClass;
    String traitSOYNumber;
    String qtlType;
    String qtlCategory;
    int    locusId;
    String locusName;
    double pValue;
    double lod;
    double r2;
    String chromosome;
    String linkageGroup;
    int    start;
    int    end;

    // set this based on start,end
    String type;

    // gwas_id
    // gwas_name
    // other_name
    // gwas_family
    // gwas_class
    // trait_SOY_number
    // QTL_type
    // QTL_category
    // locus_id
    // locus_name
    // p_value
    // LOD
    // R2
    // chromosome
    // linkage_group
    // qtl_start
    // qtl_end
    
    /**
     * Instantiate from a line from a QTL file. Do nothing if it's a comment.
     */
    public GWASFileRecord(String line) {
        String[] parts = line.split("\t");
        int i = 0;
        gwasId = Integer.parseInt(parts[i++]);
        gwasName = parts[i++];
        otherName = parts[i++];
        gwasFamily = parts[i++];
        gwasClass = parts[i++];
        traitSOYNumber = parts[i++];
        qtlType = parts[i++];
        qtlCategory = parts[i++];
        if (parts[i].length()>0) locusId = Integer.parseInt(parts[i]); i++;
        locusName = parts[i++];
        if (parts[i].length()>0) pValue = Double.parseDouble(parts[i]); i++;
        if (parts[i].length()>0) lod = Double.parseDouble(parts[i]); i++;
        if (parts[i].length()>0) r2 = Double.parseDouble(parts[i]); i++;
        chromosome = parts[i++];
        linkageGroup = parts[i++];
        // these must exist, otherwise file needs to be fixed
        start = Integer.parseInt(parts[i++]);
        end = Integer.parseInt(parts[i++]);
        // set type
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
        str += "gwasId="+gwasId;
        str += "gwasName="+gwasName;
        str += "otherName="+otherName;
        str += "gwasFamily="+gwasFamily;
        str += "gwasClass="+gwasClass;
        str += "traitSOYNumber="+traitSOYNumber;
        str += "qtlType="+qtlType;
        str += "qtlCategory="+qtlCategory;
        str += "locusId="+locusId;
        str += "locusName="+locusName;
        str += "pValue="+pValue;
        str += "lod="+lod;
        str += "r2="+r2;
        str += "chromosome="+chromosome;
        str += "linkageGroup="+linkageGroup;
        str += "start="+start;
        str += "end="+end;
        return str;
    }

}
