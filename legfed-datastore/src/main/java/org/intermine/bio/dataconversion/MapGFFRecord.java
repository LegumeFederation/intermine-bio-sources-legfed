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
 * Encapsulates a single genetic map GFF file record.
 *
 * 0                       1       2       3       4       5       6       7       8
 * phavu.G19833.gnm1.Chr01 blastn  SNP     242423  242423  .       +       .       Name=Pv_TOG905303_749;ID=1;Note=LG01 cM 0.0;alleles=T/G
 *
 * @author Sam Hokin, NCGR
 */
public class MapGFFRecord {

    String chr;
    String source;
    String type;
    int start;
    int end;
    String strand;
    String name;
    String fullname;
    String alleles;

    /**
     * Instantiate from a line from a CMap file. Do nothing if it's a header line.
     */
    public MapGFFRecord(String line) {
        if (!line.startsWith("#")) {
            try {
                // parse line
                String[] parts = line.split("\t");
                chr = parts[0];
                source = parts[1];
                type = parts[2];
                start = Integer.parseInt(parts[3]);
                end = Integer.parseInt(parts[4]);
                strand = parts[6];
                String attributeString = parts[8];
                String[] attributes = attributeString.split(";");
                for (String attribute : attributes) {
                    String[] equalparts = attribute.split("=");
                    if (equalparts[0].equals("Name")) {
                        fullname = equalparts[1];
                        name = fullname;
                        String[] underscoreparts = fullname.split("_");
                        if (underscoreparts.length==3) {
                            // strip front and back of Pv_TOG905303_749
                            name = underscoreparts[1];
                        } else if (underscoreparts.length==2) {
                            // strip back of TOG905303_749
                            name = underscoreparts[0];
                        }
                    } else if (equalparts[0].equals("alleles")) {
                        alleles = equalparts[1];
                    }
                }
            } catch (Exception ex) {
                throw new RuntimeException("Error parsing Map GFF line:\n"+line);
            }
        }
    }

    /**
     * Return true if this record actually has data.
     */
    boolean hasData() {
        return name!=null;
    }

    /**
     * Return true if record is a SNP marker.
     */
    boolean isSNP() {
        return type!=null && type.equals("SNP");
    }
}
