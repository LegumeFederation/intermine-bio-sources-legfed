package org.intermine.bio.dataconversion;

/*
 * Copyright (C) 2002-2019 FlyMine
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  See the LICENSE file for more
 * information or http://www.gnu.org/copyleft/lesser.html.
 *
 */

import org.intermine.bio.io.gff3.GFF3Record;
import org.intermine.metadata.Model;
import org.intermine.xml.full.Item;

/**
 * Handle special cases when converting LIS data store GFF3 files.
 *
 * @author Richard Smith
 * @author Sam Hokin
 */
public class LegfedGFF3RecordHandler extends GFF3RecordHandler {
    /**
     * Create a new LegfedGFF3RecordHandler object.
     * @param tgtModel the target Model
     */
    public LegfedGFF3RecordHandler(Model tgtModel) {
        super(tgtModel);
        // refsAndCollections controls references and collections that are set from the
        // Parent= attributes in the GFF3 file.
        refsAndCollections.put("Exon", "transcripts");
        refsAndCollections.put("Transcript", "gene");
        refsAndCollections.put("MRNA", "gene");
    }

    /**
     * {@inheritDoc}
     */
    public void process(GFF3Record record) {
        Item feature = getFeature();
        String clsName = feature.getClassName();

        // type             id                                                         
        // ---------------- -----------------------------------------------------------
        // five_prime_UTR   phavu.G19833.gnm2.ann1.Phvul.003G111100.1.five_prime_UTR.3 
        // CDS              phavu.G19833.gnm2.ann1.Phvul.003G111100.1.CDS.1            
        // three_prime_UTR  phavu.G19833.gnm2.ann1.Phvul.003G111100.1.three_prime_UTR.1
        // gene             phavu.G19833.gnm2.ann1.Phvul.003G111200                    
        // mRNA             phavu.G19833.gnm2.ann1.Phvul.003G111200.1                  

        String type = record.getType();
        String id = record.getId();

        // only update feature if ID is present
        if (id!=null) {
            // 0     1      2    3    4     5          6 7   8
            // phavu.G19833.gnm2.ann1.Phvul.003G111100.1.CDS.1
            // 0     1            2    3      4                         5
            // medtr.jemalong_A17.gnm5.ann1_6.exon:MtrunA17Chr1g0187771.1
            // 0     1            2    3      4
            // medtr.jemalong_A17.gnm5.ann1_6.gene:MtrunA17CPg0492171
            String[] parts = id.split("\\.");
            if (parts.length<5) {
                throw new RuntimeException("ID has too few dot-separated parts:"+id);
            }
            String gensp = parts[0];
            String strainName = parts[1];
            String assemblyVersion = parts[2];
            String annotationVersion = parts[3];
            String nameFromId = gensp+"."+parts[4];
            if (parts.length>5) nameFromId += "."+parts[5];
            if (parts.length>6) nameFromId += "."+parts[6];
            if (parts.length>7) nameFromId += "."+parts[7];
            if (parts.length>8) nameFromId += "."+parts[8];
            
            // set attributes
            feature.setAttribute("secondaryIdentifier", nameFromId);
            feature.setAttribute("assemblyVersion", assemblyVersion);
            feature.setAttribute("annotationVersion", annotationVersion);
        }
    }
}
