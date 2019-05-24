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
 * Handle special cases when converting LegFed GFF3 files.
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

        // type             id                                                          name
        // ---------------- ----------------------------------------------------------- ------------------------
        // five_prime_UTR   phavu.G19833.gnm2.ann1.Phvul.003G111100.1.five_prime_UTR.3  null
        // CDS              phavu.G19833.gnm2.ann1.Phvul.003G111100.1.CDS.1             null
        // three_prime_UTR  phavu.G19833.gnm2.ann1.Phvul.003G111100.1.three_prime_UTR.1 null
        // gene             phavu.G19833.gnm2.ann1.Phvul.003G111200                     phavu.Phvul.003G111200
        // mRNA             phavu.G19833.gnm2.ann1.Phvul.003G111200.1                   phavu.Phvul.003G111200.1

        String type = record.getType();
        String id = record.getId();
        String name = null;
        String assemblyVersion = null;
        String annotationVersion = null;
        String nameFromId = null;
        
        // 0     1      2    3    4     5          6 7   8
        // phavu.G19833.gnm2.ann1.Phvul.003G111100.1.CDS.1
        String[] parts = id.split("\\.");
        if (parts.length>2) assemblyVersion = parts[2];
        if (parts.length>3) annotationVersion = parts[3];
        if (parts.length>8) {
            nameFromId = parts[4]+"."+parts[5]+"."+parts[6]+"."+parts[7]+"."+parts[8];
        } else if (parts.length>6) {
            nameFromId = parts[4]+"."+parts[5]+"."+parts[6];
        } else if (parts.length>5) {
            nameFromId = parts[4]+"."+parts[5];
        }

        if (record.getAttributes().get("Name")!=null) {
            name = record.getAttributes().get("Name").iterator().next();
            if (name.charAt(5)=='.') name = name.substring(6);
        }

        // set attributes
        if (name!=null) {
            feature.setAttribute("secondaryIdentifier", name);
        } else if (nameFromId!=null) {
            feature.setAttribute("secondaryIdentifier", nameFromId);
        }
        if (assemblyVersion!=null) feature.setAttribute("assemblyVersion", assemblyVersion);
        if (annotationVersion!=null) feature.setAttribute("annotationVersion", annotationVersion);
    }
}
