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

import java.util.List;
import java.util.Map;
import java.util.HashMap;

import org.intermine.bio.io.gff3.GFF3Record;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreWriter;
import org.intermine.objectstore.ObjectStoreWriterFactory;
import org.intermine.xml.full.Item;
import org.intermine.xml.full.ItemFactory;
import org.intermine.xml.full.ItemHelper;

/**
 * Handle special cases when converting LIS data store GFF3 files.
 *
 * gene            ID=phavu.G19833.gnm1.ann1.Phvul.001G000100;
 * mRNA            ID=phavu.G19833.gnm1.ann1.Phvul.001G000100.1;Parent=phavu.G19833.gnm1.ann1.Phvul.001G000100;
 * CDS             ID=phavu.G19833.gnm1.ann1.Phvul.001G000100.1.v1.0.CDS.1;Parent=phavu.G19833.gnm1.ann1.Phvul.001G000100.1;
 * five_prime_UTR  ID=phavu.G19833.gnm1.ann1.Phvul.001G000100.1.v1.0.five_prime_UTR.1;Parent=phavu.G19833.gnm1.ann1.Phvul.001G000100.1;
 * three_prime_UTR ID=phavu.G19833.gnm1.ann1.Phvul.001G000100.1.v1.0.three_prime_UTR.1;Parent=phavu.G19833.gnm1.ann1.Phvul.001G000100.1;
 *
 * @author Richard Smith
 * @author Sam Hokin
 */
public class LegfedGFF3RecordHandler extends GFF3RecordHandler {

    ItemFactory itemFactory;
    Map<String, String> aliases = new HashMap<String, String>();
    Map<String, Integer> ids = new HashMap<String, Integer>();
    int nextClsId = 0;

    /**
     * Create a new LegfedGFF3RecordHandler object.
     * @param tgtModel the target Model
     */
    public LegfedGFF3RecordHandler(Model tgtModel) {
        super(tgtModel);
        // refsAndCollections controls references and collections that are set from the
        // Parent= attributes in the GFF3 file.
        refsAndCollections.put("MRNA", "gene");
        refsAndCollections.put("CDS", "mRNA");
        refsAndCollections.put("FivePrimeUTR", "mRNA");
        refsAndCollections.put("ThreePrimeUTR", "mRNA");
        refsAndCollections.put("Exon", "mRNA");
    }

    /**
     * {@inheritDoc}
     * type             id                                                         
     * ---------------- -----------------------------------------------------------
     * gene             phavu.G19833.gnm2.ann1.Phvul.003G111200                    
     * mRNA             phavu.G19833.gnm2.ann1.Phvul.003G111200.1
     * CDS              phavu.G19833.gnm2.ann1.Phvul.003G111100.1.CDS.1            
     * five_prime_UTR   phavu.G19833.gnm2.ann1.Phvul.003G111100.1.five_prime_UTR.3 
     * three_prime_UTR  phavu.G19833.gnm2.ann1.Phvul.003G111100.1.three_prime_UTR.1
     */
    public void process(GFF3Record record) {
        Item feature = getFeature();
        String className = feature.getClassName();
        String type = record.getType();
        String id = record.getId();
        Map<String,List<String>> attributesMap = record.getAttributes();

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
            String strainId = parts[1];
            String assemblyVersion = parts[2];
            String annotationVersion = parts[3];

            // set attributes
            feature.setAttribute("assemblyVersion", assemblyVersion);
            feature.setAttribute("annotationVersion", annotationVersion);
            if (attributesMap.containsKey("Name")) {
                String name = attributesMap.get("Name").get(0);
                feature.setAttribute("secondaryIdentifier", name);
            } else {
                String nameFromId = gensp+"."+parts[4];
                if (parts.length>5) nameFromId += "."+parts[5];
                if (parts.length>6) nameFromId += "."+parts[6];
                if (parts.length>7) nameFromId += "."+parts[7];
                if (parts.length>8) nameFromId += "."+parts[8];
                feature.setAttribute("secondaryIdentifier", nameFromId);
            }

            // add marker type = SNP if it is a marker with length 1
            if (type.equals("genetic_marker") && (record.getStart()-record.getEnd())==0) {
                feature.setAttribute("type", "SNP");
            }

            // handle other attributes
            for (String key : attributesMap.keySet()) {
                List<String> attributes = attributesMap.get(key);
                if (key.equals("Note")) {
                    // [50S ribosomal protein L18, chloroplastic-like [Glycine max]; IPR005484 (Ribosomal protein L18/L5); ...]
                    feature.setAttribute("description", attributes.get(0));
                } else if (type.equals("gene") && key.equals("Dbxref")) {
                    // [Gene3D:G3DSA:3.30.420.100, InterPro:IPR005484, PANTHER:PTHR12899, Pfam:PF00861, Superfamily:SSF53137]
                    for (String term : attributes) {
                        String[] pieces = term.split(":");
                        if (pieces[0].equals("Gene3D")) {
                            // need ProteinDomain Item
                        } else if (pieces[0].equals("InterPro")) {
                            String identifier = pieces[1];
			    // need ProteinDomain Item
                        } else if (pieces[0].equals("PANTHER")) {
                            String identifier = pieces[1];
                            // need ProteinDomain Item
                        } else if (pieces[0].equals("Pfam")) {
                            String identifier = pieces[1];
                            // need ProteinDomain Item
                        } else if (pieces[0].equals("Superfamily")) {
                            // need ProteinDomain Item
                        }
                    }
                } else if (type.equals("gene") && key.equals("Ontology_term")) {
                    // [GO:0003735, GO:0005622, GO:0005840, GO:0006412]
                    for (String term : attributes) {
                        if (term.startsWith("GO:")) {
                            // need OntologyTerm
                        }
                    }
                } else if (key.equals("evid_id")) {
                    // [GAR_10012494]
                    // do nothing
                }
            }
        }
    }
}
