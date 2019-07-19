package org.intermine.bio.dataconversion;

import java.util.Map;
import java.util.HashMap;

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
 * Encapsulates a single tab-delimited LIS data store info_descriptors.txt file record containing gene descriptions and ontology annotations.
 *
 * 0.0   0.1       0.2  0.3  0.4    1.0                                       1.1                               1.2
 * arahy.Tifrunner.gnm1.ann1.GU6A2U RING finger protein 5-like [Glycine max]; IPR013083 (Zinc finger, RING...); GO:0005515 (protein binding), GO:0008270 (zinc ion binding)
 *
 * @author Sam Hokin
 */
public class InfoDescriptorsRecord {
    
    String identifier;
    String description;
    Map<String,String> interpro = new HashMap<>();
    Map<String,String> go = new HashMap<>();

    /**
     * Instantiate from a line from an LIS data store annot_nfo file. Do nothing if it's a header line.
     */
    public InfoDescriptorsRecord(String line) {
        if (!line.startsWith("#")) {
            String[] parts = line.split("\t");
            String valueString = parts[1].trim();
            String[] valueParts = valueString.split("; ");
            identifier = parts[0].trim();
            if (identifier==null || identifier.trim().length()==0) {
                throw new RuntimeException("Null or empty record identifier in InfoDescriptorsRecord:\n"+line);
            }
            description = valueParts[0];                // RING finger protein 5-like [Glycine max]
            if (valueParts.length>1) {
                String interproString = valueParts[1];  // IPR000490 (Glycoside hydrolase, family 17), IPR017853 (Glycoside hydrolase, superfamily)
                int index = 0;
                while ((index=interproString.indexOf("IPR", index))>-1) {
                    String id = interproString.substring(index, index+9);
                    int paren = interproString.indexOf(')', index);
                    String desc = interproString.substring(index+11, paren);
                    interpro.put(id,desc);
                    index++;
                }
            }
            if (valueParts.length>2) {
                String goString = valueParts[2];        // GO:0005975 (carbohydrate metabolic process, foobar), GO:0005976 (carbohydrate metabolic process, barfoo)
                int index = 0;
                while ((index=goString.indexOf("GO:", index))>-1) {
                    String id = goString.substring(index, index+10);
                    int paren = goString.indexOf(')', index);
                    String desc = goString.substring(index+12, paren);
                    go.put(id,desc);
                    index++;
                }
            }
        }
    }
}
