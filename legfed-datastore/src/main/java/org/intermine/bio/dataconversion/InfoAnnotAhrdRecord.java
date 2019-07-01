package org.intermine.bio.dataconversion;

import java.util.Map;
import java.util.HashMap;

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
 * Encapsulates a single tab-delimited LIS data store annot_info_ahrd file record.
 *
 * legfed_v1_0.L_QQS5LC-consensus  glucan endo-1,3-beta-glucosidase 14-like [Glycine max];
 *             ^^^^^^^^            IPR000490 (Glycoside hydrolase, family 17), IPR017853 (Glycoside hydrolase, superfamily);
 *                                 GO:0005975 (carbohydrate metabolic process)
 *
 * legfed_v1_0 = GeneFamily.version
 * L_QQS5LC = GeneFamily.identifier, OntologyAnnotation.subject
 * glucan endo-1,3-beta-glucosidase 14-like [Glycine max] = GeneFamily.description
 * IPR000490 = ProteinDomain.primaryIdentifier, GeneFamily.proteinDomain
 * Glycoside hydrolase, family 17 = ProteinDomain.description *** handle comma inside parentheses! ***
 * GO:0005975 = OntologyTerm.identifier, OntologyAnnotation.ontologyTerm
 * carbohydrate metabolic process = OntologyTerm.description  *** handle comma inside parentheses! ***
 *
 * @author Sam Hokin
 */
public class InfoAnnotAhrdRecord {
    
    String identifier;
    String version;
    String description;
    Map<String,String> interpro = new HashMap<>();
    Map<String,String> go = new HashMap<>();

    /**
     * Instantiate from a line from an LIS data store annot_nfo file. Do nothing if it's a header line.
     */
    public InfoAnnotAhrdRecord(String line) {
        if (!line.startsWith("#")) {
            try {
                String[] parts = line.split("\t");
                String identifierString = parts[0].trim();
                String valueString = parts[1].trim();

                String[] identifierParts = identifierString.split("\\.");
                version = identifierParts[0];
                identifier = identifierParts[1].replace("-consensus","");
                
                String[] valueParts = valueString.split("; ");
                description = valueParts[0];                // glucan endo-1,3-beta-glucosidase 14-like [Glycine max]
                if (valueParts.length>1) {
                    String interproString = valueParts[1];  // IPR000490 (Glycoside hydrolase, family 17), IPR017853 (Glycoside hydrolase, superfamily)
                    int index = 0;
                    while ((index=interproString.indexOf("IPR", index))>-1) {
                        String identifier = interproString.substring(index, index+9);
                        int paren = interproString.indexOf(')', index);
                        String description = interproString.substring(index+11, paren);
                        interpro.put(identifier,description);
                        index++;
                    }
                }
                if (valueParts.length>2) {
                    String goString = valueParts[2];        // GO:0005975 (carbohydrate metabolic process, foobar), GO:0005976 (carbohydrate metabolic process, barfoo)
                    int index = 0;
                    while ((index=goString.indexOf("GO:", index))>-1) {
                        String identifier = goString.substring(index, index+10);
                        int paren = goString.indexOf(')', index);
                        String description = goString.substring(index+12, paren);
                        go.put(identifier,description);
                        index++;
                    }
                }
            } catch (Exception e) {
                throw new RuntimeException(e.toString()+"\n" +
                                           line+"\n" +
                                           "identifier="+identifier+"\n" +
                                           "version="+version+"\n" +
                                           "description="+description+"\n" +
                                           "interpro="+interpro+"\n" +
                                           "go="+go);
            }
        }

    }
}
