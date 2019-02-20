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

import java.io.BufferedReader;
import java.io.Reader;

import java.util.Map;
import java.util.HashMap;

import org.apache.log4j.Logger;

import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Item;

/**
 * Store Phenotypes and phenotype-strain associations as PhenotypeValues.
 * TaxonID 3477
 * Strain Phenotype Value
 * Strain Phenotype Value
 * Strain Phenotype Value
 *
 * Value can be text, numeric or boolean (+/-), stored in different fields as appropriate.
 *
 * @author Sam Hokin, NCGR
 */
public class PhenotypeFileConverter extends BioFileConverter {
	
    private static final Logger LOG = Logger.getLogger(PhenotypeFileConverter.class);

    // store strains and phenotypes in maps for reuse
    Map<String,Item> organismMap = new HashMap<>();
    Map<String,Item> strainMap = new HashMap<>();
    Map<String,Item> phenotypeMap = new HashMap<>();
    
    /**
     * Create a new PhenotypeFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public PhenotypeFileConverter(ItemWriter writer, Model model) {
        super(writer, model);
    }

    /**
     * {@inheritDoc}
     * Process the strain-phenotype relationships by reading in from a tab-delimited file.
     */
    @Override
    public void process(Reader reader) throws Exception {

        // don't process README files
        if (getCurrentFile().getName().contains("README")) return;

        LOG.info("Processing file "+getCurrentFile().getName()+"...");

        Item organism = null;

        BufferedReader br = new BufferedReader(reader);
	String line;
        while ((line=br.readLine())!=null) {

            if (line.startsWith("#") || line.trim().length()==0) {
                continue; // comment
            }

            String[] parts = line.split("\t");

            if (parts[0].toLowerCase().equals("taxonid")) {
                String taxonId = parts[1];
                if (organismMap.containsKey(taxonId)) {
                    organism = organismMap.get(taxonId);
                } else {
                    organism = createItem("Organism");
                    organism.setAttribute("taxonId", taxonId);
                    store(organism);
                    organismMap.put(taxonId, organism);
                }

            } else {
                String strainName = parts[0];
                String phenotypeName = parts[1];
                String value = parts[2];

                // retrieve or create the Strain item
                Item strain = null;
                if (strainMap.containsKey(strainName)) {
                    strain = strainMap.get(strainName);
                } else {
                    strain = createItem("Strain");
                    strain.setAttribute("primaryIdentifier", strainName);
                    strain.setReference("organism", organism);
                    strainMap.put(strainName, strain);
                }
            
                // retrieve or create the Phenotype item
                Item phenotype = null;
                if (phenotypeMap.containsKey(phenotypeName)) {
                    phenotype = phenotypeMap.get(phenotypeName);
                } else {
                    phenotype = createItem("Phenotype");
                    phenotype.setAttribute("primaryIdentifier", phenotypeName);
                    phenotypeMap.put(phenotypeName, phenotype);
                }
            
                // create the PhenotypeValue and associate with the strain and phenotype
                Item phenotypeValue = createItem("PhenotypeValue");
                if (value.contains(";")) {
                    // HACK: the GRIN data can have semi-colon separated values - we take the first for boolean and the mean for numeric
                    String[] measurements = value.split(";");
                    if (measurements[0].equals("+")) {
                        phenotypeValue.setAttribute("booleanValue", "true");
                    } else if (measurements[0].equals("-")) {
                        phenotypeValue.setAttribute("booleanValue", "false");
                    } else {
                        try {
                            double dval = 0.0;
                            for (int i=0; i<measurements.length; i++) {
                                dval += Double.parseDouble(measurements[i]);
                            }
                            if (dval!=0.0) dval = dval / measurements.length;
                            phenotypeValue.setAttribute("numericValue", String.valueOf(dval));
                        } catch (Exception e) {
                            phenotypeValue.setAttribute("textValue", value);
                        }
                    }
                } else if (value.equals("+") || value.equals("true")) {
                    phenotypeValue.setAttribute("booleanValue", "true");
                } else if (value.equals("-") || value.equals("false")) {
                    phenotypeValue.setAttribute("booleanValue", "false");
                } else {
                    try {
                        double dval = Double.parseDouble(value);
                        phenotypeValue.setAttribute("numericValue", value);
                    } catch (Exception e) {
                        phenotypeValue.setAttribute("textValue", value);
                    }
                }
                phenotypeValue.setReference("phenotype", phenotype);
                phenotypeValue.setReference("strain", strain);
                store(phenotypeValue);
            }
        }
        br.close();
    }

    /**
     * Store the phenotypes and strains
     */
    @Override
    public void close() throws ObjectStoreException {
        
        LOG.info("Storing "+phenotypeMap.size()+" Phenotype items...");
        store(phenotypeMap.values());
        
        LOG.info("Storing "+strainMap.size()+" Strain items...");
        store(strainMap.values());
        
    }
}
