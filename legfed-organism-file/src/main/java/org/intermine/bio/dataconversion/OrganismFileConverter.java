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

import java.util.List;
import java.util.ArrayList;

import org.apache.log4j.Logger;

import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.xml.full.Item;

/**
 * Loads organism and strain info from a tab-delimited file containing a single organism and its strains:
 *
 * organism.taxonId	385
 * organism.genus	Phaseolus
 * organism.species	vulgaris
 * organism.name	Phaseolus vulgaris
 * organism.shortname	P. vulgaris
 * organism.commonName	common bean
 * organism.description	Common bean's center of origin is Central America....
 * ##
 * strain.1.identifier	G19839	
 * strain.1.name	G 19839
 * strain.1.origin	Andes region, South America
 * strain.1.description	G 19839 is a climbing bean of Andean origin, and is the strain used for the Phaseolus vulgaris reference genome.
 *
 * @author Sam Hokin
 */
public class OrganismFileConverter extends BioFileConverter {
	
    private static final Logger LOG = Logger.getLogger(OrganismFileConverter.class);

    /**
     * Create a new OrganismFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public OrganismFileConverter(ItemWriter writer, Model model) {
        super(writer, model);
    }

    /**
     * {@inheritDoc}
     * Process the strains by reading in from a tab-delimited file.
     */
    @Override
    public void process(Reader reader) throws Exception {
        // don't process README files
        if (getCurrentFile().getName().contains("README")) return;

        LOG.info("Processing file "+getCurrentFile().getName()+"...");

        List<Item> strains = new ArrayList<>();
        
	Item organism = null;
        Item strain = null;
        int strainNum = 0; // keeps track of the current strain

        BufferedReader br = new BufferedReader(reader);
	String line;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) {
		continue; // comment
	    }
	    String[] parts = line.split("\t");
            String[] record = parts[0].split("\\.");
            // DEBUG
            System.out.println("record[0]="+record[0]);
            if (record[0].equals("organism")) {
                if (organism==null) organism = createItem("Organism");
                String attributeName = record[1];
                String attributeValue = parts[1];
                organism.setAttribute(attributeName, attributeValue);
            } else if (record[0].equals("strain")) {
                int num = Integer.parseInt(record[1]);
                if (num!=strainNum) {
                    strain = createItem("Strain");
                    strain.setReference("organism", organism);
                    strains.add(strain);
                    strainNum = num;
                }
                String attributeName = record[2];
                String attributeValue = parts[1];
                strain.setAttribute(attributeName, attributeValue);
            }
	}
        br.close();
        store(organism);
        store(strains);
    }
}
