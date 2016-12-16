package org.intermine.bio.dataconversion;

/*
 * Copyright (C) 2016 NCGR
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  See the LICENSE file for more
 * information or http://www.gnu.org/copyleft/lesser.html.
 *
 */

import java.io.BufferedReader;
import java.io.Reader;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Item;

/**
 * Store information about germplasms from tab-delimited germplasm files.
 * <pre>
 * taxonId	3917
 * primaryIdentifier	CB27
 * secondaryIdentifier	California Blackeye 27
 * description	"California Blackeye 27" (or CB27) is ideally suited to the Central Valley of California....
 * patentNumber	200000183
 * url	https://techtransfer.universityofcalifornia.edu/NCD/10164.html
 * country	United States
 * </pre>
 *
 * @author Sam Hokin
 */
public class GermplasmFileConverter extends BioFileConverter {

    private static final Logger LOG = Logger.getLogger(GermplasmFileConverter.class);

    // store organisms so we don't duplicate
    Map<String,Item> organismMap = new HashMap<String,Item>();

    /**
     * Create a new GermplasmFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public GermplasmFileConverter(ItemWriter writer, Model model) {
        super(writer, model);
    }

    /**
     * {@inheritDoc}
     * Read in each genotype file and parse the germplasms, organism, lines, markers, and genotypes
     */
    @Override
    public void process(Reader reader) throws Exception {

        // don't process README files
        if (getCurrentFile().getName().equals("README")) return;

        LOG.info("Processing Germplasm file "+getCurrentFile().getName()+"...");
        
        // ----------------------------------------------------------
        // lines, genetic markers, genotypes from the Germplasm file
        // ----------------------------------------------------------

        Item germplasm = createItem("Germplasm");
        BufferedReader br = new BufferedReader(reader);
        String line = null;
        while ((line=br.readLine())!=null) {

            String[] parts = line.split("\t");
            if (parts.length==2) {
                
                String field = parts[0].toLowerCase();
                String value = parts[1];
                
                if (field.equals("taxonid")) {
                    Item organism;
                    if (organismMap.containsKey(value)) {
                        organism = organismMap.get(value);
                    } else {
                        organism = createItem("Organism");
                        BioStoreHook.setSOTerm(this, organism, "organism", getSequenceOntologyRefId());
                        organism.setAttribute("taxonId", value);
                        store(organism);
                        organismMap.put(value, organism);
                    }
                    germplasm.setReference("organism", organism);
                }

                if (field.equals("primaryidentifier")) germplasm.setAttribute("primaryIdentifier", value);
                if (field.equals("secondaryidentifier")) germplasm.setAttribute("secondaryIdentifier", value);
                if (field.equals("description")) germplasm.setAttribute("description", value);
                if (field.equals("patentnumber")) germplasm.setAttribute("patentNumber", value);
                if (field.equals("url")) germplasm.setAttribute("url", value);
                if (field.equals("country")) germplasm.setAttribute("country", value);

            }

        }

        // store the germplasm from this file
        store(germplasm);

        // wrap up this file
        br.close();
 
    }

}
