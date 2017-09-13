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
 * Import organisms from tab-delimited germplasm files.
 *
 * <pre>
 * taxonId	3917
 * variety	CB27
 * alternateVarietyName	California Blackeye 27
 * description	California Blackeye 27 (or CB27) is ideally suited to the Central Valley of California....
 * patentNumber	200000183
 * url	https://techtransfer.universityofcalifornia.edu/NCD/10164.html
 * country	United States
 * </pre>
 *
 * @author Sam Hokin
 */
public class GermplasmFileConverter extends BioFileConverter {

    private static final Logger LOG = Logger.getLogger(GermplasmFileConverter.class);

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
        
        // the fields
        String taxonId = null; // not null
        String variety = null; // not null]
        String alternateVarietyName = null;
        String description = null;
        String patentNumber = null;
        String url = null;
        String country = null;
        
        BufferedReader br = new BufferedReader(reader);
        String line = null;
        while ((line=br.readLine())!=null) {
            String[] parts = line.split("\t");
            if (parts.length==2) {
                String field = parts[0].toLowerCase();
                String value = parts[1];
                if (field.equals("taxonid")) taxonId = value;
                if (field.equals("variety")) variety = value;
                if (field.equals("alternatevarietyname")) alternateVarietyName = value;
                if (field.equals("description")) description = value;
                if (field.equals("patentnumber")) patentNumber = value;
                if (field.equals("url")) url = value;
                if (field.equals("country")) country = value;
            }
        }

        // store the organism if taxonId and variety are given
        if (taxonId!=null && variety!=null) {
            Item organism = createItem("Organism");
            organism.setAttribute("taxonId", taxonId);
            organism.setAttribute("variety", variety);
            if (alternateVarietyName!=null) organism.setAttribute("alternateVarietyName", alternateVarietyName);
            if (description!=null) organism.setAttribute("comment", description); // use existing comment field from chado importer
            if (patentNumber!=null) organism.setAttribute("patentNumber", patentNumber);
            if (url!=null) organism.setAttribute("url", url);
            if (country!=null) organism.setAttribute("country", country);
            store(organism);
        }

        // wrap up this file
        br.close();
 
    }

}
