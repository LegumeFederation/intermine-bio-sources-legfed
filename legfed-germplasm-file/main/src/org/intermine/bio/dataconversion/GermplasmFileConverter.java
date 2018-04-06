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

import java.util.List;
import java.util.ArrayList;
import java.util.Map;
import java.util.HashMap;

import org.apache.log4j.Logger;

import org.ncgr.intermine.PubMedPublication;

import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Item;

/**
 * Import organisms from tab-delimited germplasm files.
 *
 * <pre>
 * TaxonID	3917
 * Variety	CB27
 * AlternateVarietyName	California Blackeye 27
 * Description	California Blackeye 27 (or CB27) is ideally suited to the Central Valley of California....
 * PatentNumber	200000183
 * URL          https://techtransfer.universityofcalifornia.edu/NCD/10164.html
 * Country	United States
 * PMID         123456
 * PMID         234567
 * </pre>
 *
 * @author Sam Hokin
 */
public class GermplasmFileConverter extends BioFileConverter {

    private static final Logger LOG = Logger.getLogger(GermplasmFileConverter.class);

    // store all the pubs and authors in maps so we don't duplicate
    Map<String,Item> publicationMap = new HashMap<String,Item>();
    Map<String,Item> authorMap = new HashMap<String,Item>();
    
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
        if (getCurrentFile().getName().contains("README")) return;

        LOG.info("Processing Germplasm file "+getCurrentFile().getName()+"...");

        // store the pubs for THIS organism in a list
        List<Item> publicationList = new ArrayList<Item>();
        
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
                if (field.equals("pmid")) {
                    // create or grab this publication and add it to the publicationList
                    String pubMedId = value;
                    if (publicationMap.containsKey(pubMedId)) {
                        Item publication = publicationMap.get(pubMedId);
                        publicationList.add(publication);
                    } else {
                        Item publication = createItem("Publication");
                        publication.setAttribute("pubMedId", pubMedId);
                        store(publication);
                        publicationMap.put(pubMedId, publication);
                        publicationList.add(publication);
                    }
                }
            }
        }

        // store the organism if taxonId and variety are given
        if (taxonId!=null && variety!=null) {
            Item organism = createItem("Organism");
            organism.setAttribute("taxonId", taxonId);
            organism.setAttribute("variety", variety);
            if (alternateVarietyName!=null) organism.setAttribute("alternateVarietyName", alternateVarietyName);
            if (description!=null) organism.setAttribute("description", description);
            if (patentNumber!=null) organism.setAttribute("patentNumber", patentNumber);
            if (url!=null) organism.setAttribute("url", url);
            if (country!=null) organism.setAttribute("country", country);
            for (Item publication : publicationList) {
                organism.addToCollection("publications", publication);
            }
            store(organism);
        }

        // wrap up this file
        br.close();
 
    }

}
