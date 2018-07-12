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
import org.intermine.xml.full.Item;

/**
 * Store Linkage Groups and their Genetic Maps from tab-delimited files. An optional length column may be added.
 * Number is used simply for sorting purposes, but is required.
 * 
 * Taxon ID, variety and publication info (PMID or DOI) are placed in the header.
 *
 * TaxonID        3847
 * Variety        Williams82
 * PMID           123456
 * #LinkageGroup        Number GeneticMap      [Length(cM)]
 * GmComposite2003_D1a  1      GmComposite2003 120.89
 *
 * @author Sam Hokin, NCGR
 */
public class LinkageGroupFileConverter extends BioFileConverter {
	
    private static final Logger LOG = Logger.getLogger(LinkageGroupFileConverter.class);

    // Store genetic maps and other things in maps
    Map<String,Item> geneticMapMap = new HashMap<String,Item>();
    Map<String,Item> organismMap = new HashMap<String,Item>();
    Map<String,Item> publicationMap = new HashMap<String,Item>();

    /**
     * Create a new LinkageGroupFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public LinkageGroupFileConverter(ItemWriter writer, Model model) {
        super(writer, model);
    }

    /**
     * {@inheritDoc}
     * We process the LinkageGroup-marker relationships by reading in from a tab-delimited file.
     */
    @Override
    public void process(Reader reader) throws Exception {

        // don't process README files
        if (getCurrentFile().getName().contains("README")) return;

        LOG.info("Processing file "+getCurrentFile().getName()+"...");

        // header data
        String taxonId = null;
        String variety = null;
        Item organism = null;
        int pmid = 0;
        String doi = null;
        Item publication = null;

        BufferedReader linkageGroupReader = new BufferedReader(reader);
	String line;
        while ((line=linkageGroupReader.readLine())!=null) {

            // ignore comments
            if (line.startsWith("#")) continue;

            // set organism if we're ready
            if (organism==null && taxonId!=null && variety!=null) {
                String organismKey = taxonId+"_"+variety;
                if (organismMap.containsKey(organismKey)) {
                    organism = organismMap.get(organismKey);
                } else {
                    organism = createItem("Organism");
                    organism.setAttribute("taxonId", taxonId);
                    organism.setAttribute("variety", variety);
                    store(organism);
                    organismMap.put(organismKey, organism);
                    LOG.info("Stored organism "+organismKey);
                }
            }

            // set publication if we're ready
            if (publication==null && (pmid!=0)) {
                String pubKey = String.valueOf(pmid);
                if (publicationMap.containsKey(pubKey)) {
                    publication = publicationMap.get(pubKey);
                } else {
                    publication = createItem("Publication");
                    publication.setAttribute("pubMedId", String.valueOf(pmid));
                    store(publication);
                    publicationMap.put(pubKey, publication);
                    LOG.info("Stored publication "+pubKey);
                }
            } else if (publication==null && doi!=null) {
                String pubKey = doi;
                if (publicationMap.containsKey(pubKey)) {
                    publication = publicationMap.get(pubKey);
                } else {
                    publication = createItem("Publication");
                    publication.setAttribute("doi", doi);
                    store(publication);
                    publicationMap.put(pubKey, publication);
                    LOG.info("Stored publication "+pubKey);
                }
            }
                    
            String[] parts = line.split("\t");
            String key = parts[0].trim();
            String value = parts[1].trim();

            // header stuff
            if (key.toLowerCase().equals("taxonid")) {
                taxonId = value;
            } else if (key.toLowerCase().equals("variety")) {
                variety = value;
            } else if (key.toLowerCase().equals("pmid")) {
                pmid = Integer.parseInt(value);
            } else if (key.toLowerCase().equals("doi")) {
                doi = value;

            } else {
            
                // req'd
                String lgID = parts[0];
                String number = parts[1];
                String gmID = parts[2];
                // optional
                double length = 0.0;
                if (parts.length>3) length = Double.parseDouble(parts[3]);
                
                Item lg = createItem("LinkageGroup");
                lg.setAttribute("primaryIdentifier", lgID);
                lg.setAttribute("number", number);
                
                if (length>0.0) lg.setAttribute("length", String.valueOf(length));

                Item gm = null;
                if (geneticMapMap.containsKey(gmID)) {
                    gm = geneticMapMap.get(gmID);
                } else {
                    // create and store a new genetic map
                    gm = createItem("GeneticMap");
                    gm.setAttribute("primaryIdentifier", gmID);
                    if (organism!=null) gm.setReference("organism", organism);
                    if (publication!=null) gm.addToCollection("publications", publication);
                    store(gm);
                    geneticMapMap.put(gmID, gm);
                }

                lg.setReference("geneticMap", gm);
                store(lg);
                
                gm.addToCollection("linkageGroups", lg);

            }

        }
        
        linkageGroupReader.close();
    
    }

    /**
     * Do nothing, we're storing as we build above.
     */
    @Override
    public void close() throws Exception {
    }
    
}
