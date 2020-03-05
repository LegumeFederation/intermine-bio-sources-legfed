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
 * Store GeneticMarker / QTL / Phenotyp associations
 *
 * TaxonID 3477
 * Marker QTL Phenotype
 *
 * @author Sam Hokin
 */
public class MarkerQTLFileConverter extends BioFileConverter {
	
    private static final Logger LOG = Logger.getLogger(MarkerQTLFileConverter.class);

    // store markers and linkage groups in maps for repeated use
    Map<String,Item> markerMap = new HashMap<String,Item>();
    Map<String,Item> qtlMap = new HashMap<String,Item>();
    Map<String,Item> phenotypeMap = new HashMap<String,Item>();
    Map<String,Item> organismMap = new HashMap<String,Item>();
    
    /**
     * Create a new MarkerQTLFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public MarkerQTLFileConverter(ItemWriter writer, Model model) {
        super(writer, model);
    }

    /**
     * {@inheritDoc}
     * Process the marker-linkage group relationships by reading in from a tab-delimited file.
     */
    @Override
    public void process(Reader reader) throws Exception {

        // don't process README files
        if (getCurrentFile().getName().contains("README")) return;

        LOG.info("Processing file "+getCurrentFile().getName()+"...");

        // header constants
        Item organism = null;

        BufferedReader markerReader = new BufferedReader(reader);
	String line;
        while ((line=markerReader.readLine())!=null) {

            String[] parts = line.split("\t");

            if (line.startsWith("#") || line.trim().length()==0) {

                // do nothing, comment
                
            } else if (parts[0].toLowerCase().equals("taxonid")) {
                
                String taxonId = parts[1];
                if (organismMap.containsKey(taxonId)) {
                    organism = organismMap.get(taxonId);
                } else {
                    organism = createItem("Organism");
                    organism.setAttribute("taxonId", taxonId);
                    store(organism);
                    organismMap.put(taxonId, organism);
                    LOG.info("Stored organism: "+taxonId);
                }
                
            } else {

                // bail if organism not set, otherwise merging is nightmare
                if (organism==null) {
                    LOG.error("Organism not set: supply TaxonID in header.");
                    throw new RuntimeException("Organism not set: supply TaxonID in header.");
                }
                        
                String markerID = parts[0];
                String qtlID = parts[1];
		// phenotype is optional
		String phenotypeID = null;
		if (parts.length>2) {
		    phenotypeID = parts[2];
		}
                
                // retrieve or create the marker item
                Item marker = null;
                if (markerMap.containsKey(markerID)) {
                    marker = markerMap.get(markerID);
                } else {
                    marker = createItem("GeneticMarker");
                    marker.setAttribute("primaryIdentifier", markerID);
                    marker.setReference("organism", organism);
                    markerMap.put(markerID, marker);
                    LOG.info("Storing marker "+markerID);
                }

		// retrieve or create the Phenotype item, if provided
		Item phenotype = null;
		if (phenotypeID!=null) {
		    if (phenotypeMap.containsKey(phenotypeID)) {
			phenotype = phenotypeMap.get(phenotypeID);
		    } else {
			phenotype = createItem("Phenotype");
			phenotype.setAttribute("name", phenotypeID);
			store(phenotype);
			phenotypeMap.put(phenotypeID, phenotype);
			LOG.info("Stored phenotype "+phenotypeID);
		    }
		}
                
                // retrieve or create the QTL item
                Item qtl = null;
                if (qtlMap.containsKey(qtlID)) {
                    qtl = qtlMap.get(qtlID);
                } else {
                    qtl = createItem("QTL");
                    qtl.setAttribute("primaryIdentifier", qtlID);
                    qtl.setReference("organism", organism);
		    if (phenotype!=null) qtl.setReference("phenotype", phenotype);
                    qtlMap.put(qtlID, qtl);
                    LOG.info("Added QTL "+qtlID);
                }
                
                // relate the two
                marker.addToCollection("QTLs", qtl);
                qtl.addToCollection("markers", marker);
            }
            
        }
        
        markerReader.close();
    }

    /**
     * Store the markers and linkage groups
     */
    @Override
    public void close() throws ObjectStoreException {
        LOG.info("Storing "+markerMap.size()+" GeneticMarker items...");
        store(markerMap.values());
        LOG.info("Storing "+qtlMap.size()+" QTL items...");
        store(qtlMap.values());
    }
}
