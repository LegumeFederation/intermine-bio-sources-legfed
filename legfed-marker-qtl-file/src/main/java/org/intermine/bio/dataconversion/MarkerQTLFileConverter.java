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
 * Store genetic marker / QTL associations
 *
 * Marker QTL
 *
 * @author Sam Hokin, NCGR
 */
public class MarkerQTLFileConverter extends BioFileConverter {
	
    private static final Logger LOG = Logger.getLogger(MarkerQTLFileConverter.class);

    // store markers and linkage groups in maps for repeated use
    Map<String,Item> markerMap = new HashMap<String,Item>();
    Map<String,Item> qtlMap = new HashMap<String,Item>();
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
        String taxonId = null;
        String variety = null;
        Item organism = null;

        BufferedReader markerReader = new BufferedReader(reader);
	String line;
        while ((line=markerReader.readLine())!=null) {

            // initialize organism if not set and can be
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
                    LOG.info("Storing organism: "+taxonId+" ("+variety+")");
                }
            }

            if (line.startsWith("#") || line.trim().length()==0) {

                // do nothing, comment
                
            } else if (line.startsWith("TaxonID")) {
                
                String[] parts = line.split("\t");
                taxonId = parts[1];
                
            } else if (line.startsWith("Variety")) {
                
                String[] parts = line.split("\t");
                variety = parts[1];
                
            } else {

                // bail if organism not set, otherwise merging is nightmare
                if (organism==null) {
                    LOG.error("Organism not set: taxonId="+taxonId+", variety="+variety);
                    throw new RuntimeException("Organism not set: taxonId="+taxonId+", variety="+variety);
                }
                        
                String[] parts = line.split("\t");
                String markerID = parts[0];
                String qtlID = parts[1];
                
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
                
                // retrieve or create the QTL item
                Item qtl = null;
                if (qtlMap.containsKey(qtlID)) {
                    qtl = qtlMap.get(qtlID);
                } else {
                    qtl = createItem("QTL");
                    qtl.setAttribute("primaryIdentifier", qtlID);
                    // qtl.setReference("organism", organism);
                    qtlMap.put(qtlID, qtl);
                    LOG.info("Storing QTL "+qtlID);
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
    public void close() throws Exception {

        LOG.info("Storing "+markerMap.size()+" GeneticMarker items...");
        store(markerMap.values());

        LOG.info("Storing "+qtlMap.size()+" QTL items...");
        store(qtlMap.values());
        
    }
    
}
