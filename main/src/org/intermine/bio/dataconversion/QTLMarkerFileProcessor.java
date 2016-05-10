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

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;
import java.util.Set;

import org.apache.log4j.Logger;

import org.intermine.bio.util.OrganismData;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Attribute;
import org.intermine.xml.full.Item;
import org.intermine.xml.full.Reference;

/**
 * Store the QTL-marker data and relationships from a tab-delimited file.
 *
 * @author Sam Hokin, NCGR
 */
public class QTLMarkerFileProcessor extends LegfedFileProcessor {
	
    private static final Logger LOG = Logger.getLogger(QTLMarkerFileProcessor.class);

    // global scope for use in methods
    Item organism;

    /**
     * Create a new QTLMarkerFileProcessor
     * @param legfedFileConverter the LegfedFileConverter that is controlling this processor
     */
    public QTLMarkerFileProcessor(LegfedFileConverter legfedFileConverter) {
        super(legfedFileConverter);
    }

    /**
     * {@inheritDoc}
     * We process the QTL-marker relationships by reading in from a tab-delmited file.
     */
    @Override
    public void process(Reader reader) throws ObjectStoreException {

        // ---------------------------------------------------------
        // INITIAL DATA LOADING
        // ---------------------------------------------------------

        // get the organism taxon ID; enforce only one
        OrganismData[] organisms = getLegfedFileConverter().getOrganismsToProcess().toArray(new OrganismData[0]);
        if (organisms.length>1) {
            String error = "Multiple organisms specified in data source; GeneticMarkerGFFProcessor can only process one organism at a time.";
            LOG.error(error);
            throw new RuntimeException(error);
        }
        int taxonId = organisms[0].getTaxonId();

        // create organism Item - global so can be used in populate routines
        organism = getLegfedFileConverter().createItem("Organism");
        BioStoreHook.setSOTerm(getLegfedFileConverter(), organism, "organism", getLegfedFileConverter().getSequenceOntologyRefId());
        organism.setAttribute("taxonId", String.valueOf(taxonId));
        store(organism);
        LOG.info("Created and stored organism Item for taxonId="+taxonId+".");

        // store QTLs and markers in maps
        Map<String,Item> qtlMap = new HashMap<String,Item>();
        Map<String,Item> markerMap = new HashMap<String,Item>();

        // -------------------------------------------------------------------------------------------------
        // Run through the QTL-Markers file and add the associated markers to the given QTLs
        // NOTE1: marker ZZ is a placeholder, not a real marker
        // NOTE2: given names are used as _primary_ identifiers
        // -------------------------------------------------------------------------------------------------
        
        try {

            BufferedReader qtlMarkerReader = new BufferedReader(reader);
            String qtlMarkerLine = qtlMarkerReader.readLine(); // header
            LOG.info("Reading QTL-Marker file with header:"+qtlMarkerLine);
            
            while ((qtlMarkerLine=qtlMarkerReader.readLine())!=null) {
                QTLMarkerRecord rec = new QTLMarkerRecord(qtlMarkerLine);
                if (rec.qtlName!=null && rec.markerName!=null && !rec.markerName.equals("ZZ")) {
                    // find the QTL in the map, or add it with qtlName=primaryIdentifier
                    Item qtl = null;
                    if (qtlMap.containsKey(rec.qtlName)) {
                        qtl = qtlMap.get(rec.qtlName);
                    } else {
                        qtl = getLegfedFileConverter().createItem("QTL");
                        BioStoreHook.setSOTerm(getLegfedFileConverter(), qtl, "QTL", getLegfedFileConverter().getSequenceOntologyRefId());
                        qtl.setReference("organism", organism);
                        qtl.setAttribute("primaryIdentifier", rec.qtlName);
                        qtlMap.put(rec.qtlName, qtl);
                    }
                    // find the genetic marker in the map, or add it with markerName=primaryIdentifier
                    Item marker = null;
                    if (markerMap.containsKey(rec.markerName)) {
                        marker = markerMap.get(rec.markerName);
                    } else {
                        marker = getLegfedFileConverter().createItem("GeneticMarker");
                        BioStoreHook.setSOTerm(getLegfedFileConverter(), marker, "genetic_marker", getLegfedFileConverter().getSequenceOntologyRefId());
                        marker.setReference("organism", organism);
                        marker.setAttribute("primaryIdentifier", rec.markerName);
                        markerMap.put(rec.markerName, marker);
                    }
                    // add this genetic marker to this QTL's collection
                    LOG.info("Adding QTL="+rec.qtlName+" Genetic Marker="+rec.markerName+" to associatedGeneticMarkers.");
                    qtl.addToCollection("associatedGeneticMarkers", marker);
                }
            }
            
            qtlMarkerReader.close();
            LOG.info("Done associating genetic markers with QTLs.");

        } catch (Exception ex) {

            LOG.error(ex.getMessage());

        }

        // -------------------------------------------------------------
        // store the Items
        // -------------------------------------------------------------

        for (Item qtl : qtlMap.values()) store(qtl);
        for (Item marker : markerMap.values()) store(marker);

    }
    
}
