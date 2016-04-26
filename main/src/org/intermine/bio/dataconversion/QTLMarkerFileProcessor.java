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
import java.io.File;
import java.io.FileReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
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
 * Note: this processor doesn't query the chado database.
 *
 * @author Sam Hokin, NCGR
 */
public class QTLMarkerFileProcessor extends ChadoProcessor {
	
    private static final Logger LOG = Logger.getLogger(QTLMarkerFileProcessor.class);

    // global scope for use in methods
    Item organism;

    /**
     * Create a new QTLMarkerFileProcessor
     * @param chadoDBConverter the ChadoDBConverter that is controlling this processor
     */
    public QTLMarkerFileProcessor(ChadoDBConverter chadoDBConverter) {
        super(chadoDBConverter);
    }

    /**
     * {@inheritDoc}
     * We process the chado database by reading the feature, featureloc, featurepos, featuremap, feature_relationship and featureprop tables.
     */
    @Override
    public void process(Connection connection) throws SQLException, ObjectStoreException {

        // ---------------------------------------------------------
        // INITIAL DATA LOADING
        // ---------------------------------------------------------

        // get the QTL-marker file, do nothing if it's null
        File qtlMarkerFile = getChadoDBConverter().getQtlMarkerFile();
        if (qtlMarkerFile==null) {
            LOG.error("QTL-marker file has not been set in project.xml.");
            System.exit(1);
        }
         
        // get chado organism_id from supplied taxon ID - enforce processing a single organism!
        Map<Integer,OrganismData> chadoToOrgData = getChadoDBConverter().getChadoIdToOrgDataMap();
        if (chadoToOrgData.size()>1) {
            System.err.println("ERROR - multiple chado organisms specified in data source; QTLMarkerFileProcessor can only process one organism at a time.");
            System.exit(1);
        }
        Integer organismId = 0;
        for (Integer key : chadoToOrgData.keySet()) {
            organismId = key.intValue();
        }
        OrganismData orgData = chadoToOrgData.get(organismId);
        int organism_id = organismId.intValue(); // chado organism.organism_id
        int taxonId = orgData.getTaxonId();

        // create organism Item - global so can be used in populate routines
        organism = getChadoDBConverter().createItem("Organism");
        BioStoreHook.setSOTerm(getChadoDBConverter(), organism, "organism", getChadoDBConverter().getSequenceOntologyRefId());
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

            BufferedReader qtlMarkerReader = new BufferedReader(new FileReader(qtlMarkerFile));
            String qtlMarkerLine = qtlMarkerReader.readLine(); // header
            LOG.info("Reading QTL-Marker file:"+qtlMarkerFile.getName()+" with header:"+qtlMarkerLine);
            
            while ((qtlMarkerLine=qtlMarkerReader.readLine())!=null) {
                QTLMarkerRecord rec = new QTLMarkerRecord(qtlMarkerLine);
                if (rec.qtlName!=null && rec.markerName!=null && !rec.markerName.equals("ZZ")) {
                    // find the QTL in the map, or add it with qtlName=primaryIdentifier
                    Item qtl = null;
                    if (qtlMap.containsKey(rec.qtlName)) {
                        qtl = qtlMap.get(rec.qtlName);
                    } else {
                        qtl = getChadoDBConverter().createItem("QTL");
                        BioStoreHook.setSOTerm(getChadoDBConverter(), qtl, "QTL", getChadoDBConverter().getSequenceOntologyRefId());
                        qtl.setReference("organism", organism);
                        qtl.setAttribute("primaryIdentifier", rec.qtlName);
                        qtlMap.put(rec.qtlName, qtl);
                    }
                    // find the genetic marker in the map, or add it with markerName=primaryIdentifier
                    Item marker = null;
                    if (markerMap.containsKey(rec.markerName)) {
                        marker = markerMap.get(rec.markerName);
                    } else {
                        marker = getChadoDBConverter().createItem("GeneticMarker");
                        BioStoreHook.setSOTerm(getChadoDBConverter(), marker, "genetic_marker", getChadoDBConverter().getSequenceOntologyRefId());
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

    /**
     * Store the item.
     * @param item the Item
     * @return the database id of the new Item
     * @throws ObjectStoreException if an error occurs while storing
     */
    protected Integer store(Item item) throws ObjectStoreException {
        return getChadoDBConverter().store(item);
    }
    
    /**
     * Do any extra processing that is needed before the converter starts querying features
     * @param connection the Connection
     * @throws ObjectStoreException if there is a object store problem
     * @throws SQLException if there is a database problem
     */
    protected void earlyExtraProcessing(Connection connection) throws ObjectStoreException, SQLException {
        // override in subclasses as necessary
    }

    /**
     * Do any extra processing for this database, after all other processing is done
     * @param connection the Connection
     * @param featureDataMap a map from chado feature_id to data for that feature
     * @throws ObjectStoreException if there is a problem while storing
     * @throws SQLException if there is a problem
     */
    protected void extraProcessing(Connection connection, Map<Integer, FeatureData> featureDataMap)
        throws ObjectStoreException, SQLException {
        // override in subclasses as necessary
    }

    /**
     * Perform any actions needed after all processing is finished.
     * @param connection the Connection
     * @param featureDataMap a map from chado feature_id to data for that feature
     * @throws SQLException if there is a problem
     */
    protected void finishedProcessing(Connection connection, Map<Integer, FeatureData> featureDataMap) throws SQLException {
        // override in subclasses as necessary
    }

}
