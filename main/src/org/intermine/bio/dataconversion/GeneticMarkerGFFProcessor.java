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
import java.io.IOException;
import java.math.BigDecimal;
import java.math.RoundingMode;
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
import org.intermine.xml.full.ReferenceList;

/**
 * Store the genetic marker and QTL data from a CMap tab-delimited file, along with a GFF file containing 
 * genetic markers and a QTL-markers file containing QTL-marker relations.
 *
 * Input GFF file is set in project.xml as <property name="geneticMarkerGffFile" location="/your/gff/file/path/here.gff"/>
 *
 * This processor currently makes no DB calls to the chado database, but we'll extend ChadoProcessor anyway.
 *
 * @author Sam Hokin, NCGR
 */
public class GeneticMarkerGFFProcessor extends ChadoProcessor {
	
    private static final Logger LOG = Logger.getLogger(GeneticMarkerGFFProcessor.class);

    // global objects used by populateGeneticMarker()
    Item organism;
    
    /**
     * Create a new GeneticMarkerGFFProcessor
     * @param chadoDBConverter the ChadoDBConverter that is controlling this processor
     */
    public GeneticMarkerGFFProcessor(ChadoDBConverter chadoDBConverter) {
        super(chadoDBConverter);
    }

    /**
     * {@inheritDoc}
     * We process the supplied GFF file to create GeneticMarker Items, along with Chromosome Items from the seqid column.
     */
    @Override
    public void process(Connection connection) throws SQLException, ObjectStoreException {

        // ---------------------------------------------------------
        // INITIAL DATA LOADING
        // ---------------------------------------------------------

        // get the GFF file, do nothing if it's null
        File gffFile = getChadoDBConverter().getGeneticMarkerGffFile();
        if (gffFile==null) {
            LOG.error("GFF file has not been set in project.xml.");
            System.exit(1);
        }
         
        // get chado organism_id from supplied taxon ID - enforce processing a single organism!
        Map<Integer,OrganismData> chadoToOrgData = getChadoDBConverter().getChadoIdToOrgDataMap();
        if (chadoToOrgData.size()>1) {
            LOG.error("Multiple chado organisms specified in data source; GeneticMarkerGFFProcessor can only process one organism at a time.");
            System.exit(1);
        }
        Integer organismId = 0;
        for (Integer key : chadoToOrgData.keySet()) {
            organismId = key.intValue();
        }
        int organism_id = organismId.intValue(); // chado organism.organism_id
        OrganismData orgData = chadoToOrgData.get(organismId);
        int taxonId = orgData.getTaxonId();

        // create organism Item - global so can be used in populate routines
        organism = getChadoDBConverter().createItem("Organism");
        BioStoreHook.setSOTerm(getChadoDBConverter(), organism, "organism", getChadoDBConverter().getSequenceOntologyRefId());
        organism.setAttribute("taxonId", String.valueOf(taxonId));
        store(organism);
        LOG.info("Created and stored organism Item for taxonId="+taxonId+".");

        // -------------------------------------------------------------------------------------------------------------------
        // Load the chromosomes and genetic markers from the GFF file
        // -------------------------------------------------------------------------------------------------------------------

        // store chromosomes in a map
        Map<String,Item> chromosomeMap = new HashMap<String,Item>();
        
        try {
            
            LOG.info("Reading GFF file:"+gffFile.getName());
            BufferedReader gffReader = new BufferedReader(new FileReader(gffFile));
            String gffLine = null;
            int count = 0;
            while ((gffLine=gffReader.readLine()) != null) {
                GFFRecord gff = new GFFRecord(gffLine);
                if (gff.seqid!=null && gff.attributeName!=null) {
                    // create and store this chromosome Item if not already in map; else get it from the map
                    String chrName = gff.seqid;
                    Item chromosome = null;
                    if (chromosomeMap.containsKey(chrName)) {
                        chromosome = chromosomeMap.get(chrName);
                    } else {
                        chromosome = getChadoDBConverter().createItem("Chromosome");
                        BioStoreHook.setSOTerm(getChadoDBConverter(), chromosome, "chromosome", getChadoDBConverter().getSequenceOntologyRefId());
                        chromosome.setAttribute("primaryIdentifier", chrName);
                        store(chromosome);
                        chromosomeMap.put(chrName, chromosome);
                    }
                    // create and store the genetic marker
                    Item geneticMarker = getChadoDBConverter().createItem("GeneticMarker");
                    storeGeneticMarker(chromosome, geneticMarker, gff);
                    count++;
                }
            }

            LOG.info("Read "+count+" GFF records from "+gffFile.getName());
            LOG.info("Created "+chromosomeMap.size()+" Chromosome items from "+gffFile.getName());

        } catch (Exception ex) {
            
            LOG.error(ex.getMessage());

        }

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

    /**
     * Populate a GeneticMarker Item from a GFF record
     */
    void storeGeneticMarker(Item chromosome, Item geneticMarker, GFFRecord gff) throws ObjectStoreException {
        BioStoreHook.setSOTerm(getChadoDBConverter(), geneticMarker, "genetic_marker", getChadoDBConverter().getSequenceOntologyRefId());
        geneticMarker.setReference("organism", organism);
        geneticMarker.setAttribute("primaryIdentifier", gff.attributeName);
        geneticMarker.setAttribute("type", gff.type);
        geneticMarker.setReference("chromosome", chromosome);
        geneticMarker.setAttribute("length", String.valueOf(gff.end-gff.start+1));
        Item location = getChadoDBConverter().createItem("Location");
        location.setAttribute("start", String.valueOf(gff.start));
        location.setAttribute("end", String.valueOf(gff.end));
        location.setReference("feature", geneticMarker);
        location.setReference("locatedOn", chromosome);
        store(location);
        geneticMarker.setReference("chromosomeLocation", location);
        store(geneticMarker);
    }

    /**
     * Round a double to the given number of places
     */
    public static double round(double value, int places) {
        if (places < 0) throw new IllegalArgumentException();
        BigDecimal bd = new BigDecimal(value);
        bd = bd.setScale(places, RoundingMode.HALF_UP);
        return bd.doubleValue();
    }

}
