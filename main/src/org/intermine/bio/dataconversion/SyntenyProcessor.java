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
import org.intermine.xml.full.ReferenceList;

/**
 * Read syntenic regions for two organisms from a GFF file and store them as SyntenicRegion items.
 * This is designed to use the standard (?) GFF annotation produced by DAGchainer.
 *
 * NOTE: input text files are currently HARDCODED
 *
 * @author Sam Hokin, NCGR
 */
public class SyntenyProcessor extends ChadoProcessor {
	
    private static final Logger LOG = Logger.getLogger(SyntenyProcessor.class);

    // NOTE: HARDCODED inter-organism synteny GFF FILE
    String gffFilename = "/home/intermine/data/synteny_Pv_with_Gm.gff3";

    /**
     * Create a new SyntenyProcessor
     * @param chadoDBConverter the ChadoDBConverter that is controlling this processor
     */
    public SyntenyProcessor(ChadoDBConverter chadoDBConverter) {
        super(chadoDBConverter);
    }

    /**
     * {@inheritDoc}
     * We process the GFF file by creating SytenicRegion items and storing them.
     * In order to do so, we need to access the two organisms' chromosomes from the chado database.
     */
    @Override
    public void process(Connection connection) throws SQLException, ObjectStoreException {

        // initialize our DB statement and stuff for chado queries
        Statement stmt = connection.createStatement();
        String query;
        ResultSet rs;
        
        // ---------------------------------------------------------
        // ---------------- INITIAL DATA LOADING -------------------
        // ---------------------------------------------------------
        // NOTE: HARDCODED organisms; don't know how to send in the source and target organisms from project.xml!

        // chado
        int source_organism_id = 24; // Pv
        int target_organism_id = 16; // Gm

        // InterMine
        int sourceTaxonId = 3885; // Pv
        int targetTaxonId = 3847; // Gm

        // create organism Items
        Item sourceOrganism = getChadoDBConverter().createItem("Organism");
        sourceOrganism.setAttribute("taxonId", String.valueOf(sourceTaxonId));
        Item targetOrganism = getChadoDBConverter().createItem("Organism");
        targetOrganism.setAttribute("taxonId", String.valueOf(targetTaxonId));

        // may as well store them now
        store(sourceOrganism);
        store(targetOrganism);
        
        LOG.info("Created and stored sourceOrganism, taxonId="+sourceTaxonId+" and targetOrganism, taxonId="+targetTaxonId);

        // -------------------------------------------------------------------------------------------------------------------
        // Load the GFF data into a map. Create source and target chromosome Items.
        //
        // NOTE: this is hardcoded for Pv source and Gm target chromosomes.
        // -------------------------------------------------------------------------------------------------------------------

        Map<String,GFFRecord> gffMap = new HashMap<String,GFFRecord>();
        ItemMap sourceChromosomeMap = new ItemMap();
        ItemMap targetChromosomeMap = new ItemMap();
        
        try {
            
            BufferedReader gffReader = new BufferedReader(new FileReader(gffFilename));
            LOG.info("Reading GFF file:"+gffFilename);
            String gffLine = null;
            while ((gffLine=gffReader.readLine()) != null) {
                GFFRecord gff = new GFFRecord(gffLine);
                if (gff.seqid!=null && gff.attributeName!=null) {
                    // create source chromosome, store and add to sourceChromosomeMap if needed
                    String sourceName = gff.seqid.replace("Pv","phavu.Chr"); // chado chromosome name HARDCODED
                    if (!sourceChromosomeMap.containsSecondaryIdentifier(sourceName)) {
                        // create the chromosome Item from the chado record
                        query = "SELECT * FROM feature WHERE type_id=43403 AND organism_id="+source_organism_id+" AND name='"+sourceName+"'";
                        rs = stmt.executeQuery(query);
                        if (rs.next()) {
                            ChadoFeature cf = new ChadoFeature(rs);
                            Item chromosome = getChadoDBConverter().createItem("Chromosome");
                            cf.populateBioEntity(chromosome, sourceOrganism);
                            store(chromosome);
                            sourceChromosomeMap.put(chromosome);
                            LOG.info("Stored new source chromosome "+chromosome.getAttribute("primaryIdentifier").getValue());
                        }
                        rs.close();
                    }
                    // create target chromosome, store and add to targetChromosomeMap if needed
                    String targetChromosome = gff.getTargetChromosome();
                    if (targetChromosome!=null) {
                        String targetName = targetChromosome.replace("Gm","glyma.Chr"); // chado chromosome name HARDCODED
                        if (!targetChromosomeMap.containsSecondaryIdentifier(targetName)) {
                            // create the chromosome Item from the chado record
                            query = "SELECT * FROM feature WHERE type_id=43403 AND organism_id="+target_organism_id+" AND name='"+targetName+"'";
                            rs = stmt.executeQuery(query);
                            if (rs.next()) {
                                ChadoFeature cf = new ChadoFeature(rs);
                                Item chromosome = getChadoDBConverter().createItem("Chromosome");
                                cf.populateBioEntity(chromosome, targetOrganism);
                                store(chromosome);
                                targetChromosomeMap.put(chromosome);
                                LOG.info("Stored new target chromosome "+chromosome.getAttribute("primaryIdentifier").getValue());
                            }
                            rs.close();
                        }
                    }
                    // store in map if it's a syntenic_region
                    if (gff.type.equals("syntenic_region")) {
                        gffMap.put(gff.attributeName, gff);
                        LOG.info("Syntenic region: "+gff.attributeName);
                    }
                }
            }
            gffReader.close();
            
            LOG.info("Read "+gffMap.size()+" GFF records from "+gffFilename+".");
            LOG.info("Created "+sourceChromosomeMap.size()+" source Chromosome items from "+gffFilename+".");
            LOG.info("Created "+targetChromosomeMap.size()+" target Chromosome items from "+gffFilename+".");

        } catch (Exception ex) {

            LOG.error(ex.toString());
            ex.printStackTrace();
            System.exit(1);

        }

        // ----------------------------------------------------------------------
        // Now spin through the gffMap records and store the SyntenicRegion Items
        // ----------------------------------------------------------------------

        LOG.info("Creating and storing synteny regions...");

        for (GFFRecord gff : gffMap.values())  {

            // populate and store the target region and its location
            Item targetChromosome = targetChromosomeMap.getBySecondaryIdentifier(gff.getTargetChromosome().replace("Gm","glyma.Chr")); // HARDCODED
            Item targetRegion = getChadoDBConverter().createItem("SequenceFeature");
            Item targetChromosomeLocation = getChadoDBConverter().createItem("Location");
            gff.populateDAGchainerRegion(targetRegion, targetOrganism, targetChromosome, targetChromosomeLocation);
            store(targetRegion);
            store(targetChromosomeLocation);
            
            // populate and store the source region and its location
            Item sourceChromosome = sourceChromosomeMap.getBySecondaryIdentifier(gff.seqid.replace("Pv","phavu.Chr")); // HARDCODED
            Item sourceRegion = getChadoDBConverter().createItem("SyntenicRegion");
            Item sourceChromosomeLocation = getChadoDBConverter().createItem("Location");
            gff.populateSequenceFeature(sourceRegion, sourceOrganism, sourceChromosome, sourceChromosomeLocation);
            sourceRegion.setAttribute("medianKs", String.valueOf(gff.attributeMedianKs));
            sourceRegion.setReference("targetRegion", targetRegion);
            store(sourceRegion);
            store(sourceChromosomeLocation);

        }

        LOG.info("...done.");    

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
