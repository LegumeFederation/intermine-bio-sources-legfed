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

import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.HashMap;
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
 * Store the various genetic data from the LegFed/LIS chado database.
 * These come from featuremap, featurepos, featureloc, feature_relationship and, of course, feature.
 *
 * Since this processer deals only with chado data, Items are stored in maps with Integer keys equal to
 * the chado feature.feature_id.
 *
 * @author Sam Hokin, NCGR
 */
public class GeneticProcessor extends ChadoProcessor {
	
    private static final Logger LOG = Logger.getLogger(GeneticProcessor.class);

    /**
     * Create a new GeneticProcessor
     * @param chadoDBConverter the ChadoDBConverter that is controlling this processor
     */
    public GeneticProcessor(ChadoDBConverter chadoDBConverter) {
        super(chadoDBConverter);
    }

    /**
     * {@inheritDoc}
     * We process the chado database by reading the feature, featureloc, featurepos, featuremap, feature_relationship and featureprop tables.
     */
    @Override
    public void process(Connection connection) throws SQLException, ObjectStoreException {
        
        // initialize our DB statement and other stuff
        Statement stmt = connection.createStatement();
        String query;
        ResultSet rs;
        
        // ---------------------------------------------------------
        // ---------------- INITIAL DATA LOADING -------------------
        // ---------------------------------------------------------

        // get chado organism_id from supplied taxon ID - enforce processing a single organism!
        Map<Integer,OrganismData> chadoToOrgData = getChadoDBConverter().getChadoIdToOrgDataMap();
        if (chadoToOrgData.size()>1) {
            System.err.println("ERROR - multiple chado organisms specified in data source; GeneticProcessor can only process one organism at a time.");
            System.exit(1);
        }
        Integer organismId = 0;
        for (Integer key : chadoToOrgData.keySet()) {
            organismId = key.intValue();
        }
        OrganismData org = chadoToOrgData.get(organismId);
        int organism_id = organismId.intValue();
        int taxonId = org.getTaxonId();

        // create organism Item
        Item organism = getChadoDBConverter().createItem("Organism");
        organism.setAttribute("taxonId", String.valueOf(taxonId));
        store(organism);
        LOG.info("Created and stored organism Item for taxonId="+taxonId+" and chado organism_id="+organism_id);

        
        // ---------------------------------------------------------
        // get the cvterm_id values for all our items of interest, to simplify queries below
        // ---------------------------------------------------------

        // genes
        int geneCVTermId = 0;
        rs = stmt.executeQuery("SELECT * FROM cvterm WHERE name='gene'");
        if (rs.next()) geneCVTermId = rs.getInt("cvterm_id");
        rs.close();

        // linkage groups
        int linkageGroupCVTermId = 0;
        rs = stmt.executeQuery("SELECT * FROM cvterm WHERE name='linkage_group'");
        if (rs.next()) linkageGroupCVTermId = rs.getInt("cvterm_id");
        rs.close();

        // genetic markers
        int geneticMarkerCVTermId = 0;
        rs = stmt.executeQuery("SELECT * FROM cvterm WHERE name='genetic_marker'");
        if (rs.next()) geneticMarkerCVTermId = rs.getInt("cvterm_id");
        rs.close();

        // QTLs
        int qtlCVTermId = 0;
        rs = stmt.executeQuery("SELECT * FROM cvterm WHERE name='QTL'");
        if (rs.next()) qtlCVTermId = rs.getInt("cvterm_id");
        rs.close();

        // gene families
        int geneFamilyCVTermId = 0;
        rs = stmt.executeQuery("SELECT * FROM cvterm WHERE name='gene family'");
        if (rs.next()) geneFamilyCVTermId = rs.getInt("cvterm_id");
        rs.close();

        // -----------------------------------------------------------------------------
        // load the Items into maps
        // -----------------------------------------------------------------------------
        
        // GENES
        // load the relevant genes from the feature table
        Map<Integer,Item> geneMap = new HashMap<Integer,Item>();
        query = "SELECT * FROM feature WHERE organism_id="+organism_id+" AND type_id="+geneCVTermId;
        LOG.info("executing query: "+query);
        rs = stmt.executeQuery(query);
        while (rs.next()) {
            ChadoFeature cf = new ChadoFeature(rs);
            Item gene = getChadoDBConverter().createItem("Gene");
            cf.populateBioEntity(gene, organism);
            geneMap.put(new Integer(cf.feature_id), gene);
        }
        rs.close();
        LOG.info("Created "+geneMap.size()+" Gene items.");

        // LINKAGE GROUPS
        // load the linkage groups from the feature table
        Map<Integer,Item> linkageGroupMap = new HashMap<Integer,Item>();
        query = "SELECT * FROM feature WHERE organism_id="+organism_id+" AND type_id="+linkageGroupCVTermId;
        LOG.info("executing query: "+query);
        rs = stmt.executeQuery(query);
        while (rs.next()) {
            ChadoFeature cf = new ChadoFeature(rs);
            Item linkageGroup = getChadoDBConverter().createItem("LinkageGroup");
            cf.populateBioEntity(linkageGroup, organism);
            linkageGroupMap.put(new Integer(cf.feature_id), linkageGroup);
        }
        rs.close();
        LOG.info("Created "+linkageGroupMap.size()+" LinkageGroup items.");

        // GENETIC MARKERS
        // load the genetic markers from the feature table
        Map<Integer,Item> geneticMarkerMap = new HashMap<Integer,Item>();
        query = "SELECT * FROM feature WHERE organism_id="+organism_id+" AND type_id="+geneticMarkerCVTermId;
        LOG.info("executing query: "+query);
        rs = stmt.executeQuery(query);
        while (rs.next()) {
            ChadoFeature cf = new ChadoFeature(rs);
            Item geneticMarker = getChadoDBConverter().createItem("GeneticMarker");
            cf.populateBioEntity(geneticMarker, organism);
            geneticMarkerMap.put(new Integer(cf.feature_id), geneticMarker);
        }
        rs.close();
        LOG.info("Created "+geneticMarkerMap.size()+" GeneticMarker items.");

        // QTLS
        // load the QTLs from the feature table
        Map<Integer,Item> qtlMap = new HashMap<Integer,Item>();
        query = "SELECT * FROM feature WHERE organism_id="+organism_id+" AND type_id="+qtlCVTermId;
        LOG.info("executing query: "+query);
        rs = stmt.executeQuery(query);
        while (rs.next()) {
            ChadoFeature cf = new ChadoFeature(rs);
            Item qtl = getChadoDBConverter().createItem("QTL");
            cf.populateBioEntity(qtl, organism);
            qtlMap.put(new Integer(cf.feature_id), qtl);
        }
        rs.close();
        LOG.info("Created "+qtlMap.size()+" QTL items.");

        // GENETIC MAPS
        // genetic maps are not in the feature table, so we use linkage groups to query them for this organism
        // we'll assume the units are cM, although we can query them if we want to be really pedantic
        Map<Integer,Item> geneticMapMap = new HashMap<Integer,Item>();
        query = "SELECT * FROM featuremap WHERE featuremap_id IN (" +
            "SELECT DISTINCT featuremap_id FROM featurepos WHERE feature_id IN (" +
            "SELECT feature_id FROM feature WHERE type_id="+linkageGroupCVTermId+" AND organism_id="+organism_id+"))";
        LOG.info("executing query: "+query);
        rs = stmt.executeQuery(query);
        while (rs.next()) {
            // featuremap fields
            int featuremap_id = rs.getInt("featuremap_id");
            String name = rs.getString("name");
            String description = rs.getString("description");
            // build the GeneticMap item
            Item geneticMap = getChadoDBConverter().createItem("GeneticMap");
            geneticMap.setAttribute("primaryIdentifier", name);
            geneticMap.setAttribute("description", description);
            geneticMap.setAttribute("unit", "cM");
            geneticMap.setReference("organism", organism);
            // put it in its map for further processing
            geneticMapMap.put(new Integer(featuremap_id), geneticMap);
        }
        rs.close();
        LOG.info("Created "+geneticMapMap.size()+" GeneticMap items.");

        // GENE FAMILIES
        // load the gene families from the featureprop table for our organism - distinct values with gene family type_id
        // Note: we don't have IDs for gene families, so we use their name as the key
        Map<String,Item> geneFamilyMap = new HashMap<String,Item>();
        query = "SELECT DISTINCT featureprop.value FROM featureprop,feature " +
            "WHERE featureprop.type_id="+geneFamilyCVTermId+" AND featureprop.feature_id=feature.feature_id AND feature.organism_id="+organism_id;
        LOG.info("executing query: "+query);
        rs = stmt.executeQuery(query);
        while (rs.next()) {
            String value = rs.getString("value");
            Item geneFamily = getChadoDBConverter().createItem("GeneFamily");
            geneFamily.setAttribute("primaryIdentifier", value);
            geneFamilyMap.put(value, geneFamily);
        }
        rs.close();
        LOG.info("Created "+geneFamilyMap.size()+" GeneFamily items.");

        
        //---------------------------------------------------------
        // ---------------- DATA ASSOCIATION ----------------------
        //---------------------------------------------------------

        // query the pub table to retrieve and link Publications that are referenced by our GeneticMaps
        Map<Integer,Item> publicationMap = new HashMap<Integer,Item>();
        for (Map.Entry<Integer,Item> entry : geneticMapMap.entrySet()) {
            int featuremap_id = (int) entry.getKey();
            Item geneticMap = entry.getValue();
            query = "SELECT * FROM pub WHERE pub_id IN (SELECT pub_id FROM featuremap_pub WHERE featuremap_id="+featuremap_id+")";
            rs = stmt.executeQuery(query);
            while (rs.next()) {
                // extract the relevant pub attributes
                // NOTE: should have PubMed ID, but not given in chado (generally)!
                Item publication;
                int pub_id = rs.getInt("pub_id");
                if (publicationMap.containsKey(new Integer(pub_id))) {
                    // pub already exists, so link to the existing pub
                    publication = publicationMap.get(new Integer(pub_id));
                } else {
                    // new pub, so create it
                    String title = rs.getString("title");
                    String volume = rs.getString("volume");
                    String journal = rs.getString("series_name");
                    String issue = rs.getString("issue");
                    String year = rs.getString("pyear");
                    String pages = rs.getString("pages");
                    String uniquename = rs.getString("uniquename");
                    String[] parts = uniquename.split("et al");
                    if (parts.length==1) parts = parts[0].split(",");
                    String firstAuthor = parts[0];
                    int pubMedId = 0;
                    try {
                        // search for the PubMedId, hopefully find it with just first author
                        String[] authors = new String[1];
                        authors[0] = firstAuthor;
                        pubMedId = PubMedSearch.getPubMedId(journal, Integer.parseInt(year), authors);
                    } catch (Exception ex) {
                        // do nothing
                    }
                    // build the Publication item
                    publication = getChadoDBConverter().createItem("Publication");
                    if (title!=null && title.length()>0 && !title.equals("NULL")) publication.setAttribute("title", title);
                    if (volume!=null && volume.length()>0 && !volume.equals("NULL")) publication.setAttribute("volume", volume);
                    if (journal!=null && journal.length()>0 && !journal.equals("NULL")) publication.setAttribute("journal", journal);
                    if (issue!=null && issue.length()>0 && !issue.equals("NULL")) publication.setAttribute("issue", issue);
                    if (year!=null && year.length()>0 && !year.equals("NULL")) publication.setAttribute("year", year);
                    if (pages!=null && pages.length()>0 && !pages.equals("NULL")) publication.setAttribute("pages", pages);
                    if (firstAuthor!=null && firstAuthor.length()>0) publication.setAttribute("firstAuthor", firstAuthor);
                    if (pubMedId!=0) publication.setAttribute("pubMedId", String.valueOf(pubMedId));
                    // put it in its map for further processing
                    publicationMap.put(new Integer(pub_id), publication);
                }
                // add the GeneticMap reference to the Publication (reverse-referenced)
                publication.addToCollection("bioEntities", geneticMap);
            }
            rs.close();
        }
        LOG.info("***** Done creating "+publicationMap.size()+" Publication items associated with genetic maps.");

        // query the featureprop table to associate our genes with gene families
        for (Map.Entry<Integer,Item> entry : geneMap.entrySet()) {
            int feature_id = (int) entry.getKey(); // gene
            Item gene = entry.getValue();
            query = "SELECT * FROM featureprop WHERE feature_id="+feature_id;
            rs = stmt.executeQuery(query);
            while (rs.next()) {
                int featureprop_id = rs.getInt("featureprop_id"); // not used
                int type_id = rs.getInt("type_id");               // one of three possibilities
                String value = rs.getString("value");             // the description
                // branch on type_id
                if (type_id==2125) {
                    gene.setAttribute("note", value);
                } else if (type_id==78562) {
                    gene.setAttribute("familyRepresentative", value);
                } else if (type_id==51206) {
                    // add reference to the GeneFamily to this Gene (reversed-referenced)
                    Item geneFamily = geneFamilyMap.get(value);
                    if (geneFamily!=null) {
                        gene.setReference("geneFamily", geneFamily);
                    }
                }
            }
            rs.close();
        }
        LOG.info("***** Done associating genes with gene families from featureprop.");

        // query the featurepos table for linkage groups and genetic markers associated with our genetic maps
        // Note: for linkage groups, ignore the mappos=0 start entry; use mappos>0 entry to get the length of the linkage group
        for (Map.Entry<Integer,Item> entry : geneticMapMap.entrySet()) {
            int featuremap_id = (int) entry.getKey(); // genetic map
            Item geneticMap = entry.getValue();
            // genetic markers
            query = "SELECT * FROM featurepos WHERE featuremap_id="+featuremap_id+" AND feature_id<>map_feature_id";
            LOG.info("executing query: "+query);
            rs = stmt.executeQuery(query);
            while (rs.next()) {
                int feature_id = rs.getInt("feature_id"); // genetic marker
                int map_feature_id = rs.getInt("map_feature_id"); // linkage group
                double position = rs.getDouble("mappos");
                Item geneticMarker = geneticMarkerMap.get(new Integer(feature_id));
                Item linkageGroup = linkageGroupMap.get(new Integer(map_feature_id));
                if (geneticMarker!=null) {
                    geneticMap.addToCollection("geneticMarkers", geneticMarker);
                    if (linkageGroup!=null) {
                        linkageGroup.addToCollection("geneticMarkers", geneticMarker);
                        Item linkageGroupPosition = getChadoDBConverter().createItem("LinkageGroupPosition");
                        linkageGroupPosition.setAttribute("position", String.valueOf(position));
                        linkageGroupPosition.setReference("linkageGroup", linkageGroup);
                        store(linkageGroupPosition); // store it now, no further changes
                        geneticMarker.addToCollection("linkageGroupPositions", linkageGroupPosition);
                    }
                }
            }
            rs.close();
            // linkage groups
            query = "SELECT * FROM featurepos WHERE featuremap_id="+featuremap_id+" AND feature_id=map_feature_id AND mappos>0";
            LOG.info("executing query: "+query);
            rs = stmt.executeQuery(query);
            while (rs.next()) {
                int feature_id = rs.getInt("feature_id"); // linkage group
                double length = rs.getDouble("mappos");
                Item linkageGroup = linkageGroupMap.get(new Integer(feature_id));
                if (linkageGroup!=null) {
                    linkageGroup.setAttribute("length", String.valueOf(length));
                    linkageGroup.setReference("geneticMap", geneticMap);
                }
            }
            rs.close();
        }
        LOG.info("***** Done associating genetic maps and genetic markers and linkage groups from featurepos.");
        
        // query the featureloc table for the linkage group range (begin, end) per QTL
        // NOTE 1: Ethy hack: fmin,fmax are actually begin,end in cM x 100, so divide by 100 for cM!
        for (Map.Entry<Integer,Item> entry : qtlMap.entrySet()) {
            int feature_id = (int) entry.getKey(); // QTL in featureloc
            Item qtl = entry.getValue();
            // the Big Query to get the specific linkage group and flanking markers
            query = "SELECT * FROM featureloc WHERE feature_id="+feature_id;
            LOG.info("executing query: "+query);
            rs = stmt.executeQuery(query);
            while (rs.next()) {
                int linkage_group_id = rs.getInt("srcfeature_id");
                double begin = rs.getDouble("fmin")/100.0; // Ethy hack to store cM locations as integers
                double end = rs.getDouble("fmax")/100.0;
                Item linkageGroup = linkageGroupMap.get(new Integer(linkage_group_id));
                if (linkageGroup!=null) {
                    Item linkageGroupRange = getChadoDBConverter().createItem("LinkageGroupRange");
                    linkageGroupRange.setAttribute("begin", String.valueOf(begin));
                    linkageGroupRange.setAttribute("end", String.valueOf(end));
                    linkageGroupRange.setReference("linkageGroup", linkageGroup);
                    store(linkageGroupRange); // we're done with it
                    qtl.addToCollection("linkageGroupRanges", linkageGroupRange);
                    linkageGroup.addToCollection("QTLs", qtl);
                }
            }
            rs.close();
        }
        LOG.info("***** Done populating QTLs with linkage group ranges from featureloc.");

        // 1. query the feature_relationship table for direct relations between QTLs and genetic markers
        // 2. query featureloc for genes that overlap the full range spanned by the genetic markers found in 1
        for (Map.Entry<Integer,Item> entry : qtlMap.entrySet()) {
            int subject_id = (int) entry.getKey(); // QTL in feature_relationship
            Item qtl = entry.getValue();
            int[] geneticMarkerIds = new int[3]; // there is one (nearest marker) or three (nearest, flanking low, flanking high)
            int k = 0; // index of geneticMarkerIds
            query = "SELECT * FROM feature_relationship WHERE subject_id="+subject_id;
            rs = stmt.executeQuery(query);
            while (rs.next()) {
                int object_id = rs.getInt("object_id"); // genetic marker
                int type_id = rs.getInt("type_id"); // Nearest, Flanking Low, Flanking High; we don't store this at present
                Item geneticMarker = geneticMarkerMap.get(new Integer(object_id));
                if (geneticMarker!=null) {
                    // add to QTL's associatedGeneticMarkers collection; reverse-reference will automatically associate qtl in genetic marker's QTLs collection
                    qtl.addToCollection("associatedGeneticMarkers", geneticMarker);
                    // store genetic marker for gene finding to come
                    geneticMarkerIds[k++] = object_id;
                }
            }
            rs.close(); // done collecting genetic markers

            // now query ranges of these markers to get genes that overlap
            int f1 = 200000000; // initialize minimum of range
            int f2 = 0;         // initialize maximum of range
            int srcfeature_id = 0; // chromosome
            for (int l=0; l<k; l++) {
                query = "SELECT * FROM featureloc WHERE feature_id="+geneticMarkerIds[l];
                rs = stmt.executeQuery(query);
                while (rs.next()) {
                    int fmin = rs.getInt("fmin");
                    int fmax = rs.getInt("fmax");
                    srcfeature_id = rs.getInt("srcfeature_id"); // should be same chromosome for all genetic markers
                    if (fmin<f1) f1 = fmin;
                    if (fmax>f2) f2 = fmax;
                }
                rs.close();
                if (f1<f2) {
                    query = "SELECT featureloc.feature_id FROM featureloc,feature WHERE " +
                        "featureloc.feature_id=feature.feature_id AND type_id="+geneCVTermId+" AND srcfeature_id="+srcfeature_id+" AND fmin<"+f2+" AND fmax>"+f1;
                    rs = stmt.executeQuery(query);
                    while (rs.next()) {
                        int feature_id = rs.getInt("feature_id");
                        Item gene = geneMap.get(new Integer(feature_id));
                        if (gene!=null) {
                            qtl.addToCollection("overlappingGenes", gene);
                        }
                    }
                    rs.close();
                }
            }
            
        } // QTL loop
        LOG.info("***** Done populating QTLs with genetic markers from feature_relationship and overlapping genes from featureloc.");
        
            
        // ----------------------------------------------------------------
        // we're done, so store everything
        // ----------------------------------------------------------------

        LOG.info("Storing genetic maps...");
        for (Item item : geneticMapMap.values()) {
            store(item);
        }

        LOG.info("Storing genetic map publications...");
        for (Item item : publicationMap.values()) {
            store(item);
        }

        LOG.info("Storing genes...");
        for (Item item : geneMap.values()) {
            store(item);
        }

        LOG.info("Storing gene families...");
        for (Item item : geneFamilyMap.values()) {
            store(item);
        }

        LOG.info("Storing linkage groups...");
        for (Item item : linkageGroupMap.values()) {
            store(item);
        }

        LOG.info("Storing genetic markers...");
        for (Item item : geneticMarkerMap.values()) {
            store(item);
        }
 
        LOG.info("Storing QTLs...");
        for (Item item : qtlMap.values()) {
            store(item);
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

}
