package org.intermine.bio.dataconversion;

/**
 * Copyright (C) 2015-2016 NCGR
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  See the LICENSE file for more
 * information or http://www.gnu.org/copyleft/lesser.html.
 *
 */

import java.math.BigDecimal;
import java.math.RoundingMode;
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
 * These come from feature, featuremap, featurepos, featureloc, feature_relationship, feature_stock and stock.
 *
 * Since this processer deals only with chado data, Items are stored in maps with Integer keys equal to
 * the chado feature.feature_id.
 *
 * @author Sam Hokin
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
        
        // initialize our DB statement
        Statement stmt = connection.createStatement();
        ResultSet rs; // general-purpose use
        String query; // general-purpose use
        
        // ---------------------------------------------------------
        // ---------------- INITIAL DATA LOADING -------------------
        // ---------------------------------------------------------
        
        // set the cv term IDs for our items of interest
        int linkageGroupTypeId = getTypeId(stmt, "linkage_group");
        int geneticMarkerTypeId = getTypeId(stmt, "genetic_marker");
        int qtlTypeId = getTypeId(stmt, "QTL");
        int consensusRegionTypeId = getTypeId(stmt, "consensus_region");
        int favorableAlleleSourceTypeId = getTypeId(stmt, "Favorable Allele Source");

        // store stuff in maps to avoid duplication
        Map<Integer,Item> organismMap = new HashMap<Integer,Item>();
        Map<Integer,Item> linkageGroupMap = new HashMap<Integer,Item>();
        Map<Integer,Item> geneticMarkerMap = new HashMap<Integer,Item>();
        Map<Integer,Item> qtlMap = new HashMap<Integer,Item>();
        Map<Integer,Item> geneticMapMap = new HashMap<Integer,Item>();
        Map<Integer,Item> publicationMap = new HashMap<Integer,Item>();
        
        // build the Organism map from the supplied taxon IDs, with optional variety suffixes, e.g. 3827_desi
        Map<Integer,OrganismData> chadoToOrgData = getChadoDBConverter().getChadoIdToOrgDataMap();
        for (Integer organismId : chadoToOrgData.keySet()) {
            OrganismData organismData = chadoToOrgData.get(organismId);
            String taxonId = organismData.getTaxonId();
            String variety = organismData.getVariety();
            if (variety==null) variety = OrganismData.DEFAULT_VARIETY; // need to provide a non-null default for merging
            Item organism = getChadoDBConverter().createItem("Organism");
            organism.setAttribute("taxonId", String.valueOf(taxonId));
            organism.setAttribute("variety", variety);
            store(organism);
            organismMap.put(organismId, organism);
        }
        LOG.info("Created and stored "+organismMap.size()+" organism Items.");

        // create a comma-separated list of organism IDs for WHERE organism_id IN query clauses
        String organismSQL = "(";
        boolean first = true;
        for (Integer organismId : organismMap.keySet()) {
            if (first) {
                organismSQL += organismId;
                first = false;
            } else {
                organismSQL += ","+organismId;
            }
        }
        organismSQL += ")";

        // loop over the organisms to fill the Item maps 
        for (Integer organismId : organismMap.keySet()) {
            int organism_id = organismId.intValue();
            Item organism = organismMap.get(organismId);

            // append the maps from the feature table
            linkageGroupMap.putAll(generateMap(stmt, "LinkageGroup", linkageGroupTypeId, organism_id));
	    qtlMap.putAll(generateMap(stmt, "QTL", qtlTypeId, organism_id));
            geneticMarkerMap.putAll(generateMap(stmt, "GeneticMarker", geneticMarkerTypeId, organism_id, organism));

            // genetic maps are not in the feature table, so we use linkage groups to query them for this organism
            // we'll assume the units are cM, although we can query them if we want to be really pedantic
            query = "SELECT * FROM featuremap WHERE featuremap_id IN (" +
                "SELECT DISTINCT featuremap_id FROM featurepos WHERE feature_id IN (" +
                "SELECT feature_id FROM feature WHERE organism_id="+organism_id+" AND type_id="+linkageGroupTypeId+"))";
            rs = stmt.executeQuery(query);
            while (rs.next()) {
                // featuremap fields
                int featuremap_id = rs.getInt("featuremap_id");
                String name = rs.getString("name");
                String description = rs.getString("description");
                // build the GeneticMap item
                Item geneticMap = getChadoDBConverter().createItem("GeneticMap");
                // genetic map has no SO term
                geneticMap.setAttribute("primaryIdentifier", name);
                geneticMap.setAttribute("description", description);
                geneticMap.setAttribute("unit", "cM");
                geneticMapMap.put(new Integer(featuremap_id), geneticMap);
            }
            rs.close();
        }

        //---------------------------------------------------------
        // ---------------- DATA ASSOCIATION ----------------------
        //---------------------------------------------------------

        // query the pub table to retrieve and link Publications that are referenced by our genetic maps in featuremap_pub
        for (Integer featuremapId : geneticMapMap.keySet()) {
            int featuremap_id = (int) featuremapId;
            Item geneticMap = geneticMapMap.get(featuremapId);
            ResultSet rs1 = stmt.executeQuery("SELECT * FROM pub WHERE pub_id IN (SELECT pub_id FROM featuremap_pub WHERE featuremap_id="+featuremap_id+")");
            while (rs1.next()) {
                int pub_id = rs1.getInt("pub_id");
                Integer pubId = new Integer(pub_id);
                Item publication;
                if (publicationMap.containsKey(pubId)) {
                    // pub already exists, so get from the map
                    publication = publicationMap.get(pubId);
                } else {
                    // create it, set the attributes from the ResultSet, add it to its map
                    publication = getChadoDBConverter().createItem("Publication");
                    setPublicationAttributes(publication, rs1);
                    store(publication);
                    publicationMap.put(pubId, publication);
                }
                // add this publication reference to our genetic map
                geneticMap.addToCollection("publications", publication);
            }
            rs1.close();
        }
        LOG.info("***** Done creating Publication items associated with genetic maps.");

        // query the pub table to retrieve and link Publications that are referenced by our QTLs in feature_cvterm
        for (Integer featureId : qtlMap.keySet()) {
            int feature_id = (int) featureId;
            Item qtl = qtlMap.get(featureId);
            rs = stmt.executeQuery("SELECT * FROM pub WHERE pub_id IN (SELECT pub_id FROM feature_cvterm WHERE feature_id="+feature_id+")");
            while (rs.next()) {
                int pub_id = rs.getInt("pub_id");
                Integer pubId = new Integer(pub_id);
                Item publication;
                if (publicationMap.containsKey(pubId)) {
                    // pub already exists, so get from the map
                    publication = publicationMap.get(pubId);
                } else {
                    // create it, set the attributes from the ResultSet, add it to its map
                    publication = getChadoDBConverter().createItem("Publication");
                    setPublicationAttributes(publication, rs);
                    store(publication);
                    publicationMap.put(pubId, publication);
                }
                // add this publication to our QTL.publications collection
                qtl.addToCollection("publications", publication);
            }
            rs.close();
        }
        LOG.info("***** Done creating Publication items associated with QTLs.");

        // query the feature_stock and stock tables to retrieve the favorable allele source per QTL
        for (Integer featureId : qtlMap.keySet()) {
            int feature_id = (int) featureId;
            Item qtl = qtlMap.get(featureId);
            rs = stmt.executeQuery("SELECT * FROM stock WHERE stock_id=(SELECT stock_id FROM feature_stock WHERE type_id="+favorableAlleleSourceTypeId+" AND feature_id="+feature_id+")");
            if (rs.next()) {
                String favorableAlleleSource = rs.getString("uniquename");
                qtl.setAttribute("favorableAlleleSource", favorableAlleleSource);
            }
            rs.close();
        }
        LOG.info("***** Done setting QTL.favorableAlleleSource attributes. *****");
        
        // query the featurepos table for linkage groups and genetic markers associated with our genetic maps
        // Note: for linkage groups, ignore the mappos=0 start entry; use mappos>0 entry to get the length of the linkage group
        for (Integer featuremapId : geneticMapMap.keySet()) {
            int featuremap_id = (int) featuremapId;
            Item geneticMap = geneticMapMap.get(featuremapId);
            // genetic markers
            ResultSet rs1 = stmt.executeQuery("SELECT * FROM featurepos WHERE featuremap_id="+featuremap_id+" AND feature_id<>map_feature_id");
            while (rs1.next()) {
                int feature_id = rs1.getInt("feature_id"); // genetic marker
                int map_feature_id = rs1.getInt("map_feature_id"); // linkage group
                double position = rs1.getDouble("mappos");
                Item geneticMarker = geneticMarkerMap.get(new Integer(feature_id));
                Item linkageGroup = linkageGroupMap.get(new Integer(map_feature_id));
                if (geneticMarker!=null) {
                    geneticMap.addToCollection("markers", geneticMarker);
                    if (linkageGroup!=null) {
                        linkageGroup.addToCollection("markers", geneticMarker);
                        Item linkageGroupPosition = getChadoDBConverter().createItem("LinkageGroupPosition");
                        linkageGroupPosition.setAttribute("position", String.valueOf(position));
                        linkageGroupPosition.setReference("linkageGroup", linkageGroup);
                        store(linkageGroupPosition); // store it now, no further changes
                        geneticMarker.addToCollection("linkageGroupPositions", linkageGroupPosition);
                    }
                }
            }
            rs1.close();
            // linkage groups
            ResultSet rs2 = stmt.executeQuery("SELECT * FROM featurepos WHERE featuremap_id="+featuremap_id+" AND feature_id=map_feature_id AND mappos>0");
            while (rs2.next()) {
                int feature_id = rs2.getInt("feature_id"); // linkage group
                double length = rs2.getDouble("mappos");
                Item linkageGroup = linkageGroupMap.get(new Integer(feature_id));
                if (linkageGroup!=null) {
                    linkageGroup.setAttribute("length", String.valueOf(length));
                    linkageGroup.setReference("geneticMap", geneticMap);
                }
            }
            rs2.close();
        }
        LOG.info("***** Done associating genetic maps and genetic markers and linkage groups from featurepos.");
        
        // query the featureloc table for the linkage group range (begin, end) per QTL
        // NOTE 1: Ethy hack: fmin,fmax are actually begin,end in cM x 100, so divide by 100 for cM!
        for (Integer featureId : qtlMap.keySet()) {
            int feature_id = (int) featureId;
            Item qtl = qtlMap.get(featureId);
            // the Big Query to get the specific linkage group and flanking markers
            ResultSet rs1 = stmt.executeQuery("SELECT * FROM featureloc WHERE feature_id="+feature_id);
            while (rs1.next()) {
                int linkage_group_id = rs1.getInt("srcfeature_id");
                double begin = rs1.getDouble("fmin")/100.0; // Ethy hack to store cM locations as integers
                double end = rs1.getDouble("fmax")/100.0;
                Item linkageGroup = linkageGroupMap.get(new Integer(linkage_group_id));
                if (linkageGroup!=null) {
                    Item linkageGroupRange = getChadoDBConverter().createItem("LinkageGroupRange");
                    linkageGroupRange.setAttribute("begin", String.valueOf(begin));
                    linkageGroupRange.setAttribute("end", String.valueOf(end));
                    linkageGroupRange.setAttribute("length", String.valueOf(round(end-begin,2)));
                    linkageGroupRange.setReference("linkageGroup", linkageGroup);
                    store(linkageGroupRange); // we're done with it
                    qtl.addToCollection("linkageGroupRanges", linkageGroupRange);
                    linkageGroup.addToCollection("QTLs", qtl);
                }
            }
            rs1.close();
        }
        LOG.info("***** Done populating QTLs with linkage group ranges from featureloc.");

        // Query the feature_relationship table for direct relations between QTLs and genetic markers
        for (Integer subjectId : qtlMap.keySet()) {
            int subject_id = (int) subjectId; // QTL in feature_relationship
            Item qtl = qtlMap.get(subjectId);
            // int[] geneticMarkerIds = new int[3]; // there is one (nearest marker) or three (nearest, flanking low, flanking high)
            // int k = 0; // index of geneticMarkerIds
            ResultSet rs1 = stmt.executeQuery("SELECT * FROM feature_relationship WHERE subject_id="+subject_id);
            while (rs1.next()) {
                int object_id = rs1.getInt("object_id"); // genetic marker
                int type_id = rs1.getInt("type_id"); // Nearest, Flanking Low, Flanking High; we don't store this at present
                Item geneticMarker = geneticMarkerMap.get(new Integer(object_id));
                if (geneticMarker!=null) qtl.addToCollection("markers", geneticMarker);
            }
            rs1.close(); // done collecting genetic markers
        }
        LOG.info("***** Done populating QTLs with genetic markers from feature_relationship.");
            
        // ----------------------------------------------------------------
        // we're done, so store the things that were bulk loaded
        // ----------------------------------------------------------------

        LOG.info("Storing "+linkageGroupMap.size()+" linkage groups...");
        for (Item item : linkageGroupMap.values()) store(item);

        LOG.info("Storing "+geneticMarkerMap.size()+" genetic markers...");
        for (Item item : geneticMarkerMap.values()) store(item);
 
        LOG.info("Storing "+qtlMap.size()+" QTLs...");
        for (Item item : qtlMap.values()) store(item);

        LOG.info("Storing "+geneticMapMap.size()+" genetic maps...");
        for (Item item : geneticMapMap.values()) store(item);

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
     * Set the attributes of the given publication Item from the given pub ResultSet.
     */
    void setPublicationAttributes(Item publication, ResultSet rs) throws SQLException {
        String title = rs.getString("title");
        String volume = rs.getString("volume");
        String journal = rs.getString("series_name");
        String issue = rs.getString("issue");
        String year = rs.getString("pyear");
        String pages = rs.getString("pages");
        String uniquename = rs.getString("uniquename");
        String[] parts = uniquename.split(",");
        String firstAuthor = parts[0];
        // build the Publication item
        if (title!=null && title.length()>0 && !title.equals("NULL") && !title.equals("0")) publication.setAttribute("title", title);
        if (volume!=null && volume.length()>0 && !volume.equals("NULL") && !volume.equals("0")) publication.setAttribute("volume", volume);
        if (journal!=null && journal.length()>0 && !journal.equals("NULL") && !journal.equals("0")) publication.setAttribute("journal", journal);
        if (issue!=null && issue.length()>0 && !issue.equals("NULL") && !issue.equals("0")) publication.setAttribute("issue", issue);
	if (year!=null && year.length()>0 && !year.equals("NULL")) {
	    // only set the year if it's an integer
	    try {
		int y = Integer.parseInt(year);
		publication.setAttribute("year", year);
	    } catch (Exception e) {
		LOG.info("Error parsing publication year:"+year);
	    }
	}
        if (pages!=null && pages.length()>0 && !pages.equals("NULL")) publication.setAttribute("pages", pages);
        if (firstAuthor!=null && firstAuthor.length()>0) publication.setAttribute("firstAuthor", firstAuthor);
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

    /**
     * Generate a map, keyed with feature_id, with records from the chado feature table corresponding to the given organism_id and type_id.
     *
     * @param stmt an SQL Statement object
     * @param className the IM class to be populated
     * @param typeId the CV term type_id of the chado feature
     * @param organism_id the chado organism_id for the features of interest
     */
    Map<Integer,Item> generateMap(Statement stmt, String className, int typeId, int organism_id) throws SQLException {
	String query = "SELECT * FROM feature WHERE organism_id="+organism_id+" AND type_id="+typeId;
        ResultSet rs = stmt.executeQuery(query);
        Map<Integer,Item> map = new HashMap<Integer,Item>();
	int count = 0;
        while (rs.next()) {
	    count++;
            Item item = getChadoDBConverter().createItem(className);
            ChadoFeature cf = new ChadoFeature(rs);
            cf.populateItem(item);
            map.put(new Integer(cf.feature_id), item);
        }
        rs.close();
        return map;
    }
    
    /**
     * Generate a map, keyed with feature_id, with records from the chado feature table corresponding to the given organism_id and type_id.
     * This version also sets a reference to the given organism.
     *
     * @param stmt an SQL Statement object
     * @param className the IM class to be populated
     * @param typeId the CV term type_id of the chado feature
     * @param organism_id the chado organism_id
     * @param organism the corresponding IM organism to associate with the feature
     */
    Map<Integer,Item> generateMap(Statement stmt, String className, int typeId, int organism_id, Item organism) throws SQLException {
	String query = "SELECT * FROM feature WHERE organism_id="+organism_id+" AND type_id="+typeId;
        ResultSet rs = stmt.executeQuery(query);
        Map<Integer,Item> map = new HashMap<Integer,Item>();
	int count = 0;
        while (rs.next()) {
	    count++;
            Item item = getChadoDBConverter().createItem(className);
            ChadoFeature cf = new ChadoFeature(rs);
            cf.populateItem(item, organism);
            map.put(new Integer(cf.feature_id), item);
        }
        rs.close();
        return map;
    }

    /**
     * Get the CV term type_id for a given CV term name.
     * @param stmt the database connection statement, initialized to the chado database
     * @param name the desired CV term name
     * @return the type_id
     * @throws SQLException
     */
     protected int getTypeId(Statement stmt, String name) throws SQLException {
        ResultSet rs = stmt.executeQuery("SELECT cvterm_id FROM cvterm WHERE name='"+name+"'");
        if (rs.next()) {
            int typeId = rs.getInt("cvterm_id");
            rs.close();
            return typeId;
        } else {
            throw new RuntimeException("Could not determine CV term type_id for '"+name+"'.");
        }
     }

}
