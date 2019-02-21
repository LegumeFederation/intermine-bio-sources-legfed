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
import java.util.Map;
import java.util.HashMap;
import java.util.Set;
import java.util.HashSet;

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
        
        // store stuff in maps to avoid duplication
        Map<String,Item> linkageGroupMap = new HashMap<>();
        Map<String,Item> markerMap = new HashMap<>();
        Map<String,Item> qtlMap = new HashMap<>();
        Map<String,Item> geneticMapMap = new HashMap<>();
        Map<String,Item> chromosomeMap = new HashMap<>();
        Map<Integer,Item> publicationMap = new HashMap<>();
        Map<String,Item> strainMap = new HashMap<>();

        // initialize our DB statement
        Statement stmt = connection.createStatement();
        ResultSet rs; // general-purpose use
        String query; // general-purpose use
        
        // set the cv term IDs for our items of interest
        int linkageGroupTypeId = getTypeId(stmt, "linkage_group");
        int geneticMarkerTypeId = getTypeId(stmt, "genetic_marker");
        int qtlTypeId = getTypeId(stmt, "QTL");
        int consensusRegionTypeId = getTypeId(stmt, "consensus_region");
        int favorableAlleleSourceTypeId = getTypeId(stmt, "Favorable Allele Source");
        int nearestMarkerTypeId = getTypeId(stmt, "Nearest Marker");
        int flankingMarkerLowTypeId = getTypeId(stmt, "Flanking Marker Low");
        int flankingMarkerHighTypeId = getTypeId(stmt, "Flanking Marker High");
        
        // get the desired chado organism_ids
        Set<Integer> organismIds = getChadoDBConverter().getDesiredChadoOrganismIds();

        // loop over the chado organisms to fill the Item maps 
        for (Integer organismId : organismIds) {
            int organism_id = organismId.intValue();
            Item organism = getChadoDBConverter().getOrganismItem(organismId);
            Item strain = getChadoDBConverter().getStrainItem(organismId);
            // store strains in a map since we'll get extra ones from QTL.favorableAlleleSource
            strainMap.put(getChadoDBConverter().getStrainName(organismId), strain);

            // query the featureloc table for marker chromosome locations
            query = "SELECT C.name AS cname, C.uniquename AS cuniquename, " +
                "M.name AS mname, M.uniquename AS muniquename, " +
                "featureloc.fmin AS start, featureloc.fmax AS end " +
                "FROM feature C, feature M, featureloc " +
                "WHERE C.organism_id="+organism_id+" " +
                "AND   M.type_id="+geneticMarkerTypeId+" " +
                "AND   featureloc.feature_id=M.feature_id " +
                "AND   featureloc.srcfeature_id=C.feature_id";
            rs = stmt.executeQuery(query);
            while (rs.next()) {
                String cname = rs.getString("cname");
                String cuniquename = rs.getString("cuniquename");
                String mname = rs.getString("mname");
                String muniquename = rs.getString("muniquename");
                int start = rs.getInt("start");
                int end = rs.getInt("end");
                Item chromosome;
                if (chromosomeMap.containsKey(cname)) {
                    chromosome = chromosomeMap.get(cname);
                } else {
                    // assume no markers positioned on scaffolds!
                    chromosome = createItem("Chromosome");
                    chromosome.setAttribute("primaryIdentifier", cuniquename);
                    chromosome.setAttribute("secondaryIdentifier", cname);
                    chromosome.setReference("organism", organism);
                    chromosome.setReference("strain", strain);
                    chromosomeMap.put(cname, chromosome);
                }
                Item marker;
                if (markerMap.containsKey(mname)) {
                    marker = markerMap.get(mname);
                } else {
                    marker = createItem("GeneticMarker");
                    marker.setAttribute("primaryIdentifier", mname);
                    marker.setReference("organism", organism);
                    marker.setReference("strain", strain);
                    markerMap.put(mname, marker);
                }
                // create the marker's location
                marker.setReference("chromosome", chromosome);
                Item location = createItem("Location");
                location.setAttribute("start", String.valueOf(start));
                location.setAttribute("end", String.valueOf(end));
                location.setReference("feature", marker);
                location.setReference("locatedOn", chromosome);
                store(location);
                marker.setReference("chromosomeLocation", location);
            }
            rs.close();

            // query the feature_relationship table for marker-QTL relations
            query = "SELECT DISTINCT Q.name AS qname, M.name AS mname " +
                "FROM feature_relationship, feature Q, feature M " +
                "WHERE Q.organism_id="+organism_id+" " +
                "AND feature_relationship.type_id IN ("+nearestMarkerTypeId+","+flankingMarkerLowTypeId+","+flankingMarkerHighTypeId+") " +
                "AND feature_relationship.subject_id=Q.feature_id " +
                "AND feature_relationship.object_id=M.feature_id";
            rs = stmt.executeQuery(query);
            while (rs.next()) {
                String qname = rs.getString("qname");
                String mname = rs.getString("mname");
                Item qtl;
                if (qtlMap.containsKey(qname)) {
                    qtl = qtlMap.get(qname);
                } else {
                    qtl = createItem("QTL");
                    qtl.setAttribute("primaryIdentifier", qname);
                    qtl.setReference("organism", organism);
                    qtlMap.put(qname, qtl);
                }
                Item marker;
                if (markerMap.containsKey(mname)) {
                    marker = markerMap.get(mname);
                } else {
                    marker = createItem("GeneticMarker");
                    marker.setAttribute("primaryIdentifier", mname);
                    marker.setReference("organism", organism);
                    markerMap.put(mname, marker);
                }
                qtl.addToCollection("markers", marker);
            }
            rs.close();

            // query the featurepos table for marker-linkage group positions
            query = "SELECT DISTINCT L.name AS lname, mappos, M.name AS mname " +
                "FROM featurepos, feature L, feature M " +
                "WHERE L.organism_id="+organism_id+" " +
                "AND   featurepos.map_feature_id=L.feature_id " +
                "AND   featurepos.feature_id=M.feature_id";
            rs = stmt.executeQuery(query);
            while (rs.next()) {
                String lname = rs.getString("lname");
                String mname = rs.getString("mname");
                double mappos = round(rs.getDouble("mappos"),3);
                Item linkageGroup;
                if (linkageGroupMap.containsKey(lname)) {
                    linkageGroup = linkageGroupMap.get(lname);
                } else {
                    linkageGroup = createItem("LinkageGroup");
                    linkageGroup.setAttribute("primaryIdentifier", lname);
                    linkageGroup.setReference("organism", organism);
                    linkageGroupMap.put(lname, linkageGroup);
                }
                if (lname.equals(mname)) {
                    // linkage group length is given by lname==mname and mappos>0
                    if (mappos>0.0) {
                        linkageGroup.setAttribute("length", String.valueOf(mappos));
                    }
                } else {
                    Item marker;
                    if (markerMap.containsKey(mname)) {
                        marker = markerMap.get(mname);
                    } else {
                        marker = createItem("GeneticMarker");
                        marker.setAttribute("primaryIdentifier", mname);
                        marker.setReference("organism", organism);
                        markerMap.put(mname, marker);
                    }
                    linkageGroup.addToCollection("markers", marker);
                    // set the linkage group position of the marker
                    Item linkageGroupPosition = createItem("LinkageGroupPosition");
                    linkageGroupPosition.setAttribute("position", String.valueOf(mappos));
                    linkageGroupPosition.setReference("linkageGroup", linkageGroup);
                    store(linkageGroupPosition);
                    marker.addToCollection("linkageGroupPositions", linkageGroupPosition);
                }
            }
            rs.close();

            // query the featuremap table for genetic maps that are actually associated with features
            query = "SELECT DISTINCT featuremap.name, featuremap.description " +
                "FROM featuremap,featurepos,feature " +
                "WHERE featuremap.featuremap_id=featurepos.featuremap_id " +
                "AND   featurepos.feature_id=feature.feature_id " +
                "AND   feature.organism_id="+organism_id;
            rs = stmt.executeQuery(query);
            while (rs.next()) {
                String name = rs.getString("name");
                String description = rs.getString("description");
                Item geneticMap = createItem("GeneticMap");
                geneticMap.setAttribute("primaryIdentifier", name);
                geneticMap.setAttribute("description", description);
                geneticMap.setAttribute("unit", "cM");
                geneticMap.setReference("organism", organism);
                geneticMapMap.put(name, geneticMap);
            }
            rs.close();

            // query the featurepos and featuremap tables to get the genetic maps and their linkage groups
            query = "SELECT featuremap.name AS geneticmap, feature.name AS linkagegroup " +
                "FROM feature,featurepos,featuremap " +
                "WHERE feature.feature_id=featurepos.feature_id " +
                "AND   featurepos.featuremap_id=featuremap.featuremap_id " +
                "AND   mappos>0 " +
                "AND   feature.type_id="+linkageGroupTypeId+" " +
                "AND   feature.organism_id="+organism_id;
            rs = stmt.executeQuery(query);
            while (rs.next()) {
                String geneticmap = rs.getString("geneticmap");
                String linkagegroup = rs.getString("linkagegroup");
                if (geneticMapMap.containsKey(geneticmap)) {
                    Item geneticMap = geneticMapMap.get(geneticmap);
                    Item linkageGroup;
                    if (linkageGroupMap.containsKey(linkagegroup)) {
                        linkageGroup = linkageGroupMap.get(linkagegroup);
                    } else {
                        linkageGroup = createItem("LinkageGroup");
                        linkageGroup.setAttribute("primaryIdentifier", linkagegroup);
                        linkageGroup.setReference("organism", organism);
                        linkageGroupMap.put(linkagegroup, linkageGroup);
                    }
                    linkageGroup.setReference("geneticMap", geneticMap);
                }
            }

            // query the featureloc table for QTL ranges on linkage groups, remembering to divide by 100!
            query = "SELECT L.name AS lname, Q.name AS qname, featureloc.fmin AS begin, featureloc.fmax AS end " +
                "FROM feature L, feature Q, featureloc " +
                "WHERE Q.feature_id=featureloc.feature_id " +
                "AND   L.feature_id=featureloc.srcfeature_id " +
                "AND   Q.type_id="+qtlTypeId+" " +
                "AND   L.type_id="+linkageGroupTypeId+" " +
                "AND   L.organism_id="+organism_id;
            rs = stmt.executeQuery(query);
            while (rs.next()) {
                String lname = rs.getString("lname");
                String qname = rs.getString("qname");
                double begin = rs.getDouble("begin")/100.0;
                double end = rs.getDouble("end")/100.0;
                if (linkageGroupMap.containsKey(lname)) {
                    Item linkageGroup = linkageGroupMap.get(lname);
                    Item qtl;
                    if (qtlMap.containsKey(qname)) {
                        qtl = qtlMap.get(qname);
                    } else {
                        qtl = createItem("QTL");
                        qtl.setAttribute("primaryIdentifier", qname);
                        qtl.setReference("organism", organism);
                        qtlMap.put(qname, qtl);
                    }
                    Item linkageGroupRange = createItem("LinkageGroupRange");
                    linkageGroupRange.setAttribute("begin", String.valueOf(round(begin,2)));
                    linkageGroupRange.setAttribute("end", String.valueOf(round(end,2)));
                    linkageGroupRange.setAttribute("length", String.valueOf(round(end-begin,2)));
                    linkageGroupRange.setReference("linkageGroup", linkageGroup);
                    store(linkageGroupRange);
                    qtl.addToCollection("linkageGroupRanges", linkageGroupRange);
                    linkageGroup.addToCollection("QTLs", qtl); // redundant???
                }
            }
            rs.close();

            // query the pub table to retrieve and link Publications that are referenced by our QTLs via feature_cvterm
            query = "SELECT feature.name,pub.* FROM feature,feature_cvterm,pub " +
                "WHERE feature.type_id="+qtlTypeId+" " +
                "AND feature.organism_id="+organism_id+" " +
                "AND feature.feature_id=feature_cvterm.feature_id " +
                "AND feature_cvterm.pub_id=pub.pub_id " +
                "AND volume!='NULL' AND pyear!='submitted'";
            rs = stmt.executeQuery(query);
            while (rs.next()) {
                String name = rs.getString("name");
                int pub_id = rs.getInt("pub_id");
                // only add publications that are associated with the QTLs we've pulled above
                if (qtlMap.containsKey(name)) {
                    Item qtl = qtlMap.get(name);
                    Item publication;
                    if (publicationMap.containsKey(pub_id)) {
                        publication = publicationMap.get(pub_id);
                    } else {
                        publication = createItem("Publication");
                        setPublicationAttributes(publication, rs);
                        store(publication);
                        publicationMap.put(pub_id, publication);
                    }
                    qtl.addToCollection("publications", publication);
                }
            }
            rs.close();

            // query the feature_stock table for QTL favorable allele source
            query = "SELECT feature.name,stock.name AS stockname FROM feature,feature_stock,stock " +
                "WHERE feature.feature_id=feature_stock.feature_id " +
                "AND feature.organism_id="+organism_id+" " +
                "AND feature_stock.stock_id=stock.stock_id " +
                "AND feature_stock.type_id="+favorableAlleleSourceTypeId;
            rs = stmt.executeQuery(query);
            while (rs.next()) {
                String name = rs.getString("name");
                String stockname = rs.getString("stockname");
                if (qtlMap.containsKey(name)) {
                    Item qtl = qtlMap.get(name);
                    Item stock;
                    if (strainMap.containsKey(stockname)) {
                        stock = strainMap.get(stockname);
                    } else {
                        stock = createItem("Strain");
                        stock.setAttribute("primaryIdentifier", stockname);
                        stock.setReference("organism", organism);
                        store(stock);
                        strainMap.put(stockname, stock);
                    }
                    qtl.setReference("favorableAlleleSource", stock);
                }
            }
            rs.close();
        }
        
        // query the pub table to retrieve and link Publications that are referenced by our genetic maps in featuremap_pub
        query = "SELECT featuremap.name,pub.* FROM pub,featuremap_pub,featuremap WHERE pub.pub_id=featuremap_pub.pub_id " +
            "AND featuremap_pub.featuremap_id=featuremap.featuremap_id " +
            "AND volume!='NULL' AND pyear!='submitted'";
        rs = stmt.executeQuery(query);
        while (rs.next()) {
            String name = rs.getString("name");
            int pub_id = rs.getInt("pub_id");
            // only add publications that are associated with the genetic maps we've pulled above
            if (geneticMapMap.containsKey(name)) {
                Item geneticMap = geneticMapMap.get(name);
                Item publication;
                if (publicationMap.containsKey(pub_id)) {
                    publication = publicationMap.get(pub_id);
                } else {
                    publication = createItem("Publication");
                    setPublicationAttributes(publication, rs);
                    store(publication);
                    publicationMap.put(pub_id, publication);
                }
                geneticMap.addToCollection("publications", publication);
            }
        }
        rs.close();
            
        // store the maps
        LOG.info("Storing "+chromosomeMap.size()+" chromosomes...");
        store(chromosomeMap.values());
        
        LOG.info("Storing "+linkageGroupMap.size()+" linkage groups...");
        store(linkageGroupMap.values());

        LOG.info("Storing "+markerMap.size()+" genetic markers...");
        store(markerMap.values());
 
        LOG.info("Storing "+qtlMap.size()+" QTLs...");
        store(qtlMap.values());

        LOG.info("Storing "+geneticMapMap.size()+" genetic maps...");
        store(geneticMapMap.values());
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
