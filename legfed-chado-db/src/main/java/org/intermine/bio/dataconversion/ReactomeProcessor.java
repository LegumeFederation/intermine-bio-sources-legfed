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
import org.apache.log4j.LogManager;

import org.intermine.bio.util.OrganismData;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Attribute;
import org.intermine.xml.full.Item;
import org.intermine.xml.full.Reference;

/**
 * Read plant reactome data in from a tab-delimited file, load pathways and associate them with genes.
 * We do NOT store pathways that are not associated with any of our genes.
 * The file fields are:
 * <pre>
 * identifier      name                    species           gene
 * R-CCN-1119506.1 tyrosine degradation I  Coffea canephora  Cc05_g11510
 * </pre>
 * The file is given by the parameter reactome.file in project.xml. It is likely gene_ids_by_pathway_and_species.tab from Gramene.
 *
 * NOTE: special species-specific treatment is required, as the Plant Reactome gene names do NOT match the LIS gene names by a long stretch.
 *
 * @author Sam Hokin, NCGR
 */
public class ReactomeProcessor extends ChadoProcessor {
	
    private static final Logger LOG = LogManager.getLogger(ReactomeProcessor.class);

    Map<String,Item> organismMap = new HashMap<>(); // keyed by "Genus species"

    /**
     * Create a new ReactomeProcessor
     * @param chadoDBConverter the ChadoDBConverter that is controlling this processor
     */
    public ReactomeProcessor(ChadoDBConverter chadoDBConverter) {
        super(chadoDBConverter);
    }

    /**
     * {@inheritDoc}
     * Process by spinning through the reactome file and associating pathways for our organisms with genes from chado.
     */
    @Override
    public void process(Connection connection) throws Exception {
        // Items stored in maps to avoid dupe stores
        Map<Integer,Item> geneMap = new HashMap<Integer,Item>();  // keyed by feature_id
        Map<String,Item> pathwayMap = new HashMap<String,Item>();    // keyed by identifier
        
        // initialize our DB statement
        Statement stmt = connection.createStatement();
        ResultSet rs;

        // get the desired chado organism_ids
        Set<Integer> organismIds = getChadoDBConverter().getDesiredChadoOrganismIds();

        // store desired species strings in a Set, hopefully spelled same way!
        Set<String> speciesSet = new HashSet<String>();              // "Phaseolus vulgaris"
        for (Integer organismId : organismIds) {
            OrganismData od = getChadoDBConverter().getChadoIdToOrgDataMap().get(organismId);
            String species = od.getGenus()+" "+od.getSpecies();
            speciesSet.add(species);
            Item organism = getChadoDBConverter().getOrganismItem(organismId);
            organismMap.put(species, organism);
        }

        // get the cv term ID for gene
        rs = stmt.executeQuery("SELECT cvterm_id FROM cvterm WHERE name='gene'");
        rs.next();
        int geneCVTermId = rs.getInt("cvterm_id");
        rs.close();

        // load the reactome file into a map
        String reactomeFilename = getChadoDBConverter().getReactomeFilename();
        LOG.info("Reading "+reactomeFilename);
        BufferedReader reader = new BufferedReader(new FileReader(reactomeFilename));
        String line = null;
        while ((line=reader.readLine())!=null) {
            if (!line.startsWith("#")) {
                String[] fields = line.split("\t");
                String identifier = fields[0];
                String pathwayName = fields[1];
                String species = fields[2];
                String geneName = fields[3];
                if (speciesSet.contains(species)) {
                    String featureName = null;
                    Item organism = organismMap.get(species);
                    // HACK: concoct the correct chado feature.name
                    if (geneName.contains("_")) {
                        String[] parts = geneName.split("_");
                        if (parts[0].equals("GLYMA")) {
                            featureName = "glyma.Glyma."+parts[1];
                        } else if (parts[0].equals("PHAVU")) {
                            featureName = "phavu.Phvul."+parts[1];
                            featureName = featureName.substring(0, featureName.length()-1); // drop the trailing "g"
                        } else if (parts[0].equals("MTR")) {
                            featureName = "medtr.Medtr"+parts[1];
                        } else if (parts[0].equals("C.cajan")) {
                            featureName = "cajca.C.cajan_"+parts[1];
                        } else if (parts[0].equals("TanjilG")) {
                            featureName = "lupan.Lup0"+parts[1];
                        } else if (parts[0].equals("Tp57577")) {
                            featureName = "tripr."+parts[3];
                        }
                    } else if (geneName.startsWith("Aradu")) {
                        featureName = "aradu."+geneName;
                    } else if (geneName.startsWith("Araip")) {
                        featureName = "araip."+geneName;
                    }
                    if (featureName!=null) {
                        String query = "SELECT * FROM feature WHERE type_id="+geneCVTermId+" AND name='"+featureName+"'";
                        rs = stmt.executeQuery(query);
                        if (rs.next()) {
                            // get or store the gene
                            Item gene = null;
                            int feature_id = rs.getInt("feature_id");
                            if (geneMap.containsKey(feature_id)) {
                                gene = geneMap.get(feature_id);
                            } else {
                                gene = getChadoDBConverter().createItem("Gene");
                                gene.setAttribute("chadoId", String.valueOf(feature_id));
                                gene.setAttribute("primaryIdentifier", rs.getString("uniquename"));
                                gene.setAttribute("secondaryIdentifier", rs.getString("name"));
                                gene.setReference("organism", organism);
                                geneMap.put(feature_id, gene);
                            }
                            // get or store the pathway, with reference to organism
                            Item pathway = null;
                            if (pathwayMap.containsKey(identifier)) {
                                pathway = pathwayMap.get(identifier);
                            } else {
                                pathway = getChadoDBConverter().createItem("Pathway");
                                pathway.setAttribute("identifier", identifier);
                                pathway.setAttribute("name", pathwayName);
                                pathway.setReference("organism", organism);
                                pathwayMap.put(identifier, pathway);
                            }
                            // associate gene with pathway
                            gene.addToCollection("pathways", pathway);
                        } else {
                            LOG.error(featureName+" not found in chado database.");
                        }
                        rs.close();
                    }
                }
            }
        }

        // store the maps here so that the collection get stored
        store(geneMap.values());
        store(pathwayMap.values());
    }
}
