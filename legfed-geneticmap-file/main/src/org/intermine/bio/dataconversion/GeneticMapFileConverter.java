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

import java.math.BigDecimal;
import java.math.RoundingMode;

import java.util.List;
import java.util.Map;
import java.util.HashMap;
import java.util.Set;
import java.util.HashSet;

import org.apache.log4j.Logger;

import org.ncgr.intermine.PublicationTools;
import org.ncgr.pubmed.PubMedSummary;

import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Item;

/**
 * Store the genetic marker, QTL and linkage group data for a genetic map stored in a tab-separated file. Data fields are:
 * <pre>
 * Marker	Type	Linkage group	Position (cM)	QTL	Traits associated
 * </pre>
 * The QTL name is put in the QTL.primaryIdentifier, the Traits associated are put in the QTL.secondaryIdentifier.
 *
 * A few tab-separated metadata fields precede the data, and should be in this order:
 * <pre>
 * TaxonID	3917
 * GeneticMap	Size-Drought-2019
 * PMID	12345678
 * Parents	BigTall123	SmallShort456
 * Parents	DResist234	DSuscept567
 * </pre>
 *
 * The file name gives name of the genetic map, e.g. Consensus-2019.tsv.
 *
 * The "Linkage group" column provides the numerical linkage group (chromosome, normally, but can be any LG). This is stored as the LG number.
 * The LG primary identifier is formed by prepending the genetic map name. So, for example, the metadata above would lead to a linkage group
 * with primary identifier "Size-Drought-2019_2" and secondary identifier "2".
 *
 * Markers do not get secondary identifiers. They have a unique key (primaryIdentifier,organism). Their type (e.g. SNP) is given in the second column.
 *
 * QTL and Traits associated are NOT required.
 *
 * @author Sam Hokin
 */
public class GeneticMapFileConverter extends BioFileConverter {
	
    private static final Logger LOG = Logger.getLogger(GeneticMapFileConverter.class);

    // items saved at end are stored in maps
    Map<String,Item> organismMap = new HashMap<String,Item>();          // keyed by taxonID
    Map<String,Item> geneticMapMap = new HashMap<String,Item>();        // keyed by primaryIdentifier
    Map<String,Item> germplasmMap = new HashMap<String,Item>();         // keyed by primaryIdentifier
    Map<String,Item> mappingPopulationMap = new HashMap<String,Item>(); // keyed by primaryIdentifier
    Map<String,Item> linkageGroupMap = new HashMap<String,Item>();      // keyed by primaryIdentifier
    Map<String,Item> markerMap = new HashMap<String,Item>();            // keyed by primaryIdentifier
    Map<String,Item> qtlMap = new HashMap<String,Item>();               // keyed by primaryIdentifier
    Map<String,Item> linkageGroupRangeMap = new HashMap<String,Item>(); // keyed by qtl.primaryIdentifier
    Map<String,Item> publicationMap = new HashMap<String,Item>();       // keyed by pubMedId
    Map<String,Item> authorMap = new HashMap<String,Item>();            // keyed by name
    /**
     * Create a new GeneticMapFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public GeneticMapFileConverter(ItemWriter writer, Model model) {
        super(writer, model);
    }

    /**
     * {@inheritDoc}
     * Read in the LinkageGroup file and store the linkage groups, QTLs and genetic markers along with their ranges and positions
     */
    @Override
    public void process(Reader reader) throws Exception {

        // don't process README files
        if (getCurrentFile().getName().equals("README")) return;

        LOG.info("Processing LinkageGroup file "+getCurrentFile().getName()+"...");

        // these objects are created per file
        Item organism = null;
        Item geneticMap = null;
        String geneticMapName = null;

        // -----------------------------------------------------------------------------------------------------------------------------------
        // Load genetic markers, linkage groups and QTLs from the LinkageGroup file and associate QTLs and genetic markers with linkage groups
        // -----------------------------------------------------------------------------------------------------------------------------------

        BufferedReader br = new BufferedReader(reader);
        String line = null;
        while ((line=br.readLine())!=null) {

            if (line.startsWith("TaxonID")) {

                // get the organism
                String[] parts = line.split("\t");
                String taxonID = parts[1];
                if (organismMap.containsKey(taxonID)) {
                    organism = organismMap.get(taxonID);
                } else {
                    organism = createItem("Organism");
                    BioStoreHook.setSOTerm(this, organism, "organism", getSequenceOntologyRefId());
                    organism.setAttribute("taxonId", taxonID);
                    organismMap.put(taxonID, organism);
                }
                
            } else if (line.startsWith("GeneticMap")) {

                // get the genetic map
                String[] parts = line.split("\t");
                geneticMapName = parts[1];
                if (geneticMapMap.containsKey(geneticMapName)) {
                    geneticMap = geneticMapMap.get(geneticMapName);
                } else {
                    geneticMap = createItem("GeneticMap");
                    BioStoreHook.setSOTerm(this, geneticMap, "organism", getSequenceOntologyRefId());
                    geneticMap.setAttribute("primaryIdentifier", geneticMapName);
                    geneticMap.setReference("organism", organism);
                    geneticMapMap.put(geneticMapName, geneticMap);
                }
                
            } else if (line.startsWith("PMID")) {

                // get the Publication if PMID present
                String[] parts = line.split("\t");
                String pubMedId = parts[1];
                if (publicationMap.containsKey(pubMedId)) {
                    Item publication = publicationMap.get(pubMedId);
                    geneticMap.addToCollection("publications", publication);
                } else {
                    Item publication = PublicationTools.getPublicationFromPMID(this, Integer.parseInt(pubMedId));
                    if (publication!=null) {
                        publicationMap.put(pubMedId, publication);
                        List<Item> authors = PublicationTools.getAuthorsFromPMID(this, Integer.parseInt(pubMedId));
                        for (Item author : authors) {
                            String name = author.getAttribute("name").getValue();
                            if (authorMap.containsKey(name)) {
                                Item authorToStore = authorMap.get(name);
                                publication.addToCollection("authors", authorToStore);
                            } else {
                                authorMap.put(name, author);
                                publication.addToCollection("authors", author);
                            }
                        }
                        geneticMap.addToCollection("publications", publication);
                    }
                }

            } else if (line.startsWith("Parents")) {

                // get Germplasm parents
                String[] parts = line.split("\t");
                String[] parents = new String[2];
                parents[0] = parts[1];
                parents[1] = parts[2];
                // create MappingPopulation and name for parents
                String mappingPopulationName = parents[0]+"_x_"+parents[1];
                if (mappingPopulationMap.containsKey(mappingPopulationName)) {
                    Item mappingPopulation = mappingPopulationMap.get(mappingPopulationName);
                    geneticMap.addToCollection("mappingPopulations", mappingPopulation);
                } else {
                    Item mappingPopulation = createItem("MappingPopulation");
                    mappingPopulation.setAttribute("primaryIdentifier", mappingPopulationName);
                    // add parents to collection
                    for (int i=0; i<2; i++) {
                        if (germplasmMap.containsKey(parents[i])) {
                            Item germplasm = germplasmMap.get(parents[i]);
                            mappingPopulation.addToCollection("parents", germplasm);
                        } else {
                            Item germplasm = createItem("Germplasm");
                            germplasm.setAttribute("primaryIdentifier", parents[i]);
                            germplasm.setReference("organism", organism);
                            germplasmMap.put(parents[i], germplasm);
                            mappingPopulation.addToCollection("parents", germplasm);
                        }
                    }
		    mappingPopulationMap.put(mappingPopulationName, mappingPopulation);
                    geneticMap.addToCollection("mappingPopulations", mappingPopulation);
                }

            } else {

                // looks like it's a data record
                GeneticMapRecord record = new GeneticMapRecord(line);

                // construct a fuller LinkageGroup primary identifier from genetic map name and LG number
                String lgID = geneticMapName+"_"+record.lg;
                
                // create this linkage group if new
                Item linkageGroup = null;
                if (linkageGroupMap.containsKey(lgID)) {
                    linkageGroup = linkageGroupMap.get(lgID);
                } else {
                    linkageGroup = createItem("LinkageGroup");
                    BioStoreHook.setSOTerm(this, linkageGroup, "linkage_group", getSequenceOntologyRefId());
                    linkageGroup.setAttribute("primaryIdentifier", lgID);
                    linkageGroup.setAttribute("number", String.valueOf(record.lg));
                    linkageGroup.setAttribute("length", "0.0"); // initialize with zero, will be updated by marker positions
                    linkageGroup.setReference("organism", organism);
                    linkageGroup.setReference("geneticMap", geneticMap);
                    linkageGroupMap.put(lgID, linkageGroup);
                }

                // create and store this genetic marker if new
                Item marker;
                if (markerMap.containsKey(record.marker)) {
                    marker = markerMap.get(record.marker);
                } else {
                    marker = createItem("GeneticMarker");
                    BioStoreHook.setSOTerm(this, marker, "genetic_marker", getSequenceOntologyRefId());
                    marker.setReference("organism", organism);
                    marker.setAttribute("primaryIdentifier", record.marker);
                    marker.setAttribute("type", record.type);
                    // create and store this marker's linkage group position
                    Item linkageGroupPosition = createItem("LinkageGroupPosition");
                    linkageGroupPosition.setAttribute("position", String.valueOf(record.position));
                    linkageGroupPosition.setReference("linkageGroup", linkageGroup);
                    store(linkageGroupPosition);
                    marker.addToCollection("linkageGroupPositions", linkageGroupPosition);
                    // add this marker to this genetic map and linkage group collections
                    geneticMap.addToCollection("markers", marker);
                    linkageGroup.addToCollection("markers", marker);
                    markerMap.put(record.marker, marker);
                }

                // update this linkage group's length from current marker position
                double length = Double.parseDouble(linkageGroup.getAttribute("length").getValue());
                if (record.position>length) linkageGroup.setAttribute("length", String.valueOf(record.position));
                
                if (record.qtl!=null) {
                
                    // create and add this QTL to this linkage group if new
                    Item qtl = null;
                    Item linkageGroupRange = null;
                    String key = record.qtl;
                    if (qtlMap.containsKey(key)) {
                        qtl = qtlMap.get(key);
                        linkageGroupRange = linkageGroupRangeMap.get(key);
                    } else {
                        // create new QTL 
                        qtl = createItem("QTL");
                        BioStoreHook.setSOTerm(this, qtl, "QTL", getSequenceOntologyRefId());
                        qtl.setReference("organism", organism);
                        qtl.setAttribute("primaryIdentifier", record.qtl);
                        if (record.traits!=null) {
                            qtl.setAttribute("secondaryIdentifier", record.traits);
                        }
                        // create this QTL's LinkageGroupRange, initialized at zero length
                        linkageGroupRange = createItem("LinkageGroupRange");
                        linkageGroupRange.setAttribute("begin", String.valueOf(record.position)); // start=end with first marker
                        linkageGroupRange.setAttribute("end", String.valueOf(record.position));
                        linkageGroupRange.setAttribute("length", "0.0");
                        linkageGroupRange.setReference("linkageGroup", linkageGroup);
                        // add this LinkageGroupRange to this QTLs collection
                        qtl.addToCollection("linkageGroupRanges", linkageGroupRange);
                        // store this QTL and LinkageGroupRange in their maps
                        qtlMap.put(key, qtl);
                        linkageGroupRangeMap.put(key, linkageGroupRange);
                        // add this QTL to this GeneticMap and this LinkageGroup collections
                        geneticMap.addToCollection("QTLs", qtl); // this is redundant storage, but handy all the same
                        linkageGroup.addToCollection("QTLs", qtl);
                    }
                    
                    // add this marker to this QTL's collection
                    qtl.addToCollection("associatedMarkers", marker);
                    
                    // update this QTL's linkage group range from the current marker position
                    double begin = Double.parseDouble(linkageGroupRange.getAttribute("begin").getValue());
                    double end = Double.parseDouble(linkageGroupRange.getAttribute("end").getValue());
                    if (record.position<begin) {
                        begin = record.position;
                    } else if (record.position>end) {
                        end = record.position;
                    }
                    linkageGroupRange.setAttribute("begin", String.valueOf(begin));
                    linkageGroupRange.setAttribute("end", String.valueOf(end));
                    linkageGroupRange.setAttribute("length", String.valueOf(round(end-begin,2)));

                }

            }

        }
        
        br.close();
 
    }


    /**
     * Store the items we've collected from the LinkageGroup files
     */
    @Override
    public void close() throws Exception {
        store(organismMap.values());
        store(geneticMapMap.values());
        store(germplasmMap.values());
        store(mappingPopulationMap.values());
        store(linkageGroupMap.values());
        store(markerMap.values());
        store(qtlMap.values());
        store(linkageGroupRangeMap.values());
        store(authorMap.values());
        store(publicationMap.values());
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
