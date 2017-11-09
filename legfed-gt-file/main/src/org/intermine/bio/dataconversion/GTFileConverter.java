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
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import org.ncgr.intermine.PubMedPublication;

import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Item;

/**
 * Store the organisms, markers, lines and genotype values from a MappingPopulation genotyping file.
 * <pre>
 *         line1 line2 line3 ...
 * marker1     A     a     B ...
 *   ...
 * </pre>
 *
 * A few tab-separated metadata fields precede the data, and should be in this order:
 * <pre>
 * TaxonID      3920
 * Parents	BigTall123	SmallShort456
 * PMID	12345678
 * PMID	23456789
 * Lines	Line1	Line2	Line3	...
 * </pre>
 *
 * @author Sam Hokin
 */
public class GTFileConverter extends BioFileConverter {

    static final int MAX_MARKERS = 0; // limit the number of markers per file for testing purposes; 0 to disable limiting
	
    private static final Logger LOG = Logger.getLogger(GTFileConverter.class);

    // store items at end in close() method, they may be duplicated across files
    Map<String,Item> organismMap = new HashMap<String,Item>();   // keyed by variety
    Map<String,Item> publicationMap = new HashMap<String,Item>();
    Map<String,Item> authorMap = new HashMap<String,Item>();
    Map<String,Item> lineMap = new HashMap<String,Item>();
    Map<String,Item> markerMap = new HashMap<String,Item>();

    /**
     * Create a new GTFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public GTFileConverter(ItemWriter writer, Model model) {
        super(writer, model);
    }

    /**
     * {@inheritDoc}
     * Read in each genotype file and parse the organism, lines, markers, and genotypes
     */
    @Override
    public void process(Reader reader) throws Exception {

        // don't process README files
        if (getCurrentFile().getName().equals("README")) return;

        // these objects are created per file
        String taxonId = null;
        Item mappingPopulation = null;
        Item[] lines = null;

        LOG.info("Processing GT file "+getCurrentFile().getName()+"...");
        
        // ----------------------------------------------------------
        // lines, genetic markers, genotypes from the GT file
        // ----------------------------------------------------------

        BufferedReader br = new BufferedReader(reader);

        int count = 0;
        String line = null;
        while ((line=br.readLine())!=null) {

            if (MAX_MARKERS>0 && count>MAX_MARKERS) break;  // IF TESTING: limit the number of markers
            count++; 

            if (line.startsWith("TaxonID")) {

                String[] parts = line.split("\t");
                taxonId = parts[1];

            } else if (line.startsWith("Parents")) {

                // get parent organisms
                String[] parts = line.split("\t");
                String[] parents = new String[2];
                parents[0] = parts[1];
                parents[1] = parts[2];
                // create MappingPopulation and name for parents
                mappingPopulation = createItem("MappingPopulation");
                mappingPopulation.setAttribute("primaryIdentifier", parents[0]+"_x_"+parents[1]); // this is how we know parent A and parent B
                // add parents to collection
                for (int i=0; i<2; i++) {
                    Item organism;
                    if (organismMap.containsKey(parents[i])) {
                        organism = organismMap.get(parents[i]);
                    } else {
                        organism = createItem("Organism");
                        organism.setAttribute("taxonId", taxonId);
                        organism.setAttribute("variety", parents[i]);
                        organismMap.put(parents[i], organism);
                    }
                    mappingPopulation.addToCollection("parents", organism);
                }
                
            } else if (line.startsWith("PMID")) {

                // get the Publication if PMID present
                String[] parts = line.split("\t");
                String pubMedId = parts[1];
                if (publicationMap.containsKey(pubMedId)) {
                    Item publication = publicationMap.get(pubMedId);
                    mappingPopulation.addToCollection("publications", publication);
                } else {
                    PubMedPublication pubMedPub = new PubMedPublication(this, Integer.parseInt(pubMedId));
                    // DEBUG
                    LOG.info(pubMedPub.getSummary().toString());
                    Item publication = pubMedPub.getPublication();
                    if (publication!=null) {
                        List<Item> authors = pubMedPub.getAuthors();
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
                        publicationMap.put(pubMedId, publication);
                        mappingPopulation.addToCollection("publications", publication);
                    }
                }

            } else if (line.startsWith("Lines")) {

                // get the lines
                String[] parts = line.split("\t");
                int num = parts.length - 1;    // number of genotyping lines
                lines = new Item[num];         // CB27/BB-007, etc.
                for (int i=0; i<num; i++) {
                    String lineName = parts[i+1];
                    if (lineMap.containsKey(lineName)) {
                        lines[i] = lineMap.get(lineName);
                    } else {
                        lines[i] = createItem("GenotypingLine");
                        lines[i].setAttribute("primaryIdentifier", lineName);
                        // try to parse out a number, e.g. 17 from CB27-17, for ordering purposes
                        String[] chunks = lineName.split("-");
                        if (chunks.length>1) {
                            String str = chunks[chunks.length-1]; // supports Foo-Bar-17
                            try {
                                Integer number = new Integer(Integer.parseInt(str));
                                lines[i].setAttribute("number", String.valueOf(number));
                            } catch (NumberFormatException ex) {
                                // do nothing, it's not an integer
                            }
                        }
                        lineMap.put(lineName, lines[i]);
                    }
                    lines[i].setReference("mappingPopulation", mappingPopulation);
                }

            } else {

                // add a genotyping record for a marker
                GTRecord gt = new GTRecord(line);
                Item marker;
                if (markerMap.containsKey(gt.marker)) {
                    marker = markerMap.get(gt.marker);
                } else {
                    marker = createItem("GeneticMarker");
                    marker.setAttribute("primaryIdentifier", gt.marker);
                    marker.setAttribute("type", "SNP"); // TODO: put the marker type in the file header, they don't have to be SNPs
                    markerMap.put(gt.marker, marker);
                }
                marker.addToCollection("mappingPopulations", mappingPopulation);
                
                // create and store this marker's genotypeValues; they should be same order/number as lines.
                for (int i=0; i<lines.length; i++) {
                    Item genotypeValue = createItem("GenotypeValue");
                    genotypeValue.setAttribute("value", gt.values[i]);
                    genotypeValue.setReference("line", lines[i]);
                    genotypeValue.setReference("marker", marker);
                    store(genotypeValue);
                }

            }

        } // end while reading a file

        // each file has a unique mapping population so store it here
        store(mappingPopulation);
        
        // wrap up
        br.close();
 
    }


    /**
     * Store the items stored in global maps
     */
    @Override
    public void close() throws Exception {
        store(organismMap.values());
        store(lineMap.values());
        store(markerMap.values());
        store(authorMap.values());
        store(publicationMap.values());
    }

}
