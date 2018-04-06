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
 * A few tab-separated metadata fields precede the data, and should be IN THIS ORDER:
 * <pre>
 * TaxonID      3920
 * MappingPopulation	3WayCross-2018
 * Description  Blah di blah blah blah blah blah.
 * MatrixNotes  A = Parent A, B = Parent B. Lower case: genotype calls reversed based on parental alleles.
 * Parent	BigTall123
 * Parent	SmallShort456
 * Parent	MediumTall789
 * PMID	12345678
 * PMID	23456789
 * MarkerType   SNP
 * GeneticMap   3WayCross-2018
 * Lines	Line1	Line2	Line3	...
 * </pre>
 *
 * @author Sam Hokin
 */
public class GTFileConverter extends BioFileConverter {

    static final int MAX_MARKERS = 0; // limit the number of markers per file for testing purposes; 0 to disable limiting
	
    private static final Logger LOG = Logger.getLogger(GTFileConverter.class);

    // store items at end in close() method, they may be duplicated across files
    Map<String,Item> organismMap = new HashMap<String,Item>();   // keyed by taxonId, variety (parent)
    Map<String,Item> publicationMap = new HashMap<String,Item>();
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
        if (getCurrentFile().getName().contains("README")) return;

        // these objects are created per file
        String taxonId = null;
        String markerType = null;
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

            String[] parts = line.split("\t");
            if (!line.startsWith("#") && parts.length>1) {

                String key = parts[0].trim();
                String value = parts[1].trim();

                if (key.toLowerCase().equals("taxonid")) {

                    taxonId = value;

                } else if (key.toLowerCase().equals("mappingpopulation")) {

                    mappingPopulation = createItem("MappingPopulation");
                    mappingPopulation.setAttribute("primaryIdentifier", value);

                } else if (key.toLowerCase().equals("description")) {

                    mappingPopulation.setAttribute("description", value);

                } else if (key.toLowerCase().equals("matrixnotes")) {

                    mappingPopulation.setAttribute("matrixNotes", value);
                
                } else if (key.toLowerCase().equals("geneticmap")) {

                    Item geneticMap = createItem("GeneticMap");
                    geneticMap.setAttribute("primaryIdentifier", value);
                    store(geneticMap);
                    LOG.info("Stored GeneticMap "+value);
                    mappingPopulation.setReference("geneticMap", geneticMap);

                } else if (key.toLowerCase().equals("markertype")) {

                    markerType = value;

                } else if (key.toLowerCase().equals("parent")) {

                    // add parent to mapping population's collection
                    Item organism;
                    if (organismMap.containsKey(value)) {
                        organism = organismMap.get(value);
                    } else {
                        organism = createItem("Organism");
                        organism.setAttribute("taxonId", taxonId);
                        organism.setAttribute("variety", value);
                        organismMap.put(value, organism);
                    }
                    mappingPopulation.addToCollection("parents", organism);
                    
                } else if (key.toLowerCase().equals("pmid")) {

                    // add a related publication to this mapping population
                    String pubMedId = value;
                    if (publicationMap.containsKey(pubMedId)) {
                        Item publication = publicationMap.get(pubMedId);
                        mappingPopulation.addToCollection("publications", publication);
                    } else {
                        // create a new publication
                        Item publication = createItem("Publication");
                        publication.setAttribute("pubMedId", pubMedId);
                        publicationMap.put(pubMedId, publication);
                        mappingPopulation.addToCollection("publications", publication);
                    }

                } else if (key.toLowerCase().equals("lines")) {

                    int num = parts.length - 1;    // number of genotyping lines
                    lines = new Item[num];         // CB27/BB-007, etc.
                    for (int i=0; i<num; i++) {
                        String lineName = parts[i+1];
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
                        // this line should only appear with this mapping population
                        lines[i].setReference("mappingPopulation", mappingPopulation);
                        store(lines[i]);
                    }

                } else {

                    // add a genotyping record for a marker
                    GTRecord gtr = new GTRecord(line);
                    Item marker;
                    if (markerMap.containsKey(gtr.marker)) {
                        marker = markerMap.get(gtr.marker);
                    } else {
                        marker = createItem("GeneticMarker");
                        marker.setAttribute("primaryIdentifier", gtr.marker);
                        marker.setAttribute("type", markerType);
                        markerMap.put(gtr.marker, marker);
                    }
                    marker.addToCollection("mappingPopulations", mappingPopulation);
                
                    // create and store this marker's genotypeValues; they should be same order/number as lines.
                    for (int i=0; i<lines.length; i++) {
                        Item genotypeValue = createItem("GenotypeValue");
                        genotypeValue.setAttribute("value", gtr.values[i]);
                        genotypeValue.setReference("line", lines[i]);
                        genotypeValue.setReference("marker", marker);
                        store(genotypeValue);
                    }
                    
                }

            }

        } // end while reading lines

        // each file has a unique mapping population so store it here
        store(mappingPopulation);
        
        // wrap up this file
        br.close();
 
    }


    /**
     * Store the items stored in global maps
     */
    @Override
    public void close() throws Exception {
        store(organismMap.values());
        store(markerMap.values());
        store(publicationMap.values());
    }

}
