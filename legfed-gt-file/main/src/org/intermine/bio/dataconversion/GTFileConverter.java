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

import org.ncgr.intermine.PublicationTools;
import org.ncgr.pubmed.PubMedSummary;

import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Item;

/**
 * Store the germplasms, markers, lines and genotype values from a MappingPopulation genotyping file.
 * <pre>
 *         line1 line2 line3 ...
 * marker1     A     a     B ...
 *   ...
 * </pre>
 *
 * A few tab-separated metadata fields precede the data, and should be in this order:
 * <pre>
 * TaxonID	3917
 * Parents	BigTall123	SmallShort456
 * PMID	12345678
 * Lines	Line1	Line2	Line3	...
 * </pre>
 *
 * @author Sam Hokin
 */
public class GTFileConverter extends BioFileConverter {

    static final int MAX_MARKERS = 0; // limit the number of markers per file for testing purposes; 0 to disable limiting
	
    private static final Logger LOG = Logger.getLogger(GTFileConverter.class);

    // store items at end in close() method, they may be duplicated across files
    Map<String,Item> organismMap = new HashMap<String,Item>();
    Map<String,Item> publicationMap = new HashMap<String,Item>();
    Map<String,Item> lineMap = new HashMap<String,Item>();
    Map<String,Item> markerMap = new HashMap<String,Item>();
    Map<String,Item> germplasmMap = new HashMap<String,Item>();

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
     * Read in each genotype file and parse the germplasms, organism, lines, markers, and genotypes
     */
    @Override
    public void process(Reader reader) throws Exception {

        // don't process README files
        if (getCurrentFile().getName().equals("README")) return;

        // these objects are created per file
        Item organism = null;
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

            } else if (line.startsWith("Parents")) {

                // create Germplasm parents and the MappingPopulation
                String[] parts = line.split("\t");
                String parentA = parts[1];
                String parentB = parts[2];
                Item germplasmA;
                Item germplasmB;
                if (germplasmMap.containsKey(parentA)) {
                    germplasmA = germplasmMap.get(parentA);
                } else {
                    germplasmA = createItem("Germplasm");
                    germplasmA.setAttribute("primaryIdentifier", parentA);
                    germplasmA.setReference("organism", organism);
                    germplasmMap.put(parentA, germplasmA);
                }
                if (germplasmMap.containsKey(parentB)) {
                    germplasmB = germplasmMap.get(parentB);
                } else {
                    germplasmB = createItem("Germplasm");
                    germplasmB.setAttribute("primaryIdentifier", parentB);
                    germplasmB.setReference("organism", organism);
                    germplasmMap.put(parentB, germplasmB);
                }
                mappingPopulation = createItem("MappingPopulation");
                mappingPopulation.setAttribute("primaryIdentifier", parentA+"_x_"+parentB);
                mappingPopulation.setReference("parentA", germplasmA);
                mappingPopulation.setReference("parentB", germplasmB);

            } else if (line.startsWith("PMID")) {

                // get the Publication if PMID present
                String[] parts = line.split("\t");
                String pubMedId = parts[1];
                Item publication = PublicationTools.storePublicationFromPMID(this, Integer.parseInt(pubMedId));
                if (publication!=null) mappingPopulation.setReference("publication", publication);

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
                    marker.setAttribute("type", "SNP");
                    marker.setReference("organism", organism);
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
        store(germplasmMap.values());
        store(lineMap.values());
        store(markerMap.values());
    }

}
