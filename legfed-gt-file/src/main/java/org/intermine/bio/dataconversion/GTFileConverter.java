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

import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Item;

/**
 * Store the markers, lines and genotype values from a genotyping file as well as some metadata on the genotyping study.
 *
 * <pre>
 * Lines       line1 line2 line3 ...
 * marker1     A     a     B ...
 *   ...
 * </pre>
 *
 * OR
 *
 * <pre>
 * Markers marker1 marker2 marker3 ...
 * line1   A       ...
 * line2   a       ...
 * line3   B       ...
 * </pre>
 *
 * If the last line before value data starts with "Lines" then the first format is expected; if the last line before data starts with "Markers" then
 * second (Flapjack) format is expected.
 *
 * A few tab-separated metadata fields precede the data IN THIS ORDER:
 * <pre>
 * TaxonID          3387
 * GenotypingStudy  3WayCross-2018
 * Description      Blah di blah blah blah blah blah.
 * MatrixNotes      A = Parent A, B = Parent B. Lower case: genotype calls reversed based on parental alleles.
 * MarkerType       SNP
 * PMID (optional)  12345678
 * PMID (optional)  23456789
 * Lines            line1 line2 line3 ...
 *   or
 * Markers          marker1 marker2 marker3
 * </pre>
 *
 * @author Sam Hokin
 */
public class GTFileConverter extends BioFileConverter {

    private static final Logger LOG = Logger.getLogger(GTFileConverter.class);

    private static final int MAX_ROWS = 0; // limit the number of rows read per file for testing purposes; 0 to disable limiting
	
    // store items at end in close() method, they may be duplicated across files
    Map<String,Item> publicationMap = new HashMap<String,Item>();
    Map<String,Item> markerMap = new HashMap<String,Item>();
    Map<String,Item> lineMap = new HashMap<String,Item>();
    Map<String,Item> strainMap = new HashMap<String,Item>();
    Map<String,Item> organismMap = new HashMap<String,Item>();

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
     * Read in each genotype file and parse the metadata, lines, markers, and genotypes
     */
    @Override
    public void process(Reader reader) throws Exception {

        // don't process README files
        if (getCurrentFile().getName().contains("README")) return;

        // these objects are created per file
        String markerType = null;
        Item genotypingStudy = null;

        // we'll load either an array of lines or an array of markers
        Item[] lines = new Item[0];
        Item[] markers = new Item[0];
        Item organism = null;
        boolean rowsAreLines = false;
        boolean rowsAreMarkers = false;

        LOG.info("Processing GT file "+getCurrentFile().getName()+"...");
        
        // ----------------------------------------------------------
        // lines, genetic markers, genotypes from the GT file
        // ----------------------------------------------------------

        BufferedReader br = new BufferedReader(reader);

        int count = 0;
        String row = null;
        while ((row=br.readLine())!=null) {
            if (MAX_ROWS>0 && count>MAX_ROWS) break;  // IF TESTING: limit the number of markers
            count++;

            String[] parts = row.split("\t");
            if (!row.startsWith("#") && parts.length>1) {

                String key = parts[0].trim();
                String value = parts[1].trim();

                if (key.toLowerCase().equals("taxonid")) {

                    if (organismMap.containsKey(value)) {
                        organism = organismMap.get(value);
                    } else {
                        organism = createItem("Organism");
                        organism.setAttribute("taxonId", value);
                        store(organism);
                        organismMap.put(value, organism);
                    }

                } else if (key.toLowerCase().equals("genotypingstudy")) {
                                        
                    genotypingStudy = createItem("GenotypingStudy");
                    genotypingStudy.setAttribute("primaryIdentifier", value);
                    genotypingStudy.setReference("organism", organism);

                } else if (key.toLowerCase().equals("description")) {

                    genotypingStudy.setAttribute("description", value);

                } else if (key.toLowerCase().equals("matrixnotes")) {

                    genotypingStudy.setAttribute("matrixNotes", value);
                
                } else if (key.toLowerCase().equals("markertype")) {

                    markerType = value;

                } else if (key.toLowerCase().equals("pmid")) {

                    // add a related publication to this genotyping study
                    String pubMedId = value;
                    if (publicationMap.containsKey(pubMedId)) {
                        Item publication = publicationMap.get(pubMedId);
                        genotypingStudy.addToCollection("publications", publication);
                    } else {
                        // create a new publication
                        Item publication = createItem("Publication");
                        publication.setAttribute("pubMedId", pubMedId);
                        publicationMap.put(pubMedId, publication);
                        genotypingStudy.addToCollection("publications", publication);
                    }

		} else if (key.toLowerCase().equals("parent")) {

		    // add a parent strain
		    String strainName = value;
		    if (strainMap.containsKey(strainName)) {
			Item strain = strainMap.get(strainName);
			genotypingStudy.addToCollection("parents", strain);
		    } else {
			Item strain = createItem("Strain");
			strain.setAttribute("primaryIdentifier", strainName);
			strain.setReference("organism", organism);
			store(strain);
			strainMap.put(strainName, strain);
			genotypingStudy.addToCollection("parents", strain);
		    }
			
                } else if (key.toLowerCase().equals("lines")) {

                    // We've got columns = lines and rows = markers
                    rowsAreMarkers = true;
                    LOG.info("Rows are genetic markers.");
                    System.out.println("Rows are genetic markers.");
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
			    lines[i].setReference("organism", organism);
			    store(lines[i]);
			    lineMap.put(lineName, lines[i]);
			}
			genotypingStudy.addToCollection("lines", lines[i]);
                    }

                } else if (key.toLowerCase().equals("markers")) {

                    // We've got columns = markers and rows = lines (Flapjack style)
                    rowsAreLines = true;
                    LOG.info("Rows are genotyping lines.");
                    System.out.println("Rows are genotyping lines.");
                    int num = parts.length - 1;    // number of markers
                    markers = new Item[num];       // ss1234567, etc.
                    for (int i=0; i<num; i++) {
                        String markerName = parts[i+1];
                        if (markerMap.containsKey(markerName)) {
                            markers[i] = markerMap.get(markerName);
                        } else {
                            markers[i] = createItem("GeneticMarker");
                            markers[i].setAttribute("primaryIdentifier", markerName);
                            markers[i].setAttribute("type", markerType);
                            markers[i].setReference("organism", organism);
                            store(markers[i]);
                            markerMap.put(markerName, markers[i]);
                        }
			genotypingStudy.addToCollection("markers", markers[i]);
                    }

                } else {

                    // get the genotype values
                    int num = parts.length - 1;      // number of markers or lines
                    String[] values = new String[num];
                    for (int i=0; i<num; i++) {
                        values[i] = parts[i+1];
                    }
                    
                    if (rowsAreMarkers) {
                        
                        // create/retrieve marker
                        String markerName = parts[0];
                        Item marker;
                        if (markerMap.containsKey(markerName)) {
                            marker = markerMap.get(markerName);
                        } else {
                            marker = createItem("GeneticMarker");
                            marker.setAttribute("primaryIdentifier", markerName);
                            marker.setAttribute("type", markerType);
                            marker.setReference("organism", organism);
                            store(marker);
                            markerMap.put(markerName, marker);
                        }
			genotypingStudy.addToCollection("markers", marker);
                        // create and store this marker's genotypeValues; they should be same order/number as lines.
                        for (int i=0; i<lines.length; i++) {
                            Item genotypeValue = createItem("GenotypeValue");
                            genotypeValue.setAttribute("value", values[i]);
                            genotypeValue.setReference("line", lines[i]);
                            genotypeValue.setReference("marker", marker);
                            store(genotypeValue);
                        }
                        
                    } else if (rowsAreLines) {

                        // create/retrieve line
                        String lineName = parts[0];
                        Item line;
                        if (lineMap.containsKey(lineName)) {
                            line = lineMap.get(lineName);
                        } else {
                            line = createItem("GenotypingLine");
                            line.setAttribute("primaryIdentifier", lineName);
                            line.setReference("organism", organism);
                            store(line);
                            lineMap.put(lineName, line);
                        }
			genotypingStudy.addToCollection("lines", line);
                        // create and store this lines's genotypeValues; they should be same order/number as markers.
                        for (int i=0; i<markers.length; i++) {
                            Item genotypeValue = createItem("GenotypeValue");
                            genotypeValue.setAttribute("value", values[i]);
                            genotypeValue.setReference("line", line);
                            genotypeValue.setReference("marker", markers[i]);
                            store(genotypeValue);
                        }
                        
                    }
                    
                }

            }

        } // end while reading lines/markers

        // each file has a unique genotyping study so store it here
        store(genotypingStudy);
        
        // wrap up this file
        br.close();
    }

    /**
     * Store the items stored in global maps; may have the same publication in multiple genotyping studies.
     */
    @Override
    public void close() throws Exception {
        store(publicationMap.values());
    }

}
