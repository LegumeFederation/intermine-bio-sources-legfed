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

import java.util.List;
import java.util.ArrayList;
import java.util.Map;
import java.util.HashMap;

import org.apache.log4j.Logger;

import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Item;

/**
 * Store details on SNP array markers from a tab-delimited file. Data fields are:
 *
 * TaxonID	3917
 * ArrayName	Illumina Cowpea iSelect Consortium Array
 * PMID         27775877
 * MarkerType	SNP
 * #ID	   DesignSequence   Alleles Source BeadType StepDescription  AssociatedGene_1 AssociatedGene_2 ...
 * 1_0002  TAC..[A/G]..AGA  A/G     COPA   0	    List 1 inside    AT1G53750        Phvul.010G122200 ...
 *
 * Note: genes are matched by primary identifier, so be sure that existing genes in the mine have the same style of primaryIdentifier.
 *
 * Note: Markers do not get secondary identifiers here. They are merged on primaryIdentifier+organism.
 *
 * @author Sam Hokin
 */
public class SNPMarkerFileConverter extends BioFileConverter {
	
    private static final Logger LOG = Logger.getLogger(SNPMarkerFileConverter.class);

    // items saved at end are stored in maps
    Map<String,Item> publicationMap = new HashMap<String,Item>();       // keyed by pubMedId
    Map<String,Item> authorMap = new HashMap<String,Item>();            // keyed by name
    Map<String,Item> geneMap = new HashMap<String,Item>();              // keyed by primaryIdentifier
    Map<String,Item> organismMap = new HashMap<String,Item>();          // keyed by TaxonID
    Map<String,Item> strainMap = new HashMap<String,Item>();            // keyed by primaryIdentifier

    /**
     * Create a new SNPMarkerFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public SNPMarkerFileConverter(ItemWriter writer, Model model) {
        super(writer, model);
    }

    /**
     * {@inheritDoc}
     * Read in the SNP Marker file and store the linkage groups, QTLs and genetic markers along with their ranges and positions
     */
    @Override
    public void process(Reader reader) throws Exception {

        // don't process README files
        if (getCurrentFile().getName().contains("README")) return;

        LOG.info("Processing SNP Marker file "+getCurrentFile().getName()+"...");

        // these objects are created per file
        Item publication = null;
        String arrayName = null;
        String markerType = null;

        // this file's organism
        Item organism = null;
        String taxonId = null;

        // -----------------------------------------------------------------------------------------------------------------------------------
        // Load genetic markers and associated data
        // -----------------------------------------------------------------------------------------------------------------------------------

        BufferedReader br = new BufferedReader(reader);
        String line = null;
        while ((line=br.readLine())!=null) {

            // create these markers' organism if TaxonId has been supplied
            if (organism==null && taxonId!=null) {
                if (organismMap.containsKey(taxonId)) {
                    organism = organismMap.get(taxonId);
                } else {
                    // create and store this organism
                    organism = createItem("Organism");
                    organism.setAttribute("taxonId", taxonId);
                    store(organism);
                    organismMap.put(taxonId, organism);
                    LOG.info("Stored marker organism: "+taxonId);
                }
            }

	    if (line.startsWith("#")) {

		// do nothing, comment

            } else if (line.toLowerCase().startsWith("taxonid")) {

                // this file's organism.taxonId
                String[] parts = line.split("\t");
                taxonId = parts[1];

	    } else if (line.toLowerCase().startsWith("arrayname")) {
		
                // array name is stored in a string, a marker attribute
                String[] parts = line.split("\t");
                arrayName = parts[1];

            } else if (line.toLowerCase().startsWith("markertype")) {

                // marker type is stored in a string, a marker attribute
                String[] parts = line.split("\t");
                markerType = parts[1];

            } else if (line.toLowerCase().startsWith("pmid")) {

                // get the publication
                String[] parts = line.split("\t");
                String pubMedId = parts[1];
                if (publicationMap.containsKey(pubMedId)) {
                    publication = publicationMap.get(pubMedId);
                } else {
                    publication = createItem("Publication");
                    publication.setAttribute("pubMedId", pubMedId);
                    store(publication);
                    publicationMap.put(pubMedId, publication);
                }

            } else {

                // bail if we've not specified the markers' organism
                if (organism==null) {
                    LOG.error("Marker organism not specified: taxonId="+taxonId);
                    throw new RuntimeException("Marker organism not specified: taxonId="+taxonId);
                }
                    
                // looks like it's a data record
                SNPMarkerRecord record = new SNPMarkerRecord(line);

                // create and store general stuff
                Item marker = createItem("GeneticMarker");
                marker.setReference("organism", organism);
                if (markerType!=null) marker.setAttribute("type", markerType);
                if (arrayName!=null) marker.setAttribute("arrayName", arrayName);
                if (publication!=null) marker.setReference("publication", publication);

                // add data from the marker record
                marker.setAttribute("primaryIdentifier", record.marker);
                marker.setAttribute("designSequence", record.designSequence);
                if (record.alleles!=null) marker.setAttribute("alleles", record.alleles);
                if (record.source!=null) marker.setAttribute("source", record.source);
                if (record.beadType!=null) marker.setAttribute("beadType", String.valueOf(record.beadType));
                if (record.stepDescription!=null) marker.setAttribute("stepDescription", record.stepDescription);
                if (record.associatedGenes!=null) {
                    for (int i=0; i<record.associatedGenes.length; i++) {
                        if (record.associatedGenes[i]!=null) {
                            String primaryIdentifier = record.associatedGenes[i];
                            Item gene;
                            if (geneMap.containsKey(primaryIdentifier)) {
                                gene = geneMap.get(primaryIdentifier);
                            } else {
                                gene = createItem("Gene");
                                gene.setAttribute("primaryIdentifier", primaryIdentifier);
                                store(gene);
                                geneMap.put(record.associatedGenes[i], gene);
                            }
                            marker.addToCollection("associatedGenes", gene);
                        }
                    }
                }

                // store this marker, we're done with it
                store(marker);

            }
            
        }
        
        br.close();
        
    }

}
