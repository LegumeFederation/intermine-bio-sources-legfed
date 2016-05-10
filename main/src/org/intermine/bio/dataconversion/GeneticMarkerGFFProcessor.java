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
import java.util.List;
import java.util.Map;

import org.apache.log4j.Logger;

import org.intermine.bio.util.OrganismData;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Attribute;
import org.intermine.xml.full.Item;
import org.intermine.xml.full.Reference;
import org.intermine.xml.full.ReferenceList;

/**
 * Store the genetic marker and genomic data from a GFF file.
 *
 * Input GFF file is set in project.xml as <property name="src.data.file" location="/your/gff/file/path/here.gff"/>
 * The GFF file may contain duplicates - the last one read is the one that is loaded.
 *
 * @author Sam Hokin
 */
public class GeneticMarkerGFFProcessor extends LegfedFileProcessor {
	
    private static final Logger LOG = Logger.getLogger(GeneticMarkerGFFProcessor.class);

    // global objects used by populateGeneticMarker()
    Item organism;
    
    /**
     * Create a new GeneticMarkerGFFProcessor
     * @param legfedFileConverter the LegfedFileConverter that is controlling this processor
     */
    public GeneticMarkerGFFProcessor(LegfedFileConverter legfedFileConverter) {
        super(legfedFileConverter);
    }

    /**
     * {@inheritDoc}
     * We process the supplied GFF file to create GeneticMarker Items, along with Chromosome and Supercontig Items from the seqid column.
     */
    @Override
    public void process(Reader reader) throws ObjectStoreException {

        // ---------------------------------------------------------
        // INITIAL DATA LOADING
        // ---------------------------------------------------------

        // get the organism taxon ID; enforce only one
        OrganismData[] organisms = getLegfedFileConverter().getOrganismsToProcess().toArray(new OrganismData[0]);
        if (organisms.length>1) {
            String error = "Multiple organisms specified in data source; GeneticMarkerGFFProcessor can only process one organism at a time.";
            LOG.error(error);
            throw new RuntimeException(error);
        }
        int taxonId = organisms[0].getTaxonId();

        // create organism Item - global so can be used in populate routines
        organism = getLegfedFileConverter().createItem("Organism");
        BioStoreHook.setSOTerm(getLegfedFileConverter(), organism, "organism", getLegfedFileConverter().getSequenceOntologyRefId());
        organism.setAttribute("taxonId", String.valueOf(taxonId));
        LOG.info("Created and stored organism Item for taxonId="+taxonId+".");

        // -------------------------------------------------------------------------------------------------------------------
        // Load the chromosomes, supercontigs and genetic markers from the GFF file
        // -------------------------------------------------------------------------------------------------------------------
        
        try {

            // store chromosomes in a map
            Map<String,Item> chromosomeMap = new HashMap<String,Item>();
            // and supercontigs
            Map<String,Item> supercontigMap = new HashMap<String,Item>();
            // and markers, since we may have duplicates - last one found is the one that sticks
            Map<String,Item> markerMap = new HashMap<String,Item>();
            // and marker locations, so we don't store those that belong to markers that get overwritten; also keyed by marker's name
            Map<String,Item> locationMap = new HashMap<String,Item>();
            
            BufferedReader gffReader = new BufferedReader(reader);
            String gffLine;
            LOG.info("Reading GFF file...");
            while ((gffLine = gffReader.readLine()) != null) {
                GFFRecord gff = new GFFRecord(gffLine);
                if (gff.seqid!=null && gff.attributeName!=null) {
                    Item marker = getLegfedFileConverter().createItem("GeneticMarker");
                    Item location = getLegfedFileConverter().createItem("Location");
                    // standard chromosome naming convention, e.g. "glyma.Chr01"
                    if (gff.seqid.contains("Chr")) {
                        // create and store this chromosome if not already in map; else get it from the map
                        String chrName = gff.seqid;
                        Item chromosome = null;
                        if (chromosomeMap.containsKey(chrName)) {
                            chromosome = chromosomeMap.get(chrName);
                        } else {
                            chromosome = getLegfedFileConverter().createItem("Chromosome");
                            BioStoreHook.setSOTerm(getLegfedFileConverter(), chromosome, "chromosome", getLegfedFileConverter().getSequenceOntologyRefId());
                            chromosome.setAttribute("primaryIdentifier", chrName);
                            chromosomeMap.put(chrName, chromosome);
                        }
                        // populate the chromosome location
                        location.setAttribute("start", String.valueOf(gff.start));
                        location.setAttribute("end", String.valueOf(gff.end));
                        location.setReference("feature", marker);
                        location.setReference("locatedOn", chromosome);
                        // associate the genetic marker with this chromosome/location
                        marker.setReference("chromosome", chromosome);
                        marker.setReference("chromosomeLocation", location);
                    } else {
                        // create and store this supercontig if not already in map; else get it from the map
                        String supercontigName = gff.seqid;
                        Item supercontig = null;
                        if (supercontigMap.containsKey(supercontigName)) {
                            supercontig = supercontigMap.get(supercontigName);
                        } else {
                            supercontig = getLegfedFileConverter().createItem("Supercontig");
                            BioStoreHook.setSOTerm(getLegfedFileConverter(), supercontig, "supercontig", getLegfedFileConverter().getSequenceOntologyRefId());
                            supercontig.setAttribute("primaryIdentifier", supercontigName);
                            supercontigMap.put(supercontigName, supercontig);
                        }
                        // populate the supercontig location
                        location.setAttribute("start", String.valueOf(gff.start));
                        location.setAttribute("end", String.valueOf(gff.end));
                        location.setReference("feature", marker);
                        location.setReference("locatedOn", supercontig);
                        // associate the genetic marker with this supercontig/location
                        marker.setReference("supercontig", supercontig);
                        marker.setReference("supercontigLocation", location);
                    }
                    // set other marker attributes
                    BioStoreHook.setSOTerm(getLegfedFileConverter(), marker, "genetic_marker", getLegfedFileConverter().getSequenceOntologyRefId());
                    marker.setReference("organism", organism);
                    marker.setAttribute("primaryIdentifier", gff.attributeName);
                    marker.setAttribute("type", gff.type);
                    marker.setAttribute("length", String.valueOf(gff.end-gff.start+1));
                    // add this marker and location to the maps; overwrites previous if same name
                    markerMap.put(gff.attributeName, marker);
                    locationMap.put(gff.attributeName, location);
                }
            }

            LOG.info("Created "+markerMap.size()+" distinct GeneticMarker items.");
            LOG.info("Created "+chromosomeMap.size()+" Chromosome items.");
            LOG.info("Created "+supercontigMap.size()+" Supercontig items.");

            // now store all the items
            store(organism);
            for (Item item : chromosomeMap.values()) store(item);
            for (Item item : supercontigMap.values()) store(item);
            for (Item item : locationMap.values()) store(item);
            for (Item item : markerMap.values()) store(item);

        } catch (Exception ex) {
            
            LOG.error(ex.getMessage());

        }

    }
    
}
