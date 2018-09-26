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

import java.util.Map;
import java.util.HashMap;

import org.apache.log4j.Logger;

import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.xml.full.Item;

/**
 * Store GWAS QTL/marker data from a SoyBase (or other) dump. See GWASFileRecord for individual record format.
 *
 * TaxonID	         3847
 * Variety	         Williams82
 * primaryIdentifier     KGK20170714.1
 * platformName          SoySNP50k
 * platformDetails       Illumina Infinium BeadChip
 * numberLociTested      52041
 * numberGermplasmTested 12116
 *
 * @author Sam Hokin, NCGR
 */
public class GWASFileConverter extends BioFileConverter {
	
    private static final Logger LOG = Logger.getLogger(GWASFileConverter.class);

    // store items in maps to avoid duplicates
    Map<String,Item> organismMap = new HashMap<String,Item>();
    Map<String,Item> chromosomeMap = new HashMap<String,Item>();
    Map<String,Item> markerMap = new HashMap<String,Item>();
    Map<String,Item> qtlMap = new HashMap<String,Item>();
    
    /**
     * Create a new GWASFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public GWASFileConverter(ItemWriter writer, Model model) {
        super(writer, model);
    }

    /**
     * {@inheritDoc}
     * Process the marker-chromosome relationships by reading in from a tab-delimited file.
     */
    @Override
    public void process(Reader reader) throws Exception {

        // don't process README files
        if (getCurrentFile().getName().contains("README")) return;

        LOG.info("Processing file "+getCurrentFile().getName()+"...");

        // to identify the organism
        String taxonId = null;
        String variety = null;
        Item organism = null;

        // to identify the experiment
        String primaryIdentifier = null;
        String platformName = null;
        String platformDetails = null;
        int numberLociTested = 0;
        int numberGermplasmTested = 0;
        Item gwasExperiment = null;

        BufferedReader bufferedReader = new BufferedReader(reader);
	String line;
        while ((line=bufferedReader.readLine())!=null) {

            if (!line.startsWith("#") && line.trim().length()>0) {

                String[] parts = line.split("\t");
                String lcLine = line.toLowerCase();
                
                // create the organism if we have the required stuff (and it's not already stored)
                if (organism==null && taxonId!=null && variety!=null) {
                    String organismKey = taxonId+"_"+variety;
                    if (organismMap.containsKey(organismKey)) {
                        organism = organismMap.get(organismKey);
                    } else {
                        organism = createItem("Organism");
                        organism.setAttribute("taxonId", taxonId);
                        organism.setAttribute("variety", variety);
                        store(organism);
                        LOG.info("Created and stored Organism with taxonId:variety = "+taxonId+":"+variety);
                        organismMap.put(organismKey, organism);
                    }
                }

                // create the GWAS experiment if we have the required stuff
                if (gwasExperiment==null && primaryIdentifier!=null && numberLociTested>0 && numberGermplasmTested>0) {
                    gwasExperiment = createItem("GWASExperiment");
                    gwasExperiment.setAttribute("primaryIdentifier", primaryIdentifier);
                    gwasExperiment.setAttribute("numberLociTested", String.valueOf(numberLociTested));
                    gwasExperiment.setAttribute("numberGermplasmTested", String.valueOf(numberGermplasmTested));
                    if (platformName!=null) gwasExperiment.setAttribute("platformName", platformName);
                    if (platformDetails!=null) gwasExperiment.setAttribute("platformDetails", platformDetails);
                    store(gwasExperiment);
                    LOG.info("Created and stored GWASExperiment with primaryIdentifier="+primaryIdentifier);
                }

                // header items
                if (lcLine.startsWith("taxonid")) {
                    taxonId = parts[1];
                } else if (lcLine.startsWith("variety")) {
                    variety = parts[1];
                } else if (lcLine.startsWith("primaryidentifier")) {
                    primaryIdentifier = parts[1];
                } else if (lcLine.startsWith("platformname")) {
                    platformName = parts[1];
                } else if (lcLine.startsWith("platformdetails")) {
                    platformDetails = parts[1];
                } else if (lcLine.startsWith("numberlocitested")) {
                    numberLociTested = Integer.parseInt(parts[1]);
                } else if (lcLine.startsWith("numbergermplasmtested")) {
                    numberGermplasmTested = Integer.parseInt(parts[1]);
                } else {

                    // check that we've got an organism - fatal exit if not
                    if (organism==null) {
                        String errorMsg = "Organism has not been formed for GWAS record import in file "+getCurrentFile().getName();
                        LOG.error(errorMsg);
                        throw new RuntimeException(errorMsg);
                    }

                    // check that we've got a GWAS experiment - fatal exit if not
                    if (gwasExperiment==null) {
                        String errorMsg = "GWAS experiment is not fully described:\n" +
                            "primaryIdentifier="+primaryIdentifier+"\n" +
                            "platformName="+platformName+"\n" +
                            "platformDetails="+platformDetails+"\n" +
                            "numberLociTested="+numberLociTested+"\n" +
                            "numberGermplasmTested="+numberGermplasmTested;
                        LOG.error(errorMsg);
                        throw new RuntimeException(errorMsg);
                    }

                    // process a record
                    GWASFileRecord rec = new GWASFileRecord(line);

                    //
                    // Genetic Marker, Chromosome
                    //
                    Item marker = null;
                    String markerKey = rec.locusName;
                    if (markerMap.containsKey(markerKey)) {
                        marker = markerMap.get(markerKey);
                    } else {
                        marker = createItem("GeneticMarker");
                        marker.setReference("organism", organism);
                        marker.setAttribute("primaryIdentifier", rec.locusName);
                        marker.setAttribute("type", rec.type);
                        // set the chromosome or supercontig reference
                        // HACK: assume supercontig has "scaffold" or "contig" in the name
                        boolean isSupercontig = (rec.chromosome.toLowerCase().contains("scaffold") || rec.chromosome.toLowerCase().contains("contig"));
                        Item chromosome = null;
                        if (chromosomeMap.containsKey(rec.chromosome)) {
                            chromosome = chromosomeMap.get(rec.chromosome);
                        } else {
                            // create and store this chromosome/supercontig
                            if (isSupercontig) {
                                chromosome = createItem("Supercontig");
                            } else {
                                chromosome = createItem("Chromosome");
                            }
                            chromosome.setReference("organism", organism);
                            chromosome.setAttribute("primaryIdentifier", rec.chromosome);
                            chromosomeMap.put(rec.chromosome, chromosome);
                            store(chromosome);
                            LOG.info("Created and stored chromosome/supercontig: "+rec.chromosome);
                        }
                        if (isSupercontig) {
                            marker.setReference("supercontig", chromosome);
                        } else {
                            marker.setReference("chromosome", chromosome);
                        }
                        // create and store the Location object, use + strand since it's not defined
                        Item location = createItem("Location");
                        location.setReference("locatedOn", chromosome);
                        location.setReference("feature", marker);
                        location.setAttribute("start", String.valueOf(rec.start));
                        location.setAttribute("end", String.valueOf(rec.end));
                        location.setAttribute("strand", String.valueOf(+1));
                        store(location);
                        // set the chromosomeLocation/supercontigLocation reference and store the marker
                        if (isSupercontig) {
                            marker.setReference("supercontigLocation", location);
                        } else {
                            marker.setReference("chromosomeLocation", location);
                        }                    
                        markerMap.put(markerKey, marker);
                    }

                    //
                    // QTL
                    //
                    Item qtl = null;
                    String qtlKey = rec.gwasName;
                    if (qtlMap.containsKey(qtlKey)) {
                        qtl = qtlMap.get(qtlKey);
                        LOG.info("Duplicate QTL:"+qtlKey+" with marker:"+markerKey);
                    } else {
                        qtl = createItem("QTL");
                        // qtl.setReference("organism", organism);
                        qtl.setAttribute("primaryIdentifier", rec.gwasName);
                        qtl.setAttribute("traitName", rec.gwasClass);
                        qtl.setAttribute("pValue", String.valueOf(rec.pValue));
                        qtl.setReference("gwasExperiment", gwasExperiment);
                        qtlMap.put(qtlKey, qtl);
                        LOG.info("Added QTL:"+qtlKey+" with marker:"+markerKey);
                    }
                    
                    // relate the QTL and marker
                    marker.addToCollection("QTLs", qtl);

                }
            }
        }
        
        bufferedReader.close();

    }

    /**
     * Store markers and QTLs.
     */
    @Override
    public void close() throws Exception {
        for (Item qtl : qtlMap.values()) store(qtl);
        for (Item marker : markerMap.values()) store(marker);
    }
    
}
