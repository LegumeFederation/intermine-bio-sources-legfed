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
import java.util.List;
import java.util.ArrayList;

import org.apache.log4j.Logger;

import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.xml.full.Item;

/**
 * Store QTL data and relationships from tab-delimited files organized by publication.
 * Mapping population name may be given, e.g. Sanzi_x_Vita7.
 * TaxonID and Variety are provided to establish uniqueness for merging, since the same QTL names often appear for different species/varieties.
 *
 * TaxonID     3817
 * Variety     IT97K-499-35
 * PMID        123456 [optional]
 * Journal     Intl. J. Genetics [optional]
 * Year        2007 [optional]
 * Volume      23 [optional]
 * Pages       345-357 [optional]
 * Title       A cowpea genetic map for the ages [optional]
 * DOI         10.1101/202044 [optional]
 * MappingPopulation Sanzi_x_Vita7 [optional]
 * Parent      Sanzi [optional]
 * Parent      Vita7 [optional]
 * Description A description of the purpose of this mapping population and why the parents were chosen. [optional]
 * #QTLName          TraitName
 * Seed weight 3-2   Seed weight
 *
 * @author Sam Hokin, NCGR
 */
public class QTLFileConverter extends BioFileConverter {
	
    private static final Logger LOG = Logger.getLogger(QTLFileConverter.class);

    // Typically one organism across all files
    Map<String,Item> organismMap = new HashMap<String,Item>();
    
    // Store QTLs in a map since there may be multiple pubs with the same QTLs
    Map<String,Item> qtlMap = new HashMap<String,Item>();

    // Several publications can use the same mapping population
    Map<String,Item> mappingPopulationMap = new HashMap<String,Item>();

    /**
     * Create a new QTLFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public QTLFileConverter(ItemWriter writer, Model model) {
        super(writer, model);
    }

    /**
     * {@inheritDoc}
     * We process the QTL-marker relationships by reading in from a tab-delimited file.
     */
    @Override
    public void process(Reader reader) throws Exception {

        // don't process README files
        if (getCurrentFile().getName().contains("README")) return;

        LOG.info("Processing file "+getCurrentFile().getName()+"...");

        // header values
        int taxonId = 0;
        String variety = null;

        // these are loaded to populate a publication by title, etc.
        String title = null;
        String journal = null;
        int year = 0;
        String volume = null;
        String pages = null;

        // the Items resulting from the header values
        Item organism = null;
        Item publication = null;
        List<Item> mappingPopulations = new ArrayList<Item>();
        
        BufferedReader qtlReader = new BufferedReader(reader);
	String line = null;
        while ((line=qtlReader.readLine())!=null) {

            // create and store organism (or pull it from map) if ready
            if (organism==null && taxonId>0 && variety!=null) {
                String organismKey = taxonId+"_"+variety;
                if (organismMap.containsKey(organismKey)) {
                    organism = organismMap.get(organismKey);
                } else {
                    organism = createItem("Organism");
                    organism.setAttribute("taxonId", String.valueOf(taxonId));
                    organism.setAttribute("variety", variety);
                    organismMap.put(organismKey, organism);
                    LOG.info("Storing organism: "+taxonId+" ("+variety+")");
                }
            }

            if (line.trim().length()>0 && !line.startsWith("#")) {

                String[] parts = line.split("\t");
		if (parts.length==2) {

                    String key = parts[0].trim();
                    String value = parts[1].trim();
		    
		    // header data
		    if (key.toLowerCase().equals("taxonid")) {
			taxonId = Integer.parseInt(value);
                        
		    } else if (key.toLowerCase().equals("variety")) {
			variety = value;
                        
		    } else if (key.toLowerCase().equals("pmid")) {
			int pmid = Integer.parseInt(value);
                        publication = createItem("Publication");
                        publication.setAttribute("pubMedId", String.valueOf(pmid));
                        store(publication);
                        LOG.info("Stored publication PMID="+pmid);

                    } else if (key.toLowerCase().equals("doi")) {
                        String doi = value;
                        publication = createItem("Publication");
                        publication.setAttribute("doi", doi);
                        store(publication);
                        LOG.info("Stored publication DOI="+doi);
                        
		    } else if (key.toLowerCase().equals("title")) {
			title = value;
                        
		    } else if (key.toLowerCase().equals("journal")) {
			journal = value;
                        
		    } else if (key.toLowerCase().equals("year")) {
			year = Integer.parseInt(value);
                        
		    } else if (key.toLowerCase().equals("volume")) {
			volume = value;
                        
		    } else if (key.toLowerCase().equals("pages")) {
			pages = value;
                        
		    } else if (key.toLowerCase().equals("mappingpopulation")) {

                        String mappingPopulationID = value;
                        if (mappingPopulationMap.containsKey(mappingPopulationID)) {
                            Item mappingPopulation = mappingPopulationMap.get(mappingPopulationID);
                            mappingPopulations.add(mappingPopulation);
                        } else {
                            Item mappingPopulation = createItem("MappingPopulation");
                            mappingPopulation.setAttribute("primaryIdentifier", mappingPopulationID);
                            mappingPopulationMap.put(mappingPopulationID, mappingPopulation);
                            mappingPopulations.add(mappingPopulation);
                            LOG.info("Storing mapping population:"+mappingPopulationID);
                        }

		    } else if (key.toLowerCase().equals("parent")) {

                        // be sure to put parents AFTER mapping population line!
			// if parent is present, then we assume that we only have one mapping population
			String parentVariety = value;
			String parentKey = taxonId+"_"+parentVariety;
			Item parent;
			if (organismMap.containsKey(parentKey)) {
			    parent = organismMap.get(parentKey);
			} else {
			    parent = createItem("Organism");
			    parent.setAttribute("taxonId", String.valueOf(taxonId));
			    parent.setAttribute("variety", parentVariety);
			    organismMap.put(parentKey, parent);
			    LOG.info("Storing parent: "+parentVariety);
			}
			for (Item mappingPopulation : mappingPopulations) {
			    mappingPopulation.addToCollection("parents", parent);
			}
                        LOG.info("Added parent "+parentVariety+" to mapping populations.");

                    } else if (key.toLowerCase().equals("description")) {

                        String description = value;
			for (Item mappingPopulation : mappingPopulations) {
			    mappingPopulation.setAttribute("description", description);
			}
                        LOG.info("Set description on mapping populations:"+description);
                        
		    } else {

                        // create and store the publication if we've collected at least a title
                        if (publication==null && title!=null && title.length()>0) {
                            publication = createItem("Publication");
                            publication.setAttribute("title", title);
                            if (journal!=null && journal.length()>0) publication.setAttribute("journal", journal);
                            if (year!=0) publication.setAttribute("year", String.valueOf(year));
                            if (volume!=null && volume.length()>0) publication.setAttribute("volume", volume);
                            if (pages!=null && pages.length()>0) publication.setAttribute("pages", pages);
                            store(publication);
                            LOG.info("Storing publication: "+title);
                        }

                        // associate this publication with the mapping populations
                        if (publication!=null) {
                            for (Item mappingPopulation : mappingPopulations) {
                                mappingPopulation.addToCollection("publications", publication);
                            }
                        }

			// load a QTL and trait
			String qtlName = key;
			String traitName = value;

                        // DEBUG
                        LOG.info("Storing QTL "+qtlName+":"+traitName);
                        
			// find the QTL in the map, or add it with qtlName=primaryIdentifier
			Item qtl;
			if (qtlMap.containsKey(qtlName)) {
			    qtl = qtlMap.get(qtlName);
			} else {
			    qtl = createItem("QTL");
			    qtl.setAttribute("primaryIdentifier", qtlName);
			    // qtl.setReference("organism", organism);
			    if (traitName.length()>0) qtl.setAttribute("traitName", traitName);
			    qtlMap.put(qtlName, qtl);
			}
			
			// add pub and mapping pops to this QTL's collections
			if (publication!=null) qtl.addToCollection("publications", publication);
                        for (Item mappingPopulation : mappingPopulations) {
                            qtl.addToCollection("mappingPopulations", mappingPopulation);
                        }

		    }
                }
            }
        }
        
        qtlReader.close();

    }


    /**
     * Store the items we've collected in maps across all the files
     */
    @Override
    public void close() throws Exception {
        LOG.info("Storing "+organismMap.size()+" Organism items...");
        store(organismMap.values());
        LOG.info("Storing "+mappingPopulationMap.size()+" MappingPopulation items...");
        store(mappingPopulationMap.values());
        LOG.info("Storing "+qtlMap.size()+" QTL items...");
        store(qtlMap.values());
    }
    
}
