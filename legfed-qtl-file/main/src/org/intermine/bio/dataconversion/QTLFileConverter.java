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

import org.apache.log4j.Logger;

import org.ncgr.intermine.PubMedPublication;
import org.ncgr.pubmed.PubMedSummary;

import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.xml.full.Item;

/**
 * Store QTL data and relationships from tab-delimited files organized by publication.
 * Parent names are munged into a mapping population name, e.g. Sanzi_x_Vita7 for parents Sanzi and Vita7.
 *
 * TaxonID     3817
 * Variety     IT97K-499-35
 * Parent_1    Sanzi [optional]
 * Parent_2    Vita7 [optional]
 * PMID        123456 [optional]
 * Journal     Intl. J. Genetics [optional]
 * Year        2007 [optional]
 * Volume      23 [optional]
 * Pages       345-357 [optional]
 * Title       A cowpea genetic map for the ages [optional]
 * 
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

    // Could have same author on multiple publications
    Map<String,Item> authorMap = new HashMap<String,Item>();

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
        if (getCurrentFile().getName().equals("README")) return;

        LOG.info("Processing file "+getCurrentFile().getName()+"...");

        // header values
        int taxonId = 0;
        String variety = null;
        String parent1 = null;
        String parent2 = null;
        int pmid = 0;
        String journal = null;
        int year = 0;
        String volume = null;
        String pages = null;
        String title = null;

        // the Items resulting from the header values
        Item organism = null;
        Item mappingPopulation = null;
        Item publication = null;

        BufferedReader qtlReader = new BufferedReader(reader);
	String line = null;
        while ((line=qtlReader.readLine())!=null) {
            if (line.trim().length()>0 && !line.startsWith("#")) {

                String[] parts = line.split("\t");
		if (parts.length==2) {
		    
		    // header data
		    if (parts[0].equals("TaxonID")) {
			taxonId = Integer.parseInt(parts[1].trim());
		    } else if (parts[0].equals("Variety")) {
			variety = parts[1].trim();
		    } else if (parts[0].equals("Parent_1")) {
			parent1 = parts[1].trim();
		    } else if (parts[0].equals("Parent_2")) {
			parent2 = parts[1].trim();
		    } else if (parts[0].equals("PMID")) {
			pmid = Integer.parseInt(parts[1].trim());
		    } else if (parts[0].equals("Journal")) {
			journal = parts[1].trim();
		    } else if (parts[0].equals("Year")) {
			year = Integer.parseInt(parts[1].trim());
		    } else if (parts[0].equals("Volume")) {
			volume = parts[1].trim();
		    } else if (parts[0].equals("Pages")) {
			pages = parts[1].trim();
		    } else if (parts[0].equals("Title")) {
			title = parts[1].trim();
		    } else {
			
			// create and store organism (or pull it from map)
			if (organism==null && taxonId>0 && variety!=null) {
			    String organismKey = taxonId+"_"+variety;
			    if (organismMap.containsKey(organismKey)) {
				organism = organismMap.get(organismKey);
			    } else {
				organism = createItem("Organism");
				organism.setAttribute("taxonId", String.valueOf(taxonId));
				organism.setAttribute("variety", variety);
				store(organism);
				organismMap.put(organismKey, organism);
				LOG.info("Stored organism: "+taxonId+" ("+variety+")");
			    }
			}
			
			// publication, if present, either a PMID or info columns
			if (publication==null && (pmid!=0 || (title!=null && title.trim().length()>0))) {
			    LOG.info("Creating new publication: PMID="+pmid+", Title="+title);
			    // search for PubMed summary, on either PMID or title
			    PubMedSummary pms;
			    if (pmid!=0) {
				pms = new PubMedSummary(pmid);
			    } else {
				pms = new PubMedSummary(title);
			    }
			    if (pms.id>0) {
				// it's in PubMed, so use that information
				PubMedPublication pubMedPub = new PubMedPublication(this, pms);
				publication = pubMedPub.getPublication();
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
				}
			    } else if (title!=null && title.trim().length()>0) {
				// store at least the publication title, no authors
				publication = createItem("Publication");
				publication.setAttribute("title", title.trim());
				if (journal!=null && journal.trim().length()>0) publication.setAttribute("journal", journal.trim());
				if (year!=0) publication.setAttribute("year", String.valueOf(year));
				if (volume!=null && volume.trim().length()>0) publication.setAttribute("volume", volume.trim());
				if (pages!=null && pages.trim().length()>0) publication.setAttribute("pages", pages.trim());
			    }
			    store(publication);
			    LOG.info("Stored publication: "+pmid+" "+title);
			}

                        // create and store mapping population (or pull it from map)
                        if (mappingPopulation==null && parent1!=null && parent2!=null) {
			    String mappingPopulationID = parent1+"_x_"+parent2;
			    if (mappingPopulationMap.containsKey(mappingPopulationID)) {
				mappingPopulation = mappingPopulationMap.get(mappingPopulationID);
			    } else {
				mappingPopulation = createItem("MappingPopulation");
				mappingPopulation.setAttribute("primaryIdentifier", mappingPopulationID);
				store(mappingPopulation);
				mappingPopulationMap.put(mappingPopulationID, mappingPopulation);
				LOG.info("Stored mapping population:"+mappingPopulationID);
			    }
                            if (publication!=null) mappingPopulation.addToCollection("publications", publication);
			}

			// go no further if we do not have an organism
			if (organism==null) {
			    LOG.error("No organism created for QTL references.");
			    throw new RuntimeException("No organism created for QTL references.");
			}
			
			// a QTL and trait
			String qtlName = parts[0];
			String traitName = parts[1];
			
			// find the QTL in the map, or add it with qtlName=primaryIdentifier
			Item qtl;
			if (qtlMap.containsKey(qtlName)) {
			    qtl = qtlMap.get(qtlName);
			} else {
			    qtl = createItem("QTL");
			    qtl.setAttribute("primaryIdentifier", qtlName);
			    qtl.setReference("organism", organism);
			    if (traitName.trim().length()>0) qtl.setAttribute("traitName", traitName);
			    qtlMap.put(qtlName, qtl);
			}
			
			// add to collections
			if (mappingPopulation!=null) qtl.addToCollection("mappingPopulations", mappingPopulation);
			if (publication!=null) qtl.addToCollection("publications", publication);

		    }
                }
            }
        }
        
        qtlReader.close();

    }

    /**
     * Store the items we've collected in maps and sets across all the files
     */
    @Override
    public void close() throws Exception {
        LOG.info("Storing "+qtlMap.size()+" QTL items...");
        store(qtlMap.values());
        store(authorMap.values());
    }
    
}
