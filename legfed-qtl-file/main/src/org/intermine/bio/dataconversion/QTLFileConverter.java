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
 * Store QTL data and relationships from tab-delimited files.
 *
 * QTLName Parent_1 Parent_2 TraitName [PMID] [Journal Year Volume Page Title]
 *
 * Parents are munged into a mapping population name, e.g. Sanzi_x_Vita7 for parents Sanzi and Vita7.
 *
 * @author Sam Hokin, NCGR
 */
public class QTLFileConverter extends BioFileConverter {
	
    private static final Logger LOG = Logger.getLogger(QTLFileConverter.class);

    // Store QTLs in a map since there may be multiple lines (different pubs)
    Map<String,Item> qtlMap = new HashMap<String,Item>();

    // Store mapping populations in a map, presumably far fewer than QTLs
    Map<String,Item> mappingPopulationMap = new HashMap<String,Item>();

    // Store publications in a map
    Map<String,Item> publicationMap = new HashMap<String,Item>();
    Map<String,Item> authorMap = new HashMap<String,Item>();
    
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

        BufferedReader qtlReader = new BufferedReader(reader);
	String line;
        while ((line=qtlReader.readLine())!=null) {

            if (!line.startsWith("#")) {

                QTLRecord rec = new QTLRecord(line);
                if (rec.qtlName!=null) {
                    
                    // find the QTL in the map, or add it with qtlName=primaryIdentifier
                    Item qtl = null;
                    if (qtlMap.containsKey(rec.qtlName)) {
                        qtl = qtlMap.get(rec.qtlName);
                    } else {
                        qtl = createItem("QTL");
                        qtl.setAttribute("primaryIdentifier", rec.qtlName);
                        qtl.setAttribute("traitName", rec.traitName);
                        qtlMap.put(rec.qtlName, qtl);
                    }
                    
                    // mapping population, if parents present
                    if (rec.parent1.trim().length()>0 && rec.parent2.trim().length()>0) {
                        String mappingPopulationID = rec.parent1+"_x_"+rec.parent2;
                        Item mappingPopulation = null;
                        if (mappingPopulationMap.containsKey(mappingPopulationID)) {
                            mappingPopulation = mappingPopulationMap.get(mappingPopulationID);
                        } else {
                            mappingPopulation = createItem("MappingPopulation");
                            mappingPopulation.setAttribute("primaryIdentifier", mappingPopulationID);
                            mappingPopulationMap.put(mappingPopulationID, mappingPopulation);
                        }
                        qtl.addToCollection("mappingPopulations", mappingPopulation);
                    }

                    // publication, if present, either a PMID or info columns
                    if (rec.pubPMID!=0 || (rec.pubTitle!=null && rec.pubTitle.trim().length()>0)) {
                        Item publication;
                        if (rec.pubPMID!=0 && publicationMap.containsKey(String.valueOf(rec.pubPMID))) {
                            // use PMID as key if PMID is present
                            publication = publicationMap.get(String.valueOf(rec.pubPMID));
                        } else if (rec.pubTitle!=null && rec.pubTitle.trim().length()>0 && publicationMap.containsKey(rec.pubTitle)) {
                            // use title as key if PMID is absent
                            publication = publicationMap.get(rec.pubTitle);
                        } else {
                            LOG.info("Creating new publication: PMID="+rec.pubPMID+", Title="+rec.pubTitle);
                            // new publication
                            PubMedSummary pms;
                            if (rec.pubPMID!=0) {
                                // search on PMID
                                pms = new PubMedSummary(rec.pubPMID);
                            } else {
                                // search on title
                                pms = new PubMedSummary(rec.pubTitle);
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
                                    if (rec.pubPMID!=0) {
                                        // key on PMID if PMID was present in record
                                        publicationMap.put(String.valueOf(pms.id), publication);
                                    } else {
                                        // key on title if PMID was not present in record
                                        publicationMap.put(rec.pubTitle, publication);
                                    }
                                    // add to QTL collection
                                    qtl.addToCollection("publications", publication);
                                }
                            } else if (rec.pubTitle!=null && rec.pubTitle.trim().length()>0) {
                                // store at least the publication title, no authors
                                publication = createItem("Publication");
                                publication.setAttribute("title", rec.pubTitle.trim());
                                if (rec.pubJournal!=null && rec.pubJournal.trim().length()>0) publication.setAttribute("journal", rec.pubJournal.trim());
                                if (rec.pubYear!=null && rec.pubYear.trim().length()>0) publication.setAttribute("year", rec.pubYear.trim());
                                if (rec.pubVolume!=null && rec.pubVolume.trim().length()>0) publication.setAttribute("volume", rec.pubVolume.trim());
                                if (rec.pubPages!=null && rec.pubPages.trim().length()>0) publication.setAttribute("pages", rec.pubPages.trim());
                                // key on title
                                publicationMap.put(rec.pubTitle, publication);
                                // add to QTL collection
                                qtl.addToCollection("publications", publication);
                            }
                        }
                    }
                    
                }
                
            }
        }
        
        qtlReader.close();

    }

    /**
     * Store the items we've collected in maps and sets
     */
    @Override
    public void close() throws Exception {

        LOG.info("Storing "+qtlMap.size()+" QTL items...");
        store(qtlMap.values());

        LOG.info("Storing "+mappingPopulationMap.size()+" MappingPopulation items...");
        store(mappingPopulationMap.values());

        LOG.info("Storing "+publicationMap.size()+" Publication items...");
        store(publicationMap.values());

        LOG.info("Storing "+authorMap.size()+" Author items...");
        store(authorMap.values());

    }
    
}
