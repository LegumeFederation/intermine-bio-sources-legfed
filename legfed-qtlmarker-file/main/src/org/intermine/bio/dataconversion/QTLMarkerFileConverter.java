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
import org.intermine.xml.full.Item;

/**
 * Store the QTL-marker data and relationships from tab-delimited files.
 * The files simply contain the marker name, QTL name, optional trait and optional comma-separated TO terms.
 * Preceed with mapping populations and publications if desired.
 *
 * <pre>
 * PMID     27658053
 * PMID     12345467
 * MappingPopulation     ZN016_x_Zhijiang282
 * MappingPopulation      CB27_x_Sanzi7
 * #Marker  QTL         Trait       TOTerms
 * 2_04960  Qpl.zaas-3  Pod length  TO:0002626,TO:1234567
 * 2_17765  Qpl.zaas-3  Pod length  TO:0002626
 * </pre>
 *
 * @author Sam Hokin, NCGR
 */
public class QTLMarkerFileConverter extends BioFileConverter {
	
    private static final Logger LOG = Logger.getLogger(QTLMarkerFileConverter.class);

    Map<String,Item> mappingPopulationMap = new HashMap<String,Item>();
    Map<String,Item> publicationMap = new HashMap<String,Item>();
    Map<String,Item> authorMap = new HashMap<String,Item>();

    Map<String,Item> markerMap = new HashMap<String,Item>();
    Map<String,Item> qtlMap = new HashMap<String,Item>();
    Map<String,Item> toTermMap = new HashMap<String,Item>();
    Map<String,Set<String>> toAnnotationMap = new HashMap<String,Set<String>>(); // store QTL IDs per TO ID to avoid dupes
    
    /**
     * Create a new QTLMarkerFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public QTLMarkerFileConverter(ItemWriter writer, Model model) {
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

        // -------------------------------------------------------------------------------------------------
        // Run through the QTL-Markers file and add the associated markers to the given QTLs
        // NOTE: given names are used as _primary_ identifiers
        // -------------------------------------------------------------------------------------------------

        Set<Item> qtlPubs = new HashSet<Item>(); // pubs associated with THIS set of QTLs
	Set<Item> qtlMPs = new HashSet<Item>();  // mapping populations associated with THIS set of QTLs
        
        BufferedReader qtlMarkerReader = new BufferedReader(reader);
	String line;
        while ((line=qtlMarkerReader.readLine())!=null) {

	    if (line.startsWith("#")) {

		// comment, do nothing
		
	    } else if (line.startsWith("PMID")) {
                
                // load a publication
                String[] parts = line.split("\t");
                String pubMedId = parts[1];
                if (publicationMap.containsKey(pubMedId)) {
		    LOG.info("Retrieving publication from PMID:"+pubMedId);
                    qtlPubs.add(publicationMap.get(pubMedId));
                } else {
		    LOG.info("Creating publication from PMID:"+pubMedId);
                    Item publication = PublicationTools.getPublicationFromPMID(this, Integer.parseInt(pubMedId));
                    if (publication!=null) {
                        List<Item> authors = PublicationTools.getAuthorsFromPMID(this, Integer.parseInt(pubMedId));
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
                        qtlPubs.add(publication);
                    } else {
			LOG.error("Publication is null!");
		    }
                }

	    } else if (line.startsWith("MappingPopulation")) {

		// load a mapping population
		String[] parts = line.split("\t");
		String mpId = parts[1];
		if (mappingPopulationMap.containsKey(mpId)) {
		    LOG.info("Retrieving mapping population:"+mpId);
		    qtlMPs.add(mappingPopulationMap.get(mpId));
		} else {
		    LOG.info("Creating mapping population:"+mpId);
		    Item mappingPopulation = createItem("MappingPopulation");
		    mappingPopulation.setAttribute("primaryIdentifier", mpId);
		    mappingPopulationMap.put(mpId, mappingPopulation);
		    qtlMPs.add(mappingPopulation);
		}

            } else {

                QTLMarkerRecord rec = new QTLMarkerRecord(line);
                if (rec.qtlName!=null && rec.markerName!=null) {

                    // find the QTL in the map, or add it with qtlName=primaryIdentifier
                    Item qtl = null;
                    if (qtlMap.containsKey(rec.qtlName)) {
                        qtl = qtlMap.get(rec.qtlName);
                    } else {
                        qtl = createItem("QTL");
                        qtl.setAttribute("primaryIdentifier", rec.qtlName);
                        if (rec.trait!=null) qtl.setAttribute("secondaryIdentifier", rec.trait);
			for (Item pub : qtlPubs) qtl.addToCollection("publications", pub);
			for (Item mp : qtlMPs) qtl.addToCollection("mappingPopulations", mp);
                        qtlMap.put(rec.qtlName, qtl);
                    }
                    
                    // find the genetic marker in the map, or add it with markerName=primaryIdentifier
                    Item marker = null;
                    if (markerMap.containsKey(rec.markerName)) {
                        marker = markerMap.get(rec.markerName);
                    } else {
                        marker = createItem("GeneticMarker");
                        marker.setAttribute("primaryIdentifier", rec.markerName);
                        markerMap.put(rec.markerName, marker);
                    }

                    // add this marker to this QTL's collection
                    LOG.info("Adding QTL="+rec.qtlName+" Genetic Marker="+rec.markerName+" to associatedGeneticMarkers.");
                    qtl.addToCollection("associatedGeneticMarkers", marker);

		    // TO terms
		    if (rec.toTerms!=null && rec.toTerms.length>0) {
			for (int i=0; i<rec.toTerms.length; i++) {
			    String toTermId = rec.toTerms[i];
			    Item toTerm = null;
			    if (toTermMap.containsKey(toTermId)) {
				toTerm = toTermMap.get(toTermId);
			    } else {
				toTerm = createItem("TOTerm");
				toTerm.setAttribute("identifier", toTermId);
				toTermMap.put(toTermId, toTerm);
			    }
			    // avoid dupe annotations
			    if (toAnnotationMap.containsKey(toTermId)) {
				// see if this QTL is already in
				Set<String> qtlSet = toAnnotationMap.get(toTermId);
				if (qtlSet.contains(rec.qtlName)) {
				    // do nothing
				} else {
				    Item toAnnotation = createItem("TOAnnotation");
				    toAnnotation.setReference("ontologyTerm", toTerm);
				    toAnnotation.setReference("subject", qtl);
				    store(toAnnotation);
				    qtl.addToCollection("ontologyAnnotations", toAnnotation);
				    qtlSet.add(rec.qtlName);
				}
			    } else {
				// brand new TO term in an annotation
				Item toAnnotation = createItem("TOAnnotation");
				toAnnotation.setReference("ontologyTerm", toTerm);
				toAnnotation.setReference("subject", qtl);
				store(toAnnotation);
				qtl.addToCollection("ontologyAnnotations", toAnnotation);
				Set<String> qtlSet = new HashSet<String>();
				qtlSet.add(rec.qtlName);
				toAnnotationMap.put(toTermId, qtlSet);
			    }
			}
		    }
                }

            }

        }
            
        qtlMarkerReader.close();

    }

    /**
     * Store the items we've collected from the files
     */
    @Override
    public void close() throws Exception {

	LOG.info("Storing "+mappingPopulationMap.size()+" MappingPopulation items...");
	store(mappingPopulationMap.values());

        LOG.info("Storing "+publicationMap.size()+" Publication items...");
        store(publicationMap.values());

        LOG.info("Storing "+authorMap.size()+" Author items...");
        store(authorMap.values());

        LOG.info("Storing "+markerMap.size()+" GeneticMarker items...");
        store(markerMap.values());

        LOG.info("Storing "+qtlMap.size()+" QTL items...");
        store(qtlMap.values());

	LOG.info("Storing "+toTermMap.size()+" TOTerm items...");
	store(toTermMap.values());

    }
    
}
