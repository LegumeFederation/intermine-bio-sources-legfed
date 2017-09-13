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
 * #QTL        TOTerms
 * Qpl.zaas-3  TO:0002626,TO:1234567
 * Qpl.zaas-3  TO:0002626
 * </pre>
 *
 * @author Sam Hokin, NCGR
 */
public class QTLTOFileConverter extends BioFileConverter {
	
    private static final Logger LOG = Logger.getLogger(QTLTOFileConverter.class);

    Map<String,Item> qtlMap = new HashMap<String,Item>();
    Map<String,Item> toTermMap = new HashMap<String,Item>();
    Map<String,Set<String>> toAnnotationMap = new HashMap<String,Set<String>>(); // store QTL IDs per TO ID to avoid dupes
    
    /**
     * Create a new QTLTOFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public QTLTOFileConverter(ItemWriter writer, Model model) {
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
        // Run through the QTL-TOs file and add the associated markers to the given QTLs
        // NOTE: given names are used as _primary_ identifiers
        // -------------------------------------------------------------------------------------------------

        BufferedReader qtlTOReader = new BufferedReader(reader);
	String line;
        while ((line=qtlTOReader.readLine())!=null) {

	    if (line.startsWith("#")) {

		// comment, do nothing
		
            } else {

                QTLTORecord rec = new QTLTORecord(line);
                if (rec.qtlName!=null && rec.toTerms!=null && rec.toTerms.length>0) {
                    
                    // find the QTL in the map, or add it with qtlName=primaryIdentifier
                    Item qtl = null;
                    if (qtlMap.containsKey(rec.qtlName)) {
                        qtl = qtlMap.get(rec.qtlName);
                    } else {
                        qtl = createItem("QTL");
                        qtl.setAttribute("primaryIdentifier", rec.qtlName);
                        qtlMap.put(rec.qtlName, qtl);
                    }

		    // TO terms
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
                                // add this annotation to the QTL
                                Item toAnnotation = createItem("TOAnnotation");
                                toAnnotation.setReference("ontologyTerm", toTerm);
                                toAnnotation.setReference("subject", qtl);
                                store(toAnnotation);
                                qtl.addToCollection("ontologyAnnotations", toAnnotation);
                                // add the QTL to the set of ones with this annotation stored
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
            
        qtlTOReader.close();

    }

    /**
     * Store the items we've collected in maps and sets
     */
    @Override
    public void close() throws Exception {

        LOG.info("Storing "+qtlMap.size()+" QTL items...");
        store(qtlMap.values());

	LOG.info("Storing "+toTermMap.size()+" TOTerm items...");
	store(toTermMap.values());

    }
    
}
