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
import java.util.Map;

import org.apache.log4j.Logger;

import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.xml.full.Item;

/**
 * Store the QTL-ontology term annotations from tab-delimited files. The file rows contain the QTL name and ontology term.
 *
 * #QTL        TOTerm
 * Qpl.zaas-3  TO:0002626
 * Qpl.zaas-3  PO:1234566
 *
 * @author Sam Hokin, NCGR
 */
public class QTLOntologyFileConverter extends BioFileConverter {
	
    private static final Logger LOG = Logger.getLogger(QTLOntologyFileConverter.class);

    Map<String,Item> qtlMap = new HashMap<String,Item>();
    Map<String,Item> termMap = new HashMap<String,Item>();
    
    /**
     * Create a new QTLOntologyFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public QTLOntologyFileConverter(ItemWriter writer, Model model) {
        super(writer, model);
    }

    /**
     * {@inheritDoc}
     * We process the QTL-ontology relationships by reading in from a tab-delimited file.
     */
    @Override
    public void process(Reader reader) throws Exception {

        // don't process README files
        if (getCurrentFile().getName().contains("README")) return;

        LOG.info("Processing file "+getCurrentFile().getName()+"...");

        // -------------------------------------------------------------------------------------------------
        // Run through the file and add the associated ontology terms to the QTLs.
        // -------------------------------------------------------------------------------------------------

        // header constants
        String taxonId = null;

        BufferedReader buffReader = new BufferedReader(reader);
	String line;
        while ((line=buffReader.readLine())!=null) {

            if (!line.startsWith("#") && line.trim().length()>0) {

                String[] parts = line.split("\t");
                if (parts.length==2) {
                    
                    String key = parts[0];
                    String value = parts[1];
                    
                    String qtlName = key;
                    String termID = value;
                    String termType = termID.substring(0,2); // GO, PO, TO, etc.
                    
                    // find the QTL in the map, or add it with qtlName=primaryIdentifier
                    Item qtl = null;
                    if (qtlMap.containsKey(qtlName)) {
                        qtl = qtlMap.get(qtlName);
                    } else {
                        qtl = createItem("QTL");
                        qtl.setAttribute("primaryIdentifier", qtlName);
                        qtlMap.put(qtlName, qtl);
                    }
                    
                    // find the term in the map, or add it with termID=identifier
                    Item term = null;
                    if (termMap.containsKey(termID)) {
                        term = termMap.get(termID);
                    } else {
                        // determine the type from the prefix, hopefully supported in the data model!
                        term = createItem(termType+"Term");
                        term.setAttribute("identifier", termID);
                        termMap.put(termID, term);
                    }
                    
                    // create this annotation, associate it with the term and QTL, and store it
                    Item annotation = createItem(termType+"Annotation");
                    annotation.setReference("ontologyTerm", term);
                    annotation.setReference("subject", qtl);
                    store(annotation);
                    LOG.info("Storing annotation for QTL "+qtlName+" and term "+termID);
                    
                }
		
            }
	    
        }
	
        buffReader.close();

    }

    /**
     * Store the items we've collected in maps and sets
     */
    @Override
    public void close() throws Exception {
        
        LOG.info("Storing "+qtlMap.size()+" QTL items...");
        store(qtlMap.values());
	
        LOG.info("Storing "+termMap.size()+" OntologyTerm items...");
        store(termMap.values());
        
    }
    
}
