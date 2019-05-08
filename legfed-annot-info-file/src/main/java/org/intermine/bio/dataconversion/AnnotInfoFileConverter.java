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
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Item;

/**
 * Load the information from an LIS data store info_annot file. Format is:
 *
 * pacId locusName transcriptName peptideName Pfam Panther KOG ec KO GO Best-hit-arabi-name arabi-symbol arabi-defline
 * 37170591 Phvul.001G000400 Phvul.001G000400.1 Phvul.001G000400.1.p PF00504 PTHR21649,PTHR21649:SF24 1.10.3.9 K14172 GO:0016020,GO:0009765 AT1G76570.1	Chlorophyll A-B binding family protein
 *
 * Presumption: one file per strain. (Could be several strains for an organism.)
 * NOTE: this will also link genes to proteins.
 *
 * @author Sam Hokin
 */
public class AnnotInfoFileConverter extends BioFileConverter {
	
    private static final Logger LOG = Logger.getLogger(AnnotInfoFileConverter.class);

    // store the organisms in a taxonId-keyed map since there may be multiples
    Map<String,Item> organisms = new HashMap<>();

    // store the genes and proteins in maps for the usual reasons
    Map<String,Item> genes = new HashMap<>();
    Map<String,Item> proteins = new HashMap<>();

    // a handy map to get taxon Ids from the leading piece of the file name (e.g. phavu)
    private static Map<String,String> lisOrgTaxonIds = new HashMap<String,String>() {{
            put("aradu","130453");
            put("araip","130454");
            put("vigra","157791");
            put("lotja","34305");
            put("arahy","3818");
            put("cajan","3821");
            put("cicar","3827");
            put("glyma","3847");
            put("lupan","3871");
            put("medtr","3880");
            put("phavu","3885");
            put("vigan","3914");
            put("vigun","3920");
            put("tripr","57577 ");
        }};

    /**
     * Create a new AnnotInfoFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public AnnotInfoFileConverter(ItemWriter writer, Model model) {
        super(writer, model);
    }

    /**
     * {@inheritDoc}
     * Read in the AnnotInfo file and store the linkage groups, QTLs and genetic markers along with their ranges and positions
     */
    @Override
    public void process(Reader reader) throws Exception {

        // don't process README files
        if (getCurrentFile().getName().contains("README")) return;

        LOG.info("Processing AnnotInfo file "+getCurrentFile().getName()+"...");
        
        // create and add the organism Item to its map
        Item organism;
        String taxonId = getTaxonId();
        if (organisms.containsKey(taxonId)) {
            organism = organisms.get(taxonId);
        } else {
            organism = createItem("Organism");
            organism.setAttribute("taxonId", taxonId);
            store(organism);
            organisms.put(taxonId, organism);
        }

        // one strain per file so no need to store in a map
        Item strain = createItem("Strain");
        String strainName = getStrainName();
        strain.setAttribute("primaryIdentifier", strainName);
        store(strain);
        
        // ---------------------------------------------------------------------------------------------------------------------------
        // Load and associate genes with proteins and data
        // ---------------------------------------------------------------------------------------------------------------------------

        BufferedReader br = new BufferedReader(reader);
        String line = null;
        while ((line=br.readLine())!=null) {
            // comment line
            if (line.startsWith("#")) continue;

            // String pacId;
            // String locusName;
            // String transcriptName;
            // String peptideName;
            // String[] Pfam;
            // String[] Panther;
            // String[] KOG;
            // String[] ec;
            // String[] KO;
            // String[] GO;
            // String bestHitAtName;
            // String bestHitAtSymbol;
            // String bestHitAtDefline;
            AnnotInfoRecord annot = new AnnotInfoRecord(line);
            if (annot.pacId!=null) {
                // the gene
                Item gene = null;
                if (genes.containsKey(annot.locusName)) {
                    gene = genes.get(annot.locusName);
                } else {
                    gene = createItem("Gene");
                    gene.setAttribute("primaryKey", annot.locusName);
                    genes.put(annot.locusName, gene);
                }
                // the protein
                Item protein = null;
                if (proteins.containsKey(annot.peptideName)) {
                    protein = proteins.get(annot.peptideName);
                } else {
                    protein = createItem("Protein");
                    protein.setAttribute("primaryKey", annot.peptideName);
                    protein.setReference("gene", gene);
                    proteins.put(annot.peptideName, protein);
                }
            }
        }
        br.close();
    }

    /**
     * Store the items we've collected from the AnnotInfo files
     */
    @Override
    public void close() throws ObjectStoreException {
        store(genes.values());
        store(proteins.values());
    }

    /**
     * Get the taxon ID from the current file name, e.g. phavu.G19833.gnm2.ann1.PB8d.info_annot.txt
     */
    public String getTaxonId() {
        String fileName = getCurrentFile().getName();
        String[] chunks = fileName.split(".");
        String org = chunks[0]; // e.g. phavu
        return lisOrgTaxonIds.get(org);
    }

    /**
     * Get the strain name from the current file name, e.g. phavu.G19833.gnm2.ann1.PB8d.info_annot.txt
     */
    public String getStrainName() {
        String fileName = getCurrentFile().getName();
        String[] chunks = fileName.split(".");
        return chunks[1];
    }
}
