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
 * 37170591 Phvul.001G000400 Phvul.001G000400.1 Phvul.001G000400.1.p PF00504 PTHR21649,PTHR21649:SF24 1.10.3.9 K14172 GO:0016020,GO:0009765 AT1G76570.1	Chlorophyll family protein
 *
 * Presumption: one file per strain. (Could be several strains for an organism.)
 * NOTE: this will also link genes to proteins.
 *
 * @author Sam Hokin
 */
public class AnnotInfoFileConverter extends BioFileConverter {
	
    private static final Logger LOG = Logger.getLogger(AnnotInfoFileConverter.class);

    // the file with our info
    String filename;

    // optional dataset and datasource
    String dataSourceName;
    String dataSourceUrl;
    String dataSourceDescription;
    // dataSet name is the file name
    String dataSetUrl;
    String dataSetDescription;

    // store the items in maps for the usual reasons
    Map<String,Item> genes = new HashMap<>();
    Map<String,Item> proteins = new HashMap<>();
    Map<String,Item> mRNAs = new HashMap<>();
    Map<String,Item> ontologyTerms = new HashMap<>();

    // the ontologies are created/stored only once in constructor
    Item geneOntology, pfamOntology, pantherOntology, kogOntology, ecOntology, koOntology;

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

    // set the file of interest
    public void setSrcDataFile(String filename) {
        this.filename = filename;
    }

    // set DataSource fields
    public void setDataSourceName(String name) {
        this.dataSourceName = name;
    }
    public void setDataSourceUrl(String url) {
        this.dataSourceUrl = url;
    }
    public void setDataSourceDescription(String description) {
        this.dataSourceDescription = description;
    }

    // set DataSet fields
    public void setDataSetUrl(String url) {
        this.dataSetUrl = url;
    }
    public void setDataSetDescription(String description) {
        this.dataSetDescription = description;
    }

    /**
     * Create a new AnnotInfoFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public AnnotInfoFileConverter(ItemWriter writer, Model model) {
        super(writer, model);
        // create/store the ontologies
        try {
            geneOntology = createItem("Ontology");
            geneOntology.setAttribute("name", "GO");
            geneOntology.setAttribute("url", "http://www.geneontology.org");
            store(geneOntology);

            pfamOntology = createItem("Ontology");
            pfamOntology.setAttribute("name", "Pfam");
            pfamOntology.setAttribute("url", "https://pfam.xfam.org/");
            store(pfamOntology);

            pantherOntology = createItem("Ontology");
            pantherOntology.setAttribute("name", "PANTHER");
            pantherOntology.setAttribute("url", "http://www.pantherdb.org/");
            store(pantherOntology);

            kogOntology = createItem("Ontology");
            kogOntology.setAttribute("name", "KOG");
            kogOntology.setAttribute("url", "https://genome.jgi.doe.gov/Tutorial/tutorial/kog.html");
            store(kogOntology);
            
            ecOntology = createItem("Ontology");
            ecOntology.setAttribute("name", "ENZYME");
            ecOntology.setAttribute("url", "https://enzyme.expasy.org/");
            store(ecOntology);

            koOntology = createItem("Ontology");
            koOntology.setAttribute("name", "KEGG Orthology");
            koOntology.setAttribute("url", "https://www.genome.jp/kegg/ko.html");
            store(koOntology);
        } catch (ObjectStoreException e) {
            System.err.println(e);
            System.exit(1);
        }
    }

    /**
     * {@inheritDoc}
     * Read in the AnnotInfo file and store the features.
     */
    @Override
    public void process(Reader reader) throws Exception {

        // don't process files other than the one we want
        if (!getCurrentFile().getName().equals(filename)) return;

        LOG.info("Processing info_annot file "+getCurrentFile().getName()+"...");

        // get the organism, strain, assembly and annotation from the file name
        // 0     1      2    3    4    5          6
        // phavu.G19833.gnm2.ann1.PB8d.info_annot.txt

        String[] parts = filename.split("\\.");
        String gensp = parts[0];
        String strainName = parts[1];
        String assemblyVersion = parts[2];
        String annotationVersion = parts[3];
        
        // create the Organism Item
        String taxonId = lisOrgTaxonIds.get(gensp);
        Item organism = createItem("Organism");
        organism.setAttribute("taxonId", taxonId);
        store(organism);

        // create the Strain Item
        Item strain = createItem("Strain");
        strain.setAttribute("primaryIdentifier", strainName);
        store(strain);

        // DataSource (optional)
        Item dataSource = null;
        if (dataSourceName!=null) {
            dataSource = createItem("DataSource");
            dataSource.setAttribute("name", dataSourceName);
            if (dataSourceUrl!=null) dataSource.setAttribute("url", dataSourceUrl);
            if (dataSourceDescription!=null) dataSource.setAttribute("description", dataSourceDescription);
            store(dataSource);
        }

        // DataSet (at least filename and version)
        Item dataSet = createItem("DataSet");
        dataSet.setAttribute("name", filename);
        dataSet.setAttribute("version", assemblyVersion+"."+annotationVersion);
        if (dataSetUrl!=null) dataSet.setAttribute("url", dataSetUrl);
        if (dataSetDescription!=null) dataSet.setAttribute("description", dataSetDescription);
        if (dataSource!=null) dataSet.setReference("dataSource", dataSource);
        store(dataSet);
        
        // ---------------------------------------------------------------------------------------------------------------------------
        // Load and associate genes with proteins and data
        // ---------------------------------------------------------------------------------------------------------------------------

        BufferedReader br = new BufferedReader(reader);
        String line = null;
        while ((line=br.readLine())!=null) {
            // comment line
            if (line.startsWith("#")) continue;

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
                    // Phvul.002G040500
                    gene = createItem("Gene");
                    gene.setAttribute("secondaryIdentifier", annot.locusName);
                    gene.setReference("organism", organism);
                    gene.setReference("strain", strain);
                    gene.addToCollection("dataSets", dataSet);
                    if (assemblyVersion!=null) gene.setAttribute("assemblyVersion", assemblyVersion);
                    if (annotationVersion!=null) gene.setAttribute("annotationVersion", annotationVersion);
                    genes.put(annot.locusName, gene);
                }
                // the protein
                Item protein = null;
                if (proteins.containsKey(annot.peptideName)) {
                    protein = proteins.get(annot.peptideName);
                } else {
                    // Phvul.002G040500.1.p --> Phvul.002G040500.1
                    protein = createItem("Protein");
                    protein.setAttribute("secondaryIdentifier", annot.peptideName);
                    protein.setReference("organism", organism);
                    protein.setReference("strain", strain);
                    protein.addToCollection("dataSets", dataSet);
                    if (assemblyVersion!=null) protein.setAttribute("assemblyVersion", assemblyVersion);
                    if (annotationVersion!=null) protein.setAttribute("annotationVersion", annotationVersion);
                    protein.addToCollection("genes", gene);
                    proteins.put(annot.peptideName, protein);
                }
                // the transcript = mRNA
                Item mRNA = null;
                if (mRNAs.containsKey(annot.transcriptName)) {
                    mRNA = mRNAs.get(annot.transcriptName);
                } else {
                    // Phvul.002G040500.1
                    mRNA = createItem("MRNA");
                    mRNA.setAttribute("secondaryIdentifier", annot.transcriptName);
                    mRNA.setReference("organism", organism);
                    mRNA.setReference("strain", strain);
                    mRNA.addToCollection("dataSets", dataSet);
                    if (assemblyVersion!=null) mRNA.setAttribute("assemblyVersion", assemblyVersion);
                    if (annotationVersion!=null) mRNA.setAttribute("annotationVersion", annotationVersion);
                    mRNA.setReference("gene", gene);
                    mRNA.setReference("protein", protein);
                    mRNAs.put(annot.transcriptName, mRNA);
                }
                // GO terms
                for (String identifier : annot.GO) {
                    Item goTerm = null;
                    if (ontologyTerms.containsKey(identifier)) {
                        goTerm = ontologyTerms.get(identifier);
                    } else {
                        goTerm = createItem("OntologyTerm");
                        goTerm.setAttribute("identifier", identifier);
                        goTerm.setReference("ontology", geneOntology);
                        ontologyTerms.put(identifier, goTerm);
                        Item goAnnotation = createItem("OntologyAnnotation");
                        goAnnotation.setReference("subject", gene);
                        goAnnotation.setReference("ontologyTerm", goTerm);
                        goAnnotation.addToCollection("dataSets", dataSet);
                        store(goAnnotation);
                    }
                }
                // Pfam terms
                for (String identifier : annot.Pfam) {
                    Item pfamTerm = null;
                    if (ontologyTerms.containsKey(identifier)) {
                        pfamTerm = ontologyTerms.get(identifier);
                    } else {
                        pfamTerm = createItem("OntologyTerm");
                        pfamTerm.setAttribute("identifier", identifier);
                        pfamTerm.setReference("ontology", pfamOntology);
                        ontologyTerms.put(identifier, pfamTerm);
                        Item pfamAnnotation = createItem("OntologyAnnotation");
                        pfamAnnotation.setReference("subject", protein);
                        pfamAnnotation.setReference("ontologyTerm", pfamTerm);
                        pfamAnnotation.addToCollection("dataSets", dataSet);
                        store(pfamAnnotation);
                    }
                }
                // Panther terms
                for (String identifier : annot.Panther) {
                    Item pantherTerm = null;
                    if (ontologyTerms.containsKey(identifier)) {
                        pantherTerm = ontologyTerms.get(identifier);
                    } else {
                        pantherTerm = createItem("OntologyTerm");
                        pantherTerm.setAttribute("identifier", identifier);
                        pantherTerm.setReference("ontology", pantherOntology);
                        ontologyTerms.put(identifier, pantherTerm);
                        Item pantherAnnotation = createItem("OntologyAnnotation");
                        pantherAnnotation.setReference("subject", protein);
                        pantherAnnotation.setReference("ontologyTerm", pantherTerm);
                        pantherAnnotation.addToCollection("dataSets", dataSet);
                        store(pantherAnnotation);
                    }
                }
                // KOG terms
                for (String identifier : annot.KOG) {
                    Item kogTerm = null;
                    if (ontologyTerms.containsKey(identifier)) {
                        kogTerm = ontologyTerms.get(identifier);
                    } else {
                        kogTerm = createItem("OntologyTerm");
                        kogTerm.setAttribute("identifier", identifier);
                        kogTerm.setReference("ontology", kogOntology);
                        ontologyTerms.put(identifier, kogTerm);
                        Item kogAnnotation = createItem("OntologyAnnotation");
                        kogAnnotation.setReference("subject", protein);
                        kogAnnotation.setReference("ontologyTerm", kogTerm);
                        kogAnnotation.addToCollection("dataSets", dataSet);
                        store(kogAnnotation);
                    }
                }
                // ec terms
                for (String identifier : annot.ec) {
                    Item ecTerm = null;
                    if (ontologyTerms.containsKey(identifier)) {
                        ecTerm = ontologyTerms.get(identifier);
                    } else {
                        ecTerm = createItem("OntologyTerm");
                        ecTerm.setAttribute("identifier", identifier);
                        ecTerm.setReference("ontology", ecOntology);
                        ontologyTerms.put(identifier, ecTerm);
                        Item ecAnnotation = createItem("OntologyAnnotation");
                        ecAnnotation.setReference("subject", protein);
                        ecAnnotation.setReference("ontologyTerm", ecTerm); 
                        ecAnnotation.addToCollection("dataSets", dataSet);
                       store(ecAnnotation);
                    }
                }
                // KO terms
                for (String identifier : annot.KO) {
                    Item koTerm = null;
                    if (ontologyTerms.containsKey(identifier)) {
                        koTerm = ontologyTerms.get(identifier);
                    } else {
                        koTerm = createItem("OntologyTerm");
                        koTerm.setAttribute("identifier", identifier);
                        koTerm.setReference("ontology", koOntology);
                        ontologyTerms.put(identifier, koTerm);
                        Item koAnnotation = createItem("OntologyAnnotation");
                        koAnnotation.setReference("subject", gene);
                        koAnnotation.setReference("ontologyTerm", koTerm);
                        koAnnotation.addToCollection("dataSets", dataSet);
                        store(koAnnotation);
                    }
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
        store(mRNAs.values());
        store(ontologyTerms.values());
    }
}
