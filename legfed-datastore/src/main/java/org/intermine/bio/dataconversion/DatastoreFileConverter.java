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
import java.io.File;
import java.io.FileReader;
import java.io.Reader;
import java.io.IOException;
import java.io.InputStream;

import java.util.Arrays;
import java.util.Map;
import java.util.HashMap;
import java.util.Set;
import java.util.HashSet;
import java.util.Properties;

import org.apache.log4j.Logger;

import org.intermine.dataconversion.FileConverter;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.metadata.StringUtil;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Item;

/**
 * Loads data from LIS datastore files, using the file names to detect what sort of files they are, and,
 * sometimes, what the organism is.
 *
 * @author Sam Hokin
 */
public class DatastoreFileConverter extends FileConverter {
	
    private static final Logger LOG = Logger.getLogger(DatastoreFileConverter.class);

    // store everything in maps
    Map<String,Item> organisms = new HashMap<>();
    Map<String,Item> strains = new HashMap<>();
    Map<String,Item> ontologyTerms = new HashMap<>();
    Map<String,Item> dataSets = new HashMap<>();
    Map<String,Item> geneFamilies = new HashMap<>();
    Map<String,Item> proteinDomains = new HashMap<>();
    Map<String,Item> linkageGroups = new HashMap<>();
    Map<String,Item> geneticMarkers = new HashMap<>();
    Map<String,Item> chromosomes = new HashMap<>();
    
    Set<String> ontologyAnnotations = new HashSet<>(); // concatenation of subject and term identifiers to avoid dupes

    // these are in maps of maps since we can have the same secondaryIdentifier from multiple assembly/annotation versions
    // these are merged on secondaryIdentifier,assemblyVersion,annotationVersion
    Map<String,Map<String,Item>> genes = new HashMap<String,Map<String,Item>>();
    Map<String,Map<String,Item>> proteins = new HashMap<String,Map<String,Item>>();
    Map<String,Map<String,Item>> mRNAs = new HashMap<String,Map<String,Item>>();

    // Taxon IDs of organisms (probably just one)
    Set<String> taxonIds = new HashSet<>();

    // ontologies created in constructor for random use
    Item geneOntology;
    Item pfamOntology;
    Item pantherOntology;
    Item kogOntology;
    Item ecOntology;
    Item koOntology;

    // there is only one dataSource, set in project.xml
    Item dataSource;

    // optional datasource and dataset attributes set in project.xml
    String dataSourceName;
    String dataSourceUrl;
    String dataSourceDescription;
    String dataSetUrl;
    String dataSetDescription;

    // map gensp and Genus_species to taxonId
    Map<String,String> genspTaxonId = new HashMap<>();
    Map<String,String> genusSpeciesTaxonId = new HashMap<>();

    // Set DataSource fields in project.xml
    public void setDataSourceName(String name) {
        this.dataSourceName = name;
    }
    public void setDataSourceUrl(String url) {
        this.dataSourceUrl = url;
    }
    public void setDataSourceDescription(String description) {
        this.dataSourceDescription = description;
    }

    // set DataSet fields in project.xml
    public void setDataSetUrl(String url) {
        this.dataSetUrl = url;
    }
    public void setDataSetDescription(String description) {
        this.dataSetDescription = description;
    }
        
    /**
     * Create a new DatastoreFileConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public DatastoreFileConverter(ItemWriter writer, Model model) throws ObjectStoreException {
        super(writer, model);
        // GO
        geneOntology = createItem("Ontology");
        geneOntology.setAttribute("name", "GO");
        geneOntology.setAttribute("url", "http://www.geneontology.org");
        store(geneOntology);
        // Pfam
        pfamOntology = createItem("Ontology");
        pfamOntology.setAttribute("name", "Pfam");
        pfamOntology.setAttribute("url", "https://pfam.xfam.org/");
        store(pfamOntology);
        // PANTHER
        pantherOntology = createItem("Ontology");
        pantherOntology.setAttribute("name", "PANTHER");
        pantherOntology.setAttribute("url", "http://www.pantherdb.org/");
        store(pantherOntology);
        // KOG
        kogOntology = createItem("Ontology");
        kogOntology.setAttribute("name", "KOG");
        kogOntology.setAttribute("url", "https://genome.jgi.doe.gov/Tutorial/tutorial/kog.html");
        store(kogOntology);
        // ENZYME
        ecOntology = createItem("Ontology");
        ecOntology.setAttribute("name", "ENZYME");
        ecOntology.setAttribute("url", "https://enzyme.expasy.org/");
        store(ecOntology);
        // KEGG Orthology
        koOntology = createItem("Ontology");
        koOntology.setAttribute("name", "KEGG Orthology");
        koOntology.setAttribute("url", "https://www.genome.jp/kegg/ko.html");
        store(koOntology);

        // get the organism data
        DatastoreUtils dsu = new DatastoreUtils();
        genspTaxonId = dsu.getGenspTaxonId();
        genusSpeciesTaxonId = dsu.getGenusSpeciesTaxonId();
        for (String taxonId : genspTaxonId.values()) {
            taxonIds.add(taxonId);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void process(Reader reader) throws Exception {
        // create and store the DataSource if appropriate
        if (dataSource==null && dataSourceName!=null) {
            dataSource = createItem("DataSource");
            dataSource.setAttribute("name", dataSourceName);
            if (dataSourceUrl!=null) dataSource.setAttribute("url", dataSourceUrl);
            if (dataSourceDescription!=null) dataSource.setAttribute("description", dataSourceDescription);
            store(dataSource);
        }
        // process the file
        if (getCurrentFile().getName().contains("description_")) {
            // description_Phaseolus_vulgaris.yml
            printInfoBlurb(getCurrentFile().getName());
            // processDescriptionFile(reader);
        } else if (getCurrentFile().getName().contains("strains_")) {
            // strains_Phaseolus_vulgaris.yml
            printInfoBlurb(getCurrentFile().getName());
            // processStrainsFile(reader);
        } else if (getCurrentFile().getName().endsWith(".info_annot.txt")) {
            // phavu.G19833.gnm2.ann1.PB8d.info_annot.txt
            printInfoBlurb(getCurrentFile().getName());
            // processInfoAnnotFile(reader);
        } else if (getCurrentFile().getName().endsWith(".info_annot_ahrd.tsv")) {
            // legume.genefam.fam1.M65K.info_annot_ahrd.tsv
            String fastaDirname = getCurrentFile().getParent()+"/"+getCurrentFile().getName().replace("info_annot_ahrd.tsv", "family_fasta");
            printInfoBlurb(getCurrentFile().getName());
            printInfoBlurb(fastaDirname);
            // processInfoAnnotAhrdFile(reader);
        } else if (getCurrentFile().getName().endsWith(".cmap.txt")) {
            // phavu.mixed.map1.7PMp.cmap.txt
            printInfoBlurb(getCurrentFile().getName());
            processCmapFile(reader);
        } else if (getCurrentFile().getName().endsWith(".map.gff3")) {
            // phavu.mixed.map1.7PMp.map.gff3
            printInfoBlurb(getCurrentFile().getName());
            processMapGFFFile(reader);
        } else if (getCurrentFile().getName().endsWith(".flanking_seq.fna")) {
            // phavu.mixed.map1.7PMp.flanking_seq.fna
            printInfoBlurb(getCurrentFile().getName());
            processFlankingSeqFastaFile(reader);
        } else {
            System.out.println("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx");
            System.out.println("NOT PROCESSING "+getCurrentFile().getName());
            System.out.println("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx");
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void close() throws Exception {
        System.out.println("Storing all maps...");
        // straight maps
        store(organisms.values());
        store(strains.values());
        store(ontologyTerms.values());
        store(dataSets.values());
        store(geneFamilies.values());
        store(proteinDomains.values());
        store(linkageGroups.values());
        store(geneticMarkers.values());
        store(chromosomes.values());
        // maps of maps
        for (Map<String,Item> geneMap : genes.values()) {
            store(geneMap.values());
        }
        for (Map<String,Item> proteinMap : proteins.values()) {
            store(proteinMap.values());
        }
        for (Map<String,Item> mRNAMap : mRNAs.values()) {
            store(mRNAMap.values());
        }
    }

    /**
     * Print out info about the current file or directory being processed.
     */
    static void printInfoBlurb(String blurb) {
        LOG.info("Processing file/dir "+blurb);
        System.out.println("####################################################################################################################");
        System.out.println("Processing file/dir "+blurb);
        System.out.println("####################################################################################################################");
    }

    /**
     * Get/add the DataSet Item, formed from the current filename.
     */
    Item getDataSet() {
        if (dataSets.containsKey(getCurrentFile().getName())) {
            return dataSets.get(getCurrentFile().getName());
        } else {
            Item dataSet = createItem("DataSet");
            dataSet.setAttribute("name", getCurrentFile().getName());
            if (dataSetUrl!=null) dataSet.setAttribute("url", dataSetUrl);
            if (dataSetDescription!=null) dataSet.setAttribute("description", dataSetDescription);
            if (dataSource!=null) dataSet.setReference("dataSource", dataSource);
            dataSets.put(getCurrentFile().getName(), dataSet);
            return dataSet;
        }
    }

    /**
     * Get/add the DataSet Item, with a given name (when the current filename is not appropriate).
     */
    Item getDataSet(String name) {
        if (dataSets.containsKey(name)) {
            return dataSets.get(name);
        } else {
            Item dataSet = createItem("DataSet");
            dataSet.setAttribute("name", name);
            if (dataSetUrl!=null) dataSet.setAttribute("url", dataSetUrl);
            if (dataSetDescription!=null) dataSet.setAttribute("description", dataSetDescription);
            if (dataSource!=null) dataSet.setReference("dataSource", dataSource);
            dataSets.put(name, dataSet);
            return dataSet;
        }
    }

    /**
     * Get the taxon ID for a string array like ["something","Phaseolus","vulgaris"]
     */
    String getTaxonId(String[] dashparts) {
        String genus = dashparts[1];
        String species = dashparts[2];
        String taxonId = genusSpeciesTaxonId.get(genus+"_"+species);
        if (taxonId==null) {
            throw new RuntimeException("Taxon ID not available for "+genus+"_"+species);
        }
        return taxonId;
    }

    /**
     * Get the taxon ID for a gensp string like "phavu".
     */
    String getTaxonId(String gensp) {
        String taxonId = genspTaxonId.get(gensp);
        if (taxonId==null) {
            throw new RuntimeException("Taxon ID not available for "+gensp);
        }
        return taxonId;
    }

    /**
     * Get/add the organism Item associated with the given taxon ID.
     */
    Item getOrganism(String taxonId) {
        Item organism;
        if (organisms.containsKey(taxonId)) {
            organism = organisms.get(taxonId);
        } else {
            organism = createItem("Organism");
            organism.setAttribute("taxonId", taxonId);
            organisms.put(taxonId, organism);
        }
        return organism;
    }

    /**
     * Get/add the strain Item associated with the given strain name.
     * Sets the organism reference if created.
     */
    Item getStrain(String strainName, Item organism) {
        Item strain;
        if (strains.containsKey(strainName)) {
            strain = strains.get(strainName);
        } else {
            strain = createItem("Strain");
            strain.setAttribute("primaryIdentifier", strainName);
            strain.setReference("organism", organism);
            strains.put(strainName, strain);
        }
        return strain;
    }

    /**
     * Process an organism description file, which is in YAML format:
     *
     * description_Phaseolus_vulgaris.yml
     *
     * %YAML 1.2
     * ##### Phaseolus vulgaris	
     * organism.taxid:	3885
     * organism.genus:	Phaseolus
     * organism.species:	vulgaris
     * organism.abbrev:	phavu
     * organism.commonName:	common bean
     * organism.description:	Common bean was likely domesticated independently both in Central America and in the Andes....
     */
    void processDescriptionFile(Reader reader) throws IOException {
        // get the organism
        String[] dotparts = getCurrentFile().getName().split("\\.");
        String[] dashparts = dotparts[0].split("_");
        Item organism = getOrganism(getTaxonId(dashparts));
        String genus = null;   // for forming Organism.name
        String species = null; // for forming Organism.name
        // now load the attributes
        BufferedReader br = new BufferedReader(reader);
        String line = null;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("%")) continue;
            if (line.startsWith("#")) continue;
            String[] parts = line.split("\t");
            if (parts.length>1) {
                String attributeName = parts[0].replace("organism.","").replace(":","");
                     // munge
                if (attributeName.equals("taxid")) attributeName = "taxonId";
                if (attributeName.equals("abbrev")) attributeName = "abbreviation";
                String attributeValue = parts[1].trim();
                if (attributeValue.length()>0) {
                    organism.setAttribute(attributeName, attributeValue);
                    if (attributeName.equals("genus")) genus = attributeValue;
                    if (attributeName.equals("species")) species = attributeValue;
                }
            }
        }
        br.close();
        if (genus!=null && species!=null) organism.setAttribute("name", genus+" "+species);
    }

    /**
     * Process a strains file, which is in YAML format:
     *
     * strains_Phaseolus_vulgaris.yml
     *
     * %YAML 1.2
     * #####
     * strain.identifier:	G19833
     * strain.accession:	
     * strain.name:	G19833
     * strain.origin:	Peru
     * strain.description:	Andean landrace G19833 was selected for genome sequencing partly due to its resistance to numerous diseases...
     * #####
     * strain.identifier:	BAT93
     * strain.accession:	PI 633451
     * strain.name:	BAT93
     * strain.origin:	CIAT
     * strain.description:	Accession BAT93 is a Mesomarican line that has been used in numerous breeding projects and trait-mapping studies.
     */
    void processStrainsFile(Reader reader) throws IOException {
        // get the organism
        String[] dotparts = getCurrentFile().getName().split("\\.");
        String[] dashparts = dotparts[0].split("_");
        Item organism = getOrganism(getTaxonId(dashparts));
        // spin through the strain sections
        BufferedReader br = new BufferedReader(reader);
        String strainName = null; // keep track of current strain name
        String line = null;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#") || line.startsWith("%") || line.startsWith("#####")) {
                continue;
            } else {
                // strain attributes
                String[] parts = line.split("\t");
                if (parts.length>1) {
                    String attributeName = parts[0].replace("strain.","").replace(":","");
                    String attributeValue = parts[1].trim();
                    if (attributeName.equals("identifier")) {
                        strainName = attributeValue;
                    }
                    Item strain = getStrain(strainName, organism);
                    // special case identifier --> primaryIdentifier (Annotatable)
                    if (attributeName.equals("identifier")) attributeName = "primaryIdentifier";
                    if (attributeValue.length()>0) {
                        strain.setAttribute(attributeName, attributeValue);
                    }
                }
            }
        }
        br.close();
    }

    /**
     * Process an info_annot.txt file which contains relationships between genes, transcripts, proteins and ontology terms.
     * This will also link genes to proteins. The file name starts with gensp.strainName.assemblyVersion.annotationVersion.
     * 0     1      2    3    4    5          6
     * phavu.G19833.gnm2.ann1.PB8d.info_annot.txt
     *
     * pacId locusName transcriptName peptideName Pfam Panther KOG ec KO GO Best-hit-arabi-name arabi-symbol arabi-defline
     * 37170591 Phvul.001G000400 Phvul.001G000400.1 Phvul.001G000400.1.p PF00504 PTHR21649,PTHR21649:SF24 1.10.3.9 K14172 GO:0016020,GO:0009765 AT1G76570.1 Chlorophyll family protein
     */
    void processInfoAnnotFile(Reader reader) throws IOException, ObjectStoreException {
        // get the dot-separated file information
        String[] dotparts = getCurrentFile().getName().split("\\.");
        String gensp = dotparts[0];
        String strainName = dotparts[1];
        String assemblyVersion = dotparts[2];
        String annotationVersion = dotparts[3];
        String versionKey = assemblyVersion+"."+annotationVersion;
        // get the dataSet, organism and strain
        Item dataSet = getDataSet();
        dataSet.setAttribute("version", versionKey);
        Item organism = getOrganism(getTaxonId(gensp));
        Item strain = getStrain(strainName, organism);
        // spin through the file
        BufferedReader br = new BufferedReader(reader);
        String line = null;
        while ((line=br.readLine())!=null) {
            // comment line
            if (line.startsWith("#")) continue;
            // String bestHitAtName;
            // String bestHitAtSymbol;
            // String bestHitAtDefline;
            InfoAnnotRecord record = new InfoAnnotRecord(line);
            if (record.pacId!=null) {
                // the gene
                Item gene = null;
                String geneIdentifier = record.locusName;
                if (genes.containsKey(geneIdentifier)) {
                    Map<String,Item> geneMap = genes.get(geneIdentifier);
                    if (geneMap.containsKey(versionKey)) {
                        gene = geneMap.get(versionKey);
                    } else {
                        // Phvul.002G040500
                        gene = createItem("Gene");
                        gene.setAttribute("secondaryIdentifier", geneIdentifier);
                        gene.setReference("organism", organism);
                        gene.setReference("strain", strain);
                        gene.addToCollection("dataSets", dataSet);
                        gene.setAttribute("assemblyVersion", assemblyVersion);
                        gene.setAttribute("annotationVersion", annotationVersion);
                        geneMap.put(versionKey, gene);
                    }
                } else {
                    // Phvul.002G040500
                    gene = createItem("Gene");
                    gene.setAttribute("secondaryIdentifier", geneIdentifier);
                    gene.setReference("organism", organism);
                    gene.setReference("strain", strain);
                    gene.addToCollection("dataSets", dataSet);
                    gene.setAttribute("assemblyVersion", assemblyVersion);
                    gene.setAttribute("annotationVersion", annotationVersion);
                    Map<String,Item> geneMap = new HashMap<>();
                    geneMap.put(versionKey, gene);
                    genes.put(geneIdentifier, geneMap);
                }
                // the protein
                Item protein = null;
                String proteinIdentifier = record.peptideName;
                if (proteins.containsKey(proteinIdentifier)) {
                    Map<String,Item> proteinMap = proteins.get(proteinIdentifier);
                    if (proteinMap.containsKey(versionKey)) {
                        protein = proteinMap.get(versionKey);
                    } else {
                        // Phvul.002G040500.1.p --> Phvul.002G040500.1
                        protein = createItem("Protein");
                        protein.setAttribute("secondaryIdentifier", proteinIdentifier);
                        protein.setReference("organism", organism);
                        protein.setReference("strain", strain);
                        protein.addToCollection("dataSets", dataSet);
                        protein.setAttribute("assemblyVersion", assemblyVersion);
                        protein.setAttribute("annotationVersion", annotationVersion);
                        protein.addToCollection("genes", gene);
                        proteinMap.put(versionKey, protein);
                    }
                } else {
                    // Phvul.002G040500.1.p --> Phvul.002G040500.1
                    protein = createItem("Protein");
                    protein.setAttribute("secondaryIdentifier", proteinIdentifier);
                    protein.setReference("organism", organism);
                    protein.setReference("strain", strain);
                    protein.addToCollection("dataSets", dataSet);
                    protein.setAttribute("assemblyVersion", assemblyVersion);
                    protein.setAttribute("annotationVersion", annotationVersion);
                    protein.addToCollection("genes", gene);
                    Map<String,Item> proteinMap = new HashMap<>();
                    proteinMap.put(versionKey, protein);
                    proteins.put(proteinIdentifier, proteinMap);
                }
                // the transcript = mRNA
                Item mRNA = null;
                String mRNAIdentifier = record.transcriptName;
                if (mRNAs.containsKey(mRNAIdentifier)) {
                    Map<String,Item> mRNAMap = mRNAs.get(mRNAIdentifier);
                    if (mRNAMap.containsKey(versionKey)) {
                        mRNA = mRNAMap.get(versionKey);
                    } else {
                        // Phvul.002G040500.1
                        mRNA = createItem("MRNA");
                        mRNA.setAttribute("secondaryIdentifier", mRNAIdentifier);
                        mRNA.setReference("organism", organism);
                        mRNA.setReference("strain", strain);
                        mRNA.addToCollection("dataSets", dataSet);
                        mRNA.setAttribute("assemblyVersion", assemblyVersion);
                        mRNA.setAttribute("annotationVersion", annotationVersion);
                        mRNA.setReference("gene", gene);
                        mRNA.setReference("protein", protein);
                        mRNAMap.put(versionKey, mRNA);
                    }
                } else {
                    // Phvul.002G040500.1
                    mRNA = createItem("MRNA");
                    mRNA.setAttribute("secondaryIdentifier", mRNAIdentifier);
                    mRNA.setReference("organism", organism);
                    mRNA.setReference("strain", strain);
                    mRNA.addToCollection("dataSets", dataSet);
                    mRNA.setAttribute("assemblyVersion", assemblyVersion);
                    mRNA.setAttribute("annotationVersion", annotationVersion);
                    mRNA.setReference("gene", gene);
                    mRNA.setReference("protein", protein);
                    Map<String,Item> mRNAMap = new HashMap<>();
                    mRNAMap.put(versionKey, mRNA);
                    mRNAs.put(mRNAIdentifier, mRNAMap);
                }
                // GO terms
                for (String identifier : record.GO) {
                    Item goTerm = null;
                    if (ontologyTerms.containsKey(identifier)) {
                        goTerm = ontologyTerms.get(identifier);
                    } else {
                        goTerm = createItem("OntologyTerm");
                        goTerm.setAttribute("identifier", identifier);
                        goTerm.setReference("ontology", geneOntology);
                        ontologyTerms.put(identifier, goTerm);
                    }
                    String annotKey = identifier+"_"+versionKey+"_"+geneIdentifier;
                    if (!ontologyAnnotations.contains(annotKey)) {
                        Item goAnnotation = createItem("OntologyAnnotation");
                        goAnnotation.setReference("subject", gene);
                        goAnnotation.setReference("ontologyTerm", goTerm);
                        goAnnotation.addToCollection("dataSets", dataSet);
                        store(goAnnotation);
                        ontologyAnnotations.add(annotKey);
                    }
                }
                // Pfam terms
                for (String identifier : record.Pfam) {
                    Item pfamTerm = null;
                    if (ontologyTerms.containsKey(identifier)) {
                        pfamTerm = ontologyTerms.get(identifier);
                    } else {
                        pfamTerm = createItem("OntologyTerm");
                        pfamTerm.setAttribute("identifier", identifier);
                        pfamTerm.setReference("ontology", pfamOntology);
                        ontologyTerms.put(identifier, pfamTerm);
                    }
                    String annotKey = identifier+"_"+versionKey+"_"+proteinIdentifier;
                    if (!ontologyAnnotations.contains(annotKey)) {
                        Item pfamAnnotation = createItem("OntologyAnnotation");
                        pfamAnnotation.setReference("subject", protein);
                        pfamAnnotation.setReference("ontologyTerm", pfamTerm);
                        pfamAnnotation.addToCollection("dataSets", dataSet);
                        store(pfamAnnotation);
                        ontologyAnnotations.add(annotKey);
                    }
                }
                // Panther terms
                for (String identifier : record.Panther) {
                    Item pantherTerm = null;
                    if (ontologyTerms.containsKey(identifier)) {
                        pantherTerm = ontologyTerms.get(identifier);
                    } else {
                        pantherTerm = createItem("OntologyTerm");
                        pantherTerm.setAttribute("identifier", identifier);
                        pantherTerm.setReference("ontology", pantherOntology);
                        ontologyTerms.put(identifier, pantherTerm);
                    }
                    String annotKey = identifier+"_"+versionKey+"_"+proteinIdentifier;
                    if (!ontologyAnnotations.contains(annotKey)) {
                        Item pantherAnnotation = createItem("OntologyAnnotation");
                        pantherAnnotation.setReference("subject", protein);
                        pantherAnnotation.setReference("ontologyTerm", pantherTerm);
                        pantherAnnotation.addToCollection("dataSets", dataSet);
                        store(pantherAnnotation);
                        ontologyAnnotations.add(annotKey);
                    }                }
                // KOG terms
                for (String identifier : record.KOG) {
                    Item kogTerm = null;
                    if (ontologyTerms.containsKey(identifier)) {
                        kogTerm = ontologyTerms.get(identifier);
                    } else {
                        kogTerm = createItem("OntologyTerm");
                        kogTerm.setAttribute("identifier", identifier);
                        kogTerm.setReference("ontology", kogOntology);
                        ontologyTerms.put(identifier, kogTerm);
                    }
                    String annotKey = identifier+"_"+versionKey+"_"+proteinIdentifier;
                    if (!ontologyAnnotations.contains(annotKey)) {
                        Item kogAnnotation = createItem("OntologyAnnotation");
                        kogAnnotation.setReference("subject", protein);
                        kogAnnotation.setReference("ontologyTerm", kogTerm);
                        kogAnnotation.addToCollection("dataSets", dataSet);
                        store(kogAnnotation);
                        ontologyAnnotations.add(annotKey);
                    }
                }
                // ec terms
                for (String identifier : record.ec) {
                    Item ecTerm = null;
                    if (ontologyTerms.containsKey(identifier)) {
                        ecTerm = ontologyTerms.get(identifier);
                    } else {
                        ecTerm = createItem("OntologyTerm");
                        ecTerm.setAttribute("identifier", identifier);
                        ecTerm.setReference("ontology", ecOntology);
                        ontologyTerms.put(identifier, ecTerm);
                    }
                    String annotKey = identifier+"_"+versionKey+"_"+proteinIdentifier;
                    if (!ontologyAnnotations.contains(annotKey)) {
                        Item ecAnnotation = createItem("OntologyAnnotation");
                        ecAnnotation.setReference("subject", protein);
                        ecAnnotation.setReference("ontologyTerm", ecTerm); 
                        ecAnnotation.addToCollection("dataSets", dataSet);
                        store(ecAnnotation);
                        ontologyAnnotations.add(annotKey);
                    }
                }
                // KO terms
                for (String identifier : record.KO) {
                    Item koTerm = null;
                    if (ontologyTerms.containsKey(identifier)) {
                        koTerm = ontologyTerms.get(identifier);
                    } else {
                        koTerm = createItem("OntologyTerm");
                        koTerm.setAttribute("identifier", identifier);
                        koTerm.setReference("ontology", koOntology);
                        ontologyTerms.put(identifier, koTerm);
                    }
                    String annotKey = identifier+"_"+versionKey+"_"+geneIdentifier;
                    if (!ontologyAnnotations.contains(annotKey)) {
                        Item koAnnotation = createItem("OntologyAnnotation");
                        koAnnotation.setReference("subject", gene);
                        koAnnotation.setReference("ontologyTerm", koTerm);
                        koAnnotation.addToCollection("dataSets", dataSet);
                        store(koAnnotation);
                        ontologyAnnotations.add(annotKey);
                    }
                }
            }
        }
        br.close();
    }

    /**
     * Process an info_annot_ahrd.tsv file which contains gene families and semi-colon separated groups of ontology terms.
     * 0         1       2    3    4               5
     * z1-legfed.genefam.fam1.M65K.info_annot_ahrd.tsv
     * z1-legfed.genefam.fam1.M65K.family_fasta
     *
     * legfed_v1_0.L_H6Q7T0-consensus type I inositol-1,4,5-trisphosphate 5-phosphatase;
     *             ^^^^^^^^           IPR005135 (Endonuclease/exonuclease/phosphatase), IPR015943 (WD40/YVTN repeat-like-containing domain);
     *                                GO:0005515 (protein binding), GO:0046856 (phosphatidylinositol dephosphorylation)
     *
     * zzzzz.genefam.fam1.M65K.family_fasta/L_H6Q7T0
     *                                      ^^^^^^^^
     * >lupan.Lup015831.1
     * ----CASFAKLT--TLSPHWIGNNSFSSRRGGSSPLTATRRVSLPIRASSYSDELVQTAK
     * TIASPGRGILAIDESNATCGKRLASIGLDNTEVNRQAYRQLLLTTPGLGEYISGAILFEE
     * ...
     * >phavu.Phvul.007G033800.1
     * -----------------------------------TFSPRRVSLPIRASSYQQELVQTAK
     * SIASPGRGILAIDESNATCGKRLASIGLDNTEVNRQAYRQLLLTTPGLGEYISGAILFEE
     * ...
     */
    void processInfoAnnotAhrdFile(Reader reader) throws IOException, ObjectStoreException {
        // get the dataSet version from the filename
        String[] dotparts = getCurrentFile().getName().split("\\.");
        String[] dashparts = getCurrentFile().getName().split("-"); // z1- may be prepended to run after other files
        String dataSetName = getCurrentFile().getName();
        if (dashparts.length>1) dataSetName = dashparts[1]; // hopefully no other dashes!
        String dataSetVersion = dotparts[2];
        Item dataSet = getDataSet(dataSetName);
        dataSet.setAttribute("version", dataSetVersion);
        // get the FASTA directory
        String fastaDirname = getCurrentFile().getParent()+"/"+getCurrentFile().getName().replace("info_annot_ahrd.tsv", "family_fasta");
        // spin through the lines
        BufferedReader br = new BufferedReader(reader);
        String line = null;
        while ((line=br.readLine())!=null) {
            // comment line
            if (line.startsWith("#")) continue;
            // parse record and create items
            InfoAnnotAhrdRecord record = new InfoAnnotAhrdRecord(line);
            Item geneFamily = createItem("GeneFamily");
            geneFamily.setAttribute("identifier", record.identifier);
            geneFamily.setAttribute("version", record.version);
            geneFamily.setAttribute("description", record.description);
            geneFamily.setReference("dataSet", dataSet);
            geneFamilies.put(record.identifier, geneFamily);
            // GO terms
            for (String identifier : record.go.keySet()) {
                String description = record.go.get(identifier);
                Item goTerm = null;
                if (ontologyTerms.containsKey(identifier)) {
                    goTerm = ontologyTerms.get(identifier);
                } else {
                    goTerm = createItem("OntologyTerm");
                    goTerm.setAttribute("identifier", identifier);
                    goTerm.setAttribute("description", description);
                    goTerm.setReference("ontology", geneOntology);
                    ontologyTerms.put(identifier, goTerm);
                }
                Item goAnnotation = createItem("OntologyAnnotation");
                goAnnotation.setReference("subject", geneFamily);
                goAnnotation.setReference("ontologyTerm", goTerm);
                goAnnotation.addToCollection("dataSets", dataSet);
                store(goAnnotation);
            }
            // interpro domains
            for (String identifier : record.interpro.keySet()) {
                String description = record.interpro.get(identifier);
                Item proteinDomain = null;
                if (proteinDomains.containsKey(identifier)) {
                    proteinDomain = proteinDomains.get(identifier);
                } else {
                    proteinDomain = createItem("ProteinDomain");
                    proteinDomain.setAttribute("primaryIdentifier", identifier);
                    proteinDomain.setAttribute("description", description);
                    proteinDomains.put(identifier, proteinDomain);
                }
                proteinDomain.addToCollection("geneFamilies", geneFamily);
            }
            // load the gene family FASTA if present to link proteins of the desired organism
            String fastaFilename = fastaDirname+"/"+record.identifier;
            File fastaFile = new File(fastaFilename);
            if (fastaFile.exists()) {
                BufferedReader fbr = new BufferedReader(new FileReader(fastaFile));
                String fline = null;
                while ((fline=fbr.readLine())!=null) {
                    if (fline.startsWith(">")) {
                        String genspName = fline.substring(1);
                        String[] parts = genspName.split("\\.");
                        String gensp = parts[0];
                        String name = parts[1];
                        for (int i=2; i<parts.length; i++) name += "."+parts[i];
                        String taxonId = genspTaxonId.get(gensp);
                        if (taxonId!=null && taxonIds.contains(taxonId)) {
                            if (proteins.containsKey(name)) {
                                // update proteins we already have (hopefully most of them)
                                for (Item protein : proteins.get(name).values()) {
                                    protein.setReference("geneFamily", geneFamily);
                                    protein.addToCollection("dataSets", dataSet);
                                }
                            } else {
                                // create protein without assemblyVersion, annotationVersion if not already loaded (boo)
                                Item protein = createItem("Protein");
                                protein.setAttribute("secondaryIdentifier", name);
                                protein.setReference("geneFamily", geneFamily);
                                protein.addToCollection("dataSets", dataSet);
                                Map<String,Item> proteinMap = new HashMap<>();
                                proteinMap.put(dataSetVersion, protein);
                                proteins.put(name, proteinMap);
                            }
                        }
                    }
                }
                fbr.close();
            }
        }
        br.close();
    }

    /**
     * Process a genetic Cmap file.
     *
     * 0     1     2    3    4    5
     * phavu.mixed.map1.7PMp.cmap.txt
     *
     * map_acc                map_name map_start map_stop feature_acc                feature_name feature_aliases feature_start feature_stop feature_type_acc is_landmark
     * PvCookUCDavis2009_Pv09 Pv09     0         66.1     PvCook...Pv09_Pv_TOG913042 Pv_TOG913042 TOG913042       66.1          66.1         SNP              0
     * PvCookUCDavis2009_Pv09 Pv09     0         66.1     PvCook...Pv09_Pv_TOG899751 Pv_TOG899751 TOG899751       44.4          44.4         SNP              0
     */
    void processCmapFile(Reader reader) throws IOException, ObjectStoreException {
        // get the identifiers from the filename
        String[] dotparts = getCurrentFile().getName().split("\\.");
        String dataSetName = getCurrentFile().getName();
        String gensp = dotparts[0];
        String parents = dotparts[1];
        String dataSetVersion = dotparts[2];
        // create Items
        Item dataSet = getDataSet(dataSetName);
        dataSet.setAttribute("version", dataSetVersion);
        Item organism = getOrganism(getTaxonId(gensp));
        // spin through the lines
        BufferedReader br = new BufferedReader(reader);
        String line = null;
        while ((line=br.readLine())!=null) {
            CMapRecord record = new CMapRecord(line);
            if (!record.hasData()) continue;
            Item lg = getLinkageGroup(record.map_acc);
            lg.setAttribute("secondaryIdentifier", record.map_name);
            lg.setAttribute("length", String.valueOf(record.map_stop));
            lg.setReference("organism", organism);
            lg.setReference("dataSet", dataSet);
            String fullname = record.feature_name;
            String name = fullname;
            if (record.feature_aliases!=null) name = record.feature_aliases;
            Item gm = getGeneticMarker(name);
            gm.setAttribute("type", String.valueOf(record.feature_type_acc));
            gm.setReference("organism", organism);
            gm.addToCollection("dataSets", dataSet);
            Item lgp = createItem("LinkageGroupPosition");
            lgp.setAttribute("position", String.valueOf(record.feature_start));
            lgp.setReference("linkageGroup", lg);
            store(lgp);
            gm.addToCollection("linkageGroupPositions", lgp);
            lg.addToCollection("markers", gm);
        }
        br.close();
    }

    /**
     * Process a flanking sequence FASTA file -- store the flanking sequence as the genetic marker's "sequence".
     *
     * 0     1     2    3    4            5
     * phavu.mixed.map1.7PMp.flanking_seq.fna
     *
     * >TOG905303_749
     * GACACGTAACTGAAATTTCACYATCCCTTACTGTTTCTCAAATCTTAGGATGAAAATAGATGCAATTGCTGGAAGTTCCCAATATTTTCTTTTRCACCTATGTACTCAGGTAATTTTATAT
     */
    void processFlankingSeqFastaFile(Reader reader) throws IOException, ObjectStoreException {
        // get the identifiers from the filename
        String[] dotparts = getCurrentFile().getName().split("\\.");
        String dataSetName = getCurrentFile().getName();
        String gensp = dotparts[0];
        String parents = dotparts[1];
        String dataSetVersion = dotparts[2];
        // create Items
        Item dataSet = getDataSet(dataSetName);
        dataSet.setAttribute("version", dataSetVersion);
        Item organism = getOrganism(getTaxonId(gensp));
        // spin through the lines
        BufferedReader br = new BufferedReader(reader);
        String line = null;
        while ((line=br.readLine())!=null) {
            if (!line.startsWith(">")) continue;
            String fullname = line.substring(1);
            String name = fullname;
            String[] underscoreparts = fullname.split("_");
            if (underscoreparts.length==3) {
                // strip front and back of Pv_TOG905303_749
                name = underscoreparts[1];
            } else if (underscoreparts.length==2) {
                // strip back of TOG905303_749
                name = underscoreparts[0];
            }
            Item gm = getGeneticMarker(name);
            if (!name.equals(fullname)) gm.setAttribute("secondaryIdentifier", fullname);
            String residues = br.readLine();
            Item sequence = createItem("Sequence");
            sequence.setAttribute("residues", residues);
            sequence.setAttribute("length", String.valueOf(residues.length()));
            store(sequence);
            gm.setReference("sequence", sequence);
        }
        br.close();
    }

    /**
     * Process a genetic map GFF file.
     *
     * 0     1     2    3    4   5
     * phavu.mixed.map1.7PMp.map.gff3
     *
     * ##gff-version 3
     * ##date Sat Jul  8 11:12:08 2017
     * ##source gbrowse gbgff gff3 dumper, from https://legumeinfo.org/genomes/gbrowse/Pv1.0 with post-processing
     * phavu.G19833.gnm1.Chr01 blastn  SNP     242423  242423  .       +       .       Name=Pv_TOG905303_749;ID=1;Note=LG01 cM 0.0;alleles=T/G
     */
    void processMapGFFFile(Reader reader) throws IOException, ObjectStoreException {
        // get the identifiers from the filename
        String[] dotparts = getCurrentFile().getName().split("\\.");
        String dataSetName = getCurrentFile().getName();
        String gensp = dotparts[0];
        String parents = dotparts[1];
        String dataSetVersion = dotparts[2];
        // create Items
        Item dataSet = getDataSet(dataSetName);
        dataSet.setAttribute("version", dataSetVersion);
        Item organism = getOrganism(getTaxonId(gensp));
        // spin through the lines
        BufferedReader br = new BufferedReader(reader);
        String line = null;
        while ((line=br.readLine())!=null) {
            MapGFFRecord record = new MapGFFRecord(line);
            if (!record.hasData()) continue;
            Item chr = getChromosomeOrSupercontig(record.chr);
            chr.setReference("organism", organism);
            chr.addToCollection("dataSets", dataSet);
            Item gm = getGeneticMarker(record.name);
            gm.setReference("organism", organism);
            gm.addToCollection("dataSets", dataSet);
            if (!record.name.equals(record.fullname)) gm.setAttribute("secondaryIdentifier", record.fullname);
            gm.setAttribute("type", record.type);
            gm.setAttribute("length", String.valueOf(record.end-record.start+1));
            if (record.alleles!=null) gm.setAttribute("alleles", record.alleles);
            Item location = createItem("Location");
            location.setAttribute("strand", record.strand);
            location.setAttribute("start", String.valueOf(record.start));
            location.setAttribute("end", String.valueOf(record.end));
            location.setReference("feature", gm);
            location.addToCollection("dataSets", dataSet);
            location.setReference("locatedOn", chr);
            store(location);
            if (chr.getClassName().equals("Supercontig")) {
                gm.setReference("supercontig", chr);
                gm.setReference("supercontigLocation", location);
            } else {
                gm.setReference("chromosome", chr);
                gm.setReference("chromosomeLocation", location);
            }
        }
        br.close();
    }

    /**
     * Get an old LinkageGroup or create a new one
     */
    public Item getLinkageGroup(String primaryIdentifier) {
        if (linkageGroups.containsKey(primaryIdentifier)) {
            return linkageGroups.get(primaryIdentifier);
        } else {
            Item lg = createItem("LinkageGroup");
            lg.setAttribute("primaryIdentifier", primaryIdentifier);
            linkageGroups.put(primaryIdentifier, lg);
            return lg;
        }
    }

    /**
     * Get an old GeneticMarker or create a new one
     */
    public Item getGeneticMarker(String primaryIdentifier) {
        if (geneticMarkers.containsKey(primaryIdentifier)) {
            return geneticMarkers.get(primaryIdentifier);
        } else {
            Item gm = createItem("GeneticMarker");
            gm.setAttribute("primaryIdentifier", primaryIdentifier);
            geneticMarkers.put(primaryIdentifier, gm);
            return gm;
        }
    }

    /**
     * Get an existing Chromosome/Supercontig record or create a new one.
     */
    public Item getChromosomeOrSupercontig(String primaryIdentifier) {
        if (chromosomes.containsKey(primaryIdentifier)) {
            return chromosomes.get(primaryIdentifier);
        } else {
            Item chr;
            // HACK: change the className to "Supercontig" if identifier contains "scaffold" or ends in "sc" etc.
            String[] dotparts = primaryIdentifier.split("\\.");
            String lastPart = dotparts[dotparts.length-1];
            if (primaryIdentifier.toLowerCase().contains("scaffold") || lastPart.contains("sc")) {
                chr = createItem("Supercontig");
            } else {
                chr = createItem("Chromosome");
            }
            chr.setAttribute("primaryIdentifier", primaryIdentifier);
            chromosomes.put(primaryIdentifier, chr);
            return chr;
        }
    }
}