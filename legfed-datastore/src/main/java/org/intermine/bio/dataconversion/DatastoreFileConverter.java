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
import java.io.IOException;

import java.util.Arrays;
import java.util.Map;
import java.util.HashMap;
import java.util.Set;
import java.util.HashSet;

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
    Map<String,Item> genes = new HashMap<>();
    Map<String,Item> proteins = new HashMap<>();
    Map<String,Item> mRNAs = new HashMap<>();
    Map<String,Item> ontologyTerms = new HashMap<>();
    Map<String,Item> dataSets = new HashMap<>();
    Map<String,Item> geneFamilies = new HashMap<>();
    Map<String,Item> proteinDomains = new HashMap<>();

    // Taxon IDs of organisms of interest, set in project.xml
    Set<String> taxonIds;

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

    // map gensp to taxon Ids for associating dir/file names with organisms
    private static Map<String,String> genspTaxonIds = new HashMap<String,String>() {{
            put("aradu","130453");
            put("arahy","3818");
            put("araip","130454");
            put("cajan","3821");
            put("cicar","3827");
            put("glyma","3847");
            put("lotja","34305");
            put("lupan","3871");
            put("medtr","3880");
            put("phavu","3885");
            put("tripr","57577 ");
            put("vigan","3914");
            put("vigra","157791");
            put("vigun","3920");
        }};

    // map Genus_species to taxon Ids for associating dir/file names with organisms
    private static Map<String,String> genusSpeciesTaxonIds = new HashMap<String,String>() {{
            put("Arachis_duranensis","130453");
            put("Arachis_hypogaea","3818");
            put("Arachis_ipaensis","130454");
            put("Glycine_max","3847");
            put("Medicago_truncatula","3880");
            put("Phaseolus_vulgaris","3885");
        }};

    /**
     * Set the list of taxonIds that should be imported in cases where multiple organism data are present (like gene families).
     *
     * @param taxonIds a space-separated list of taxonIds
     */
    public void setOrganisms(String taxonIdString) {
        this.taxonIds = new HashSet<String>(Arrays.asList(StringUtil.split(taxonIdString, " ")));
        System.out.println("Will store for organisms:"+taxonIdString);
    }

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
        if (getCurrentFile().getName().startsWith("description_")) {
            // Phaseolus_vulgaris/about_this_collection/description_Phaseolus_vulgaris.yml
            printInfo();
            processDescriptionFile(reader);
        } else if (getCurrentFile().getName().startsWith("strains_")) {
            // Phaseolus_vulgaris/about_this_collection/strains_Phaseolus_vulgaris.yml
            printInfo();
            processStrainsFile(reader);
        } else if (getCurrentFile().getName().endsWith(".info_annot.txt")) {
            // Phaseolus_vulgaris/G19833.gnm2.ann1.PB8d/phavu.G19833.gnm2.ann1.PB8d.info_annot.txt
            printInfo();
            processInfoAnnotFile(reader);
        } else if (getCurrentFile().getName().endsWith(".info_annot_ahrd.tsv")) {
            // Gene_families/legume.genefam.fam1.M65K/legume.genefam.fam1.M65K.info_annot_ahrd.tsv
            printInfo();
            processInfoAnnotAhrdFile(reader);
        } else if (getCurrentFile().getParent().endsWith(".family_fasta")) {
            // Gene_families/legume.genefam.fam1.M65K/legume.genefam.fam1.M65K.family_fasta
            processFamilyFastaFile(reader);
        }
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public void close() throws Exception {
        System.out.println("Storing all maps...");
        store(organisms.values());
        store(strains.values());
        store(genes.values());
        store(proteins.values());
        store(mRNAs.values());
        store(ontologyTerms.values());
        store(dataSets.values());
        store(geneFamilies.values());
        store(proteinDomains.values());
    }

    /**
     * Print out info about the current file being processed.
     */
    void printInfo() {
        LOG.info("Processing file "+getCurrentFile().getName());
        System.out.println("#################################################################");
        System.out.println("Processing file "+getCurrentFile().getName());
        System.out.println("#################################################################");
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
        String taxonId = genusSpeciesTaxonIds.get(genus+"_"+species);
        if (taxonId==null) {
            throw new RuntimeException("Taxon ID not available for "+genus+"_"+species);
        }
        return taxonId;
    }

    /**
     * Get the taxon ID for a gensp string like "phavu".
     */
    String getTaxonId(String gensp) {
        String taxonId = genspTaxonIds.get(gensp);
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
     * organism.taxonId:	3885
     * organism.genus:	Phaseolus
     * organism.species:	vulgaris
     * organism.name:	Phaseolus vulgaris
     * organism.commonName:	common bean
     * organism.description:	Common bean was likely domesticated independently both in Central America and in the Andes....
     */
    void processDescriptionFile(Reader reader) throws IOException {
        // get the organism
        String[] dotparts = getCurrentFile().getName().split("\\.");
        String[] dashparts = dotparts[0].split("_");
        Item organism = getOrganism(getTaxonId(dashparts));
        // now load the attributes
        BufferedReader br = new BufferedReader(reader);
        String line = null;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("%")) continue;
            if (line.startsWith("#")) continue;
            String[] parts = line.split("\t");
            if (parts.length>1) {
                String attributeName = parts[0].replace("organism.","").replace(":","");
                String attributeValue = parts[1].trim();
                organism.setAttribute(attributeName, attributeValue);
            }
        }
        br.close();
    }

    /**
     * Process a strains file, which is in YAML format:
     *
     * strains_Phaseolus_vulgaris.yml
     *
     * %YAML 1.2
     * ##### G19833	
     * strain.identifier:	G19833
     * strain.accession:	
     * strain.name:	G19833
     * strain.origin:	Peru
     * strain.description:	Andean landrace G19833 was selected for genome sequencing partly due to its resistance to numerous diseases...
     * ##### BAT93	
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
            if (line.startsWith("%")) continue;
            if (line.startsWith("#####")) {
                // new strain
                strainName = line.substring(6).trim();
                Item strain = getStrain(strainName, organism);
            } else if (line.startsWith("#")) {
                // comment
                continue;
            } else {
                // strain attributes
                Item strain = getStrain(strainName, organism);
                String[] parts = line.split("\t");
                if (parts.length>1) {
                    String attributeName = parts[0].replace("strain.","").replace(":","");
                    String attributeValue = parts[1].trim();
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
        // get the dataSet, organism and strain
        Item dataSet = getDataSet();
        dataSet.setAttribute("version", assemblyVersion+"."+annotationVersion);
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
                if (genes.containsKey(record.locusName)) {
                    gene = genes.get(record.locusName);
                } else {
                    // Phvul.002G040500
                    gene = createItem("Gene");
                    gene.setAttribute("secondaryIdentifier", record.locusName);
                    gene.setReference("organism", organism);
                    gene.setReference("strain", strain);
                    gene.addToCollection("dataSets", dataSet);
                    gene.setAttribute("assemblyVersion", assemblyVersion);
                    gene.setAttribute("annotationVersion", annotationVersion);
                    genes.put(record.locusName, gene);
                }
                // the protein
                Item protein = null;
                if (proteins.containsKey(record.peptideName)) {
                    protein = proteins.get(record.peptideName);
                } else {
                    // Phvul.002G040500.1.p --> Phvul.002G040500.1
                    protein = createItem("Protein");
                    protein.setAttribute("secondaryIdentifier", record.peptideName);
                    protein.setReference("organism", organism);
                    protein.setReference("strain", strain);
                    protein.addToCollection("dataSets", dataSet);
                    protein.setAttribute("assemblyVersion", assemblyVersion);
                    protein.setAttribute("annotationVersion", annotationVersion);
                    protein.addToCollection("genes", gene);
                    proteins.put(record.peptideName, protein);
                }
                // the transcript = mRNA
                Item mRNA = null;
                if (mRNAs.containsKey(record.transcriptName)) {
                    mRNA = mRNAs.get(record.transcriptName);
                } else {
                    // Phvul.002G040500.1
                    mRNA = createItem("MRNA");
                    mRNA.setAttribute("secondaryIdentifier", record.transcriptName);
                    mRNA.setReference("organism", organism);
                    mRNA.setReference("strain", strain);
                    mRNA.addToCollection("dataSets", dataSet);
                    mRNA.setAttribute("assemblyVersion", assemblyVersion);
                    mRNA.setAttribute("annotationVersion", annotationVersion);
                    mRNA.setReference("gene", gene);
                    mRNA.setReference("protein", protein);
                    mRNAs.put(record.transcriptName, mRNA);
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
                        Item goAnnotation = createItem("OntologyAnnotation");
                        goAnnotation.setReference("subject", gene);
                        goAnnotation.setReference("ontologyTerm", goTerm);
                        goAnnotation.addToCollection("dataSets", dataSet);
                        store(goAnnotation);
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
                        Item pfamAnnotation = createItem("OntologyAnnotation");
                        pfamAnnotation.setReference("subject", protein);
                        pfamAnnotation.setReference("ontologyTerm", pfamTerm);
                        pfamAnnotation.addToCollection("dataSets", dataSet);
                        store(pfamAnnotation);
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
                        Item pantherAnnotation = createItem("OntologyAnnotation");
                        pantherAnnotation.setReference("subject", protein);
                        pantherAnnotation.setReference("ontologyTerm", pantherTerm);
                        pantherAnnotation.addToCollection("dataSets", dataSet);
                        store(pantherAnnotation);
                    }
                }
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
                        Item kogAnnotation = createItem("OntologyAnnotation");
                        kogAnnotation.setReference("subject", protein);
                        kogAnnotation.setReference("ontologyTerm", kogTerm);
                        kogAnnotation.addToCollection("dataSets", dataSet);
                        store(kogAnnotation);
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
                        Item ecAnnotation = createItem("OntologyAnnotation");
                        ecAnnotation.setReference("subject", protein);
                        ecAnnotation.setReference("ontologyTerm", ecTerm); 
                        ecAnnotation.addToCollection("dataSets", dataSet);
                       store(ecAnnotation);
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
     * Process an info_annot_ahrd.tsv file which contains gene families and semi-colon separated groups of ontology terms.
     * 0      1       2    3    4               5
     * legume.genefam.fam1.M65K.info_annot_ahrd.tsv
     *
     * legfed_v1_0.L_H6Q7T0-consensus type I inositol-1,4,5-trisphosphate 5-phosphatase;
     *                                IPR005135 (Endonuclease/exonuclease/phosphatase), IPR015943 (WD40/YVTN repeat-like-containing domain);
     *                                GO:0005515 (protein binding), GO:0046856 (phosphatidylinositol dephosphorylation)
     */
    void processInfoAnnotAhrdFile(Reader reader) throws IOException, ObjectStoreException {
        // get the dataSet version from the filename
        String[] dotparts = getCurrentFile().getName().split("\\.");
        String version = dotparts[2];
        Item dataSet = getDataSet();
        dataSet.setAttribute("version", version);
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
        }
        br.close();
    }

    /**
     * Process a single family_fasta file which contains the protein associations with a particular gene family,
     * along with the sequence match, which we do not store.
     *
     * The gene family is given by the name of the file. Only process the organisms specified in project.xml.
     *                                        N-2                                   N-1
     * Gene_families/legume.genefam.fam1.M65K/legume.genefam.fam1.M65K.family_fasta/L_1KSDD1
     *
     * >lupan.Lup015831.1
     * ----CASFAKLT--TLSPHWIGNNSFSSRRGGSSPLTATRRVSLPIRASSYSDELVQTAK
     * TIASPGRGILAIDESNATCGKRLASIGLDNTEVNRQAYRQLLLTTPGLGEYISGAILFEE
     * ...
     * >phavu.Phvul.007G033800.1
     * -----------------------------------TFSPRRVSLPIRASSYQQELVQTAK
     * SIASPGRGILAIDESNATCGKRLASIGLDNTEVNRQAYRQLLLTTPGLGEYISGAILFEE
     * ...
     */
    void processFamilyFastaFile(Reader reader) throws IOException, ObjectStoreException {
        // the DataSet is given by the parent directory
        String[] slashParts = getCurrentFile().getPath().split("/");
        String geneFamilyIdentifier = slashParts[slashParts.length-1];
        String dataSetName = slashParts[slashParts.length-2];
        String[] dotParts = dataSetName.split("\\."); // legume.genefam.fam1.M65K.family_fasta
        String dataSetVersion = dotParts[2];

        Item dataSet = getDataSet(dataSetName);
        dataSet.setAttribute("version", dataSetVersion);
        
        // the gene family is the file name
        Item geneFamily;
        if (geneFamilies.containsKey(getCurrentFile().getName())) {
            geneFamily = geneFamilies.get(getCurrentFile().getName());
        } else {
            geneFamily = createItem("GeneFamily");
            geneFamily.setAttribute("identifier", getCurrentFile().getName());
            geneFamilies.put(getCurrentFile().getName(), geneFamily);
        }

        // proteins are in FASTA headers
        BufferedReader br = new BufferedReader(reader);
        String line = null;
        while ((line=br.readLine())!=null) {
            if (line.startsWith(">")) {
                String genspName = line.substring(1);
                String[] parts = genspName.split("\\.");
                String gensp = parts[0];
                String name = parts[1];
                for (int i=2; i<parts.length; i++) name += "."+parts[i];
                String taxonId = genspTaxonIds.get(gensp);
                if (taxonId!=null && taxonIds.contains(taxonId)) {
                    Item protein = null;
                    if (proteins.containsKey(name)) {
                        protein = proteins.get(name);
                    } else {
                        protein = createItem("Protein");
                        protein.setAttribute("secondaryIdentifier", name);
                        protein.setReference("geneFamily", geneFamily);
                        proteins.put(name, protein);
                    }
                }
            }
        }
        br.close();
    }
}
