package org.intermine.bio.dataconversion;

/*
 * Copyright (C) 2002-2016 FlyMine, Legume Federation
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
import java.util.List;
import java.util.ArrayList;
import java.util.Map;
import java.util.HashMap;

import org.apache.log4j.Logger;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Item;

/**
 * DataConverter to create ExpressionSource, ExpressionSample and ExpressionValue items from a single set of datastore expression files.
 * ICPL87119.gnm1.ann1.expr.KEY4
 * ├── cajca.ICPL87119.gnm1.ann1.KEY4.exprSamples.tsv
 * ├── cajca.ICPL87119.gnm1.ann1.KEY4.exprSource.tsv
 * └── cajca.ICPL87119.gnm1.ann1.KEY4.genesSamplesTpm.tsv
 *
 * cajca.ICPL87119.gnm1.ann1.KEY4.exprSource.tsv
 * ---------------------------------------------
 * #DATASOURCE	
 * ##	
 * ##Name, descriptive (max 250 char)	
 * NAME	Gene expression atlas of pigeonpea Asha(ICPL87119)
 * #Shortname: max 64 char	
 * SHORTNAME	Pigeonpea gene expression atlas
 * #Origin: ex. SRA, GEO, A labname, etc.	
 * ORIGIN	SRA
 * #Description: our description if needed	
 * DESCRIPTION	To be able to link the genome sequence information to the phenotype, especially...
 * #NCBI BioProj acc, PRJNA num	
 * BIOPROJ_ACC	PRJNA354681
 * #NCBI BioProj Title	
 * BIOPROJ_TITLE	Gene expression atlas of pigeonpea
 * #NCBI BioProj description/ abstract(?)	
 * BIOPROJ_DESCRIPTION	Pigeonpea (Cajanus cajan) is an important grain legume of the semi-arid tropics, mainly used...
 * #NCBI SRP number	
 * SRA_PROJ_ACC	SRP097728
 * # GEO series if exists	
 * GEO_SERIES	
 * #Associated publication at NCBI	
 * BIOPROJ_PUBLICATION	Pazhamala, L. T., Purohit, S., Saxena, R. K., Garg, V., Krishnamurthy, L., Verdier, J., & Varshney, R. K. (2017). Gene expression...
 * #Link to Pub	
 * PUB_LINK	https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5429002/
 * #Link to full publication	
 * PUB_FULLLINK	https://academic.oup.com/jxb/article/68/8/2037/3051749/Gene-expression-atlas-of-pigeonpea-and-its
 
 * cajca.ICPL87119.gnm1.ann1.KEY4.exprSamples.tsv
 * ----------------------------------------------
 * #SAMPLE																
 * #																
 * #Hints:																
 * #sample display name	Lookup key for data column (=SRR no)	unique name for LIS Chado	Descriptive	Exp design: factors-rep combinations	tissue/organ/plant part			Chado spelling	sub-spp, cv-type	strain, line, cultivar, genotype	Other_attributes, place-holder, add extra columns	SRR(not in recc)	SAMN number	SRS number	PRJN number	SRP number
 * sample_name	key	sample_uniquename	description	treatment	tissue	dev_stage	age	organism	infraspecies	cultivar	other	sra_run	biosample_accession	sra_accession	bioproject_accession	sra_study
 * Mature seed at reprod (SRR5199304)	SRR5199304	Mature seed at reprod (SRR5199304)	Mature seed at Reproductive stage (SRR5199304)	Mature seed at reprod	Mature seed	Reproductive stage	Cajanus cajan	ICPL87119	Asha(ICPL87119)		SRR5199304	SAMN06264156	SRS1937936	PRJNA354681	SRP097728
 *
 * cajca.ICPL87119.gnm1.ann1.KEY4.genesSamplesTpm.tsv
 * --------------------------------------------------
 * geneID	SRR5199304	SRR5199305	SRR5199306	SRR5199307	...
 * cajca.ICPL87119.gnm1.ann1.C.cajan_00002	0	0	0	0	...
 *
 * @author Sam Hokin
 */
public class ExpressionFileConverter extends BioFileConverter {
    
    private static final Logger LOG = Logger.getLogger(ExpressionFileConverter.class);

    // DataSource is set in project.xml; URL and description are optional since may already exist from other loads.
    Item dataSource;
    String dataSourceName;
    String dataSourceUrl;
    String dataSourceDescription;

    // this particular directory creates a DataSet and ExpressionSource and other stuff for a single experiment
    Item dataSet;
    Item expressionSource;
    Item bioProject;
    Item pub;
    Item organism;
    Item strain;

    // Item maps
    Map<String,Item> samples = new HashMap<>();
    Map<String,Item> genes = new HashMap<>();

    // for non-static utility methods
    DatastoreUtils dsu;

    /**
     * Constructor.
     *
     * @param writer the ItemWriter used to handle the resultant items
     * @param model the Model
     * @throws ObjectStoreException os
     */
    public ExpressionFileConverter(ItemWriter writer, Model model) throws ObjectStoreException {
        super(writer, model);
        // for non-static utility methods
        dsu = new DatastoreUtils();
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

    /**
     * Called for each file found.
     *
     * {@inheritDoc}
     */
    public void process(Reader reader) throws Exception {
        if (dataSource==null) {
            // set defaults for LIS if not given
            if (dataSourceName==null) dataSourceName = DatastoreUtils.DEFAULT_DATASOURCE_NAME;
            if (dataSourceName.equals(DatastoreUtils.DEFAULT_DATASOURCE_NAME)) {
                dataSourceUrl = DatastoreUtils.DEFAULT_DATASOURCE_URL;
                dataSourceDescription = DatastoreUtils.DEFAULT_DATASOURCE_DESCRIPTION;
            }
            // create the DataSource once
            dataSource = createItem("DataSource");
            dataSource.setAttribute("name", dataSourceName);
            if (dataSourceUrl!=null) dataSource.setAttribute("url", dataSourceUrl);
            if (dataSourceDescription!=null) dataSource.setAttribute("description", dataSourceDescription);
        }
        if (getCurrentFile().getName().contains("README")) {
            return;
        } else if (getCurrentFile().getName().endsWith("exprSource.tsv")) {
            // cajca.ICPL87119.gnm1.ann1.KEY4.exprSource.tsv
            if (organism==null) createOrganismAndStrain();
            processSource(reader);
        } else if (getCurrentFile().getName().endsWith("exprSamples.tsv")) {
            // cajca.ICPL87119.gnm1.ann1.KEY4.exprSamples.tsv
            if (organism==null) createOrganismAndStrain();
            processSamples(reader);
        } else if (getCurrentFile().getName().endsWith("genesSamplesTpm.tsv")) {
            // cajca.ICPL87119.gnm1.ann1.KEY4.genesSamplesTpm.tsv
            if (organism==null) createOrganismAndStrain();
            processExpression(reader);
        }
    }

    /**
     * Create the organism and strain Items from the current filename.
     */
    void createOrganismAndStrain() {
        // get the organism and strain from this filename
        String[] pieces = getCurrentFile().getName().split("\\.");
        String gensp = pieces[0];
        String strainId = pieces[1];
        String genomeVersion = pieces[2];
        String annotationVersion = pieces[3];
        organism = createItem("Organism");
        organism.setAttribute("taxonId", dsu.getTaxonId(gensp));
        organism.setAttribute("genus", dsu.getGenus(gensp));
        organism.setAttribute("species", dsu.getSpecies(gensp));
        strain = createItem("Strain");
        strain.setAttribute("identifier", strainId);
        strain.setReference("organism", organism);
    }

    /**
     * Process the datasource meta data file and put the info into a DataSet.
     *
     * cajca.ICPL87119.gnm1.ann1.KEY4.exprSource.tsv
     *
     * NAME	Gene expression atlas of pigeonpea Asha(ICPL87119)
     * SHORTNAME	Pigeonpea gene expression atlas
     * ORIGIN	SRA
     * DESCRIPTION	To be able to link the genome sequence information to the phenotype, especially...
     * BIOPROJ_ACC	PRJNA354681
     * BIOPROJ_TITLE	Gene expression atlas of pigeonpea
     * BIOPROJ_DESCRIPTION	Pigeonpea (Cajanus cajan) is an important grain legume of the semi-arid tropics, mainly used...
     * SRA_PROJ_ACC	SRP097728
     * GEO_SERIES	
     * BIOPROJ_PUBLICATION	Pazhamala, L. T., Purohit, S., Saxena, R. K., Garg, V., Krishnamurthy, L., Verdier, J., & Varshney, R. K. (2017). Gene expression...
     * PUB_LINK	https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5429002/
     * PUB_FULLLINK	https://academic.oup.com/jxb/article/68/8/2037/3051749/Gene-expression-atlas-of-pigeonpea-and-its
     */
    void processSource(Reader reader) throws IOException, ObjectStoreException {
        dataSet = createItem("DataSet");
        if (dataSource!=null) dataSet.setReference("dataSource", dataSource);
        expressionSource = createItem("ExpressionSource");
        expressionSource.setAttribute("unit", "TPM"); // NOTE: assume TPM
        expressionSource.setReference("dataSet", dataSet);
        // we'll favor the pubLink, which is likely a free PMC article
        String pubLink = null;
        String pubFullLink = null;
        // create BioProject only if BIOPROJ lines exist
        bioProject = null;
        // spin through the file
        BufferedReader br = new BufferedReader(reader);
        String line = null;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) continue; // comment or blank
            String[] parts = line.split("\t");
            if (parts.length==2) {
                switch(parts[0]) {
                case "NAME" :
                    dataSet.setAttribute("name", parts[1]);
                    expressionSource.setAttribute("primaryIdentifier", parts[1]);
                    break;
                case "SHORTNAME" :
                    dataSet.setAttribute("shortName", parts[1]);
                    break;
                case "ORIGIN" :
                    dataSet.setAttribute("origin", parts[1]);
                    break;
                case "DESCRIPTION" :
                    dataSet.setAttribute("description", parts[1]);
                    break;
                case "GEO_SERIES" :
                    dataSet.setAttribute("geoSeries", parts[1]);
                    break;
                case "SRA_PROJ_ACC" :
                    dataSet.setAttribute("sra", parts[1]);
                    break;
                case "PUB_LINK" :
                    pubLink = parts[1];
                    break;
                case "PUB_FULLLINK" :
                    pubFullLink = parts[1];
                    break;
                case "BIOPROJ_ACC" :
                    if (bioProject==null) bioProject = createItem("BioProject");
                    bioProject.setAttribute("accession", parts[1]);
                    break;
                case "BIOPROJ_TITLE" :
                    if (bioProject==null) bioProject = createItem("BioProject");
                    bioProject.setAttribute("title", parts[1]);
                    break;
                case "BIOPROJ_DESCRIPTION" :
                    if (bioProject==null) bioProject = createItem("BioProject");
                    bioProject.setAttribute("description", parts[1]);
                    break;
                case "BIOPROJ_PUBLICATION" :
                    // do nothing
                    break;
                default :
                    // log existence of unknown field
                    LOG.info(getCurrentFile().getName()+" contains unknown field:"+line);
                }
            }
        }
        if (pubLink!=null || pubFullLink!=null) {
            Item pub = createItem("Publication");
            if (pubLink!=null) {
                pub.setAttribute("url", pubLink);
            } else if (pubFullLink!=null) {
                pub.setAttribute("url", pubFullLink);
            }
            dataSet.setReference("publication", pub);
            store(pub);
        }
        if (bioProject!=null) {
            dataSet.setReference("bioProject", bioProject);
        }
    }

    /**
     * Process the file which describes the samples.
     *
     * cajca.ICPL87119.gnm1.ann1.KEY4.exprSamples.tsv
     * 
     * sample_name                        key        sample_uniquename                  description                                    treatment             tissue
     * Mature seed at reprod (SRR5199304) SRR5199304 Mature seed at reprod (SRR5199304) Mature seed at Reproductive stage (SRR5199304) Mature seed at reprod Mature seed 
     *
     * dev_stage          age      organism      infraspecies cultivar        other sra_run    biosample_accession sra_accession bioproject_accession sra_study
     * Reproductive stage          Cajanus cajan ICPL87119    Asha(ICPL87119)       SRR5199304 SAMN06264156        SRS1937936    PRJNA354681          SRP097728
     */
    void processSamples(Reader reader) throws IOException {
        String[] colnames = null;
        int num = 0;
        BufferedReader br = new BufferedReader(reader);
        String line = null;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) continue; // comment or blank
            String[] parts = line.split("\t");
            // header contains colnames
            if (line.startsWith("sample_name")) {
                colnames = parts;
            } else if (colnames!=null) {
                if (parts.length!=colnames.length) {
                    System.err.println(line);
                    throw new RuntimeException("ERROR: colnames.length="+colnames.length+" but line has "+parts.length+" fields!");
                }
                String key = null;
                Item sample = createItem("ExpressionSample");
                num++;
                sample.setAttribute("num", String.valueOf(num));
                for (int i=0; i<colnames.length; i++) {
                    switch(colnames[i]) {
                    case "sample_name" :
                        sample.setAttribute("name", parts[i]);
                        break;
                    case "key" :
                        key = parts[i];
                        sample.setAttribute("primaryIdentifier", parts[i]);
                        break;
                    case "sample_uniquename" :
                        break;
                    case "description" :
                        sample.setAttribute("description", parts[i]);
                        break;
                    case "treatment" :
                        break;
                    case "tissue" :
                        break;
                    case "dev_stage" :
                        break;
                    case "age" :
                        break;
                    case "organism" :
                        // we use the filename for organism
                        break;
                    case "infraspecies" :
                        // we use the filename for strain
                        break;
                    case "cultivar" :
                        // we use the filename for strain
                        break;
                    case "other" :
                        break;
                    case "sra_run" :
                        break;
                    case "biosample_accession" :
                        // bioSample accession is all we need, everything else is on NCBI
                        sample.setAttribute("bioSample", parts[i]);
                        break;
                    case "sra_accession" :
                        break;
                    case "bioproject_accession" :
                        break;
                    case "sra_study" :
                        break;
                    }
                }
                if (key!=null && sample!=null) {
                    samples.put(key, sample);
                }
            }
        }
    }

    /**
     * Process a gene expression file. Each gene-sample entry creates an ExpressionValue.
     *
     * cajca.ICPL87119.gnm1.ann1.KEY4.genesSamplesTpm.tsv
     * --------------------------------------------------
     * geneID	SRR5199304	SRR5199305	SRR5199306	SRR5199307	...
     * cajca.ICPL87119.gnm1.ann1.C.cajan_00002	0	0	0	0	...
     */
    void processExpression(Reader reader) throws IOException, ObjectStoreException {
        List<Item> sampleList = new ArrayList<>();
        BufferedReader br = new BufferedReader(reader);
        String line = null;
        while ((line=br.readLine())!=null) {
            if (line.startsWith("#") || line.trim().length()==0) continue; // comment or blank
            String[] parts = line.split("\t");
            if (parts[0].equals("geneID")) {
                // header line gives samples in order so add to list
                for (int i=1; i<parts.length; i++) {
                    String sampleId = parts[i];
                    if (samples.containsKey(sampleId)) {
                        sampleList.add(samples.get(sampleId));
                    } else {
                        Item sample = createItem("ExpressionSample");
                        sample.setAttribute("primaryIdentifier", sampleId);
                        samples.put(sampleId, sample);
                        sampleList.add(sample);
                    }
                }
            } else {
                // a gene line
                String geneId = parts[0];
                Item gene = createItem("Gene");
                gene.setAttribute("primaryIdentifier", geneId);
                genes.put(geneId, gene);
                for (int i=1; i<parts.length; i++) {
                    double value = Double.parseDouble(parts[i]);
                    Item sample = sampleList.get(i-1);
                    Item expressionValue = createItem("ExpressionValue");
                    expressionValue.setAttribute("value", String.valueOf(value));
                    expressionValue.setReference("gene", gene);
                    expressionValue.setReference("sample", sample);
                    store(expressionValue); // this is a one-off so store right away
                }
            }
        }
    }

    /**
     * Store all Items.
     */
    public void close() throws ObjectStoreException {
        // associate the samples with the source and dataSet and organism and strain
        for (Item sample : samples.values()) {
            sample.setReference("source", expressionSource);
            sample.setReference("dataSet", dataSet);
            sample.setReference("organism", organism);
            sample.setReference("strain", strain);
        }
        // store everything
        store(organism);
        store(strain);
        store(dataSource);
        store(dataSet);
        store(bioProject);
        store(expressionSource);
	store(genes.values());
        store(samples.values());
    }
}
