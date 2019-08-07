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
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.xml.full.Item;

/**
 * Read synteny blocks for two organisms from a DAGchainer synteny GFF file and store them as SyntenyBlock items,
 * each related to two SyntenicRegions. This is designed to use the GFF annotation produced by DAGchainer.
 *
 * ##gff-version 3
 * ##date Mon May  2 14:03:32 2016
 * ##source gbrowse gbgff gff3 dumper
 * #SourceTaxonID	130453
 * #SourceStrain	V14167
 * #TargetTaxonID	130454
 * #TargetStrain	K30076
 * Aradu.A01	DAGchainer	syntenic_region	63347	127681	173.0	-	.	Name=Araip.B01.1;median_Ks=0.0430;Target=Araip.B01:17125379..17229197
 *
 * Source and target strains are given by the taxon ID and strain names in the header.
 *
 * @author Sam Hokin
 */
public class SyntenyGFFConverter extends BioFileConverter {
	
    private static final Logger LOG = Logger.getLogger(SyntenyGFFConverter.class);

    static final String RANGE_SEPARATOR = "-";
    static final String REGION_SEPARATOR = "|";

    // these maps prevent duplicate stores
    Map<String,Item> strainMap = new HashMap<>();
    Map<String,Item> organismMap = new HashMap<>();
    Map<String,Item> chromosomeMap = new HashMap<>();

    // use this map to prevent storing duplicate synteny blocks with regions swapped
    Map<String,String> syntenyBlocks = new HashMap<>();

    // map gensp to taxon Id
    Map<String,String> genspTaxonId = new HashMap<>();
        
    /**
     * Create a new SyntenyGFFConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public SyntenyGFFConverter(ItemWriter writer, Model model) {
        super(writer, model);

        // get the organism data
        DatastoreUtils dsu = new DatastoreUtils();
        genspTaxonId = dsu.getGenspTaxonId();
    }

    /**
     * {@inheritDoc}
     * We process each GFF file by creating SyntenyBlock and SyntenicRegion items and storing them.
     */
    @Override
    public void process(Reader reader) throws Exception {

        LOG.info("Processing Synteny file "+getCurrentFile().getName()+"...");
        System.out.println("##### Processing Synteny file "+getCurrentFile().getName()+" #####");

        String sourceTaxonId = null;
        String sourceStrainName = null;
        String targetTaxonId = null;
        String targetStrainName = null;

        Item sourceOrganism = null;
	Item sourceStrain = null;
        Item targetOrganism = null;
	Item targetStrain = null;
        
        // -------------------------------------------------------------------------------------------------------
        // Load the GFF data into a map. Add new chromosomes to chromosome map, keyed by primaryIdentifier.
        // -------------------------------------------------------------------------------------------------------

        BufferedReader gffReader = new BufferedReader(reader);
        String line = null;
        while ((line=gffReader.readLine())!=null) {
            if (line.startsWith("#SourceTaxonID")) {
                String[] parts = line.split("\t");
                sourceTaxonId = parts[1];
                sourceOrganism = getOrganismItem(sourceTaxonId);
            } else if (line.startsWith("#SourceStrain")) {
                if (sourceOrganism==null) throw new RuntimeException("#SourceTaxonID must be placed before #SourceStrain in synteny GFF file.");
                String[] parts = line.split("\t");
                sourceStrainName = parts[1];
                sourceStrain = getStrainItem(sourceStrainName, sourceOrganism);
            } else if (line.startsWith("#TargetTaxonID")) {
                String[] parts = line.split("\t");
                targetTaxonId = parts[1];
                targetOrganism = getOrganismItem(targetTaxonId);
            } else if (line.startsWith("#TargetStrain")) {
                if (targetOrganism==null) throw new RuntimeException("#TargetTaxonID must be placed before #TargetStrain in synteny GFF file.");
                String[] parts = line.split("\t");
                targetStrainName = parts[1];
                targetStrain = getStrainItem(targetStrainName, targetOrganism);
            } else if (line.startsWith("#") || line.trim().length()==0) {
                // do nothing, true comment or blank line
            } else {
                // load the GFF line
                GFF3Record gff = new GFF3Record(line);
                if (gff.getType().equals("syntenic_region")) {
                    String sourceChrName = getSourceChromosomeName(gff);
                    String targetChrName = getTargetChromosomeName(gff);
                    if (targetChrName==null) {
                        throw new RuntimeException("GFF syntenic_region record is missing target attribute:"+line);
                    }
                    // ignore this record if it's a scaffold
                    if (targetChrName.toLowerCase().contains("scaffold") || targetChrName.toLowerCase().contains("superscaf")) {
                        continue;
                    }
                    Item sourceChromosome;
                    if (chromosomeMap.containsKey(sourceChrName)) {
                        sourceChromosome = chromosomeMap.get(sourceChrName);
                    } else {
                        // create and store the source chromosome and add to the chromosome map
                        sourceChromosome = createItem("Chromosome");
                        sourceChromosome.setAttribute("primaryIdentifier", sourceChrName);
                        sourceChromosome.setReference("organism", sourceOrganism);
                        sourceChromosome.setReference("strain", sourceStrain);
                        store(sourceChromosome);
                        chromosomeMap.put(sourceChrName, sourceChromosome);
                    }
                    Item targetChromosome;
                    if (chromosomeMap.containsKey(targetChrName)) {
                        targetChromosome = chromosomeMap.get(targetChrName);
                    } else {
                        // create and store the target chromosome and add to the chromosome map
                        targetChromosome = createItem("Chromosome");
                        targetChromosome.setAttribute("primaryIdentifier", targetChrName);
                        targetChromosome.setReference("organism", targetOrganism);
                        targetChromosome.setReference("strain", targetStrain);
                        store(targetChromosome);
                        chromosomeMap.put(targetChrName, targetChromosome);
                    }
                    
                    // populate the source region and its location
                    Item sourceRegion = createItem("SyntenicRegion");
                    Item sourceChromosomeLocation = createItem("Location");
                    populateSourceRegion(sourceRegion, gff, sourceOrganism, sourceStrain, sourceChromosome, sourceChromosomeLocation);
                    String sourceIdentifier = getSourceRegionName(gff);
                    
                    // populate the target region and its location
                    Item targetRegion = createItem("SyntenicRegion");
                    Item targetChromosomeLocation = createItem("Location");
                    populateTargetRegion(targetRegion, gff, targetOrganism, targetStrain, targetChromosome, targetChromosomeLocation);
                    String targetIdentifier = getTargetRegionName(gff);
                    
                    // only continue if we haven't stored a synteny block with these regions
                    if (syntenyBlocks.containsKey(sourceIdentifier) && syntenyBlocks.get(sourceIdentifier).equals(targetIdentifier) ||
                        syntenyBlocks.containsKey(targetIdentifier) && syntenyBlocks.get(targetIdentifier).equals(sourceIdentifier)) {
                        
                        // do nothing
                        
                    } else {
                        
                        // store the regions in the map for future non-duplication
                        syntenyBlocks.put(sourceIdentifier, targetIdentifier);
                        
                        // get the medianKs value for this block
                        Map<String, List<String>> attributes = gff.getAttributes();
                        String medianKs = attributes.get("median_Ks").get(0);
                        
                        // associate the two regions with this synteny block
                        Item syntenyBlock = createItem("SyntenyBlock");
                        syntenyBlock.setAttribute("primaryIdentifier", getSyntenyBlockName(gff));
                        syntenyBlock.setAttribute("medianKs", medianKs);
                        syntenyBlock.addToCollection("syntenicRegions", sourceRegion);
                        syntenyBlock.addToCollection("syntenicRegions", targetRegion);
                        store(syntenyBlock);
                        
                        // associate the block with the regions and store them
                        sourceRegion.setReference("syntenyBlock", syntenyBlock);
                        store(sourceRegion);
                        store(sourceChromosomeLocation);
                        targetRegion.setReference("syntenyBlock", syntenyBlock);
                        store(targetRegion);
                        store(targetChromosomeLocation);
                    }
                }
            }
        }
	gffReader.close();
    }

    /**
     * Get/add the organism Item associated with the given taxon ID.
     */
    Item getOrganismItem(String taxonId) throws ObjectStoreException {
        Item organism;
        if (organismMap.containsKey(taxonId)) {
            organism = organismMap.get(taxonId);
        } else {
            organism = createItem("Organism");
            organism.setAttribute("taxonId", taxonId);
            store(organism);
            organismMap.put(taxonId, organism);
        }
        return organism;
    }


    /**
     * Get/add the strain Item associated with the given strain name.
     * Sets the organism reference if created.
     */
    Item getStrainItem(String strainName, Item organism) throws ObjectStoreException {
        Item strain;
        if (strainMap.containsKey(strainName)) {
            strain = strainMap.get(strainName);
        } else {
            strain = createItem("Strain");
            strain.setAttribute("primaryIdentifier", strainName);
            strain.setReference("organism", organism);
            store(strain);
            strainMap.put(strainName, strain);
        }
        return strain;
    }

    /**
     * Populate the attributes of a source SyntenicRegion with a GFF3Record's data.
     *
     * @param syntenicRegion the SyntenicRegion Item
     * @param gff the GFF3Record holding the data
     * @param chromosome the source Chromosome Item
     * @param chromosomeLocation the source Location Item to be filled in
     */
    void populateSourceRegion(Item syntenicRegion, GFF3Record gff, Item organism, Item strain, Item chromosome, Item chromosomeLocation) {
	syntenicRegion.setAttribute("primaryIdentifier", getSourceRegionName(gff));
	syntenicRegion.setAttribute("length", String.valueOf(getSourceEnd(gff)-getSourceStart(gff)+1));
	syntenicRegion.setAttribute("score", String.valueOf(gff.getScore()));
	syntenicRegion.setReference("organism", organism);
	syntenicRegion.setReference("strain", strain);
	syntenicRegion.setReference("chromosome", chromosome);
	syntenicRegion.setReference("chromosomeLocation", chromosomeLocation);
	chromosomeLocation.setAttribute("start", String.valueOf(getSourceStart(gff)));
	chromosomeLocation.setAttribute("end", String.valueOf(getSourceEnd(gff)));
	chromosomeLocation.setAttribute("strand", String.valueOf(gff.getStrand()));
	chromosomeLocation.setReference("feature", syntenicRegion);
	chromosomeLocation.setReference("locatedOn", chromosome);
    }

    /**
     * Populate the attributes of a target SyntenicRegion with a GFF3Record's DAGchainer attributes data; Organism, Chromosome and ChromosomeLocation Items must be passed in as well.
     *
     * @param syntenicRegion the SyntenicRegion Item
     * @param gff the GFF3Record holding the data
     * @param chromosome the target Chromosome Item
     * @param chromosomeLocation the target Location Item to be filled in
     */
    void populateTargetRegion(Item syntenicRegion, GFF3Record gff, Item organism, Item strain, Item chromosome, Item chromosomeLocation) {
	syntenicRegion.setAttribute("primaryIdentifier", getTargetRegionName(gff));
	syntenicRegion.setAttribute("length", String.valueOf(getTargetEnd(gff)-getTargetStart(gff)+1));
	syntenicRegion.setAttribute("score", String.valueOf(gff.getScore()));
	syntenicRegion.setReference("organism", organism);
	syntenicRegion.setReference("strain", strain);
	syntenicRegion.setReference("chromosome", chromosome);
	syntenicRegion.setReference("chromosomeLocation", chromosomeLocation);
	chromosomeLocation.setAttribute("start", String.valueOf(getTargetStart(gff)));
	chromosomeLocation.setAttribute("end", String.valueOf(getTargetEnd(gff)));
	if (getTargetStrand(gff)!='\0') chromosomeLocation.setAttribute("strand", String.valueOf(getTargetStrand(gff)));
	chromosomeLocation.setReference("feature", syntenicRegion);
	chromosomeLocation.setReference("locatedOn", chromosome);
    }

    /**
     * Return the DAGchainer target strand from a DAGchainer Name attribute, checking for ' ' instead of '+' since GFF3Record converts a plus to space.
     */
    char getTargetStrand(GFF3Record gff) {
	String name = gff.getNames().get(0);
	char endChar = name.charAt(name.length()-1);
	if (endChar==' ' || endChar=='+') {
	    return '+';
	} else if (endChar=='-') {
	    return '-';
	} else {
	    return '\0';
	}
    }

    /**
     * Return the source chromosome nanme - opportunity for tweaks here
     */
    String getSourceChromosomeName(GFF3Record gff) {
	return gff.getSequenceID();
    }

    /**
     * Return the DAGchainer target chromosome from a DAGchainer GFF3Record
     * Target=Araip.B01:17125379..17229197
     */
    String getTargetChromosomeName(GFF3Record gff) {
        String[] chunks = gff.getTarget().split(":");
        return chunks[0];
    }
    
    /**
     * Return the source sequence start from a GFF3 record - just echoes GFF3Record.getStart() for clarity
     */
    int getSourceStart(GFF3Record gff) {
	return gff.getStart();
    }

    /**
     * Return the source sequence end from a GFF3 record - just echoes GFF3Record.getEnd() for clarity
     */
    int getSourceEnd(GFF3Record gff) {
	return gff.getEnd();
    }

    
    /**
     * Return the target sequence start from a DAGchainer GFF3Record
     * Target=Araip.B01:17125379..17229197
     */
    int getTargetStart(GFF3Record gff) {
	String[] chunks = gff.getTarget().split(":");
	String range = chunks[1];
	String[] pieces = range.split("\\.\\.");
	return Integer.parseInt(pieces[0]);
    }

    /**
     * Return the target sequence end from a DAGchainer GFF3Record
     */
    int getTargetEnd(GFF3Record gff) {
	String[] chunks = gff.getTarget().split(":");
	String range = chunks[1];
	String[] pieces = range.split("\\.\\.");
	return Integer.parseInt(pieces[1]);
    }

    /**
     * Return a source syntenic region primary identifier from a GFF record; same as GBrowse standard
     */
    String getSourceRegionName(GFF3Record gff) {
	String sourceChrName = getSourceChromosomeName(gff);
	int sourceStart = getSourceStart(gff);
	int sourceEnd = getSourceEnd(gff);
	return sourceChrName+":"+sourceStart+RANGE_SEPARATOR+sourceEnd;
    }

    /**
     * Return a target syntenic region primary identifier from a GFF record; same as GBrowse standard
     */
    String getTargetRegionName(GFF3Record gff) {
	String targetChrName = getTargetChromosomeName(gff);
	int targetStart = getTargetStart(gff);
	int targetEnd = getTargetEnd(gff);
	return targetChrName+":"+targetStart+RANGE_SEPARATOR+targetEnd;
    }

    /**
     * Return a synteny block primary identifier formed from the source and target names with a separator
     */
    String getSyntenyBlockName(GFF3Record gff) {
	return getSourceRegionName(gff)+REGION_SEPARATOR+getTargetRegionName(gff);
    }
}
