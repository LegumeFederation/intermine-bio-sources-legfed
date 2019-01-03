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

import org.intermine.bio.io.gff3.GFF3Record;
import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.xml.full.Item;

/**
 * Read synteny blocks for two organisms from a DAGchainer synteny GFF file and store them as SyntenyBlock items,
 * each related to two SyntenicRegions. This is designed to use the GFF annotation produced by DAGchainer.
 *
 * The DAGchainer GFF file lines must look like this:
 *
 * phavu.Chr01 DAGchainer syntenic_region 125452 912158 2665.5 - . Name=Pv01.Gm14.2.+;Parent=19;ID=20;Target=glyma.Chr14:48062215..48932270;median_Ks=0.3559
 *
 * That is:
 * - the target sequence must be given by a Target attribute.
 * - the Name attribute must be a unique identifier.
 * - the median_Ks attribute must be spelled that way.
 * - the source ID and target ID must match the primaryIdentifier of chromosomes in the production database.
 * - the record type must be "syntenic_region."
 * - the target strand, + or -, may optionally be given in the last character of the Name attribute, as in this example. Otherwise, strand isn't recorded for the target.
 * The Parent and ID attributes are not used.
 *
 * In addition, the source and target organisms/strains are identified by comments at the TOP of the file, such as:
 *
 * #SourceTaxonID  130453
 * #SourceStrain  V14167
 * #TargetTaxonID  130454
 * #TargetStrain  K30076
 *
 * @author Sam Hokin, NCGR
 */
public class SyntenyGFFConverter extends BioFileConverter {
	
    private static final Logger LOG = Logger.getLogger(SyntenyGFFConverter.class);

    static final String RANGE_SEPARATOR = "-";
    static final String REGION_SEPARATOR = "|";

    // these maps prevent duplicate stores
    Map<String,Item> organismMap = new HashMap<String,Item>();
    Map<String,Item> strainMap = new HashMap<String,Item>();
    Map<String,Item> chromosomeMap = new HashMap<String,Item>();

    // use this map to prevent storing duplicate synteny blocks with regions swapped
    Map<String,String> syntenyBlockMap = new HashMap<String,String>();
        
    /**
     * Create a new SyntenyGFFConverter
     * @param writer the ItemWriter to write out new items
     * @param model the data model
     */
    public SyntenyGFFConverter(ItemWriter writer, Model model) {
        super(writer, model);
    }

    /**
     * {@inheritDoc}
     * We process each GFF file by creating SyntenyBlock and SyntenicRegion items and storing them.
     */
    @Override
    public void process(Reader reader) throws Exception {

        LOG.info("Processing Synteny file "+getCurrentFile().getName()+"...");

        // constants for each file
        String sourceTaxonId = null;
        String targetTaxonId = null;
        String sourceStrainName = null;
        String targetStrainName = null;
        Item sourceOrganism = null;
        Item targetOrganism = null;
	Item sourceStrain = null;
	Item targetStrain = null;
        
        // -------------------------------------------------------------------------------------------------------
        // Load the GFF data into a map. Add new chromosomes to chromosome map, keyed by primaryIdentifier.
        // -------------------------------------------------------------------------------------------------------

        BufferedReader gffReader = new BufferedReader(reader);
        String line = null;
        while ((line=gffReader.readLine()) != null) {

            // create and store organism and strain Items if info has been loaded
	    if (sourceOrganism==null && sourceTaxonId!=null) {
                if (organismMap.containsKey(sourceTaxonId)) {
                    sourceOrganism = organismMap.get(sourceTaxonId);
                } else {
                    sourceOrganism = createItem("Organism");
                    sourceOrganism.setAttribute("taxonId", sourceTaxonId);
                    store(sourceOrganism);
                    organismMap.put(sourceTaxonId, sourceOrganism);
                    LOG.info("Created source organism: "+sourceTaxonId);
                }
	    }
	    if (targetOrganism==null && targetTaxonId!=null) {
                if (organismMap.containsKey(targetTaxonId)) {
                    targetOrganism = organismMap.get(targetTaxonId);
                } else {
                    targetOrganism = createItem("Organism");
                    targetOrganism.setAttribute("taxonId", targetTaxonId);
                    store(targetOrganism);
                    organismMap.put(targetTaxonId, targetOrganism);
                    LOG.info("Created target organism: "+targetTaxonId);
                }
            }
	    if (sourceStrain==null && sourceStrainName!=null) {
                if (strainMap.containsKey(sourceStrainName)) {
                    sourceStrain = strainMap.get(sourceStrainName);
                } else {
                    sourceStrain = createItem("Strain");
                    sourceStrain.setAttribute("primaryIdentifier", sourceStrainName);
                    store(sourceStrain);
                    strainMap.put(sourceStrainName, sourceStrain);
                    LOG.info("Created source strain: "+sourceStrainName);
                }
	    }
	    
	    if (targetStrain==null && targetStrainName!=null) {
		if (strainMap.containsKey(targetStrainName)) {
		    targetStrain = strainMap.get(targetStrainName);
		} else {
		    targetStrain = createItem("Strain");
		    targetStrain.setAttribute("primaryIdentifier", targetStrainName);
		    store(targetStrain);
		    strainMap.put(targetStrainName, targetStrain);
		    LOG.info("Created target strain: "+targetStrainName);
		}
	    }

	    if (line.startsWith("#SourceTaxonID")) {

		String[] parts = line.split("\t");
		try {
		    sourceTaxonId = parts[1];
		    LOG.info("sourceTaxonId="+sourceTaxonId);
		} catch (Exception e) {
		    System.err.println("INPUT ERROR:"+line);
		    throw new RuntimeException(e);
		}

	    } else if (line.startsWith("#SourceStrain")) {
                
		String[] parts = line.split("\t");
		try {
		    sourceStrainName = parts[1];
		    LOG.info("sourceStrain="+sourceStrainName);
		} catch (Exception e) {
		    System.err.println("INPUT ERROR:"+line);
		    throw new RuntimeException(e);
		}

	    } else if (line.startsWith("#TargetTaxonID")) {

		String[] parts = line.split("\t");
		try {
		    targetTaxonId = parts[1];
		    LOG.info("targetTaxonId="+targetTaxonId);
		} catch (Exception e) {
		    System.err.println("INPUT ERROR:"+line);
		    throw new RuntimeException(e);
		}

	    } else if (line.startsWith("#TargetStrain")) {
                
		String[] parts = line.split("\t");
		try {
		    targetStrainName = parts[1];
		    LOG.info("targetStrain="+targetStrainName);
		} catch (Exception e) {
		    System.err.println("INPUT ERROR:"+line);
		    throw new RuntimeException(e);
		}

	    } else if (line.startsWith("#") || line.trim().length()==0) {
                
		// do nothing, comment or line skip

	    } else {

		// bail if we don't have organisms or strains instantiated
		if (sourceOrganism==null) {
		    LOG.error("Source organism not established:"+sourceTaxonId);
		    throw new RuntimeException("Source organism not established:"+sourceTaxonId);
		}
		if (targetOrganism==null) {
		    LOG.error("Target organism not established:"+targetTaxonId);
		    throw new RuntimeException("Target organism not established:"+targetTaxonId);
		}
		if (sourceStrain==null) {
		    LOG.error("Source strain not established:"+sourceStrainName);
		    throw new RuntimeException("Source strain not established:"+sourceStrainName);
		}
		if (targetStrain==null) {
		    LOG.error("Target strain not established:"+targetStrainName);
		    throw new RuntimeException("Target strain not established:"+targetStrainName);
		}
		
		// load the GFF line
		GFF3Record gff = new GFF3Record(line);
		if (gff.getType().equals("syntenic_region")) {
		    String sourceChrName = getSourceChromosomeName(gff);
		    String targetChrName = getTargetChromosomeName(gff);
		    if (targetChrName==null) {
			throw new RuntimeException("GFF syntenic_region record is missing Target= attribute:"+line);
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
			LOG.info("Created new source chromosome:"+sourceChrName);
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
			LOG.info("Created new target chromosome:"+targetChrName);
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
		    if (syntenyBlockMap.containsKey(sourceIdentifier) && syntenyBlockMap.get(sourceIdentifier).equals(targetIdentifier) ||
			syntenyBlockMap.containsKey(targetIdentifier) && syntenyBlockMap.get(targetIdentifier).equals(sourceIdentifier)) {
                        
			// do nothing
                        
		    } else {

			// store the regions in the map for future non-duplication
			syntenyBlockMap.put(sourceIdentifier, targetIdentifier);
                    
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
     */
    String getTargetChromosomeName(GFF3Record gff) {
	if (gff.getTarget()==null) {
	    return null;
	} else {
	    String[] chunks = gff.getTarget().split(":");
	    return chunks[0];
	}
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
