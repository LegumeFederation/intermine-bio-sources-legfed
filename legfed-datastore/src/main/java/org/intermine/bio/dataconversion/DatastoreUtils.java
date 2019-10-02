package org.intermine.bio.dataconversion;

import java.io.InputStream;
import java.io.IOException;

import java.util.Map;
import java.util.HashMap;
import java.util.Properties;

/**
 * Utility methods for data store.
 *
 * @author Sam Hokin
 */
public class DatastoreUtils {

    private static final String MINE_PROP_FILE = "organism_config.properties";

    private Map<String,String> taxonIdGenus = new HashMap<>();
    private Map<String,String> taxonIdSpecies = new HashMap<>();
    private Map<String,String> genspTaxonId = new HashMap<>();
    private Map<String,String> genusSpeciesTaxonId = new HashMap<>();

    public DatastoreUtils() {
        // map taxonId to genus, species
        // 0     1    2
        // taxon.3702.genus=Arabidopsis
        // taxon.3702.species=thaliana
        // taxon.3702.uniprot=ARATH
        Properties props = new Properties();
        try {
            InputStream propsResource = DatastoreUtils.class.getClassLoader().getResourceAsStream(MINE_PROP_FILE);
            if (propsResource == null) {
                System.err.println("Did not find organism properties file:"+MINE_PROP_FILE);
                System.exit(1);
            }
            props.load(propsResource);
            for (Object obj : props.keySet()) {
                String key = (String) obj;
                String value = props.getProperty(key);
                String[] dotparts = key.split("\\.");
                if (dotparts[0].equals("taxon")) {
                    String taxonId = dotparts[1];
                    if (dotparts[2].equals("genus")) {
                        taxonIdGenus.put(taxonId, value);
                    } else if (dotparts[2].equals("species")) {
                        taxonIdSpecies.put(taxonId, value);
                    }
                }
            }
        } catch (IOException e) {
            throw new RuntimeException("Problem loading properties from:"+MINE_PROP_FILE, e);
        }
        // map gensp, Genus_species to taxonId
        for (String taxonId : taxonIdGenus.keySet()) {
            String genus = taxonIdGenus.get(taxonId);
            String species = taxonIdSpecies.get(taxonId);
            String gensp = genus.substring(0,3).toLowerCase()+species.substring(0,2).toLowerCase();
            String genusSpecies = genus+"_"+species;
            genspTaxonId.put(gensp, taxonId);
            genusSpeciesTaxonId.put(genusSpecies, taxonId);
        }
    }

    /** 
     * get the taxonId -> Genus map
     */
    public Map<String,String> getTaxonIdGenus() {
        return taxonIdGenus;
    }

    /**
     * get the taxonId -> species map
     */
    public Map<String,String> getTaxonIdSpecies() {
        return taxonIdSpecies;
    }

    /**
     * get the gensp -> taxonId map
     */
    public Map<String,String> getGenspTaxonId() {
        return genspTaxonId;
    }

    /**
     * get the genus_species -> taxonId map
     */
    public Map<String,String> getGenusSpeciesTaxonId() {
        return genusSpeciesTaxonId;
    }

    /**
     * Determine whether the given primaryIdentifier is for a Supercontig (based on some string content).
     */
    static public boolean isSupercontig(String primaryIdentifier) {
        String lc = primaryIdentifier.toLowerCase();
	if (lc.contains("scaffold")
            || lc.contains("contig")
            || lc.contains("pilon")
            || primaryIdentifier.contains("Aipa")
            || primaryIdentifier.contains("Adur"))
	    return true;
	// tricky ones
	String[] parts = primaryIdentifier.split("\\.");
	// 0     1           2    3
	// cicar.CDCFrontier.gnm1.C11044140
	if (parts[3].length()==9 && parts[3].charAt(0)=='C') return true;
	// 0     1   2    3
	// glyma.Lee.gnm1.sc119
	if (parts[1].equals("Lee") && parts[3].startsWith("sc")) return true;
	// 0     1        2    3
	// glyso.PI483463.gnm1.sc255
	if (parts[0].equals("glyso") && parts[3].startsWith("sc")) return true;
	// 0     1            2    3
	// medtr.jemalong_A17.gnm5.MtrunA17Chr0c01
	if (parts[3].contains("Chr0c")) return true;
	// it's a chromosome!
	return false;
    }
}
