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

    // defaults for LIS datasource
    public static final String DEFAULT_DATASOURCE_NAME = "Legume Information System";
    public static final String DEFAULT_DATASOURCE_URL = "https://legumeinfo.org/";
    public static final String DEFAULT_DATASOURCE_DESCRIPTION =
        "A collaborative, community resource to facilitate crop improvement by integrating genetic, genomic, and trait data across legume species.";
        
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
     * Determine whether the given primaryIdentifier is for a Supercontig (based on some string content).
     */
    static public boolean isSupercontig(String primaryIdentifier) {
        String lc = primaryIdentifier.toLowerCase();
	if (lc.contains("scaffold")
            || lc.contains("contig")
            || lc.contains("pilon")
            || primaryIdentifier.contains("Aipa")
            || primaryIdentifier.contains("Adur")) {
	    return true;
        }
	// tricky ones
	String[] parts = primaryIdentifier.split("\\.");
        if (parts.length>=4) {
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
        }
	// it's a chromosome!
	return false;
    }

    /**
     * Get the taxon ID for a gensp string like "phavu".
     */
    public String getTaxonId(String gensp) {
        String taxonId = genspTaxonId.get(gensp);
        if (taxonId==null) {
            throw new RuntimeException("Taxon ID not available for "+gensp);
        }
        return taxonId;
    }

    /**
     * Get the Genus for a gensp string like "phavu".
     */
    public String getGenus(String gensp) {
        String taxonId = getTaxonId(gensp);
        if (taxonIdGenus.containsKey(taxonId)) {
            return taxonIdGenus.get(taxonId);
        } else {
            throw new RuntimeException("Genus not available for taxon ID "+taxonId);
        }
    }
    
    /**
     * Get the species for a gensp string like "phavu".
     */
    public String getSpecies(String gensp) {
        String taxonId = genspTaxonId.get(gensp);
        if (taxonIdSpecies.containsKey(taxonId)) {
            return taxonIdSpecies.get(taxonId);
        } else {
            throw new RuntimeException("Species not available for taxon ID "+taxonId);
        }
    }
}
