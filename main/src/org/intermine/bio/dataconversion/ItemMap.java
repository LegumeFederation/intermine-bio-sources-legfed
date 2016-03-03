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

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.intermine.xml.full.Item;
import org.intermine.xml.full.Reference;
import org.intermine.xml.full.ReferenceList;

/**
 * Extends HashMap to add some handy methods for finding specific Items given (hopefully unique) values.
 * Key is Item refId, value is Item itself.
 *
 * @author Sam Hokin
 */
public class ItemMap extends HashMap<String,Item> {

    /**
     * Put an Item into this Map using its identifier as key
     */
    public void put(Item item) {
        put(item.getIdentifier(), item);
    }

    /**
     * Get an Item from this map by its primaryIdentifier; return null if not found.
     */
    public Item getByPrimaryIdentifier(String primaryIdentifier) {
        for (Item item : this.values()) {
            String primaryId = item.getAttribute("primaryIdentifier").getValue();
            if (primaryId.equals(primaryIdentifier)) return item;
        }
        return null;
    }

    /**
     * Get an Item from this map by its secondaryIdentifier; return null if not found.
     */
    public Item getBySecondaryIdentifier(String secondaryIdentifier) {
        for (Item item : this.values()) {
            String secondaryId = item.getAttribute("secondaryIdentifier").getValue();
            if (secondaryId.equals(secondaryIdentifier)) return item;
        }
        return null;
    }
    
    /**
     * Get an Item from this map by its chadoFeatureId; return null if not found.
     */
    public Item getByChadoFeatureId(int feature_id) {
        for (Item item : this.values()) {
            int chadoFeatureId = Integer.parseInt(item.getAttribute("chadoFeatureId").getValue());
            if (chadoFeatureId==feature_id) return item;
        }
        return null;
    }

    /**
     * Return true if this map contains an Item with the given primary identifier
     */
    public boolean containsPrimaryIdentifier(String primaryIdentifier) {
        for (Item item : this.values()) {
            String primaryId = item.getAttribute("primaryIdentifier").getValue();
            if (primaryId.equals(primaryIdentifier)) return true;
        }
        return false;
    }

    /**
     * Return true if this map contains an Item with the given secondary identifier
     */
    public boolean containsSecondaryIdentifier(String secondaryIdentifier) {
        for (Item item : this.values()) {
            String secondaryId = item.getAttribute("secondaryIdentifier").getValue();
            if (secondaryId.equals(secondaryIdentifier)) return true;
        }
        return false;
    }

    /**
     * Return a List of Items from this map corresponding to the named collection of the given Item.
     * Returns an empty List if no items are found.
     */
    public List<Item> getForCollection(Item item, String collectionName) {
        ReferenceList refList = item.getCollection(collectionName);
        List<Item> itemList = new ArrayList<Item>();
        if (refList!=null) {
            for (String refId : refList.getRefIds()) {
                Item thisItem = this.get(refId);
                if (thisItem!=null) itemList.add(thisItem);
            }
        }
        return itemList;
    }

    /**
     * Return an Item from this map corresponding to the named reference of the given Item.
     * Returns null if none is found.
     */
    public Item getForReference(Item item, String referenceName) {
        Item refItem = null;
        Reference ref = item.getReference(referenceName);
        if (ref!=null) {
            String refId = ref.getRefId();
            refItem = this.get(refId);
        }
        return refItem;
    }

    
}
