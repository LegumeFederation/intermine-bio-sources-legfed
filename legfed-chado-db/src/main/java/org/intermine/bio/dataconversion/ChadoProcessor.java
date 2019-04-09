package org.intermine.bio.dataconversion;

/*
 * Copyright (C) 2002-2015 FlyMine
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  See the LICENSE file for more
 * information or http://www.gnu.org/copyleft/lesser.html.
 *
 */

import java.util.Collection;
import java.sql.Connection;

import org.intermine.dataconversion.ItemWriter;
import org.intermine.metadata.Model;
import org.intermine.objectstore.ObjectStoreException;
import org.intermine.sql.Database;
import org.intermine.xml.full.Attribute;
import org.intermine.xml.full.Item;
import org.intermine.xml.full.ItemHelper;
import org.intermine.xml.full.Reference;
    
/**
 * A processor for a chado module.  See http://www.gmod.org/wiki/index.php/Chado#Modules for
 * documentation about the possible modules.
 *
 * @author Kim Rutherford
 */
public abstract class ChadoProcessor {
    
    private final ChadoDBConverter chadoDBConverter;

    /**
     * Create a new ChadoModuleProcessor object.
     * @param chadoDBConverter the converter that created this Processor
     */
    public ChadoProcessor(ChadoDBConverter chadoDBConverter) {
        this.chadoDBConverter = chadoDBConverter;
    }

    /**
     * Return the ChadoDBConverter that was passed to the constructor.
     * @return the chadoDBConverter
     */
    public ChadoDBConverter getChadoDBConverter() {
        return chadoDBConverter;
    }

    /**
     * Return the database to read from
     * @return the database
     */
    public Database getDatabase() {
        return chadoDBConverter.getDatabase();
    }

    /**
     * Return an ItemWriter used to handle the resultant Items
     * @return the writer
     */
    public ItemWriter getItemWriter() {
        return chadoDBConverter.getItemWriter();
    }

    /**
     * Return the Model of the target database.
     * @return the model
     */
    public Model getModel() {
        return chadoDBConverter.getModel();
    }

    /**
     * Do the processing for this module - called by ChadoDBConverter.
     * @param connection the database connection to chado
     * @throws Exception if the is a problem while processing
     */
    public abstract void process(Connection connection) throws Exception;

    /**
     * Set an attribute in an Item by creating an Attribute object and storing it.
     * @param intermineObjectId the intermine object ID of the item to create this attribute for.
     * @param attributeName the attribute name
     * @param value the value to set
     * @throws ObjectStoreException if there is a problem while storing
     */
    protected void setAttribute(Integer intermineObjectId, String attributeName, String value)
        throws ObjectStoreException {
        Attribute att = new Attribute();
        att.setName(attributeName);
        att.setValue(value);
        chadoDBConverter.store(att, intermineObjectId);
    }

    /**
     * Store an item.
     * @param item the Item
     * @return the database id of the new Item
     * @throws ObjectStoreException if an error occurs while storing
     */
    protected Integer store(Item item) throws ObjectStoreException {
        return chadoDBConverter.store(item);
    }

    /**
     * Store a single XML Reference
     * @param reference the Reference to store
     * @param itemId the InterMine ID of the Item that holds this Reference
     * @throws ObjectStoreException if an error occurs in storing
     */
    public void store(Reference reference, Integer itemId) throws ObjectStoreException {
	chadoDBConverter.store(reference, itemId);
    }

    /**
     * Store a single org.intermine.model.fulldata.ReferenceList
     * @param referenceList the list to store
     * @param itemId the InterMine ID of the Item that holds this ReferenceList
     * @throws ObjectStoreException if an error occurs in storing
     */
    public void store(org.intermine.model.fulldata.ReferenceList referenceList, Integer itemId) throws ObjectStoreException {
        chadoDBConverter.getItemWriter().store(referenceList, itemId);
    }

    /**
     * Store a single org.intermine.xml.full.ReferenceList
     * @param referenceList the list to store
     * @param itemId the InterMine ID of the Item that holds this ReferenceList
     * @throws ObjectStoreException if an error occurs in storing
     */
    public void store(org.intermine.xml.full.ReferenceList referenceList, Integer itemId) throws ObjectStoreException {
        chadoDBConverter.store(referenceList, itemId);
    }

    /**
     * Store the items.
     * @param items the Item Collection
     * @throws ObjectStoreException if an error occurs while storing
     */
    protected void store(Collection<Item> items) throws ObjectStoreException {
        for (Item item : items) store(item);
    }

    /**
     * Create an Item.
     * @param className the class of the Item to create
     * @return the Item
     * @throws ObjectStoreException
     */
    protected Item createItem(String className) {
        return chadoDBConverter.createItem(className);
    }
}
