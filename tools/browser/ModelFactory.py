# (C) Copyright 2009-2016 ECMWF.
# 
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
# In applying this licence, ECMWF does not waive the privileges and immunities 
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.


#===============================================================================
# ModelFactory
# A class responsible for the instantiation of objects from the right OOPS's model
#===============================================================================
class ModelFactory:
    __factories = {}

    #===========================================================================
    # createPlot
    # Instatiate a new plot
    #===========================================================================
    def createPlot(id):
        return ModelFactory.__factories[id].createPlot()

    #===========================================================================
    # createFdb
    # Instantiate a new states manager
    #===========================================================================
    def createFdb(id, dir):
        return ModelFactory.__factories[id].createFdb(dir)

    #===========================================================================
    # createOdb
    # Instantiate a new Odb
    #===========================================================================
    def createOdb(id, file):
        return ModelFactory.__factories[id].createOdb(file)
    
    #===========================================================================
    # title
    # the title
    #===========================================================================
    def title(id):
        return ModelFactory.__factories[id].title()

    #===========================================================================
    # addFactory
    # Adds a new factory in the builder
    #===========================================================================
    def addFactory(factory):
        ModelFactory.__factories[factory.id()] = factory

    #===========================================================================
    # List of static methods
    #===========================================================================
    addFactory = staticmethod(addFactory)
    createPlot = staticmethod(createPlot)
    createFdb = staticmethod(createFdb)
    createOdb = staticmethod(createOdb)
    title = staticmethod(title)
