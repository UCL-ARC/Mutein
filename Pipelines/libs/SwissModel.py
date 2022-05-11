"""
RSA 11/5/22
------------------------
Class to manage interactions and data from Swiss Model, the homology structure tool

"""
import os
import subprocess


class SwissModel:
    def __init__(self, accession):                
        self.accession = accession

    def searchForStructures(self):
        #https://swissmodel.expasy.org/repository/uniprot/P46531.json
        html = "https://swissmodel.expasy.org/repository/uniprot/" + self.accession.upper() + ".json"