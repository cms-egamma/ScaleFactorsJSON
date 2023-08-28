import correctionlib.schemav2 as schema
from correctionlib.schemav2 import CorrectionSet,Correction,Category,CategoryItem
from JSONTools import *
import pandas as pd
import numpy as np
class MakePhotonJSONs:
    def __init__(self,years):
        print('Initialized')
        self.corrs={}
        for year_i in years:
            self.corrs[year_i]=[]
        self.addedID=False
        self.addedHLT=False
        self.addedCSEV=False
        self.addedPixVeto=False
    def ID(self,ID_dict):
        if self.addedID==True:
            print("Added ID already??? Will not add another. Skipping.")
            return 1
        self.addedID=True
        for year in ID_dict:
            yearname=list(year.keys())[0]
            wps=list(year[yearname].keys())
            #print(year[yearname])
            corr = Correction.parse_obj(
            {
                "version": 2,
                "name": "Photon-ID-SF",
                "description": f"These are the photon ID Scale Factors (nominal, up or down). They are available for the cut-based and MVA IDs. They are dependent on the transverse momenta and pseudorapidity of the photon. ",
                "inputs": [
                    {"name": "year","type": "string", "description": "year/scenario: example 2016preVFP, 2017 etc"},
                    {"name": "ValType","type": "string", "description": "sf/sfup/sfdown (sfup = sf + syst, sfdown = sf - syst) "},
                    {"name": "WorkingPoint","type": "string", "description": "Working Point of choice : Loose, Medium etc."},
                    {"type": "real", "name": "eta", "description": "supercluster eta"},
                    {"name": "pt", "type": "real", "description": "photon pT"},
                ],
                "output": {"name": "weight", "type": "real", "description": "value of scale factor (nominal, up or down)"},
                "data": schema.Category.parse_obj({
                    "nodetype": "category",
                    "input": "year",
                    "content": [
                        schema.CategoryItem.parse_obj({"key":yearname,
                                                       "value": SFyearwise(files=year[yearname],names=wps)})]
                })
            })
            self.corrs[yearname].append(corr)
            
    def CSEV(self,CSEV_dict):
        if self.addedCSEV==True:
            print("Added CSEV already??? Will not add another. Skipping.")
            return 1
        self.addedCSEV=True
        for year in CSEV_dict:
            yearname=list(year.keys())[0]
            wps=list(year[yearname].keys())
            #print(year[yearname])
            corr = Correction.parse_obj(
            {
                "version": 2,
                "name": "Photon-CSEV-SF",
                "description": f"These are the Photon CSEV Scale Factors (nominal, up or down). They are available for the cut-based and MVA IDs. For each working point of choice, they are dependent on the photon pseudorapidity and R9: Possible bin choices: ['EBInc','EBHighR9','EBLowR9','EEInc','EEHighR9','EELowR9']",
                "inputs": [
                    {"name": "year","type": "string", "description": "year/scenario: example 2016preVFP, 2017 etc"},
                    {"name": "ValType","type": "string", "description": "sf/sfup/sfdown (sfup = sf + syst, sfdown = sf - syst) "},
                    {"name": "WorkingPoint","type": "string", "description": "Working Point of choice : Loose, Medium etc."},
                    {"type": "string", "name": "CSEVBin", "description": "Possible bin choices ['EBInc','EBHighR9','EBLowR9','EEInc','EEHighR9','EELowR9']"},
            ],
                "output": {"name": "weight", "type": "real", "description": "value of scale factor (nominal, up or down)"},
                "data": schema.Category.parse_obj({
                    "nodetype": "category",
                    "input": "year",
                    "content": [
                        schema.CategoryItem.parse_obj({"key":yearname,
                                                       "value": CSEVSFyearwise(files=year[yearname],names=wps)})]
                })
            })
            self.corrs[yearname].append(corr)
            
    def PixVeto(self,PixVeto_dict):
        if self.addedPixVeto==True:
            print("Added PixVeto already??? Will not add another. Skipping.")
            return 1
        self.addedPixVeto=True
        for year in PixVeto_dict:
            yearname=list(year.keys())[0]
            wps=list(year[yearname].keys())
            #print(year[yearname])
            corr = Correction.parse_obj(
                {
                    "version": 2,
                    "name": "Photon-PixVeto-SF",
                    "description": f"These are the Photon Pixel Veto Scale Factors (nominal, up or down). They are available for the cut-based and MVA IDs. For each working point of choice, they are dependent on the photon pseudorapidity and R9: Possible bin choices: ['EBInc','EBHighR9','EBLowR9','EEInc','EEHighR9','EELowR9']",
                    "inputs": [
                        {"name": "year","type": "string", "description": "year/scenario: example 2016preVFP, 2017 etc"},
                        {"name": "ValType","type": "string", "description": "sf/sfup/sfdown (sfup = sf + syst, sfdown = sf - syst) "},
                        {"name": "WorkingPoint","type": "string", "description": "Working Point of choice : Loose, Medium etc."},
                        {"type": "string", "name": "HasPixBin", "description": "Possible choices: ['EBInc','EBHighR9','EBLowR9','EEInc','EEHighR9','EELowR9']"},
                    ],
                    "output": {"name": "weight", "type": "real", "description": "value of scale factor (nominal, up or down)"},
                    "data": schema.Category.parse_obj({
                        "nodetype": "category",
                        "input": "year",
                        "content": [
                        schema.CategoryItem.parse_obj({"key":yearname,
                                                       "value": HasPixSFyearwise(files=year[yearname],names=wps)})]
                    })
                })
            self.corrs[yearname].append(corr)
    def HLT(self,HLT_dict):
        if self.addedHLT==True:
            print("Added HLT already??? Will not add another. Skipping.")
            return 1
        self.addedHLT=True
        for year in HLT_dict:
            yearname=list(year.keys())[0]
            wps=list(year[yearname].keys())
            #print(year[yearname])
            corr = Correction.parse_obj(
            {
                "version": 2,
                "name": "Photon-HLT-SF",
                "description": f"These are the photon ID Scale Factors (nominal, up or down). They are available for the cut-based and MVA IDs. They are dependent on the transverse momenta and pseudorapidity of the photon. ",
                "inputs": [
                    {"name": "year","type": "string", "description": "year/scenario: example 2016preVFP, 2017 etc"},
                    {"name": "ValType","type": "string", "description": "sf/sfup/sfdown (sfup = sf + syst, sfdown = sf - syst) "},
                    {"name": "WorkingPoint","type": "string", "description": "Working Point of choice : Loose, Medium etc."},
                    {"type": "real", "name": "eta", "description": "supercluster eta"},
                    {"name": "pt", "type": "real", "description": "photon pT"},
                ],
                "output": {"name": "weight", "type": "real", "description": "value of scale factor (nominal, up or down)"},
                "data": schema.Category.parse_obj({
                    "nodetype": "category",
                    "input": "year",
                    "content": [
                        schema.CategoryItem.parse_obj({"key":yearname,
                                                       "value": SFyearwise(files=year[yearname],names=wps)})]
                })
            })
            self.corrs[yearname].append(corr)
    def SaveJSONs(self,desc,years,era,nameJSON):
        for year in years:
            cset = CorrectionSet(schema_version=2, corrections=self.corrs[year],description=desc)
            with open(year+'_'+era+'/'+nameJSON, "w") as fout:
                fout.write(cset.json(exclude_unset=True, indent=4))

   

