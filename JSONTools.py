# Technical part below
from ROOT import TFile, TH1F, TCanvas, TString
import correctionlib.schemav2 as schema
from correctionlib.schemav2 import Correction
corrs=[]
#Function to extract SFs from EGamma standard root files
def getSFs(fn="filename",IsSF="sf",sfhist="EGamma_SF2D"):
    tf = TFile(fn)
    fo = tf.Get(sfhist)
    Xbins=fo.GetNbinsX()
    Ybins=fo.GetNbinsY()
    X=[None]*(Xbins+1)
    Y=[None]*(Ybins+1)
    values=[]
    errors=[]
    for i in range(1,Xbins + 1):
        X[i-1]=(fo.GetXaxis().GetBinLowEdge(i))
    X[Xbins]=fo.GetXaxis().GetBinUpEdge(Xbins)
    for j in range(1,Ybins + 1):
        Y[j-1]=(fo.GetYaxis().GetBinLowEdge(j))
    Y[Ybins]=fo.GetYaxis().GetBinUpEdge(Ybins)
    for i in range(1,Xbins + 1):
        for j in range(1,Ybins + 1):
            values.append(fo.GetBinContent(i,j))
            errors.append(fo.GetBinError(i,j))
    if IsSF=="sf":
        valSFs=schema.MultiBinning.parse_obj({
            "inputs":["eta","pt"],
            "nodetype": "multibinning",
            "edges": [
                X,
                Y,
            ],
            "content": values,
            "flow": 'error',
        })
        return valSFs
    if IsSF=="syst":
        valerrors=schema.MultiBinning.parse_obj({
            "inputs":["eta","pt"],
            "nodetype": "multibinning",
            "edges": [
                X,
                Y,
            ],
            "content": errors,
            "flow": 'error',
        })
        return valerrors
    
def SFyearwise(files=[],names=[],valtypes=["sf","syst"]):
    output = schema.Category.parse_obj({
                "nodetype": "category",
                "input": "ValType",
                "content":[
                    schema.CategoryItem.parse_obj({
                        "key": val, 
                        "value": schema.Category.parse_obj({
                            "nodetype": "category",
                            "input": "WorkingPoint",
                            "content":[
                                schema.CategoryItem.parse_obj({
                                        "key": name, 
                                        "value": getSFs(fn=files[name],IsSF=val)})
                                for name in names
                                ],
                            })
                    })
                    for val in valtypes
                ],
    })
    return output


def CSEVSFs(files,name,i,IsSF="sf"):
    file=TFile(files[name])
    hist=file.Get(name+"ID/SF_CSEV_"+name+"ID")
    SF=0
    if IsSF=="sf":
        SF=hist.GetBinContent(i)
    else:
        SF=hist.GetBinError(i)
    return SF

def CSEVSFyearwise(files=[],names=[],valtypes=["sf","syst"]):
    binlist=['EBInc','EBHighR9','EBLowR9','EEInc','EEHighR9','EELowR9']
    output = schema.Category.parse_obj({
                "nodetype": "category",
                "input": "ValType",
                "content":[
                    schema.CategoryItem.parse_obj({
                        "key": val, 
                        "value": schema.Category.parse_obj({
                            "nodetype": "category",
                            "input": "WorkingPoint",
                            "content":[
                                schema.CategoryItem.parse_obj({
                                        "key": name, 
                                        "value": schema.Category.parse_obj({
                                            "nodetype": "category",
                                            "input": "CSEVBin",
                                            "content": [schema.CategoryItem.parse_obj({"key": binlist[i-1], 
                                                                                       "value": CSEVSFs(files,name,i,val)})
                                                        for i in range(1,7)
                                                       ],
                                        })
                                })
                                for name in names
                                ],
                            })
                    })
                    for val in valtypes
                ],
    })
    return output
def HasPixSFs(files,name,i,IsSF="sf"):
    file=TFile(files[name])
    hist=file.Get(name+"ID/SF_HasPix_"+name+"ID")
    SF=0
    if IsSF=="sf":
        SF=hist.GetBinContent(i)
    else:
        SF=hist.GetBinError(i)
    return SF

def HasPixSFyearwise(files=[],names=[],valtypes=["sf","syst"]):
    binlist=['EBInc','EBHighR9','EBLowR9','EEInc','EEHighR9','EELowR9']
    output = schema.Category.parse_obj({
                "nodetype": "category",
                "input": "ValType",
                "content":[
                    schema.CategoryItem.parse_obj({
                        "key": val, 
                        "value": schema.Category.parse_obj({
                            "nodetype": "category",
                            "input": "WorkingPoint",
                            "content":[
                                schema.CategoryItem.parse_obj({
                                        "key": name, 
                                        "value": schema.Category.parse_obj({
                                            "nodetype": "category",
                                            "input": "HasPixBin",
                                            "content": [schema.CategoryItem.parse_obj({"key": binlist[i-1], 
                                                                                       "value": HasPixSFs(files,name,i,val)})
                                                        for i in range(1,7)
                                                       ],
                                        })
                                })
                                for name in names
                                ],
                            })
                    })
                    for val in valtypes
                ],
    })
    return output
