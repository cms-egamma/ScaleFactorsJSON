# Technical part below
from ROOT import TFile, TH1F, TCanvas, TString
import correctionlib.schemav2 as schema
from correctionlib.schemav2 import Correction
from math import fabs
corrs=[]
#Function to extract SFs from EGamma standard root files
def getSFs(fn="filename",IsSF="sf",sfhist="EGamma_SF2D"):
    tf = TFile(fn)
    fo = tf.Get(sfhist)
    Xbins=fo.GetNbinsX()
    Ybins=fo.GetNbinsY()
    X=[None]*(Xbins+1)
    
    Y=[None]*(Ybins+1)
    #print("Y: "+str(Y))
    values=[]
    errors=[]
    for i in range(1,Xbins + 1):
        X[i-1]=(fo.GetXaxis().GetBinLowEdge(i)) if fo.GetXaxis().GetBinLowEdge(i)!=-2.5 else float('-inf')
    X[Xbins]=fo.GetXaxis().GetBinUpEdge(Xbins) if fo.GetXaxis().GetBinUpEdge(Xbins)!=2.5 else float('inf')
    print("X: "+str(X))
    for j in range(1,Ybins + 1):
        Y[j-1]=(fo.GetYaxis().GetBinLowEdge(j))
    Y[Ybins]=fo.GetYaxis().GetBinUpEdge(Ybins) if fo.GetYaxis().GetBinUpEdge(Ybins)!=500 else float('inf')
    print("Y: "+str(Y))
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
    if IsSF=="sfup":
        valerrorsup=schema.MultiBinning.parse_obj({
            "inputs":["eta","pt"],
            "nodetype": "multibinning",
            "edges": [
                X,
                Y,
            ],
            "content": [m + n for m,n in zip(values,errors)],
            "flow": 'error',
        })
        return valerrorsup
    if IsSF=="sfdown":
        valerrorsdown=schema.MultiBinning.parse_obj({
            "inputs":["eta","pt"],
            "nodetype": "multibinning",
            "edges": [
                X,
                Y,
            ],
            "content": [m - n for m,n in zip(values,errors)],
            "flow": 'error',
        })
        return valerrorsdown
    
def SFyearwise(files=[],names=[],valtypes=["sf","sfup","sfdown"]):
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
        if IsSF=="sfup":
            SF=hist.GetBinContent(i)+hist.GetBinError(i)
        if IsSF=="sfdown":
            SF=hist.GetBinContent(i)-hist.GetBinError(i)
    return SF

def CSEVSFyearwise(files=[],names=[],valtypes=["sf","sfup","sfdown"]):
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
        if IsSF=="sfup":
            SF=hist.GetBinContent(i)+hist.GetBinError(i)
        if IsSF=="sfdown":
            SF=hist.GetBinContent(i)-hist.GetBinError(i)
    return SF

def HasPixSFyearwise(files=[],names=[],valtypes=["sf","sfup","sfdown"]):
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
