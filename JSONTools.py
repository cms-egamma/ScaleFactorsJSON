# Technical part below
from ROOT import TFile, TH1F, TCanvas, TString
import correctionlib.schemav2 as schema
from correctionlib.schemav2 import Correction
from math import fabs
import pandas as pd
import numpy as np

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

def getScaleUnc(year,val,etabin,gain):
    
    SandSdict={
        "2016preVFP" :{ "EBin":{12 : 0.05, 6 : 0.5, 1 : 2.0}, "EBout":{12 : 0.075, 6 : 0.5, 1 : 2.0}, "EE":{12 : 0.15,  6 : 1.0,  1 : 3.0} },
        "2016postVFP":{ "EBin":{12 : 0.05, 6 : 0.5, 1 : 2.0}, "EBout":{12 : 0.075, 6 : 0.5, 1 : 2.0}, "EE":{12 : 0.15,  6 : 1.0,  1 : 3.0} },
        "2017"        :{ "EBin":{12 : 0.05, 6 : 0.5, 1 : 1.0}, "EBout":{12 : 0.05,  6 : 0.5, 1 : 1.0}, "EE":{12 : 0.1,   6 : 0.5,  1 : 2.0} },
        "2018"        :{ "EBin":{12 : 0.05, 6 : 0.5, 1 : 1.0}, "EBout":{12 : 0.05,  6 : 0.5, 1 : 1.0}, "EE":{12 : 0.125, 6 : 0.75, 1 : 2.0} }
    }
    
    if val=="scaleup":
        return (1+ SandSdict[year][etabin][gain]/100.)
    elif val=="scaledown": 
        return (1- SandSdict[year][etabin][gain]/100.)
    else:
        return (SandSdict[year][etabin][gain]/100.)

def getScaleUncbins(year,val,gain):
    etabins=[float('-inf'),-1.442, -1, 0, 1, 1.442, float('inf')]
    return schema.MultiBinning.parse_obj({
        "inputs":["eta"],
        "nodetype": "multibinning",
        "edges": [etabins],
        "content": [getScaleUnc(year,val,etabin,gain) for etabin in ['EE','EBout','EBin','EBin','EBout','EE']] ,
        "flow": 'error'})

def ScaleUnc(year,valtypes=["scaleup","scaledown","scaleunc"]):
    gains=[12,6,1]
    output = schema.Category.parse_obj({
                "nodetype": "category",
                "input": "ValType",
                "content":[
                    schema.CategoryItem.parse_obj({
                        "key": val, 
                        "value":schema.Category.parse_obj({
                            "nodetype": "category",
                            "input": "gain",
                            "content":[schema.CategoryItem.parse_obj({"key": gain, 
                                                                      "value": getScaleUncbins(year,val,gain)})
                                       for gain in gains
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

def createDictFromPandas(df):
        if (df.index.nlevels==1):
            return df.to_dict()
        dict_f = {}
        for level in df.index.levels[0]:
            if (level in df.index):
                dict_f[level] = createDictFromPandas(df.xs([level]))
        return dict_f

def getSSdictRun2(year):
    
    url_dict = {
        "2016preVFP": "https://raw.githubusercontent.com/cms-data/EgammaAnalysis-ElectronTools/master/ScalesSmearings/Run2016_UltraLegacy_preVFP_RunFineEtaR9Gain_scales.dat",
        "2016postVFP": "https://raw.githubusercontent.com/cms-data/EgammaAnalysis-ElectronTools/master/ScalesSmearings/Run2016_UltraLegacy_postVFP_RunFineEtaR9Gain_scales.dat",
        "2017": "https://raw.githubusercontent.com/cms-data/EgammaAnalysis-ElectronTools/master/ScalesSmearings/Run2017_24Feb2020_runEtaR9Gain_v2_scales.dat",
        "2018": "https://raw.githubusercontent.com/cms-data/EgammaAnalysis-ElectronTools/master/ScalesSmearings/Run2018_29Sep2020_RunFineEtaR9Gain_scales.dat",
        "Prompt2022FG":"https://akapoordocs.web.cern.ch/akapoordocs/step4closure_Prompt2022FG_28_06_2023_v0_scales.dat"
    }
    
    print(url_dict[year])
    data = np.genfromtxt(url_dict[year])
    df = pd.DataFrame(np.squeeze(data), columns=['run_bin_low', 'run_bin_high', 'eta_bin_low', 'eta_bin_high',
                                                 'r9_bin_low',"r9_bin_high","et_bin_low","et_bin_high",
                                                 "gain_seed","total_correction","total_uncertainty"])
    #df = df.iloc[1:]
    #df = df.iloc[:-1]
    
    #temporary hack to fix eta boundaries
    df.loc[(df['eta_bin_low'] >= 1.56) & (df['eta_bin_low'] <= 1.57), 'eta_bin_low'] = 1.566
    df.loc[(df['eta_bin_high'] >= 1.56) & (df['eta_bin_high'] <= 1.57), 'eta_bin_high'] = 1.566
    
    df.loc[(df['eta_bin_low'] >= 1.44) & (df['eta_bin_low'] <= 1.45), 'eta_bin_low'] = 1.442
    df.loc[(df['eta_bin_high'] >= 1.44) & (df['eta_bin_high'] <= 1.45), 'eta_bin_high'] = 1.442
    
    
    mydf=df.set_index(['run_bin_low','run_bin_high','eta_bin_low',
                       'eta_bin_high','r9_bin_low',"r9_bin_high","et_bin_low","et_bin_high","gain_seed"])
    newdf=createDictFromPandas(mydf)
    return newdf

def get_SS_runbins_and_SS_dict(year):
    newdf=getSSdictRun2(year)
    runbins=[]
    runbinsa=list(newdf.keys())
    for runbina in runbinsa:
        runbins.append([runbina,list(newdf[runbina].keys())[0]])
    return runbins,newdf

def returnsmearingdf(url):
    import pandas as pd
    from io import StringIO
    import requests

    # Fetch the content from the URL
    response = requests.get(url)
    content = response.text
    
    # Skip the header line and create a StringIO object
    data_io = StringIO(content)
    next(data_io)  # Skip the header line

    # Define column names for the DataFrame
    column_names = ['category', 'Emean', 'err_Emean', 'rho', 'err_rho', 'phi', 'err_phi']

    # Read the data into a DataFrame
    df = pd.read_csv(data_io, sep='\t', header=None, names=column_names)

    # Split the 'category' column into four columns
    df[['abseta_min', 'abseta_max', 'R9_min', 'R9_max']] = df['category'].str.split('-|_', expand=True)[[1, 2, 4, 5]]

    # Convert the new columns to numeric
    df[['abseta_min', 'abseta_max', 'R9_min', 'R9_max']] = df[['abseta_min', 'abseta_max', 'R9_min', 'R9_max']].astype(float)

    df.loc[(df['abseta_min'] >= 1.56) & (df['abseta_min'] <= 1.57), 'abseta_min'] = 1.566
    df.loc[(df['abseta_max'] >= 1.56) & (df['abseta_max'] <= 1.57), 'abseta_max'] = 1.566
    
    df.loc[(df['abseta_min'] >= 1.44) & (df['abseta_min'] <= 1.45), 'abseta_min'] = 1.442
    df.loc[(df['abseta_max'] >= 1.44) & (df['abseta_max'] <= 1.45), 'abseta_max'] = 1.442
    # Drop the original 'category' column
    df.drop('category', axis=1, inplace=True)

    # Print the resulting DataFrame
    return df

def getSmearingdictRun2(year):
    
    url_dict = {
        "Prompt2022FG":"https://akapoordocs.web.cern.ch/akapoordocs/step2_Prompt2022FG_26_06_2023_v0_smearings.dat"
    }
    
    return returnsmearingdf(url_dict[year])

