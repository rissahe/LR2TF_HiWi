import numpy as np
import pandas as pd
from importlib.resources import files



#Will not work until package directory is created, please enter regulon csv manually.
def load_dorothea_regulon (organism):
    if (organism == "human"):
        dorothea_regulon_human = files("LR2TF_py.data").joinpath("human_dorothea_reg.csv").pd.read_csv(index_col=0)
        regulon =  dorothea_regulon_human.loc[dorothea_regulon_human["confidence"].isin(["A","B", "C", "D"])]
        regulon = pd.DataFrame.rename(regulon, columns={"source" : "tf"})

    elif (organism == "mouse"):
        dorothea_regulon_mouse = files("LR2TF_py.data").joinpath("mouse_dorothea_reg.csv").pd.read_csv(index_col=0)
        regulon =  dorothea_regulon_mouse.loc[dorothea_regulon_mouse["confidence"].isin(["A","B", "C", "D"])]
        regulon = pd.DataFrame.rename(regulon, columns={"source" : "tf"})
        
    else:
        print("Only human and mouse regulons can be loaded by default!")
    return regulon


def validate_input_arguments (arguments_list):
    if arguments_list["out_path"] is None:
        print("Please provide an output path")
    elif arguments_list["out_path"][-1] != "/":
        arguments_list["out_path"] = arguments_list["out_path"] + "/"

    if arguments_list["celltype"] is None:
        print("Please provide the name of the metadata field containing cell type annotations")

    if arguments_list["condition"] is None:
        print("Please provide the name of the metadata field containing condition annotations")

    if arguments_list["organism"] is None:
        arguments_list["organism"] = "human"

    if arguments_list["comparison_list"] is None:
        arguments_list["comparison_list"] = np.nan

    if arguments_list["logfc"] is None:
        arguments_list["logfc"] = 0.0

    if arguments_list ["pval"] is None:
        arguments_list["pval"] = 0.05

    if arguments_list["reg"] is None:
        arguments_list["reg"] = load_dorothea_regulon(arguments_list["organism"])

    elif isinstance(arguments_list["reg"], str):
        arguments_list["reg"] = pd.read_csv(arguments_list["reg"], index_col=0)
        arguments_list["reg"] = pd.DataFrame.rename(arguments_list["reg"], columns={"source" : "tf"})

    if not "tf" in arguments_list["reg"] and "target" in arguments_list["reg"] and "weight"in arguments_list["reg"]:
        raise Exception("Not all necessary columns found in regulon table! Please make sure that the regulon has the columns source, target and weight!")
    
    return(arguments_list)


#Please edit the following argument list to suit your own data.

#Enter your regulon csv path for reg

arguments = {"out_path" : None, "celltype" : None, "condition" : None, "organism" : None, "comparison_list" : None, "logfc" : None, "pval" : None, "reg" : "/home/larissa/Documents/LR2TF_HiWi/human_dorothea_reg.csv"}

arguments_list = validate_input_arguments(arguments)

print(arguments_list["reg"])