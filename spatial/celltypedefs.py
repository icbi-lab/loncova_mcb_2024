import pandas as pd


class ctdef:
    def __init__(self, granularity, def_table):
        self.granularity = granularity
        self.def_table = def_table

        use_lumped = True if (self.granularity == "lumped") else False

        # Cell type definitions
        celltype_def = pd.read_csv(self.def_table, sep="\t")

        self.cellTypeMap = {}
        self.cellTypeColorMap = {}
        self.cellTypeNameMap = {}
        self.cellTypeClassMap = {}
        self.cellTypeClassColors = {}
        self.cellTypeLumpedMap = {}
        self.cellTypeExclude = {}

        self.immuneCellIdx = []
        self.tumorCellIdx = []

        for idx, row in celltype_def.iterrows():
            if use_lumped:
                self.cellTypeLumpedMap[row["celltype"]] = row["celltype_lumped"]
                self.cellTypeMap[row["celltype_lumped"]] = row["number_lumped"]
                self.cellTypeColorMap[row["celltype_lumped"]] = row["color_lumped"]
                self.cellTypeNameMap[row["celltype_lumped"]] = row["cname_lumped"]
                self.cellTypeClassMap[row["celltype_lumped"]] = row["celltype_class"]
                if row["celltype_class"] in ["Myelocytes", "Macorphages", "Lymphocytes", "Immunocytes"]:
                    self.immuneCellIdx.append(row["number_lumped"])
                if row["celltype_class"] in ["Tumorcells"]:
                    self.tumorCellIdx.append(row["number_lumped"])

                self.cellTypeExclude[row["celltype_lumped"]] = row["exclude"]

            else:
                self.cellTypeMap[row["celltype"]] = row["number"]
                self.cellTypeColorMap[row["celltype"]] = row["color"]
                self.cellTypeNameMap[row["celltype"]] = row["cname"]
                self.cellTypeClassMap[row["celltype"]] = row["celltype_class"]
                if row["celltype_class"] in ["Myelocytes", "Macorphages", "Lymphocytes", "Immunocytes"]:
                    self.immuneCellIdx.append(row["number"])
                if row["celltype_class"] in ["Tumorcells"]:
                    self.tumorCellIdx.append(row["number"])

                self.cellTypeExclude[row["celltype"]] = row["exclude"]

            self.cellTypeClassColors[row["celltype_class"]] = row["color_celltype_class"]

        # make unique
        self.immuneCellIdx = list(set(self.immuneCellIdx))
        self.tumorCellIdx = list(set(self.tumorCellIdx))

        # map cell types to the cell type call color
        self.cellTypeClassColorMap = {
            n: self.cellTypeClassColors[self.cellTypeClassMap[n]] for n in self.cellTypeClassMap
        }

        # In[4]:

        # map cell type colors to cell type names used in manuscript
        self.cellTypeNameColorMap = {self.cellTypeNameMap[n]: self.cellTypeColorMap[n] for n in self.cellTypeMap}

        # In[5]:

        # map cell type names as used in manuscript to cell type numbers (01, 02, 03 ...) as used in manuscript
        self.cellTypeNameColorMapIdx = {
            self.cellTypeNameMap[n]: str(self.cellTypeMap[n] + 1).zfill(2) for n in self.cellTypeMap
        }

        # In[6]:

        # map cell type numbers (01, 02, 03 ...) as used in manuscript to cell type colors used in manuscript
        self.cellTypeIdxColorMap = {
            self.cellTypeNameColorMapIdx[x]: self.cellTypeNameColorMap[x] for x in self.cellTypeNameColorMap
        }

        # In[7]:

        # make pandas series
        self.cellTypeIdxColorMap_s = pd.Series(self.cellTypeIdxColorMap)

        # In []:

        # map cell type numbers (01, 02, 03 ...) as used in manuscript to cell type colors used in manuscript
        self.cellTypeClassIdxColorMap = {
            str(list(self.cellTypeColorMap.keys()).index(n) + 1).zfill(2): self.cellTypeClassColors[
                self.cellTypeClassMap[n]
            ]
            for n in self.cellTypeClassMap
        }
        self.cellTypeDummyIdxColorMap = {
            str(list(self.cellTypeColorMap.keys()).index(n) + 1).zfill(2): "white" for n in self.cellTypeClassMap
        }

        self.cellTypeClassIdxColorMap_s = pd.Series(self.cellTypeClassIdxColorMap)
        self.cellTypeDummyIdxColorMap_s = pd.Series(self.cellTypeDummyIdxColorMap)

        # In[8]:

        # make pandas series
        self.cellTypeColorMap_s = pd.Series(self.cellTypeNameColorMap)

        # In[9]:

        # swap key/value of custom_cellTypeMapOrder (custom order number is key and name is value now)
        self.inv_cellTypeMapOrder = {self.cellTypeMap[x]: x for x in self.cellTypeMap}

        # get custom sorted cell type numbers as used in manuscript
        self.cellTypeMapOrder_idx = [
            str(i + 1).zfill(2) for i in list(dict(sorted(self.inv_cellTypeMapOrder.items())).keys())
        ]

        # In[10]:
        # swap key/value of cellTypeNameMap (cellTypeNameMap name is key and data celltype name is value now)
        self.inv_cellTypeNameMap = {self.cellTypeNameMap[x]: x for x in self.cellTypeNameMap}
