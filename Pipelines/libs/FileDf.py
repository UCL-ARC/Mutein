"""
RSA 4/5/22
------------------------
Class to manage reading and writing files
Specifically upper and lower casing

"""

import pandas as pd


class FileDf:
    def __init__(self, file_path, sep=",", header=True, cols=[]):
        self.file_path = file_path
        self.sep = sep
        self.header = header
        self.cols = cols

    def openDataFrame(self):
        print("FileDf opening", self.file_path)
        if self.header:
            df = pd.read_csv(self.file_path, sep=self.sep)
        else:
            df = pd.read_csv(self.file_path, sep=self.sep, header=None)

        if self.cols != []:
            df.columns = self.cols
        return df

    def getLines(self):
        with open(self.file_path) as fr:
            lines = fr.readlines()
            return lines
        return []


##################################################################################
class FileDic:
    def __init__(self, file_path, dic):
        self.file_path = file_path
        self.dic = dic
    
    def add(self,key,val,fill_last=[]):
        if key not in self.dic:
            self.dic[key] = []
        self.dic[key].append(val)    
    
    def fillLast(self,fill_last):
        for key in fill_last:
            if key in self.dic:            
                self.dic[key].append(self.dic[key][-1])    

    def saveAsDf(self, sep=",", header=True):
        len0 = 0
        for key,val in self.dic.items():
            if len0 == 0:
                len0 = len(val)
            else:
                if len(val) != len0:
                    print(len0,key,len(val))

        df = pd.DataFrame.from_dict(self.dic)
        df.to_csv(self.file_path, index=False, sep=sep, header=header)
        return df
