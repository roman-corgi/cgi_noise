# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 17:22:28 2023

@author: david
"""
import pandas as pd
import xlwings as xw
import os

class MyExcelApp:
    def __init__(self):
        self.excel = xw.App(visible=False)

    def __enter__(self):
        return self.excel

    def __exit__(self, exc, value, traceback):
        # Handle your app-specific exceptions (exc) here
        self.excel.quit()
        return True   
        # ^ return True only if you intend to catch all errors in here.
        # Otherwise, leave as is and use try... except on the outside.

class MyExcelWorkbook:
    def __init__(self, xlapp, bookname):
        self.workbook = xlapp.books.open(bookname)

    def __enter__(self):
        return self.workbook

    def __exit__(self, exc, value, traceback):
        # Handle your workbook specific exceptions (exc) here
        # self.workbook.save()   # depends what you want to do here
        self.workbook.close()
        return True   
        # ^ return True only if you intend to catch all errors in here.
        # Otherwise, leave as is and use try... except on the outside.   

class loadXLcol:
    
    """Class loadXLcol opens Excel and converts to pandas dataframe
    
    Parameters ---
    filename: str (e.g. C: User\Desktop\CGI_Perf\EBcsvData\Photometry\DET_CBE_200318.xlsx)
        
    Attributes ---
    fullfile:str (same as filename parameter above)
    xlname: str (aka os.path.basename(filename))(e.g. DET_CBE_200318.csv)
    prefix:str (type of data (DET, CGPERF, THPT, etc.) ) 
   """     
        
    
    
    def __init__(self, filename, nrows):        
        self.fullfile = filename
        self.xlname = os.path.basename(filename)
        self.prefix = ((os.path.basename(filename)).split("_"))[0]      
       
        with MyExcelApp() as app:
            with MyExcelWorkbook(app, filename) as wb:
                
                sheet1 = wb.sheets[0]
                # sheet1 = wbs.sheets[0]
                range = sheet1.range('A1')
                self.colNo = range.current_region.last_cell.column            
                self.ScenData = sheet1.range((1,self.colNo),(nrows,self.colNo)).value
                self.RowNames = sheet1.range((1,1),(nrows,1)).value
       
        
    @property
    def df(self): #getDataFrame(self):
        dataf = pd.DataFrame(data = self.ScenData, index=self.RowNames, columns=['Latest'])
        return dataf 