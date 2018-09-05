#! /usr/bin/env python
# -*- coding: utf-8 -*-

from openpyxl import Workbook
from openpyxl.styles import Font, Fill, Border, Side
from string import ascii_uppercase 

wb = Workbook()
wb.guess_types = True

def format_table(work_sheet, input_file, title):
    
    work_sheet["A1"].value = title
    work_sheet["A1"].font = Font(name = "Helvetica", size = 12, bold = True)

    col_len = {} 
    italic_ind = [] 
    with open(input_file, 'r') as hin:
        header = hin.readline().rstrip('\n').split('\t')
        for i in range(len(header)):
            work_sheet[ascii_uppercase[i] + '3'].value = header[i]
            work_sheet[ascii_uppercase[i] + '3'].font = Font(name = "Helvetica", size = 11, bold = True)        
            work_sheet[ascii_uppercase[i] + '3'].border = Border(top = Side(border_style = "medium", color = "FF000000"),
                                                                 bottom = Side(border_style = "medium", color = "FF000000"))

            if header[i] in ["Gene", "Associated splicing factor mutation"]: italic_ind.append(i)

            col = work_sheet[ascii_uppercase[i] + '3'].column
            if col not in col_len: col_len[col] = 0
            col_len[col] = max(col_len[col], len(header[i]))


        row = 0
        for line in hin:
            F = line.rstrip('\n').split('\t')
            for i in range(len(F)):
                work_sheet[ascii_uppercase[i] + str(row + 4)].value = F[i]
                if i in italic_ind:
                    work_sheet[ascii_uppercase[i] + str(row + 4)].font = Font(name = "Helvetica", size = 11, italic = True)
                else:
                    work_sheet[ascii_uppercase[i] + str(row + 4)].font = Font(name = "Helvetica", size = 11)

                col = work_sheet[ascii_uppercase[i] + str(row + 4)].column
                if col not in col_len: col_len[col] = 0
                col_len[col] = max(col_len[col], len(F[i]))

            row = row + 1
        
    def as_text(value): return str(value) if value is not None else ""

    for i in col_len:
        work_sheet.column_dimensions[i].width = float(col_len[i]) * 1.2


ws1 = wb.active
ws1.title = "S14_SAV_list"
format_table(ws1, "../output/table/stable.raw.txt", "Supplementary Table 14.  List of splicing associated variant identified by SAVNet")

wb.save("../output/table/Supplementary_Table.xlsx")



