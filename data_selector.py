# -*- coding: cp1252 -*-
"""
Data selector : a pop-up window which lets the user choose whoch data use for multiple scripts of the nanowire vibrations project

Author :
    Clément Chardin
    clement.chardin@nottingham.ac.uk
"""
import os
from tkinter import *
from tkinter import filedialog

class DataSelector(object):
    def __init__(self, description):
        super(DataSelector, self).__init__()
        self.directory = 'None'
        self.description = description
        self.setup_ui()

    def GetDirectory(self):
        self.directory = filedialog.askdirectory()
        self.text.set(self.directory)

    def end(self):
        self.window.destroy()

    def setup_ui(self):
        self.window = Tk()
        self.window.title("Directory Selector")

        self.description_label = Label(self.window, text=self.description)

        self.label = Label(self.window, text="Current directory :")

        self.text = StringVar()
        self.text.set(self.directory)
        self.directory_label = Label(self.window, textvariable=self.text)

        self.getdir_button = Button(self.window, text="Choose directory", command=self.GetDirectory)

        self.sens_label = Label(self.window, text="Sens (10 ou 01) :")
        self.sens_var = StringVar(value='10')
        self.sens_entry = Entry(self.window, textvariable=self.sens_var)

        self.AR_label = Label(self.window, text="AR (A, R, both or None) :")
        self.AR_var = StringVar(value='None')
        self.AR_entry = Entry(self.window, textvariable=self.AR_var)

        self.savefigs_var = BooleanVar(value=True)
        self.savefigs_check = Checkbutton(self.window, text="Save figures", variable=self.savefigs_var)

        self.savefiles_var = BooleanVar(value=True)
        self.savefiles_check = Checkbutton(self.window, text="Save files", variable=self.savefiles_var)

        self.end_button = Button(self.window, text="Use this directory", command=self.end)

        self.description_label.pack()
        self.label.pack()
        self.directory_label.pack()
        self.getdir_button.pack()
        self.sens_label.pack()
        self.sens_entry.pack()
        self.AR_label.pack()
        self.AR_entry.pack()
        self.savefigs_check.pack()
        self.savefiles_check.pack()
        self.end_button.pack()
        self.window.geometry("400x400")

    def select_directory(self):
        self.window.mainloop()
