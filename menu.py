from tkinter import *
from eos_information import *
import numpy as np


#ecut plot executers
def ecut_diamond():
    exec(open("ecut_energies_diamond.py").read())


def ecut_betasn():
    exec(open("ecut_energies_betasn.py").read())


def ecut_difference():
    exec(open("ecut_energies_difference.py").read())

#kpoint plot executers
def kpoint_diamond():
    exec(open("kpoint_energies_diamond.py").read())


def kpoint_betasn():
    exec(open("kpoint_energies_betasn.py").read())


def kpoint_difference():
    exec(open("kpoint_energies_difference.py").read())


#eos plot executers
def eos_diamond():
    exec(open("eos_energies_diamond.py").read())


def eos_betasn():
    exec(open("eos_energies_betasn.py").read())


def eos_twophase():
    exec(open("eos_energies_joint.py").read())

#band plot executers
def bands_diamond():
    exec(open("bands_diamond.py").read())


def bands_betasn():
    exec(open("bands_betasn.py").read())


def wiki_format():
    exec(open("wiki_format_information.py").read())
    notice = Tk()
    notice_text = Label(notice, text='Wiki format information is located in the command line, copy information to pate on wiki from there')
    notice_text.pack()

def info_win():
    #information window
    info = Tk()

    diamond_label = Label(info, text='Diamond')
    diamond_label.grid(row=0, column=1)
    vol_0_diamond_lab = Label(info, text=vol_0_diamond)
    vol_0_diamond_lab.grid(row=1, column=1)
    k_0_diamond_lab = Label(info, text=k_0_diamond)
    k_0_diamond_lab.grid(row=2, column=1)
    k_0_prime_diamond_lab = Label(info, text=k_0_prime_diamond)
    k_0_prime_diamond_lab.grid(row=3, column=1)
    vol_t_diamond_lab = Label(info, text=vol_t_diamond)
    vol_t_diamond_lab.grid(row=4, column=1)
    betasn_label = Label(info, text='BetaSn')
    betasn_label.grid(row=0, column=2)
    vol_0_betasn_lab = Label(info, text=vol_0_betasn)
    vol_0_betasn_lab.grid(row=1, column=2)
    k_0_betasn_lab = Label(info, text=k_0_betasn)
    k_0_betasn_lab.grid(row=2, column=2)
    k_0_prime_betasn_lab = Label(info, text=k_0_prime_betasn)
    k_0_prime_betasn_lab.grid(row=3, column=2)
    vol_t_betasn_lab = Label(info, text=vol_t_betasn)
    vol_t_betasn_lab.grid(row=4, column=2)
    celldm_ratio_betasn_lab = Label(info, text=celldm_ratio_betasn)
    celldm_ratio_betasn_lab.grid(row=5, column=2)
    t_pres_lab = Label(info, text=t_pres)
    t_pres_lab.grid(row=6, column=1, columnspan=2)

    infolab1 = Label(info, text='vol_0')
    infolab1.grid(row=1, column=0)
    infolab2 = Label(info, text='K_{0}')
    infolab2.grid(row=2, column=0)
    infolab3 = Label(info, text='K_0_prime')
    infolab3.grid(row=3, column=0)
    infolab4 = Label(info, text='vol_t')
    infolab4.grid(row=4, column=0)
    infolab5 = Label(info, text='c/a')
    infolab5.grid(row=5, column=0)
    infolab6 = Label(info, text='Transition Pressure')
    infolab6.grid(row=6, column=0)
    copy_button = Button(info, text='wiki format', command=wiki_format)
    copy_button.grid(row=7,column=0,columnspan=3)
root = Tk()
#main menu
label1 = Label(root, text='QEInformation')
label1.grid()
#ecut menu
label2 = Label(root, text='Ecut Plots', height=3, width=20)
label2.grid(row=1, column=0)

button1  = Button(root, text='diamond', command=ecut_diamond)
button1.grid(row=2, column=0)

button2  = Button(root, text='betasn', command=ecut_betasn)
button2.grid(row=3, column=0)

button3  = Button(root, text='difference', command=ecut_difference)
button3.grid(row=4, column=0)
#kpoint menu
label3 = Label(root, text='Kpoint Plots',height=3, width=20)
label3.grid(row=1, column=1)

button4  = Button(root, text='diamond', command=kpoint_diamond)
button4.grid(row=2, column=1)

button5  = Button(root, text='betasn', command=kpoint_betasn)
button5.grid(row=3, column=1)

button6  = Button(root, text='difference', command=kpoint_difference)
button6.grid(row=4, column=1)
#eos menu
label4 = Label(root, text='EOS Plots',height=3, width=20)
label4.grid(row=1, column=2)

button7  = Button(root, text='diamond', command=eos_diamond)
button7.grid(row=2, column=2)

button8  = Button(root, text='betasn', command=eos_betasn)
button8.grid(row=3, column=2)

button9  = Button(root, text='two-phase', command=eos_twophase)
button9.grid(row=4, column=2)
#bands menu
label5 = Label(root, text='Bands Plots',height=3, width=20)
label5.grid(row=1, column=3)

button10 = Button(root, text='diamond', command=bands_diamond)
button10.grid(row=2, column=3)

button11 = Button(root, text='betasn', command=bands_betasn)
button11.grid(row=3, column=3)

#information
label16 = Label(root, text='Information',height=3, width=20)
label16.grid(row=1, column=4)
button12 = Button(root, text='Show Info', command=info_win)
button12.grid(row=2, column=4)

root.mainloop()
