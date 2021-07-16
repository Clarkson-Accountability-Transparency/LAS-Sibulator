import tkinter as tk

from sibulator import *

class SibulatorGUI():
    def __init__(self, master):
        self.master = master
        self.frame = tk.Frame(master)
        self.create_widgets()

    def create_widgets(self):
        # create labels and entries
        lbl_known_sample_path = tk.Label(text='Path to known sample .csv file')
        ent_known_sample_path = tk.Entry(width=75)
        self.ent_known_sample_path = ent_known_sample_path

        lbl_allele_freq_data = tk.Label(text='Allele frequency database to use (nist or strider)')
        ent_allele_freq_data = tk.Entry(width=75)
        self.ent_allele_freq_data = ent_allele_freq_data

        lbl_subpop = tk.Label(text='Subpopulation of allele frequency database to use (see documentation based on chosen database)')
        ent_subpop = tk.Entry(width=75)
        self.ent_subpop = ent_subpop

        lbl_num_samples = tk.Label(text='Number of sibling samples to generate')
        ent_num_samples = tk.Entry(width=75)
        self.ent_num_samples = ent_num_samples

        lbl_dest_path = tk.Label(text='Name of file you want to save generated samples to')
        ent_dest_path = tk.Entry(width=75)
        self.ent_dest_path = ent_dest_path

        # display all widgets
        lbl_known_sample_path.pack()
        ent_known_sample_path.pack()

        lbl_allele_freq_data.pack()
        ent_allele_freq_data.pack()

        lbl_subpop.pack()
        ent_subpop.pack()

        lbl_num_samples.pack()
        ent_num_samples.pack()

        lbl_dest_path.pack()
        ent_dest_path.pack()

        btn_run = tk.Button(text='Run', command=self.run_simulation)
        btn_run.pack()

    def run_simulation(self):
        # read in data from user
        known_sample_path = self.ent_known_sample_path.get()
        allele_freq_data = self.ent_allele_freq_data.get()
        subpop = self.ent_subpop.get()
        num_samples = int(self.ent_num_samples.get())
        dest_path = self.ent_dest_path.get()

        # run simulation
        # print(known_sample_path, allele_freq_data, subpop, num_samples, dest_path)
        run_simulation(known_sample_path, allele_freq_data, subpop, num_samples, dest_path)

        # display done and wait for user to close
        self.frame.destroy()
        self.frame = tk.Frame(self.master)
        lbl_done = tk.Label(text='Done!', height=10)
        btn_quit = tk.Button(text='Quit', command=self.master.destroy)

        lbl_done.pack()
        btn_quit.pack()
        self.frame.pack()

window = tk.Tk()
app = SibulatorGUI(window)
window.mainloop()

