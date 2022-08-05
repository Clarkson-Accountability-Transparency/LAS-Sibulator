'''
Implements sibling simulator for creation of profiles related to a given known profile
Next steps:
    - Support reading in different file types for known sample
    - Calibrate replacement probabilities of 0, 1, 2 to population prevalence?
'''
import pandas as pd
import numpy as np

import os
import argparse
from bs4 import BeautifulSoup





class Sibulator:
    def __init__(self, *args, **kwargs):
        self.load_known_sample(kwargs.get('knownSampleFile', None))
        self.subpop = kwargs.get('subpop', None)
        self.alleleFreqDict = {} 

    def load_known_sample(self, knownSampleFile):
        '''
        Reads in known sample from file and puts information into a python dict with following format
        {[site_name]: [allele1, allele2], ... }
        '''
        if knownSampleFile is None:
            return
        self.knownSampleFile = knownSampleFile

        file_format = self.knownSampleFile.split('.')[-1]

        if file_format == 'txt':
            # assume text file is formatted as necessary python dict
            with open(self.knownSampleFile, 'r') as f:
                known_sample = eval(f.read()) #TODO this isn't safe

        elif file_format == 'csv':
            data = pd.read_csv(self.knownSampleFile)
   
            # account for vertical format
            if data.shape[1] < 3:
                data.columns = data.iloc[0].values # assume there is some identifier and actual column names in first row
                data = data.drop(0).reset_index(drop=True).transpose() # convert to horizontal orientation
                data.columns = data.iloc[0].values # first row is the allele location
                data.drop(data.index[0], inplace=True) # get rid of first location
    
            data = data.iloc[0] # assume first row is desired known profile
    
            self.knownSample = {}
            for col in ['CaseNumber', 'Sample']:
                if col in data.index:
                    data.drop(col, inplace=True)
            for loc in data.index:
                allele = data[loc]
                key = loc.split('.')[0].replace('_', ' ') # standardize key values (contributor profiles use space, database uses _)
                self.knownSample[key] = self.knownSample.get(key, []) + [allele]
        else:
            raise NotImplementedError('File type: .{} not yet supported'.format(file_format))
        
        return self.knownSample

    def load_allele_frequencies(self, *vargs):
        pass

    def generate_sibling(self, subpop):
        ''' Returns a sibling'''
        rand_sample = {'numReplaced':0}
        for loci, allele in self.knownSample.items():
            # at each loci, choose to either replae 0, 1, or 2 alleles
            #print(self.alleleFreqDict[loci].keys())
            a1, a2 = allele
            test = np.random.random()
            if test < 0.25:
                # replace none
                new_alleles = [a1, a2]
            elif test < 0.75:
                # replace one
                p = self.alleleFreqDict[loci][subpop]['frequency']
                p = p / np.sum(p) # softmax b/c frequencies don't always add to exactly 1
                a3 = np.random.choice(a=self.alleleFreqDict[loci][subpop]['alleles'], size=1, p=p)[0]
                if np.random.random() < 0.5:
                    # replace 1st
                    new_alleles = [a3, a2]
                else:
                    # replace 2nd
                    new_alleles = [a1, a3]
                rand_sample['numReplaced'] += 1
            else:
                # replace both
                p = self.alleleFreqDict[loci][subpop]['frequency']
                p = p / np.sum(p) # softmax b/c frequencies don't always add to exactly 1
                a3, a4 = np.random.choice(a=self.alleleFreqDict[loci][subpop]['alleles'], size=2, p=p)
                new_alleles = [a3, a4]
                rand_sample['numReplaced'] += 2
    
            rand_sample[loci] = new_alleles
        
        return rand_sample

    def generate_sibling_samples(self,  numSiblings):
        '''
        Generates possible sibling DNA profiles given known sample and user-given specifications
        '''
        #simulated_samples = []
        #for i in range(num_siblings):
        #    rand_sample = self.generate_sibling() generate_one_sample()
        #    simulated_samples.append(rand_sample)
    
        return list([self.generate_sibling(self.subpop) for i in range(numSiblings)])


class NistSibulator(Sibulator):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._NIST_DIR         = 'NIST'
        self._AFAM_FREQ_FILE   = 'NIST Fusion AfAm_Amended2017[1].csv'
        self._ASIAN_FREQ_FILE  = 'NIST Fusion Asian_Amended2017[1].csv'
        self._CAUC_FREQ_FILE   = 'NIST Fusion Cauc_Amended2017[1].csv'
        self._HISP_FREQ_FILE   = 'NIST Fusion Hisp_Amended2017[1].csv'
        self.load_allele_frequencies()

    def load_allele_frequencies(self, *vargs):
        num_vargs = len(vargs)
        if num_vargs == 5:
            self._NIST_DIR = vargs[0]
            self._AFAM_FREQ_FILE = vargs[1]
            self._ASIAN_FREQ_FILE = vargs[2]
            self._CAUC_FREQ_FILE = vargs[3]
            self._HISP_FREQ_FILE = vargs[4]
        if num_vargs != 5 and num_vargs != 0:
            raise ArgumentError('Expected either no arguments or 5 (Dir, afam, asian, cauc, hispanic)') 
             
        # read data in
        afam_data  = pd.read_csv(os.path.join(self._NIST_DIR, self._AFAM_FREQ_FILE))
        asian_data = pd.read_csv(os.path.join(self._NIST_DIR, self._ASIAN_FREQ_FILE)) 
        cauc_data  = pd.read_csv(os.path.join(self._NIST_DIR, self._CAUC_FREQ_FILE))
        hisp_data  = pd.read_csv(os.path.join(self._NIST_DIR, self._HISP_FREQ_FILE)) 
        data = {'AfAm': afam_data, 'Asian': asian_data, 'Cauc': cauc_data, 'Hisp': hisp_data}
        N = []
        for df in data.values():
            df_n = df[df['Allele'] == 'N']
            N.append(int(df_n.values[0][1]))

            # drop allele from index
            df.drop(df_n.index[0], axis=0, inplace=True)

        # convert data to dictionaries in expected format, also combine to create a total population frequency
        i = 0
        total_allele_counts = {}
        for subpop, df in data.items():
            for col in df.drop('Allele', axis=1).columns:
                key = col.replace('_', ' ') # standardize key values (contributor profiles use space, database uses _)
                # column is the location name ex. D3S1358
                self.alleleFreqDict[key] = {**self.alleleFreqDict.get(key, {}), **{subpop: {'alleles': df['Allele'].astype('float'), 'frequency': df[col].astype('float')}}}

                # update total counts of each allele
                n = N[i]
                for j, allele in enumerate(df['Allele']):
                    # iterate over Series
                    allele_freq = df.iloc[j][col] # allele frequency at current location
                    allele_count = allele_freq * n
                    temp_dict = total_allele_counts.get(key, {}) # get allele counts associated w/ location
                    temp_dict[allele] = temp_dict.get(allele, 0) + allele_count  # update counts
                    total_allele_counts[key] = temp_dict # replace counts
            i += 1

        # add total population allele frequencies to data
        for loc, counts in total_allele_counts.items():
            for k, v in counts.items():
                total_allele_counts[loc][k] = v / sum(N) # normalize counts by total number of profiles
            self.alleleFreqDict[loc]['total'] = {'alleles': [float(k) for k in list(counts.keys())], 
                                              'frequency': [float(f) for f in list(counts.values())]}


        return self.alleleFreqDict

class StriderSibulator(Sibulator):
    # Reading the data inside the xml file to a variable under the name data
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._FREQ_FILE = 'STRidER_frequencies_2019-08-02.xml'
        self.load_allele_frequencies()

    def load_allele_frequencies(self, *vargs):
        '''
        Loads allele frequency data as seen from STRidER (https://strider.online/frequencies)
        '''
        if len(vargs)> 0:
            self._FREQ_FILE = vargs[0]
        with open(self._FREQ_FILE, 'r') as f:
            data = f.read()

        # Passing the stored data inside the beautifulsoup parser, storing the returned object 
        bs_data = BeautifulSoup(data, "xml")

        # turn xml into dict for use
        marker_data = bs_data.find_all('marker')
        # print(len(marker_data)) # number of markers 
        for md in marker_data:
            loc_name = md.find('name').text # get loci name
            self.alleleFreqDict[loc_name] = {}
            
            alleles_list = md.find('alleles').text.split(', ') # list of alleles at this site
            origins = md.find_all('origin') # frequency data for different regions/subpopulations
            for orig in origins:
                orig_name = orig.get('name') # name of region
                frequency_dict = {}
                allele_freqs = orig.find_all('frequency') # get allele frequencies for this region
                # add frequency data to dict for all that appear in this specific region
                for allele in allele_freqs:
                    al = allele.get('allele')
                    freq = float(allele.text)
                    frequency_dict[al] = float(allele.text)

                # fill in missing alleles with frequency 0
                for al in alleles_list:
                    frequency_dict[al] = frequency_dict.get(al, 0)
                    
                al_list = list(frequency_dict.keys()) # get allele names in list
                freq_list = list(frequency_dict.values()) # get allele frequencies in list
                # print(sum(freq_list)) # for verifying sum to approx 1
           
                #TODO: This seems wrong. 
                self.alleleFreqDict[loc_name][orig_name] = {'alleles': al_list} # add info to returned dict
                self.alleleFreqDict[loc_name][orig_name] = {'frequency': freq_list} # add info to returned dict

        return self.alleleFreqDict
    
       

def run_simulation(knownSampleFile, useStrider, subpop, num_siblings, dest_path):
    # generate samples
    if useStrider:
        sibGen = StriderSibulator(knownSampleFile = knownSampleFile, subpop = subpop)
    else:
        sibGen = NistSibulator(knownSampleFile = knownSampleFile, subpop = subpop)
    
    
    sibling_samples = sibGen.generate_sibling_samples(num_siblings)

    # write output to file
    output = {}
    col_rename = {}
    for i, sample in enumerate([sibGen.knownSample] + sibling_samples): # output first row is known sample
        for k, v in sample.items():
            if k=='numReplaced':
                output[k] = v
                continue
     
            output[k] = output.get(k, []) + [float(v[0])] # bit hacky that some are single item lists
            output[k + '.1'] = output.get(k + '.1', []) + [float(v[1])] # bit hacky that some are single item lists
            if i == 0:
                col_rename[k + '.1'] = k # for renaming columns from key.1 to key to preserve same format as input
    output = pd.DataFrame(output)
    output = output.rename(col_rename, axis=1)
    output.to_csv(dest_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = 'Evaluate a single trained model.')
    parser.add_argument('--known_sample_file', dest='known_sample_file', required=True,
                        help='path to known sample profile for which siblings are based on')
    parser.add_argument('--use-strider', dest='useStrider', action='store_true', help='use Strider, default is NIST')
    parser.add_argument('--subpop', dest='subpop', required=True,
                        help='The population of interest to sample from')
    parser.add_argument('--num_siblings', dest='num_siblings',type=int, required=True,
                        help='number of sibling samples to generate')
    parser.add_argument('--dest_path', dest='dest_path', required=True,
                        help='file to output generated samples to')
    args = parser.parse_args()

    run_simulation(args.known_sample_file, args.useStrider, args.subpop, args.num_siblings, args.dest_path)

    
