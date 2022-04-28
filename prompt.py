import os

class Program:
    def __init__(self):
        self.settings_path=[]
        self.settings={'deltav': 0,'lag_max': 0, 'spectrum_cutoff': 0}
        self.data_path=[]
            
    def start(self):
        self.settings_get_path()        
    
    # Find where the settings file is.
    def settings_get_path(self):
        # Checks current directory for a settings file. If one is found, goes 
        # on to load the settings.
        print('Scanning %s for a settings.txt file.' % os.getcwd())
        if 'settings.txt' in os.listdir():
            print('Settings file found.')
            self.settings_path=os.getcwd()
            
            self.settings_get()
        # Otherwise, a prompt comes up to provide a valid path 
        else:
            print('Settings file not found. Please provide the directory of the settings file.')
            os.chdir(input('Settings directory path: '))
            print()
            self.settings_get_path()
    
    def settings_get(self):
        print('\nLoading settings.')
        with open('settings.txt') as the_file:
            lines=the_file.readlines()
        for line in lines:
            setting_values=line.split('#')[0].replace(' ','').replace('\t','').split('=')
            self.settings[setting_values[0]]=setting_values[1]
        self.settings_get_answer()
        
    def settings_get_answer(self):
        print('Review settings? [y/n]')
        answer=input()
        print()
        if answer=='y':
            self.settings_info()
            self.settings_review()
        elif answer=='n':
            self.data_get_path()
        else:
            print('Input not recognized. Please try again.')
            self.settings_get_answer()
        
    def settings_info(self):
        print('The following commands are recognized:')
        print('[disp] shows the current setting values.')
        print('Matching a setting key directly will allow changing its value.')
        print('[cont] will move on.')
        print('[help] shows this message again.')        
        
    def settings_review(self):
        answer=input('Please input a command: ')
        if answer=='disp':
            for key in self.settings:
                print('%s: %s' % (key,self.settings[key]))
            self.settings_review()
        elif answer in self.settings:
            new_val=input('Insert new value for %s: ' % answer)
            try:
                float(new_val)
                print('%s changed to %s.' % (answer,new_val))
                self.settings[answer]=new_val
            except ValueError:
                print('Entered value is not a number. Returning.')
            self.settings_review()        
        elif answer=='cont':
            if self.data_path==[]:
                self.data_get_path()
            else:
                self.data_ready()
        elif answer=='help':
            self.settings_info()
            self.settings_review()
        else:
            print('Input not recognized. Please try again.')
            self.settings_review()
        
    def data_get_path(self):
        print('Please provide the path of the data to be processed.')
        self.data_path=input('Data path: ')
        os.chdir(self.data_path)
        fits_num=0
        for file in os.listdir():
            if '.fits' in file:
                fits_num+=1
        if fits_num==0:
            print('No .fits files found. Returning.\n')
            self.data_get_path()
        else:
            print('Found %i .fits files. Is this the current directory? [y/n]' % fits_num)
            answer=input()
            print()
            if answer=='y':
                print('Data path set to %s.' % self.data_path)
                self.data_ready()
            elif answer=='n':
                self.data_get_path()
            else:
                print('Input not recognized. Please try again.')
                self.data_get_path()
    
    def data_ready(self):
        print('Start data processing [start], review settings [settings], or change data path [data]?')
        answer=input('Command: ')        
        if answer=='settings':
            self.settings_info()
            self.settings_review()
        elif answer=='data':
            self.data_get_path()
        elif answer=='start':
            print('Commencing!')
        
        else:
            print('Input not recognized. Please try again.')
            self.data_ready()
            
