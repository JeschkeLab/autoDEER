"""
This is a class designed to make the saving of data and meta data as simple as possible. 


2022 Hugo Karas
"""


import os
import h5py as h
import logging

data_dir = None

log = logging.getLogger('core.Saving')

class save_file:

    def create_file(self,file_name):
        """ This function create a file with this file name. If the file does not exist it will return an error"""
        if not os.path.isabs(file_name): # path is not an absolute path
            if self.data_dir == None:
                raise RuntimeError("Data directory must be set unless an absolute path is given")
            else:
                split_name = os.path.split(file_name)
                if split_name[1] == []:
                    raise RuntimeError("A file must be given not a directory")
                else:
                    split_name_no_ext = os.path.splitext(split_name[1])[0]
                    file_path = self.data_dir + split_name_no_ext
                    self.file = h.File(file_path.encode(encoding="ascii"),'w')
                    return self.file

        else:
            file_path = file_name
            self.file = h.File(file_path.encode(encoding="ascii"),'w')
            return self.file

    def open_file(self,file_name):
        """This function will open a file with the given name/path. If one does not exist it will create one."""
        if not os.path.isabs(file_name): # path is not an absolute path
            if self.data_dir == None:
                raise RuntimeError("Data directory must be set unless an absolute path is given")
            else:
                split_name = os.path.split(file_name)
                if split_name[1] == []:
                    raise RuntimeError("A file must be given not a directory")
                else:
                    split_name_no_ext = os.path.splitext(split_name[1])[0]
                    file_path = self.data_dir + split_name_no_ext
                    self.file = h.File(file_path.encode(encoding="ascii"),'a')
                    return self.file

        else:
            file_path = file_name
            self.file = h.File(file_path.encode(encoding="ascii"),'r+')
            return self.file

    def set_data_directory(self,directory):

        if  not os.path.isdir(directory):
            print("Can't find directory")
            answer = input(f"Would you like to create a directory at {directory}: (Y/N)")
            if (answer == "Y") or (answer == 'y'):
                try: 
                    os.mkdir(directory)
                except:
                    raise RuntimeError("Can't make directory")
                else:
                    self.data_dir = directory              
            else:
                RuntimeError("Directory Error")
        else:
            self.data_dir = directory

    def save_spectrometer_config(self,api):
        # This function is designed to save the spectrometer configuration.
        # For a spectrometer built with multiple hardware components this will 
        # need to be run once for each component.

        file = self.file
        meta = api.meta

        grp = file.create_group("Spectrometer")
        
        for key in meta:
            data = meta[key]
            grp.attrs.create(key,data)
        pass

    def save_experiment_meta(self,exp):
        # This doesn't work yet, no such thing as experimental metadata yet
        # This is a function that saves an experiment to this save file
        file = self.file

        list_grps = file.keys()

        if exp.type in list_grps:
            cur_grp = file[exp.type]
        else:
            cur_grp = file.create_group(exp.type)
        

    def save_experimental_data(self,data,grp_name:str,dset_name:str=None,meta:dict=None) -> h.Dataset:

        file = self.file
        list_grps = file.keys()

        if grp_name not in list_grps:
            cur_grp = file.create_group(grp_name)
        else:
            cur_grp = file[grp_name]

        if dset_name == None:
            dset_name = "raw_data"

        dset_name = self._find_unique_name(cur_grp,dset_name) # Check if dset name already exists, add a suffix if needed

        if dset_name == 1: #ERROR: We are no unique group can be found.
            return 0

        dset_grp = cur_grp.create_group(dset_name)
        # dset = cur_grp.create_dataset(dset_name,data=[data.time,data.data]) # This attempts to place it all in one dataset
        dset_grp.create_dataset('time',data=data.time)
        dset_grp.create_dataset('data',data=data.data)

        if meta != None:
            for key in meta:
                cur_grp.attrs.create(key,meta[key])

        return dset_grp
        

    def autosave_experiment_data(self,exp) -> None:
        # This does not work yet
        file = self.file
        list_grps = file.keys()
        if exp.type in list_grps:
            cur_grp = file[exp.type]
        else:
            print("Metadata must have been saved before the autosaving data")
        
        if "autosave" not in cur_grp.keys():
            cur_grp.create_group("autosave")

        pass

    def save_attribute(self,data,name,grp_name) -> None:
        file = self.file
        list_grps = file.keys()

        if grp_name not in list_grps:
            cur_grp = file.create_group(grp_name)
        else:
            cur_grp = file[grp_name]
        
        cur_grp.attrs.create(name,data)

    def _find_unique_name(self,cur_grp,dset_name):

        if dset_name in cur_grp.keys():
            i=1
            tmp_dset_name = dset_name
            while (tmp_dset_name in cur_grp.keys()) and i <10:
                tmp_dset_name = dset_name +f'({i})'
                i = i + 1
            
            if i<10:
                log.warning(f'Can\'t find a unique group, using {tmp_dset_name} instead')
                return tmp_dset_name
            else:
                log.error('Can\'t find a unique group')
                return 1 
        else:
            return dset_name
