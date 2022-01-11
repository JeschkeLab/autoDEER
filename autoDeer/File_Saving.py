import os
import h5py as h

data_dir = None

def create_file(file_name):
    global data_dir
    if not os.path.isabs(file_name): # path is not an absolute path
        if data_dir == None:
            raise RuntimeError("Data directory must be set unless an absoluet path is given")
        else:
            split_name = os.path.split(file_name)
            if split_name[1] == []:
                raise RuntimeError("A file must be given not a directory")
            else:
                split_name_no_ext = os.path.splitext(split_name[1])[0]
                file_path = data_dir + split_name_no_ext
                file = h.file(file_path,'w')
                return file

    else:
        file_path = file_name
        file = h.file(file_path,'w')
        return file

def set_data_directory(directory):

    if  not os.path.isdir(directory):
        print("Can't find directory")
        answer = input(f"Would you like to create a directory at {directory}: (Y/N)")
        if (answer == "Y") or (answer == 'y'):
            try: 
                os.mkdir(directory)
            except:
                raise RuntimeError("Can't make directory")
            else:
                global data_dir
                data_dir = directory              
        else:
            RuntimeError("Directory Error")
    else:
        #global data_dir
        data_dir = directory

