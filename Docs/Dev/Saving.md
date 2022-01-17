# A standardised data format for multi-step EPR experiments

As EPR experiments move beyond single pulse sequences and become a train of sub-experiments with data processing between, it become critical that there is a single way of saving such data. Here I present a comprehensive, hierarchical data format for this task.

The criteria for this format are:
- [ ] The spectrometer's specifications must be saved
- [ ] Each sub-experiment must be saved
        - [ ] The raw data
        - [ ] Any processed data
        - [ ] The pulse sequence
        - [ ] The sequence input parameters

## Implementation:
The proposed implementation for this is to use an open standard, scientific data format, such as HDF5. HDF5 is a widely used data format for the saving of data-sets and metadata, and is even the backbone and basis of MATLAB's ".mat" file format. One of the main benefits HDF5 has over older more traditional format is that includes data-set chunking, this is the idea that a dataset it broken into smaller chunks to reduce the file size and ease compression. This was particularly a problem in astronomy due to their large file sizes, but in EPR where our file sizes are currently still small is of less use.

Whilst the HDF5 file will be produced by Python it is a completely open standard and subsequently can be opened by many other programming languages such as MATLAB as well as graphical interfaces such as HDFview. This means that while an HDF5 file is not naively human readable and can quickly and easily be read by a human. It should also be simple and easy to write a small extension to both "eprload" and "deerload" to read these files. 

### The structure of a HDF5 file
A HDF5 file is formed of 3 main components: groups, datasets and attributes.
#### Groups
The group is what makes this a hierarchical data format, to the extent where one could even consider this to be a saved folder not a single file. Every file starts with the root group "/". In this we can create new groups, and inside more sub groups etc…
#### Datasets
Datasets are used for the saving of large array. They are in many ways similar to a numpy array, in that they have a size,shape,data type and num bytes attributes. They can also contain attributes
#### Attributes
Attributes are the way of adding metadata to a HDF5 file, they are used for either single values, strings or small arrays. An attribute can be added to either a group of a dataset. There is no ability to partially read an attribute. You must always read all of it.

### Our implementation
One of other changes that would be beneficial to add is saving as we go. Currently if the software or spectrometer crashes all unsaved data is lost. We propose that the save file is constantly generated in the tmp folder until the end of the experiment. It is only at the end of the experiment that the user specifies a name and location and the file is then simply moved to this location. We could also extend this further such that during a very long experiment the data is autosaved once an hour. This would have a few benefits: Firstly, if the software crashes it can be recovered. Secondly, if there is a fundamental change in the sample environment (such as a change in temperature), then the only a short period of data is lost. 

#### Current Limitations
- Only 1 type of experiment may feature i