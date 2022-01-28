# Guide for development


## Packaging 
### Installing the package in editable/develop mode
When actively devloping, it is best that the package is installed through pip in an editable enviroment. This allows you to edit the source code whilst still being able to use the standard import method.

#### Steps:
1) Navigate to the directory of your source code. This should be where your local git repository is stored.
2) `  python -m pip install -e .`



### Building a new release package
Whenever a new milestone has been reached and there is a complete functioning set of code, a complete package should be created. This should **only** be done from the main branch and **never** from the dev branch. So a new package can only be created after a merge request from dev to main, that has by definition passed all the tests. 

#### Required Files
The following files are needed for a a build to be possible
- `pyproject.toml`
- `setup.cfg`
- `requirments.txt`
- `MANIFEST.in`
- `__init__.py` in **every** module/submodule directory

#### Code:
