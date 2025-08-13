import importlib.metadata
import toml
import importlib.util
import os.path
MODULE_DIR = importlib.util.find_spec('pyepr').submodule_search_locations[0]

# with open(os.path.join(MODULE_DIR,"../pyproject.toml"), "r") as f:
#     config = toml.load(f)
#     __version__ = config["tool"]["poetry"]["version"]

__version__ = importlib.metadata.version('autoDEER')

try: 
    import git
    import os

    branch = git.Repo(os.path.dirname(__file__),search_parent_directories=True).active_branch.name
    commit = git.Repo(os.path.dirname(__file__),search_parent_directories=True).head.commit.hexsha

    if branch == 'main':
        __version__ = __version__
    else:
        __version__ = f"{__version__}-{branch}-{commit[:7]}"

except:
    __version__ = __version__
