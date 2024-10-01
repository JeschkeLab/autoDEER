__version__ = "0.9"

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
