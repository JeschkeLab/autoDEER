# Packaging and Building Guide for autoDEER


This document provides instructions on how to package and build the autoDEER software for distribution. It covers the necessary steps to create source distributions and binary wheels, ensuring compatibility across different platforms. To enable a smooth build process, a docker container should be used.

## Prerequisites
Before starting the build process, ensure you have the following prerequisites installed:
- Docker: Install Docker from [here](https://www.docker.com/get-started).
- Git: Install Git from [here](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git).

## Cloning the Repository
First, clone the autoDEER repository from GitHub:

```bash
git clone www.github.com/jeschkelab/autoDEER.git
cd autoDEER
```

## Building with the Docker Container
Build the Docker container using the provided Dockerfile:
```bash
docker docker build --platform linux/amd64 --build-arg ENV=development -t autodeer:dev .
```
or interatively
```bash
docker run -it -v $(pwd):/app autodeer:dev                                       
```
This command creates a Docker image named `autodeer-build` that contains all the necessary dependencies for building autoDEER.
## Running the Build Process
Run the Docker container and mount the autoDEER source code directory:
```bash
docker run --rm -v $(pwd):/app autodeer-build bash -c "
    cd /app &&
    pip install --upgrade build &&
    python -m build
"
```
This command mounts the current directory (the autoDEER source code) to `/app` in the Docker container, installs the build tool, and runs the build process.

## Building locally without Docker
The build process is first done with PyInstaller to create a standalone executable.

```bash
pip install pyinstaller
pyinstaller --onefile --name autodeer_gui autodeer/gui/main.py
```
This will create a standalone executable named `autodeer_gui` in the `dist` directory.
