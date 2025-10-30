# Use AlmaLinux 9 as base image
ARG ENV=production PYTHON_VERSION=3.12
FROM --platform=linux/amd64 almalinux:9

# Set environment variables
# ENV PYTHON_VERSION=3.12
ENV PYTHONUNBUFFERED=1
ENV DISPLAY=:1
ENV PYENV_ROOT="$HOME/.pyenv"
ENV PATH="$PYENV_ROOT/shims:$PYENV_ROOT/bin:$PATH"

# Set root password (for demonstration purposes)
RUN usermod root -p root

# Update system and install development tools
RUN dnf update -y && \
    dnf groupinstall "Development Tools" -y && \
    dnf install --allowerasing -y \
    gcc \
    gcc-c++ \
    openssl-devel \
    wget \
    make \
    sqlite-devel \
    readline-devel \
    tk-devel \
    python3-devel \
    epel-release \
    git \
    curl \
    cairo-devel \
    pkg-config

RUN /usr/bin/crb enable
# Install Qt6 and X11 dependencies for PyQt6
RUN dnf install -y \
    qt6-qtbase-devel \
    qt6-qttools-devel \
    mesa-libGL-devel \
    xcb-util \
    xcb-util-wm \
    xcb-util-image \
    xcb-util-keysyms \
    xcb-util-renderutil \
    libxkbcommon \
    libxkbcommon-x11 \
    libX11 \
    libXext \
    libXrender \
    fontconfig \
    dbus-x11 && \
    dnf clean all

# Create a working directory for applications
WORKDIR /app

# Create a non-root user for running applications
RUN useradd -m -s /bin/bash pyqtuser && \
    chown -R pyqtuser:pyqtuser /app

# Install pyenv
RUN git clone https://github.com/pyenv/pyenv.git $PYENV_ROOT && \
    cd $PYENV_ROOT && \
    src/configure && \
    make -C src

# Set up pyenv for the non-root user
RUN echo 'export PYENV_ROOT="$HOME/.pyenv"' >> /home/pyqtuser/.bashrc && \
    echo 'export PATH="$PYENV_ROOT/bin:$PATH"' >> /home/pyqtuser/.bashrc && \
    echo 'eval "$(pyenv init -)"' >> /home/pyqtuser/.bashrc

# RUN echo 'export PYENV_ROOT="$HOME/.pyenv"' >> ~/.bash_profile && \ 
#     echo '[[ -d $PYENV_ROOT/bin ]] && export PATH="$PYENV_ROOT/bin:$PATH"' >> ~/.bash_profile && \
#     echo 'eval "$(pyenv init - bash)"' >> ~/.bash_profile

# Install Python using pyenv
RUN pyenv install ${PYTHON_VERSION} && \
    pyenv global ${PYTHON_VERSION} && \
    pyenv rehash

# Switch to non-root user
USER pyqtuser

# Default command
RUN if [ "$ENV" = "production" ]; then \
    COPY --chown=pyqtuser:pyqtuser . /app; \
    fi

CMD ["/bin/bash", "pip install --upgrade pip setuptools wheel && pip install -e ."]
