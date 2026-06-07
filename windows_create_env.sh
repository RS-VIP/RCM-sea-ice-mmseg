conda create -n env_mmseg python=3.10 -y
conda activate env_mmseg

conda install -c conda-forge numpy matplotlib scipy packaging prettytable xarray h5netcdf scikit-learn scikit-image tqdm joblib pyyaml pre-commit ipywidgets pandas requests -y
python -m pip install --upgrade pip setuptools wheel

pip install torch torchvision torchmetrics --index-url https://download.pytorch.org/whl/cu118
# check: python -c "import torch; print(torch.__version__); print(torch.version.cuda)"

pip install wandb icecream termcolor torch-summary rasterio lxml geopandas shapely ftfy

pip install -U openmim
mim install mmengine


#====================== INSTALL MMCV LIBRARY =========================


#============ for windows
# 1- 
# Python: 3.10
# PyTorch: 2.7.1+cu118
# CUDA Toolkit: 11.8
# GPU arch: 6.1 for GTX 1070
# Visual Studio: VS 2022 with MSVC v142 / 14.29 toolset

# 2- Install NVIDIA CUDA Toolkit 11.8.
# check: nvcc --version     #  release 11.8 

# 3- Install Visual Studio compiler
# In Visual Studio Installer, install:
# - Desktop development with C++
# - MSVC v142 - VS 2019 C++ x64/x86 build tools
# - Windows 10 SDK
# - C++ CMake tools

# 4- Open the correct terminal
# x64 Native Tools Command Prompt for VS 2022
call "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Auxiliary\Build\vcvars64.bat" -vcvars_ver=14.29

conda activate env_mmseg
cd ..
git clone https://github.com/open-mmlab/mmcv.git
cd /mmcv
# checkout the correct version compatible with mmsegmentation
git fetch --tags
git checkout v2.1.0

# remove build folder if reinstaling 
rmdir /s /q build

# check compiler
where cl
cl
# mmcv around 19.29

# check GPU capability

# 5- Set build variables
set CUDA_HOME=C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v11.8
set CUDA_PATH=C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v11.8
set PATH=C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v11.8\bin;%PATH%

set DISTUTILS_USE_SDK=1
# GPU capability found in step 4
set TORCH_CUDA_ARCH_LIST=6.1  
set MMCV_WITH_OPS=1
set MAX_JOBS=1

# 6- Install MMCV:
python setup.py build_ext
python setup.py develop

# 7- Check installation:
python .dev_scripts/check_installation.py
#============ ===========

#============ for linux
# install the package from source code:
# https://mmcv.readthedocs.io/en/latest/get_started/build.html (branch 2.x)
#============ =========

#====================== ==================== =========================

# INSTALL MMSEGMENTATION LIBRARY
cd d:/Temp/sea-ice-mmseg
pip install -v -e .

# INSTALL MAGIC LIBRARY
# 1- Have the magic_lib folder
# 2- Code edits in: 
#   CMakeLists.txt:
#     change cmake_minimum_required(VERSION 3.4) to cmake_minimum_required(VERSION 3.15)
#     delete the lines: set(PYBIND11_PYTHON_VERSION 3.8 CACHE STRING "") and set(PYBIND11_PYTHON_VERSION 3.8)
#   magic_core/CMakeLists.txt
#     change cmake_minimum_required(VERSION 3.4) to cmake_minimum_required(VERSION 3.15)
#   pyproject.toml
#     change the contents to:
#       [build-system]
#       requires = ["setuptools", "wheel", "cmake", "ninja"]
#       build-backend = "setuptools.build_meta"
# 3- Go to root directory of magic_lib and run
pip install -e .