# Code to run on a fresh Ubuntu machine to get it to run this stuff
# Run on new server boot

sudo apt-get update

# Install basics
sudo apt install python
sudo apt install libfontconfig1
sudo apt-get install libgomp1
sudo apt-get install gcc

# Install a sane shell
sudo apt install zsh
sudo chsh -s /bin/zsh
sudo chsh -s /usr/bin/zsh ubuntu

sh -c "$(curl -fsSL https://raw.github.com/ohmyzsh/ohmyzsh/master/tools/install.sh)"

# Install conda
# Currently using sha256 hash 1314b90489f154602fd794accfc90446111514a5a72fe1f71ab83e07de9504a7
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh

# Add conda to path
export PATH="$PATH:home/ubuntu/miniconda3/bin"

# Conda stuff
conda create --name qosf
conda activate qosf

conda install numpy
conda install matplotlib
conda install python=3.7.1
conda install psi4 -c psi4
python -m pip install --user openfermion  # Roy went to sleep here
python -m pip install --user openfermionpyscf
pip install qiskit
