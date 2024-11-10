## FastCCI

## installation
```shell
conda create -n fastcci python=3.11
conda activate fastcci
cd ~/Softwares/
mkdir poetry_fastcci
cd poetry_fastcci/
python3 -m venv .
./bin/pip install -U pip setuptools
./bin/pip install poetry
ln -s /net/mulan/home/siyuh/Softwares/poetry_fastcci/bin/poetry /net/mulan/home/siyuh/anaconda3/envs/fastcci/bin/
cd ~/Projects/FastCCI/code/FastCCI/

poetry init
poetry add numpy@^1.24
poetry add pandas
poetry add scipy
poetry add scanpy
poetry add loguru

poetry add ipykernel -G dev
poetry add matplotlib -G dev
poetry add seaborn -G dev
```