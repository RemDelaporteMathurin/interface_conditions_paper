# interface_conditions_paper


Scripts for the article [Rémi Delaporte-Mathurin et al 2021 Nucl. Fusion 61 036038](https://iopscience.iop.org/article/10.1088/1741-4326/abd95f/meta)

1. Run a FEniCS docker container

```
docker run -ti -v $(pwd):/home/fenics/shared --name fenics quay.io/fenicsproject/stable:latest
```
2. Install FESTIM 0.7.1

```
pip install git+https://github.com/RemDelaporteMathurin/FESTIM@0.7.1
```
