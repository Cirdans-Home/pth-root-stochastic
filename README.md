# Stochastic pth root approximation of a stochastic matrix: A Riemannian optimization approach

Algorithms for approximating a stochastic pth root of a stochastic matrix using Riemannian optimization. 

## Manopt

The code uses the optimization algorithms on Riemannian manifold from the 
[Manopt library](https://www.manopt.org). The library is included as a
git submodule. The code can be obtained from this repository by doing
```
git clone --recurse-submodules git@github.com:Cirdans-Home/approximate-embedded-markov.git
```
Older versions of gith may need to do:
```
git clone --recursive git@github.com:Cirdans-Home/approximate-embedded-markov.git
```
If none of this works
```
git clone git@github.com:Cirdans-Home/approximate-embedded-markov.git
cd approximate-embedded-markov
git submodule update --init --recursive
```
## Contributors

- F. Durastante [:link:](https://fdurastante.github.io/)
- B. Meini [:link:](https://people.dm.unipi.it/meini/)

## License

In accordance with the Manopt library, the new code is released under the
GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007 (GPL3).
