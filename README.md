# A Riemannian Approach to approximately embed Markov Chains

Algorithms for the approximation of an embedding for Markov chains.

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

## License

In accordance with the Manopt library, the new code is released under the
GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007 (GPL3).
