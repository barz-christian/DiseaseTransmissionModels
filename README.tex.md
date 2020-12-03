Disease transition modeling in stan
================
creation date: 29-07-2020, update date: 2020-12-03

![\\begin{aligned}\\frac{S}{dt} = & -\\beta S \\frac{I}{N} \\\\
\\frac{dI}{dt} = & \\beta S\\frac{I}{N} - \\gamma I \\\\ \\frac{dR}{dt}
= & \\gamma I
\\end{aligned}](https://latex.codecogs.com/png.latex?%5Cbegin%7Baligned%7D%5Cfrac%7BS%7D%7Bdt%7D%20%3D%20%26%20-%5Cbeta%20S%20%5Cfrac%7BI%7D%7BN%7D%20%5C%5C%20%5Cfrac%7BdI%7D%7Bdt%7D%20%3D%20%26%20%5Cbeta%20S%5Cfrac%7BI%7D%7BN%7D%20-%20%5Cgamma%20I%20%5C%5C%20%5Cfrac%7BdR%7D%7Bdt%7D%20%3D%20%26%20%5Cgamma%20I%20%5Cend%7Baligned%7D
"\\begin{aligned}\\frac{S}{dt} = & -\\beta S \\frac{I}{N} \\\\ \\frac{dI}{dt} = & \\beta S\\frac{I}{N} - \\gamma I \\\\ \\frac{dR}{dt} = & \\gamma I \\end{aligned}")

where :

  - ![S(t)](https://latex.codecogs.com/png.latex?S%28t%29 "S(t)") is the
    number of people susceptible to becoming infected (no immunity),
  - ![I(t)](https://latex.codecogs.com/png.latex?I%28t%29 "I(t)") is the
    number of people currently infected (and infectious),
  - ![R(t)](https://latex.codecogs.com/png.latex?R%28t%29 "R(t)") is the
    number of recovered people (we assume they remain immune
    indefinitely),
  - ![N](https://latex.codecogs.com/png.latex?N "N") is the total number
    of the population,
  - ![\\beta](https://latex.codecogs.com/png.latex?%5Cbeta "\\beta") is
    the constant rate of infectious contact between people,
  - ![\\gamma](https://latex.codecogs.com/png.latex?%5Cgamma "\\gamma")
    the constant recovery rate of infected individuals.
