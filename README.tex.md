Disease transition modeling in stan
================
creation date: 29-07-2020, update date: 2020-12-03

\(\begin{aligned}\frac{S}{dt} = & -\beta S \frac{I}{N} \\ \frac{dI}{dt} = & \beta S\frac{I}{N} - \gamma I \\ \frac{dR}{dt} = & \gamma I \end{aligned}\)

where :

  - \(S(t)\) is the number of people susceptible to becoming infected
    (no immunity),
  - \(I(t)\) is the number of people currently infected (and
    infectious),
  - \(R(t)\) is the number of recovered people (we assume they remain
    immune indefinitely),
  - \(N\) is the total number of the population,
  - \(\beta\) is the constant rate of infectious contact between people,
  - \(\gamma\) the constant recovery rate of infected individuals.
