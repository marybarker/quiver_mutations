## Exchange number stratification

* Begin by running "Generate QPs from Triangulation" and exporting the list of QPs
* Using the page developer tools, paste in the QP list: `var qps = ...`
* Load the initial QP into the tool:
```
updateQPFromJSON(qps[0])
var qp = makeQP(QPglobalEdges, QPglobalNodes.map((i, idx) => idx), QPglobalFrozenNodes, QPglobalPotential, 'fromThing')
```
* Run the random search. Arguments are:
 * Base QP
 * Expected exchange number (potentials matching this exchange number will be displayed separately in the output)
 * maximum exchange number: stop searching for quivers after this many are found for a potential
 * Expected quivers (potentials producing this expected quiver set will be displayed separately in the output)
 * maxCycleLength for each term
 * minPotentialTerms
 * maxPotentialTerms
 * required terms to include
 * Number of potentials to generate
```
potentialRandomSearch(qp, 5, Infinity, [], 4, 1, 30, [], 1000000);
```

## Finding valid potentials

```
potentialStructuredTest(max)
```

By default this will try everything from (1, 1, 1) up to (max, max, max). The initial lines of the function can be edited to test a more specific family.

## Files

* `r123_1mil_2.json`: A test of exchange number stratification on (1, 2, 3) following the procedure above
* `structuredTest_8.json`: Finding potentials for (1, 1, 1) up to (8, 8, 8)
* `structuredTest_1_2_a_up_to_30.json`: finding potentials for (1, 2, a), a <= 30
* `structuredTest_1_a_a2_up_to_a_eq_5.json`: finding potentials for (1, a, a^2), a<=5
* `structuredTest_1_a_2a_up_to_a_eq_13.json`: finding potentials for (1, a, 2a), a <= 13