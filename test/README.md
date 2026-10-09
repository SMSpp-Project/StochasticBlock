# test

The testers of the `StochasticBlock` module, which need nothing but the
module and the core SMS++ library.

- `StochasticBlock_test` exercises the core `StochasticBlock` machinery:
  wrapping a Block into a stochastic problem, setting the scenario data,
  re-solving under random modifications of the scenarios, and checking the
  results.

- `test_discrete` validates `DiscreteScenarioSet`, the class that manages a
  discrete set of scenarios: scenario storage and access, the random and the
  baseline pool selection, the configuration and serialization patterns, and
  the handling of invalid inputs.

- `MultiStageDiscreteScenarioSet_unit_test` checks the scenario tree of
  `MultiStageDiscreteScenarioSet` through its view API, the independence of
  the views and of their clones, the joint probabilities of the leaves and
  the serialization round trip.

- `IndependentMultiStageScenarioGenerator_unit_test` checks the same view API
  on `IndependentMultiStageScenarioGenerator`, whose stages do not depend on
  the history.

The reduction of a `DiscreteScenarioSet` to its representatives is asked of
a Solver that no dependency of this module provides: the heuristics of
`ScenarioReductionSolver` are tested in that module, and the reductions of
the instances of a model in the suites of the umbrella.

All of them are built by the provided `makefile` (or via CMake from the
umbrella, where each is registered as a separate `ctest` labelled
`StochasticBlock`). Run them as `./<name>`, with no argument.


## Authors

- **Rafael Durbano Lobato**  
  Dipartimento di Informatica  
  Università di Pisa

- **Benoît Tran**  
  Dipartimento di Informatica  
  Università di Pisa

- **Donato Meoli**  
  Dipartimento di Informatica  
  Università di Pisa


## License

This code is provided free of charge under the [GNU Lesser General Public
License version 3.0](https://opensource.org/licenses/lgpl-3.0.html),
see the [LICENSE](../LICENSE) file for details.
