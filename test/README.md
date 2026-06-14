# test

Two testers for the `StochasticBlock` module.

- `StochasticBlock_test` exercises the core `StochasticBlock` machinery:
  wrapping a Block into a stochastic problem, setting the scenario data,
  re-solving under random modifications of the scenarios, and checking the
  results.

- `test_discrete` validates `DiscreteScenarioSet`, the class that manages a
  discrete set of scenarios: scenario storage and access, pool selection,
  optimization-based scenario reduction, the configuration and serialization
  patterns, and the handling of invalid inputs. The reduction step relies on
  a `MILPSolver` backend when one is available.

Both are built by the provided `makefile` (or via CMake from the umbrella,
where each is registered as a separate `ctest`). Run them as
`./StochasticBlock_test` and `./test_discrete`.


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
