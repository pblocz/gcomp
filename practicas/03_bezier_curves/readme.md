# Bezier curves

This assigment consists on different implementation of
bezier curves computation algorithms.

## Install

Use `setup.py` to install the module `p3bezier` and compile
the extension `p3cbezier`.

## Structure

- The module `p3bezier.policies` contains the actual
  implementations of the algorithms.

- The class `p3bezier.Bezier` can be instantiated with a `Policy`
  to compute different curves from control polygons.

- `p3bezier.tester` contains a testting framework to test pieces
  of algorithms.


## Executables

- `p3bezier-tester`: tests the algorithms and compares them
  showing graphs of the algorithm execution.

- `p3bezier-plot`: interactive add dots to the control polygon
  and plot its respective bezier curve.
