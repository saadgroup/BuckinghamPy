# BuckinghamPy

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/saadgroup/BuckinghamPi/master?filepath=buckinghampy-gui.ipynb) (GUI App)

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/mk-95/BuckinghamPi/master?filepath=examples.ipynb) (Script)

Python code that implements the Buckingham-Pi theorem and return all sets of dimensionless groups

## Installation
---
Clone the package from the github repository into the current directory
```buildoutcfg
git clone https://github.com/saadgroup/BuckinghamPi.git . 
```
Use `pip` tool to install the package in the active python evironment
```buildoutcfg
pip install .
```
## Example
Now you can import the module and use it as follows
```buildoutcfg
from buckinghampy import BuckinghamPi

Example = BuckinghamPi()
Example.add_variable(name='u', expression='l/t')
Example.add_variable(name='rho', expression='m/(l**3)')
Example.add_variable(name='mu', expression='m/(t*l)')
Example.add_variable(name='dx', expression='l')
Example.add_variable(name='dt', expression='t', explicit=True)

Example.generate_pi_terms()

Example.print_all()
```
or you can import the graphic user interface only in a Jupyter cell
```buildoutcfg
from buckinghampy import BuckinghamPiGui

GUI=BuckinghamPiGui()
```

---
## See Also

* [Documentation](https://htmlpreview.github.io/?https://github.com/mk-95/BuckinghamPi/blob/master/doc/buckinghampi.m.html)
--- 
