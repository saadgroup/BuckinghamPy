# BuckinghamPi
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/mk-95/BuckinghamPi/master?filepath=examples.ipynb)

Python code that implement the Buckingham-Pi theorem for different variables and return all possible dimensionless pi terms

## Installation
---
clone the package into a directory
```buildoutcfg
git clone https://github.com/mk-95/BuckinghamPi.git . 
python setup.py install
```

Now you can import the module and use it as follow
```buildoutcfg
from buckinghampi import BuckinghamPi

Examp = BuckinghamPi()
Examp.add_variable(name='u', expression='l/t')
Examp.add_variable(name='rho', expression='m/(l**3)')
Examp.add_variable(name='mu', expression='m/(t*l)')
Examp.add_variable(name='dx', expression='l')
Examp.add_variable(name='dt', expression='t', select=True)

Examp.generate_pi_terms()

for space in Examp.pi_terms:
    print('-------------------------')
    for term in space:
        print(term)
```

---
## See Also

* [Documentation](https://htmlpreview.github.io/?https://github.com/mk-95/BuckinghamPi/blob/master/doc/buckinghampi.m.html)
--- 
