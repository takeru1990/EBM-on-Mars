`test.c` is a C code made for a study of my master's thesis. It was last modified on 2015/01/06, according to file's metadata.

# Compile & Execute
```
$ gcc test.c -lm
$ ./a.out 
```
Then you would hopefuly get progress info in console during execution.
```
The 0 th loop begins.
The 1 th loop begins.
...
The 75 th loop begins.
```

When finished, you would get 7 files below with result info in console.
```
0.00756976(air), 0.0899085(ice), 0.0325217(rego), 149.484(T_sub), 76 loops, with no bug
```
* M_flat.dat
* T_flat.dat
* M_nega.dat
* T_nega.dat
* M_posi.dat
* T_posi.dat
* fCO2.dat

The following Abstract is from the thesis, which is available on [my blog](https://lookbackmargin.blog/2020/05/11/mars-co2-system-and-obliquity/).

# Abstract
The climate of Mars is dominated by CO2. This is because 95% of the atmosphere of Mars consists of CO2, and the ice-caps and surface regolith layers are large reservoirs of CO2. Exchange of CO2 among these three CO2 reservoirs on the Martian surface should therefore play a great role in determining the global climate of Mars.

On the other hand, Laskar et al. (2003) showed theoretically that orbital elements, especially obliquity and eccentricity, has changed dynamically on Mars for the last 10 million years. Because both the obliquity and eccentricity are critical parameters to determine distribution and seasonality of insolation, Martian climate system may have changed dynamically owing to the changes in these orbital parameters.

Analyses of behaviors of the Martian climate coupled with surface CO2 system have been studied by Nakamura and Tajika (2001, 2002, 2003). They focused on multiplicity of the Martian climate states and the evolution of the climate system of Mars, but did not investigate the details of behaviors of the Martian CO2 system against the orbital parameters such as obliquity and eccentricity, and possible relation to periglacial landforms observed on the Martian surface today.

In this study, we investigate the behaviors of the Martian surface CO2 system (partition of CO2 among the three CO2 reservoirs ??? atmosphere, ice caps, and regolith) against changes in obliquity and eccentricity. We found that the Martian atmosphere becomes very thin (down to ~0.3 mbar) when obliquity becomes lower because permanent ice cap develops significantly under low obliquity. An increase of obliquity would result in increases of the reservoir sizes of the atmosphere and regolith. However, we also found that the reservoir sizes of the Martian surface CO2 system would become constant against the obliquity changes when obliquity is larger than some critical value (~40??). This is because permanent ice cap cannot be stable under high obliquity and seasonal ice cap does not change a net CO2 budget through a year. Formation of CO2 ice on the polar-facing slope of crater is shown, and it may explain possible periglacial features as seen in the interiors of Martial craters located in the middle latitude.
