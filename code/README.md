Code
=============================


To recreate the figures used in the manuscript run the workflow.R file contained within the Multistrain folder.

```r
	setwd("Multistrain")
	source("workflow.R")
```


The *in silico* dilution methods were run in R (version 3.4.3), using the package deSolve (version 1.20). Additional packages used for plotting include lhs (version 0.14), fields (version 9.0) and viridis (version 0.4.0).


###References

R Core Team. 2015. R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.

Carnell R. 2016. lhs: Latin Hypercube Samples. R package version 0.14. https://CRAN.R-project.org/package=lhs

Soetaert K., Petzoldt T., Setzer R. W. 2010. Solving Differential Equations in R: Package deSolve. Journal of Statistical Software, 33(9), 1--25. [doi:10.18637/jss.v033.i09](http://dx.doi.org/10.18637/jss.v033.i09)

Nychka D., Furrer R., Paige J., Sain S. 2015. “fields: Tools for spatial data.” doi: 10.5065/D6W957CT (URL: http://doi.org/10.5065/D6W957CT), R package version 9.0

Garnier S. 2017. viridis: Default Color Maps from 'matplotlib'. R package version 0.4.0. https://CRAN.R-project.org/package=viridis
