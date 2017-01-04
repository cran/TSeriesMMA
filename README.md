# Multiscale Multifractal Analysis
Multiscale multifractal analysis (MMA) (Gierałtowski et al., 2012) is a time series analysis method, designed to describe scaling properties of fluctuations within the signal analyzed. The main result of this procedure is the so called Hurst surface h(q,s) , which is a dependence of the local Hurst exponent h (fluctuation scaling exponent) on the multifractal parameter q (Kantelhardt et al., 2002) and the scale of observation s (data window width). 
The package contains main function mma() to generate the Hurst Surface for the given input signal.


## Usage

First we read data from an external file

```{r}
signal<-read.table('/home/data2.txt')
```

Now use function MMA to generate the hurst surface plot for the dataframe passed
Note that except the dataframe, all parameters are optional. Column name of the column to be analysed must be specified if the dataframe has column names .

```{r}
library('TSeriesMMA')
TSeriesMMA::mma(smin=10, smax=600, qmin=-5, qmax=5, data=signal, col=1)
```


