# Summarize values of a statistic

Calculates the estimate (mean or median), the confidence interval of the
estimate with the specified confidence level and the variance of the
data provided.

## Usage

``` r
summarizeStat(values, summary_method, conf_level = 0.95)
```

## Arguments

- values:

  Numeric, integer or logical vector.

- summary_method:

  Character, the measure of central tendency ("mean" or "median") to be
  used.

- conf_level:

  Numeric, confidence level of the interval (default: 0.95).

## Value

A data frame with 1 row and 4 columns:

- estimate:

  Numeric, the central tendency (mean or median) of 'values'.

- var:

  Numeric, the variance of 'values'.

- lwr:

  Numeric, the lower bound of the confidence interval.

- upr:

  Numeric, the upper bound of the confidence interval.
