# Print and Summary Methods for DNS Objects

Compact displays for the objects the fitting, estimation and simulation
functions return.

## Usage

``` r
# S3 method for class 'ns_fit'
print(x, ...)

# S3 method for class 'dns_fit'
print(x, ...)

# S3 method for class 'dns_fit'
summary(object, ...)

# S3 method for class 'summary.dns_fit'
print(x, ...)

# S3 method for class 'dns_dynamics'
print(x, ...)

# S3 method for class 'dns_simulation'
print(x, ...)

# S3 method for class 'dns_state_space'
print(x, ...)

# S3 method for class 'dns_state_space'
summary(object, ...)

# S3 method for class 'summary.dns_state_space'
print(x, ...)
```

## Arguments

- x:

  An object returned by the corresponding dnsr function.

- ...:

  Ignored, present for consistency with the generic.

## Value

`x`, invisibly. Called for the display it produces.

## Details

These are classes layered onto the lists the package already returned:
every element is still present under the same name, and
[`is.list()`](https://rdrr.io/r/base/list.html) and `$` access behave as
before. Only the default display changes.
