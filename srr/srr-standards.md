## [General Standards](https://stats-devguide.ropensci.org/standards.html#general-standards)

### Documentation

- [ ] **G1.0** *Statistical Software should list at least one primary reference from published academic literature.* 
- [ ] **G1.1** *Statistical Software should document whether the algorithm(s) it implements are:* - *The first implementation of a novel algorithm*; or - *The first implementation within **R** of an algorithm which has previously been implemented in other languages or contexts*; or - *An improvement on other implementations of similar algorithms in **R***. 
- [ ] **G1.2** *Statistical Software should include a* Life Cycle Statement *describing current and anticipated future states of development.* 

#### Statistical Terminology

- [ ] **G1.3** *All statistical terminology should be clarified and unambiguously defined.* 

#### Function-level Documentation

- [ ] **G1.4** *Software should use [`roxygen2`](https://roxygen2.r-lib.org/) to document all functions.*
 - [ ] **G1.4a** *All internal (non-exported) functions should also be documented in standard [`roxygen2`](https://roxygen2.r-lib.org/) format, along with a final `@noRd` tag to suppress automatic generation of `.Rd` files or [`@keywords internal`](https://roxygen2.r-lib.org/reference/tags-index-crossref.html?q=keywords%20internal#null) if documentation is still desired.* 

#### Supplementary Documentation

- [ ] **G1.5** *Software should include all code necessary to reproduce results which form the basis of performance claims made in associated publications.* 
- [ ] **G1.6** *Software should include code necessary to compare performance claims with alternative implementations in other R packages.* 

### Input Structures

#### Uni-variate (Vector) Input

- [ ] **G2.0** *Implement assertions on lengths of inputs, particularly through asserting that inputs expected to be single- or multi-valued are indeed so.*
 - [ ] **G2.0a** Provide explicit secondary documentation of any expectations on lengths of inputs
- [ ] **G2.1** *Implement assertions on types of inputs (see the initial point on nomenclature above).*
 - [ ] **G2.1a** *Provide explicit secondary documentation of expectations on data types of all vector inputs.*
- [ ] **G2.2** *Appropriately prohibit or restrict submission of multivariate input to parameters expected to be univariate.*
- [ ] **G2.3** *For univariate character input:*
 - [ ] **G2.3a** *Use `match.arg()` or equivalent where applicable to only permit expected values.*
 - [ ] **G2.3b** *Either: use `tolower()` or equivalent to ensure input of character parameters is not case dependent; or explicitly document that parameters are strictly case-sensitive.*
- [ ] **G2.4** *Provide appropriate mechanisms to convert between different data types, potentially including:*
 - [ ] **G2.4a** *explicit conversion to `integer` via `as.integer()`*
 - [ ] **G2.4b** *explicit conversion to continuous via `as.numeric()`*
 - [ ] **G2.4c** *explicit conversion to character via `as.character()` (and not `paste` or `paste0`)*
 - [ ] **G2.4d** *explicit conversion to factor via `as.factor()`*
 - [ ] **G2.4e** *explicit conversion from factor via `as...()` functions*
- [ ] **G2.5** *Where inputs are expected to be of `factor` type, secondary documentation should explicitly state whether these should be `ordered` or not, and those inputs should provide appropriate error or other routines to ensure inputs follow these expectations.* 
- [ ] **G2.6** *Software which accepts one-dimensional input should ensure values are appropriately pre-processed regardless of class structures.* 

#### Tabular Input

- [ ] **G2.7** *Software should accept as input as many of the above standard tabular forms as possible, including extension to domain-specific forms.* 
- [ ] **G2.8** *Software should provide appropriate conversion or dispatch routines as part of initial pre-processing to ensure that all other sub-functions of a package receive inputs of a single defined class or type.*
- [ ] **G2.9** *Software should issue diagnostic messages for type conversion in which information is lost (such as conversion of variables from factor to character; standardisation of variable names; or removal of meta-data such as those associated with [`sf`-format](https://r-spatial.github.io/sf/) data) or added (such as insertion of variable or column names where none were provided).* 
- [ ] **G2.10** *Software should ensure that extraction or filtering of single columns from tabular inputs should not presume any particular default behaviour, and should ensure all column-extraction operations behave consistently regardless of the class of tabular data used as input.* 
- [ ] **G2.11** *Software should ensure that `data.frame`-like tabular objects which have columns which do not themselves have standard class attributes (typically, `vector`) are appropriately processed, and do not error without reason. This behaviour should be tested. Again, columns created by the [`units` package](https://github.com/r-quantities/units/) provide a good test case.*
- [ ] **G2.12** *Software should ensure that `data.frame`-like tabular objects which have list columns should ensure that those columns are appropriately pre-processed either through being removed, converted to equivalent vector columns where appropriate, or some other appropriate treatment such as an informative error. This behaviour should be tested.* 

#### Missing or Undefined Values

- [ ] **G2.13** *Statistical Software should implement appropriate checks for missing data as part of initial pre-processing prior to passing data to analytic algorithms.*
- [ ] **G2.14** *Where possible, all functions should provide options for users to specify how to handle missing (`NA`) data, with options minimally including:*
 - [ ] **G2.14a** *error on missing data*
 - [ ] **G2.14b** *ignore missing data with default warnings or messages issued*
 - [ ] **G2.14c** *replace missing data with appropriately imputed values*
- [ ] **G2.15** *Functions should never assume non-missingness, and should never pass data with potential missing values to any base routines with default `na.rm = FALSE`-type parameters (such as [`mean()`](https://stat.ethz.ch/R-manual/R-devel/library/base/html/mean.html), [`sd()`](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/sd.html) or [`cor()`](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/cor.html)).*
- [ ] **G2.16** *All functions should also provide options to handle undefined values (e.g., `NaN`, `Inf` and `-Inf`), including potentially ignoring or removing such values.* 

### Algorithms

- [ ] **G3.0** *Statistical software should never compare floating point numbers for equality. All numeric equality comparisons should either ensure that they are made between integers, or use appropriate tolerances for approximate equality.* 
- [ ] **G3.1** *Statistical software which relies on covariance calculations should enable users to choose between different algorithms for calculating covariances, and should not rely solely on covariances from the `stats::cov` function.*
 - [ ] **G3.1a** *The ability to use arbitrarily specified covariance methods should be documented (typically in examples or vignettes).* 

### Output Structures

- [ ] **G4.0** *Statistical Software which enables outputs to be written to local files should parse parameters specifying file names to ensure appropriate file suffixes are automatically generated where not provided.* 

### Testing 

#### Test Data Sets

- [ ] **G5.0** *Where applicable or practicable, tests should use standard data sets with known properties (for example, the [NIST Standard Reference Datasets](https://www.itl.nist.gov/div898/strd/), or data sets provided by other widely-used R packages).*
- [ ] **G5.1** *Data sets created within, and used to test, a package should be exported (or otherwise made generally available) so that users can confirm tests and run examples.* 

#### Responses to Unexpected Input

- [ ] **G5.2** *Appropriate error and warning behaviour of all functions should be explicitly demonstrated through tests. In particular,*
 - [ ] **G5.2a** *Every message produced within R code by `stop()`, `warning()`, `message()`, or equivalent should be unique*
 - [ ] **G5.2b** *Explicit tests should demonstrate conditions which trigger every one of those messages, and should compare the result with expected values.*
- [ ] **G5.3** *For functions which are expected to return objects containing no missing (`NA`) or undefined (`NaN`, `Inf`) values, the absence of any such values in return objects should be explicitly tested.* 

#### Algorithm Tests

- [ ] **G5.4** **Correctness tests** *to test that statistical algorithms produce expected results to some fixed test data sets (potentially through comparisons using binding frameworks such as [RStata](https://github.com/lbraglia/RStata)).*
 - [ ] **G5.4a** *For new methods, it can be difficult to separate out correctness of the method from the correctness of the implementation, as there may not be reference for comparison. In this case, testing may be implemented against simple, trivial cases or against multiple implementations such as an initial R implementation compared with results from a C/C++ implementation.*
 - [ ] **G5.4b** *For new implementations of existing methods, correctness tests should include tests against previous implementations. Such testing may explicitly call those implementations in testing, preferably from fixed-versions of other software, or use stored outputs from those where that is not possible.*
 - [ ] **G5.4c** *Where applicable, stored values may be drawn from published paper outputs when applicable and where code from original implementations is not available*
- [ ] **G5.5** *Correctness tests should be run with a fixed random seed*
- [ ] **G5.6** **Parameter recovery tests** *to test that the implementation produce expected results given data with known properties. For instance, a linear regression algorithm should return expected coefficient values for a simulated data set generated from a linear model.*
 - [ ] **G5.6a** *Parameter recovery tests should generally be expected to succeed within a defined tolerance rather than recovering exact values.*
 - [ ] **G5.6b** *Parameter recovery tests should be run with multiple random seeds when either data simulation or the algorithm contains a random component. (When long-running, such tests may be part of an extended, rather than regular, test suite; see G5.10-4.12, below).* 
- [ ] **G5.7** **Algorithm performance tests** *to test that implementation performs as expected as properties of data change. For instance, a test may show that parameters approach correct estimates within tolerance as data size increases, or that convergence times decrease for higher convergence thresholds.*
- [ ] **G5.8** **Edge condition tests** *to test that these conditions produce expected behaviour such as clear warnings or errors when confronted with data with extreme properties including but not limited to:*
 - [ ] **G5.8a** *Zero-length data*
 - [ ] **G5.8b** *Data of unsupported types (e.g., character or complex numbers in for functions designed only for numeric data)*
 - [ ] **G5.8c** *Data with all-`NA` fields or columns or all identical fields or columns*
 - [ ] **G5.8d** *Data outside the scope of the algorithm (for example, data with more fields (columns) than observations (rows) for some regression algorithms)*
- [ ] **G5.9** **Noise susceptibility tests** *Packages should test for expected stochastic behaviour, such as through the following conditions:*
 - [ ] **G5.9a** *Adding trivial noise (for example, at the scale of `.Machine$double.eps`) to data does not meaningfully change results*
 - [ ] **G5.9b** *Running under different random seeds or initial conditions does not meaningfully change results* 

#### Extended tests

- [ ] **G5.10** *Extended tests should included and run under a common framework with other tests but be switched on by flags such as as a `<MYPKG>_EXTENDED_TESTS="true"` environment variable.* - The extended tests can be then run automatically by GitHub Actions for example by adding the following to the `env` section of the workflow: 
- [ ] **G5.11** *Where extended tests require large data sets or other assets, these should be provided for downloading and fetched as part of the testing workflow.*
 - [ ] **G5.11a** *When any downloads of additional data necessary for extended tests fail, the tests themselves should not fail, rather be skipped and implicitly succeed with an appropriate diagnostic message.*
- [ ] **G5.12** *Any conditions necessary to run extended tests such as platform requirements, memory, expected runtime, and artefacts produced that may need manual inspection, should be described in developer documentation such as a `CONTRIBUTING.md` or `tests/README.md` file.*


---

## [Regression and Supervised Learning Standards](https://stats-devguide.ropensci.org/standards.html#regression-and-supervised-learning)

### Input data structures and validation

- [ ] **RE1.0** *Regression Software should enable models to be specified via a formula interface, unless reasons for not doing so are explicitly documented.*
- [ ] **RE1.1** *Regression Software should document how formula interfaces are converted to matrix representations of input data.* 
- [ ] **RE1.2** *Regression Software should document expected format (types or classes) for inputting predictor variables, including descriptions of types or classes which are not accepted.* 
- [ ] **RE1.3** *Regression Software which passes or otherwise transforms aspects of input data onto output structures should ensure that those output structures retain all relevant aspects of input data, notably including row and column names, and potentially information from other `attributes()`.*
 - [ ] **RE1.3a** *Where otherwise relevant information is not transferred, this should be explicitly documented.* 
- [ ] **RE1.4** *Regression Software should document any assumptions made with regard to input data; for example distributional assumptions, or assumptions that predictor data have mean values of zero. Implications of violations of these assumptions should be both documented and tested.* 

### Pre-processing and Variable Transformation

- [ ] **RE2.0** *Regression Software should document any transformations applied to input data, for example conversion of label-values to `factor`, and should provide ways to explicitly avoid any default transformations (with error or warning conditions where appropriate).*
- [ ] **RE2.1** *Regression Software should implement explicit parameters controlling the processing of missing values, ideally distinguishing `NA` or `NaN` values from `Inf` values (for example, through use of `na.omit()` and related functions from the `stats` package).* 
- [ ] **RE2.2** *Regression Software should provide different options for processing missing values in predictor and response data. For example, it should be possible to fit a model with no missing predictor data in order to generate values for all associated response points, even where submitted response values may be missing.*
- [ ] **RE2.3** *Where applicable, Regression Software should enable data to be centred (for example, through converting to zero-mean equivalent values; or to z-scores) or offset (for example, to zero-intercept equivalent values) via additional parameters, with the effects of any such parameters clearly documented and tested.*
- [ ] **RE2.4** *Regression Software should implement pre-processing routines to identify whether aspects of input data are perfectly collinear, notably including:*
 - [ ] **RE2.4a** *Perfect collinearity among predictor variables*
 - [ ] **RE2.4b** *Perfect collinearity between independent and dependent variables* 

### Algorithms

- [ ] **RE3.0** *Issue appropriate warnings or other diagnostic messages for models which fail to converge.*
- [ ] **RE3.1** *Enable such messages to be optionally suppressed, yet should ensure that the resultant model object nevertheless includes sufficient data to identify lack of convergence.*
- [ ] **RE3.2** *Ensure that convergence thresholds have sensible default values, demonstrated through explicit documentation.*
- [ ] **RE3.3** *Allow explicit setting of convergence thresholds, unless reasons against doing so are explicitly documented.* 

### Return Results

- [ ] **RE4.0** *Regression Software should return some form of "model" object, generally through using or modifying existing class structures for model objects (such as `lm`, `glm`, or model objects from other packages), or creating a new class of model objects.*
- [ ] **RE4.1** *Regression Software may enable an ability to generate a model object without actually fitting values. This may be useful for controlling batch processing of computationally intensive fitting algorithms.* 

#### Accessor Methods

- [ ] **RE4.2** *Model coefficients (via `coef()` / `coefficients()`)*
- [ ] **RE4.3** *Confidence intervals on those coefficients (via `confint()`)*
- [ ] **RE4.4** *The specification of the model, generally as a formula (via `formula()`)*
- [ ] **RE4.5** *Numbers of observations submitted to model (via `nobs()`)*
- [ ] **RE4.6** *The variance-covariance matrix of the model parameters (via `vcov()`)*
- [ ] **RE4.7** *Where appropriate, convergence statistics* 
- [ ] **RE4.8** *Response variables, and associated "metadata" where applicable.*
- [ ] **RE4.9** *Modelled values of response variables.*
- [ ] **RE4.10** *Model Residuals, including sufficient documentation to enable interpretation of residuals, and to enable users to submit residuals to their own tests.*
- [ ] **RE4.11** *Goodness-of-fit and other statistics associated such as effect sizes with model coefficients.*
- [ ] **RE4.12** *Where appropriate, functions used to transform input data, and associated inverse transform functions.* 
- [ ] **RE4.13** *Predictor variables, and associated "metadata" where applicable.* 

#### Prediction, Extrapolation, and Forecasting

- [ ] **RE4.14** *Where possible, values should also be provided for extrapolation or forecast *errors*.*
- [ ] **RE4.15** *Sufficient documentation and/or testing should be provided to demonstrate that forecast errors, confidence intervals, or equivalent values increase with forecast horizons.* 
- [ ] **RE4.16** *Regression Software which models distinct responses for different categorical groups should include the ability to submit new groups to `predict()` methods.* 

#### Reporting Return Results

- [ ] **RE4.17** *Model objects returned by Regression Software should implement or appropriately extend a default `print` method which provides an on-screen summary of model (input) parameters and (output) coefficients.*
- [ ] **RE4.18** *Regression Software may also implement `summary` methods for model objects, and in particular should implement distinct `summary` methods for any cases in which calculation of summary statistics is computationally non-trivial (for example, for bootstrapped estimates of confidence intervals).* 

### Documentation

- [ ] **RE5.0** *Scaling relationships between sizes of input data (numbers of observations, with potential extension to numbers of variables/columns) and speed of algorithm.* 

### Visualization

- [ ] **RE6.0** *Model objects returned by Regression Software (see* **RE4***) should have default `plot` methods, either through explicit implementation, extension of methods for existing model objects, or through ensuring default methods work appropriately.*
- [ ] **RE6.1** *Where the default `plot` method is **NOT** a generic `plot` method dispatched on the class of return objects (that is, through an S3-type `plot.<myclass>` function or equivalent), that method dispatch (or equivalent) should nevertheless exist in order to explicitly direct users to the appropriate function.*
- [ ] **RE6.2** *The default `plot` method should produce a plot of the `fitted` values of the model, with optional visualisation of confidence intervals or equivalent.* 
- [ ] **RE6.3** *Where a model object is used to generate a forecast (for example, through a `predict()` method), the default `plot` method should provide clear visual distinction between modelled (interpolated) and forecast (extrapolated) values.* 

### Testing

#### Input Data

- [ ] **RE7.0** *Tests with noiseless, exact relationships between predictor (independent) data.*
 - [ ] **RE7.0a** In particular, these tests should confirm ability to reject perfectly noiseless input data.
- [ ] **RE7.1** *Tests with noiseless, exact relationships between predictor (independent) and response (dependent) data.*
 - [ ] **RE7.1a** *In particular, these tests should confirm that model fitting is at least as fast or (preferably) faster than testing with equivalent noisy data (see RE2.4b).* 

#### Return Results

- [ ] **RE7.2** Demonstrate that output objects retain aspects of input data such as row or case names (see **RE1.3**).
- [ ] **RE7.3** Demonstrate and test expected behaviour when objects returned from regression software are submitted to the accessor methods of **RE4.2**--**RE4.7**.
- [ ] **RE7.4** Extending directly from **RE4.15**, where appropriate, tests should demonstrate and confirm that forecast errors, confidence intervals, or equivalent values increase with forecast horizons.
