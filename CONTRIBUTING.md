# Contributing to the prolfqua package

## Code of conduct

Please note that this package is released with a [Contributor Code of
Conduct](https://fgcz.github.io/prolfqua/CODE_OF_CONDUCT.html). By
contributing to the package, you agree to abide by its terms.

## How you can contribute

There are several ways you can contribute to this project. If you want
to know more about why and how to contribute to open source projects
like this one, see this [Open Source
Guide](https://opensource.guide/how-to-contribute/).

### Share the work

Think prolfqua is useful? Let others discover it, by telling them in
person, via Twitter or a blog post.

Using prolfqua for a paper you are writing? Retrieve all citation
information by running `citation("prolfqua")` in the R console.

### Ask a question

Using prolfqua and got stuck? Browse the
[documentation](https://fgcz.github.io/prolfqua/) to see if you can find
a solution. Still stuck? Post your question as an [issue on
GitHub](https://github.com/fgcz/prolfqua/issues).

Want to ask a question in private? Contact the package maintainer by
[mail](mailto:wew@fgcz.ethz.ch).

### Propose an idea

Have an idea for a new prolfqua feature? Take a look at the
[documentation](https://fgcz.github.io/prolfqua/) and [issue
list](https://github.com/fgcz/prolfqua/issues) to see if it isn’t
included or suggested yet. If not, suggest your idea as an [issue on
GitHub](https://github.com/fgcz/prolfqua/issues). While we can’t promise
to implement your idea, it helps to:

- Explain in detail how it would work.
- Keep the scope as narrow as possible.

See below if you want to contribute code for your idea as well.

### Get support for your metabolic cart

Do you want your metabolic cart to be supported by the prolfqua package?
Consider sending me example raw data files via
[mail](mailto:wew@fgcz.ethz.ch).

### Report a bug

Using prolfqua and discovered a bug? That’s annoying! Don’t let others
have the same experience and report it as an [issue on
GitHub](https://github.com/fgcz/prolfqua/issues) so we can fix it. A
good bug report makes it easier for us to do so, so please include:

- Any details about your local setup that might be helpful in
  troubleshooting.
- Detailed steps to reproduce the bug.

### Improve the documentation

Noticed a typo on the website? Think a function could use a better
example? Good documentation makes all the difference, so your help to
improve it is very welcome!

#### The website

[This website](https://fgcz.github.io/prolfqua/) is generated with
[`pkgdown`](http://pkgdown.r-lib.org/). That means we don’t have to
write any html: content is pulled together from documentation in the
code, vignettes,
[Markdown](https://guides.github.com/features/mastering-markdown/)
files, the package `DESCRIPTION` and `_pkgdown.yml` settings. If you
know your way around `pkgdown`, you can [propose a file
change](https://help.github.com/articles/editing-files-in-another-user-s-repository/)
to improve documentation. If not, [report an
issue](https://github.com/fgcz/prolfqua/issues) and we can point you in
the right direction.

#### Function documentation

Functions are described as comments near their code and translated to
documentation using [`roxygen2`](https://klutometis.github.io/roxygen/).
If you want to improve a function description:

1.  Go to `R/` directory in the [code
    repository](https://github.com/fgcz/prolfqua).
2.  Look for the file with the name of the function.
3.  [Propose a file
    change](https://help.github.com/articles/editing-files-in-another-user-s-repository/)
    to update the function documentation in the roxygen comments
    (starting with `#'`).

### Contribute code

Care to fix bugs or implement new functionality for prolfqua? Awesome!
Have a look at the [issue list](https://github.com/fgcz/prolfqua/issues)
and leave a comment on the things you want to work on. See also the
development guidelines below.

## Development guidelines

We try to follow the [GitHub
flow](https://guides.github.com/introduction/flow/) for development.

1.  Fork [this repo](https://github.com/fgcz/prolfqua) and clone it to
    your computer. To learn more about this process, see [this
    guide](https://guides.github.com/activities/forking/).
2.  If you have forked and cloned the project before and it has been a
    while since you worked on it, [pull changes from the original
    repo](https://help.github.com/articles/merging-an-upstream-repository-into-your-fork/)
    to your clone by using `git pull upstream master`.
3.  Open the RStudio project file (`.Rproj`).
4.  Make your changes:
    - Write your code.
    - Test your code (bonus points for adding unit tests).
    - Document your code (see function documentation above).
    - Check your code with
      [`devtools::check()`](https://devtools.r-lib.org/reference/check.html)
      and aim for 0 errors and warnings.
5.  Commit and push your changes.
6.  Submit a [pull
    request](https://guides.github.com/activities/forking/#making-a-pull-request).
